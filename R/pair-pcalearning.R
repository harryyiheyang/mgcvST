# PCAlearning saddlepoint path for the pair tests of large score dimension q.
#
# Each score covariance H_j is projected onto an r-dimensional orthonormal
# basis B of symmetric q x q matrices learned from training genes. With
# c_j = B' vech_w(H_j) (weighted vech, Frobenius isometric), the first two pair
# traces tr(H_i H_j) = c_i' c_j and tr((H_i H_j)^2) are those of the projected
# matrices, which the level-2 trace table of B gives by contraction of the
# monomials of c_i and c_j. The leading singular values of the pair spectrum
# come from a shared basis V (k leading eigenvectors of the training matrices):
# R_j = chol(V' H_j V / scale_j) and svd(R_i R_j'); the remainder is a
# Satterthwaite node matching the two traces. The matrices H_j come from the
# fit: reconstructed from the sparse INLA state, or built from the working
# model of an mgcv fit.

# Per-gene variance scales of the first-stage fit that define training strata.
# INLA: NB V = 1 / mu + 1 / theta, Poisson (and other families) V = 1 / mu, and
# mu_bar is the mean of the fitted mean, stored at estimation (a fit without it
# was estimated before it was stored and is refused). An mgcv fit is stratified
# by the field scale only.
.mgcvst_pca_scales <- function(fit) {
  G <- length(fit$feature_id)
  if (!identical(fit$estimator, "INLA")) {
    sigma_g2 <- as.numeric(fit$dispersion) / as.numeric(fit$lambda)
    return(data.frame(nb = rep(FALSE, G), theta = rep(NA_real_, G),
                      mu_bar = rep(NA_real_, G), sigma_g2 = sigma_g2,
                      sigma_e2 = rep(1, G), tau = 1 / sigma_g2))
  }
  nb <- fit$diagnostics$family_used == "negative_binomial"
  theta <- vapply(fit$family_parameters, function(x) {
    if (length(x)) x[1L] else NA_real_
  }, numeric(1L))
  mu_bar <- fit$mu_bar
  if (!is.numeric(mu_bar) || length(mu_bar) != G) {
    stop("The fit does not store mu_bar, the mean of each fitted mean; re-run ",
         "inlaST.estimate().", call. = FALSE)
  }
  mu_bar <- as.numeric(mu_bar)
  sigma_g2 <- as.numeric(fit$dispersion) /
    as.numeric(fit$smoothing_parameters[, fit$score_sparse$sp_index])
  sigma_e2 <- ifelse(nb, 1 + mu_bar / theta, 1)
  data.frame(nb = nb, theta = theta, mu_bar = mu_bar, sigma_g2 = sigma_g2,
             sigma_e2 = sigma_e2, tau = sigma_e2 / sigma_g2)
}

# 10 quantile bins of log sigma_g2 by 10 variance bins: non-NB genes form
# bin 1 and NB genes 9 quantile bins of log sigma_e2. Quantiles use `universe`.
.mgcvst_pca_cells <- function(scales, universe) {
  cell <- rep(NA_integer_, nrow(scales))
  lg <- log(scales$sigma_g2[universe])
  g <- findInterval(lg, stats::quantile(lg, seq(0, 1, 0.1)),
                    rightmost.closed = TRUE, all.inside = TRUE)
  e <- rep(1L, length(universe))
  nb <- scales$nb[universe]
  if (any(nb)) {
    le <- log(scales$sigma_e2[universe][nb])
    e[nb] <- 1L + findInterval(le, stats::quantile(le, seq(0, 1, length.out = 10)),
                               rightmost.closed = TRUE, all.inside = TRUE)
  }
  cell[universe] <- (g - 1L) * 10L + e
  cell
}

# Draw n_per_cell genes per cell from `eligible` with the current RNG state.
# Quota left by sparse cells is reallocated proportionally to cell sizes, with
# largest remainders, until n_per_cell * 100 genes (or every eligible gene)
# are drawn.
.mgcvst_pca_draw <- function(cell, eligible, n_per_cell) {
  total <- min(as.integer(n_per_cell) * 100L, length(eligible))
  N <- tabulate(cell[eligible], 100L)
  quota <- pmin(as.integer(n_per_cell), N)
  deficit <- total - sum(quota)
  while (deficit > 0) {
    room <- N - quota
    w <- ifelse(room > 0, N, 0)
    x <- deficit * w / sum(w)
    add <- pmin(floor(x), room)
    rem <- deficit - sum(add)
    if (rem > 0) {
      o <- order(-(x - floor(x)) * (room - add > 0))[seq_len(rem)]
      add[o] <- add[o] + 1L
    }
    add <- pmin(add, room)
    quota <- quota + add
    deficit <- total - sum(quota)
  }
  drawn <- integer()
  for (k in which(quota > 0)) {
    cand <- eligible[cell[eligible] == k]
    drawn <- c(drawn, cand[sample.int(length(cand), quota[k])])
  }
  drawn
}

# Stratified training genes; the caller's random-number state is restored.
.mgcvst_pca_training <- function(scales, universe, n_per_cell, seed) {
  cell <- .mgcvst_pca_cells(scales, universe)
  old_exists <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (old_exists) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit(if (old_exists) {
    assign(".Random.seed", old_seed, envir = .GlobalEnv)
  } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  set.seed(as.integer(seed))
  list(train = sort(.mgcvst_pca_draw(cell, universe, n_per_cell)), cell = cell)
}

# How the PCAlearning route obtains the curvature matrices of the used genes.
# `materialize(ids, B, V, pack)` returns, for genes `ids` in blocks of 32: the
# score coordinates a_j (q x n), the coefficients c_j = B' vech_w(H_j) when `B`
# is given, ||H_j||_F^2, the scale max|H_j|, with `pack = TRUE` the float32
# weighted-vech H_j and, when the shared basis `V` is given, the factors R_j
# (k^2 x n). Failed genes carry `error`.
.mgcvst_pca_producer <- function(fit, used, basis, threads) {
  if (identical(fit$score_backend, "sparse")) {
    fit <- .inlast_sparse_prepare(fit)
    g <- fit$score_sparse
    q <- as.integer(basis$rank)
    produce <- function(ids, B, V, pack) {
      units <- .inlast_sparse_units(fit, ids, threads = threads)
      unit_error <- vapply(units, function(z) {
        if (is.null(z$error)) NA_character_ else as.character(z$error)
      }, character(1L))
      good <- which(is.na(unit_error))
      out <- list(error = unit_error)
      if (length(good)) {
        z <- mgcvst_inla_sparse_materialize_pca_cpp(
          units[good], g$cache$general_Q, as.numeric(g$constraint),
          basis$coordinate, basis$basis, B, rep(pack, length(good)),
          threads, g$cache$prepared, V
        )
        out$good <- good
        out$z <- z
        out$error[good] <- z$error
      }
      out
    }
  } else {
    builder <- .mgcvst_state_builder(fit, used, threads)
    q <- .mgcvst_state_width(fit)
    produce <- function(ids, B, V, pack) {
      states <- builder$build(ids)
      state_error <- vapply(states, function(z) {
        if (is.null(z$error)) NA_character_ else as.character(z$error)
      }, character(1L))
      good <- which(is.na(state_error))
      out <- list(error = state_error)
      if (length(good)) {
        z <- mgcvst_pca_dense_cpp(
          lapply(states[good], `[[`, "M"),
          do.call(cbind, lapply(states[good], `[[`, "a")), B,
          rep(pack, length(good)), V, threads
        )
        out$good <- good
        out$z <- z
        out$error[good] <- z$error
      }
      out
    }
  }
  materialize <- function(ids, B = NULL, V = NULL, pack = FALSE) {
    n <- length(ids)
    r <- if (is.null(B)) 0L else ncol(B)
    k <- if (is.null(V)) 0L else ncol(V)
    A <- matrix(NA_real_, q, n)
    C <- matrix(NA_real_, n, r)
    R <- matrix(NA_real_, k * k, n)
    fro2 <- scale <- rep(NA_real_, n)
    error <- rep(NA_character_, n)
    packed <- vector("list", n)
    for (first in if (n) seq.int(1L, n, by = 32L) else integer()) {
      kk <- first:min(n, first + 31L)
      out <- produce(ids[kk], B, V, pack)
      error[kk] <- out$error
      if (is.null(out$good)) next
      kg <- kk[out$good]
      z <- out$z
      A[, kg] <- z$a
      if (r) C[kg, ] <- z$C
      fro2[kg] <- z$fro2
      scale[kg] <- z$scale
      if (k) R[, kg] <- z$R
      if (pack) packed[kg] <- z$packed
    }
    list(A = A, C = C, fro2 = fro2, error = error, packed = packed, R = R,
         scale = scale)
  }
  list(width = q, materialize = materialize)
}

## ---- BASIS ----
# Orthonormal basis B (weighted vech, q (q + 1) / 2 x rank) of
# V = [tau_j vech_w(H_j)] from the Gram eigen decomposition,
# B = V R with R = U_r diag(1 / sqrt(lambda_r)); the shared basis V (q x k) from
# the leading eigenvectors of sum_j H_j / scale_j over the training genes; and
# the factors R_j of the training genes, taken from their packed matrices so
# that the training genes are not materialized twice. The achievable rank is the
# number of materialized training genes with Gram eigenvalue above 1e-10 of the
# largest.
.mgcvst_pca_basis <- function(producer, train, tau, rank, k, threads) {
  t0 <- proc.time()[["elapsed"]]
  z <- producer$materialize(train, pack = TRUE)
  ok <- is.na(z$error)
  q <- producer$width
  t1 <- proc.time()[["elapsed"]]
  gram <- mgcvst_pca_gram_cpp(z$packed[ok], tau[ok], threads)
  t2 <- proc.time()[["elapsed"]]
  eg <- eigen(gram, symmetric = TRUE)
  positive <- if (sum(ok)) sum(eg$values > 1e-10 * eg$values[1L]) else 0L
  if (rank > positive) {
    stop("rank = ", rank, " exceeds the achievable PCAlearning rank ", positive,
         " (", sum(ok), " training genes materialized, ", positive,
         " Gram eigenvalues above 1e-10 of the largest); use a smaller rank, or ",
         "moments = \"exact\" for so few genes.")
  }
  rotation <- sweep(eg$vectors[, seq_len(rank), drop = FALSE], 2L,
                    sqrt(eg$values[seq_len(rank)]), "/")
  t3 <- proc.time()[["elapsed"]]
  B <- mgcvst_pca_basis_cpp(z$packed[ok], tau[ok], rotation, threads)
  t4 <- proc.time()[["elapsed"]]
  S <- mgcvst_pca_packed_sum_cpp(z$packed[ok], q)
  kv <- as.integer(min(k, q))
  V <- eigen(S, symmetric = TRUE)$vectors[, seq_len(kv), drop = FALSE]
  proj <- mgcvst_pca_packed_project_cpp(z$packed[ok], V, threads)
  t5 <- proc.time()[["elapsed"]]
  list(B = B, V = V, V_sha = digest::digest(V, algo = "sha256"), k = kv,
       train = train[ok], tau = tau[ok], values = eg$values,
       rotation = rotation, A = z$A[, ok, drop = FALSE],
       C = (gram %*% rotation) / tau[ok], fro2 = z$fro2[ok],
       R = proj$R, scale = proj$scale, error = proj$error,
       elapsed = c(materialize = t1 - t0, gram = t2 - t1, eigen = t3 - t2,
                   basis = t4 - t3, shared = t5 - t4))
}

# Peak bytes of mgcvst_pca_tables_cpp: the basis matrices (double), the
# r (r + 1) / 2 products B_a B_b with a <= b (float), the float panel of the
# Gram step (block 32) and the double Gram of the r^2 products.
.mgcvst_pca_table_bytes <- function(q, r) {
  8 * q^2 * r + 4 * q^2 * (r * (r + 1) / 2) + 4 * q * 32 * r^2 + 8 * r^4
}

# Atomic RDS commit in a checkpoint directory.
.mgcvst_pca_save <- function(object, path, file) {
  target <- file.path(path, file)
  tmp <- tempfile("pca-", tmpdir = path, fileext = ".tmp")
  on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
  saveRDS(object, tmp, compress = FALSE)
  if (file.exists(target) || !file.rename(tmp, target)) {
    stop("Could not commit the PCAlearning checkpoint file ", file, ".")
  }
  invisible(NULL)
}

# Completed projection blocks: feature_index, A, C, fro2, error, R and scale per
# gene.
.mgcvst_pca_blocks_read <- function(path) {
  files <- sort(list.files(path, "^pca-projection-[0-9]{6}\\.rds$"))
  lapply(files, function(file) {
    z <- readRDS(file.path(path, file))
    body <- z[c("feature_index", "A", "C", "fro2", "error", "R", "scale")]
    if (!identical(z$checksum, digest::digest(body, algo = "sha256"))) {
      stop("The PCAlearning checkpoint block is damaged: ", file, ".")
    }
    body
  })
}

# Defaults of the PCAlearning controls of the pair tests, in one place: the
# rank r of the trace tables and the number k of leading singular values taken
# from the shared basis. On the full-rank MAGIC basis (q = 1961) rank 20 with
# k = 50 missed the accuracy limits at p = 1e-12 and 1e-20, and rank 30 with
# k = 80 met them.
.mgcvst_pca_defaults <- list(rank = 30L, n_per_cell = 3L, seed = 1L, k = 80L)

# Validated PCAlearning controls as integers. The pair tests call this before
# they build any basis, so a bad control fails at once.
.mgcvst_pca_check_args <- function(rank, n_per_cell, seed, k = NULL) {
  for (x in list(list("rank", rank), list("n_per_cell", n_per_cell))) {
    v <- x[[2L]]
    if (!is.numeric(v) || length(v) != 1L || !is.finite(v) || v < 1 ||
        v != floor(v) || v > .Machine$integer.max) {
      stop(x[[1L]], " must be one positive integer.", call. = FALSE)
    }
  }
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
      seed < 0 || seed != floor(seed) || seed > .Machine$integer.max) {
    stop("seed must be one non-negative integer.", call. = FALSE)
  }
  if (!is.null(k) && (!is.numeric(k) || length(k) != 1L || !is.finite(k) ||
                      k < 1 || k != floor(k) || k > .Machine$integer.max)) {
    stop("k must be NULL or one positive integer.", call. = FALSE)
  }
  list(rank = as.integer(rank), n_per_cell = as.integer(n_per_cell),
       seed = as.integer(seed), k = if (is.null(k)) NULL else as.integer(k))
}

# Everything the pair kernels need of the genes `used`: the training genes and
# the PCAlearning basis B (learned, or read from the checkpoint directory
# `path`), the shared basis V, the coefficients C, the factors R_g of V, the
# scales, the score coordinates A, the degree-2 monomials K2 and the trace
# table T2. Genes without a usable state (`failed_gene`) and genes that the
# Stage 1 test did not select (`skipped`, neither trained on nor materialized)
# carry neutral values. The pair stage and the validation scripts both start from this
# object.
.mgcvst_pca_prepare <- function(fit, used, basis, q, rank, n_per_cell, seed, k,
                                threads, path = NULL, verbose = FALSE,
                                started = proc.time()[["elapsed"]]) {
  k_requested <- k
  basis_file <- if (is.null(path)) NULL else file.path(path, "pca-basis.rds")
  degenerate <- .mgcvst_stage1_unselected(fit)
  skipped <- degenerate[used]
  available <- .mgcvst_feature_available(fit) & !degenerate
  basis_resumed <- !is.null(basis_file) && file.exists(basis_file)
  if (basis_resumed) {
    saved <- readRDS(basis_file)
    if (!identical(saved$checksum, digest::digest(saved[c("sampled", "learned")],
                                                  algo = "sha256"))) {
      stop("The PCAlearning checkpoint basis is damaged: ", basis_file, ".")
    }
    sampled <- saved$sampled
    learned <- saved$learned
    t_sample <- 0
    if (verbose) message("Resumed the rank-", rank, " PCAlearning basis.")
  } else {
    scales <- .mgcvst_pca_scales(fit)
    sampled <- .mgcvst_pca_training(scales, which(available), n_per_cell, seed)
    sampled$scales <- scales
    t_sample <- proc.time()[["elapsed"]] - started
    if (verbose) message("Sampled ", length(sampled$train),
                         " PCAlearning training genes.")
  }
  # The producer serves the training genes as well as the genes of the pairs.
  producer <- .mgcvst_pca_producer(
    fit, sort(unique(c(used[!skipped], sampled$train))), basis, threads)
  if (!basis_resumed) {
    learned <- .mgcvst_pca_basis(producer, sampled$train,
                                 sampled$scales$tau[sampled$train], rank,
                                 k_requested, threads)
    if (verbose) message("Learned the rank-", rank, " PCAlearning basis and the ",
                         "shared basis of k = ", learned$k, ".")
    if (!is.null(path)) {
      body <- list(sampled = sampled, learned = learned)
      body$checksum <- digest::digest(body, algo = "sha256")
      .mgcvst_pca_save(body, path, "pca-basis.rds")
    }
  }
  scales <- sampled$scales
  kv <- learned$k

  # Training genes reuse their first materialization; the others are projected
  # in checkpoint blocks of 256 genes.
  t0 <- proc.time()[["elapsed"]]
  rest <- setdiff(used[!skipped], learned$train)
  blocks <- if (is.null(path)) list() else .mgcvst_pca_blocks_read(path)
  done <- unlist(lapply(blocks, `[[`, "feature_index"))
  resumed_genes <- sum(rest %in% done)
  missing <- setdiff(rest, done)
  next_block <- if (is.null(path)) 0L else max(0L, as.integer(substr(
    list.files(path, "^pca-projection-[0-9]{6}\\.rds$"), 16L, 21L)))
  for (first in if (length(missing)) seq.int(1L, length(missing), by = 256L) else integer()) {
    ids <- missing[first:min(length(missing), first + 255L)]
    z <- producer$materialize(ids, B = learned$B, V = learned$V)
    block <- list(feature_index = ids, A = z$A, C = z$C, fro2 = z$fro2,
                  error = z$error, R = z$R, scale = z$scale)
    if (!is.null(path)) {
      next_block <- next_block + 1L
      .mgcvst_pca_save(c(block, list(checksum = digest::digest(block, algo = "sha256"))),
                       path, sprintf("pca-projection-%06d.rds", next_block))
    }
    blocks[[length(blocks) + 1L]] <- block
    if (verbose) message("Projected ", min(length(missing), first + 255L), " of ",
                         length(missing), " PCAlearning genes.")
  }
  proj <- list(
    ids = unlist(lapply(blocks, `[[`, "feature_index")),
    A = do.call(cbind, c(list(matrix(0, q, 0L)), lapply(blocks, `[[`, "A"))),
    C = do.call(rbind, c(list(matrix(0, 0L, rank)), lapply(blocks, `[[`, "C"))),
    fro2 = unlist(lapply(blocks, `[[`, "fro2")),
    error = unlist(lapply(blocks, `[[`, "error")),
    R = do.call(cbind, c(list(matrix(0, kv * kv, 0L)), lapply(blocks, `[[`, "R"))),
    scale = unlist(lapply(blocks, `[[`, "scale"))
  )
  rm(blocks)
  kt <- match(used, learned$train)
  kr <- match(used, proj$ids)
  tr <- !is.na(kt)
  other <- which(!tr & !skipped)
  A <- matrix(NA_real_, q, length(used))
  C <- matrix(NA_real_, length(used), ncol(learned$B))
  R <- matrix(NA_real_, kv * kv, length(used))
  scale <- rep(NA_real_, length(used))
  fro2 <- rep(NA_real_, length(used))
  gene_error <- rep(NA_character_, length(used))
  A[, tr] <- learned$A[, kt[tr]]
  A[, other] <- proj$A[, kr[other]]
  C[tr, ] <- learned$C[kt[tr], ]
  C[other, ] <- proj$C[kr[other], ]
  R[, tr] <- learned$R[, kt[tr]]
  R[, other] <- proj$R[, kr[other]]
  scale[tr] <- learned$scale[kt[tr]]
  scale[other] <- proj$scale[kr[other]]
  fro2[tr] <- learned$fro2[kt[tr]]
  fro2[other] <- proj$fro2[kr[other]]
  gene_error[tr] <- learned$error[kt[tr]]
  gene_error[other] <- proj$error[kr[other]]
  learned$A <- NULL
  learned$R <- NULL
  rm(proj)
  t_project <- proc.time()[["elapsed"]] - t0

  # A gene without a usable state, or one that the Stage 1 test did not select,
  # enters the kernels with neutral values; its pairs are marked below (status 3 and
  # status 4) and carry no kernel result.
  failed_gene <- !is.na(gene_error)
  C_out <- C
  dead <- failed_gene | skipped
  if (any(dead)) {
    A[, dead] <- 0
    C[dead, ] <- 0
    R[, dead] <- 0
    scale[dead] <- 1
  }
  # The basis is orthonormal: t_1 = c_i' c_j needs no level-1 table.
  basis_check <- max(abs(crossprod(learned$B) - diag(ncol(learned$B))))
  tables <- mgcvst_pca_tables_cpp(learned$B, q, threads)
  K2 <- mgcvst_pca_monomials_cpp(C)
  T2 <- tables$Tsym2
  learned$B <- NULL
  if (verbose) message("Computed PCAlearning trace tables in ",
                       round(tables$timing[["total"]], 2), " s.")
  list(A = A, C = C, C_out = C_out, K2 = K2, T2 = T2, R = R, scale = scale,
       used = used, k = kv, V = learned$V, V_sha = learned$V_sha, fro2 = fro2,
       gene_error = gene_error, failed_gene = failed_gene, skipped = skipped,
       trained = tr,
       sampled = sampled, scales = scales, learned = learned, tables = tables,
       basis_check = basis_check, basis_resumed = basis_resumed,
       resumed_genes = resumed_genes, projected_genes = length(missing),
       t_sample = t_sample, t_project = t_project)
}

# Saddlepoint log p-values for the requested pairs from the PCAlearning basis,
# streamed to raw Parquet shards. `index` is NULL for every pair of the
# available genes (gene blocks are generated and scored on the fly) or a
# two-column matrix of available feature indices i < j, evaluated by (i, j).
# A sparse INLA fit needs its observation `basis`; an mgcv fit has none.
.mgcvst_pair_pcalearning <- function(fit, index, threads, chunk_size, verbose,
                                     basis = NULL,
                                     rank = .mgcvst_pca_defaults$rank,
                                     n_per_cell = .mgcvst_pca_defaults$n_per_cell,
                                     seed = .mgcvst_pca_defaults$seed,
                                     k = .mgcvst_pca_defaults$k,
                                     checkpoint_dir = NULL, resume = TRUE,
                                     route = NULL) {
  started <- proc.time()[["elapsed"]]
  controls <- .mgcvst_pca_check_args(rank, n_per_cell, seed, k)
  rank <- controls$rank
  n_per_cell <- controls$n_per_cell
  seed <- controls$seed
  k_requested <- if (is.null(controls$k)) .mgcvst_pca_defaults$k else controls$k
  chunk_size <- .mgcvst_check_chunk_size(chunk_size)
  if (is.null(chunk_size)) stop("chunk_size must be one positive integer.")
  threads <- as.integer(threads)
  sparse <- identical(fit$score_backend, "sparse")
  if (sparse && is.null(basis)) stop("A sparse INLA fit needs its observation basis.")
  q <- .mgcvst_state_width(fit, basis)
  table_bytes <- .mgcvst_pca_table_bytes(q, rank)
  available_memory <- .mgcvst_memory_probe()$available
  if (is.finite(available_memory) && table_bytes > available_memory) {
    stop(sprintf(paste0("PCAlearning trace tables for rank = %d and q = %d need ",
                        "about %.1f GB, above the %.1f GB of available memory; ",
                        "use a smaller rank."),
                 rank, q, table_bytes / 1024^3, available_memory / 1024^3))
  }
  if (sparse) fit <- .inlast_sparse_prepare(fit)

  available <- .mgcvst_feature_available(fit)
  all_pairs <- is.null(index)
  if (all_pairs) {
    used <- which(available)
    if (length(used) < 2L) stop("At least two available features are required.")
  } else {
    if (!is.matrix(index) || ncol(index) != 2L || !nrow(index) ||
        anyNA(index) || any(index[, 1L] >= index[, 2L]) ||
        !all(available[index])) {
      stop("index must contain two-column indices i < j of available features.")
    }
    used <- sort(unique(as.vector(index)))
  }

  contract_early <- .mgcvst_contract("pcalearning")
  path <- NULL
  if (!is.null(checkpoint_dir)) {
    signature <- list(version = 3L, method = "spa_pcalearning",
                      contract = contract_early,
                      fit = .mgcvst_pair_signature(fit, basis), rank = rank,
                      n_per_cell = n_per_cell, seed = seed, k = k_requested,
                      levels = 2L, q = q)
    manifest <- file.path(checkpoint_dir, "manifest.rds")
    if (resume && file.exists(manifest)) {
      old <- readRDS(manifest)$signature
      if (is.null(old$method)) {
        stop("The checkpoint in ", checkpoint_dir, " holds the score states of ",
             "the exact route; use a new checkpoint_dir for the PCAlearning ",
             "route.", call. = FALSE)
      }
      if (!identical(old$contract$calibration_contract,
                     contract_early$calibration_contract)) {
        stop("The checkpoint in ", checkpoint_dir, " was written under an ",
             "earlier algorithm contract (this version writes ",
             contract_early$calibration_contract, "); use a new checkpoint_dir.")
      }
      if (!identical(old, signature)) {
        stop("The checkpoint in ", checkpoint_dir, " was written for a different ",
             "fit, score basis, route, rank, k, n_per_cell, or seed; use a new ",
             "checkpoint_dir.")
      }
    }
    path <- .mgcvst_store_open(checkpoint_dir, signature, fit$feature_id,
                               resume = resume)$path
    .mgcvst_route_save(path, route)
  }
  root <- if (is.null(path)) tempfile("mgcvst-pairs-") else path
  # Pair directories of another algorithm contract are refused before any
  # work; the directory of this run is opened once the shared basis exists.
  .mgcvst_pairs_refuse_stale(root, contract_early)
  degenerate <- .mgcvst_stage1_unselected(fit)
  universe <- if (all_pairs) {
    list(all = TRUE, used = used, n_feature = length(fit$feature_id),
         degenerate = which(degenerate[used]))
  } else list(index = index, degenerate = which(degenerate[used]))

  prep <- .mgcvst_pca_prepare(fit, used, basis, q, rank, n_per_cell, seed,
                              k_requested, threads, path, verbose, started)
  A <- prep$A; C <- prep$C; C_out <- prep$C_out; K2 <- prep$K2; T2 <- prep$T2
  R <- prep$R; scale <- prep$scale; kv <- prep$k; fro2 <- prep$fro2
  gene_error <- prep$gene_error; tr <- prep$trained; sampled <- prep$sampled
  skipped <- prep$skipped
  scales <- prep$scales; learned <- prep$learned; tables <- prep$tables
  t_sample <- prep$t_sample; t_project <- prep$t_project
  preparation_elapsed <- proc.time()[["elapsed"]] - started

  # Compact rows of one block of pairs: local positions (li, lj) in `used`
  # and the kernel output matrix. A pair with a gene whose projection failed
  # has status 3 and no p-value. Otherwise a pair with an unselected gene has
  # status 4, p = 1 and no score. A pair the kernel could not evaluate keeps its
  # status 1 or 2 and has p = 1: two-sided and both one-sided log p are 0.
  ok <- is.na(gene_error)
  pcols <- c("log_p_two_sided", "log_p_positive", "log_p_negative")
  pair_frame <- function(li, lj, out) {
    frame <- .mgcvst_pairs_frame(
      used[li], used[lj], out[, "U"], out[, "logp_two_sided"],
      out[, "logp_positive"], out[, "logp_negative"],
      remainder_kind = out[, "remainder_kind"], status = out[, "status"]
    )
    invalid <- frame$status != .mgcvst_pair_status[["ok"]]
    frame[invalid, pcols] <- 0
    bad <- !ok[li] | !ok[lj]
    skip <- !bad & (skipped[li] | skipped[lj])
    frame$status[skip] <- .mgcvst_pair_status[["degenerate"]]
    frame$remainder_kind[skip] <- 0L
    frame$score[skip] <- NA_real_
    frame[skip, pcols] <- 0
    frame$status[bad] <- .mgcvst_pair_status[["feature"]]
    frame$remainder_kind[bad] <- 0L
    frame[bad, c("score", pcols)] <- NA_real_
    frame
  }

  contract <- .mgcvst_contract("pcalearning", k = kv, basis_sha = learned$V_sha)
  pair_dir <- .mgcvst_pairs_open(root, universe, contract, resume)
  # The block schedule depends on chunk_size, so it is fixed when the pair
  # directory is first used and a resumed run follows the stored one, whatever
  # chunk_size it is given.
  n_schedule <- if (all_pairs) length(used) else nrow(index)
  schedule_file <- file.path(pair_dir, "schedule.rds")
  if (file.exists(schedule_file)) {
    stored <- .mgcvst_schedule_load(schedule_file)
    if (!identical(stored$kind, "pcalearning") || !identical(stored$n, n_schedule) ||
        !identical(stored$all_pairs, all_pairs)) {
      stop("The stored pair schedule is damaged or incompatible.")
    }
    if (verbose && !identical(stored$chunk_size, as.integer(chunk_size))) {
      message("Resumed with the stored chunk_size = ", stored$chunk_size, ".")
    }
    chunk_size <- stored$chunk_size
  } else {
    chunk_size <- as.integer(min(chunk_size, .Machine$integer.max))
    .mgcvst_schedule_save(schedule_file, list(
      kind = "pcalearning", all_pairs = all_pairs, n = n_schedule,
      chunk_size = chunk_size))
  }

  t0 <- proc.time()[["elapsed"]]
  chunks <- 0L
  resumed_pairs <- 0
  nodes_above <- 0
  shard_files <- character()
  shard_rows <- integer()
  n_used <- length(used)
  if (all_pairs) {
    per_gene <- n_used - seq_len(n_used)
    first <- 1L
    id <- 0L
    while (first < n_used) {
      last <- first
      total <- per_gene[first]
      while (last + 1L < n_used && total + per_gene[last + 1L] <= chunk_size) {
        last <- last + 1L
        total <- total + per_gene[last]
      }
      id <- id + 1L
      left <- first:last
      len <- n_used - left
      li <- rep(left, len)
      lj <- sequence(len, from = left + 1L)
      if (.mgcvst_shard_complete(pair_dir, id, used[li], used[lj])) {
        resumed_pairs <- resumed_pairs + length(li)
      } else {
        out <- mgcvst_pca_spa_block_cpp(A, C, K2, T2, R, scale, first, last, threads)
        nodes_above <- nodes_above + attr(out, "nodes_above_leading")
        .mgcvst_write_parquet(pair_frame(li, lj, out),
                              .mgcvst_shard_file(pair_dir, id))
      }
      shard_files <- c(shard_files, .mgcvst_shard_file(pair_dir, id))
      shard_rows <- c(shard_rows, length(li))
      chunks <- chunks + 1L
      first <- last + 1L
    }
  } else {
    local <- matrix(match(index, used), ncol = 2L)
    n <- nrow(index)
    for (first in seq.int(1L, n, by = chunk_size)) {
      z <- first:min(n, first + chunk_size - 1L)
      if (.mgcvst_shard_complete(pair_dir, first, index[z, 1L], index[z, 2L])) {
        resumed_pairs <- resumed_pairs + length(z)
      } else {
        out <- mgcvst_pca_spa_pairs_cpp(A, C, K2, T2, R, scale, local[z, 1L],
                                        local[z, 2L], threads)
        nodes_above <- nodes_above + attr(out, "nodes_above_leading")
        .mgcvst_write_parquet(pair_frame(local[z, 1L], local[z, 2L], out),
                              .mgcvst_shard_file(pair_dir, first))
      }
      shard_files <- c(shard_files, .mgcvst_shard_file(pair_dir, first))
      shard_rows <- c(shard_rows, length(z))
      chunks <- chunks + 1L
    }
  }
  pair_elapsed <- proc.time()[["elapsed"]] - t0

  C <- C_out
  e2 <- fro2 - rowSums(C^2)
  genes <- data.frame(
    feature_id = fit$feature_id[used], feature_index = used,
    training = tr, cell = sampled$cell[used],
    sigma_g2 = scales$sigma_g2[used], sigma_e2 = scales$sigma_e2[used],
    tau = scales$tau[used], fro2 = fro2, e2 = e2,
    e2_relative = e2 / fro2, error_message = gene_error, degenerate = skipped,
    stringsAsFactors = FALSE
  )
  dimnames(C) <- list(fit$feature_id[used], paste0("c", seq_len(ncol(C))))
  training <- data.frame(
    feature_id = fit$feature_id[learned$train], feature_index = learned$train,
    cell = sampled$cell[learned$train], tau = learned$tau,
    stringsAsFactors = FALSE
  )
  failed <- genes[!is.na(genes$error_message), c("feature_id", "error_message")]
  names(failed) <- c("feature_id", "error")
  rownames(failed) <- NULL
  list(
    pair_dir = pair_dir, shards = shard_files, rows = shard_rows,
    n_pairs = sum(shard_rows), temporary = is.null(path), failed = failed,
    elapsed = pair_elapsed,
    metadata = list(
      preparation_backend = if (sparse) "sparse" else "model_native",
      pair_schedule = if (all_pairs) "pcalearning_gene_blocks" else
        "pcalearning_pair_list",
      chunks = chunks, chunk_size = chunk_size, resumed_pairs = resumed_pairs,
      pair_dir = pair_dir, k = kv, basis_sha = learned$V_sha, q = q,
      nodes_above_leading = nodes_above,
      degenerate_genes = sum(skipped),
      preparation_elapsed = preparation_elapsed, contract = contract,
      pca_learning = list(
        rank = rank, n_per_cell = n_per_cell,
        seed = seed, k = kv, q = q, training = training,
        checkpoint = list(path = path, basis_resumed = prep$basis_resumed,
                          resumed_genes = prep$resumed_genes,
                          projected_genes = prep$projected_genes),
        basis_orthonormality = prep$basis_check,
        gram_values = learned$values, rotation = learned$rotation,
        coefficients = C, genes = genes,
        table_timing = tables$timing, table_memory = tables$memory,
        elapsed = c(sample = t_sample, learned$elapsed, project = t_project,
                    tables = tables$timing[["total"]], pairs = pair_elapsed)
      )
    )
  )
}
