# PCAlearning approximate Liu path for sparse INLA pair tests.
#
# Each score covariance H_j is projected onto an r-dimensional orthonormal
# basis B of symmetric q x q matrices learned from training genes. With
# c_j = B' vech_w(H_j) (weighted vech, Frobenius isometric), the pair traces
# tr((H_i H_j)^s), s = 1, ..., 4, are
# replaced by the traces of the projected matrices, which the trace tables of
# B give by contraction of monomials of c_i and c_j.

# Per-gene variance scales of the first-stage fit that define training strata.
# NB: V = 1 / mu + 1 / theta; Poisson (and other families): V = 1 / mu.
.mgcvst_pca_scales <- function(fit, threads = 1L) {
  G <- length(fit$feature_id)
  nb <- fit$diagnostics$family_used == "negative_binomial"
  theta <- vapply(fit$family_parameters, function(x) {
    if (length(x)) x[1L] else NA_real_
  }, numeric(1L))
  mu_bar <- numeric(G)
  for (first in seq.int(1L, G, by = 500L)) {
    ids <- first:min(first + 499L, G)
    it <- ifelse(nb[ids], 1 / theta[ids], 0)
    V <- .inlast_working_state(fit, ids, threads = threads)$variance
    D <- sweep(V, 2L, it, "-")
    mu_bar[ids] <- colMeans(1 / D)
  }
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


# Materialize genes `ids` in blocks of 32: score coordinates a_j (q x n),
# coefficients c_j = B' vech_w(H_j) when `B` is given, ||H_j||_F^2 and, with
# `pack = TRUE`, float32 weighted-vech H_j. Failed genes carry `error`.
.mgcvst_pca_materialize <- function(fit, ids, basis, threads, B = NULL,
                                    pack = FALSE) {
  g <- fit$score_sparse
  n <- length(ids)
  A <- matrix(NA_real_, basis$rank, n)
  C <- matrix(NA_real_, n, if (is.null(B)) 0L else ncol(B))
  fro2 <- rep(NA_real_, n)
  error <- rep(NA_character_, n)
  packed <- vector("list", n)
  for (first in if (n) seq.int(1L, n, by = 32L) else integer()) {
    k <- first:min(n, first + 31L)
    units <- .inlast_sparse_units(fit, ids[k], threads = threads)
    unit_error <- vapply(units, function(z) {
      if (is.null(z$error)) NA_character_ else as.character(z$error)
    }, character(1L))
    error[k] <- unit_error
    good <- which(is.na(unit_error))
    if (!length(good)) next
    z <- mgcvst_inla_sparse_materialize_pca_cpp(
      units[good], g$cache$general_Q, as.numeric(g$constraint),
      basis$coordinate, basis$basis, B, rep(pack, length(good)),
      threads, g$cache$prepared
    )
    kg <- k[good]
    A[, kg] <- z$a
    if (!is.null(B)) C[kg, ] <- z$C
    fro2[kg] <- z$fro2
    error[kg] <- z$error
    if (pack) packed[kg] <- z$packed
  }
  list(A = A, C = C, fro2 = fro2, error = error, packed = packed)
}

# Orthonormal basis B (weighted vech, q (q + 1) / 2 x rank) of
# V = [tau_j vech_w(H_j)] from the Gram eigen decomposition,
# B = V R with R = U_r diag(1 / sqrt(lambda_r)). Training genes are
# materialized once; their coefficients follow from V'B = G R, so that
# c_j = (G R)_j / tau_j. The achievable rank is the number of materialized
# training genes with Gram eigenvalue above 1e-10 of the largest.
.mgcvst_pca_basis <- function(fit, train, tau, basis, rank, threads) {
  t0 <- proc.time()[["elapsed"]]
  z <- .mgcvst_pca_materialize(fit, train, basis, threads, pack = TRUE)
  ok <- is.na(z$error)
  t1 <- proc.time()[["elapsed"]]
  gram <- mgcvst_pca_gram_cpp(z$packed[ok], tau[ok], threads)
  t2 <- proc.time()[["elapsed"]]
  eg <- eigen(gram, symmetric = TRUE)
  positive <- if (sum(ok)) sum(eg$values > 1e-10 * eg$values[1L]) else 0L
  if (rank > positive) {
    stop("rank = ", rank, " exceeds the achievable PCAlearning rank ", positive,
         " (", sum(ok), " training genes materialized, ", positive,
         " Gram eigenvalues above 1e-10 of the largest).")
  }
  rotation <- sweep(eg$vectors[, seq_len(rank), drop = FALSE], 2L,
                    sqrt(eg$values[seq_len(rank)]), "/")
  t3 <- proc.time()[["elapsed"]]
  B <- mgcvst_pca_basis_cpp(z$packed[ok], tau[ok], rotation, threads)
  t4 <- proc.time()[["elapsed"]]
  list(B = B, train = train[ok], tau = tau[ok], values = eg$values,
       rotation = rotation, A = z$A[, ok, drop = FALSE],
       C = (gram %*% rotation) / tau[ok], fro2 = z$fro2[ok],
       elapsed = c(materialize = t1 - t0, gram = t2 - t1, eigen = t3 - t2,
                   basis = t4 - t3))
}

# Peak bytes of mgcvst_pca_tables_cpp: the larger of the level-3 stage
# (float B, Y = B_a B_b, N) and the level-4 stage (float B, N, Y sums,
# W panel, double Grams G3 and G4); d2 = r (r + 1) / 2.
.mgcvst_pca_table_bytes <- function(q, r) {
  d2 <- r * (r + 1) / 2
  level3 <- 4 * q^2 * (r + r^2 + d2 * r)
  level4 <- 4 * q^2 * (r + d2 * r + d2) + 4 * 32 * q * d2^2 +
    8 * d2^4 + 8 * (d2 * r)^2
  max(level3, level4)
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

# Completed projection blocks: feature_index, A, C, fro2, error per gene.
.mgcvst_pca_blocks_read <- function(path) {
  files <- sort(list.files(path, "^pca-projection-[0-9]{6}\\.rds$"))
  lapply(files, function(file) {
    z <- readRDS(file.path(path, file))
    body <- z[c("feature_index", "A", "C", "fro2", "error")]
    if (!identical(z$checksum, digest::digest(body, algo = "sha256"))) {
      stop("The PCAlearning checkpoint block is damaged: ", file, ".")
    }
    body
  })
}

# Defaults of the PCAlearning controls of inlaST.test(), in one place.
.mgcvst_pca_defaults <- list(rank = 10L, n_per_cell = 3L, seed = 1L)

# Approximate Liu log p-values for the requested pairs from the PCAlearning
# basis, streamed to raw Parquet shards. `index` is NULL for every pair of the
# available genes (gene blocks are generated and scored on the fly) or a
# two-column matrix of available feature indices i < j, evaluated by (i, j).
.mgcvst_pair_pcalearning <- function(fit, index, threads, chunk_size, verbose,
                                     basis, rank = .mgcvst_pca_defaults$rank,
                                     n_per_cell = .mgcvst_pca_defaults$n_per_cell,
                                     seed = .mgcvst_pca_defaults$seed,
                                     checkpoint_dir = NULL, resume = TRUE) {
  started <- proc.time()[["elapsed"]]
  for (x in list(list("rank", rank), list("n_per_cell", n_per_cell),
                 list("chunk_size", chunk_size))) {
    v <- x[[2L]]
    if (!is.numeric(v) || length(v) != 1L || !is.finite(v) || v < 1 ||
        v != floor(v)) stop(x[[1L]], " must be one positive integer.")
  }
  rank <- as.integer(rank)
  n_per_cell <- as.integer(n_per_cell)
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
      seed < 0 || seed != floor(seed) || seed > .Machine$integer.max) {
    stop("seed must be one non-negative integer.")
  }
  seed <- as.integer(seed)
  threads <- as.integer(threads)
  table_bytes <- .mgcvst_pca_table_bytes(basis$rank, rank)
  available_memory <- .mgcvst_memory_probe()$available
  if (is.finite(available_memory) && table_bytes > available_memory) {
    stop(sprintf(paste0("PCAlearning trace tables for rank = %d and q = %d need ",
                        "about %.1f GB, above the %.1f GB of available memory; ",
                        "use a smaller rank."),
                 rank, basis$rank, table_bytes / 1024^3,
                 available_memory / 1024^3))
  }
  fit <- .inlast_sparse_prepare(fit)

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

  contract <- .mgcvst_contract("pcalearning")
  path <- NULL
  if (!is.null(checkpoint_dir)) {
    signature <- list(version = 2L, method = "pca_learning",
                      contract = contract,
                      fit = .mgcvst_pair_signature(fit, basis), rank = rank,
                      n_per_cell = n_per_cell, seed = seed,
                      q = basis$rank)
    manifest <- file.path(checkpoint_dir, "manifest.rds")
    if (resume && file.exists(manifest)) {
      old <- readRDS(manifest)$signature
      if (!identical(old$version, 2L) || !identical(old$method, "pca_learning") ||
          !identical(old$contract$calibration_contract,
                     contract$calibration_contract)) {
        stop("The checkpoint in ", checkpoint_dir, " was written under an ",
             "earlier PCAlearning algorithm contract (this version writes ",
             contract$calibration_contract, "); use a new checkpoint_dir.")
      }
      if (!identical(old, signature)) {
        stop("The checkpoint in ", checkpoint_dir, " was written for a different ",
             "fit, score basis, rank, n_per_cell, or seed; use a new checkpoint_dir.")
      }
    }
    path <- .mgcvst_store_open(checkpoint_dir, signature, fit$feature_id,
                               resume = resume)$path
  }
  root <- if (is.null(path)) tempfile("mgcvst-pairs-") else path
  universe <- if (all_pairs) {
    list(all = TRUE, used = used, n_feature = length(fit$feature_id))
  } else list(index = index)
  pair_dir <- .mgcvst_pairs_open(root, universe, contract, resume)

  basis_file <- if (is.null(path)) NULL else file.path(path, "pca-basis.rds")
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
    universe_genes <- which(available)
    scales <- .mgcvst_pca_scales(fit, threads = threads)
    sampled <- .mgcvst_pca_training(scales, universe_genes, n_per_cell, seed)
    sampled$scales <- scales
    t_sample <- proc.time()[["elapsed"]] - started
    if (verbose) message("Sampled ", length(sampled$train),
                         " PCAlearning training genes.")
    learned <- .mgcvst_pca_basis(fit, sampled$train,
                                 scales$tau[sampled$train], basis, rank, threads)
    if (verbose) message("Learned the rank-", rank, " PCAlearning basis.")
    if (!is.null(path)) {
      body <- list(sampled = sampled, learned = learned)
      body$checksum <- digest::digest(body, algo = "sha256")
      .mgcvst_pca_save(body, path, "pca-basis.rds")
    }
  }
  scales <- sampled$scales

  # Training genes reuse their first materialization; the others are projected
  # in checkpoint blocks of 256 genes.
  t0 <- proc.time()[["elapsed"]]
  rest <- setdiff(used, learned$train)
  blocks <- if (is.null(path)) list() else .mgcvst_pca_blocks_read(path)
  done <- unlist(lapply(blocks, `[[`, "feature_index"))
  resumed_genes <- sum(rest %in% done)
  missing <- setdiff(rest, done)
  next_block <- if (is.null(path)) 0L else max(0L, as.integer(substr(
    list.files(path, "^pca-projection-[0-9]{6}\\.rds$"), 16L, 21L)))
  for (first in if (length(missing)) seq.int(1L, length(missing), by = 256L) else integer()) {
    ids <- missing[first:min(length(missing), first + 255L)]
    z <- .mgcvst_pca_materialize(fit, ids, basis, threads, B = learned$B)
    block <- list(feature_index = ids, A = z$A, C = z$C, fro2 = z$fro2,
                  error = z$error)
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
    A = do.call(cbind, c(list(matrix(0, basis$rank, 0L)), lapply(blocks, `[[`, "A"))),
    C = do.call(rbind, c(list(matrix(0, 0L, rank)), lapply(blocks, `[[`, "C"))),
    fro2 = unlist(lapply(blocks, `[[`, "fro2")),
    error = unlist(lapply(blocks, `[[`, "error"))
  )
  rm(blocks)
  kt <- match(used, learned$train)
  kr <- match(used, proj$ids)
  tr <- !is.na(kt)
  A <- matrix(NA_real_, basis$rank, length(used))
  C <- matrix(NA_real_, length(used), ncol(learned$B))
  A[, tr] <- learned$A[, kt[tr]]
  A[, !tr] <- proj$A[, kr[!tr]]
  C[tr, ] <- learned$C[kt[tr], ]
  C[!tr, ] <- proj$C[kr[!tr], ]
  fro2 <- ifelse(tr, learned$fro2[kt], proj$fro2[kr])
  gene_error <- ifelse(tr, NA_character_, proj$error[kr])
  learned$A <- NULL
  rm(proj)
  t_project <- proc.time()[["elapsed"]] - t0

  tables <- mgcvst_pca_tables_cpp(learned$B, basis$rank, threads)
  learned$B <- NULL
  preparation_elapsed <- proc.time()[["elapsed"]] - started
  if (verbose) message("Computed PCAlearning trace tables in ",
                       round(tables$timing[["total"]], 2), " s.")

  # Compact rows of one block of pairs: local positions (li, lj) in `used`
  # and the kernel output matrix. A pair with a gene whose projection failed
  # has status 3; a non-finite log p-value has status 2.
  ok <- is.na(gene_error)
  pair_frame <- function(li, lj, out) {
    frame <- .mgcvst_pairs_frame(
      used[li], used[lj], out[, "U"], out[, "logp_two_sided"],
      out[, "logp_positive"], out[, "logp_negative"]
    )
    invalid <- !is.finite(frame$log_p_two_sided)
    frame$status[invalid] <- .mgcvst_pair_status[["p_value"]]
    frame[invalid, c("log_p_two_sided", "log_p_positive", "log_p_negative")] <-
      NA_real_
    bad <- !ok[li] | !ok[lj]
    frame$status[bad] <- .mgcvst_pair_status[["feature"]]
    frame[bad, c("score", "log_p_two_sided", "log_p_positive",
                 "log_p_negative")] <- NA_real_
    frame
  }

  t0 <- proc.time()[["elapsed"]]
  chunks <- 0L
  resumed_pairs <- 0
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
        out <- mgcvst_pca_pairs_block_cpp(A, C, tables, first, last, threads,
                                          moments = FALSE)
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
        out <- mgcvst_pca_pairs_cpp(A, C, tables, local[z, 1L], local[z, 2L],
                                    threads, moments = FALSE)
        .mgcvst_write_parquet(pair_frame(local[z, 1L], local[z, 2L], out),
                              .mgcvst_shard_file(pair_dir, first))
      }
      shard_files <- c(shard_files, .mgcvst_shard_file(pair_dir, first))
      shard_rows <- c(shard_rows, length(z))
      chunks <- chunks + 1L
    }
  }
  pair_elapsed <- proc.time()[["elapsed"]] - t0

  e2 <- fro2 - rowSums(C^2)
  genes <- data.frame(
    feature_id = fit$feature_id[used], feature_index = used,
    training = tr, cell = sampled$cell[used],
    sigma_g2 = scales$sigma_g2[used], sigma_e2 = scales$sigma_e2[used],
    tau = scales$tau[used], fro2 = fro2, e2 = e2,
    e2_relative = e2 / fro2, error_message = gene_error,
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
      liu_approximation = "pca_learning", preparation_backend = "sparse",
      pair_schedule = if (all_pairs) "pcalearning_gene_blocks" else
        "pcalearning_pair_list",
      chunks = chunks, resumed_pairs = resumed_pairs, pair_dir = pair_dir,
      preparation_elapsed = preparation_elapsed, contract = contract,
      pca_learning = list(
        rank = rank, n_per_cell = n_per_cell,
        seed = seed, q = basis$rank, training = training,
        checkpoint = list(path = path, basis_resumed = basis_resumed,
                          resumed_genes = resumed_genes,
                          projected_genes = length(missing)),
        gram_values = learned$values, rotation = learned$rotation,
        coefficients = C, genes = genes,
        table_timing = tables$timing, table_memory = tables$memory,
        elapsed = c(sample = t_sample, learned$elapsed, project = t_project,
                    tables = tables$timing[["total"]], pairs = pair_elapsed)
      )
    )
  )
}
