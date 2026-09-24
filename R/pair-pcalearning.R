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
.mgcvst_pca_scales <- function(fit) {
  G <- length(fit$feature_id)
  nb <- fit$diagnostics$family_used == "negative_binomial"
  theta <- vapply(fit$family_parameters, function(x) {
    if (length(x)) x[1L] else NA_real_
  }, numeric(1L))
  mu_bar <- numeric(G)
  for (first in seq.int(1L, G, by = 500L)) {
    ids <- first:min(first + 499L, G)
    it <- ifelse(nb[ids], 1 / theta[ids], 0)
    D <- sweep(fit$working_variance[, ids, drop = FALSE], 2L, it, "-")
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
# c_j = (G R)_j / tau_j.
.mgcvst_pca_basis <- function(fit, train, tau, basis, rank, threads) {
  t0 <- proc.time()[["elapsed"]]
  z <- .mgcvst_pca_materialize(fit, train, basis, threads, pack = TRUE)
  ok <- is.na(z$error)
  t1 <- proc.time()[["elapsed"]]
  gram <- mgcvst_pca_gram_cpp(z$packed[ok], tau[ok], threads)
  t2 <- proc.time()[["elapsed"]]
  eg <- eigen(gram, symmetric = TRUE)
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

# Approximate Liu log p-values for all requested pairs from the PCAlearning
# basis. An all-pairs universe streams gene blocks; other pair lists are
# evaluated by (i, j).
.mgcvst_pair_pcalearning <- function(fit, index, pair_index, threads,
                                     chunk_size, verbose, basis, rank = 10L,
                                     n_per_cell = 3L, seed = 1L) {
  started <- proc.time()[["elapsed"]]
  threads <- as.integer(threads)
  fit <- .inlast_sparse_prepare(fit)
  universe <- which(.mgcvst_feature_available(fit))
  scales <- .mgcvst_pca_scales(fit)
  sampled <- .mgcvst_pca_training(scales, universe, n_per_cell, seed)
  t_sample <- proc.time()[["elapsed"]] - started
  if (verbose) message("Sampled ", length(sampled$train),
                       " PCAlearning training genes.")

  learned <- .mgcvst_pca_basis(fit, sampled$train, scales$tau[sampled$train],
                               basis, as.integer(rank), threads)
  if (verbose) message("Learned the rank-", rank, " PCAlearning basis.")

  # Training genes reuse their first materialization; the others are projected.
  used <- sort(unique(as.vector(index)))
  t0 <- proc.time()[["elapsed"]]
  rest <- setdiff(used, learned$train)
  proj <- .mgcvst_pca_materialize(fit, rest, basis, threads, B = learned$B)
  kt <- match(used, learned$train)
  kr <- match(used, rest)
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

  n_used <- length(used)
  local <- matrix(match(index, used), ncol = 2L)
  lo <- pmin(local[, 1L], local[, 2L])
  hi <- pmax(local[, 1L], local[, 2L])
  key <- (lo - 1) * n_used + hi
  n <- nrow(index)
  score <- rep(NA_real_, n)
  log_p <- matrix(NA_real_, n, 3L)
  columns <- c("logp_two_sided", "logp_positive", "logp_negative")
  t0 <- proc.time()[["elapsed"]]
  chunks <- 0L
  all_pairs <- n == n_used * (n_used - 1) / 2 && all(lo < hi) &&
    !anyDuplicated(key)
  if (all_pairs) {
    target <- order(key)
    per_gene <- n_used - seq_len(n_used)
    done <- 0
    first <- 1L
    while (first < n_used) {
      last <- first
      total <- per_gene[first]
      while (last + 1L < n_used && total + per_gene[last + 1L] <= chunk_size) {
        last <- last + 1L
        total <- total + per_gene[last]
      }
      out <- mgcvst_pca_pairs_block_cpp(A, C, tables, first, last, threads,
                                        moments = FALSE)
      rows <- target[done + seq_len(nrow(out))]
      score[rows] <- out[, "U"]
      log_p[rows, ] <- out[, columns]
      done <- done + nrow(out)
      chunks <- chunks + 1L
      first <- last + 1L
    }
  } else {
    for (first in seq.int(1L, n, by = chunk_size)) {
      z <- first:min(n, first + chunk_size - 1L)
      out <- mgcvst_pca_pairs_cpp(A, C, tables, local[z, 1L], local[z, 2L],
                                  threads, moments = FALSE)
      score[z] <- out[, "U"]
      log_p[z, ] <- out[, columns]
      chunks <- chunks + 1L
    }
  }
  pair_elapsed <- proc.time()[["elapsed"]] - t0

  error_message <- rep(NA_character_, n)
  error_message[!is.finite(log_p[, 1L])] <-
    "PCAlearning Liu calibration returned an invalid p-value."
  ok <- is.na(gene_error)
  for (j in which(!ok[local[, 1L]] | !ok[local[, 2L]])) {
    k <- local[j, ][!ok[local[j, ]]]
    error_message[j] <- paste(paste0(fit$feature_id[used[k]], ": ",
                                     gene_error[k]), collapse = " | ")
  }
  bad <- !is.na(error_message)
  score[bad] <- NA_real_
  log_p[bad, ] <- NA_real_

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
  result <- data.frame(
    pair_index = pair_index, score = score, information = NA_real_,
    effective_rank = NA_real_, p_value = exp(log_p[, 1L]),
    log_p_two_sided = log_p[, 1L], log_p_positive = log_p[, 2L],
    log_p_negative = log_p[, 3L], error_message = error_message,
    stringsAsFactors = FALSE
  )
  list(
    result = result, elapsed = pair_elapsed,
    metadata = list(
      approximate = "PCAlearning", preparation_backend = "sparse",
      pair_schedule = if (all_pairs) "pcalearning_gene_blocks" else
        "pcalearning_pair_list",
      chunks = chunks, preparation_elapsed = preparation_elapsed,
      pca_learning = list(
        rank = as.integer(rank), n_per_cell = as.integer(n_per_cell),
        seed = as.integer(seed), q = basis$rank, training = training,
        gram_values = learned$values, rotation = learned$rotation,
        coefficients = C, genes = genes,
        table_timing = tables$timing, table_memory = tables$memory,
        elapsed = c(sample = t_sample, learned$elapsed, project = t_project,
                    tables = tables$timing[["total"]], pairs = pair_elapsed)
      )
    )
  )
}
