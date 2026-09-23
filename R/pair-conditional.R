# Conditional Cauchy pair testing for the sparse INLA score backend.
.mgcvst_conditional_test <- function(fit, pairs, q.value, threads,
                                    chunk_size, checkpoint_dir, resume,
                                    call) {
  if (!inherits(fit, "mgcvST_model_fit") ||
      !identical(fit$estimator, "INLA") ||
      !identical(fit$score_backend, "sparse")) {
    stop("pairwise_method = 'conditional' requires a sparse INLA fit.")
  }
  if (!is.numeric(q.value) || length(q.value) != 1L ||
      !is.finite(q.value) || q.value <= 0 || q.value > 1) {
    stop("q.value must be one finite value in (0, 1].")
  }
  if (is.null(threads)) threads <- 1L
  if (!is.numeric(threads) || length(threads) != 1L ||
      !is.finite(threads) || threads < 1 || threads != floor(threads)) {
    stop("threads must be one positive integer.")
  }
  threads <- as.integer(threads)
  if (is.null(chunk_size)) chunk_size <- max(16L, min(32L, threads * 4L))
  if (!is.numeric(chunk_size) || length(chunk_size) != 1L ||
      !is.finite(chunk_size) || chunk_size < 1 ||
      chunk_size != floor(chunk_size) || chunk_size > 64L) {
    stop("For conditional pairs, chunk_size must be 1 to 64 genes per batch.")
  }
  chunk_size <- as.integer(chunk_size)
  if (!is.logical(resume) || length(resume) != 1L || is.na(resume)) {
    stop("resume must be TRUE or FALSE.")
  }

  .mgcvst_thread_limit()
  RhpcBLASctl::blas_set_num_threads(1L)
  available <- .mgcvst_feature_available(fit)
  if (is.null(pairs)) {
    used <- which(available)
    if (length(used) < 2L) stop("At least two available INLA genes are required.")
    index <- mgcvst_conditional_all_pairs_cpp(length(used))
    universe <- "all_available_gene_pairs"
  } else {
    requested <- .mgcvst_pair_index(pairs, fit$feature_id)
    if (any(!available[requested])) {
      stop("Every requested conditional pair must have two available INLA genes.")
    }
    used <- sort(unique(as.vector(requested)))
    index <- matrix(match(requested, used), ncol = 2L)
    universe <- "explicit_requested_gene_pairs"
  }
  ids <- fit$feature_id[used]
  if (anyDuplicated(ids)) stop("Conditional feature IDs must be unique.")
  fit <- .inlast_sparse_prepare(fit)
  basis <- .inlast_sparse_observation_basis(fit, coverage = 0.995)
  q <- basis$rank
  G <- length(used)
  A <- matrix(NA_real_, q, G, dimnames = list(NULL, ids))
  starts <- seq.int(1L, G, by = 32L)
  for (first in starts) {
    rows <- first:min(G, first + 31L)
    states <- .inlast_sparse_batch(fit, used[rows], threads,
                                   score_only = TRUE)
    if (length(states) != length(rows)) {
      stop("Sparse INLA score preparation returned the wrong gene count.")
    }
    for (k in seq_along(rows)) {
      state <- states[[k]]
      if (!is.null(state$error)) {
        stop("Conditional score preparation failed for ", ids[rows[k]],
             ": ", state$error)
      }
      A[, rows[k]] <- as.numeric(crossprod(basis$coordinate, state$a))
    }
  }
  if (any(!is.finite(A))) stop("Conditional score coordinates are non-finite.")
  S <- crossprod(A)

  geometry_hash <- digest::digest(list(
    fit$score_sparse$A, fit$score_sparse$Q,
    fit$score_sparse$constraint, fit$geometry$nuisance_design,
    basis$coordinate, basis$basis
  ), algo = "sha256")
  working_hash <- vapply(used, function(j) {
    digest::digest(fit$working_variance[, j], algo = "sha256")
  }, character(1L))
  signature <- digest::digest(list(
    version = 1L, ids = ids, A = A, geometry_hash = geometry_hash,
    working_hash = working_hash,
    smoothing = fit$smoothing_parameters[used, , drop = FALSE],
    dispersion = fit$dispersion[used]
  ), algo = "sha256")

  temporary <- is.null(checkpoint_dir)
  if (temporary) checkpoint_dir <- tempfile("mgcvst-conditional-", tmpdir = tempdir())
  if (!is.character(checkpoint_dir) || length(checkpoint_dir) != 1L ||
      is.na(checkpoint_dir) || !nzchar(checkpoint_dir)) {
    stop("checkpoint_dir must be NULL or one non-empty directory path.")
  }
  if (!dir.exists(checkpoint_dir)) {
    if (!dir.create(checkpoint_dir, recursive = TRUE)) {
      stop("Could not create the conditional checkpoint directory.")
    }
  }
  checkpoint_dir <- normalizePath(checkpoint_dir, winslash = "/", mustWork = TRUE)
  if (temporary) on.exit({
    root <- normalizePath(tempdir(), winslash = "/", mustWork = TRUE)
    if (startsWith(checkpoint_dir, paste0(root, "/"))) {
      unlink(checkpoint_dir, recursive = TRUE)
    }
  }, add = TRUE)
  manifest_path <- file.path(checkpoint_dir, "manifest.rds")
  if (file.exists(manifest_path)) {
    if (!resume) stop("A conditional checkpoint already exists; use a new directory.")
    manifest <- readRDS(manifest_path)
    if (!identical(manifest$signature, signature) ||
        !identical(manifest$feature_id, ids)) {
      stop("The conditional checkpoint belongs to a different fit or gene set.")
    }
  } else {
    existing <- list.files(checkpoint_dir, all.files = TRUE, no.. = TRUE)
    if (length(existing)) stop("The conditional checkpoint directory is not empty.")
    saveRDS(list(signature = signature, feature_id = ids, rank = q),
            paste0(manifest_path, ".pending"))
    if (!file.rename(paste0(manifest_path, ".pending"), manifest_path)) {
      stop("Could not finalize the conditional checkpoint manifest.")
    }
  }

  path <- file.path(checkpoint_dir, sprintf("row-%05d.bin", seq_len(G)))
  done <- paste0(path, ".done")
  complete <- vapply(seq_len(G), function(k) {
    file.exists(done[k]) && file.exists(path[k]) &&
      identical(file.info(path[k])$size, as.numeric(8 * G)) &&
      identical(readLines(done[k], warn = FALSE), signature)
  }, logical(1L))
  missing <- which(!complete)
  if (length(missing)) for (first in seq.int(1L, length(missing), by = chunk_size)) {
    rows <- missing[first:min(length(missing), first + chunk_size - 1L)]
    units <- .inlast_sparse_units(fit, used[rows], threads = threads)
    if (length(units) != length(rows)) {
      stop("Sparse INLA unit preparation returned the wrong gene count.")
    }
    for (k in seq_along(rows)) if (!is.null(units[[k]]$error)) {
      stop("Conditional covariance preparation failed for ", ids[rows[k]],
           ": ", units[[k]]$error)
    }
    states <- .inlast_sparse_materialize_reduced(fit, units, basis,
                                                 threads = threads)
    if (length(states) != length(rows)) {
      stop("Sparse INLA materialization returned the wrong gene count.")
    }
    for (k in seq_along(rows)) if (!is.null(states[[k]]$error)) {
      stop("Conditional covariance materialization failed for ", ids[rows[k]],
           ": ", states[[k]]$error)
    }
    Vb <- mgcvst_conditional_variance_rows_cpp(
      A, lapply(states, `[[`, "M"), threads = threads, block_size = 512L
    )
    if (any(!is.finite(Vb)) || any(Vb <= 0)) {
      stop("A conditional variance is non-finite or non-positive.")
    }
    for (k in seq_along(rows)) {
      row <- rows[k]
      pending <- paste0(path[row], ".pending")
      if (file.exists(pending)) unlink(pending)
      con <- file(pending, "wb")
      writeBin(as.numeric(Vb[k, ]), con, size = 8L)
      close(con)
      if (file.exists(path[row])) unlink(path[row])
      if (!file.rename(pending, path[row])) {
        stop("Could not finalize conditional variance row ", row, ".")
      }
      done_pending <- paste0(done[row], ".pending")
      if (file.exists(done_pending)) unlink(done_pending)
      writeLines(signature, done_pending)
      if (file.exists(done[row])) unlink(done[row])
      if (!file.rename(done_pending, done[row])) {
        stop("Could not finalize conditional done marker ", row, ".")
      }
    }
    rm(units, states, Vb)
  }

  V <- matrix(NA_real_, G, G)
  for (k in seq_len(G)) {
    con <- file(path[k], "rb")
    V[k, ] <- readBin(con, numeric(), n = G, size = 8L)
    close(con)
  }
  if (any(!is.finite(V)) || any(V <= 0)) {
    stop("The conditional variance checkpoint contains invalid rows.")
  }
  pairs_out <- mgcvst_conditional_pairs_cpp(S, V, index)
  lp <- pairs_out$log_p
  if (anyNA(lp) || any(lp > 0)) stop("Cauchy calibration returned invalid log p-values.")
  m <- length(lp)
  Hm <- digamma(m + 1) - digamma(1)
  ord <- order(lp)
  raw <- lp[ord] + log(m) + log(Hm) - log(seq_len(m))
  sorted_by <- pmin(0, rev(cummin(rev(raw))))
  lby <- numeric(m)
  lby[ord] <- sorted_by
  results <- data.frame(
    feature1 = ids[index[, 1L]], feature2 = ids[index[, 2L]],
    S = pairs_out$score, p = exp(lp), log_p = lp,
    p_BY = exp(lby), log_p_BY = lby,
    BY_reject = lby <= log(q.value),
    stringsAsFactors = FALSE
  )
  discoveries <- sum(results$BY_reject)
  structure(list(
    results = results,
    threshold = list(q_value = q.value, FDR = TRUE,
                     adjustment_method = "BY",
                     raw_p_threshold = if (discoveries) max(results$p[results$BY_reject]) else NA_real_),
    discoveries = list(pairs_requested = m, pairs_tested = m,
                       pairs_with_p_value = m, pairs_discovered = discoveries,
                       pairs_discovered_positive = NA_integer_,
                       pairs_discovered_negative = NA_integer_,
                       pairs_highlighted = 0L, pairs_retained = discoveries),
    pair_contract = universe,
    test_definition = "conditional_normal_directions_cauchy_combination",
    calibration = "conditional_cauchy",
    timing = list(threads = threads, score_rank = q,
                  score_genes = G, variance_rows_reused = sum(complete),
                  checkpoint_dir = if (temporary) NULL else checkpoint_dir),
    call = call
  ), class = "mgcvST_test")
}
