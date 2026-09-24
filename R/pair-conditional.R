# Benjamini-Yekutieli step-up on natural-log p-values; returns log adjusted
# p-values, so tails below the double range keep their ordering and decisions.
.mgcvst_log_by <- function(lp) {
  m <- length(lp)
  Hm <- digamma(m + 1) - digamma(1)
  ord <- order(lp)
  raw <- lp[ord] + log(m) + log(Hm) - log(seq_len(m))
  out <- numeric(m)
  out[ord] <- pmin(0, rev(cummin(rev(raw))))
  out
}

# Conditional Cauchy pair testing for the sparse INLA score backend.
.mgcvst_conditional_test <- function(fit, pairs, q.value, FDR, method, threads,
                                    chunk_size, checkpoint_dir, resume,
                                    conditional_precision, call) {
  if (!inherits(fit, "mgcvST_model_fit") ||
      !identical(fit$estimator, "INLA") ||
      !identical(fit$score_backend, "sparse")) {
    stop("pairwise_method = 'conditional_cauchy' requires a sparse INLA fit.")
  }
  if (!is.numeric(q.value) || length(q.value) != 1L ||
      !is.finite(q.value) || q.value <= 0 || q.value > 1) {
    stop("q.value must be one finite value in (0, 1].")
  }
  if (!is.logical(FDR) || length(FDR) != 1L || is.na(FDR)) {
    stop("FDR must be TRUE or FALSE.")
  }
  if (!is.character(method) || length(method) != 1L || is.na(method) ||
      !(method %in% stats::p.adjust.methods)) {
    stop("method must be one of stats::p.adjust.methods.")
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
  if (!is.character(conditional_precision) ||
      length(conditional_precision) != 1L ||
      !conditional_precision %in% c("double", "float32")) {
    stop("conditional_precision must be 'double' or 'float32'.")
  }
  float32 <- identical(conditional_precision, "float32")

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
  t_basis <- proc.time()[["elapsed"]]
  basis <- .inlast_sparse_observation_basis(fit, coverage = 0.995)
  basis_elapsed <- proc.time()[["elapsed"]] - t_basis
  q <- basis$rank
  G <- length(used)
  A <- crossprod(basis$coordinate, fit$score_a[, used, drop = FALSE])
  dimnames(A) <- list(NULL, ids)
  score_elapsed <- 0
  if (any(!is.finite(A))) stop("Conditional score coordinates are non-finite.")
  S <- crossprod(A)

  existing_manifest <- !is.null(checkpoint_dir) && dir.exists(checkpoint_dir) &&
    file.exists(file.path(checkpoint_dir, "manifest.rds"))
  checkpoint_manifest <- if (existing_manifest) {
    readRDS(file.path(checkpoint_dir, "manifest.rds"))
  } else NULL
  signature_data <- list(
    version = 3L, ids = ids,
    compact = .inlast_compact_signature(fit, used, basis),
    score_rank = q, sp_index = fit$score_sparse$sp_index,
    nuisance_precision = .inlast_sparse_nuisance_precision(fit, used),
    coverage = 0.995, conditional_precision = conditional_precision
  )
  signature <- digest::digest(signature_data, algo = "sha256")

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
    if (!identical(checkpoint_manifest$signature, signature) ||
        !identical(checkpoint_manifest$feature_id, ids)) {
      stop("The conditional checkpoint belongs to a different fit or gene set.")
    }
  } else {
    existing <- list.files(checkpoint_dir, all.files = TRUE, no.. = TRUE)
    if (length(existing)) stop("The conditional checkpoint directory is not empty.")
    saveRDS(list(version = 3L, signature = signature, feature_id = ids,
                 rank = q),
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
  materialize_elapsed <- 0
  variance_elapsed <- 0
  if (length(missing)) for (first in seq.int(1L, length(missing), by = chunk_size)) {
    rows <- missing[first:min(length(missing), first + chunk_size - 1L)]
    t_units <- proc.time()[["elapsed"]]
    units <- .inlast_sparse_units(fit, used[rows], threads = threads)
    score_elapsed <- score_elapsed + proc.time()[["elapsed"]] - t_units
    if (length(units) != length(rows)) {
      stop("Sparse INLA unit preparation returned the wrong gene count.")
    }
    for (k in seq_along(rows)) if (!is.null(units[[k]]$error)) {
      stop("Conditional covariance preparation failed for ", ids[rows[k]],
           ": ", units[[k]]$error)
    }
    t_materialize <- proc.time()[["elapsed"]]
    states <- .inlast_sparse_materialize_reduced(fit, units, basis,
                                                 threads = threads)
    materialize_elapsed <- materialize_elapsed +
      proc.time()[["elapsed"]] - t_materialize
    if (length(states) != length(rows)) {
      stop("Sparse INLA materialization returned the wrong gene count.")
    }
    for (k in seq_along(rows)) if (!is.null(states[[k]]$error)) {
      stop("Conditional covariance materialization failed for ", ids[rows[k]],
           ": ", states[[k]]$error)
    }
    t_variance <- proc.time()[["elapsed"]]
    Vb <- mgcvst_conditional_variance_rows_cpp(
      A, lapply(states, `[[`, "M"), threads = threads, block_size = 512L,
      float32 = float32
    )
    variance_elapsed <- variance_elapsed + proc.time()[["elapsed"]] - t_variance
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
  s <- pairs_out$score
  lp12 <- log(2) + stats::pnorm(abs(s) / sqrt(V[index]), lower.tail = FALSE,
                                log.p = TRUE)
  lp21 <- log(2) + stats::pnorm(abs(s) / sqrt(V[index[, 2:1, drop = FALSE]]),
                                lower.tail = FALSE, log.p = TRUE)
  results <- data.frame(
    feature1 = ids[index[, 1L]], feature2 = ids[index[, 2L]],
    signed_score = s, statistic = s^2,
    p_two_sided = exp(lp), log_p_two_sided = lp,
    p_1_given_2 = exp(lp12), log_p_1_given_2 = lp12,
    p_2_given_1 = exp(lp21), log_p_2_given_1 = lp21,
    stringsAsFactors = FALSE
  )
  m <- nrow(results)
  # Multiple testing is optional; BY stays on the log scale so that tails
  # below the double range keep their ordering and decisions.
  la <- if (!FDR) lp else if (method == "BY") .mgcvst_log_by(lp) else
    log(stats::p.adjust(exp(lp), method))
  results$p_adjusted <- exp(la)
  results$log_p_adjusted <- la
  results$discovered <- la <= log(q.value)
  discoveries <- sum(results$discovered)
  structure(list(
    results = results,
    threshold = list(q_value = q.value, FDR = FDR,
                     adjustment_method = if (FDR) method else "none",
                     raw_p_threshold = if (discoveries)
                       max(results$p_two_sided[results$discovered]) else NA_real_),
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
                  conditional_precision = conditional_precision,
                  basis_elapsed = basis_elapsed,
                  score_unit_elapsed = score_elapsed,
                  reduced_materialize_elapsed = materialize_elapsed,
                  variance_elapsed = variance_elapsed,
                  checkpoint_dir = if (temporary) NULL else checkpoint_dir),
    call = call
  ), class = "mgcvST_test")
}
