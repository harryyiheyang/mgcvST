# Sparse INLA downstream evaluation uses one OpenMP layer in the manager.
.mgcvst_inla_downstream <- function(fit) {
  identical(fit$estimator, "INLA")
}

.mgcvst_inla_require_sparse <- function(fit) {
  if (!identical(fit$score_backend, "sparse")) {
    stop("INLA downstream tests require the single-global sparse score ",
         "geometry built by inlaST.set()/inlaST.estimate(); the dense INLA ",
         "score no longer exists.")
  }
  invisible(NULL)
}

.mgcvst_inla_pair_chunk_size <- function(fit, memory_bytes = NULL,
                                         basis = NULL) {
  if (is.null(basis)) basis <- .inlast_sparse_observation_basis(fit)
  if (is.null(memory_bytes)) memory_bytes <- .mgcvst_inla_memory_plan(
    fit, basis, pairs = 1L, threads = 1L)$cache_bytes
  r <- basis$rank
  as.integer(max(1L, min(128L, floor(memory_bytes / (2 * 8 * r^2)))))
}

.mgcvst_inla_serial_backend <- function(BPPARAM) {
  if (!inherits(BPPARAM, "SerialParam")) {
    stop(
      "Sparse INLA downstream tests use C++ OpenMP; BPPARAM must be SerialParam()."
    )
  }
  invisible(NULL)
}

.mgcvst_inla_test_pairs <- function(fit, index, pair_index, threads,
                                    chunk_size, verbose, coverage = 0.995,
                                    full_rank = FALSE, basis = NULL,
                                    cache_bytes = NULL, checkpoint_dir = NULL,
                                    resume = TRUE, liu_approximation = "exact",
                                    rank = 10L, n_per_cell = 3L, seed = 1L) {
  fit <- .inlast_sparse_prepare(fit)
  if (is.null(basis)) basis <- .inlast_sparse_observation_basis(
    fit, coverage = coverage, full_rank = full_rank
  )
  evaluated <- if (liu_approximation == "pca_learning") {
    .mgcvst_pair_pcalearning(
      fit, index, pair_index, threads, chunk_size, verbose, basis = basis,
      rank = rank, n_per_cell = n_per_cell, seed = seed,
      checkpoint_dir = checkpoint_dir, resume = resume
    )
  } else {
    .mgcvst_pair_pipeline(
      fit, index, pair_index, threads, chunk_size, verbose, basis = basis,
      cache_bytes = cache_bytes, checkpoint_dir = checkpoint_dir, resume = resume
    )
  }
  out <- evaluated$result
  names(out)[names(out) == "score"] <- "signed_score"
  names(out)[names(out) == "p_value"] <- "p_two_sided"
  valid <- is.finite(out$p_two_sided) & out$p_two_sided >= 0 &
    out$p_two_sided <= 1
  out$p_positive <- out$p_negative <- NA_real_
  out$p_positive[valid] <- ifelse(out$signed_score[valid] >= 0,
    out$p_two_sided[valid] / 2, 1 - out$p_two_sided[valid] / 2)
  out$p_negative[valid] <- ifelse(out$signed_score[valid] <= 0,
    out$p_two_sided[valid] / 2, 1 - out$p_two_sided[valid] / 2)
  attr(out, "inla_pairwise") <- c(list(
    q = ncol(fit$score_sparse$Q), r = basis$rank,
    target_coverage = basis$coverage, kept_coverage = basis$kept,
    tail = basis$tail,
    basis = "constrained_observation_kernel_A_Qg_inverse_At",
    unit_cache = "score_state_shards"
  ), evaluated$metadata)
  list(result = out, elapsed = evaluated$elapsed)
}
.mgcvst_inla_wgcna_scores <- function(fit, used, threads, verbose) {
  group <- "global"
  A <- fit$score_a[, used, drop = FALSE]
  colnames(A) <- fit$feature_id[used]
  if (any(!is.finite(A))) {
    stop("The fit does not retain valid sparse INLA score vectors for: ",
         paste(fit$feature_id[used[!is.finite(colSums(A))]], collapse = ", "), ".")
  }
  width <- stats::setNames(nrow(A), "global")
  normalization <- as.integer(fit$score_sparse$normalization)
  if (verbose) {
    message("Constructed scores for ", length(used),
            " features from the saved sparse INLA scores.")
  }
  list(
    A = A, group = group, width = width,
    normalization = normalization,
    feature_id = fit$feature_id[used]
  )
}
