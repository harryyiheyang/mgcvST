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
                                    resume = TRUE, approximate = FALSE,
                                    n_ref = 100L,
                                    ref_method = c("random", "score", "hyper"),
                                    ref_seed = 1L, ref_tol = 1e-6,
                                    diagnostic_pairs = 0L) {
  fit <- .inlast_sparse_prepare(fit)
  if (is.null(basis)) basis <- .inlast_sparse_observation_basis(
    fit, coverage = coverage, full_rank = full_rank
  )
  evaluate <- if (approximate) .mgcvst_pair_approximate else .mgcvst_pair_pipeline
  args <- list(
    fit, index, pair_index, threads, chunk_size, verbose, basis = basis,
    cache_bytes = cache_bytes, checkpoint_dir = checkpoint_dir, resume = resume
  )
  if (approximate) args <- c(args, list(n_ref = n_ref, ref_method = ref_method,
    ref_seed = ref_seed, ref_tol = ref_tol,
    diagnostic_pairs = diagnostic_pairs))
  evaluated <- do.call(evaluate, args)
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
  fit <- .inlast_sparse_prepare(fit)
  blocks <- split(seq_along(used), ceiling(seq_along(used) / 32L))
  states <- vector("list", length(used))
  for (rows in blocks) {
    states[rows] <- .inlast_sparse_batch(
      fit, features = used[rows], threads = threads, score_only = TRUE
    )
  }
  if (length(states) != length(used)) {
    stop("The sparse INLA batch returned an incompatible feature count.")
  }
  failed <- vapply(states, function(z) {
    !is.null(z$error) && length(z$error) == 1L && !is.na(z$error) && nzchar(z$error)
  }, logical(1L))
  if (any(failed)) {
    stop("Sparse INLA score construction failed: ", paste(
      paste0(fit$feature_id[used[failed]], ": ",
             vapply(states[failed], `[[`, character(1L), "error")),
      collapse = " | "
    ))
  }
  coordinate_width <- length(states[[1L]]$a)
  normalization <- as.integer(states[[1L]]$normalization)
  width <- stats::setNames(coordinate_width, "global")
  A <- do.call(cbind, lapply(states, `[[`, "a"))
  colnames(A) <- fit$feature_id[used]
  if (nrow(A) != coordinate_width || any(!is.finite(A))) {
    stop("The sparse INLA batch returned invalid score coordinates.")
  }
  if (verbose) {
    message("Constructed scores for ", length(used), " features with C++ OpenMP.")
  }
  list(
    A = A, group = group, width = width,
    normalization = normalization,
    feature_id = fit$feature_id[used]
  )
}
