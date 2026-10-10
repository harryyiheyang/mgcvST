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

# The pair route of inlaST.test(): PCAlearning-approximated Liu moments on the
# constrained observation-kernel basis, streamed to raw shards.
.mgcvst_inla_test_pairs <- function(fit, index, threads, chunk_size, verbose,
                                    basis = NULL,
                                    rank = .mgcvst_pca_defaults$rank,
                                    n_per_cell = .mgcvst_pca_defaults$n_per_cell,
                                    seed = .mgcvst_pca_defaults$seed,
                                    checkpoint_dir = NULL, resume = TRUE) {
  fit <- .inlast_sparse_prepare(fit)
  if (is.null(basis)) basis <- .inlast_sparse_observation_basis(fit)
  evaluated <- .mgcvst_pair_pcalearning(
    fit, index, threads, chunk_size, verbose, basis = basis,
    rank = rank, n_per_cell = n_per_cell, seed = seed,
    checkpoint_dir = checkpoint_dir, resume = resume
  )
  evaluated$metadata <- c(list(
    q = ncol(fit$score_sparse$Q), r = basis$rank,
    target_coverage = basis$coverage, kept_coverage = basis$kept,
    tail = basis$tail,
    basis = "constrained_observation_kernel_A_Qg_inverse_At",
    unit_cache = "score_state_shards"
  ), evaluated$metadata)
  evaluated
}

.mgcvst_inla_wgcna_scores <- function(fit, used, threads, verbose) {
  group <- "global"
  fit <- .inlast_sparse_prepare(fit)
  basis <- .inlast_sparse_observation_basis(fit)
  A <- crossprod(basis$coordinate, fit$score_a[, used, drop = FALSE])
  colnames(A) <- fit$feature_id[used]
  if (any(!is.finite(A))) {
    stop("The fit does not retain valid sparse INLA score vectors for: ",
         paste(fit$feature_id[used[!is.finite(colSums(A))]], collapse = ", "), ".")
  }
  width <- stats::setNames(nrow(A), "global")
  normalization <- basis$rank
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
