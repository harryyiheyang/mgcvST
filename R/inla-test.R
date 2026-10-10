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

.mgcvst_inla_wgcna_scores <- function(fit, used, threads, verbose) {
  group <- "global"
  fit <- .inlast_sparse_prepare(fit)
  basis <- .inlast_check_basis(fit, .inlast_sparse_observation_basis(fit))
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
