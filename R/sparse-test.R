.mgcvst_wgcna_scores <- function(fit, used, threads, verbose) {
  group <- "global"
  fit <- .spdest_sparse_prepare(fit)
  basis <- .spdest_check_basis(fit, .spdest_sparse_observation_basis(fit))
  A <- crossprod(basis$coordinate, fit$score_a[, used, drop = FALSE])
  colnames(A) <- fit$feature_id[used]
  if (any(!is.finite(A))) {
    stop("The fit does not retain valid sparse score vectors for: ",
         paste(fit$feature_id[used[!is.finite(colSums(A))]], collapse = ", "), ".")
  }
  width <- stats::setNames(nrow(A), "global")
  normalization <- basis$rank
  if (verbose) {
    message("Constructed scores for ", length(used),
            " features from the saved sparse scores.")
  }
  list(
    A = A, group = group, width = width,
    normalization = normalization,
    feature_id = fit$feature_id[used]
  )
}
