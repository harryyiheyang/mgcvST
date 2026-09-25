# Factor one positive-definite fixed-kappa SPDE covariance.
.mgcvst_spde_factor <- function(B, Q, scale) {
  sqrt(scale) * .mgcvst_spde_factor_base(B, Q)
}

.mgcvst_spde_factor_base <- function(B, Q) {
  Q <- Matrix::forceSymmetric(Matrix::Matrix(Q, sparse = TRUE))
  R <- Matrix::chol(Q)
  Rinv <- Matrix::solve(R, Matrix::Diagonal(nrow(Q)))
  .magic_mm(B, as.matrix(Rinv))
}

# Per-test geometry only. Cache conditions without moving failures out of
# the original per-pair tryCatch or changing feature-validation precedence.
.mgcvst_model_fixed_factors <- function(fit) {
  if (identical(fit$score_backend, "sparse")) {
    return(vector("list", length(fit$geometry$smooth)))
  }
  lapply(fit$geometry$smooth, function(s) {
    if (s$fixed || is.null(s$score_component) || length(s$penalties) != 1L ||
        length(s$sp_index) != 1L) return(NULL)
    tryCatch(.mgcvst_spde_factor_base(s$B, s$penalties[[1L]]),
             error = function(e) e)
  })
}
