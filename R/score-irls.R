# Reconstruct and validate the training linear-predictor matrix.
.gam_training_lpmatrix <- function(fit) {
  X <- as.matrix(mgcv::predict.gam(fit, type = "lpmatrix"))
  if (nrow(X) != length(fit$linear.predictors) || any(!is.finite(X))) {
    stop("Could not reconstruct a finite training lpmatrix from fit.")
  }
  X
}

