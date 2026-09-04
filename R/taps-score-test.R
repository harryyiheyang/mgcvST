# Package-local fitted-model TAPS, adapted from the source credited in inst/COPYRIGHTS.
# The spectrum and calibration arithmetic live in marginal-taps.R / marginal-api.R.
taps_score_test <- function(fit, test.component = 1L, null.tol = 1e-10,
                            method = "davies", max_eps = 1e-8, max_iter = 1e5,
                            n_threads = 1L, lpmatrix = NULL) {
  if (!inherits(fit, "gam")) stop("fit must be a fitted gam object.")
  .working_family_id(fit$family$family)
  method <- match.arg(method, c("davies", "liu"))
  if (length(test.component) != 1L || !is.finite(test.component) ||
      test.component != as.integer(test.component) || test.component < 1L ||
      test.component > length(fit$smooth)) {
    stop("test.component must identify a fitted smooth term.")
  }
  if (!isTRUE(n_threads == 1L)) stop("Inline marginal tests require n_threads = 1.")
  X <- lpmatrix
  if (is.null(X)) X <- fit$.taps_score_X
  if (is.null(X)) X <- .gam_training_lpmatrix(fit)
  if (!is.matrix(X) || !is.numeric(X) || any(!is.finite(X)) ||
      nrow(X) != NROW(fit$linear.predictors) || ncol(X) != length(fit$coefficients)) {
    stop("The training lpmatrix must be finite and aligned with the fitted rows and coefficients.")
  }
  if (!is.null(colnames(X)) && !is.null(names(fit$coefficients)) &&
      !identical(colnames(X), names(fit$coefficients))) {
    stop("The training lpmatrix columns must follow the fitted coefficient order.")
  }
  geometry <- .mgcvst_marginal_geometry(fit, X, test.component)
  z <- .mgcvst_marginal_spectrum(fit, geometry, null.tol)
  if (!length(z$lambda) || !is.finite(z$statistic) || any(!is.finite(z$lambda))) {
    stop("The marginal score has no finite positive mixture spectrum.")
  }
  if (method == "davies") {
    result <- .mgcvst_marginal_davies(z, "liu", max_eps, max_iter)
    p <- result$p_value
    method <- result$method_used
  } else {
    p <- .mgcvst_marginal_liu(z$statistic, .mgcvst_marginal_moments(z$lambda))
  }
  data.frame(smooth.term = z$smooth.term, smooth.pvalue = p, method = method)
}
