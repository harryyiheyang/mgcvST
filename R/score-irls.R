# Map an mgcv family label to the supported working-model identifier.
.working_family_id <- function(family_name) {
  x <- tolower(family_name)
  if (grepl("negative binomial", x, fixed = TRUE)) return("negative_binomial")
  if (grepl("^tweedie", x)) return("tweedie")
  if (x %in% c("binomial", "quasibinomial")) return(x)
  if (x %in% c("poisson", "quasipoisson", "gaussian", "gamma")) return(x)
  stop(
    "mgcvST IRLS extraction currently supports Gaussian, Poisson, ",
    "quasi-Poisson, negative-binomial, Tweedie, binomial, ",
    "quasi-binomial, and Gamma fits."
  )
}

# Reconstruct and validate the training linear-predictor matrix.
.gam_training_lpmatrix <- function(fit) {
  X <- as.matrix(mgcv::predict.gam(fit, type = "lpmatrix"))
  if (nrow(X) != length(fit$linear.predictors) || any(!is.finite(X))) {
    stop("Could not reconstruct a finite training lpmatrix from fit.")
  }
  X
}

# Extract one full-rank smooth basis, penalty, and smoothing parameter.
.gam_single_smooth <- function(fit, smooth_index, L) {
  if (length(fit$smooth) != 1L) {
    stop(
      "The current IRLS pair API requires exactly one penalized smooth per fit; ",
      "additional nuisance smooth covariances are not yet implemented."
    )
  }
  smooth_index <- as.integer(smooth_index)
  if (length(smooth_index) != 1L || is.na(smooth_index) || smooth_index != 1L) {
    stop("smooth_index must identify the fit's single smooth term.")
  }
  s <- fit$smooth[[smooth_index]]
  if (length(s$S) != 1L || is.null(s$first.sp) || is.null(s$last.sp) ||
      s$first.sp != s$last.sp) {
    stop("The tested smooth must have exactly one full-rank penalty and one smoothing parameter.")
  }
  sp_values <- fit$sp
  if (!is.null(fit$full.sp) && length(fit$full.sp) >= s$last.sp) {
    sp_values <- fit$full.sp
  }
  sp <- as.numeric(sp_values[s$first.sp])
  if (length(sp) != 1L || !is.finite(sp) || sp <= 0) {
    stop("The tested smooth has no positive finite fitted smoothing parameter.")
  }
  if (inherits(s, "spdePC.smooth")) {
    B <- s$score_basis
    Q <- s$pc_score_Q
    if (is.null(B) || is.null(Q)) {
      stop("The spdePC smooth does not retain its full-coordinate score geometry.")
    }
    B <- .as_numeric_matrix(B, "spdePC score basis")
    Q <- Matrix::Matrix(Q, sparse = TRUE)
    if (nrow(B) != nrow(L) || ncol(B) != nrow(Q) || nrow(Q) != ncol(Q)) {
      stop("The spdePC full-coordinate score geometry is incompatible with the fit.")
    }
    score.psd <- TRUE
    idx <- seq.int(s$first.para, s$last.para)
  } else {
    idx <- seq.int(s$first.para, s$last.para)
    B <- L[, idx, drop = FALSE]
    Q <- Matrix::Matrix(s$S[[1L]], sparse = TRUE)
    score.psd <- FALSE
  }
  list(
    label = s$label, B = B, Q = Q, sp = sp, columns = idx,
    score_precision_psd = score.psd
  )
}

#' Extract the final IRLS working model from an mgcv fit
#'
#' Uses the generic final-PIRLS identities
#' \deqn{z=\eta+(y-\mu)g'(\mu),\qquad
#' W_i=w_i^{prior}/\{V(\mu_i)g'(\mu_i)^2\},}
#' and returns the diagonal working variance `phi / W`. For binomial matrix
#' responses, `mgcv` stores proportions in `fit$y` and trials in the prior
#' weights, so the same expression applies.
#'
#' @param fit A converged `mgcv::gam` fit.
#' @return A list with working response, working variance, offset, family
#'   metadata, and convergence information.
#' @export
rkhs_extract_working_model <- function(fit) {
  if (!inherits(fit, "gam")) stop("fit must inherit from class 'gam'.")
  if (isFALSE(fit$converged)) stop("fit did not converge.")
  if (!is.null(fit$mgcv.conv) && isFALSE(fit$mgcv.conv)) {
    stop("fit smoothing-parameter optimization did not converge.")
  }
  family_id <- .working_family_id(fit$family$family)

  eta <- as.numeric(fit$linear.predictors)
  mu <- as.numeric(fit$fitted.values)
  y <- as.numeric(fit$y)
  n <- length(eta)
  if (length(mu) != n || length(y) != n) {
    stop("fit response, fitted mean, and linear predictor have incompatible lengths.")
  }
  prior <- fit$prior.weights
  if (is.null(prior)) prior <- rep(1, n)
  prior <- as.numeric(prior)
  mu_eta <- as.numeric(fit$family$mu.eta(eta))
  variance <- as.numeric(fit$family$variance(mu))
  if (length(prior) != n || length(mu_eta) != n || length(variance) != n ||
      any(!is.finite(c(eta, mu, y, prior, mu_eta, variance))) ||
      any(prior <= 0) || any(mu_eta == 0) || any(variance <= 0)) {
    stop("The final IRLS response, derivative, variance, or prior weights are invalid.")
  }

  g_prime <- 1 / mu_eta
  z <- eta + (y - mu) * g_prime
  W <- prior / (variance * g_prime^2)
  phi <- fit$sig2
  if (is.null(phi)) phi <- 1
  phi <- as.numeric(phi)
  if (length(phi) != 1L || !is.finite(phi) || phi <= 0 ||
      any(!is.finite(z)) || any(!is.finite(W)) || any(W <= 0)) {
    stop("The final IRLS working response, weights, or dispersion are invalid.")
  }
  offset <- fit$offset
  if (is.null(offset)) offset <- rep(0, n)
  offset <- as.numeric(offset)
  if (length(offset) != n || any(!is.finite(offset))) {
    stop("The fitted offset is invalid.")
  }

  theta <- NULL
  if (is.function(fit$family$getTheta)) theta <- fit$family$getTheta(TRUE)
  list(
    pseudo_response = z,
    working_weight = W / phi,
    working_variance = phi / W,
    dispersion = phi,
    offset = offset,
    working_error = z - offset,
    family = family_id,
    family_label = fit$family$family,
    family_parameters = theta,
    converged = fit$converged
  )
}

#' Conditional IRLS RKHS covariance score for two mgcv fits
#'
#' Extracts each fit's final working Gaussian model, requires one shared
#' full-rank smooth penalty `Q`, and calls [rkhs_covariance_score()]. The fitted
#' smoothing parameters enter only as feature-specific marginal field scales
#' `phi / sp`; the unscaled `Q` and coefficient coordinates must agree between
#' features.
#'
#' This is a conditional working-Gaussian calibration. It does not account for
#' estimation of smoothing, dispersion, negative-binomial shape, or Tweedie
#' power parameters. Binomial small-sample and sparse-count applications may
#' require parametric-bootstrap calibration.
#'
#' @param fit1,fit2 Converged `mgcv::gam` fits on the same rows
#'   and with the same one-smooth setup.
#' @param smooth_index1,smooth_index2 Index of the shared full-rank smooth in
#'   each fit. The current prototype requires it to be the sole smooth.
#' @param method Calibration method passed to [rkhs_covariance_score()].
#' @param geometry_tol Relative tolerance for comparing `B` and unscaled `Q`.
#' @return An `rkhs_covariance_score` object augmented with IRLS metadata.
#' @export
rkhs_covariance_score_irls <- function(
    fit1, fit2, smooth_index1 = 1L, smooth_index2 = 1L,
    method = c("liu", "davies"), geometry_tol = 1e-9) {
  method <- match.arg(method)
  geometry_tol <- as.numeric(geometry_tol)
  if (length(geometry_tol) != 1L || !is.finite(geometry_tol) ||
      geometry_tol <= 0) {
    stop("geometry_tol must be one positive finite value.")
  }
  W1 <- rkhs_extract_working_model(fit1)
  W2 <- rkhs_extract_working_model(fit2)
  L1 <- .gam_training_lpmatrix(fit1)
  L2 <- .gam_training_lpmatrix(fit2)
  S1 <- .gam_single_smooth(fit1, smooth_index1, L1)
  S2 <- .gam_single_smooth(fit2, smooth_index2, L2)

  if (nrow(S1$B) != nrow(S2$B) || ncol(S1$B) != ncol(S2$B)) {
    stop("The two fits do not share the same observation and smooth dimensions.")
  }
  B_scale <- max(1, max(abs(S1$B)), max(abs(S2$B)))
  if (max(abs(S1$B - S2$B)) > geometry_tol * B_scale) {
    stop("The two fits do not share the same smooth basis and row ordering.")
  }
  Q1 <- as.matrix(S1$Q)
  Q2 <- as.matrix(S2$Q)
  if (!all(dim(Q1) == dim(Q2))) stop("The two smooth penalties have different dimensions.")
  Q_scale <- max(1, max(abs(Q1)), max(abs(Q2)))
  if (max(abs(Q1 - Q2)) > geometry_tol * Q_scale) {
    stop("The two fits do not share the same fixed-kappa RKHS precision Q.")
  }
  if (!identical(S1$score_precision_psd, S2$score_precision_psd)) {
    stop("The two fits do not use the same score-precision type.")
  }

  smooth_cols1 <- unique(unlist(lapply(
    fit1$smooth, function(s) seq.int(s$first.para, s$last.para)
  )))
  smooth_cols2 <- unique(unlist(lapply(
    fit2$smooth, function(s) seq.int(s$first.para, s$last.para)
  )))
  X1 <- L1[, setdiff(seq_len(ncol(L1)), smooth_cols1), drop = FALSE]
  X2 <- L2[, setdiff(seq_len(ncol(L2)), smooth_cols2), drop = FALSE]

  op1 <- rkhs_score_operator(
    S1$B, S1$Q, W1$working_variance, X1,
    field_scale = W1$dispersion / S1$sp,
    allow_psd = S1$score_precision_psd
  )
  op2 <- rkhs_score_operator(
    S2$B, S2$Q, W2$working_variance, X2,
    field_scale = W2$dispersion / S2$sp,
    allow_psd = S2$score_precision_psd
  )
  ans <- rkhs_covariance_score(
    W1$working_error, W2$working_error, op1, op2,
    method = method
  )
  ans$working_family1 <- W1$family
  ans$working_family2 <- W2$family
  ans$smooth_label1 <- S1$label
  ans$smooth_label2 <- S2$label
  ans$field_scale1 <- W1$dispersion / S1$sp
  ans$field_scale2 <- W2$dispersion / S2$sp
  ans$calibration_scope <- "conditional_irls_working_gaussian"
  ans$working_model1 <- W1
  ans$working_model2 <- W2
  ans
}
