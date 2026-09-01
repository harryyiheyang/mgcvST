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
