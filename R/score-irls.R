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
