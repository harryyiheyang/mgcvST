.spatial_crash <- function() {
  structure(list(message = "The inla program call crashed.", call = quote(INLA::inla())),
    class = c("inlaCrashError", "error", "condition"))
}

test_that("native crash fallback preserves null nuisance states and raw precision units", {
  n <- 8L
  X <- cbind(1, seq(-1, 1, length.out = n))
  A <- Matrix::Matrix(diag(n), sparse = TRUE)
  Z <- Matrix::Matrix(model.matrix(~ factor(rep(1:2, 4)) - 1), sparse = TRUE)
  spec <- list(fixed = list(X = X), geometry_sp_length = 2L,
    random = list(
      list(name = "global", target = TRUE, A = A, precision_scale = 4,
        projection = matrix(0, n, n - 1L), sp_index = 1L),
      list(name = "group", target = FALSE, A = Z, sp_index = 2L)))
  for (family in c("gaussian", "poisson", "negative_binomial")) {
    null <- list(converged = TRUE, dispersion = if (family == "gaussian") .3 else 1,
      family = family, family_parameters = if (family == "negative_binomial") 2 else numeric(),
      fixed_mode = c(.2, .4), fixed_mean = c(.2, .4), random_mode = list(group = c(-.1, .1)),
      tau = c(group = 3), tau_internal = c(group = 3), precision_scale = c(group = 1),
      lambda = c(group = if (family == "gaussian") .9 else 3),
      constraint_residual = c(group = NA_real_),
      constraint_residual_uncorrected = c(group = NA_real_),
      observation_spatial_mean = c(group = NA_real_), estimation = list())
    offset <- seq(-.2, .2, length.out = n)
    y <- 1:8
    z <- spdeST:::.inlast_spatial_fallback(.spatial_crash(), null, spec, y, offset)
    expect_identical(z$fixed_mode, null$fixed_mode)
    expect_identical(z$random_mode$group, null$random_mode$group)
    expect_identical(z$random_mode$global, numeric(n))
    expect_identical(z$coefficients$global, numeric(n - 1L))
    expect_identical(z$dispersion, null$dispersion)
    expect_identical(z$family_parameters, null$family_parameters)
    expect_equal(z$tau, c(global = 1e8, group = 3), tolerance = 0)
    expect_equal(z$tau_internal, c(global = 2.5e7, group = 3), tolerance = 0)
    expect_equal(z$lambda, z$dispersion * z$tau, tolerance = 1e-15)
    expect_equal(z$smoothing_parameters, unname(z$lambda), tolerance = 0)
    eta <- offset + as.numeric(X %*% null$fixed_mode + Z %*% null$random_mode$group)
    expect_equal(z$eta, eta, tolerance = 1e-15)
    mu <- if (family == "gaussian") eta else exp(eta)
    expect_equal(z$mu, mu, tolerance = 1e-15)
    V <- if (family == "gaussian") rep(.3, n) else 1 / mu + if (family == "negative_binomial") .5 else 0
    expect_equal(z$working_variance, V, tolerance = 1e-15)
    E <- if (family == "gaussian") y - offset else eta + (y - mu) / mu - offset
    expect_equal(z$working_error, E, tolerance = 1e-15)
    expect_identical(z$spatial_fallback$error_message, conditionMessage(.spatial_crash()))
    expect_identical(z$spatial_fallback$native_converged, FALSE)
    expect_true(is.na(z$log_marginal_likelihood))
  }
})

test_that("fallback leaves input errors and unsuccessful null fits unchanged", {
  e <- simpleError("Invalid input")
  expect_identical(spdeST:::.inlast_spatial_fallback(e, list(converged = TRUE),
    NULL, NULL, NULL), e)
  e <- .spatial_crash()
  for (null in list(simpleError("null failed"), list(converged = FALSE))) {
    expect_identical(spdeST:::.inlast_spatial_fallback(e, null, NULL, NULL, NULL), e)
  }
})
