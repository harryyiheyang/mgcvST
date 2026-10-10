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
    z <- mgcvST:::.inlast_spatial_fallback(.spatial_crash(), null, spec, y, offset)
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
    expect_true(is.na(z$mode_status))
  }
})

test_that("fallback leaves input errors and unsuccessful null fits unchanged", {
  e <- simpleError("Invalid input")
  expect_identical(mgcvST:::.inlast_spatial_fallback(e, list(converged = TRUE),
    NULL, NULL, NULL), e)
  e <- .spatial_crash()
  for (null in list(simpleError("null failed"), list(converged = FALSE))) {
    expect_identical(mgcvST:::.inlast_spatial_fallback(e, null, NULL, NULL, NULL), e)
  }
})

test_that("public native crash recovery reuses null fits without replacing p-values", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  skip_if_not_installed("geometry")
  set.seed(9024L)
  xy <- as.matrix(expand.grid(x = seq(0, 1, length.out = 4), y = seq(0, 1, length.out = 4)))
  mesh <- list(loc = xy, graph = list(tv = geometry::delaunayn(xy)))
  d <- data.frame(x = runif(90, .03, .97), y = runif(90, .03, .97),
    cov = seq(-1, 1, length.out = 90), group = factor(rep(1:3, 30)),
    offset0 = seq(-.2, .2, length.out = 90))
  model <- inlaST.set(response ~ cov + offset(offset0) + s(group, bs = "re"),
    d, family = mgcv::nb(theta = 2), mesh = mesh, kappa = .7,
    coordinates = c("x", "y"), precision_scale = "raw",
    control = list(fixed_precision = c(2, 3), poisson_screen_phi = 0))
  Y <- rbind(a = rnbinom(90, mu = exp(1 + d$cov + d$offset0), size = 2),
             b = rnbinom(90, mu = exp(1 - d$cov + d$offset0), size = 2))
  extra <- rbind(seq(-.1, .1, length.out = 90), seq(.1, -.1, length.out = 90))
  calls <- new.env(parent = emptyenv())
  calls$null <- list()
  calls$spatial <- 0L
  native <- mgcvST:::.inlast_fit_feature
  testthat::local_mocked_bindings(.inlast_fit_feature = function(spec, y, offset, control, diagnostics) {
    if (any(vapply(spec$random, function(z) isTRUE(z$target), logical(1L)))) {
      calls$spatial <- calls$spatial + 1L
      stop(.spatial_crash())
    }
    z <- native(spec, y, offset, control, diagnostics)
    calls$null[[length(calls$null) + 1L]] <- z
    z
  }, .package = "mgcvST")
  fit <- inlaST.estimate(Y, model, offset = extra, retain_smooth = TRUE,
    diagnostics = TRUE, BPPARAM = BiocParallel::SerialParam(), spatial = "all")
  expect_length(calls$null, 2L)
  expect_identical(calls$spatial, 2L)
  expect_true(all(fit$diagnostics$converged))
  expect_true(all(fit$diagnostics$spatial_fallback))
  expect_identical(fit$diagnostics$outer_convergence, rep("null_zero_spatial_fallback", 2))
  expect_identical(fit$diagnostics$error_class, rep("inlaCrashError", 2))
  expect_true(all(is.na(fit$diagnostics$criterion)))
  expect_equal(fit$target_coefficients, fit$target_coefficients * 0, tolerance = 0)
  expect_true(all(fit$smooth_coefficients$global == 0))
  expect_equal(unname(fit$lambda), rep(1e8, 2), tolerance = 0)
  expect_equal(fit$offset, sweep(extra, 2, model$offset, "+"), tolerance = 0)
  expect_true(all(is.finite(fit$score_a)))
  expect_true(all(is.finite(fit$diagnostics$marginal_p_value)))
  S <- mgcvST:::.inlast_null_spec(model$inla_spec)
  expected <- mgcvST:::.inlast_null_marginal(fit$feature_id, fit$score_sparse,
    model$geometry$nuisance_design, calls$null, S,
    vapply(calls$null, `[[`, numeric(1), "dispersion"),
    do.call(rbind, lapply(calls$null, `[[`, "smoothing_parameters")),
    1:2, chunk_size = 2L, threads = 1L)
  expect_equal(fit$diagnostics$marginal_p_value, expected$p_value, tolerance = 0)
  pair <- inlaST.test(fit, adjust = "none", rank = 2L,
    pairs = matrix(c("a", "b"), ncol = 2L), threads = 1L, moments = "exact")
  expect_true(all(is.finite(pair$results$log_p_two_sided)))
  expect_true(all(is.finite(pair$results$score)))
  for (j in 1:2) {
    expected_nuisance <- c(calls$null[[j]]$fixed_mode, unlist(calls$null[[j]]$random_mode))
    expect_equal(unname(fit$nuisance_coefficients[, j]), unname(expected_nuisance), tolerance = 0)
    expect_equal(fit$smoothing_parameters[j, 2], calls$null[[j]]$smoothing_parameters[2], tolerance = 0)
    expect_equal(unname(fit$inla_diagnostics[[j]]$tau[1]), 1e8, tolerance = 0)
  }
  expect_error(inlaST.estimate(Y + .1, model, spatial = "all"), "non-negative integers")
  expect_identical(calls$spatial, 2L)
})
