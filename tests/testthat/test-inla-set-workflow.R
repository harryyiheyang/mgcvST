.inlast_set_workflow_fixture <- function(n = 42L, seed = 1601L) {
  skip_if_not_installed("INLA")
  skip_if_not_installed("geometry")
  set.seed(seed)
  vertices <- as.matrix(expand.grid(
    u = seq(0, 1, length.out = 4L),
    v = seq(0, 1, length.out = 4L)
  ))
  mesh <- list(loc = vertices, graph = list(tv = geometry::delaunayn(vertices)))
  data <- data.frame(
    u = runif(n, .03, .97), v = runif(n, .03, .97),
    z = seq(-1, 1, length.out = n),
    offset0 = seq(-.12, .12, length.out = n)
  )
  basis <- spde_basis(
    mesh, as.matrix(data[c("u", "v")]), kappa = .7,
    project_intercept = TRUE
  )
  basis$component <- basis$score.component <- "global"
  eta <- .4 + .25 * data$z + data$offset0 +
    .2 * sin(2 * pi * data$u) - .15 * cos(2 * pi * data$v)
  y <- eta + rnorm(n, sd = .35)
  Y <- rbind(feature_a = y, feature_b = y + rnorm(n, sd = .08))
  list(data = data, basis = basis, mesh = mesh, Y = Y)
}

.inlast_expect_same_set_geometry <- function(actual, expected,
                                              tolerance = 1e-10) {
  expect_identical(actual$setting, expected$setting)
  expect_identical(actual$components, expected$components)
  expect_identical(actual$inla_spec$family, expected$inla_spec$family)
  expect_equal(actual$offset, expected$offset, tolerance = tolerance)
  expect_equal(unname(actual$L), unname(expected$L), tolerance = tolerance)
  expect_equal(actual$geometry$X, expected$geometry$X, tolerance = tolerance)
  expect_identical(names(actual$inla_spec$random),
                   names(expected$inla_spec$random))
  for (j in seq_along(actual$inla_spec$random)) {
    a <- actual$inla_spec$random[[j]]
    b <- expected$inla_spec$random[[j]]
    expect_identical(a$name, b$name)
    expect_identical(a$target, b$target)
    expect_equal(as.matrix(a$A), as.matrix(b$A), tolerance = tolerance)
    expect_equal(as.matrix(a$Q), as.matrix(b$Q), tolerance = tolerance)
    expect_equal(a$constraint, b$constraint, tolerance = tolerance)
    expect_equal(a$projection, b$projection, tolerance = tolerance)
  }
}

.inlast_set_workflow_spatial_mean <- function(fit, feature = 1L) {
  vapply(fit$score_components, function(component) {
    j <- fit$geometry$target[[component]]
    B <- fit$geometry$smooth[[j]]$B
    coefficient <- fit$smooth_coefficients[[component]][feature, ]
    mean(as.numeric(B %*% coefficient))
  }, numeric(1L))
}

test_that("inlaST.set accepts basis, complete-formula and frozen-G workflows", {
  skip_on_cran()
  f <- .inlast_set_workflow_fixture()
  d <- f$data
  basis <- f$basis
  stored <- list(fixed_precision = 1.25, gaussian_precision = 4,
                 fixed_effect_precision = 0)

  by_basis <- inlaST.set(
    response ~ z + offset(offset0), d, basis,
    family = gaussian(), coordinates = c("u", "v"), control = stored
  )
  complete <- response ~ z + offset(offset0) +
    s(u, v, bs = "spde", xt = basis)
  # The third positional family is retained for mgcvST.set-style calls.
  by_formula <- inlaST.set(complete, d, gaussian(), control = stored)
  prepared <- mgcvST.set(complete, d, gaussian())
  by_G <- inlaST.set(G = prepared$G, control = stored)

  for (model in list(by_basis, by_formula, by_G)) {
    expect_s3_class(model, "inlaST_model")
    expect_true(model$mean_constraint_active)
    expect_identical(model$mean_constraint, "observation")
    expect_true(model$inla_spec$mean_constraint_active)
    expect_identical(model$setting, "global")
    expect_identical(model$components, "global")
    expect_identical(model$inla_control$fixed_precision, 1.25)
    expect_identical(model$inla_control$gaussian_precision, 4)
  }
  .inlast_expect_same_set_geometry(by_formula, by_basis)
  .inlast_expect_same_set_geometry(by_G, by_basis)

  local <- spde_basis(
    f$mesh, as.matrix(d[c("u", "v")]), kappa = 5,
    project_intercept = TRUE
  )
  local$component <- local$score.component <- "local"
  two_component <- inlaST.set(
    response ~ z + offset(offset0), d,
    list(global = basis, local = local), family = gaussian(),
    setting = "global_local", coordinates = c("u", "v"), control = stored
  )
  expect_identical(two_component$setting, "global_local")
  expect_identical(two_component$components, c("global", "local"))
  expect_true(all(vapply(two_component$inla_spec$random[1:2], function(block) {
    expected <- as.numeric(crossprod(block$A, rep(1 / nrow(d), nrow(d))))
    isTRUE(block$target) && isTRUE(all.equal(block$constraint, expected,
                                            tolerance = 1e-12)) &&
      max(abs(colMeans(as.matrix(block$A %*% block$projection)))) < 1e-10
  }, logical(1L))))

  d$u_local <- d$u
  d$v_local <- d$v
  complete_two <- response ~ z + offset(offset0) +
    s(u, v, bs = "spde", xt = basis) +
    s(u_local, v_local, bs = "spde", xt = local)
  inferred_two <- inlaST.set(complete_two, d, gaussian(), control = stored)
  prepared_two <- mgcvST.set(complete_two, d, gaussian())
  frozen_two <- inlaST.set(G = prepared_two$G, control = stored)
  expect_identical(inferred_two$setting, "global_local")
  expect_identical(frozen_two$setting, "global_local")
  .inlast_expect_same_set_geometry(inferred_two, two_component)
  .inlast_expect_same_set_geometry(frozen_two, two_component)

  expect_error(
    inlaST.set(complete, d, gaussian(), setting = "global_local",
               control = stored),
    "setting|component|global_local"
  )
  expect_error(
    inlaST.set(G = prepared$G, setting = "global_local", control = stored),
    "setting|component|global_local"
  )

  uncentred <- spde_basis(
    f$mesh, as.matrix(d[c("u", "v")]), kappa = .7,
    project_intercept = FALSE
  )
  bad_formula <- response ~ z + offset(offset0) +
    s(u, v, bs = "spde", xt = uncentred)
  bad_G <- mgcvST.set(bad_formula, d, gaussian())$G
  expect_error(inlaST.set(G = bad_G), "mean|cent|formula|projection")
})

test_that("diagnostic covariance is optional and leaves the native score unchanged", {
  skip_on_cran()
  f <- .inlast_set_workflow_fixture()
  model <- inlaST.set(
    response ~ z + offset(offset0), f$data, f$basis, family = gaussian(),
    coordinates = c("u", "v"),
    control = list(fixed_precision = 1.25, gaussian_precision = 4)
  )
  diagnostic <- inlaST.estimate(f$Y, model, diagnostics = TRUE)
  expect_true(all(vapply(diagnostic$expected_nuisance_covariance, is.matrix, logical(1))))
  testthat::local_mocked_bindings(
    .inlast_expected_covariance = function(...) stop("Unrequested diagnostic solve"),
    .package = "mgcvST"
  )
  ordinary <- inlaST.estimate(f$Y, model)
  expect_true(all(ordinary$diagnostics$converged))
  expect_true(all(vapply(ordinary$expected_nuisance_covariance, is.null, logical(1))))
  for (field in c("nuisance_covariance", "working_error", "working_variance")) {
    expect_equal(ordinary[[field]], diagnostic[[field]], tolerance = 1e-9)
  }
  a <- inlaST.test(ordinary, pairs = matrix(c(1, 2), 1), calibration = "davies")
  b <- inlaST.test(diagnostic, pairs = matrix(c(1, 2), 1), calibration = "davies")
  for (field in c("signed_score", "information", "p_two_sided")) {
    expect_equal(a$results[[field]], b$results[[field]], tolerance = 1e-9)
  }
})

test_that("stored INLA controls are inherited and explicit controls override", {
  skip_on_cran()
  f <- .inlast_set_workflow_fixture(seed = 1602L)
  stored <- list(fixed_precision = 1.25, gaussian_precision = 4,
                 fixed_effect_precision = 0,
                 precision_prior = list(
                   prior = "normal", param = c(0, 1 / 9), initial = -.5
                 ))
  model <- inlaST.set(
    response ~ z + offset(offset0), f$data, f$basis,
    family = gaussian(), coordinates = c("u", "v"), control = stored,
    score_backend = "dense"
  )
  inherited <- inlaST.estimate(
    f$Y, model, retain_smooth = TRUE, retain_marginal = TRUE,
    BPPARAM = BiocParallel::SerialParam(), control = list(),
    marginal_args = list(method = "liu")
  )
  overridden <- inlaST.estimate(
    f$Y[1L, , drop = FALSE], model, retain_smooth = TRUE,
    BPPARAM = BiocParallel::SerialParam(),
    control = list(
      fixed_precision = 2,
      precision_prior = list(prior = "flat", param = numeric(), initial = 1)
    ),
    marginal_args = list(method = "liu")
  )

  expect_equal(unname(inherited$dispersion), rep(.25, 2L), tolerance = 1e-12)
  expect_equal(unname(overridden$dispersion), .25, tolerance = 1e-12)
  expect_equal(unname(inherited$lambda), rep(.25 * 1.25, 2L),
               tolerance = 1e-10)
  expect_equal(unname(overridden$lambda), .25 * 2, tolerance = 1e-10)
  expect_lt(max(abs(inherited$observation_spatial_mean)), 1e-10)
  expect_lt(max(abs(overridden$observation_spatial_mean)), 1e-10)
  expect_lt(max(abs(.inlast_set_workflow_spatial_mean(inherited))), 1e-10)
  expect_true(all(vapply(inherited$nuisance_covariance, is.matrix, logical(1L))))
  expect_identical(inherited$geometry$nuisance_projection,
                   "conditional_INLA_block")
  expect_identical(inherited$score_backend, "dense")
  expect_identical(
    inherited$estimation$control$precision_prior,
    list(prior = "normal", param = c(0, 1 / 9), initial = -.5)
  )
  expect_identical(
    overridden$estimation$control$precision_prior,
    list(prior = "flat", param = numeric(), initial = 1)
  )

  pair <- matrix(c(1L, 2L), nrow = 1L)
  direct_pair <- mgcvST.test(
    inherited, pairs = pair, calibration = "liu",
    BPPARAM = BiocParallel::SerialParam()
  )
  named_pair <- inlaST.test(
    inherited, pairs = pair, calibration = "liu",
    BPPARAM = BiocParallel::SerialParam()
  )
  expect_numerically_equivalent_test(named_pair, direct_pair)
  direct_marginal <- mgcvST.marginal(
    inherited, features = 2:1, calibration = "liu",
    BPPARAM = BiocParallel::SerialParam()
  )
  named_marginal <- inlaST.marginal(
    inherited, features = 2:1, calibration = "liu",
    BPPARAM = BiocParallel::SerialParam()
  )
  expect_identical(named_marginal, direct_marginal)
  expect_error(
    inlaST.set(
      response ~ z, f$data, f$basis, family = gaussian(),
      coordinates = c("u", "v"), control = list(not_a_control = 1)
    ),
    "Unknown|unknown|control"
  )
})

test_that("a fixed NB family does not partially match nb_size_prior", {
  skip_on_cran()
  f <- .inlast_set_workflow_fixture(seed = 1604L)
  prior <- list(prior = "normal", param = c(0, 1 / 9), initial = .2)
  model <- inlaST.set(
    response ~ z + offset(offset0), f$data, f$basis,
    family = mgcv::nb(theta = 4), coordinates = c("u", "v"),
    control = list(nb_size_prior = prior)
  )

  expect_identical(model$inla_control[["nb_size", exact = TRUE]], 4)
  expect_identical(model$inla_control[["nb_size_prior", exact = TRUE]], prior)
})

test_that("set control and native geometry survive serialization and SOCK", {
  skip_on_cran()
  f <- .inlast_set_workflow_fixture(n = 46L, seed = 1603L)
  basis <- f$basis
  complete <- response ~ z + offset(offset0) +
    s(u, v, bs = "spde", xt = basis)
  stored <- list(fixed_precision = 1.4, gaussian_precision = 5,
                 fixed_effect_precision = 0)
  model <- inlaST.set(complete, f$data, gaussian(), control = stored)
  restored <- unserialize(serialize(model, NULL))

  expect_identical(restored$inla_control, model$inla_control)
  expect_identical(restored$score_backend, model$score_backend)
  .inlast_expect_same_set_geometry(restored, model, tolerance = 0)
  serial <- inlaST.estimate(
    f$Y, model, retain_smooth = TRUE,
    BPPARAM = BiocParallel::SerialParam(), control = list(),
    marginal_args = list(method = "liu")
  )
  bp <- BiocParallel::SnowParam(2L, type = "SOCK")
  parallel <- tryCatch(
    inlaST.estimate(
      f$Y, restored, retain_smooth = TRUE, BPPARAM = bp,
      control = list(), chunk_size = 1L,
      marginal_args = list(method = "liu")
    ),
    finally = BiocParallel::bpstop(bp)
  )

  expect_equal(serial$working_error, parallel$working_error, tolerance = 2e-5)
  expect_equal(serial$working_variance, parallel$working_variance,
               tolerance = 2e-5)
  expect_equal(serial$nuisance_covariance, parallel$nuisance_covariance,
               tolerance = 2e-5)
  expect_equal(serial$lambda, parallel$lambda, tolerance = 2e-5)
  expect_lt(max(abs(parallel$observation_spatial_mean)), 1e-10)
  expect_lt(max(abs(.inlast_set_workflow_spatial_mean(parallel, 1L))), 1e-10)
  expect_lt(max(abs(.inlast_set_workflow_spatial_mean(parallel, 2L))), 1e-10)
  expect_identical(parallel$geometry$nuisance_projection,
                   "conditional_INLA_block")
})
