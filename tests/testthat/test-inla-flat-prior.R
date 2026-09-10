.inlast_flat_fixture <- function(n = 56L, seed = 1201L) {
  skip_if_not_installed("INLA")
  skip_if_not_installed("geometry")
  set.seed(seed)
  vertices <- as.matrix(expand.grid(
    x = seq(0, 1, length.out = 4L),
    y = seq(0, 1, length.out = 4L)
  ))
  mesh <- list(loc = vertices, graph = list(tv = geometry::delaunayn(vertices)))
  data <- data.frame(
    x = runif(n, 0.02, 0.98), y = runif(n, 0.02, 0.98),
    z = seq(-1, 1, length.out = n), exposure = runif(n, 0.8, 1.4)
  )
  data$offset0 <- log(data$exposure)
  basis <- spde_basis(
    mesh, as.matrix(data[c("x", "y")]), kappa = 0.7,
    project_intercept = TRUE
  )
  list(data = data, basis = basis)
}

test_that("flat log-precision prior supports custom starts and retains mean zero", {
  skip_on_cran()
  f <- .inlast_flat_fixture()
  d <- f$data
  set.seed(1202L)
  y <- 0.3 + 0.25 * d$z + d$offset0 +
    0.18 * sin(2 * pi * d$x) + rnorm(nrow(d), sd = 0.35)
  model <- inlaST.set(
    response ~ z + offset(offset0), d, f$basis, family = gaussian()
  )
  engine <- mgcvST:::.inlast_fit_feature(
    model$inla_spec, y, offset = model$offset,
    control = list(
      precision_prior = list(
        prior = "flat", param = numeric(), initial = -0.5
      ),
      gaussian_precision = 1 / 0.35^2,
      keep_fit = TRUE
    )
  )

  expect_true(engine$converged)
  expect_true(all(is.finite(engine$tau_internal)))
  expect_lt(max(abs(engine$constraint_residual), na.rm = TRUE), 1e-12)
  expect_lt(max(abs(engine$observation_spatial_mean), na.rm = TRUE), 1e-12)
  expect_match(engine$estimation$prior_semantics, "improper flat")
  expect_true(engine$estimation$hyperpriors$any_improper)
  expect_identical(
    engine$estimation$hyperpriors$latent_precision$type,
    "improper_flat_log_hyperparameter"
  )
  expect_false(engine$estimation$hyperpriors$latent_precision$proper)
  expect_identical(engine$estimation$hyperpriors$observation$type, "fixed")
  expect_true(engine$estimation$hyper_mode_diagnostics$finite)
  expect_identical(
    engine$estimation$hyper_mode_diagnostics$optimizer_status, 0L
  )
  expect_match(
    engine$estimation$hyper_mode_diagnostics$boundary_check,
    "multi-start sensitivity"
  )
})

test_that("flat prior controls reject parameters and non-finite starts", {
  fields <- c("precision_prior", "gaussian_precision_prior", "nb_size_prior")
  for (field in fields) {
    bad_param <- list()
    bad_param[[field]] <- list(prior = "flat", param = 0, initial = 0)
    expect_error(mgcvST:::.inlast_control(bad_param), "flat")

    bad_initial <- list()
    bad_initial[[field]] <- list(
      prior = "flat", param = numeric(), initial = Inf
    )
    expect_error(mgcvST:::.inlast_control(bad_initial), "finite")
  }
})
