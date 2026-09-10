test_that("INLA control merges prior objects atomically and preserves NULL", {
  base <- list(
    precision_prior = list(
      prior = "flat", param = numeric(), initial = -3
    ),
    nb_size = 7,
    verbose = FALSE,
    control.inla = list(tolerance = 1e-3, h = 0.004)
  )
  override <- list(
    precision_prior = list(prior = "normal", param = c(-1.25, 0.4)),
    nb_size = NULL,
    control.inla = list(tolerance = 2e-4)
  )
  merged <- mgcvST:::.inlast_merge_control(base, override)
  expect_identical(
    merged$precision_prior,
    list(prior = "normal", param = c(-1.25, 0.4))
  )
  expect_true("nb_size" %in% names(merged))
  expect_null(merged$nb_size)
  expect_identical(
    merged$control.inla,
    list(tolerance = 2e-4, h = 0.004)
  )

  checked <- mgcvST:::.inlast_control(merged)
  expect_identical(
    checked$precision_prior,
    list(prior = "normal", param = c(-1.25, 0.4), initial = 0)
  )
  expect_error(
    mgcvST:::.inlast_merge_control(
      list(), structure(list(FALSE, TRUE), names = c("verbose", "verbose"))
    ),
    "unique"
  )
  expect_error(mgcvST:::.inlast_merge_control(list(), list(FALSE)), "names")
})

test_that("native INLA tuning is validated and EB Gaussian strategy is forced", {
  control <- mgcvST:::.inlast_control(list(
    control.inla = list(
      tolerance = 1e-4, tolerance.f = 5e-5,
      strategy = "gaussian", int.strategy = "eb"
    )
  ))
  expect_identical(
    control$control.inla,
    list(
      tolerance = 1e-4, tolerance.f = 5e-5,
      strategy = "gaussian", int.strategy = "eb"
    )
  )
  expect_error(mgcvST:::.inlast_control(list(
    control.inla = list(strategy = "simplified.laplace")
  )), "gaussian")
  expect_error(mgcvST:::.inlast_control(list(
    control.inla = list(int.strategy = "grid")
  )), "eb")
  expect_error(mgcvST:::.inlast_control(list(
    control.inla = list(reordering = "default")
  )), "Unsupported")
  expect_error(mgcvST:::.inlast_control(list(
    control.inla = list(tolerance = -1)
  )), "must be positive")
  expect_error(mgcvST:::.inlast_control(list(
    control.inla = structure(list(1, 2), names = c("h", "h"))
  )), "unique")
})

test_that("registered and custom INLA prior specifications are validated", {
  checked <- mgcvST:::.inlast_control(list(
    precision_prior = list(
      prior = "normal", param = c(-1.5, 0.25), initial = 1.2
    ),
    gaussian_precision_prior = list(
      prior = "loggamma", param = c(2, 0.4)
    ),
    nb_size_prior = list(prior = "flat")
  ))
  expect_identical(
    checked$precision_prior,
    list(prior = "normal", param = c(-1.5, 0.25), initial = 1.2)
  )
  expect_identical(
    checked$gaussian_precision_prior,
    list(prior = "loggamma", param = c(2, 0.4), initial = 0)
  )
  expect_identical(
    checked$nb_size_prior,
    list(prior = "flat", param = numeric(), initial = 0)
  )

  expression_prior <- paste(
    "expression:",
    "mu = 0.3;",
    "prec = 0.5;",
    "return(0.5*log(prec/(2*pi))-0.5*prec*(log_precision-mu)^2);"
  )
  expression_control <- mgcvST:::.inlast_control(list(
    precision_prior = list(prior = expression_prior, initial = -0.2)
  ))
  expect_identical(
    expression_control$precision_prior,
    list(prior = expression_prior, param = numeric(), initial = -0.2)
  )

  expect_error(mgcvST:::.inlast_control(list(
    precision_prior = list(prior = "normal")
  )), "explicit param")
  expect_error(mgcvST:::.inlast_control(list(
    precision_prior = list(prior = "normal", param = c(0, 0))
  )), "positive precision")
  expect_error(mgcvST:::.inlast_control(list(
    precision_prior = list(prior = "loggamma", param = c(1, -0.1))
  )), "positive shape and rate")
  expect_error(mgcvST:::.inlast_control(list(
    precision_prior = list(prior = "not-an-INLA-prior", param = numeric())
  )), "Unknown INLA")
  expect_error(mgcvST:::.inlast_control(list(
    precision_prior = list(prior = "flat", fixed = TRUE)
  )), "fixed_precision")
})

test_that("prior metadata records the actual family, parameters, and starts", {
  control <- mgcvST:::.inlast_control(list(
    precision_prior = list(
      prior = "loggamma", param = c(1.7, 0.08), initial = -0.7
    ),
    gaussian_precision_prior = list(
      prior = "gaussian", param = c(1.2, 0.25), initial = 0.4
    )
  ))
  metadata <- mgcvST:::.inlast_prior_metadata(
    control, family = "gaussian", fixed_precision = NULL
  )
  expect_identical(metadata$latent_precision$type, "INLA_hyperprior")
  expect_identical(metadata$latent_precision$prior, "loggamma")
  expect_equal(metadata$latent_precision$param, c(1.7, 0.08))
  expect_identical(metadata$latent_precision$initial, -0.7)
  expect_true(metadata$latent_precision$proper)
  expect_identical(metadata$observation$prior, "gaussian")
  expect_equal(metadata$observation$param, c(1.2, 0.25))
  expect_identical(metadata$observation$initial, 0.4)
  expect_match(metadata$observation$statement, "N\\(1.2, 2\\^2\\)")
  expect_false(metadata$any_improper)

  nb_metadata <- mgcvST:::.inlast_prior_metadata(
    mgcvST:::.inlast_control(), family = "negative_binomial",
    fixed_precision = NULL
  )
  expect_true(nb_metadata$any_improper)
  expect_identical(
    nb_metadata$latent_precision$type,
    "improper_flat_log_hyperparameter"
  )
  expect_identical(
    nb_metadata$observation$type,
    "improper_flat_log_hyperparameter"
  )
})

.inlast_custom_prior_fixture <- function(n = 42L, seed = 1811L) {
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
    z = seq(-1, 1, length.out = n), offset0 = runif(n, -0.1, 0.1)
  )
  basis <- spde_basis(
    mesh, as.matrix(data[c("x", "y")]), kappa = 1.1,
    project_intercept = TRUE
  )
  signal <- 0.3 + 0.2 * data$z + data$offset0 +
    0.12 * sin(2 * pi * data$x)
  list(data = data, basis = basis,
       y = signal + rnorm(n, sd = 0.35))
}

test_that("INLA accepts validated loggamma and expression precision priors", {
  skip_on_cran()
  f <- .inlast_custom_prior_fixture()
  model <- inlaST.set(
    response ~ z + offset(offset0), f$data, f$basis, family = gaussian()
  )
  expression_prior <- paste(
    "expression:",
    "mu = 0.3;",
    "prec = 0.5;",
    "return(0.5*log(prec/(2*pi))-0.5*prec*(log_precision-mu)^2);"
  )
  controls <- list(
    loggamma = list(
      prior = "loggamma", param = c(1, 0.1), initial = 0
    ),
    expression = list(prior = expression_prior, initial = 0)
  )
  for (name in names(controls)) {
    engine <- mgcvST:::.inlast_fit_feature(
      model$inla_spec, f$y, offset = model$offset,
      control = list(
        precision_prior = controls[[name]],
        gaussian_precision = 1 / 0.35^2,
        control.inla = list(tolerance = 1e-4)
      )
    )
    expect_true(engine$converged, info = name)
    expect_true(is.finite(engine$tau[["global"]]), info = name)
    expect_true(
      abs(engine$constraint_residual[["global"]]) < 1e-12,
      info = name
    )
    expect_identical(
      engine$estimation$hyperpriors$latent_precision$prior,
      mgcvST:::.inlast_control(list(
        precision_prior = controls[[name]]
      ))$precision_prior$prior
    )
  }
})
