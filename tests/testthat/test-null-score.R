test_that("null-first BAM score retains a calibration cache", {
  f <- st_fixture(nuisance = TRUE)
  fit <- mgcvST.estimate(
    f$Y[1L, , drop = FALSE], f$model, retain_marginal = TRUE,
    marginal_args = list(method = "liu"), BPPARAM = BiocParallel::SerialParam()
  )
  expect_true(is.finite(fit$diagnostics$marginal_p_value[1L]))
  expect_identical(fit$marginal_data$version, 2L)
  out <- mgcvST.marginal(fit, calibration = "liu")
  expect_true(is.finite(out$p_value[1L]))
})

test_that("external GAM setup uses its fixed null design", {
  f <- st_fixture()
  fit <- mgcvST.estimate(
    f$Y[1L, , drop = FALSE], f$G,
    marginal_args = list(method = "liu"), BPPARAM = BiocParallel::SerialParam()
  )
  expect_true(is.finite(fit$diagnostics$marginal_p_value[1L]))
  expect_null(fit$gam)
})

test_that("Poisson prescreen routes independent null and full BAM fits", {
  control <- mgcv::gam.control(nthreads = 1L)
  control$poisson_screen_phi <- 1e6
  f <- st_fixture(n = 60L, nuisance = TRUE)
  for (setup in list(f$G, f$model)) {
    fit <- mgcvST.estimate(
      f$Y[1L, , drop = FALSE], setup, control = control,
      marginal_args = list(method = "liu"),
      BPPARAM = BiocParallel::SerialParam()
    )
    expect_identical(fit$diagnostics$family_used, "quasipoisson")
    expect_true(is.finite(fit$diagnostics$marginal_p_value[1L]))
    expect_true(is.finite(fit$dispersion[1L]))
  }
})

test_that("model-set null BAM includes a feature offset", {
  f <- st_fixture(n = 60L, nuisance = TRUE)
  model <- mgcvST.set(
    response ~ offset(offset0) + s(z, k = 5) +
      s(x, y, bs = "spde", xt = f$basis),
    f$data, family = f$model$G$family
  )
  extra_offset <- seq(-0.2, 0.2, length.out = ncol(f$Y))
  control <- mgcv::gam.control(nthreads = 1L)
  control$poisson_screen_phi <- 0
  fit <- mgcvST.estimate(
    f$Y[1L, , drop = FALSE], model, offset = extra_offset,
    control = control, marginal_args = list(method = "liu"),
    BPPARAM = BiocParallel::SerialParam()
  )
  data <- model$null_data
  data[[model$null_response]] <- as.numeric(f$Y[1L, ])
  null_fit <- mgcv::bam(
    model$null_formula, data = data, family = model$G$family,
    offset = extra_offset, method = "fREML", discrete = TRUE, nthreads = 1L,
    control = mgcv::gam.control(nthreads = 1L)
  )
  G <- structure(model$G, null_spec = list(
    formula = model$null_formula, data = model$null_data,
    response = model$null_response, X0 = model$null_X
  ))
  target <- which(vapply(
    G$smooth, function(s) identical(s$score.component, "global"), logical(1L)
  ))
  setup <- mgcvST:::.mgcvst_null_score_setup(G, target, attr(G, "null_spec"))
  expected <- mgcvST:::.mgcvst_null_score_test(null_fit, setup, method = "liu")
  expect_equal(fit$diagnostics$marginal_p_value, expected$smooth.pvalue,
               tolerance = 1e-8)
})
