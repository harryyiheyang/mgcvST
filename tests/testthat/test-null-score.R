test_that("null-first BAM score retains a calibration cache", {
  f <- st_fixture(nuisance = TRUE)
  fit <- mgcvST.estimate(
    f$Y[1L, , drop = FALSE], f$model, retain_marginal = TRUE,
    BPPARAM = BiocParallel::SerialParam()
  )
  expect_true(is.finite(fit$diagnostics$marginal_p_value[1L]))
  expect_identical(fit$marginal_data$version, 2L)
  out <- mgcvST.marginal(fit)
  expect_true(is.finite(out$p_value[1L]))
})

test_that("external GAM setup uses its fixed null design", {
  f <- st_fixture()
  fit <- mgcvST.estimate(
    f$Y[1L, , drop = FALSE], f$G,
    BPPARAM = BiocParallel::SerialParam()
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
    control = control,
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
  expected <- mgcvST:::.mgcvst_null_score_test(null_fit, setup)
  expect_equal(fit$diagnostics$marginal_p_value, expected$smooth.pvalue,
               tolerance = 1e-8)
})

test_that("single-coefficient null scores agree with direct GAM fits", {
  for (family in list(stats::gaussian(), stats::quasipoisson(), mgcv::nb())) {
    f <- st_fixture(n = 60L, family = family)
    basis <- f$basis
    model <- mgcvST.set(
      response ~ offset(offset0) + s(x, y, bs = "spde", xt = basis),
      f$data, family = family
    )
    control <- mgcv::gam.control(nthreads = 1L)
    control$poisson_screen_phi <- 0
    data <- model$null_data
    data[[model$null_response]] <- as.numeric(f$Y[1L, ])
    direct <- mgcv::gam(model$null_formula, data = data,
                        family = unserialize(serialize(family, NULL)),
                        method = "REML", control = mgcv::gam.control(nthreads = 1L))
    setup <- mgcvST:::.mgcvst_null_score_setup(model$G, 1L, list(
      formula = model$null_formula, data = model$null_data,
      response = model$null_response, X0 = model$null_X
    ))
    expected <- mgcvST:::.mgcvst_null_score_test(direct, setup)
    for (G in list(model, model$G)) {
      fit <- mgcvST.estimate(
        f$Y[1L, , drop = FALSE], G, control = control, retain_marginal = TRUE,
        BPPARAM = BiocParallel::SerialParam()
      )
      expect_true(is.finite(fit$diagnostics$marginal_p_value))
      expect_true(all(is.na(fit$diagnostics$error_message)))
      expect_equal(fit$diagnostics$marginal_p_value, expected$smooth.pvalue,
                   tolerance = 1e-8)
      retained <- mgcvST.marginal(fit)
      expect_equal(retained$p_value, expected$smooth.pvalue, tolerance = 1e-8)
    }
  }
})

test_that("one slope without an intercept retains per-feature null offsets", {
  f <- st_fixture(n = 60L)
  basis <- f$basis
  model <- mgcvST.set(
    response ~ 0 + z + offset(offset0) + s(x, y, bs = "spde", xt = basis),
    f$data, family = mgcv::nb()
  )
  Y <- f$Y[1:2, , drop = FALSE]
  offset <- matrix(seq(-0.2, 0.2, length.out = length(Y)), nrow(Y))
  control <- mgcv::gam.control(nthreads = 1L)
  control$poisson_screen_phi <- 0
  fit <- mgcvST.estimate(
    Y, model, offset = offset, control = control,
    BPPARAM = BiocParallel::SerialParam()
  )
  expect_true(all(is.finite(fit$diagnostics$marginal_p_value)))
  setup <- mgcvST:::.mgcvst_null_score_setup(model$G, 1L, list(
    formula = model$null_formula, data = model$null_data,
    response = model$null_response, X0 = model$null_X
  ))
  for (j in seq_len(nrow(Y))) {
    data <- model$null_data
    data[[model$null_response]] <- as.numeric(Y[j, ])
    direct <- mgcv::gam(model$null_formula, data = data, family = mgcv::nb(),
                        offset = offset[j, ], method = "REML",
                        control = mgcv::gam.control(nthreads = 1L))
    expected <- mgcvST:::.mgcvst_null_score_test(direct, setup)
    expect_equal(fit$diagnostics$marginal_p_value[j], expected$smooth.pvalue,
                 tolerance = 1e-8)
  }
})

test_that("single-coefficient null fitting is available on SOCK workers", {
  skip_on_cran()
  f <- st_fixture(n = 60L)
  basis <- f$basis
  model <- mgcvST.set(
    response ~ offset(offset0) + s(x, y, bs = "spde", xt = basis),
    f$data, family = mgcv::nb()
  )
  Y <- f$Y[1:2, , drop = FALSE]
  serial <- mgcvST.estimate(Y, model$G,
                            BPPARAM = BiocParallel::SerialParam(), chunk_size = 1L)
  workers <- BiocParallel::SnowParam(2L, type = "SOCK", progressbar = FALSE)
  on.exit(BiocParallel::bpstop(workers), add = TRUE)
  parallel <- mgcvST.estimate(Y, model$G,
                              BPPARAM = workers, chunk_size = 1L)
  expect_true(all(is.finite(parallel$diagnostics$marginal_p_value)))
  expect_true(all(is.na(parallel$diagnostics$error_message)))
  expect_equal(parallel$diagnostics$marginal_p_value,
               serial$diagnostics$marginal_p_value, tolerance = 1e-10)
})

test_that("parametric null fitting preserves weights and rejects unsupported AR", {
  d <- data.frame(response = rep(c(1, 2, 4, 3, 7), 8),
                   off = seq(-0.2, 0.2, length.out = 40))
  setup <- list(spec = list(formula = response ~ offset(off)), single_parametric = TRUE)
  offset <- seq(0.1, -0.1, length.out = nrow(d))
  weights <- rep(c(1, 2), length.out = nrow(d))
  control <- mgcv::gam.control(nthreads = 1L)
  for (family in list(stats::gaussian(), stats::quasipoisson(), mgcv::nb())) {
    family_raw <- serialize(family, NULL)
    fit <- mgcvST:::.mgcvst_fit_null(setup, d, unserialize(family_raw), offset, control,
      list(weights = weights, gamma = 1.2, chunk.size = 10L))
    expected <- mgcv::gam(response ~ offset(off), data = d,
      family = unserialize(family_raw), offset = offset,
      weights = weights, gamma = 1.2, control = control, method = "REML")
    expect_equal(stats::coef(fit), stats::coef(expected), tolerance = 1e-12)
    expect_equal(fit$Vp, expected$Vp, tolerance = 1e-12)
    expect_equal(fit$offset, expected$offset, tolerance = 1e-12)
    expect_equal(fit$prior.weights, weights)
  }
  expect_error(mgcvST:::.mgcvst_fit_null(setup, d, stats::gaussian(), offset,
    control, list(rho = 0.5)), "cannot preserve BAM's nonzero rho")
})
