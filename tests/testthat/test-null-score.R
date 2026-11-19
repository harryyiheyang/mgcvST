test_that("external GAM setup uses its fixed null design", {
  f <- st_fixture()
  fit <- spdeST:::mgcvST.estimate_spatial(
    mgcvST.estimate(
    f$Y[1L, , drop = FALSE], mgcvST.set(G = f$G),
    BPPARAM = BiocParallel::SerialParam()),
    f$Y[1L, , drop = FALSE], features = "all", BPPARAM = BiocParallel::SerialParam())
  expect_true(is.finite(fit$diagnostics$marginal_p_value[1L]))
  expect_null(fit$gam)
})

test_that("Poisson prescreen routes independent null and full BAM fits", {
  control <- mgcv::gam.control(nthreads = 1L)
  local_poisson_route(TRUE)
  f <- st_fixture(n = 60L)
  for (setup in list(mgcvST.set(G = f$G), f$model)) {
    fit <- spdeST:::mgcvST.estimate_spatial(
    mgcvST.estimate(
      f$Y[1L, , drop = FALSE], setup, control = control,
      BPPARAM = BiocParallel::SerialParam()),
    f$Y[1L, , drop = FALSE], features = "all", BPPARAM = BiocParallel::SerialParam())
    expect_identical(fit$diagnostics$family_used, "poisson")
    expect_true(is.finite(fit$diagnostics$marginal_p_value[1L]))
    expect_true(is.finite(fit$dispersion[1L]))
  }
})

