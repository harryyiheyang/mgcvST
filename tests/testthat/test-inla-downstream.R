test_that("INLA-named downstream wrappers preserve delegated validation", {
  bad <- structure(list(), class = "inlaST_fit")
  expect_error(inlaST.test(bad, pairs = matrix(c(1L, 2L), nrow = 1L)),
               "score-test engine|feature geometry")
  expect_error(inlaST.marginal(bad, calibration = "liu"),
               "retain_marginal")
})

test_that("INLA-named downstream wrappers are exact delegates", {
  f <- st_fixture(n = 45L, family = gaussian())
  fit <- mgcvST.estimate(
    f$Y[1:2, ], f$model, retain_marginal = TRUE,
    marginal_args = list(method = "liu"),
    BPPARAM = BiocParallel::SerialParam()
  )
  class(fit) <- c("inlaST_fit", class(fit))
  pair <- matrix(c(1L, 2L), nrow = 1L)
  direct_test <- mgcvST.test(fit, pairs = pair, calibration = "liu",
                             BPPARAM = BiocParallel::SerialParam())
  named_test <- inlaST.test(fit, pairs = pair, calibration = "liu",
                            BPPARAM = BiocParallel::SerialParam())
  expect_numerically_equivalent_test(named_test, direct_test)

  direct_marginal <- mgcvST.marginal(
    fit, features = 2:1, calibration = "liu",
    BPPARAM = BiocParallel::SerialParam()
  )
  named_marginal <- inlaST.marginal(
    fit, features = 2:1, calibration = "liu",
    BPPARAM = BiocParallel::SerialParam()
  )
  expect_identical(named_marginal, direct_marginal)
})
