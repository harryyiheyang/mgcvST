test_that("public Liu approximation preserves the score and adjustment contract", {
  f <- st_fixture(n = 60L)
  pairs <- rbind(c(1L, 2L), c(1L, 3L), c(2L, 3L))
  for (model in list(f$model, f$G)) {
    fit <- mgcvST.estimate(f$Y, model, diagnostics = FALSE)
    exact <- mgcvST.test(fit, pairs = pairs, method = "BY", threads = 2L)
    approx <- mgcvST.test(fit, pairs = pairs, method = "BY", threads = 2L,
                         approximate = TRUE, n_ref = 3L, ref_tol = 1e-12)
    expect_equal(approx$results$signed_score, exact$results$signed_score,
                 tolerance = 1e-10)
    expect_equal(approx$results$p_two_sided, exact$results$p_two_sided,
                 tolerance = 1e-6)
    expect_equal(approx$results$p_adjusted,
                 p.adjust(approx$results$p_two_sided, "BY"))
    expect_true(approx$timing$pair_pipeline$approximate)
    expect_null(approx$internal)
    expect_error(mgcvST.test(fit, pairs = pairs, approximate = TRUE,
                            calibration = "davies"), "requires calibration")
  }
})
