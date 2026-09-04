test_that("package-local TAPS validates and reuses its training matrix", {
  f <- st_fixture()
  fit <- mgcv::gam(G = f$G, method = "REML")
  L <- predict(fit, type = "lpmatrix")
  expected <- mgcvST:::taps_score_test(fit, method = "liu")
  trace("predict.gam", where = asNamespace("mgcv"), print = FALSE,
        tracer = quote(stop("unexpected predictor call")))
  on.exit(untrace("predict.gam", where = asNamespace("mgcv")), add = TRUE)
  expect_identical(mgcvST:::taps_score_test(fit, method = "liu", lpmatrix = L), expected)
  expect_error(mgcvST:::taps_score_test(fit, lpmatrix = L[-1, ]), "aligned")
  expect_error(mgcvST:::taps_score_test(fit, lpmatrix = L[, ncol(L):1]), "coefficient order")
  expect_error(mgcvST:::taps_score_test(fit, lpmatrix = L * NA_real_), "finite")
})

test_that("local Davies preserves success and all numerical failures use Liu", {
  f <- st_fixture()
  fit <- mgcv::gam(G = f$G, method = "REML")
  liu <- mgcvST:::taps_score_test(fit, method = "liu")
  for (result in list(list(ifault = 1L, Qq = .4), list(ifault = 0L, Qq = NA_real_),
                      list(ifault = 0L, Qq = 0), list(ifault = 0L, Qq = 1.1),
                      list(ifault = 0L, Qq = numeric()), list(Qq = .4))) {
    local({
      local_mocked_bindings(davies = function(...) result, .package = "CompQuadForm")
      expect_identical(mgcvST:::taps_score_test(fit), liu)
    })
  }
  local({
    local_mocked_bindings(davies = function(...) list(ifault = 0L, Qq = .321),
                          .package = "CompQuadForm")
    out <- mgcvST:::taps_score_test(fit)
    expect_identical(out$smooth.pvalue, .321)
    expect_identical(out$method, "davies")
  })
  local({
    local_mocked_bindings(davies = function(...) stop("integration failed"),
                          .package = "CompQuadForm")
    expect_identical(mgcvST:::taps_score_test(fit), liu)
    expect_identical(mgcvST:::taps_score_test(fit, method = "liu"), liu)
  })
})
