test_that("package-local TAPS validates and reuses its training matrix", {
  f <- st_fixture()
  fit <- mgcv::gam(G = f$G, method = "REML")
  L <- predict(fit, type = "lpmatrix")
  expected <- mgcvST:::taps_score_test(fit)
  trace("predict.gam", where = asNamespace("mgcv"), print = FALSE,
        tracer = quote(stop("unexpected predictor call")))
  on.exit(untrace("predict.gam", where = asNamespace("mgcv")), add = TRUE)
  expect_identical(mgcvST:::taps_score_test(fit, lpmatrix = L), expected)
  expect_error(mgcvST:::taps_score_test(fit, lpmatrix = L[-1, ]), "aligned")
  expect_error(mgcvST:::taps_score_test(fit, lpmatrix = L[, ncol(L):1]), "coefficient order")
  expect_error(mgcvST:::taps_score_test(fit, lpmatrix = L * NA_real_), "finite")
})

test_that("local Davies keeps p-values in (0, 1] and failures use the saddlepoint", {
  f <- st_fixture()
  fit <- mgcv::gam(G = f$G, method = "REML")
  reference <- mgcvST:::taps_score_test(fit)
  spec <- attr(reference, "marginal_spectrum")
  sp <- mgcvST:::.mgcvst_marginal_saddlepoint(spec$statistic, spec$lambda)
  expect_true(is.finite(sp))
  expect_saddlepoint <- function(out) {
    expect_identical(out$smooth.pvalue, sp)
    expect_identical(out$method, "saddlepoint")
    expect_identical(attr(out, "marginal_spectrum"), spec)
  }
  for (result in list(list(ifault = 0L, Qq = NA_real_), list(ifault = 0L, Qq = 0),
                      list(ifault = 0L, Qq = 1.1), list(ifault = 0L, Qq = numeric()),
                      list(ifault = 1L, Qq = -1e-14))) {
    local({
      local_mocked_bindings(davies = function(...) result, .package = "CompQuadForm")
      expect_saddlepoint(mgcvST:::taps_score_test(fit))
    })
  }
  for (result in list(list(ifault = 0L, Qq = .321), list(ifault = 1L, Qq = .4),
                      list(Qq = .4))) {
    local({
      local_mocked_bindings(davies = function(...) result, .package = "CompQuadForm")
      out <- mgcvST:::taps_score_test(fit)
      expect_identical(out$smooth.pvalue, result$Qq)
      expect_identical(out$method, "davies")
      expect_identical(attr(out, "marginal_spectrum"), spec)
    })
  }
  local({
    local_mocked_bindings(davies = function(...) stop("integration failed"),
                          .package = "CompQuadForm")
    expect_saddlepoint(mgcvST:::taps_score_test(fit))
  })
  expect_error(mgcvST:::taps_score_test(fit, method = "liu"), "should be")
})
