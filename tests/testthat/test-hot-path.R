test_that("legacy switches exactly reproduce existing fitting outputs", {
  old <- old_st()
  f <- st_fixture()
  marginal_callback <- function(...) list(smooth.pvalue = .375)
  for (G in list(f$G, f$model)) {
    a <- old$mgcvST.estimate(f$Y, G, marginal_test = marginal_callback, retain_smooth = TRUE)
    b <- mgcvST.estimate(f$Y, G, marginal = TRUE, diagnostics = TRUE,
                         marginal_test = marginal_callback, retain_smooth = TRUE)
    expect_identical(strip_elapsed(b), strip_elapsed(a))
    fast <- mgcvST.estimate(f$Y, G, marginal = FALSE, diagnostics = FALSE,
                            marginal_test = function(...) stop("must not run"))
    expect_identical(fast$working_error, a$working_error)
    expect_identical(fast$working_variance, a$working_variance)
    expect_identical(fast$lambda, a$lambda)
    expect_true(all(is.na(fast$diagnostics$marginal_p_value)))
    expect_null(fast$marginal_data)
  }
})

test_that("diagnostics FALSE does not invoke summary.gam", {
  f <- st_fixture()
  trace("summary.gam", where = asNamespace("mgcv"),
        tracer = quote(stop("summary was called")), print = FALSE)
  on.exit(untrace("summary.gam", where = asNamespace("mgcv")), add = TRUE)
  for (G in list(f$G, f$model)) {
    fit <- mgcvST.estimate(f$Y, G, diagnostics = FALSE, marginal = FALSE)
    expect_true(all(is.finite(fit$working_error)))
    if (!is.null(fit$diagnostics$wood_p_value)) expect_true(all(is.na(fit$diagnostics$wood_p_value)))
  }
})

test_that("model score hot path and all pair output fields remain identical", {
  old <- old_st()
  f <- st_fixture(nuisance = TRUE)
  fit <- mgcvST.estimate(f$Y, f$model)
  pairs <- rbind(c(1L,2L), c(3L,1L), c(2L,3L))
  duplicated <- rbind(pairs, c(2L,1L))
  expect_error(old$mgcvST.test(fit,pairs=duplicated), "duplicated tests")
  expect_error(mgcvST.test(fit,pairs=duplicated), "duplicated tests")
  for (cal in c("liu", "davies")) {
    if (cal == "davies") skip_if_not_installed("CompQuadForm")
    for (chunk in c(1L, 100L)) {
      a <- old$mgcvST.test(fit, pairs = pairs, highlight = matrix(c(1L,3L),1), calibration = cal, chunk_size = chunk)
      b <- mgcvST.test(fit, pairs = pairs, highlight = matrix(c(1L,3L),1), calibration = cal, chunk_size = chunk)
      expect_identical(strip_elapsed(b), strip_elapsed(a))
    }
  }
  for (i in 1:3) expect_identical(mgcvST:::.mgcvst_model_score_state(fit,i), old$.mgcvst_model_score_state(fit,i))
  fit$smoothing_parameters[2,1] <- -1
  a <- old$mgcvST.test(fit, pairs = pairs)
  b <- mgcvST.test(fit, pairs = pairs)
  expect_identical(strip_elapsed(a), strip_elapsed(b))
})

test_that("ordinary Liu engine is unchanged", {
  old <- old_st()
  for (pc in c(FALSE, TRUE)) {
    f <- st_fixture(pc = pc)
    fit <- mgcvST.estimate(f$Y, f$G)
    pairs <- t(combn(1:3,2))
    expect_identical(strip_elapsed(mgcvST.test(fit,pairs=pairs)),
                     strip_elapsed(old$mgcvST.test(fit,pairs=pairs)))
  }
})

test_that("new switches reject non-logical values", {
  expect_error(mgcvST.estimate(NULL, NULL, marginal = NA), "marginal must")
  expect_error(mgcvST.estimate(NULL, NULL, diagnostics = 0), "diagnostics must")
  expect_error(mgcvST.estimate(NULL, NULL, retain_marginal = 1), "retain_marginal must")
})
