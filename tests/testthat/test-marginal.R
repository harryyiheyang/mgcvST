test_that("retained marginal data is opt-in, compact and survives serialization", {
  f <- st_fixture()
  fit <- mgcvST.estimate(f$Y, f$G, retain_marginal = TRUE, chunk_size = 1L, spatial = "all")
  expect_identical(fit$marginal_data$version, 2L)
  expect_length(fit$marginal_data$state, 3L)
  expect_false(any(vapply(fit$marginal_data$state, inherits, logical(1), what = "gam")))
  expect_true(all(vapply(fit$marginal_data$state, function(z) {
    is.numeric(z$marginal_cache$statistic) &&
      length(z$marginal_cache$lambda) > 0L
  }, logical(1L))))
  a <- mgcvST.marginal(fit)
  b <- mgcvST.marginal(unserialize(serialize(fit, NULL)), features = c(3L,1L))
  expect_identical(unname(b$p_value), unname(a$p_value[c(3,1)]))
  expect_true(all(is.finite(a$p_value)))
  expect_error(mgcvST.marginal(mgcvST.estimate(f$Y,f$G, spatial = "all")), "retain_marginal")

  calls <- 0L
  original <- mgcvST:::.mgcvst_marginal_spectrum
  testthat::local_mocked_bindings(
    .mgcvst_marginal_spectrum = function(...) {
      calls <<- calls + 1L
      original(...)
    },
    .package = "mgcvST"
  )
  payload <- list(index = 1L, state = fit$marginal_data$state[1L], version = 2L)
  cached <- mgcvST:::.mgcvst_marginal_chunk(
    payload, fit$marginal_data$geometry, 1e-10, 1e-8, 1e5, 1L
  )
  expect_identical(calls, 0L)
  expect_equal(cached[[1L]]$statistic,
               fit$marginal_data$state[[1L]]$marginal_cache$statistic)
  recomputed <- mgcvST:::.mgcvst_marginal_chunk(
    payload, fit$marginal_data$geometry, 1e-9, 1e-8, 1e5, 1L
  )
  expect_identical(calls, 0L)
  expect_true(is.finite(recomputed[[1L]]$statistic))
})

test_that("custom marginal callbacks do not populate the built-in spectrum cache", {
  callback <- function(fit, test.component, n_threads) {
    data.frame(smooth.pvalue = 0.25, method = "custom")
  }
  z <- mgcvST:::.mgcvst_marginal_score(list(), callback, list())
  expect_null(z$cache)
})

test_that("null-first TAPS is finite and its retained calibration is stable", {
  for (fam in c("gaussian", "nb")) {
    f <- st_fixture(family = if (fam == "gaussian") gaussian() else mgcv::nb())
    fit <- mgcvST.estimate(f$Y, f$G, retain_marginal = TRUE, spatial = "all")
    got <- mgcvST.marginal(fit)
    expect_true(all(is.finite(fit$diagnostics$marginal_p_value)))
    expect_equal(got$p_value, fit$diagnostics$marginal_p_value,
                 tolerance = 1e-8)
    expect_identical(got$method_used, fit$diagnostics$marginal_method)
    expect_identical(got$fallback_used, fit$diagnostics$marginal_fallback)
    expect_true(all(is.na(got$error_message)))
  }
  f <- st_fixture(nuisance = TRUE)
  fit <- mgcvST.estimate(f$Y, f$model, spatial = "all")
  expect_true(all(is.finite(fit$diagnostics$marginal_p_value)))
  expect_true(all(fit$diagnostics$marginal_requested_method == "davies"))
  expect_true(all(fit$diagnostics$marginal_method %in% c("davies", "saddlepoint")))
})

test_that("a Davies failure uses the saddlepoint and sets the fallback flag", {
  z <- list(statistic = 100, lambda = c(1, 2, 3))
  sp <- mgcvST:::.mgcvst_marginal_saddlepoint(z$statistic, z$lambda)
  expect_true(is.finite(sp) && sp > 0 && sp < 1)
  failures <- list(list(Qq = 0, ifault = 0L), list(Qq = -1e-12, ifault = 1L),
                   list(Qq = NA_real_, ifault = 0L), list(Qq = NaN, ifault = 0L),
                   list(Qq = Inf, ifault = 0L), list(Qq = 1.1, ifault = 0L),
                   list(Qq = numeric(), ifault = 0L), list(ifault = 0L))
  for (result in failures) {
    local({
      testthat::local_mocked_bindings(davies = function(...) result,
                                      .package = "CompQuadForm")
      out <- mgcvST:::.mgcvst_marginal_davies(z, 1e-8, 1e5)
      expect_identical(out$p_value, sp)
      expect_identical(out$method_used, "saddlepoint")
      expect_true(out$fallback_used)
      expect_match(out$fallback_reason, "Davies numerical failure")
      expect_identical(out$davies_ifault, result$ifault)
      expect_true(is.na(out$error_message))
    })
  }
  local({
    testthat::local_mocked_bindings(davies = function(...) stop("integration failed"),
                                    .package = "CompQuadForm")
    out <- mgcvST:::.mgcvst_marginal_davies(z, 1e-8, 1e5)
    expect_identical(out$p_value, sp)
    expect_identical(out$method_used, "saddlepoint")
    expect_true(out$fallback_used)
    expect_identical(out$fallback_reason, "integration failed")
    expect_identical(out$davies_ifault, NA_integer_)
  })
})

test_that("a Davies p-value in (0, 1] is kept whatever ifault reports", {
  z <- list(statistic = 100, lambda = c(1, 2, 3))
  kept <- list(list(Qq = 0.4, ifault = 1L), list(Qq = 0.4, ifault = 2L),
               list(Qq = 0.4), list(Qq = 1, ifault = 0L),
               list(Qq = 1e-300, ifault = 1L))
  for (result in kept) {
    local({
      testthat::local_mocked_bindings(davies = function(...) result,
                                      .package = "CompQuadForm")
      out <- mgcvST:::.mgcvst_marginal_davies(z, 1e-8, 1e5)
      expect_identical(out$p_value, result$Qq)
      expect_identical(out$method_used, "davies")
      expect_false(out$fallback_used)
      expect_true(is.na(out$fallback_reason))
      expect_identical(out$davies_ifault,
                       if (is.null(result$ifault)) NA_integer_ else result$ifault)
    })
  }
})

test_that("estimation and retained recalibration report the method actually used", {
  f <- st_fixture()
  local({
    testthat::local_mocked_bindings(davies = function(...) list(Qq = 0.4, ifault = 1L),
                                    .package = "CompQuadForm")
    fit <- mgcvST.estimate(f$Y, f$G, BPPARAM = BiocParallel::SerialParam(), spatial = "all")
    expect_identical(unname(fit$diagnostics$marginal_p_value), rep(0.4, 3L))
    expect_identical(fit$diagnostics$marginal_method, rep("davies", 3L))
    expect_identical(fit$diagnostics$marginal_fallback, rep(FALSE, 3L))
  })
  local({
    testthat::local_mocked_bindings(davies = function(...) list(Qq = 0, ifault = 0L),
                                    .package = "CompQuadForm")
    for (setup in list(f$G, f$model)) {
      fit <- mgcvST.estimate(f$Y, setup, retain_marginal = TRUE,
                             BPPARAM = BiocParallel::SerialParam(), spatial = "all")
      expect_identical(fit$diagnostics$marginal_requested_method, rep("davies", 3L))
      expect_identical(fit$diagnostics$marginal_method, rep("saddlepoint", 3L))
      expect_identical(fit$diagnostics$marginal_fallback, rep(TRUE, 3L))
      expect_true(all(is.finite(fit$diagnostics$marginal_p_value)))
      spectrum <- lapply(fit$marginal_data$state, `[[`, "marginal_cache")
      expect_equal(unname(fit$diagnostics$marginal_p_value),
                   unname(vapply(spectrum, function(z) mgcvST:::.mgcvst_marginal_saddlepoint(
                     z$statistic, z$lambda), numeric(1L))), tolerance = 1e-12)
      got <- mgcvST.marginal(fit)
      expect_identical(got$method_requested, rep("davies", 3L))
      expect_identical(got$method_used, rep("saddlepoint", 3L))
      expect_identical(got$fallback_used, rep(TRUE, 3L))
      expect_identical(got$davies_ifault, rep(0L, 3L))
      expect_identical(unname(got$p_value), unname(fit$diagnostics$marginal_p_value))
    }
  })
})

test_that("Liu is not available for the marginal test", {
  f <- st_fixture()
  for (setup in list(f$G, f$model)) {
    expect_error(mgcvST.estimate(f$Y, setup, marginal_args = list(method = "liu"), spatial = "all"),
                 "method = \"liu\" is not available")
  }
  gam_fit <- mgcv::gam(G = f$G, method = "REML")
  expect_error(mgcvST:::taps_score_test(gam_fit, method = "liu"), "should be")
  expect_error(mgcvST:::.mgcvst_null_score_test(NULL, NULL, method = "liu"), "should be")
  fit <- mgcvST.estimate(f$Y, f$G, retain_marginal = TRUE, spatial = "all")
  expect_error(mgcvST.marginal(fit, calibration = "liu"), "unused argument")
  expect_error(mgcvST.marginal(fit, fallback = "saddlepoint"), "unused argument")
  expect_false(exists("mgcvst_marginal_liu_moments_cpp", asNamespace("mgcvST")))
  expect_false(exists(".mgcvst_marginal_liu", asNamespace("mgcvST")))
  explicit <- mgcvST.estimate(f$Y, f$G, marginal_args = list(method = "davies"), spatial = "all")
  expect_identical(explicit$diagnostics$marginal_p_value,
                   fit$diagnostics$marginal_p_value)
})

test_that("the saddlepoint fallback agrees with Davies where Davies is accurate", {
  skip_if_not_installed("CompQuadForm")
  set.seed(11)
  lambda <- sort(rexp(300)^2, decreasing = TRUE)
  mu <- sum(lambda); sd <- sqrt(2 * sum(lambda^2))
  for (q in mu + sd * c(-1, 0.5, 3, 6, 10)) {
    d <- CompQuadForm::davies(q, lambda, lim = 1e6, acc = 1e-14)
    expect_identical(d$ifault, 0L)
    sp <- mgcvST:::.mgcvst_marginal_saddlepoint(q, lambda)
    expect_lt(abs(log10(sp) - log10(d$Qq)), 0.05)
  }
})

test_that("the saddlepoint fallback stays finite beyond machine precision", {
  lambda <- c(5, 2, rep(0.5, 50))
  q <- 1200
  sp <- mgcvST:::.mgcvst_marginal_saddlepoint(q, lambda)
  expect_true(is.finite(sp) && sp > 0 && sp < 1e-40)
  # exponential tail rate 1/(2 * max(lambda)) (Chen and Lumley 2019, Theorem 1)
  sp2 <- mgcvST:::.mgcvst_marginal_saddlepoint(q + 10, lambda)
  expect_lt(abs((log(sp) - log(sp2)) / 10 - 1 / (2 * max(lambda))), 0.01)
})

test_that("the saddlepoint uses its limit at the mean and returns 1 for q <= 0", {
  lambda <- c(5, 2, rep(0.5, 50))
  mu <- sum(lambda)
  rho3 <- 8 * sum(lambda^3) / (2 * sum(lambda^2))^1.5
  centre <- mgcvST:::.mgcvst_marginal_saddlepoint(mu, lambda)
  expect_equal(centre, 0.5 - rho3 / (6 * sqrt(2 * pi)), tolerance = 1e-12)
  for (f in c(0.999, 1.001)) {
    expect_lt(abs(mgcvST:::.mgcvst_marginal_saddlepoint(mu * f, lambda) - centre), 0.01)
  }
  expect_identical(mgcvST:::.mgcvst_marginal_saddlepoint(0, lambda), 1)
  expect_identical(mgcvST:::.mgcvst_marginal_saddlepoint(-2, lambda), 1)
  expect_true(is.na(mgcvST:::.mgcvst_marginal_saddlepoint(3, numeric(0))))
})

test_that("Snow workers use retained state and chunk caches", {
  skip_on_cran()
  f <- st_fixture(nuisance=TRUE)
  bp <- BiocParallel::SnowParam(workers=2,type="SOCK",progressbar=FALSE)
  on.exit(BiocParallel::bpstop(bp),add=TRUE)
  Y <- f$Y[rep(1:3, 2),,drop=FALSE]
  rownames(Y) <- paste0("snow", seq_len(nrow(Y)))
  fit <- mgcvST.estimate(Y,f$model,retain_marginal=TRUE,BPPARAM=bp,chunk_size=1, spatial = "all")
  expect_true(all(vapply(fit$nuisance_covariance, is.matrix, logical(1L))))
  pairs <- t(combn(1:3,2))
  serial <- mgcvST.test(fit,pairs=pairs, moments = "exact")
  blocks <- mgcvST.test(fit,pairs=pairs,chunk_size=1, moments = "exact")
  expect_identical(serial$results,blocks$results)
  a <- mgcvST.marginal(fit)
  b <- mgcvST.marginal(fit,BPPARAM=bp,chunk_size=1)
  expect_identical(a,b)
})
