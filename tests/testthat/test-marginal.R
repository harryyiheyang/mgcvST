test_that("marginal spectral powers agree exactly with R at all thread counts", {
  set.seed(10)
  xs <- lapply(c(1, 5, 31, 200), function(n) exp(runif(n, -50, 50)))
  powers <- lapply(xs, function(x) cbind(x, x^2, x^3, x^4))
  ref <- do.call(rbind, lapply(xs, mgcvST:::.mgcvst_marginal_moments))
  expect_identical(mgcvST:::mgcvst_marginal_liu_moments_cpp(powers, 1L), ref)
  if (.Platform$OS.type == "windows") {
    expect_identical(mgcvST:::mgcvst_marginal_liu_moments_cpp(powers, 2L), ref)
  }
})

test_that("retained marginal data is opt-in, compact and survives serialization", {
  f <- st_fixture()
  fit <- mgcvST.estimate(f$Y, f$G, retain_marginal = TRUE, chunk_size = 1L)
  expect_identical(fit$marginal_data$version, 2L)
  expect_length(fit$marginal_data$state, 3L)
  expect_false(any(vapply(fit$marginal_data$state, inherits, logical(1), what = "gam")))
  expect_true(all(vapply(fit$marginal_data$state, function(z) {
    is.numeric(z$marginal_cache$statistic) &&
      length(z$marginal_cache$lambda) > 0L
  }, logical(1L))))
  a <- mgcvST.marginal(fit, calibration = "liu")
  b <- mgcvST.marginal(unserialize(serialize(fit, NULL)), calibration = "liu", features = c(3L,1L))
  expect_identical(unname(b$p_value), unname(a$p_value[c(3,1)]))
  expect_true(all(is.finite(a$p_value)))
  expect_error(mgcvST.marginal(mgcvST.estimate(f$Y,f$G)), "retain_marginal")

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
    payload, fit$marginal_data$geometry, "liu", "none",
    1e-10, 1e-8, 1e5, 1L
  )
  expect_identical(calls, 0L)
  expect_equal(cached[[1L]]$statistic,
               fit$marginal_data$state[[1L]]$marginal_cache$statistic)
  recomputed <- mgcvST:::.mgcvst_marginal_chunk(
    payload, fit$marginal_data$geometry, "liu", "none",
    1e-9, 1e-8, 1e5, 1L
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
    for (pc in c(FALSE, TRUE)) {
      f <- st_fixture(family = if (fam == "gaussian") gaussian() else mgcv::nb(), pc = pc)
      fit <- mgcvST.estimate(f$Y, f$G, retain_marginal = TRUE,
                             marginal_args = list(method = "liu"))
      got <- mgcvST.marginal(fit, calibration = "liu")
      expect_true(all(is.finite(fit$diagnostics$marginal_p_value)))
      expect_equal(got$p_value, fit$diagnostics$marginal_p_value,
                   tolerance = 1e-8)
      expect_true(all(is.na(got$error_message)))
    }
  }
  f <- st_fixture(nuisance = TRUE)
  fit <- mgcvST.estimate(f$Y, f$model, marginal_args = list(method = "liu"))
  expect_true(all(is.finite(fit$diagnostics$marginal_p_value)))
})

test_that("Davies failures never switch calibration without explicit consent", {
  skip_if_not_installed("CompQuadForm")
  z <- list(statistic = 100, lambda = c(1,2,3))
  testthat::local_mocked_bindings(davies = function(...) list(Qq=0,ifault=1L), .package="CompQuadForm")
  none <- mgcvST:::.mgcvst_marginal_davies(z,"none",1e-8,1e5)
  yes <- mgcvST:::.mgcvst_marginal_davies(z,"saddlepoint",1e-8,1e5)
  expect_true(is.na(none$p_value))
  expect_identical(none$method_used,"davies")
  expect_false(none$fallback_used)
  expect_true(yes$fallback_used)
  expect_identical(yes$method_used,"saddlepoint")
  expect_true(is.finite(yes$p_value))
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

test_that("Snow workers use retained state and chunk caches", {
  skip_on_cran()
  f <- st_fixture(nuisance=TRUE)
  bp <- BiocParallel::SnowParam(workers=2,type="SOCK",progressbar=FALSE)
  on.exit(BiocParallel::bpstop(bp),add=TRUE)
  Y <- f$Y[rep(1:3, 2),,drop=FALSE]
  rownames(Y) <- paste0("snow", seq_len(nrow(Y)))
  fit <- mgcvST.estimate(Y,f$model,retain_marginal=TRUE,BPPARAM=bp,chunk_size=1)
  expect_true(all(vapply(fit$nuisance_covariance, is.matrix, logical(1L))))
  pairs <- t(combn(1:3,2))
  serial <- mgcvST.test(fit,pairs=pairs)
  snow <- mgcvST.test(fit,pairs=pairs,BPPARAM=bp,chunk_size=1)
  expect_identical(serial$results,snow$results)
  for (cal in c("liu","davies")) {
    a <- mgcvST.marginal(fit,calibration=cal)
    b <- mgcvST.marginal(fit,calibration=cal,BPPARAM=bp,chunk_size=1)
    expect_identical(a,b)
  }
})
