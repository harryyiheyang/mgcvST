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
  expect_length(fit$marginal_data$geometry, 1L)
  expect_length(fit$marginal_data$state, 3L)
  expect_false(any(vapply(fit$marginal_data$state, inherits, logical(1), what = "gam")))
  a <- mgcvST.marginal(fit, calibration = "liu")
  b <- mgcvST.marginal(unserialize(serialize(fit, NULL)), calibration = "liu", features = c(3L,1L))
  expect_identical(unname(b$p_value), unname(a$p_value[c(3,1)]))
  expect_true(all(is.finite(a$p_value)))
  expect_error(mgcvST.marginal(mgcvST.estimate(f$Y,f$G)), "retain_marginal")
})

test_that("native marginal evaluation matches pinned upstream TAPS", {
  repo <- Sys.getenv("MGCVST_TAPS_BASELINE")
  skip_if(!dir.exists(file.path(repo, "R")), "Set MGCVST_TAPS_BASELINE for upstream equivalence")
  skip_if_not_installed("mgcv.taps")
  upstream <- new.env(parent = asNamespace("mgcv.taps"))
  for (name in c("basicfunction.R", "gammfast.R", "extract_pseudo_response.R", "taps_score_test.R")) {
    sys.source(file.path(repo,"R",name), upstream)
  }
  for (fam in list(stats::gaussian(), stats::poisson(), mgcv::nb())) {
    for (pc in c(FALSE, TRUE)) {
      f <- st_fixture(family = fam, pc = pc)
      fit <- mgcvST.estimate(f$Y, f$G, retain_marginal = TRUE)
      reference <- numeric(3)
      for (k in 1:3) {
        G <- f$G
        G$y <- as.numeric(f$Y[k,])
        G$mf[[attr(G$terms,"response")]] <- G$y
        G$family <- unserialize(serialize(f$G$family,NULL))
        gam <- mgcv::gam(G=G, method="REML", control=mgcv::gam.control(nthreads=1L))
        reference[k] <- upstream$taps_score_test(gam, method="liu")$smooth.pvalue
      }
      got <- mgcvST.marginal(fit,calibration="liu")
      expect_identical(got$p_value, reference)
      expect_true(all(is.na(got$error_message)))
      if (requireNamespace("CompQuadForm",quietly=TRUE)) {
        # The final feature GAM is still in scope. Compare successful Davies
        # calibration, without hiding numerical failures with Liu.
        davies <- mgcvST.marginal(fit,features=3L,calibration="davies")
        if (is.na(davies$error_message)) {
          expect_identical(davies$p_value,
            upstream$taps_score_test(gam,method="davies")$smooth.pvalue)
        }
      }
    }
  }
  f <- st_fixture(nuisance=TRUE)
  fit <- mgcvST.estimate(f$Y,f$model,retain_marginal=TRUE)
  G <- f$model$G
  G$y <- as.numeric(f$Y[1,])
  G$mf[[attr(G$terms,"response")]] <- G$y
  gam <- mgcv::gam(G=G,method="REML",control=mgcv::gam.control(nthreads=1L))
  component <- fit$marginal_data$geometry[[1]]$test_component
  expect_identical(mgcvST.marginal(fit,features=1L,calibration="liu")$p_value,
                   upstream$taps_score_test(gam,test.component=component,method="liu")$smooth.pvalue)
})

test_that("Davies failures never switch calibration without explicit consent", {
  skip_if_not_installed("CompQuadForm")
  z <- list(statistic = 100, lambda = c(1,2,3))
  testthat::local_mocked_bindings(davies = function(...) list(Qq=0,ifault=1L), .package="CompQuadForm")
  none <- mgcvST:::.mgcvst_marginal_davies(z,"none",1e-8,1e5)
  yes <- mgcvST:::.mgcvst_marginal_davies(z,"liu",1e-8,1e5)
  expect_true(is.na(none$p_value))
  expect_identical(none$method_used,"davies")
  expect_false(none$fallback_used)
  expect_true(yes$fallback_used)
  expect_identical(yes$method_used,"liu")
  expect_true(is.finite(yes$p_value))
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
