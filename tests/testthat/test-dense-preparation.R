test_that("dense batches agree with observation-space precision matrices", {
  set.seed(914)
  n <- 18L
  q <- 4L
  T0 <- matrix(rnorm(n * q), n, q)
  E <- matrix(rnorm(n * 3L), n, 3L)
  V <- matrix(runif(n * 3L, 0.4, 2), n, 3L)
  scale <- c(0.2, 1, 3)
  for (X in list(matrix(numeric(), n, 0L), cbind(1, seq_len(n) / n),
                 cbind(rep(1, n), rep(0, n)))) {
    ans <- mgcvST:::mgcvst_dense_score_batch_cpp(T0, V, E, scale, X, list(), 1L)
    for (i in seq_len(ncol(E))) {
      F <- sqrt(scale[i]) * T0
      W <- solve(diag(V[, i]) + tcrossprod(F))
      P <- W
      if (ncol(X)) {
        WX <- W %*% X
        P <- W - WX %*% CppMatrix::matrixGeneralizedInverse(crossprod(X, WX)) %*%
          t(WX)
      }
      expect_null(ans[[i]]$error)
      expect_null(dim(ans[[i]]$a))
      expect_equal(ans[[i]]$a, drop(crossprod(F, P %*% E[, i])), tolerance = 1e-10)
      expect_equal(ans[[i]]$H, crossprod(F, P %*% F), tolerance = 1e-10)
    }
  }
  X <- cbind(1, seq_len(n) / n)
  vp <- list(diag(c(0.03, 0.1)), diag(c(0.01, 0.05)), diag(c(0.1, 0.02)))
  ans <- mgcvST:::mgcvst_dense_score_batch_cpp(T0, V, E, scale, X, vp, 1L)
  for (i in seq_len(ncol(E))) {
    F <- sqrt(scale[i]) * T0
    W <- solve(diag(V[, i]) + tcrossprod(F))
    WX <- W %*% X
    P <- W - WX %*% vp[[i]] %*% t(WX)
    expect_equal(ans[[i]]$a, drop(crossprod(F, P %*% E[, i])), tolerance = 1e-10)
    expect_equal(ans[[i]]$H, crossprod(F, P %*% F), tolerance = 1e-10)
  }
})

test_that("dense preparation isolates failed features and preserves thread results", {
  skip_on_cran()
  set.seed(915)
  T0 <- matrix(rnorm(60), 15L, 4L)
  X <- matrix(1, 15L, 1L)
  E <- matrix(rnorm(90), 15L, 6L)
  V <- matrix(1, 15L, 6L)
  scale <- rep(1, 6L)
  V[1, 2] <- 0
  E[2, 3] <- NA_real_
  scale[4] <- -1
  vp <- rep(list(matrix(0.05, 1L, 1L)), 6L)
  vp[[5]] <- matrix(1, 2L, 2L)
  one <- mgcvST:::mgcvst_dense_score_batch_cpp(T0, V, E, scale, X, vp, 1L)
  two <- mgcvST:::mgcvst_dense_score_batch_cpp(T0, V, E, scale, X, vp, 2L)
  expect_equal(two, one, tolerance = 1e-12)
  expect_null(one[[1]]$error)
  expect_null(one[[6]]$error)
  for (i in 2:5) expect_type(one[[i]]$error, "character")
  expect_error(mgcvST:::mgcvst_dense_score_batch_cpp(
    T0, V[-1, ], E, scale, X, vp, 1L), "dimensions")
})

test_that("score_only batches match the full a vector for both nuisance modes", {
  set.seed(916)
  n <- 18L
  q <- 4L
  T0 <- matrix(rnorm(n * q), n, q)
  E <- matrix(rnorm(n * 3L), n, 3L)
  V <- matrix(runif(n * 3L, 0.4, 2), n, 3L)
  scale <- c(0.2, 1, 3)
  X <- cbind(1, seq_len(n) / n)
  full <- mgcvST:::mgcvst_dense_score_batch_cpp(T0, V, E, scale, X, list(), 1L)
  only <- mgcvST:::mgcvst_dense_score_batch_cpp(
    T0, V, E, scale, X, list(), 1L, score_only = TRUE
  )
  for (i in seq_len(ncol(E))) {
    expect_null(only[[i]]$error)
    expect_null(only[[i]]$H)
    expect_equal(only[[i]]$a, full[[i]]$a, tolerance = 1e-12)
  }

  vp <- list(diag(c(0.03, 0.1)), diag(c(0.01, 0.05)), diag(c(0.1, 0.02)))
  full2 <- mgcvST:::mgcvst_dense_score_batch_cpp(T0, V, E, scale, X, vp, 1L)
  only2 <- mgcvST:::mgcvst_dense_score_batch_cpp(
    T0, V, E, scale, X, vp, 1L, score_only = TRUE
  )
  for (i in seq_len(ncol(E))) {
    expect_null(only2[[i]]$error)
    expect_null(only2[[i]]$H)
    expect_equal(only2[[i]]$a, full2[[i]]$a, tolerance = 1e-12)
  }
})

test_that("model preparation keeps the conditional nuisance covariance", {
  f <- st_fixture(nuisance = TRUE)
  fit <- mgcvST.estimate(f$Y, f$model, diagnostics = FALSE,
                         BPPARAM = BiocParallel::SerialParam())
  fit$.mgcvst_fixed_factors <- mgcvST:::.mgcvst_model_fixed_factors(fit)
  ids <- c(3L, 1L, 2L)
  native <- mgcvST:::.mgcvst_model_dense_preparation(fit, ids)
  expect_type(native, "list")
  ans <- mgcvST:::mgcvst_dense_score_batch_cpp(
    native$T0, fit$working_variance[, ids], fit$working_error[, ids],
    fit$dispersion[ids] / fit$smoothing_parameters[ids, native$sp_index],
    native$X, fit$nuisance_covariance[ids], 1L)
  for (k in seq_along(ids)) {
    expect_true(is.numeric(ans[[k]]$a))
    expect_true(is.matrix(ans[[k]]$H))
  }
  old <- fit
  old$nuisance_covariance <- NULL
  expect_null(mgcvST:::.mgcvst_model_dense_preparation(old, ids))
  pairs <- rbind(c(1L, 2L), c(2L, 3L))
  out <- mgcvST.test(fit, pairs = pairs, threads = 1L)
  expect_identical(out$timing$preparation_backend, "C++ OpenMP")
  expect_identical(out$timing$preparation_threads, 1L)
  expect_gte(out$timing$summary_elapsed, 0)
  expect_equal(out$timing$elapsed,
               out$timing$summary_elapsed + out$timing$pair_elapsed)
})

test_that("a feature without a usable nuisance covariance fails at estimation, no fallback", {
  f <- st_fixture(nuisance = TRUE)
  original <- mgcvST:::.mgcvst_nuisance_state
  # .mgcvst_estimate_model() runs every per-feature fit through a worker
  # bundle (.mgcvst_worker_bundle()) that rebinds each copied function's
  # environment to an isolated bundle env; an ordinary closure over a local
  # counter would lose that counter once the mock's environment is
  # reassigned. Stash both the call counter and the true original function
  # as attributes of the mock closure itself, which survive
  # environment(mock) <- bundle, and read them back via sys.function().
  options(.mgcvst_test_nuisance_calls = 0L)
  on.exit(options(.mgcvst_test_nuisance_calls = NULL), add = TRUE)
  mock <- function(fit, geometry, cache) {
    real <- attr(sys.function(), "real")
    n <- getOption(".mgcvst_test_nuisance_calls", 0L) + 1L
    options(.mgcvst_test_nuisance_calls = n)
    if (n == 2L) return(list(error = "rank-deficient fit (rank 3 of 4)"))
    real(fit, geometry, cache)
  }
  attr(mock, "real") <- original
  testthat::local_mocked_bindings(
    .mgcvst_nuisance_state = mock,
    .package = "mgcvST"
  )
  fit <- mgcvST.estimate(f$Y, f$model, diagnostics = FALSE,
                         BPPARAM = BiocParallel::SerialParam())
  expect_identical(fit$feature_id[2L], "response2")
  expect_true(grepl(
    "nuisance covariance unavailable: rank-deficient fit",
    fit$diagnostics$error_message[2L], fixed = TRUE
  ))
  available <- mgcvST:::.mgcvst_feature_available(fit)
  expect_false(available[2L])
  expect_true(all(available[-2L]))

  pairs <- rbind(c(1L, 2L), c(2L, 3L), c(1L, 3L))
  out <- mgcvST.test(fit, pairs = pairs, threads = 1L)
  bad <- out$results$feature1 == "response2" | out$results$feature2 == "response2"
  expect_true(any(bad) && !all(bad))
  expect_true(all(is.na(out$results$p_two_sided[bad])))
  expect_true(all(!is.na(out$results$error_message[bad])))
  expect_false(any(out$results$discovered[bad]))
  expect_true(all(is.finite(out$results$p_two_sided[!bad])))
})
