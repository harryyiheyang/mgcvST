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
    ref <- mgcvST:::.mgcvst_model_score_state(fit, ids[k])
    expect_equal(ans[[k]]$a, ref$a, tolerance = 1e-10)
    expect_equal(ans[[k]]$H, unname(ref$M), tolerance = 1e-10)
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
