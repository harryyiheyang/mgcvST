test_that("packed dense score states restore exactly", {
  M <- crossprod(matrix(seq_len(20), 5L, 4L))
  state <- list(a = as.numeric(seq_len(4)), M = M, width = c(global = 4L))
  restored <- mgcvST:::.mgcvst_unpack_score_state(
    mgcvST:::.mgcvst_pack_score_state(state)
  )
  expect_identical(restored, state)
})

test_that("a raw gam setup and a prepared model share one estimation and test path", {
  f <- st_fixture()
  pairs <- rbind(c(1L, 2L), c(1L, 3L), c(2L, 3L))
  for (design in list(f$G, f$model)) {
    fit <- mgcvST.estimate(f$Y, design, diagnostics = FALSE,
                          BPPARAM = BiocParallel::SerialParam(), spatial = "all")
    expect_s3_class(fit, "mgcvST_model_fit")
    expect_identical(fit$test_engine, "single_model")
    expect_true(is.list(fit$geometry$smooth))
    one <- mgcvST.test(fit, pairs = pairs, chunk_size = 1L, moments = "exact")
    block <- mgcvST.test(fit, pairs = pairs, chunk_size = 100L, moments = "exact")
    expect_equal(one$results, block$results, tolerance = 1e-12)
    expect_true(all(is.finite(one$results$log_p_two_sided)))
  }
})

test_that("model pair states are constructed once per unique feature", {
  f <- st_fixture(nuisance = TRUE)
  fit <- mgcvST.estimate(
    f$Y, f$model, BPPARAM = BiocParallel::SerialParam(), diagnostics = FALSE, spatial = "all"
  )
  pairs <- rbind(c(1L, 2L), c(1L, 3L), c(2L, 3L))
  count <- new.env(parent = emptyenv())
  count$features <- 0L
  original_batch <- mgcvST:::mgcvst_dense_score_batch_cpp
  testthat::local_mocked_bindings(
    mgcvst_dense_score_batch_cpp = function(T0, variance, error, scale, X,
                                            nuisance, threads) {
      count$features <- count$features + ncol(variance)
      original_batch(T0, variance, error, scale, X, nuisance, threads)
    }, .package = "mgcvST"
  )
  ans <- mgcvST.test(fit, pairs = pairs, chunk_size = 1L, moments = "exact")
  expect_identical(count$features, 3L)
  expect_true(all(is.finite(ans$results$log_p_two_sided)))
})
