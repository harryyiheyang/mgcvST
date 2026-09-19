test_that("packed dense score states restore exactly", {
  M <- crossprod(matrix(seq_len(20), 5L, 4L))
  state <- list(a = as.numeric(seq_len(4)), M = M, width = c(global = 4L))
  restored <- mgcvST:::.mgcvst_unpack_score_state(
    mgcvST:::.mgcvst_pack_score_state(state)
  )
  expect_identical(restored, state)
})

test_that("dense pair caches are shared with local Snow workers", {
  skip_on_cran()
  f <- st_fixture()
  bp <- BiocParallel::SnowParam(workers = 2L, type = "SOCK", progressbar = FALSE)
  on.exit(BiocParallel::bpstop(bp), add = TRUE)
  pairs <- rbind(c(1L, 2L), c(1L, 3L), c(2L, 3L))
  for (design in list(f$G, f$model)) {
    fit <- mgcvST.estimate(f$Y, design, diagnostics = FALSE,
                          BPPARAM = BiocParallel::SerialParam())
    for (calibration in c("liu", "davies")) {
      serial <- mgcvST.test(fit, pairs = pairs, calibration = calibration,
                            chunk_size = 1L, BPPARAM = BiocParallel::SerialParam())
      parallel <- mgcvST.test(fit, pairs = pairs, calibration = calibration,
                              chunk_size = 1L, BPPARAM = bp)
      expect_equal(parallel$results, serial$results, tolerance = 1e-12)
    }
  }
})

test_that("model pair states are constructed once per unique feature", {
  f <- st_fixture(nuisance = TRUE)
  fit <- mgcvST.estimate(
    f$Y, f$model, BPPARAM = BiocParallel::SerialParam(), diagnostics = FALSE
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
  ans <- mgcvST.test(
    fit, pairs = pairs, calibration = "liu", chunk_size = 1L,
    BPPARAM = BiocParallel::SerialParam()
  )
  expect_identical(count$features, 3L)
  expect_true(all(is.finite(ans$results$p_two_sided)))
})

test_that("packed model pair evaluation preserves direct pair results", {
  f <- st_fixture(nuisance = TRUE)
  fit <- mgcvST.estimate(
    f$Y, f$model, BPPARAM = BiocParallel::SerialParam(), diagnostics = FALSE
  )
  pairs <- rbind(c(1L, 2L), c(1L, 3L), c(2L, 3L))
  fit$.mgcvst_fixed_factors <- mgcvST:::.mgcvst_model_fixed_factors(fit)
  for (calibration in c("liu", "davies")) {
    if (calibration == "davies") skip_if_not_installed("CompQuadForm")
    expected <- lapply(seq_len(nrow(pairs)), function(k) {
      mgcvST:::.mgcvst_model_pair_single(
        fit, pairs[k, 1L], pairs[k, 2L], calibration
      )
    })
    ans <- mgcvST.test(
      fit, pairs = pairs, calibration = calibration, chunk_size = 1L,
      BPPARAM = BiocParallel::SerialParam()
    )
    expect_equal(ans$results$signed_score,
                 vapply(expected, `[[`, numeric(1L), "score"), tolerance = 1e-10)
    expect_equal(ans$results$p_two_sided,
                 vapply(expected, `[[`, numeric(1L), "p_two_sided"),
                 tolerance = 1e-12)
  }
})
