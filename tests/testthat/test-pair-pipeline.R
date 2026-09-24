.pair_pipeline_fit <- function() {
  ids <- c("g1", "g2", "g3")
  list(
    feature_id = ids, test_engine = "single_model",
    estimator = "INLA", score_backend = "sparse",
    working_error = matrix(0, 4L, 3L),
    working_variance = matrix(1, 4L, 3L),
    dispersion = rep(1, 3L), lambda = rep(1, 3L),
    smoothing_parameters = matrix(1, 3L, 1L),
    nuisance_covariance = list(),
    geometry = list(nuisance_design = matrix(numeric(), 4L, 0L)),
    score_sparse = list(A = Matrix::Diagonal(4L, 3L),
                        Q = Matrix::Diagonal(3L),
                        constraint = rep(1, 3L), sp_index = 1L)
  )
}

test_that("bounded pair blocks preserve self, duplicate, and reversed pairs", {
  fit <- .pair_pipeline_fit()
  basis <- list(coordinate = matrix(c(1, 0, 0, 0, 1, 0), 3L, 2L),
                basis = matrix(c(1, 0, 0, 0, 1, 0), 3L, 2L),
                rank = 2L, coverage = 0.995)
  pairs <- rbind(c(1L, 1L), c(1L, 2L), c(2L, 1L), c(1L, 2L),
                 c(2L, 2L), c(2L, 3L), c(3L, 2L), c(3L, 3L), c(1L, 3L))
  original <- mgcvST:::.mgcvst_liu_pairs
  seen <- list()
  testthat::local_mocked_bindings(
    .inlast_sparse_prepare = function(fit) fit,
    .inlast_sparse_units = function(fit, features, threads = 1L)
      lapply(features, function(i) list(feature = i)),
    .inlast_sparse_materialize_reduced = function(fit, units, basis, threads = 1L)
      lapply(units, function(z) {
        i <- z$feature
        list(a = c(i, i + 0.25), M = diag(c(i + 0.5, i + 1)), width = 2L)
      }),
    .mgcvst_liu_pairs = function(index, pair_index, feature_id, summaries,
                                 threads, chunk_size, verbose) {
      seen[[length(seen) + 1L]] <<- c(rows = nrow(index),
                                      features = length(summaries$used))
      original(index, pair_index, feature_id, summaries, threads,
               chunk_size, verbose)
    },
    .package = "mgcvST"
  )
  full <- mgcvST:::.mgcvst_pair_pipeline(
    fit, pairs, seq_len(nrow(pairs)), threads = 1L, chunk_size = 9L,
    verbose = FALSE, basis = basis, cache_bytes = 100000
  )
  seen <- list()
  bounded <- mgcvST:::.mgcvst_pair_pipeline(
    fit, pairs, seq_len(nrow(pairs)), threads = 1L, chunk_size = 9L,
    verbose = FALSE, basis = basis, cache_bytes = 4200
  )
  expect_identical(bounded$result$pair_index, seq_len(nrow(pairs)))
  expect_equal(bounded$result, full$result, tolerance = 1e-12)
  expect_gt(length(seen), 1L)
  expect_true(all(vapply(seen, `[[`, integer(1L), "features") <= 2L))
  expect_true(all(vapply(seen, `[[`, integer(1L), "rows") <= 9L))
})

test_that("shared preparation keeps dense native and fallback state contracts", {
  fit <- list(
    working_error = matrix(0, 3L, 2L),
    working_variance = matrix(1, 3L, 2L),
    dispersion = c(1, 2), smoothing_parameters = matrix(c(2, 4), 2L),
    nuisance_covariance = list(matrix(1), matrix(1)),
    geometry = list(X = matrix(1, 3L, 1L))
  )
  seen <- new.env(parent = emptyenv())
  seen$scale <- NULL
  native_batch <- function(T0, variance, error, scale, X, nuisance, threads) {
    seen$scale <- scale
    lapply(seq_along(scale), function(i) list(a = c(i, i + 1),
                                              H = diag(c(i, i + 1))))
  }
  fallback <- function(fit, feature) {
    list(a = c(feature, feature + 1), M = diag(2), width = 2L)
  }
  testthat::local_mocked_bindings(
    mgcvst_dense_score_batch_cpp = native_batch,
    .mgcvst_model_score_state = fallback,
    .package = "mgcvST"
  )
  native <- list(T0 = matrix(1, 3L, 2L), X = matrix(1, 3L, 1L),
                 sp_index = 1L, width = c(global = 2L))
  model <- mgcvST:::.mgcvst_pair_build_batch(
    fit, 1:2, 2L, NULL, "model_native", native = native
  )
  expect_equal(seen$scale, c(0.5, 0.5))
  expect_identical(model[[1L]]$a, c(1, 2))
  expect_equal(model[[2L]]$M, diag(c(2, 3)))
  expect_identical(model[[1L]]$width, c(global = 2L))

  legacy <- mgcvST:::.mgcvst_pair_build_batch(
    fit, 1:2, 2L, NULL, "legacy_native",
    T0 = matrix(1, 3L, 2L), field_scale = c(0.25, 0.75)
  )
  expect_equal(seen$scale, c(0.25, 0.75))
  expect_equal(legacy[[1L]]$M, diag(c(1, 2)))
  fallback_states <- mgcvST:::.mgcvst_pair_build_batch(
    fit, 1:2, 1L, NULL, "model_fallback"
  )
  expect_equal(fallback_states[[2L]]$M, diag(2))
})
