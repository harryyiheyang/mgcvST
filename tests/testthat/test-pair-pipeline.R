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

test_that("shared preparation keeps dense native state contracts", {
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
  testthat::local_mocked_bindings(
    mgcvst_dense_score_batch_cpp = native_batch,
    .package = "mgcvST"
  )
  native <- list(T0 = matrix(1, 3L, 2L), X = matrix(1, 3L, 1L),
                 sp_index = 1L, width = c(global = 2L))
  model <- mgcvST:::.mgcvst_pair_build_batch(
    fit, 1:2, 2L, mode = "model_native", native = native
  )
  expect_equal(seen$scale, c(0.5, 0.5))
  expect_identical(model[[1L]]$a, c(1, 2))
  expect_equal(model[[2L]]$M, diag(c(2, 3)))
  expect_identical(model[[1L]]$width, c(global = 2L))

  legacy <- mgcvST:::.mgcvst_pair_build_batch(
    fit, 1:2, 2L, mode = "legacy_native",
    T0 = matrix(1, 3L, 2L), field_scale = c(0.25, 0.75)
  )
  expect_equal(seen$scale, c(0.25, 0.75))
  expect_equal(legacy[[1L]]$M, diag(c(1, 2)))
})

test_that("the fused C++ Liu pair kernel matches the old trace-powers + R Liu path", {
  set.seed(20260924)
  q <- 6L
  K <- 8L
  H <- lapply(seq_len(K), function(k) {
    z <- matrix(rnorm(q * q), q, q)
    crossprod(z) + diag(q) * 0.1
  })
  a <- matrix(rnorm(q * K), q, K)

  idx <- which(upper.tri(matrix(0, K, K)), arr.ind = TRUE)
  left <- idx[, 2L]
  right <- idx[, 1L]
  ord <- order(left)
  left <- left[ord]
  right <- right[ord]

  new <- mgcvST:::mgcvst_pair_liu_cpp(H, a, left, right, threads = 1L)

  pairs <- cbind(left, right)
  old_score <- colSums(a[, left, drop = FALSE] * a[, right, drop = FALSE])
  old_moments <- mgcvST:::mgcvst_pair_trace_powers_cpp(H, pairs, maxPower = 4L, threads = 1L)
  old_liu <- mgcvST:::.liu_squared_score_moments(
    abs(old_score), old_moments[, 1L], old_moments[, 2L],
    old_moments[, 3L], old_moments[, 4L]
  )
  old_information <- old_moments[, 1L]
  old_effective_rank <- old_moments[, 1L]^2 / old_moments[, 2L]

  expect_equal(new$score, old_score, tolerance = 1e-12)
  expect_equal(new$information, old_information, tolerance = 1e-12)
  expect_equal(new$effective_rank, old_effective_rank, tolerance = 1e-12)
  finite_p <- old_liu$p_value > 1e-300
  expect_equal(
    exp(new$log_p_two_sided)[finite_p], old_liu$p_value[finite_p],
    tolerance = 1e-10
  )

  # A constructed strong pair whose old p underflows to 0 has a finite
  # log_p_two_sided well below log(1e-300).
  strong_H <- lapply(1:2, function(k) diag(rep(100, q)))
  strong_a <- cbind(rep(60, q), rep(60, q))
  strong <- mgcvST:::mgcvst_pair_liu_cpp(
    strong_H, strong_a, 1L, 2L, threads = 1L
  )
  strong_moments <- mgcvST:::mgcvst_pair_trace_powers_cpp(
    strong_H, matrix(c(1L, 2L), nrow = 1L), maxPower = 4L, threads = 1L
  )
  strong_old <- mgcvST:::.liu_squared_score_moments(
    abs(sum(strong_a[, 1L] * strong_a[, 2L])), strong_moments[, 1L],
    strong_moments[, 2L], strong_moments[, 3L], strong_moments[, 4L]
  )
  expect_equal(strong_old$p_value, 0)
  expect_true(is.finite(strong$log_p_two_sided))
  expect_lt(strong$log_p_two_sided, log(1e-300))
})
