.pair_pipeline_fit <- function(p = 4L) {
  ids <- paste0("g", seq_len(p))
  list(
    feature_id = ids, test_engine = "single_model",
    estimator = "mgcv", score_backend = "dense",
    working_error = matrix(0, 4L, p),
    working_variance = matrix(1, 4L, p),
    dispersion = rep(1, p), lambda = rep(1, p),
    smoothing_parameters = matrix(1, p, 1L),
    nuisance_covariance = list(),
    geometry = list(
      target = c(global = 1L),
      smooth = list(list(B = matrix(c(1, 0, 0, 0, 0, 1, 0, 0), 4L, 2L)))
    )
  )
}

# The pipeline builds states through these three steps; the fixture replaces
# the dense score kernel with deterministic states.
.pair_pipeline_mocks <- function(env = parent.frame()) {
  testthat::local_mocked_bindings(
    .mgcvst_model_fixed_factors = function(fit) list(NULL),
    .mgcvst_model_dense_preparation = function(fit, features) {
      list(T0 = NULL, X = matrix(numeric(), 4L, 0L), sp_index = 1L,
           width = c(global = 2L))
    },
    .mgcvst_pair_build_batch = function(fit, ids, threads, native) {
      lapply(ids, function(i) {
        list(a = c(i, i + 0.25), M = diag(c(i + 0.5, i + 1)), width = 2L)
      })
    },
    .package = "mgcvST", .env = env
  )
}

test_that("bounded pair blocks stream every pair once whatever the state budget", {
  .pair_pipeline_mocks()
  fit <- .pair_pipeline_fit()
  pairs <- t(utils::combn(4L, 2L))
  original <- mgcvST:::.mgcvst_liu_pairs
  seen <- list()
  testthat::local_mocked_bindings(
    .mgcvst_liu_pairs = function(index, active, states, threads) {
      seen[[length(seen) + 1L]] <<- c(rows = nrow(index), features = length(active))
      original(index, active, states, threads)
    },
    .package = "mgcvST"
  )
  run <- function(index, cache_bytes, chunk_size = 9L) {
    z <- mgcvST:::.mgcvst_pair_pipeline(fit, index, threads = 1L,
      chunk_size = chunk_size, verbose = FALSE, cache_bytes = cache_bytes)
    out <- do.call(rbind, lapply(z$shards, mgcvST:::.mgcvst_read_shard))
    expect_identical(sum(z$rows), nrow(out))
    out[order(out$i, out$j), ]
  }
  for (index in list(NULL, cbind(pairs[, 1L], pairs[, 2L]))) {
    seen <- list()
    full <- run(index, 1e5)
    seen <- list()
    bounded <- run(index, 4200)
    expect_identical(nrow(bounded), nrow(pairs))
    expect_equal(bounded, full, tolerance = 1e-12)
    expect_identical(unname(as.matrix(bounded[, c("i", "j")])), unname(pairs))
    expect_gt(length(seen), 1L)
    expect_true(all(vapply(seen, `[[`, numeric(1L), "features") <= 2L))
    expect_true(all(full$status == 0L))
    expect_true(all(is.finite(full$log_p_two_sided)))
  }
  small <- run(NULL, 1e5, chunk_size = 2L)
  expect_equal(small, run(NULL, 1e5), tolerance = 1e-12)
  expect_true(all(vapply(seen, `[[`, numeric(1L), "rows") <= 9L))
})

test_that("a pair of a failed state is returned with status 3, never a p-value", {
  .pair_pipeline_mocks()
  testthat::local_mocked_bindings(
    .mgcvst_pair_build_batch = function(fit, ids, threads, native) {
      lapply(ids, function(i) {
        if (i == 2L) list(error = "state failed") else
          list(a = c(i, i + 0.25), M = diag(c(i + 0.5, i + 1)), width = 2L)
      })
    }, .package = "mgcvST"
  )
  z <- mgcvST:::.mgcvst_pair_pipeline(.pair_pipeline_fit(3L), NULL, 1L, 10L, FALSE)
  out <- mgcvST:::.mgcvst_read_shard(z$shards[1L])
  bad <- out$i == 2L | out$j == 2L
  expect_true(all(out$status[bad] == 3L))
  expect_true(all(is.na(out$score[bad]) & is.na(out$log_p_two_sided[bad])))
  expect_identical(out$status[!bad], 0L)
  expect_identical(z$failed, data.frame(feature_id = "g2", error = "state failed"))
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
  model <- mgcvST:::.mgcvst_pair_build_batch(fit, 1:2, 2L, native)
  expect_equal(seen$scale, c(0.5, 0.5))
  expect_identical(model[[1L]]$a, c(1, 2))
  expect_equal(model[[2L]]$M, diag(c(2, 3)))
  expect_identical(model[[1L]]$width, c(global = 2L))
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
