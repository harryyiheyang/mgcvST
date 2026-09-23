test_that("holdout sampling excludes reference endpoints and preserves RNG", {
  index <- rbind(c(1L, 2L), c(2L, 3L), c(3L, 4L),
                 c(4L, 5L), c(5L, 3L), c(3L, 4L))
  set.seed(19)
  old <- .Random.seed
  rows <- mgcvST:::.mgcvst_landmark_holdout(index, 2L, 3L, 0L)
  expect_identical(.Random.seed, old)
  expect_identical(rows,
                   mgcvST:::.mgcvst_landmark_holdout(index, 2L, 3L, 0L))
  expect_length(rows, 3L)
  expect_false(any(index[rows, ] == 2L))
  expect_identical(mgcvST:::.mgcvst_landmark_holdout(index, 2L, 0L, 0L),
                   integer())
  expect_identical(mgcvST:::.mgcvst_landmark_holdout(
    matrix(integer(), 0L, 2L), 2L, 3L, 0L), integer())
})

test_that("shared Liu log tails keep central underflow information", {
  liu <- mgcvST:::.liu_squared_score_moments(400, 21, 273, 4161, 65793)
  expect_identical(liu$p_value, 0)
  expect_lt(mgcvST:::.mgcvst_liu_logp(liu), -308)
  extreme <- mgcvST:::.liu_squared_score_moments(2500, 21, 273, 4161, 65793)
  expect_identical(extreme$p_value, 0)
  expect_true(is.finite(mgcvST:::.mgcvst_liu_logp(extreme)))
  expect_identical(mgcvST:::.mgcvst_liu_logp(
    list(transformed = 2000, df = 0.13, ncp = 1)), -Inf)
  noncentral <- list(transformed = 5, df = 3, ncp = 2)
  expect_equal(mgcvST:::.mgcvst_liu_logp(noncentral),
               stats::pchisq(5, df = 3, ncp = 2,
                             lower.tail = FALSE, log.p = TRUE) / log(10))
  expect_true(is.na(mgcvST:::.mgcvst_liu_logp(
    list(transformed = NaN, df = 3, ncp = 0))))
  expect_identical(mgcvST:::.mgcvst_liu_logp(
    list(transformed = numeric(), df = numeric(), ncp = numeric())), numeric())
})

test_that("exact holdout reads only double shards with bounded resident states", {
  path <- tempfile("mgcvst-landmark-exact-")
  on.exit(unlink(path, recursive = TRUE), add = TRUE)
  store <- mgcvST:::.mgcvst_store_open(
    path, signature = list(test = "landmark"),
    feature_ids = paste0("g", seq_len(6L)), storage = "double"
  )
  M <- list(
    diag(c(1, 2, 3)),
    crossprod(matrix(c(1, 1, 0, 0, 2, 1, 1, 0, 2), 3L)) + diag(3) / 2,
    tcrossprod(matrix(c(1, 0, 2, 1, 2, 1, 0, 1, 1), 3L)) + diag(3) / 3,
    diag(c(2, 1, 4)),
    matrix(0, 3L, 3L)
  )
  a <- list(c(1, 2, 3), c(2, -1, 1), c(-1, 1, 2),
            c(20, 0, 0), c(0, 0, 0))
  for (i in seq_along(M)) {
    mgcvST:::.mgcvst_store_write(store, i, list(a = a[[i]], M = M[[i]]))
  }
  mgcvST:::.mgcvst_store_write(store, 6L, list(error = "fit failed"))
  state_bytes <- max(vapply(seq_len(4L), function(i)
    as.numeric(object.size(mgcvST:::.mgcvst_store_read(store, i))), numeric(1L)))
  index <- rbind(c(4L, 4L), c(3L, 1L), c(2L, 4L), c(1L, 2L),
                 c(4L, 3L), c(2L, 3L), c(1L, 4L), c(1L, 2L),
                 c(5L, 1L), c(6L, 1L))

  testthat::local_mocked_bindings(
    .mgcvst_pair_build_batch = function(...) stop("score builder was called"),
    .mgcvst_model_score_state = function(...) stop("score builder was called"),
    .package = "mgcvST"
  )
  small <- mgcvST:::.mgcvst_landmark_exact(
    store, index, threads = 2L, chunk_size = 100L,
    cache_bytes = 2.1 * state_bytes
  )
  large <- mgcvST:::.mgcvst_landmark_exact(
    store, index, threads = 1L, chunk_size = 100L,
    cache_bytes = 10 * state_bytes
  )
  expect_equal(small$moments, large$moments)
  expect_equal(small$score, large$score)
  expect_equal(small$log10_p, large$log10_p)
  expect_identical(small$builds, 0L)
  expect_identical(large$builds, 0L)
  expect_gt(small$chunks, 1L)
  expect_lte(small$cache_peak_bytes, 2.1 * state_bytes)
  expect_match(small$error_message[9L], "non-finite or non-positive")
  expect_match(small$error_message[10L], "fit failed")
  expect_true(is.finite(small$log10_p[1L]))
  expect_lt(small$log10_p[1L], -308)
  expect_identical(small$p_value[1L], 0)

  good <- seq_len(8L)
  expected <- mgcvST:::mgcvst_pair_trace_powers_cpp(
    M[seq_len(4L)], index[good, , drop = FALSE], 4L, 1L
  )
  expect_equal(small$moments[good, ], expected, tolerance = 1e-12)
  expect_equal(small$score[good], vapply(good, function(k)
    sum(a[[index[k, 1L]]] * a[[index[k, 2L]]]), numeric(1L)))
  expect_equal(small$p_value[good], vapply(good, function(k)
    mgcvST:::.liu_squared_score_moments(abs(small$score[k]),
      small$moments[k, 1L], small$moments[k, 2L],
      small$moments[k, 3L], small$moments[k, 4L])$p_value, numeric(1L)))
  expect_error(mgcvST:::.mgcvst_landmark_exact(
    store, index[2L, , drop = FALSE], 1L, 1L, cache_bytes = state_bytes),
    "cannot hold")
  empty <- mgcvST:::.mgcvst_landmark_exact(
    store, matrix(integer(), 0L, 2L), 1L, 1L, state_bytes
  )
  expect_identical(dim(empty$moments), c(0L, 4L))
  expect_identical(empty$builds, 0L)
  expect_identical(empty$chunks, 0L)
})

test_that("unrepresentable log tail does not invalidate a valid zero p-value", {
  path <- tempfile("mgcvst-landmark-log-boundary-")
  on.exit(unlink(path, recursive = TRUE), add = TRUE)
  store <- mgcvST:::.mgcvst_store_open(
    path, signature = list(test = "log-boundary"), feature_ids = "g"
  )
  mgcvST:::.mgcvst_store_write(store, 1L,
    list(a = c(50, 0, 0), M = diag(c(2, 1, 4))))
  testthat::local_mocked_bindings(
    .mgcvst_liu_logp = function(liu) rep(-Inf, length(liu$p_value)),
    .package = "mgcvST"
  )
  exact <- mgcvST:::.mgcvst_landmark_exact(
    store, matrix(c(1L, 1L), 1L), 1L, 1L, 1e6
  )
  expect_identical(exact$p_value, 0)
  expect_identical(exact$log10_p, -Inf)
  expect_true(is.na(exact$error_message))
})
