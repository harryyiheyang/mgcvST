.approx_fixture <- function() {
  u <- 1 / sqrt(2)
  v <- c(1, 0)
  w <- c(0, 1)
  plus <- c(u, u)
  minus <- c(u, -u)
  M <- list(diag(2L), tcrossprod(v), tcrossprod(w),
            tcrossprod(plus), tcrossprod(minus))
  fit <- list(
    feature_id = paste0("g", seq_along(M)),
    test_engine = "single_model", estimator = "INLA",
    score_backend = "sparse",
    working_error = matrix(0, 4L, length(M)),
    working_variance = matrix(1, 4L, length(M)),
    dispersion = rep(1, length(M)), lambda = rep(1, length(M)),
    smoothing_parameters = matrix(1, length(M), 1L),
    geometry = list(nuisance_design = matrix(numeric(), 4L, 0L)),
    score_sparse = list(A = Matrix::Diagonal(4L, 3L),
                        Q = Matrix::Diagonal(3L),
                        constraint = rep(1, 3L), sp_index = 1L)
  )
  basis <- list(coordinate = matrix(c(1, 0, 0, 0, 1, 0), 3L, 2L),
                basis = matrix(c(1, 0, 0, 0, 1, 0), 3L, 2L),
                rank = 2L, coverage = 0.995)
  list(fit = fit, basis = basis, M = M)
}

test_that("signed reference reconstruction handles noncommuting PSD matrices", {
  f <- .approx_fixture()
  calls <- new.env(parent = emptyenv())
  calls$built <- integer()
  calls$trace <- 0L
  original_trace <- mgcvST:::mgcvst_pair_trace_powers_cpp
  testthat::local_mocked_bindings(
    .inlast_sparse_prepare = function(fit) fit,
    .inlast_sparse_units = function(fit, features, threads = 1L) {
      calls$built <- c(calls$built, features)
      lapply(features, function(i) list(feature = i))
    },
    .inlast_sparse_materialize_reduced = function(fit, units, basis,
                                                   threads = 1L) {
      lapply(units, function(z) list(a = c(z$feature, z$feature + 0.25),
                                     M = f$M[[z$feature]], width = 2L))
    },
    mgcvst_pair_trace_powers_cpp = function(...) {
      calls$trace <- calls$trace + 1L
      original_trace(...)
    },
    .package = "mgcvST"
  )
  pairs <- t(utils::combn(seq_along(f$M), 2L))
  path <- tempfile("mgcvst-approx-")
  on.exit(unlink(path, recursive = TRUE), add = TRUE)
  run <- function() mgcvST:::.mgcvst_pair_approximate(
    f$fit, pairs, seq_len(nrow(pairs)), threads = 2L,
    chunk_size = 4L, verbose = FALSE, basis = f$basis,
    n_ref = 5L, ref_method = "random", ref_seed = 42L,
    cache_bytes = 100000L, checkpoint_dir = path, ref_tol = 1e-12
  )
  first <- run()
  expect_identical(sort(calls$built), seq_along(f$M))
  expect_identical(first$metadata$builds, 5L)
  expect_gt(first$metadata$reference_negative[2L], 0L)
  expect_gt(first$metadata$reference_negative[3L], 0L)
  expect_identical(first$result$pair_index, seq_len(nrow(pairs)))
  for (k in seq_len(nrow(pairs))) {
    i <- pairs[k, 1L]
    j <- pairs[k, 2L]
    score <- sum(c(i, i + 0.25) * c(j, j + 0.25))
    expected <- mgcvST:::rkhs_score_calibrate(
      score, f$M[[i]], f$M[[j]], method = "liu"
    )
    expect_equal(first$result$score[k], score, tolerance = 1e-12)
    if (expected$information <= 1e-10) {
      expect_true(is.na(first$result$information[k]))
      expect_true(is.na(first$result$p_value[k]))
    } else {
      expect_equal(first$result$information[k], expected$information,
                   tolerance = 1e-8)
      expect_equal(first$result$p_value[k], expected$p_two_sided,
                   tolerance = 1e-7)
    }
  }
  trace_before_resume <- calls$trace
  resumed <- run()
  expect_identical(resumed$metadata$builds, 0L)
  expect_identical(sort(calls$built), seq_along(f$M))
  expect_identical(calls$trace, trace_before_resume)
  expect_equal(resumed$result, first$result, tolerance = 1e-12)

  partial_path <- tempfile("mgcvst-approx-partial-")
  on.exit(unlink(partial_path, recursive = TRUE), add = TRUE)
  partial <- function(directory) mgcvST:::.mgcvst_pair_approximate(
    f$fit, pairs, seq_len(nrow(pairs)), threads = 2L,
    chunk_size = 3L, verbose = FALSE, basis = f$basis,
    n_ref = 3L, ref_method = "random", ref_seed = 0L,
    cache_bytes = 100000L, checkpoint_dir = directory
  )
  smaller <- partial(partial_path)
  expect_identical(smaller$metadata$builds, 5L)
  expect_identical(smaller$metadata$reference_count, 3L)
  expect_identical(smaller$result$pair_index, seq_len(nrow(pairs)))
  expect_true(all(is.finite(smaller$result$p_value) |
                    !is.na(smaller$result$error_message)))
  trace_before_resume <- calls$trace
  smaller_resumed <- partial(partial_path)
  expect_identical(smaller_resumed$metadata$builds, 0L)
  expect_identical(calls$trace, trace_before_resume)
  expect_equal(smaller_resumed$result, smaller$result, tolerance = 1e-12)

  memory_only <- partial(NULL)
  expect_identical(memory_only$metadata$builds, 5L)
  expect_null(memory_only$metadata$path)
})

test_that("reference selection preserves the caller RNG", {
  f <- .approx_fixture()
  set.seed(123L)
  before <- .Random.seed
  used <- seq_along(f$M)
  for (method in "random") {
    selected <- mgcvST:::.mgcvst_landmark_select(
      f$fit, used, 3L, method, 7L
    )
    expect_length(selected, 3L)
    expect_true(all(selected %in% used))
    expect_identical(.Random.seed, before)
  }
})

test_that("approximate preparation stops when detected working memory is exhausted", {
  f <- .approx_fixture()
  testthat::local_mocked_bindings(
    .inlast_sparse_prepare = function(fit) fit,
    .mgcvst_memory_probe = function(...) list(available = 0),
    .package = "mgcvST"
  )
  expect_error(mgcvST:::.mgcvst_pair_approximate(
    f$fit, matrix(c(1L, 2L), nrow = 1L), 1L, threads = 1L,
    chunk_size = 1L, verbose = FALSE, basis = f$basis,
    n_ref = 1L, cache_bytes = 100000L
  ), "Insufficient memory for one approximate score-state batch")
})
