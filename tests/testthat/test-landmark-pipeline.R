test_that("score-selected CUR constructs units once and diagnoses held-out pairs", {
  M <- list(diag(c(1, 2)), matrix(c(2, .3, .3, 1), 2),
            diag(c(.8, 1.5)), matrix(c(1.2, -.2, -.2, .9), 2),
            diag(c(2.1, .7)), diag(c(.6, .8)))
  G <- length(M)
  fit <- list(feature_id = paste0("g", seq_len(G)), estimator = "INLA",
    test_engine = "single_model", score_backend = "sparse",
    working_error = matrix(0, 4L, G), working_variance = matrix(1, 4L, G),
    dispersion = rep(1, G), smoothing_parameters = matrix(1, G, 1L),
    geometry = list(nuisance_design = matrix(numeric(), 4L, 0L)),
    score_sparse = list(Q = Matrix::Diagonal(3L)))
  basis <- list(coordinate = rbind(diag(2L), c(0, 0)), rank = 2L)
  built <- materialized <- integer()
  testthat::local_mocked_bindings(
    .inlast_sparse_prepare = function(fit) fit,
    .inlast_sparse_units = function(fit, features, threads = 1L) {
      built <<- c(built, features)
      lapply(features, function(i) list(id = i, a = c(i + 1, i + .25, 0)))
    },
    .inlast_sparse_materialize_reduced = function(fit, units, basis, threads = 1L) {
      lapply(units, function(z) {
        materialized <<- c(materialized, z$id)
        list(a = z$a[1:2], M = M[[z$id]], width = 2L)
      })
    }, .package = "mgcvST")
  path <- tempfile("mgcvst-cur-score-")
  on.exit(unlink(path, recursive = TRUE), add = TRUE)
  pairs <- t(utils::combn(seq_len(G), 2L))
  run <- function(chunk_size = 4L) mgcvST:::.mgcvst_pair_approximate(
    fit, pairs, seq_len(nrow(pairs)), 2L, chunk_size, FALSE, basis,
    n_ref = 2L, ref_method = "score", ref_seed = 27L,
    cache_bytes = 1000000, checkpoint_dir = path,
    diagnostic_pairs = 10000L)
  z <- run()
  expect_equal(sort(built), seq_len(G))
  expect_equal(sort(materialized), seq_len(G))
  expect_equal(z$metadata$unit_builds, G)
  expect_equal(z$metadata$builds, G)
  expect_equal(length(z$metadata$diagnostics$pair_index), 6L)
  expect_equal(z$metadata$diagnostics$error$pairs, rep(6L, 5L))
  expect_equal(z$metadata$diagnostics$exact_builds, 0L)
  for (k in seq_len(nrow(pairs))) {
    i <- pairs[k, 1L]; j <- pairs[k, 2L]
    U <- sum(c(i + 1, i + .25) * c(j + 1, j + .25))
    expect_equal(z$result$score[k], U, tolerance = 1e-10)
    expect_true(is.finite(z$result$p_value[k]) &&
                z$result$p_value[k] >= 0 && z$result$p_value[k] <= 1)
  }
  resumed <- run(3L)
  expect_identical(resumed$metadata$builds, 0L)
  expect_identical(resumed$metadata$unit_builds, 0L)
  expect_identical(resumed$metadata$resumed_pairs, nrow(pairs))
  expect_equal(sort(built), seq_len(G))
  expect_equal(sort(materialized), seq_len(G))
  expect_equal(resumed$result, z$result)
})
