.inla_openmp_fit <- function(p = 4L) {
  ids <- paste0("g", seq_len(p))
  structure(list(
    feature_id = ids,
    working_error = matrix(0, 5L, p, dimnames = list(NULL, ids)),
    working_variance = matrix(1, 5L, p, dimnames = list(NULL, ids)),
    dispersion = stats::setNames(rep(1, p), ids),
    lambda = stats::setNames(rep(1, p), ids),
    smoothing_parameters = matrix(1, p, 1L,
      dimnames = list(ids, "global")),
    diagnostics = data.frame(error_message = rep(NA_character_, p)),
    geometry = list(
      target = c(global = 1L), nuisance_design = matrix(numeric(), 5L, 0L),
      smooth = list(list(score_component = "global"))
    ),
    score_sparse = list(sp_index = 1L, Q = Matrix::Diagonal(10L)),
    score_backend = "sparse",
    test_engine = "single_model", score_components = "global",
    estimator = "INLA"
  ), class = c("inlaST_fit", "mgcvST_model_fit", "mgcvST_fit", "mgcvST"))
}

test_that("sparse INLA pair routing is bounded, OpenMP-only and exact Liu", {
  fit <- .inla_openmp_fit()
  calls <- new.env(parent = emptyenv())
  calls$features <- list()
  calls$materialized <- integer()
  units <- function(fit, features, threads = 1L) {
    calls$features[[length(calls$features) + 1L]] <- features
    lapply(features, function(i) list(
      feature = i, a = c(i, i + 0.5), error = NULL
    ))
  }
  materialize <- function(fit, units, basis, threads = 1L) {
    lapply(units, function(z) list(
      a = z$a,
      M = diag(c(0.6 + z$feature / 10, 1.1 + z$feature / 20)),
      error = NULL, padding = raw(1024L)
    ))
  }
  reduced <- function(fit, units, basis, threads = 1L) {
    calls$materialized <- c(calls$materialized, vapply(units, `[[`, numeric(1L), "feature"))
    materialize(fit, units, basis, threads)
  }
  basis <- function(fit, coverage = 0.995, full_rank = FALSE) {
    list(coordinate = diag(2L), basis = diag(2L), rank = 2L,
      coverage = coverage, kept = 1, tail = 0)
  }
  testthat::local_mocked_bindings(
    .inlast_sparse_prepare = function(fit) fit,
    .inlast_sparse_units = units,
    .inlast_sparse_materialize_reduced = reduced,
    .inlast_sparse_observation_basis = basis,
    .package = "mgcvST")
  pairs <- rbind(c(1L, 2L), c(1L, 3L), c(1L, 2L))
  out <- mgcvST:::.mgcvst_inla_test_pairs(
    fit, pairs, seq_len(nrow(pairs)), threads = 2L,
    chunk_size = 1L, verbose = FALSE
  )
  z <- out$result

  expect_identical(sort(unlist(calls$features)), seq_len(3L))
  expect_true(all(table(unlist(calls$features)) == 1L))
  expect_equal(sort(calls$materialized), as.numeric(seq_len(3L)))
  expect_equal(attr(z, "inla_pairwise")$cache_hits, 3L)
  expect_equal(mgcvST:::.mgcvst_inla_pair_chunk_size(fit), 128L)
  for (k in seq_len(nrow(pairs))) {
    i <- pairs[k, 1L]
    j <- pairs[k, 2L]
    a1 <- c(i, i + 0.5)
    a2 <- c(j, j + 0.5)
    M1 <- diag(c(0.6 + i / 10, 1.1 + i / 20))
    M2 <- diag(c(0.6 + j / 10, 1.1 + j / 20))
    expected <- mgcvST:::rkhs_score_calibrate(
      sum(a1 * a2), M1, M2, method = "liu"
    )
    expect_equal(z$signed_score[k], sum(a1 * a2), tolerance = 1e-12)
    expect_equal(z$information[k], expected$information, tolerance = 1e-12)
    expect_equal(z$p_two_sided[k], expected$p_two_sided, tolerance = 1e-12)
  }

  grouped <- mgcvST:::.mgcvst_inla_test_pairs(
    fit, pairs, seq_len(nrow(pairs)), threads = 2L,
    chunk_size = 2L, verbose = FALSE
  )$result
  expect_identical(grouped$pair_index, z$pair_index)
  expect_equal(grouped$signed_score, z$signed_score, tolerance = 1e-12)
  expect_equal(grouped$information, z$information, tolerance = 1e-12)
  expect_equal(grouped$p_two_sided, z$p_two_sided, tolerance = 1e-12)

  calls$basis <- 0L
  counted_basis <- function(fit, coverage = 0.995, full_rank = FALSE) {
    calls$basis <- calls$basis + 1L
    basis(fit, coverage, full_rank)
  }
  testthat::local_mocked_bindings(
    .inlast_sparse_observation_basis = counted_basis,
    .package = "mgcvST")
  # mgcvST.test() no longer accepts inlaST.estimate() fits (Task E1); the
  # basis-caching path it used to exercise is reached through inlaST.test()
  # with liu_approximation = "pca_learning" instead.
  public <- inlaST.test(
    fit, pairwise_method = "score_liu", liu_approximation = "pca_learning",
    pairs = pairs[1:2, , drop = FALSE], calibration = "liu",
    BPPARAM = BiocParallel::SerialParam()
  )
  expect_identical(calls$basis, 1L)
  expect_identical(public$timing$inla_projection$r, 2L)
  expect_identical(public$timing$inla_projection$unit_cache, "score_state_shards")

  evicted <- mgcvST:::.mgcvst_inla_test_pairs(
    fit, pairs, seq_len(nrow(pairs)), threads = 1L,
    chunk_size = 1L, verbose = FALSE, cache_bytes = 1
  )$result
  expect_equal(attr(evicted, "inla_pairwise")$cache_misses, 6L)
  expect_equal(attr(evicted, "inla_pairwise")$builds, 3L)

  expect_error(
    mgcvST.test(fit, pairs = pairs, calibration = "davies"),
    "mgcvST.test\\(\\) does not accept inlaST.estimate\\(\\) fits; use inlaST.test\\(\\)."
  )
  snow <- BiocParallel::SnowParam(2L, type = "SOCK", progressbar = FALSE)
  expect_error(
    mgcvST.test(fit, pairs = pairs, calibration = "liu", BPPARAM = snow),
    "mgcvST.test\\(\\) does not accept inlaST.estimate\\(\\) fits; use inlaST.test\\(\\)."
  )
})
