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
  units <- function(fit, features, threads = 1L) {
    calls$features[[length(calls$features) + 1L]] <- features
    lapply(features, function(i) list(
      feature = i, a = c(i, i + 0.5), error = NULL
    ))
  }
  materialize <- function(fit, units, threads = 1L) {
    lapply(units, function(z) list(
      a = z$a,
      M = diag(c(0.6 + z$feature / 10, 1.1 + z$feature / 20)),
      error = NULL
    ))
  }
  testthat::local_mocked_bindings(
    .inlast_sparse_prepare = function(fit) fit,
    .inlast_sparse_units = units,
    .inlast_sparse_materialize = materialize,
    .package = "mgcvST")
  pairs <- rbind(c(1L, 2L), c(3L, 4L), c(1L, 4L))
  z <- mgcvST:::.mgcvst_inla_test_pairs(
    fit, pairs, seq_len(nrow(pairs)), threads = 2L,
    chunk_size = 1L, verbose = FALSE
  )$result

  expect_identical(sort(unlist(calls$features)), seq_len(4L))
  expect_true(all(table(unlist(calls$features)) == 1L))
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

  expect_error(
    mgcvST.test(fit, pairs = pairs, calibration = "davies"),
    "calibration = 'liu' only"
  )
  snow <- BiocParallel::SnowParam(2L, type = "SOCK", progressbar = FALSE)
  expect_error(
    mgcvST.test(fit, pairs = pairs, calibration = "liu", BPPARAM = snow),
    "must be SerialParam"
  )
})

test_that("sparse INLA WGCNA scores retain redundant coordinates and q-1 scale", {
  fit <- .inla_openmp_fit(3L)
  seen_threads <- NULL
  batch <- function(fit, features, threads = 1L, score_only = FALSE,
                    null_target = FALSE) {
    seen_threads <<- threads
    lapply(features, function(i) list(
      a = c(i, i + 1, -2 * i - 1), M = NULL, error = NULL,
      normalization = 2L, width = c(global = 3L), backend = "test"
    ))
  }
  testthat::local_mocked_bindings(
    .inlast_sparse_prepare = function(fit) fit,
    .inlast_sparse_batch = batch,
    .package = "mgcvST")
  z <- mgcvST:::.mgcvst_inla_wgcna_scores(
    fit, c(3L, 1L, 2L), threads = 3L, verbose = FALSE
  )
  expect_identical(seen_threads, 3L)
  expect_identical(z$width, c(global = 3L))
  expect_identical(z$normalization, 2L)
  expect_identical(z$feature_id, c("g3", "g1", "g2"))
  expect_equal(crossprod(z$A) / z$normalization,
    crossprod(z$A) / 2, tolerance = 0)
})
