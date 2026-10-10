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

test_that("mgcvST.test() rejects sparse INLA fits for every pair universe", {
  fit <- .inla_openmp_fit()
  pairs <- rbind(c(1L, 2L), c(1L, 3L), c(1L, 4L))
  message <- "does not accept inlaST.estimate() fits; use inlaST.test()."
  expect_error(mgcvST.test(fit, pairs = pairs, moments = "exact"), message, fixed = TRUE)
  expect_error(mgcvST.test(fit, moments = "exact"), message, fixed = TRUE)
  expect_error(mgcvST.test(fit, pairs = pairs, adjust = "BH", threads = 2L, moments = "exact"),
               message, fixed = TRUE)
  expect_false(any(c("BPPARAM", "calibration", "...") %in% names(formals(mgcvST.test))))
  expect_false(any(c("BPPARAM", "calibration", "approximate_test", "...") %in%
                     names(formals(inlaST.test))))
})
