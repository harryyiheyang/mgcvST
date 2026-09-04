test_that("mandatory marginal preserves existing fitting outputs", {
  old <- old_st()
  f <- st_fixture()
  marginal_callback <- function(...) list(smooth.pvalue = .375)
  for (G in list(f$G, f$model)) {
    a <- old$mgcvST.estimate(f$Y, G, marginal_test = marginal_callback, retain_smooth = TRUE)
    b <- mgcvST.estimate(f$Y, G, diagnostics = TRUE,
                         marginal_test = marginal_callback, retain_smooth = TRUE)
    expect_identical(strip_elapsed(b), strip_elapsed(a))
    fast <- mgcvST.estimate(f$Y, G, diagnostics = FALSE,
                            marginal_test = marginal_callback)
    expect_identical(fast$working_error, a$working_error)
    expect_identical(fast$working_variance, a$working_variance)
    expect_identical(fast$lambda, a$lambda)
    expect_true(all(fast$diagnostics$marginal_p_value == .375))
    expect_null(fast$marginal_data)
  }
})

test_that("diagnostics FALSE does not invoke summary.gam", {
  f <- st_fixture()
  trace("summary.gam", where = asNamespace("mgcv"),
        tracer = quote(stop("summary was called")), print = FALSE)
  on.exit(untrace("summary.gam", where = asNamespace("mgcv")), add = TRUE)
  for (G in list(f$G, f$model)) {
    fit <- mgcvST.estimate(f$Y, G, diagnostics = FALSE)
    expect_true(all(is.finite(fit$working_error)))
    if (!is.null(fit$diagnostics$wood_p_value)) expect_true(all(is.na(fit$diagnostics$wood_p_value)))
  }
})

test_that("model batches reuse one exact formal training lpmatrix", {
  f <- st_fixture(nuisance = TRUE)
  original <- mgcvST:::.gam_training_lpmatrix
  old_options <- options(mgcvST.test_lpmatrix = original,
                         mgcvST.test_lpmatrix_count = 0L)
  on.exit(options(old_options), add = TRUE)
  testthat::local_mocked_bindings(
    .gam_training_lpmatrix = function(fit) {
      options(mgcvST.test_lpmatrix_count =
                getOption("mgcvST.test_lpmatrix_count") + 1L)
      getOption("mgcvST.test_lpmatrix")(fit)
    },
    .package = "mgcvST"
  )
  fit <- mgcvST.estimate(
    f$Y, f$model, BPPARAM = BiocParallel::SerialParam(), chunk_size = 1L,
    diagnostics = FALSE
  )
  expect_identical(getOption("mgcvST.test_lpmatrix_count"), 1L)
  expect_true(all(is.finite(fit$working_error)))
})

test_that("unknown prediction methods disable exact-lpmatrix reuse", {
  f <- st_fixture(nuisance = TRUE)
  original_lpmatrix <- mgcvST:::.gam_training_lpmatrix
  original_predictor <- mgcvST:::Predict.matrix.spde.smooth
  old_options <- options(
    mgcvST.test_lpmatrix = original_lpmatrix,
    mgcvST.test_lpmatrix_count = 0L,
    mgcvST.test_predictor = original_predictor
  )
  on.exit(options(old_options), add = TRUE)
  testthat::local_mocked_bindings(
    Predict.matrix.spde.smooth = function(object, data) {
      getOption("mgcvST.test_predictor")(object, data)
    },
    .gam_training_lpmatrix = function(fit) {
      options(mgcvST.test_lpmatrix_count =
                getOption("mgcvST.test_lpmatrix_count") + 1L)
      getOption("mgcvST.test_lpmatrix")(fit)
    },
    .package = "mgcvST"
  )
  fit <- mgcvST.estimate(
    f$Y, f$model, BPPARAM = BiocParallel::SerialParam(), chunk_size = 1L,
    diagnostics = FALSE
  )
  expect_identical(getOption("mgcvST.test_lpmatrix_count"), nrow(f$Y))
  expect_true(all(is.finite(fit$working_error)))
  expect_true(all(vapply(fit$nuisance_covariance, is.null, logical(1L))))
  expect_s3_class(mgcvST:::.mgcvst_model_operator(fit, 1L)$operator,
                  "rkhs_score_operator")
})

test_that("custom worker initialization disables shared prediction geometry", {
  f <- st_fixture(nuisance = TRUE)
  fit <- mgcvST.estimate(
    f$Y, f$model, BPPARAM = BiocParallel::SerialParam(), chunk_size = 1L,
    worker_init = function() invisible(NULL),
    diagnostics = FALSE
  )
  expect_true(all(vapply(fit$nuisance_covariance, is.null, logical(1L))))
  expect_true(all(is.finite(fit$working_error)))
})

test_that("Vp model projection is numerically equivalent with unchanged contract", {
  old <- old_st()
  f <- st_fixture(nuisance = TRUE)
  fit <- mgcvST.estimate(f$Y, f$model)
  pairs <- rbind(c(1L,2L), c(3L,1L), c(2L,3L))
  duplicated <- rbind(pairs, c(2L,1L))
  expect_error(old$mgcvST.test(fit,pairs=duplicated), "duplicated tests")
  expect_error(mgcvST.test(fit,pairs=duplicated), "duplicated tests")
  for (cal in c("liu", "davies")) {
    if (cal == "davies") skip_if_not_installed("CompQuadForm")
    for (chunk in c(1L, 100L)) {
      a <- old$mgcvST.test(fit, pairs = pairs, highlight = matrix(c(1L,3L),1), calibration = cal, chunk_size = chunk)
      b <- mgcvST.test(fit, pairs = pairs, highlight = matrix(c(1L,3L),1), calibration = cal, chunk_size = chunk)
      expect_identical(a$call,b$call)
      expect_numerically_equivalent_test(a, b)
    }
  }
  for (i in 1:3) {
    a <- old$.mgcvst_model_score_state(fit, i)
    b <- mgcvST:::.mgcvst_model_score_state(fit, i)
    expect_equal(a$a, b$a, tolerance = 1e-10)
    expect_equal(a$M, b$M, tolerance = 1e-10)
    expect_identical(a$width, b$width)
  }
  fit$smoothing_parameters[2,1] <- -1
  a <- old$mgcvST.test(fit, pairs = pairs)
  b <- mgcvST.test(fit, pairs = pairs)
  expect_numerically_equivalent_test(a, b)
})

test_that("conditional nuisance state is compact, shared and CppMatrix-backed", {
  f <- st_fixture(nuisance = TRUE)
  fit <- mgcvST.estimate(f$Y, f$model, diagnostics = FALSE)
  LN <- fit$geometry$nuisance_design
  blocks <- fit$nuisance_covariance
  expect_true(is.matrix(LN))
  expect_identical(length(blocks), nrow(f$Y))
  expect_true(all(vapply(blocks, is.matrix, logical(1L))))
  expect_true(all(vapply(blocks, function(x) identical(dim(x), rep(ncol(LN), 2L)), logical(1L))))
  expect_false(any(vapply(blocks[-1L], identical, logical(1L), y = blocks[[1L]])))
  expect_null(fit$gam)
  original_multiply <- CppMatrix::matrixMultiply
  calls <- 0L
  testthat::local_mocked_bindings(
    matrixMultiply = function(...) {
      calls <<- calls + 1L
      original_multiply(...)
    },
    .package = "CppMatrix"
  )
  state <- mgcvST:::.mgcvst_model_score_state(fit, 1L)
  expect_gt(calls, 0L)
  expect_true(all(is.finite(c(state$a, state$M))))
  expect_false(any(c("%*%", "crossprod", "tcrossprod") %in%
                   all.names(body(mgcvST:::.mgcvst_model_apply_P))))
})

test_that("ordinary overall low-rank smooths share the same Vp machinery", {
  f <- st_fixture()
  model <- model.set(
    response ~ offset(offset0) + z + s(x, k = 6) + s(y, k = 6),
    f$data, f$basis, family = mgcv::nb()
  )
  fit <- mgcvST.estimate(f$Y, model, diagnostics = FALSE)
  expect_identical(fit$geometry$nuisance_projection, "conditional_Vp_block")
  expect_true(all(vapply(fit$nuisance_covariance, is.matrix, logical(1L))))
  testthat::local_mocked_bindings(
    matrixEigen = function(...) stop("nuisance penalty was eigendecomposed"),
    .package = "CppMatrix"
  )
  state <- mgcvST:::.mgcvst_model_score_state(fit, 1L)
  expect_true(all(is.finite(c(state$a, state$M))))
})

test_that("ordinary Liu engine is unchanged", {
  old <- old_st()
  for (pc in c(FALSE, TRUE)) {
    f <- st_fixture(pc = pc)
    fit <- mgcvST.estimate(f$Y, f$G)
    pairs <- t(combn(1:3,2))
    expect_identical(strip_elapsed(mgcvST.test(fit,pairs=pairs)),
                     strip_elapsed(old$mgcvST.test(fit,pairs=pairs)))
  }
})

test_that("new switches reject non-logical values", {
  expect_false("marginal" %in% names(formals(mgcvST.estimate)))
  expect_error(mgcvST.estimate(NULL, NULL, marginal_test = NULL,
                               marginal_args = list(), marginal = NA), "always runs")
  expect_error(mgcvST.estimate(NULL, NULL, diagnostics = 0), "diagnostics must")
  expect_error(mgcvST.estimate(NULL, NULL, retain_marginal = 1), "retain_marginal must")
})
