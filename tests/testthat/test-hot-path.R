test_that("mandatory marginal preserves existing fitting outputs", {
  old <- old_st()
  f <- st_fixture()
  marginal_callback <- function(...) list(smooth.pvalue = .375)
  for (G in list(f$model)) {
    a <- old$mgcvST.estimate(f$Y, G, marginal_test = marginal_callback, retain_smooth = TRUE, spatial = "all")
    b <- mgcvST.estimate(f$Y, G, diagnostics = TRUE,
                         marginal_test = marginal_callback, retain_smooth = TRUE, spatial = "all")
    expect_identical(strip_elapsed(b), strip_elapsed(a))
    fast <- mgcvST.estimate(f$Y, G, diagnostics = FALSE,
                            marginal_test = marginal_callback, spatial = "all")
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
    fit <- mgcvST.estimate(f$Y, G, diagnostics = FALSE, spatial = "all")
    expect_true(all(is.finite(fit$working_error)))
    if (!is.null(fit$diagnostics$wood_p_value)) expect_true(all(is.na(fit$diagnostics$wood_p_value)))
  }
})

test_that("model batches reuse the prepared GAM design without lpmatrix", {
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
    diagnostics = FALSE, spatial = "all"
  )
  expect_identical(getOption("mgcvST.test_lpmatrix_count"), 0L)
  expect_true(all(is.finite(fit$working_error)))
})

test_that("prepared GAM designs do not depend on later prediction methods", {
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
    diagnostics = FALSE, spatial = "all"
  )
  expect_identical(getOption("mgcvST.test_lpmatrix_count"), 0L)
  expect_true(all(is.finite(fit$working_error)))
  expect_true(all(vapply(fit$nuisance_covariance, is.matrix, logical(1L))))
})

test_that("custom worker initialization retains shared prediction geometry", {
  f <- st_fixture(nuisance = TRUE)
  fit <- mgcvST.estimate(
    f$Y, f$model, BPPARAM = BiocParallel::SerialParam(), chunk_size = 1L,
    worker_init = function() invisible(NULL), diagnostics = FALSE, spatial = "all"
  )
  expect_true(all(vapply(fit$nuisance_covariance, is.matrix, logical(1L))))
  expect_true(all(is.finite(fit$working_error)))
})

test_that("pair universes are validated and duplicated tests are rejected", {
  f <- st_fixture(nuisance = TRUE)
  fit <- mgcvST.estimate(f$Y, f$model, spatial = "all")
  pairs <- rbind(c(1L, 2L), c(3L, 1L), c(2L, 3L))
  expect_error(mgcvST.test(fit, pairs = rbind(pairs, c(2L, 1L)), moments = "exact"),
               "duplicated tests")
  expect_error(mgcvST.test(fit, pairs = rbind(c(1L, 1L)), moments = "exact"), "two different features")
  expect_error(mgcvST.test(fit, pairs = matrix(c("a", "b"), 1L), moments = "exact"), "unknown feature IDs")
  for (chunk in c(1L, 100L)) {
    out <- mgcvST.test(fit, pairs = pairs, chunk_size = chunk, moments = "exact")
    expect_identical(nrow(out$results), 3L)
    expect_true(all(out$results$i < out$results$j))
    expect_true(all(is.finite(out$results$log_p_two_sided)))
  }
  fit$smoothing_parameters[2, 1] <- -1
  bad <- mgcvST.test(fit, pairs = pairs, moments = "exact")
  expect_true(all(bad$results$status[bad$results$i == 2L | bad$results$j == 2L] == 3L))
  expect_identical(bad$failed$feature_id, "response2")
})

test_that("conditional nuisance state is compact and shared", {
  f <- st_fixture(nuisance = TRUE)
  fit <- mgcvST.estimate(
    f$Y, f$model, diagnostics = FALSE,
    BPPARAM = BiocParallel::SerialParam(), spatial = "all"
  )
  LN <- fit$geometry$nuisance_design
  blocks <- fit$nuisance_covariance
  expect_true(is.matrix(LN))
  expect_identical(length(blocks), nrow(f$Y))
  expect_true(all(vapply(blocks, is.matrix, logical(1L))))
  expect_true(all(vapply(blocks, function(x) identical(dim(x), rep(ncol(LN), 2L)), logical(1L))))
  expect_false(any(vapply(blocks[-1L], identical, logical(1L), y = blocks[[1L]])))
  expect_null(fit$gam)
})

test_that("ordinary overall low-rank smooths share the same Vp machinery", {
  f <- st_fixture()
  data <- f$data
  data$w <- sin(seq_len(nrow(data)) / 5)
  model <- model.set(
    response ~ offset(offset0) + s(z, k = 6) + s(w, k = 6),
    data, f$basis, family = mgcv::nb()
  )
  fit <- mgcvST.estimate(
    f$Y, model, diagnostics = FALSE,
    BPPARAM = BiocParallel::SerialParam(), spatial = "all"
  )
  expect_identical(fit$geometry$nuisance_projection, "conditional_Vp_block")
  expect_true(all(vapply(fit$nuisance_covariance, is.matrix, logical(1L))))
})

test_that("new switches reject non-logical values", {
  expect_false("marginal" %in% names(formals(mgcvST.estimate)))
  expect_error(mgcvST.estimate(NULL, NULL, marginal_test = NULL,
                               marginal_args = list(), marginal = NA, spatial = "all"), "always runs")
  expect_error(mgcvST.estimate(NULL, NULL, diagnostics = 0, spatial = "all"), "diagnostics must")
  expect_error(mgcvST.estimate(NULL, NULL, retain_marginal = 1, spatial = "all"), "retain_marginal must")
})

test_that("a raw gam setup and the prepared model give identical estimates and tests", {
  f <- st_fixture()
  sp <- BiocParallel::SerialParam()
  from_G <- suppressWarnings(mgcvST.estimate(f$Y, f$G, BPPARAM = sp, spatial = "all"))
  from_model <- suppressWarnings(mgcvST.estimate(f$Y, f$model, BPPARAM = sp, spatial = "all"))
  # Only the stored setup object differs; every estimate is the same.
  same <- function(x) {
    x <- strip_elapsed(x)
    x$model <- NULL
    x$signature <- x$estimation_context <- NULL
    x
  }
  expect_identical(same(from_G), same(from_model))
  expect_identical(mgcvST.test(from_G, moments = "exact")$results, mgcvST.test(from_model, moments = "exact")$results)
})
