.spdest_inla_fixture <- function(n = 120L, seed = 2402L) {
  skip_if_not_installed("INLA")
  skip_if_not_installed("geometry")
  set.seed(seed)
  vertices <- as.matrix(expand.grid(x = seq(0, 1, length.out = 4L),
                                    y = seq(0, 1, length.out = 4L)))
  mesh <- list(loc = vertices, graph = list(tv = geometry::delaunayn(vertices)))
  data <- data.frame(x = runif(n, 0.02, 0.98), y = runif(n, 0.02, 0.98),
                     z = seq(-1, 1, length.out = n), exposure = runif(n, 0.8, 1.4))
  data$offset0 <- log(data$exposure)
  basis <- spde_basis(mesh, as.matrix(data[c("x", "y")]), kappa = 0.7,
                      project_intercept = TRUE)
  list(data = data, basis = basis)
}

test_that("spdeST.scan and spdeST.estimate run the mgcv implementation in two calls", {
  f <- st_fixture(n = 200L)
  set.seed(9060L)
  mu <- exp(1 + f$data$offset0 + .5 * f$data$z)
  Y <- t(vapply(c(.8, 2, 6), function(s) stats::rnbinom(200L, mu = mu, size = s),
                numeric(200L)))
  Y <- rbind(Y, f$Y[1L, ])
  rownames(Y) <- paste0("g", 1:4)
  withr::local_options(mgcvST.progress = FALSE)
  sp <- BiocParallel::SerialParam()
  model <- f$model
  model$estimation <- NULL  # defaults when the model predates the settings
  scan <- spdeST.scan(Y, model, BPPARAM = sp)
  expect_s3_class(scan, "spdeST_scan")
  expect_identical(scan$stage1$null_engine, "spde.scan")
  expect_false(any(scan$diagnostics$spatial_fitted))
  direct <- mgcvST.estimate(Y, model, BPPARAM = sp, spatial = "none")
  expect_identical(scan$diagnostics$marginal_p_value, direct$diagnostics$marginal_p_value)

  fit <- spdeST.estimate(scan, Y, features = c("g2", "g4"), BPPARAM = sp)
  expect_s3_class(fit, "spdeST_fit")
  expect_identical(fit$diagnostics$spatial_fitted, c(FALSE, TRUE, FALSE, TRUE))
  both <- mgcvST.estimate_spatial(direct, Y, features = c("g2", "g4"), BPPARAM = sp)
  expect_identical(fit$lambda, both$lambda)
  expect_identical(fit$working_error, both$working_error)
  more <- spdeST.estimate(fit, Y, features = "g1", BPPARAM = sp)
  expect_identical(more$diagnostics$spatial_fitted, c(TRUE, TRUE, FALSE, TRUE))

  seen <- NULL
  testthat::local_mocked_bindings(
    mgcvST.test = function(fitmgcvST, ...) { seen <<- list(...); "mgcv test" },
    .package = "mgcvST")
  expect_identical(spdeST.test(fit, pairs = matrix(c(2L, 4L), 1L), moments = "exact"),
                   "mgcv test")
  expect_identical(seen$moments, "exact")
  expect_error(spdeST.test(scan), "spdeST.estimate")
  expect_error(spdeST.estimate(direct, Y), "spdeST.scan")
  expect_error(spdeST.scan(Y, list()), "mgcvST.set\\(\\) or inlaST.set")
})

test_that("mgcvST.set fixes the estimation settings with the model", {
  f <- st_fixture(n = 60L)
  model <- mgcvST.set(response ~ offset(offset0) + z +
                        s(x, y, bs = "spde", xt = f$basis), f$data,
                      bam_args = list(gamma = 1.2), retain_smooth = TRUE)
  expect_identical(model$estimation$bam_args, list(gamma = 1.2))
  expect_true(model$estimation$retain_smooth)
  expect_error(mgcvST.set(response ~ z, f$data, bam_args = list(1)), "named list")
  expect_error(mgcvST.set(response ~ z, f$data, diagnostics = NA), "diagnostics must")
  calls <- NULL
  testthat::local_mocked_bindings(
    mgcvST.estimate = function(...) { calls <<- list(...); structure(list(), class = "mgcvST_model_fit") },
    .package = "mgcvST")
  spdeST.scan(f$Y, model, chunk_size = 2L)
  expect_identical(calls$gamma, 1.2)
  expect_true(calls$retain_smooth)
  expect_identical(calls$spatial, "none")
  expect_identical(calls$chunk_size, 2L)
  expect_null(calls$BPPARAM)
})

test_that("spdeST.scan and spdeST.estimate run the INLA implementation in two calls", {
  skip_if_not_installed("INLA")
  f <- .spdest_inla_fixture()
  d <- f$data
  n <- nrow(d)
  set.seed(9061L)
  mu <- exp(1 + d$offset0 + .3 * d$z)
  Y <- t(vapply(c(1, 4, 10), function(s) stats::rnbinom(n, mu = mu, size = s),
                numeric(n)))
  rownames(Y) <- paste0("g", 1:3)
  model <- inlaST.set(response ~ z + offset(offset0), d, f$basis,
                      family = mgcv::nb(), retain_smooth = TRUE)
  expect_true(model$estimation$retain_smooth)
  sp <- BiocParallel::SerialParam()
  scan <- spdeST.scan(Y, model, BPPARAM = sp)
  expect_s3_class(scan, c("spdeST_scan", "inlaST_fit"))
  expect_identical(scan$stage1$null_engine, "spde.scan")
  direct <- inlaST.estimate(Y, model, BPPARAM = sp, spatial = "none",
                            retain_smooth = TRUE)
  expect_equal(scan$diagnostics$marginal_p_value,
               direct$diagnostics$marginal_p_value, tolerance = 1e-12)
  fit <- spdeST.estimate(scan, Y, features = "g2", BPPARAM = sp)
  expect_s3_class(fit, "spdeST_fit")
  expect_identical(fit$diagnostics$spatial_fitted, c(FALSE, TRUE, FALSE))
  seen <- NULL
  testthat::local_mocked_bindings(
    inlaST.test = function(fitinlaST, ...) { seen <<- list(...); "inla test" },
    .package = "mgcvST")
  expect_identical(spdeST.test(fit, moments = "exact"), "inla test")
  expect_identical(seen$moments, "exact")
})
