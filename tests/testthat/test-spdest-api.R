.spdest_inla_fixture <- function(n = 120L, seed = 2402L) {
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

test_that("spdeST.scan and spdeST.estimate take the model of inlaST.set in two calls", {
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
  direct <- inlaST.estimate(Y, model, BPPARAM = sp,
                            retain_smooth = TRUE)
  expect_equal(scan$diagnostics$marginal_p_value,
               direct$diagnostics$marginal_p_value, tolerance = 1e-12)
  fit <- spdeST.estimate(scan, Y, features = "g2", BPPARAM = sp)
  expect_s3_class(fit, "spdeST_fit")
  expect_identical(fit$diagnostics$spatial_fitted, c(FALSE, TRUE, FALSE))
  more <- spdeST.estimate(fit, Y, features = "g1", BPPARAM = sp)
  expect_identical(more$diagnostics$spatial_fitted, c(TRUE, TRUE, FALSE))
  seen <- NULL
  testthat::local_mocked_bindings(
    inlaST.test = function(fitinlaST, ...) { seen <<- list(...); "pair test" },
    .package = "spdeST")
  expect_identical(spdeST.test(fit, moments = "exact"), "pair test")
  expect_identical(seen$moments, "exact")
  expect_error(spdeST.test(scan), "spdeST.estimate")
  expect_error(spdeST.estimate(direct, Y), "spdeST.scan")
  expect_error(spdeST.scan(Y, list()), "inlaST.set")
})
