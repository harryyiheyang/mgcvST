test_that("the removed threshold control is refused by both entries", {
  expect_error(mgcvST:::.mgcvst_refuse_screen_control(list(poisson_screen_phi = 1.1)),
               "removed")
  expect_error(mgcvST:::.mgcvst_refuse_screen_control(list(poisson_screen_phi = NULL)),
               "removed")
  expect_silent(mgcvST:::.mgcvst_refuse_screen_control(list(nthreads = 1L)))
  skip_if_not_installed("INLA")
  expect_error(mgcvST:::.inlast_control(list(poisson_screen_phi = 0)), "removed")
  expect_null(mgcvST:::.inlast_control()$poisson_screen_phi)
})

test_that("a feature is routed exactly when its NB2 fit is the Poisson fit", {
  Y <- rbind(below = c(rep(c(0, 2), 49), 1, 1),
    between = rep(c(0, 2), 50), above = rep(c(0, 4), 50))
  X <- matrix(1, ncol(Y), 1L)
  z <- mgcvST:::.mgcvst_prescreen_route(Y, X, numeric(ncol(Y)), TRUE)
  expect_equal(unname(z$phi), c(98, 100, 200) / 99, tolerance = 1e-8)
  # score = (sum (y - mu)^2 - sum y) / 2: -1, 0 and 100; 0 is the end point.
  expect_equal(unname(z$score), c(-1, 0, 100), tolerance = 1e-8)
  expect_identical(unname(z$poisson), c(TRUE, TRUE, FALSE))
  sparse <- mgcvST:::.mgcvst_prescreen_route(Matrix::Matrix(Y, sparse = TRUE), X,
    NULL, TRUE, threads = 2L)
  expect_equal(sparse$score, z$score)
  no_x <- mgcvST:::.mgcvst_prescreen_route(Y, matrix(numeric(), ncol(Y), 0L),
    numeric(ncol(Y)), TRUE)
  expect_equal(unname(no_x$score), c(-1, 0, 150), tolerance = 1e-8)
  inactive <- mgcvST:::.mgcvst_prescreen_route(Y, X, numeric(ncol(Y)), FALSE)
  expect_identical(inactive$poisson, rep(FALSE, 3))
  expect_true(all(is.na(inactive$phi)))
})

test_that("public mgcv and INLA entries apply the same routing", {
  skip_if_not_installed("INLA")
  f <- st_fixture(n = 100L)
  f$data$offset0 <- 0
  model <- mgcvST.set(response ~ offset(offset0) + s(x, y, bs = "spde", xt = f$basis),
    f$data, family = mgcv::nb())
  mesh <- list(loc = f$basis$mesh_vertices, graph = list(tv = f$basis$mesh_triangles))
  model_inla <- inlaST.set(response ~ offset(offset0), f$data,
    family = mgcv::nb(), mesh = mesh, kappa = .7,
    coordinates = c("x", "y"))
  Y <- rbind(below = c(rep(c(0, 2), 49), 1, 1), above = rep(c(0, 4), 50))
  route <- mgcvST:::.mgcvst_prescreen_route
  captured <- new.env(parent = emptyenv())
  testthat::local_mocked_bindings(.mgcvst_prescreen_route = function(Y, X, offset, active, threads = 1L) {
    captured$result <- route(Y, X, offset, active, threads)
    stop("prescreen captured before model fitting")
  }, .package = "mgcvST")
  for (setup in list(model, model$G)) {
    expect_error(mgcvST.estimate(Y, setup, control = mgcv::gam.control(),
      BPPARAM = BiocParallel::SerialParam(), spatial = "all"),
      "prescreen captured before model fitting")
    expect_identical(unname(captured$result$poisson), c(TRUE, FALSE))
  }
  expect_error(inlaST.estimate(Y, model_inla,
    BPPARAM = BiocParallel::SerialParam(), spatial = "all"),
    "prescreen captured before model fitting")
  expect_identical(unname(captured$result$poisson), c(TRUE, FALSE))
  control <- mgcv::gam.control()
  control$poisson_screen_phi <- 1.1
  expect_error(mgcvST.estimate(Y, model, control = control,
    BPPARAM = BiocParallel::SerialParam()), "removed")
})
