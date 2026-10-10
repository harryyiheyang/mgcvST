test_that("both control paths share the default and preserve explicit thresholds", {
  expect_identical(mgcvST:::.mgcvst_prescreen_threshold(NULL), 1.01)
  expect_identical(mgcvST:::.inlast_control()$poisson_screen_phi, 1.01)
  for (value in c(0, 1.005, 1.1)) {
    expected <- if (value == 0) NULL else value
    expect_identical(mgcvST:::.mgcvst_prescreen_threshold(value), expected)
    control <- mgcvST:::.inlast_control(list(poisson_screen_phi = value))
    expect_identical(control$poisson_screen_phi, value)
    expect_identical(mgcvST:::.mgcvst_prescreen_threshold(control$poisson_screen_phi), expected)
  }
  control <- mgcvST:::.inlast_control(list(poisson_screen_phi = NULL))
  expect_identical(mgcvST:::.mgcvst_prescreen_threshold(control$poisson_screen_phi), 1.01)
})

test_that("actual Poisson GLMs route at the tighter default", {
  Y <- rbind(below = c(rep(c(0, 2), 49), 1, 1),
    between = rep(c(0, 2), 50), above = rep(c(0, 4), 50))
  X <- matrix(1, ncol(Y), 1L)
  threshold <- mgcvST:::.mgcvst_prescreen_threshold(NULL)
  z <- mgcvST:::.mgcvst_prescreen_route(Y, X, numeric(ncol(Y)), threshold, TRUE)
  expect_equal(unname(z$phi), c(98, 100, 200) / 99, tolerance = 1e-8)
  expect_identical(unname(z$poisson), c(TRUE, FALSE, FALSE))
  old <- mgcvST:::.mgcvst_prescreen_route(Y, X, numeric(ncol(Y)), 1.1, TRUE)
  expect_identical(unname(old$poisson), c(TRUE, TRUE, FALSE))
  off <- mgcvST:::.mgcvst_prescreen_route(Y, X, numeric(ncol(Y)),
    mgcvST:::.mgcvst_prescreen_threshold(0), TRUE)
  expect_identical(off, list(phi = rep(NA_real_, 3), poisson = rep(FALSE, 3)))
  inactive <- mgcvST:::.mgcvst_prescreen_route(Y, X, numeric(ncol(Y)), threshold, FALSE)
  expect_identical(inactive, off)
})

test_that("the exact boundary is inclusive without changing phi", {
  Y <- matrix(c(rep(0, 97), 3, 1, 1), nrow = 1L)
  X <- matrix(numeric(), ncol(Y), 0L)
  z <- mgcvST:::.mgcvst_prescreen_route(Y, X, NULL, 1.01, TRUE)
  expect_identical(unname(z$phi), 1.01)
  expect_identical(unname(z$poisson), TRUE)
  below <- mgcvST:::.mgcvst_prescreen_route(Y, X, NULL, 1.01 - 1e-12, TRUE)
  expect_identical(below$phi, z$phi)
  expect_identical(unname(below$poisson), FALSE)
})

test_that("disabled screening is unchanged by the default value", {
  control <- mgcvST:::.inlast_control(list(poisson_screen_phi = 0))
  threshold <- mgcvST:::.mgcvst_prescreen_threshold(control$poisson_screen_phi)
  testthat::local_mocked_bindings(.mgcvst_prescreen_default = 1.1,
    .mgcvst_prescreen_dispersion = function(...) stop("disabled screen must not fit a GLM"),
    .package = "mgcvST")
  expect_identical(mgcvST:::.inlast_control(list(poisson_screen_phi = 0)), control)
  expect_identical(mgcvST:::.mgcvst_prescreen_threshold(0), threshold)
  expect_identical(mgcvST:::.mgcvst_prescreen_route(matrix(1, 2L, 10L), NULL, NULL,
    threshold, TRUE), list(phi = rep(NA_real_, 2), poisson = rep(FALSE, 2)))
})

test_that("public mgcv and INLA entries apply the same routing controls", {
  skip_if_not_installed("INLA")
  f <- st_fixture(n = 100L)
  f$data$offset0 <- 0
  model <- mgcvST.set(response ~ offset(offset0) + s(x, y, bs = "spde", xt = f$basis),
    f$data, family = mgcv::nb())
  mesh <- list(loc = f$basis$mesh_vertices, graph = list(tv = f$basis$mesh_triangles))
  model_inla <- inlaST.set(response ~ offset(offset0), f$data,
    family = mgcv::nb(), mesh = mesh, kappa = .7,
    coordinates = c("x", "y"))
  Y <- rbind(below = c(rep(c(0, 2), 49), 1, 1), between = rep(c(0, 2), 50))
  route <- mgcvST:::.mgcvst_prescreen_route
  captured <- new.env(parent = emptyenv())
  testthat::local_mocked_bindings(.mgcvst_prescreen_route = function(Y, X, offset, threshold, active) {
    captured$result <- route(Y, X, offset, threshold, active)
    captured$threshold <- threshold
    stop("prescreen captured before model fitting")
  }, .package = "mgcvST")
  controls <- list(list(), list(poisson_screen_phi = NULL),
    list(poisson_screen_phi = 1.1), list(poisson_screen_phi = 0))
  for (k in seq_along(controls)) {
    expected <- if (k == 3L) c(TRUE, TRUE) else if (k == 4L) c(FALSE, FALSE) else c(TRUE, FALSE)
    threshold <- if (k == 3L) 1.1 else if (k == 4L) NULL else 1.01
    control <- utils::modifyList(mgcv::gam.control(), controls[[k]], keep.null = TRUE)
    for (setup in list(model, model$G)) {
      expect_error(mgcvST.estimate(Y, setup, control = control,
        BPPARAM = BiocParallel::SerialParam(), spatial = "all"), "prescreen captured before model fitting")
      expect_identical(captured$threshold, threshold)
      expect_identical(unname(captured$result$poisson), expected)
    }
    expect_error(inlaST.estimate(Y, model_inla, control = controls[[k]],
      BPPARAM = BiocParallel::SerialParam(), spatial = "all"), "prescreen captured before model fitting")
    expect_identical(captured$threshold, threshold)
    expect_identical(unname(captured$result$poisson), expected)
  }
})
