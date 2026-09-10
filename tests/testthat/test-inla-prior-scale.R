prior_scale_fixture <- function(kappa = .7) {
  vertices <- as.matrix(expand.grid(x = seq(0, 1, length.out = 5),
                                  y = seq(0, 1, length.out = 5)))
  mesh <- list(loc = vertices, graph = list(tv = geometry::delaunayn(vertices)))
  set.seed(1382)
  d <- data.frame(x = runif(50, .01, .99), y = runif(50, .01, .99))
  basis <- spde_basis(mesh, as.matrix(d), kappa = kappa, project_intercept = TRUE)
  list(d = d, basis = basis)
}

test_that("observation precision prior has unit mean constrained variance", {
  skip_if_not_installed("INLA")
  for (kappa in c(.1, .7, 6)) {
    z <- prior_scale_fixture(kappa)
    model <- inlaST.set(response ~ 1, z$d, z$basis, family = gaussian(),
                        precision_scale = "observation")
    b <- model$inla_spec$random[[1]]
    B <- as.matrix(b$A %*% b$projection)
    Q <- crossprod(b$projection, b$Q %*% b$projection)
    G <- B %*% solve(as.matrix(Q), t(B)) / b$precision_scale
    expect_equal(mean(diag(G)), 1, tolerance = 1e-9)
    expect_lt(max(abs(colMeans(B))), 1e-12)
    expect_identical(model$precision_scale, "observation")
  }
  z <- prior_scale_fixture()
  raw <- inlaST.set(response ~ 1, z$d, z$basis, family = gaussian())
  expect_identical(raw$inla_spec$random[[1]]$precision_scale, 1)
  expect_error(inlaST.set(response ~ 1, z$d, z$basis,
                          precision_scale = "unknown"), "arg")
})

test_that("prior scaling preserves a fixed original-precision model", {
  skip_if_not_installed("INLA")
  z <- prior_scale_fixture()
  raw <- inlaST.set(response ~ 1, z$d, z$basis, family = gaussian())
  scaled <- inlaST.set(response ~ 1, z$d, z$basis, family = gaussian(),
                       precision_scale = "observation")
  y <- .4 + sin(z$d$x * 4) + rnorm(nrow(z$d), sd = .25)
  control <- list(fixed_precision = .3, gaussian_precision = 16)
  a <- mgcvST:::.inlast_fit_feature(raw$inla_spec, y, raw$offset, control)
  b <- mgcvST:::.inlast_fit_feature(scaled$inla_spec, y, scaled$offset, control)
  expect_equal(unname(a$tau), .3, tolerance = 1e-12)
  expect_equal(unname(b$tau), .3, tolerance = 1e-12)
  expect_equal(a$eta, b$eta, tolerance = 1e-5)
  expect_equal(a$nuisance_covariance, b$nuisance_covariance, tolerance = 1e-5)
  expect_equal(unname(b$tau_internal * b$precision_scale), .3, tolerance = 1e-12)
  expect_lt(max(abs(b$observation_spatial_mean)), 1e-10)
})
