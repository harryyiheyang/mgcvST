.spdest_set_workflow_fixture <- function(n = 42L, seed = 1601L) {
  skip_if_not_installed("geometry")
  set.seed(seed)
  vertices <- as.matrix(expand.grid(
    u = seq(0, 1, length.out = 4L),
    v = seq(0, 1, length.out = 4L)
  ))
  mesh <- list(loc = vertices, graph = list(tv = geometry::delaunayn(vertices)))
  data <- data.frame(
    u = runif(n, .03, .97), v = runif(n, .03, .97),
    z = seq(-1, 1, length.out = n),
    offset0 = seq(-.12, .12, length.out = n)
  )
  basis <- spde_basis(
    mesh, as.matrix(data[c("u", "v")]), kappa = .7,
    project_intercept = TRUE
  )
  basis$component <- basis$score.component <- "global"
  eta <- .4 + .25 * data$z + data$offset0 +
    .2 * sin(2 * pi * data$u) - .15 * cos(2 * pi * data$v)
  y <- eta + rnorm(n, sd = .35)
  Y <- rbind(feature_a = y, feature_b = y + rnorm(n, sd = .08))
  list(data = data, basis = basis, mesh = mesh, Y = Y)
}

.spdest_expect_same_set_geometry <- function(actual, expected,
                                              tolerance = 1e-10) {
  expect_identical(actual$setting, expected$setting)
  expect_identical(actual$components, expected$components)
  expect_identical(actual$spde_spec$family, expected$spde_spec$family)
  expect_equal(actual$offset, expected$offset, tolerance = tolerance)
  expect_identical(dim(actual$L), dim(expected$L))
  expect_equal(as.numeric(actual$L), as.numeric(expected$L), tolerance = tolerance)
  expect_identical(dim(actual$geometry$X), dim(expected$geometry$X))
  expect_equal(as.numeric(actual$geometry$X), as.numeric(expected$geometry$X),
               tolerance = tolerance)
  expect_identical(names(actual$spde_spec$random),
                   names(expected$spde_spec$random))
  for (j in seq_along(actual$spde_spec$random)) {
    a <- actual$spde_spec$random[[j]]
    b <- expected$spde_spec$random[[j]]
    expect_identical(a$name, b$name)
    expect_identical(a$target, b$target)
    expect_equal(as.matrix(a$A), as.matrix(b$A), tolerance = tolerance)
    expect_equal(as.matrix(a$Q), as.matrix(b$Q), tolerance = tolerance)
    expect_equal(a$constraint, b$constraint, tolerance = tolerance)
    expect_equal(a$projection, b$projection, tolerance = tolerance)
  }
}

test_that("spdeST.set accepts basis, complete-formula and frozen-G workflows", {
  skip_on_cran()
  f <- .spdest_set_workflow_fixture()
  d <- f$data
  basis <- f$basis
  stored <- list(spde = list(tol = 1e-9))

  by_basis <- spdeST.set(
    response ~ z + offset(offset0), d, basis,
    family = poisson(), coordinates = c("u", "v"), control = stored
  )
  complete <- response ~ z + offset(offset0) +
    s(u, v, bs = "spde", xt = basis)
  # The family may be the third positional argument.
  by_formula <- spdeST.set(complete, d, poisson(), control = stored)
  prepared <- spdeST:::.mgcvst_set_prepare(complete, d, poisson(), .allow_poisson = TRUE)
  by_G <- spdeST.set(G = prepared$G, control = stored)

  for (model in list(by_basis, by_formula, by_G)) {
    expect_s3_class(model, "spdeST_model")
    expect_true(model$mean_constraint_active)
    expect_identical(model$mean_constraint, "observation")
    expect_true(model$spde_spec$mean_constraint_active)
    expect_identical(model$setting, "global")
    expect_identical(model$components, "global")
    expect_identical(model$spde_control$spde$tol, 1e-9)
  }
  .spdest_expect_same_set_geometry(by_formula, by_basis)
  .spdest_expect_same_set_geometry(by_G, by_basis)

  # The public setup contract accepts one global spatial process only.
  expect_error(
    spdeST.set(complete, d, poisson(), setting = "global_local",
               control = stored),
    "setting must be \"global\""
  )
  expect_error(
    spdeST.set(G = prepared$G, setting = "global_local", control = stored),
    "setting must be \"global\""
  )
  expect_error(
    model.set(response ~ z + offset(offset0), d, basis,
              family = poisson(), setting = "global_local",
              coordinates = c("u", "v")),
    "setting must be \"global\""
  )
  # A second SPDE term is a duplicate, and the component tag must be global.
  second_basis <- spde_basis(
    f$mesh, as.matrix(d[c("u", "v")]), kappa = 5, project_intercept = TRUE
  )
  d$u_local <- d$u
  d$v_local <- d$v
  complete_two <- response ~ z + offset(offset0) +
    s(u, v, bs = "spde", xt = basis) +
    s(u_local, v_local, bs = "spde", xt = second_basis)
  expect_error(
    spdeST.set(complete_two, d, poisson(), control = stored),
    "exactly one spatial SPDE term"
  )
  local_basis <- second_basis
  local_basis$component <- local_basis$score.component <- "local"
  expect_error(
    spdeST.set(response ~ z + offset(offset0) +
                 s(u, v, bs = "spde", xt = local_basis),
               d, poisson(), control = stored),
    "xt\\$component must be 'global'"
  )

  uncentred <- spde_basis(
    f$mesh, as.matrix(d[c("u", "v")]), kappa = .7,
    project_intercept = FALSE
  )
  bad_formula <- response ~ z + offset(offset0) +
    s(u, v, bs = "spde", xt = uncentred)
  bad_G <- spdeST:::.mgcvst_set_prepare(bad_formula, d, poisson(), .allow_poisson = TRUE)$G
  expect_error(spdeST.set(G = bad_G), "mean|cent|formula|projection")
})

test_that("native mesh setup preserves the sparse 2D contract and exposes 3D geometry", {
  skip_on_cran()
  skip_if_not_installed("geometry")
  f <- .spdest_set_workflow_fixture(n = 28L, seed = 1605L)
  control <- list(spde = list(tol = 1e-9))
  legacy <- spdeST.set(
    response ~ z + offset(offset0), f$data, f$basis, family = poisson(),
    coordinates = c("u", "v"), control = control
  )
  native <- spdeST.set(
    response ~ z + offset(offset0), f$data, family = poisson(), mesh = f$mesh,
    kappa = .7, coordinates = c("u", "v"), control = control
  )
  expect_s3_class(native, "spdeST_native_model")
  expect_true(native$native)
  expect_equal(as.matrix(native$spde_spec$random[[1L]]$A),
               as.matrix(legacy$spde_spec$random[[1L]]$A), tolerance = 1e-12)
  expect_equal(as.matrix(native$spde_spec$random[[1L]]$Q),
               as.matrix(legacy$spde_spec$random[[1L]]$Q), tolerance = 1e-12)
  expect_equal(native$spde_spec$random[[1L]]$constraint,
               legacy$spde_spec$random[[1L]]$constraint, tolerance = 1e-12)
  native_score <- spdeST:::.spdest_sparse_score_geometry(native)
  expect_identical(native_score$normalization,
                   ncol(native$spde_spec$random[[1L]]$Q) - 1L)

  skip_if_not_installed("fmesher")
  xyz <- as.matrix(expand.grid(
    x = seq(0, 1, length.out = 3L), y = seq(0, 1, length.out = 3L),
    z = seq(0, 1, length.out = 3L)
  ))
  mesh3 <- fmesher::fm_mesh_3d(loc = xyz, tv = geometry::delaunayn(xyz))
  d3 <- data.frame(x = runif(18L, .05, .95), y = runif(18L, .05, .95),
                   z = runif(18L, .05, .95), a = rnorm(18L))
  native3 <- spdeST.set(
    response ~ a, d3, family = poisson(), mesh = mesh3, kappa = .7,
    coordinates = c("x", "y", "z"), control = control
  )
  expect_s3_class(native3, "spdeST_native_model")
  expect_identical(native3$spde$dim, 3L)
  expect_identical(ncol(native3$spde_spec$random[[1L]]$A), nrow(xyz))
  native3_score <- spdeST:::.spdest_sparse_score_geometry(native3)
  expect_identical(native3_score$normalization, nrow(xyz) - 1L)
})

test_that("a family outside the model of spde.estimate() is refused at setup", {
  skip_on_cran()
  f <- .spdest_set_workflow_fixture(seed = 1604L)
  expect_error(spdeST.set(response ~ z + offset(offset0), f$data, f$basis,
                          family = mgcv::nb(theta = 4), coordinates = c("u", "v")),
               "fixed theta")
  expect_error(spdeST.set(response ~ z + offset(offset0), f$data, family = mgcv::nb(theta = 4),
                          mesh = f$mesh, kappa = .7, coordinates = c("u", "v")),
               "fixed theta")
  expect_error(spdeST.set(response ~ z + offset(offset0), f$data, f$basis,
                          family = gaussian(), coordinates = c("u", "v")),
               "Poisson|negative-binomial")
  expect_error(spdeST.set(response ~ z, f$data, f$basis, coordinates = c("u", "v"),
                          control = list(fixed_precision = 2)), "Unknown control field")
})
