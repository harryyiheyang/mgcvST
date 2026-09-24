.inlast_set_workflow_fixture <- function(n = 42L, seed = 1601L) {
  skip_if_not_installed("INLA")
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

.inlast_expect_same_set_geometry <- function(actual, expected,
                                              tolerance = 1e-10) {
  expect_identical(actual$setting, expected$setting)
  expect_identical(actual$components, expected$components)
  expect_identical(actual$inla_spec$family, expected$inla_spec$family)
  expect_equal(actual$offset, expected$offset, tolerance = tolerance)
  expect_identical(dim(actual$L), dim(expected$L))
  expect_equal(as.numeric(actual$L), as.numeric(expected$L), tolerance = tolerance)
  expect_identical(dim(actual$geometry$X), dim(expected$geometry$X))
  expect_equal(as.numeric(actual$geometry$X), as.numeric(expected$geometry$X),
               tolerance = tolerance)
  expect_identical(names(actual$inla_spec$random),
                   names(expected$inla_spec$random))
  for (j in seq_along(actual$inla_spec$random)) {
    a <- actual$inla_spec$random[[j]]
    b <- expected$inla_spec$random[[j]]
    expect_identical(a$name, b$name)
    expect_identical(a$target, b$target)
    expect_equal(as.matrix(a$A), as.matrix(b$A), tolerance = tolerance)
    expect_equal(as.matrix(a$Q), as.matrix(b$Q), tolerance = tolerance)
    expect_equal(a$constraint, b$constraint, tolerance = tolerance)
    expect_equal(a$projection, b$projection, tolerance = tolerance)
  }
}

test_that("inlaST.set accepts basis, complete-formula and frozen-G workflows", {
  skip_on_cran()
  f <- .inlast_set_workflow_fixture()
  d <- f$data
  basis <- f$basis
  stored <- list(fixed_precision = 1.25, gaussian_precision = 4)

  by_basis <- inlaST.set(
    response ~ z + offset(offset0), d, basis,
    family = gaussian(), coordinates = c("u", "v"), control = stored
  )
  complete <- response ~ z + offset(offset0) +
    s(u, v, bs = "spde", xt = basis)
  # The third positional family is retained for mgcvST.set-style calls.
  by_formula <- inlaST.set(complete, d, gaussian(), control = stored)
  prepared <- mgcvST.set(complete, d, gaussian())
  by_G <- inlaST.set(G = prepared$G, control = stored)

  for (model in list(by_basis, by_formula, by_G)) {
    expect_s3_class(model, "inlaST_model")
    expect_true(model$mean_constraint_active)
    expect_identical(model$mean_constraint, "observation")
    expect_true(model$inla_spec$mean_constraint_active)
    expect_identical(model$setting, "global")
    expect_identical(model$components, "global")
    expect_identical(model$inla_control$fixed_precision, 1.25)
    expect_identical(model$inla_control$gaussian_precision, 4)
  }
  .inlast_expect_same_set_geometry(by_formula, by_basis)
  .inlast_expect_same_set_geometry(by_G, by_basis)

  # The public setup contract accepts one global spatial process only.
  expect_error(
    inlaST.set(complete, d, gaussian(), setting = "global_local",
               control = stored),
    "setting must be \"global\""
  )
  expect_error(
    inlaST.set(G = prepared$G, setting = "global_local", control = stored),
    "setting must be \"global\""
  )
  expect_error(
    model.set(response ~ z + offset(offset0), d, basis,
              family = gaussian(), setting = "global_local",
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
    inlaST.set(complete_two, d, gaussian(), control = stored),
    "exactly one spatial SPDE term"
  )
  local_basis <- second_basis
  local_basis$component <- local_basis$score.component <- "local"
  expect_error(
    inlaST.set(response ~ z + offset(offset0) +
                 s(u, v, bs = "spde", xt = local_basis),
               d, gaussian(), control = stored),
    "xt\\$component must be 'global'"
  )

  uncentred <- spde_basis(
    f$mesh, as.matrix(d[c("u", "v")]), kappa = .7,
    project_intercept = FALSE
  )
  bad_formula <- response ~ z + offset(offset0) +
    s(u, v, bs = "spde", xt = uncentred)
  bad_G <- mgcvST.set(bad_formula, d, gaussian())$G
  expect_error(inlaST.set(G = bad_G), "mean|cent|formula|projection")
})

test_that("native mesh setup preserves the sparse 2D contract and exposes 3D geometry", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  skip_if_not_installed("geometry")
  f <- .inlast_set_workflow_fixture(n = 28L, seed = 1605L)
  control <- list(fixed_precision = 1.25, gaussian_precision = 4)
  legacy <- inlaST.set(
    response ~ z + offset(offset0), f$data, f$basis, family = gaussian(),
    coordinates = c("u", "v"), control = control
  )
  native <- inlaST.set(
    response ~ z + offset(offset0), f$data, family = gaussian(), mesh = f$mesh,
    kappa = .7, coordinates = c("u", "v"), control = control
  )
  expect_s3_class(native, "inlaST_native_model")
  expect_true(native$native)
  expect_equal(as.matrix(native$inla_spec$random[[1L]]$A),
               as.matrix(legacy$inla_spec$random[[1L]]$A), tolerance = 1e-12)
  expect_equal(as.matrix(native$inla_spec$random[[1L]]$Q),
               as.matrix(legacy$inla_spec$random[[1L]]$Q), tolerance = 1e-12)
  expect_equal(native$inla_spec$random[[1L]]$constraint,
               legacy$inla_spec$random[[1L]]$constraint, tolerance = 1e-12)
  native_score <- mgcvST:::.inlast_sparse_score_geometry(native)
  expect_identical(native_score$normalization,
                   ncol(native$inla_spec$random[[1L]]$Q) - 1L)

  skip_if_not_installed("fmesher")
  xyz <- as.matrix(expand.grid(
    x = seq(0, 1, length.out = 3L), y = seq(0, 1, length.out = 3L),
    z = seq(0, 1, length.out = 3L)
  ))
  mesh3 <- fmesher::fm_mesh_3d(loc = xyz, tv = geometry::delaunayn(xyz))
  d3 <- data.frame(x = runif(18L, .05, .95), y = runif(18L, .05, .95),
                   z = runif(18L, .05, .95), a = rnorm(18L))
  native3 <- inlaST.set(
    response ~ a, d3, family = gaussian(), mesh = mesh3, kappa = .7,
    coordinates = c("x", "y", "z"), control = control
  )
  expect_s3_class(native3, "inlaST_native_model")
  expect_identical(native3$spde$dim, 3L)
  expect_identical(ncol(native3$inla_spec$random[[1L]]$A), nrow(xyz))
  native3_score <- mgcvST:::.inlast_sparse_score_geometry(native3)
  expect_identical(native3_score$normalization, nrow(xyz) - 1L)
})

test_that("a fixed NB family does not partially match nb_size_prior", {
  skip_on_cran()
  f <- .inlast_set_workflow_fixture(seed = 1604L)
  prior <- list(prior = "normal", param = c(0, 1 / 9), initial = .2)
  model <- inlaST.set(
    response ~ z + offset(offset0), f$data, f$basis,
    family = mgcv::nb(theta = 4), coordinates = c("u", "v"),
    control = list(nb_size_prior = prior)
  )

  expect_identical(model$inla_control[["nb_size", exact = TRUE]], 4)
  expect_identical(model$inla_control[["nb_size_prior", exact = TRUE]], prior)
})

