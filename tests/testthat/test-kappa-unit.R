.kappa_unit_mesh2d <- function(scale = 1) {
  vertices <- as.matrix(expand.grid(
    x = seq(0, 4, length.out = 6L), y = seq(0, 3, length.out = 5L)
  ))
  tv <- geometry::delaunayn(vertices)
  list(loc = vertices * scale, graph = list(tv = tv))
}

.kappa_unit_data2d <- function(n = 48L, seed = 3101L, scale = 1) {
  set.seed(seed)
  d <- data.frame(
    x = runif(n, 0.3, 3.7), y = runif(n, 0.2, 2.6),
    z = seq(-1, 1, length.out = n)
  )
  d$x <- d$x * scale
  d$y <- d$y * scale
  d
}

.kappa_unit_span <- function(d, columns) {
  apply(as.matrix(d[columns]), 2L, function(z) max(z) - min(z))
}

# Covariance at the observations of the field conditioned on g'u = 0.
.kappa_unit_kernel <- function(A, Q, g) {
  Z <- qr.Q(qr(matrix(g, ncol = 1L)), complete = TRUE)[, -1L, drop = FALSE]
  B <- as.matrix(A %*% Z)
  QZ <- crossprod(Z, as.matrix(Q %*% Z))
  B %*% solve((QZ + t(QZ)) / 2, t(B))
}

.kappa_unit_expect_proportional <- function(K1, K2, tolerance = 1e-8) {
  ratio <- sum(K1 * K2) / sum(K1 * K1)
  expect_gt(ratio, 0)
  expect_lt(max(abs(K2 - ratio * K1)) / max(abs(K2)), tolerance)
  ratio
}

test_that("kappa = NULL is an error because mgcvST never estimates kappa", {
  skip_if_not_installed("geometry")
  mesh <- .kappa_unit_mesh2d()
  d <- .kappa_unit_data2d()
  loc <- as.matrix(d[c("x", "y")])
  expect_identical(formals(spde_basis)$kappa, 0.05)
  expect_identical(formals(inlaST.set)$kappa, 0.05)
  expect_error(spde_basis(mesh, loc, kappa = NULL), "never estimates kappa")
  expect_error(inlaST.set(kappa = NULL), "never estimates kappa")
  expect_error(
    inlaST.set(response ~ z, d, family = gaussian(), mesh = mesh,
               kappa = NULL, coordinates = c("x", "y")),
    "never estimates kappa"
  )
  expect_error(spde_basis(mesh, loc, kappa = -1), "unit-scale")
  expect_error(spde_basis(mesh, loc, kappa = c(0.05, 0.1)), "unit-scale")
  basis <- spde_basis(mesh, loc)
  expect_null(basis$penalty)
  expect_error(
    inlaST.set(response ~ z, d, basis, family = gaussian(), kappa = 0.1),
    "applies only to the native mesh setup"
  )
  stale <- basis
  stale$kappa_internal <- NULL
  expect_error(model.set(response ~ z, d, stale, family = gaussian()),
               "rebuild it with spde_basis")
})

test_that("native fmesher meshes store kappa_unit and kappa_internal = kappa_unit / L", {
  skip_if_not_installed("INLA")
  skip_if_not_installed("fmesher")
  skip_if_not_installed("geometry")
  d <- .kappa_unit_data2d()
  mesh <- fmesher::fm_rcdt_2d_inla(loc = .kappa_unit_mesh2d()$loc)
  model <- inlaST.set(response ~ z, d, family = gaussian(), mesh = mesh,
                      coordinates = c("x", "y"))
  span <- .kappa_unit_span(d, c("x", "y"))
  L <- max(span)
  expect_identical(model$kappa_unit, 0.05)
  expect_equal(model$unit_length, L)
  expect_equal(model$coordinate_span, span)
  expect_equal(model$kappa_internal, 0.05 / L)
  expect_equal(model$spde$range_unit, sqrt(8) / 0.05)
  fem <- mgcvST:::.spde_basis_fem(mgcvST:::.spde_basis_mesh(mesh))
  k <- 0.05 / L
  expect_equal(as.matrix(model$inla_spec$random[[1L]]$Q),
               as.matrix(k^4 * fem$M0 + 2 * k^2 * fem$M1 + fem$M2),
               tolerance = 1e-12)
  expect_output(print(model), "kappa \\(unit scale, fixed\\): 0.05")
  expect_output(print(model), "unit length L")

  set.seed(3102L)
  xyz <- as.matrix(expand.grid(x = seq(0, 2, length.out = 3L),
                               y = seq(0, 1, length.out = 3L),
                               z = seq(0, 3, length.out = 3L)))
  mesh3 <- fmesher::fm_mesh_3d(loc = xyz, tv = geometry::delaunayn(xyz))
  d3 <- data.frame(x = runif(30L, 0.1, 1.9), y = runif(30L, 0.1, 0.9),
                   z = runif(30L, 0.2, 2.5), a = rnorm(30L))
  model3 <- inlaST.set(response ~ a, d3, family = gaussian(), mesh = mesh3,
                       coordinates = c("x", "y", "z"))
  span3 <- .kappa_unit_span(d3, c("x", "y", "z"))
  L3 <- max(span3)
  expect_identical(model3$kappa_unit, 0.05)
  expect_equal(model3$unit_length, L3)
  expect_equal(model3$coordinate_span, span3)
  expect_equal(model3$kappa_internal, 0.05 / L3)
  expect_equal(model3$spde$range_unit, 2 / 0.05)
  fem3 <- fmesher::fm_fem(mesh3, order = 2L)
  k3 <- 0.05 / L3
  expect_equal(as.matrix(model3$inla_spec$random[[1L]]$Q),
               as.matrix(k3^4 * fem3$c0 + 2 * k3^2 * fem3$g1 + fem3$g2),
               tolerance = 1e-12)
})

test_that("the mgcv path uses kappa_internal = kappa_unit * transform_scale / L", {
  skip_if_not_installed("geometry")
  d <- .kappa_unit_data2d()
  loc <- as.matrix(d[c("x", "y")])
  raw <- .kappa_unit_mesh2d()
  center <- c(2, 1.5)
  scale <- 4
  normalized <- structure(list(
    mesh = list(loc = sweep(raw$loc, 2L, center, "-") / scale,
                graph = raw$graph),
    transform = list(center = center, scale = scale)
  ), class = "spde_mesh")
  span <- .kappa_unit_span(d, c("x", "y"))
  L <- max(span)

  a <- spde_basis(raw, loc, kappa = 0.3)
  b <- spde_basis(normalized, loc, kappa = 0.3)
  expect_identical(b$kappa_unit, 0.3)
  expect_equal(b$unit_length, L)
  expect_equal(unname(b$coordinate_span), unname(span))
  expect_equal(b$kappa_internal, 0.3 * scale / L)
  expect_equal(a$kappa_internal, 0.3 / L)
  fem <- mgcvST:::.spde_basis_fem(mgcvST:::.spde_basis_mesh(normalized))
  k <- 0.3 * scale / L
  Qraw <- as.matrix(k^4 * fem$M0 + 2 * k^2 * fem$M1 + fem$M2)
  expect_equal(b$Q, crossprod(b$projection, Qraw %*% b$projection),
               tolerance = 1e-10)
  # Normalizing coordinates by s multiplies the 2D precision by s^2 only.
  expect_equal(b$B, a$B, tolerance = 1e-12)
  expect_equal(b$Q, scale^2 * a$Q, tolerance = 1e-10)
  expect_output(print(b), "kappa \\(unit scale, fixed\\): 0.3")

  model <- model.set(response ~ z, d, b, family = gaussian())
  expect_identical(model$kappa_unit, 0.3)
  expect_equal(model$kappa_internal, 0.3 * scale / L)
  expect_equal(model$unit_length, L)
  expect_output(print(model), "unit length L")
  prepared <- mgcvST.set(response ~ z + s(x, y, bs = "spde", xt = b), d,
                         gaussian())
  expect_equal(prepared$kappa_internal, 0.3 * scale / L)
  expect_equal(prepared$G$smooth[[1L]]$kappa_internal, 0.3 * scale / L)
  expect_null(prepared$G$smooth[[1L]]$kappa.estimated)
})

test_that("rescaling mm to um leaves the kernel and marginal scores unchanged", {
  skip_if_not_installed("geometry")
  mm <- .kappa_unit_data2d()
  um <- .kappa_unit_data2d(scale = 1000)
  mesh_mm <- .kappa_unit_mesh2d()
  mesh_um <- .kappa_unit_mesh2d(scale = 1000)

  # mgcv path: the projected SPDE covariance B Q^-1 B'.
  b_mm <- spde_basis(mesh_mm, as.matrix(mm[c("x", "y")]))
  b_um <- spde_basis(mesh_um, as.matrix(um[c("x", "y")]))
  expect_equal(b_um$unit_length / b_mm$unit_length, 1000)
  expect_equal(b_um$kappa_internal * 1000, b_mm$kappa_internal)
  ratio <- .kappa_unit_expect_proportional(
    b_mm$B %*% solve(b_mm$Q, t(b_mm$B)), b_um$B %*% solve(b_um$Q, t(b_um$B))
  )
  expect_equal(ratio, 1e6, tolerance = 1e-8)

  # Native INLA path: constrained kernel and sparse marginal score.
  skip_if_not_installed("INLA")
  m_mm <- inlaST.set(response ~ z, mm, family = gaussian(), mesh = mesh_mm,
                     coordinates = c("x", "y"))
  m_um <- inlaST.set(response ~ z, um, family = gaussian(), mesh = mesh_um,
                     coordinates = c("x", "y"))
  expect_identical(m_um$kappa_unit, m_mm$kappa_unit)
  r_mm <- m_mm$inla_spec$random[[1L]]
  r_um <- m_um$inla_spec$random[[1L]]
  expect_equal(as.matrix(r_um$A), as.matrix(r_mm$A), tolerance = 1e-12)
  ratio <- .kappa_unit_expect_proportional(
    .kappa_unit_kernel(r_mm$A, r_mm$Q, r_mm$constraint),
    .kappa_unit_kernel(r_um$A, r_um$Q, r_um$constraint)
  )
  expect_equal(ratio, 1e6, tolerance = 1e-8)

  set.seed(3103L)
  n <- nrow(mm)
  E <- cbind(rnorm(n), sin(mm$x) + rnorm(n, sd = 0.2))
  V <- cbind(rep(1, n), runif(n, 0.5, 2))
  X <- m_mm$geometry$nuisance_design
  s_mm <- mgcvST:::.inlast_sparse_null_batch(
    mgcvST:::.inlast_sparse_score_geometry(m_mm), X, E, V, NULL)
  s_um <- mgcvST:::.inlast_sparse_null_batch(
    mgcvST:::.inlast_sparse_score_geometry(m_um), X, E, V, NULL)
  for (j in seq_len(ncol(E))) {
    expect_equal(s_um[[j]]$statistic, s_mm[[j]]$statistic, tolerance = 1e-8)
    expect_equal(s_um[[j]]$moments, s_mm[[j]]$moments, tolerance = 1e-8)
    expect_equal(
      mgcvST:::.mgcvst_marginal_liu(s_um[[j]]$statistic, s_um[[j]]$moments),
      mgcvST:::.mgcvst_marginal_liu(s_mm[[j]]$statistic, s_mm[[j]]$moments),
      tolerance = 1e-8
    )
  }

  skip_on_cran()
  eta <- 0.3 + 0.4 * mm$z + 0.5 * sin(mm$x) * cos(mm$y)
  Y <- rbind(a = eta + rnorm(n, sd = 0.3), b = eta + rnorm(n, sd = 0.3))
  # Fixed precisions keep tau * Q identical: Q_um = 1000^-2 Q_mm in 2D.
  fit_mm <- inlaST.estimate(
    Y, m_mm, BPPARAM = BiocParallel::SerialParam(),
    control = list(fixed_precision = 2, gaussian_precision = 1 / 0.09)
  )
  fit_um <- inlaST.estimate(
    Y, m_um, BPPARAM = BiocParallel::SerialParam(),
    control = list(fixed_precision = 2e6, gaussian_precision = 1 / 0.09)
  )
  expect_identical(fit_um$kappa_unit, 0.05)
  expect_equal(fit_um$unit_length / fit_mm$unit_length, 1000)
  expect_true(all(is.finite(fit_mm$diagnostics$marginal_p_value)))
  expect_equal(fit_um$diagnostics$marginal_p_value,
               fit_mm$diagnostics$marginal_p_value, tolerance = 1e-6)
  expect_equal(fit_um$score_a, fit_mm$score_a, tolerance = 1e-6)
  pairs <- matrix(c("a", "b"), 1L)
  t_mm <- inlaST.test(fit_mm, pairs = pairs, approximate_test = FALSE)$result
  t_um <- inlaST.test(fit_um, pairs = pairs, approximate_test = FALSE)$result
  # The exact pair path packs normalized states in fp16.
  expect_equal(t_um$score, t_mm$score, tolerance = 1e-3)
  expect_equal(t_um$mlog10p, t_mm$mlog10p, tolerance = 1e-3)
})
