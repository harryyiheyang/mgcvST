.native_iid_mesh3d <- function() {
  xyz <- as.matrix(expand.grid(
    x = c(0, 1), y = c(0, 1), z = c(0, 1)
  ))
  fmesher::fm_mesh_3d(loc = xyz, tv = geometry::delaunayn(xyz))
}

.native_iid_data <- function(n = 36L, levels = 3L, seed = 2501L) {
  set.seed(seed)
  xyz <- matrix(runif(3L * n), ncol = 3L)
  xyz <- xyz / pmax(1, rowSums(xyz) / .92)
  slide <- factor(
    rep(paste0("slide", seq_len(levels)), length.out = n),
    levels = c("unused", paste0("slide", seq_len(levels)))
  )
  cell <- factor(rep(letters[1:5], length.out = n))
  batch <- factor(rep(paste0("batch", 1:4), length.out = n))
  data.frame(
    x = xyz[, 1L], y = xyz[, 2L], z = xyz[, 3L],
    offset0 = seq(-.2, .2, length.out = n), age = seq(-1, 1, length.out = n),
    cell = cell, slide = slide, batch = batch
  )
}

test_that("native iid design preserves multiple random-effect terms in order", {
  d <- .native_iid_data(n = 18L, levels = 3L)
  f <- response ~ cell + age + offset(offset0) +
    s(slide, bs = "re") + s(batch, bs = "re")
  design <- mgcvST:::.inlast_native_design(f, d)
  d$response <- 0
  fixed <- stats::model.frame(response ~ cell + age + offset(offset0), d,
                              na.action = stats::na.fail)

  expect_identical(colnames(design$X),
                   colnames(stats::model.matrix(stats::terms(fixed), fixed)))
  expect_equal(design$X,
               stats::model.matrix(stats::terms(fixed), fixed))
  expect_equal(design$offset, stats::model.offset(fixed))
  expect_identical(vapply(design$nuisance, `[[`, character(1L), "type"),
                   c("re", "re"))
  expect_identical(design$nuisance[[1L]]$levels, paste0("slide", 1:3))
  expect_identical(design$nuisance[[2L]]$levels, paste0("batch", 1:4))
  expect_identical(colnames(design$nuisance[[1L]]$Z),
                   paste0("slide:slide", 1:3))
  expect_equal(rowSums(as.matrix(design$nuisance[[1L]]$Z)), rep(1, nrow(d)))
  expect_equal(rowSums(as.matrix(design$nuisance[[2L]]$Z)), rep(1, nrow(d)))
})

test_that("native iid maps fixed and random nuisance columns in combined order", {
  skip_if_not_installed("fmesher")
  skip_if_not_installed("geometry")
  d <- .native_iid_data(n = 186L, levels = 93L)
  d$slide <- factor(rep(paste0("slide", seq_len(93L)), length.out = nrow(d)),
                    levels = c("unused", paste0("slide", seq_len(93L))))
  model1 <- inlaST.set(
    response ~ s(slide, bs = "re"), d, family = gaussian(),
    mesh = .native_iid_mesh3d(), kappa = .8, coordinates = c("x", "y", "z"),
    control = list(fixed_precision = c(1.1, 1.7), gaussian_precision = 4)
  )
  expect_identical(dim(model1$inla_spec$nuisance_design), c(nrow(d), 94L))
  expect_identical(ncol(model1$inla_spec$fixed$X), 1L)
  expect_identical(ncol(model1$inla_spec$random[[2L]]$A), 93L)

  model6 <- inlaST.set(
    response ~ cell + age + s(slide, bs = "re") + s(batch, bs = "re"),
    d, family = gaussian(),
    mesh = .native_iid_mesh3d(), kappa = .8, coordinates = c("x", "y", "z"),
    control = list(fixed_precision = c(1.1, 1.7, 2.2), gaussian_precision = 4)
  )
  spec <- model6$inla_spec
  p <- ncol(spec$fixed$X)
  q <- ncol(spec$random[[1L]]$A)
  expect_identical(p, 6L)
  expect_identical(dim(spec$nuisance_design), c(nrow(d), 103L))
  expect_identical(vapply(spec$random, `[[`, integer(1L), "sp_index"), 1:3)
  expect_identical(vapply(spec$random, `[[`, character(1L), "name"),
                   c("global", "slide", "batch"))
  expect_identical(spec$geometry_sp_length, 3L)
  expect_identical(
    spec$nuisance_index,
    c(seq_len(p), p + q + seq_len(93L), p + q + 93L + seq_len(4L))
  )
  expect_equal(spec$nuisance_design,
               cbind(spec$fixed$X, as.matrix(spec$random[[2L]]$A),
                     as.matrix(spec$random[[3L]]$A)))
  expect_identical(vapply(spec$nuisance_map, `[[`, integer(1L), "full_column"),
                   spec$nuisance_index)
  expect_identical(vapply(spec$nuisance_map[seq_len(p)], `[[`, character(1L), "source"),
                   rep("fixed", p))
  expect_identical(vapply(spec$nuisance_map[-seq_len(p)], `[[`, character(1L), "source"),
                   rep("random", 97L))
  expect_identical(
    vapply(spec$nuisance_map[-seq_len(p)], `[[`, integer(1L), "block"),
    c(rep(2L, 93L), rep(3L, 4L))
  )
  expect_true(mgcvST:::.inlast_sparse_score_capability(model6)$eligible)
})

test_that("native iid modes follow the nuisance coefficient map", {
  spec <- list(
    nuisance_design = matrix(0, 2L, 4L,
      dimnames = list(NULL, c("fixed1", "fixed2", "iid1", "iid2"))),
    nuisance_map = list(
      list(source = "fixed", block = NA_integer_, index = 1L),
      list(source = "fixed", block = NA_integer_, index = 2L),
      list(source = "random", block = 2L, index = 1L),
      list(source = "random", block = 2L, index = 2L)
    )
  )
  fit <- list(
    fixed_mode = c(0.5, -0.25),
    random_mode = list(global = c(9, 8), slide = c(0.1, 0.3))
  )
  expect_equal(
    mgcvST:::.inlast_nuisance_mode(fit, spec),
    c(fixed1 = 0.5, fixed2 = -0.25, iid1 = 0.1, iid2 = 0.3)
  )
})

test_that("native iid modes survive inlaST estimate compaction", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  skip_if_not_installed("fmesher")
  skip_if_not_installed("geometry")
  d <- .native_iid_data(n = 48L, levels = 3L, seed = 4911L)
  effect <- c(-0.4, 0.1, 0.35)[match(d$slide, paste0("slide", 1:3))]
  set.seed(4912L)
  y <- 0.3 + effect + d$offset0 + rnorm(nrow(d), sd = 0.25)
  model <- inlaST.set(
    response ~ offset(offset0) + s(slide, bs = "re"), d,
    family = gaussian(), mesh = .native_iid_mesh3d(), kappa = 0.8,
    coordinates = c("x", "y", "z")
  )
  control <- list(
    fixed_precision = c(2, 3), gaussian_precision = 16,
    num_threads = 1L
  )
  fit <- inlaST.estimate(
    matrix(y, nrow = 1L, dimnames = list("gene", NULL)), model,
    BPPARAM = BiocParallel::SerialParam(), control = control, threads = 1L
  )
  direct <- mgcvST:::.inlast_fit_feature(
    model$inla_spec, y, offset = model$offset, control = control
  )
  expected <- unname(c(direct$fixed_mode, direct$random_mode[[2L]]))
  expect_equal(as.numeric(fit$nuisance_coefficients[, 1L]), expected,
               tolerance = 1e-8)
  working <- mgcvST:::.inlast_working_state(fit, 1L, threads = 1L)
  expect_equal(as.numeric(working$eta), direct$eta, tolerance = 1e-8)
})

test_that("native iid rejects unsupported grouping structures", {
  d <- .native_iid_data()
  d$numeric_group <- seq_len(nrow(d))
  d$slide2 <- factor(rep(c("a", "b"), length.out = nrow(d)))
  expect_error(mgcvST:::.inlast_native_design(
    response ~ s(slide, bs = "rw1"), d), "nuisance term")
  expect_error(mgcvST:::.inlast_native_design(
    response ~ s(slide, bs = "rw2"), d), "nuisance term")
  expect_error(mgcvST:::.inlast_native_design(
    response ~ s(numeric_group, bs = "re"), d), "factor or character")
  expect_error(mgcvST:::.inlast_native_design(
    response ~ s(slide, bs = "re", by = age), d), "does not support by")
  expect_length(mgcvST:::.inlast_native_design(
    response ~ s(slide, bs = "re") + s(slide2, bs = "re"), d
  )$nuisance, 2L)
  ignored <- mgcvST:::.inlast_native_design(
    response ~ s(slide, bs = "re", k = 2, xt = list(dummy = TRUE),
                 id = "shared", fx = TRUE), d
  )
  expect_identical(ignored$nuisance[[1L]]$levels, paste0("slide", 1:3))
  repeated_re <- mgcvST:::.inlast_native_design(
    response ~ s(slide, bs = "re", k = 2) +
      s(slide, bs = "re", k = 3, id = "second"), d
  )
  expect_identical(
    vapply(repeated_re$nuisance, `[[`, character(1L), "name"),
    c("slide", "slide.1")
  )
  expect_equal(repeated_re$nuisance[[1L]]$Z,
               repeated_re$nuisance[[2L]]$Z)
  expect_error(mgcvST:::.inlast_native_design(
    response ~ s(age, bs = "gp"), d), "full-rank-penalty")
  expect_error(mgcvST:::.inlast_native_design(
    response ~ s(age, bs = "gp", m = -2, by = offset0), d),
    "does not support by")
  expect_error(mgcvST:::.inlast_native_design(
    response ~ s(age, bs = "gp", m = -2, fx = TRUE), d),
    "full-rank-penalty")
  repeated <- mgcvST:::.inlast_native_design(
    response ~ s(age, bs = "gp", k = 8, m = -2) +
      s(age, bs = "gp", k = 10, m = c(-2, .3)), d
  )
  expect_identical(vapply(repeated$nuisance, `[[`, character(1L), "name"),
                   c("s(age)", "s(age).1"))
  expect_identical(anyDuplicated(vapply(
    repeated$nuisance, `[[`, character(1L), "label"
  )), 0L)
})

test_that("the GP bridge projects the intercept and whitens the penalty", {
  d <- .native_iid_data(n = 40L, seed = 2510L)
  d$response <- 0
  f <- response ~ s(age, bs = "gp", k = 10, m = -2)
  split <- mgcv::interpret.gam(f)
  mf <- stats::model.frame(split$fake.formula, d, na.action = stats::na.fail)
  raw <- mgcv::smoothCon(
    split$smooth.spec[[1L]], mf, absorb.cons = FALSE, scale.penalty = TRUE
  )[[1L]]
  bridge <- mgcvST:::.inlast_native_gp_bridge(
    split$smooth.spec[[1L]], mf, "s(age, bs = 'gp')"
  )[[1L]]

  g <- colMeans(raw$X)
  K <- qr.Q(qr(matrix(g, ncol = 1L)), complete = TRUE)[, -1L, drop = FALSE]
  S <- crossprod(K, raw$S[[1L]] %*% K)
  S <- (S + t(S)) / 2
  R <- chol(S)
  Ri <- backsolve(R, diag(ncol(R)))

  expect_identical(ncol(bridge$Z), ncol(raw$X) - 1L)
  expect_lt(max(abs(colMeans(bridge$Z))), 1e-12)
  expect_equal(bridge$projection, K, tolerance = 1e-12)
  expect_equal(bridge$centred_penalty, S, tolerance = 1e-12)
  expect_equal(unname(bridge$Z), unname(raw$X %*% K %*% Ri),
               tolerance = 1e-12)
  expect_equal(crossprod(Ri, S %*% Ri), diag(ncol(R)), tolerance = 1e-11)
})

test_that("native setup registers re and whitened GP blocks in formula order", {
  skip_if_not_installed("fmesher")
  skip_if_not_installed("geometry")
  d <- .native_iid_data(n = 40L, seed = 2511L)
  model <- inlaST.set(
    response ~ s(slide, bs = "re") +
      s(age, bs = "gp", k = 10, m = -2),
    d, family = gaussian(), mesh = .native_iid_mesh3d(), kappa = .8,
    coordinates = c("x", "y", "z")
  )
  random <- model$inla_spec$random

  expect_identical(vapply(random, `[[`, character(1L), "name"),
                   c("global", "slide", "s(age)"))
  expect_identical(vapply(random, `[[`, integer(1L), "sp_index"), 1:3)
  expect_identical(random[[2L]]$nuisance_type, "re")
  expect_identical(random[[3L]]$nuisance_type, "gp")
  expect_equal(as.matrix(random[[2L]]$Q), diag(ncol(random[[2L]]$A)))
  expect_equal(as.matrix(random[[3L]]$Q), diag(ncol(random[[3L]]$A)))
  expect_lt(max(abs(Matrix::colMeans(random[[3L]]$A))), 1e-12)
  expect_equal(
    model$inla_spec$nuisance_design,
    cbind(model$inla_spec$fixed$X, as.matrix(random[[2L]]$A),
          as.matrix(random[[3L]]$A)),
    tolerance = 1e-12
  )
})

test_that("a fixed GP precision reproduces the mgcv REML smooth", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  skip_if_not_installed("fmesher")
  skip_if_not_installed("geometry")
  d <- .native_iid_data(n = 60L, seed = 2512L)
  phi <- .12
  eta <- .3 + sin(2 * pi * d$age) +
    c(-.2, .1, .25)[match(d$slide, paste0("slide", 1:3))] + d$offset0
  set.seed(2513L)
  d$response <- eta + rnorm(nrow(d), sd = sqrt(phi))
  f <- response ~ offset(offset0) + s(slide, bs = "re") +
    s(age, bs = "gp", k = 10, m = -2)
  reference <- mgcv::gam(f, data = d, method = "REML", scale = phi)
  model <- inlaST.set(
    f, d, family = gaussian(), mesh = .native_iid_mesh3d(), kappa = .8,
    coordinates = c("x", "y", "z")
  )
  engine <- mgcvST:::.inlast_fit_feature(
    model$inla_spec, d$response, offset = model$offset,
    control = list(
      fixed_precision = c(1e10, unname(reference$sp) / phi),
      gaussian_precision = 1 / phi
    )
  )
  gp_inla <- as.numeric(
    model$inla_spec$random[[3L]]$A %*% engine$random_mode[[3L]]
  )
  gp_mgcv <- as.numeric(predict(reference, type = "terms")[, "s(age)"])

  expect_equal(engine$smoothing_parameters[2:3], unname(reference$sp),
               tolerance = 1e-10)
  expect_equal(gp_inla, gp_mgcv, tolerance = 2e-4)
})

test_that("estimated re and GP smooths agree with mgcv REML", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  skip_if_not_installed("geometry")
  set.seed(2520L)
  n <- 160L
  vertices <- as.matrix(expand.grid(
    x = seq(0, 1, length.out = 5L), y = seq(0, 1, length.out = 5L)
  ))
  mesh <- list(loc = vertices,
               graph = list(tv = geometry::delaunayn(vertices)))
  d <- data.frame(
    x = runif(n, .03, .97), y = runif(n, .03, .97),
    age = runif(n, -1, 1),
    slide = factor(rep(paste0("s", 1:4), length.out = n)),
    offset0 = seq(-.1, .1, length.out = n)
  )
  basis <- spde_basis(
    mesh, as.matrix(d[c("x", "y")]), kappa = 1.1,
    project_intercept = TRUE
  )
  basis$component <- basis$score.component <- "global"
  phi <- .12
  slide_effect <- c(-.25, .05, .12, .3)[as.integer(d$slide)]
  eta <- .4 + .45 * sin(2 * pi * d$x) * cos(2 * pi * d$y) +
    .65 * sin(pi * d$age) + slide_effect + d$offset0
  d$response <- eta + rnorm(n, sd = sqrt(phi))

  reference <- mgcv::gam(
    response ~ offset(offset0) +
      s(x, y, bs = "spde", xt = basis) +
      s(slide, bs = "re") + s(age, bs = "gp", k = 10, m = -2),
    data = d, method = "REML", scale = phi
  )
  model <- inlaST.set(
    response ~ offset(offset0) +
      s(slide, bs = "re") + s(age, bs = "gp", k = 10, m = -2),
    d, family = gaussian(), mesh = mesh, kappa = 1.1,
    coordinates = c("x", "y")
  )
  engine <- mgcvST:::.inlast_fit_feature(
    model$inla_spec, d$response, offset = model$offset,
    control = list(gaussian_precision = 1 / phi)
  )
  gp_inla <- as.numeric(
    model$inla_spec$random[[3L]]$A %*% engine$random_mode[[3L]]
  )
  gp_mgcv <- as.numeric(predict(reference, type = "terms")[, "s(age)"])

  # The native spatial block uses raw FEM scaling while mgcv reports its
  # rescaled SPDE multiplier. The re and GP penalties have identical scaling.
  expect_lt(max(abs(log(
    engine$smoothing_parameters[2:3] / unname(reference$sp[2:3])
  ))), 5e-3)
  expect_equal(gp_inla, gp_mgcv, tolerance = 1e-3)
})

