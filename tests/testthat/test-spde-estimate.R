estimate_fixture <- function(n = 400L, k = 9L, seed = 11L) {
  skip_if_not_installed("geometry")
  set.seed(seed)
  xy <- as.matrix(expand.grid(x = seq(0, 1, length.out = k), y = seq(0, 1, length.out = k)))
  mesh <- list(loc = xy, graph = list(tv = geometry::delaunayn(xy)))
  coords <- cbind(runif(n, .01, .99), runif(n, .01, .99))
  field <- .spde_scan_field(spde_basis(mesh, coords, kappa = 3))
  A <- field$A
  Q <- field$Q
  m <- ncol(A)
  X <- cbind(1, rnorm(n))
  off <- log(runif(n, .5, 2))
  L <- t(chol(solve(as.matrix(Q))))
  sim <- function(size, amp) {
    u <- as.numeric(A %*% (L %*% rnorm(m)))
    u <- amp * (u - mean(u)) / stats::sd(u)
    mu <- exp(off + 1 + .3 * X[, 2L] + u)
    if (is.finite(size)) stats::rnbinom(n, mu = mu, size = size) else stats::rpois(n, mu)
  }
  # The third feature is underdispersed (binomial): its fit is the Poisson
  # end point alpha = 0.
  under <- stats::rbinom(n, 40L, pmin(exp(off + 1 + .3 * X[, 2L]) / 40, .95))
  Y <- rbind(nb2 = sim(2, .6), nb8 = sim(8, .6), under = under)
  list(Y = Y, X = X, A = A, Q = Q, off = off, n = n,
       constraint = as.numeric(Matrix::crossprod(A, rep(1 / n, n))))
}

# The same model in mgcv: the field on the null space of the constraint, its
# penalty N'QN, REML.
estimate_reference <- function(f, y, family) {
  N <- MASS::Null(f$constraint)
  Zc <- as.matrix(f$A %*% N)
  S <- as.matrix(Matrix::crossprod(N, f$Q %*% N))
  x2 <- f$X[, 2L]
  off <- f$off
  mgcv::gam(y ~ x2 + Zc + offset(off), paraPen = list(Zc = list(S)), family = family,
            method = "REML",
            control = mgcv::gam.control(epsilon = 1e-10, newton = list(conv.tol = 1e-10)))
}

test_that("spde.estimate is the REML fit of mgcv::gam", {
  skip_if_not_installed("MASS")
  f <- estimate_fixture()
  fit <- spde.estimate(f$Y, f$X, f$A, f$Q, offset = f$off)
  expect_true(all(fit$converged))
  expect_true(all(is.na(fit$error)))
  eta <- function(g) f$off + as.numeric(f$X %*% attr(fit, "coefficients")[, g]) +
    as.numeric(f$A %*% attr(fit, "field")[, g])
  expect_equal(as.numeric(f$constraint %*% attr(fit, "field")), rep(0, 3L), tolerance = 1e-10)
  for (g in 1:2) {
    ref <- estimate_reference(f, f$Y[g, ], mgcv::nb())
    poisson <- estimate_reference(f, f$Y[g, ], stats::poisson())
    expect_identical(fit$family[g], "negative_binomial")
    expect_equal(fit$tau[g], unname(ref$sp), tolerance = 1e-3)
    expect_equal(fit$size[g], ref$family$getTheta(TRUE), tolerance = 1e-3)
    expect_equal(eta(g), as.numeric(ref$linear.predictors), tolerance = 1e-4)
    # The criterion is -REML up to a constant shared by the two families.
    expect_equal(fit$laml[g] - fit$laml_poisson[g],
                 as.numeric(poisson$gcv.ubre - ref$gcv.ubre), tolerance = 1e-4)
  }
  ref <- suppressWarnings(estimate_reference(f, f$Y[3L, ], stats::poisson()))
  expect_identical(fit$family[3L], "poisson")
  expect_identical(fit$size[3L], Inf)
  expect_equal(fit$laml[3L], fit$laml_poisson[3L], tolerance = 1e-8)
  expect_equal(eta(3L), as.numeric(ref$linear.predictors), tolerance = 1e-4)

  only <- spde.estimate(f$Y, f$X, f$A, f$Q, offset = f$off, poisson = c(TRUE, FALSE, FALSE))
  expect_identical(only$family, c("poisson", "negative_binomial", "poisson"))
  expect_identical(only$tau[2:3], fit$tau[2:3])
  threaded <- spde.estimate(f$Y, f$X, f$A, f$Q, offset = f$off, threads = 2L)
  expect_identical(threaded$tau, fit$tau)
  expect_identical(attr(threaded, "field"), attr(fit, "field"))
})

test_that("spde.estimate reports a feature it cannot fit and keeps the others", {
  f <- estimate_fixture(n = 200L, k = 6L)
  Y <- rbind(f$Y[1L, ], zero = 0)
  fit <- spde.estimate(Y, f$X, f$A, f$Q, offset = f$off)
  expect_identical(fit$error[2L], "All counts are zero.")
  expect_true(is.na(fit$tau[2L]) && is.na(fit$family[2L]) && !fit$converged[2L])
  expect_true(fit$converged[1L])
  # A level with no counts: its coefficient has no finite estimate, the fitted
  # means of the level vanish and the direction is aliased.
  rare <- as.numeric(seq_len(f$n) <= 5L)
  y <- f$Y[1L, ]
  y[rare == 1] <- 0
  fit <- spde.estimate(rbind(y), cbind(f$X, rare), f$A, f$Q, offset = f$off)
  expect_true(is.na(fit$error) && is.finite(fit$tau) && is.finite(fit$laml))
  expect_error(spde.estimate(f$Y, f$X[-1L, ], f$A, f$Q), "one column per row")
  expect_error(spde.estimate(f$Y, f$X, f$A, f$Q, constraint = 0), "constraint must")
})

test_that("engine = 'spde' runs the three stages without INLA", {
  skip_if_not_installed("geometry")
  set.seed(2402L)
  n <- 300L
  vertices <- as.matrix(expand.grid(x = seq(0, 1, length.out = 6L),
                                    y = seq(0, 1, length.out = 6L)))
  mesh <- list(loc = vertices, graph = list(tv = geometry::delaunayn(vertices)))
  d <- data.frame(x = runif(n, 0.02, 0.98), y = runif(n, 0.02, 0.98),
                  z = seq(-1, 1, length.out = n), offset0 = log(runif(n, 0.8, 1.4)))
  basis <- spde_basis(mesh, as.matrix(d[c("x", "y")]), kappa = 2, project_intercept = TRUE)
  wave <- sin(2 * pi * d$x) * cos(2 * pi * d$y)
  mu <- exp(1 + d$offset0 + .3 * d$z)
  Y <- rbind(g1 = stats::rnbinom(n, mu = mu * exp(.8 * wave), size = 3),
             g2 = stats::rnbinom(n, mu = mu * exp(.8 * wave), size = 6),
             g3 = stats::rbinom(n, 40L, pmin(mu * exp(-.8 * wave) / 40, .95)),
             g4 = stats::rnbinom(n, mu = mu, size = 2))
  expect_error(inlaST.set(response ~ z + offset(offset0), d, basis, family = mgcv::nb(),
                          control = list(engine = "mgcv")), "'inla' or 'spde'")
  expect_error(inlaST.set(response ~ z + offset(offset0), d, basis, family = mgcv::nb(),
                          control = list(engine = "spde", spde = list(tolerance = 1))),
               "spde.estimate\\(\\) controls")
  model <- inlaST.set(response ~ z + offset(offset0), d, basis, family = mgcv::nb(),
                      control = list(engine = "spde"), diagnostics = TRUE)
  withr::local_options(mgcvST.progress = FALSE)
  sp <- BiocParallel::SerialParam()
  scan <- spdeST.scan(Y, model, BPPARAM = sp)
  fit <- spdeST.estimate(scan, Y, features = c("g1", "g2", "g3"), BPPARAM = sp)
  expect_identical(fit$diagnostics$spatial_fitted, c(TRUE, TRUE, TRUE, FALSE))
  expect_true(all(fit$diagnostics$converged[1:3]))
  expect_identical(fit$diagnostics$family_used[1:3],
                   c("negative_binomial", "negative_binomial", "poisson"))
  expect_identical(fit$feature_family[1:3], fit$diagnostics$family_used[1:3])
  expect_length(fit$family_parameters[[3L]], 0L)

  # The record of a feature is the direct spde.estimate() fit of the block.
  block <- model$inla_spec$random[[1L]]
  direct <- spde.estimate(Y[1:3, ], model$inla_spec$fixed$X, block$A, block$Q,
                          offset = model$offset)
  expect_equal(unname(fit$smoothing_parameters[1:3, block$sp_index]), direct$tau, tolerance = 1e-12)
  expect_equal(unlist(fit$family_parameters[1:2], use.names = FALSE), direct$size[1:2],
               tolerance = 1e-12)
  expect_equal(unname(fit$target_coefficients[, 1:3]), unname(attr(direct, "field")),
               tolerance = 1e-12)
  expect_identical(fit$inla_diagnostics[[1L]]$estimation$engine, "spde.estimate")

  skip_if_not_installed("arrow")
  res <- spdeST.test(fit, moments = "exact")
  expect_equal(nrow(res$results), 3L)
  expect_true(all(is.finite(res$results$p_value)))
  r <- res$results
  expect_gt(r$score[r$i == 1 & r$j == 2], 0)
  expect_lt(r$score[r$i == 1 & r$j == 3], 0)
})
