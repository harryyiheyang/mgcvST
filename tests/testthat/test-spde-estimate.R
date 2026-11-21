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
  # The returned point is a stationary point of the criterion: the gradient
  # is reported, and no point of a grid around it is higher.
  expect_true(all(abs(fit$grad_log_tau[1:2]) < .1 & abs(fit$grad_log_alpha[1:2]) < .1))
  expect_true(is.na(fit$grad_log_alpha[3L]))
  for (g in 1:2) {
    grid <- expand.grid(r = log(fit$tau[g]) + seq(-.5, .5, length.out = 11L),
                        t = -log(fit$size[g]) + seq(-.5, .5, length.out = 11L))
    L <- spde.laml(f$Y[g, ], f$X, f$A, f$Q, grid$r, exp(grid$t), offset = f$off)
    expect_lt(max(L) - fit$laml[g], 1e-3)
    expect_equal(L[61L], fit$laml[g], tolerance = 1e-6)
  }
  expect_equal(spde.laml(f$Y[1L, ], f$X, f$A, f$Q, log(only$tau[1L]), 0, offset = f$off),
               only$laml[1L], tolerance = 1e-6)
  expect_error(spde.laml(f$Y[1L, ], f$X, f$A, f$Q, 0, -1), "non-negative")
  # The start of the size search does not move the optimum.
  started <- spde.estimate(f$Y, f$X, f$A, f$Q, offset = f$off, size_start = c(.5, 200, 3))
  expect_identical(started$family, fit$family)
  expect_equal(started$tau, fit$tau, tolerance = 1e-5)
  expect_equal(started$size[1:2], fit$size[1:2], tolerance = 1e-5)
  expect_equal(started$laml, fit$laml, tolerance = 1e-8)
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

test_that("the three stages run and a feature's record is its spde.estimate() fit", {
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
  expect_error(spdeST.set(response ~ z + offset(offset0), d, basis, family = mgcv::nb(),
                          control = list(engine = "inla")), "Unknown control field")
  expect_error(spdeST.set(response ~ z + offset(offset0), d, basis, family = mgcv::nb(),
                          control = list(spde = list(tolerance = 1))),
               "spde.estimate\\(\\) controls")
  model <- spdeST.set(response ~ z + offset(offset0), d, basis, family = mgcv::nb(),
                      diagnostics = TRUE)
  withr::local_options(spdeST.progress = FALSE)
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
  block <- model$spde_spec$random[[1L]]
  direct <- spde.estimate(Y[1:3, ], model$spde_spec$fixed$X, block$A, block$Q,
                          offset = model$offset)
  expect_equal(unname(fit$smoothing_parameters[1:3, block$sp_index]), direct$tau, tolerance = 1e-5)
  expect_equal(unlist(fit$family_parameters[1:2], use.names = FALSE), direct$size[1:2],
               tolerance = 1e-5)
  expect_equal(unname(fit$target_coefficients[, 1:3]), unname(attr(direct, "field")),
               tolerance = 1e-4)
  expect_identical(fit$fit_diagnostics[[1L]]$estimation$engine, "spde.estimate")
  # The size of the null model is the start of the search and is kept.
  expect_equal(unname(fit$diagnostics$null_size[1:2]),
               vapply(fit$null_state[1:2], `[[`, numeric(1L), "family_parameters"),
               ignore_attr = TRUE)
  expect_equal(fit$diagnostics$spatial_size[1:3], c(direct$size[1:2], Inf), tolerance = 1e-4)
  expect_true(is.na(fit$diagnostics$spatial_size[4L]))
  expect_true(all(is.finite(fit$diagnostics$gradient_log_tau[1:3])))
  expect_true(all(is.finite(fit$diagnostics$gradient_log_alpha[1:2])) &&
                is.na(fit$diagnostics$gradient_log_alpha[3L]))
  expect_identical(fit$diagnostics$criterion_name[1L], "spde.estimate Laplace marginal likelihood")

  skip_if_not_installed("arrow")
  res <- spdeST.test(fit, moments = "exact")
  expect_equal(nrow(res$results), 3L)
  expect_true(all(is.finite(res$results$p_value)))
  r <- res$results
  expect_gt(r$score[r$i == 1 & r$j == 2], 0)
  expect_lt(r$score[r$i == 1 & r$j == 3], 0)
})

test_that("a model outside the scope of spde.estimate() is refused and chunks follow the cap", {
  expect_identical(.spdest_control(), list(spde = list()))
  expect_error(.spdest_control(list(control.inla = list(tolerance = 1e-6))), "Unknown control field")
  expect_error(.spde_estimate_check(list(family = "gaussian", random = list(list(kind = "spde")))),
               "Poisson or negative-binomial")
  expect_error(.spde_estimate_check(list(family = "poisson", random = list())), "one SPDE field")
  expect_identical(lengths(.mgcvst_feature_chunks(1:400, NULL, BiocParallel::SerialParam(),
                                                  checkpoint = TRUE, cap = 160L)),
                   c(`1` = 160L, `2` = 160L, `3` = 80L))
})

test_that("the criterion holds when the field can nearly take up the intercept", {
  f <- estimate_fixture()
  skip_if_not_installed("MASS")
  # A precision that hardly penalises a constant field: without the constraint
  # the intercept and the constant would be one direction.
  one <- rep(1, ncol(f$Q))
  q1 <- as.numeric(f$Q %*% one)
  Q <- as.matrix(f$Q) - (1 - 1e-10) * tcrossprod(q1) / sum(q1)
  Q <- methods::as(Matrix::Matrix(0.5 * (Q + t(Q)), sparse = TRUE), "generalMatrix")
  X <- cbind(f$X, as.numeric(seq_len(f$n) %% 3L == 0L))
  # The same criterion by dense algebra on the null space of the constraint.
  N <- MASS::Null(f$constraint)
  Z <- cbind(X, as.matrix(f$A %*% N))
  p <- ncol(X)
  reference <- function(y, log_tau, alpha) {
    P <- matrix(0, ncol(Z), ncol(Z))
    P[-seq_len(p), -seq_len(p)] <- exp(log_tau) * crossprod(N, as.matrix(Q %*% N))
    loglik <- function(mu) if (alpha == 0) sum(stats::dpois(y, mu, log = TRUE)) else
      sum(stats::dnbinom(y, mu = mu, size = 1 / alpha, log = TRUE))
    objective <- function(th) loglik(exp(f$off + as.numeric(Z %*% th))) - 0.5 * sum(th * (P %*% th))
    th <- c(log(mean(y + 0.1)) - mean(f$off), rep(0, ncol(Z) - 1L))
    for (it in seq_len(200L)) {
      mu <- exp(f$off + as.numeric(Z %*% th))
      shrink <- 1 / (1 + alpha * mu)
      H <- crossprod(Z, mu * (1 + alpha * y) * shrink^2 * Z) + P
      g <- crossprod(Z, (y - mu) * shrink) - P %*% th
      step <- as.numeric(solve(H, g))
      if (0.5 * sum(g * step) < 1e-12) break
      scale <- 1
      while (objective(th + scale * step) < objective(th)) scale <- scale / 2
      th <- th + scale * step
    }
    objective(th) + 0.5 * (ncol(f$A) - 1) * log_tau -
      0.5 * as.numeric(determinant(H)$modulus) - 0.5 * log(sum(f$constraint^2))
  }
  for (log_tau in c(-2, 0, 2)) for (alpha in c(0, 0.3)) {
    expect_equal(spde.laml(f$Y[1L, ], X, f$A, Q, log_tau, alpha, offset = f$off,
                           constraint = f$constraint),
                 reference(f$Y[1L, ], log_tau, alpha), tolerance = 1e-8)
  }
})

test_that("a factor level without counts leaves the criterion of the model without it", {
  f <- estimate_fixture()
  level <- seq_len(f$n) %% 3L
  X <- cbind(1, as.numeric(level == 1L), as.numeric(level == 2L))
  y <- f$Y[1L, ]
  y[level != 0L] <- 0
  keep <- level == 0L
  for (alpha in c(0, 0.3)) {
    log_tau <- seq(-6, 4, by = 2)
    full <- spde.laml(y, X, f$A, f$Q, log_tau, alpha, offset = f$off,
                      constraint = f$constraint)
    # Each point on its own, not warm-started from its neighbour.
    cold <- vapply(log_tau, function(x) spde.laml(
      y, X, f$A, f$Q, x, alpha, offset = f$off, constraint = f$constraint), 0)
    reduced <- spde.laml(y[keep], X[keep, 1L, drop = FALSE], f$A[keep, , drop = FALSE],
                         f$Q, log_tau, alpha, offset = f$off[keep],
                         constraint = f$constraint)
    expect_equal(full, reduced, tolerance = 1e-8)
    expect_equal(cold, reduced, tolerance = 1e-8)
  }
  Y <- methods::as(Matrix::Matrix(rbind(a = y), sparse = TRUE), "CsparseMatrix")
  fit <- spde.estimate(Y, X, f$A, f$Q, offset = f$off, constraint = f$constraint)
  expect_true(fit$converged)
  beta <- attr(fit, "coefficients")[, 1L]
  eta <- f$off + as.numeric(X %*% beta + f$A %*% attr(fit, "field")[, 1L])
  expect_lte(max(exp(eta[!keep])), 1e-11 * (1 + 1e-8))
})
