scan_data <- function(seed, n = 400L, features = 12L) {
  set.seed(seed)
  d <- data.frame(group = factor(sample(letters[1:3], n, TRUE)),
                  x = stats::rnorm(n))
  X <- stats::model.matrix(~ group + x, d)
  offset <- log(stats::runif(n, .5, 2))
  size <- rep(c(.3, 1, 5, 50), length.out = features)
  Y <- t(vapply(seq_len(features), function(j) {
    beta <- c(stats::runif(1, -2, 2), stats::rnorm(3, 0, .3))
    stats::rnbinom(n, mu = exp(offset + X %*% beta), size = size[j])
  }, numeric(n)))
  list(X = X, offset = offset, Y = Y)
}

scan_kernel <- function(data, poisson = rep(FALSE, nrow(data$Y)), ...) {
  Yt <- methods::as(methods::as(methods::as(
    Matrix::Matrix(t(data$Y), sparse = TRUE), "CsparseMatrix"),
    "generalMatrix"), "dMatrix")
  mgcvST:::mgcvst_nb_null_batch_cpp(Yt, data$X, matrix(data$offset),
                                    poisson, ...)
}

test_that("the batched NB null kernel reproduces glm.nb", {
  skip_if_not_installed("MASS")
  data <- scan_data(9045L)
  z <- scan_kernel(data, threads = 2L)
  expect_true(all(is.na(z$error)))
  expect_true(all(z$converged))
  for (j in seq_len(nrow(data$Y))) {
    ref <- tryCatch(suppressWarnings(MASS::glm.nb(
      data$Y[j, ] ~ data$X - 1 + offset(data$offset),
      control = stats::glm.control(epsilon = 1e-12, maxit = 200))),
      error = function(e) NULL)
    if (is.null(ref)) {
      # No overdispersion: glm.nb fails, and the fit is the Poisson end point.
      pois <- stats::glm(data$Y[j, ] ~ data$X - 1 + offset(data$offset),
                         family = stats::poisson())
      expect_lte(z$score0[j], 0)
      expect_identical(z$alpha[j], 0)
      expect_identical(z$size[j], Inf)
      expect_equal(z$beta[, j], unname(stats::coef(pois)), tolerance = 1e-6)
      next
    }
    expect_equal(z$beta[, j], unname(stats::coef(ref)), tolerance = 1e-6)
    expect_equal(z$size[j], ref$theta, tolerance = 1e-5)
    expect_gt(z$score0[j], 0)
    expect_equal(z$loglik[j], as.numeric(stats::logLik(ref)), tolerance = 1e-8)
  }
})

test_that("the Poisson route and a fixed size reproduce glm", {
  data <- scan_data(9046L, features = 4L)
  z <- scan_kernel(data, poisson = rep(TRUE, 4L))
  for (j in 1:4) {
    ref <- stats::glm(data$Y[j, ] ~ data$X - 1 + offset(data$offset),
                      family = stats::poisson(),
                      control = stats::glm.control(epsilon = 1e-12))
    expect_equal(z$beta[, j], unname(stats::coef(ref)), tolerance = 1e-6)
    expect_equal(z$loglik[j], as.numeric(stats::logLik(ref)), tolerance = 1e-8)
    expect_identical(z$size[j], Inf)
    mu <- stats::fitted(ref)
    expect_equal(z$score0[j], sum((data$Y[j, ] - mu)^2 - data$Y[j, ]) / 2,
                 tolerance = 1e-8)
    expect_equal(z$phi[j], sum((data$Y[j, ] - mu)^2 / mu) / ref$df.residual,
                 tolerance = 1e-8)
  }
  skip_if_not_installed("MASS")
  z <- scan_kernel(data, alpha_fixed = 1 / 2)
  for (j in 1:4) {
    ref <- stats::glm(data$Y[j, ] ~ data$X - 1 + offset(data$offset),
                      family = MASS::negative.binomial(2),
                      control = stats::glm.control(epsilon = 1e-12))
    expect_equal(z$beta[, j], unname(stats::coef(ref)), tolerance = 1e-6)
    expect_identical(z$size[j], 2)
    expect_equal(z$loglik[j], as.numeric(stats::logLik(ref)) + 0, tolerance = 1e-8)
  }
})

test_that("a feature that cannot be fitted reports an error alone", {
  data <- scan_data(9047L, features = 3L)
  data$X <- cbind(data$X, dup = data$X[, "x"])
  z <- scan_kernel(data)
  expect_true(all(!is.na(z$error)))
  expect_true(all(is.na(z$size)))
})

test_that("large counts use the exact tail sums", {
  skip_if_not_installed("MASS")
  set.seed(9049L)
  n <- 300L
  X <- cbind(1, stats::rnorm(n))
  y <- stats::rnbinom(n, mu = exp(6 + .5 * X[, 2]), size = 4)
  expect_gt(max(y), 1000)
  data <- list(X = X, offset = numeric(n), Y = matrix(y, 1L))
  z <- scan_kernel(data)
  ref <- MASS::glm.nb(y ~ X - 1, control = stats::glm.control(epsilon = 1e-12))
  expect_equal(z$size, ref$theta, tolerance = 1e-6)
  expect_equal(z$loglik, as.numeric(stats::logLik(ref)), tolerance = 1e-8)
  expect_equal(z$beta[, 1L], unname(stats::coef(ref)), tolerance = 1e-7)
})

test_that("the end point is a maximum when the overdispersion score is not positive", {
  set.seed(9050L)
  n <- 500L
  X <- cbind(1, stats::rnorm(n))
  y <- stats::rbinom(n, 20, .1)  # underdispersed counts
  data <- list(X = X, offset = numeric(n), Y = matrix(y, 1L))
  z <- scan_kernel(data)
  expect_lt(z$score0, 0)
  expect_identical(z$alpha, 0)
  for (a in c(1e-6, 1e-3, .1)) {
    fixed <- scan_kernel(data, alpha_fixed = a)
    expect_lt(fixed$loglik, z$loglik)
  }
})

scan_field <- function(n = 300L, seed = 9052L) {
  skip_if_not_installed("geometry")
  set.seed(seed)
  xy <- as.matrix(expand.grid(x = seq(0, 1, length.out = 6),
                              y = seq(0, 1, length.out = 6)))
  mesh <- list(loc = xy, graph = list(tv = geometry::delaunayn(xy)))
  coords <- cbind(stats::runif(n, .01, .99), stats::runif(n, .01, .99))
  basis <- spde_basis(mesh, coords, kappa = .8)
  field <- mgcvST:::.spde_scan_field(basis)
  list(basis = basis, A = field$A, Q = field$Q, coords = coords)
}

# Dense reference of the GLM score test: T = r'P0 A Q^-1 A' P0 r and the
# mixture weights eig(Q^-1/2 A' P0 A Q^-1/2), P0 = W - WX(X'WX)^-1X'W.
scan_reference <- function(y, X, offset, beta, alpha, A, Q) {
  mu <- exp(drop(X %*% beta) + offset)
  w <- mu / (1 + alpha * mu)
  r <- (y - mu) / mu
  W <- diag(w)
  P0 <- W - W %*% X %*% solve(crossprod(X, w * X), t(X) %*% W)
  A <- as.matrix(A)
  S <- solve(as.matrix(Q))
  stat <- drop(t(r) %*% P0 %*% A %*% S %*% t(A) %*% P0 %*% r)
  R <- chol(S)
  lam <- eigen(R %*% t(A) %*% P0 %*% A %*% t(R), symmetric = TRUE)$values
  lam <- lam[lam > 1e-12 * max(lam)]
  list(statistic = stat,
       p_value = CompQuadForm::davies(stat, lam, acc = 1e-8, lim = 1e5)$Qq)
}

test_that("spde.scan is the GLM score test of the raw SPDE field", {
  field <- scan_field()
  n <- nrow(field$coords)
  X <- cbind(1, stats::rnorm(n))
  offset <- log(stats::runif(n, .5, 2))
  u <- as.numeric(field$A %*% as.numeric(t(chol(solve(as.matrix(field$Q)))) %*%
                                   stats::rnorm(ncol(field$A))))
  u <- .6 * (u - mean(u)) / stats::sd(u)
  eta <- offset + .5 + .2 * X[, 2]
  Y <- rbind(stats::rnbinom(n, mu = exp(eta), size = 2),
             stats::rnbinom(n, mu = exp(eta + u), size = 2),
             stats::rpois(n, exp(eta)),
             stats::rnbinom(n, mu = exp(eta + .5 * u), size = 8))
  rownames(Y) <- paste0("g", 1:4)
  out <- spde.scan(Y, X, field$A, field$Q, offset = offset,
                   poisson = c(FALSE, FALSE, TRUE, FALSE), threads = 2L,
                   chunk_size = 3L)
  expect_identical(out$feature_id, rownames(Y))
  expect_identical(out$family[3], "poisson")
  expect_true(all(out$converged))
  beta <- attr(out, "coefficients")
  for (j in 1:4) {
    alpha <- if (out$family[j] == "poisson") 0 else 1 / out$size[j]
    ref <- scan_reference(Y[j, ], X, offset, beta[, j], alpha, field$A, field$Q)
    expect_equal(out$statistic[j], ref$statistic, tolerance = 1e-8)
    expect_equal(out$p_value[j], ref$p_value, tolerance = 1e-6)
  }
  expect_lt(out$p_value[2], 1e-3)
  # The scale of Q does not change the p-values.
  scaled <- spde.scan(Y, X, field$A, 7 * field$Q, offset = offset,
                      poisson = c(FALSE, FALSE, TRUE, FALSE))
  expect_equal(scaled$p_value, out$p_value, tolerance = 1e-8)
  fixed <- spde.scan(Y[1, , drop = FALSE], X, field$A, field$Q, offset = offset,
                     size = 3)
  expect_identical(fixed$size, 3)
})

test_that("mgcvST.estimate and inlaST.estimate share the spde.scan Stage 1", {
  f <- st_fixture(n = 200L)
  set.seed(9051L)
  mu <- exp(1 + f$data$offset0 + .5 * f$data$z)
  Y <- t(vapply(c(.8, 2, 6), function(s) stats::rnbinom(200L, mu = mu, size = s),
                numeric(200L)))
  Y <- rbind(Y, f$Y[2L, ])
  target <- which(vapply(f$model$G$smooth, function(s)
    identical(s$score.component, "global"), logical(1L)))
  expect_true(mgcvST:::.mgcvst_scan_eligible(
    f$model$G, target, f$model$null_formula, f$model$G$family, NULL, list(), list()))
  expect_false(mgcvST:::.mgcvst_scan_eligible(
    f$model$G, target, f$model$null_formula, mgcv::nb(theta = 3), NULL, list(), list()))
  sp <- BiocParallel::SerialParam()
  withr::local_options(mgcvST.progress = FALSE)
  fit <- mgcvST.estimate(Y, f$model, BPPARAM = sp, spatial = "none")
  expect_identical(fit$stage1$null_engine, "spde.scan")

  # The Stage 1 of inlaST.estimate() on the same field and null design.
  raw <- mgcvST:::.inlast_prepare_basis(f$basis, f$basis$coordinates)$raw
  X0 <- as.matrix(f$model$null_X)
  spec <- list(random = list(list(A = raw$A, Q = raw$Q)))
  null_spec <- list(fixed = list(X = X0), random = list(),
                    family = "negative_binomial")
  expect_true(mgcvST:::.inlast_scan_eligible(null_spec))
  payload <- list(index = 1:4, feature_id = as.character(1:4),
                  Y = Y, extra_offset = NULL,
                  poisson = fit$diagnostics$family_used != "negative_binomial")
  inla <- mgcvST:::.inlast_scan_chunk(payload, spec, null_spec,
                                      f$data$offset0, list())
  p_inla <- vapply(inla, function(z) z$marginal$p_value, numeric(1L))
  expect_equal(fit$diagnostics$marginal_p_value, p_inla, tolerance = 1e-5)

  withr::local_options(mgcvST.scan_null = FALSE)
  expect_false(mgcvST:::.inlast_scan_eligible(null_spec))
  gam <- suppressWarnings(mgcvST.estimate(Y, f$model, BPPARAM = sp, spatial = "none"))
  expect_identical(gam$stage1$null_engine, "gam")
  # The GAM null uses the intercept-projected field: close, not equal.
  expect_equal(fit$diagnostics$marginal_p_value,
               gam$diagnostics$marginal_p_value, tolerance = .05)
})
