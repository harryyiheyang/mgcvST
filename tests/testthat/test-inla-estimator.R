.inlast_fixture <- function(n = 48L, seed = 804L) {
  skip_if_not_installed("INLA")
  skip_if_not_installed("geometry")
  set.seed(seed)
  vertices <- as.matrix(expand.grid(
    x = seq(0, 1, length.out = 4L),
    y = seq(0, 1, length.out = 4L)
  ))
  mesh <- list(loc = vertices, graph = list(tv = geometry::delaunayn(vertices)))
  data <- data.frame(
    x = runif(n, 0.02, 0.98),
    y = runif(n, 0.02, 0.98),
    z = seq(-1, 1, length.out = n),
    exposure = runif(n, 0.8, 1.4)
  )
  data$offset0 <- log(data$exposure)
  basis <- spde_basis(
    mesh, as.matrix(data[c("x", "y")]), kappa = 0.7,
    project_intercept = TRUE
  )
  # The estimator must rebuild A and the observation constraint at the data
  # supplied to inlaST.set(), rather than silently using basis$coordinates.
  data$x <- data$x + 1e-5
  data$y <- data$y - 1e-5
  list(data = data, basis = basis, mesh = mesh)
}

test_that("INLA defaults estimated precision and NB size hyperparameters to flat", {
  control <- mgcvST:::.inlast_control()
  fields <- c("precision_prior", "gaussian_precision_prior", "nb_size_prior")
  for (field in fields) {
    expect_identical(control[[field]]$prior, "flat")
    expect_length(control[[field]]$param, 0L)
    expect_identical(control[[field]]$initial, 0)

    bad <- list()
    bad[[field]] <- list(
      prior = "not-an-INLA-prior", param = numeric(), initial = 0
    )
    expect_error(mgcvST:::.inlast_control(bad), "Unknown INLA")

    flat <- list()
    flat[[field]] <- list(
      prior = "flat", param = numeric(), initial = -1.25
    )
    flat_control <- mgcvST:::.inlast_control(flat)
    expect_identical(flat_control[[field]]$prior, "flat")
    expect_length(flat_control[[field]]$param, 0L)
    expect_identical(flat_control[[field]]$initial, -1.25)

    bad_flat <- list()
    bad_flat[[field]] <- list(
      prior = "flat", param = 0, initial = 0
    )
    expect_error(mgcvST:::.inlast_control(bad_flat), "flat")
  }
  normal_nb <- mgcvST:::.inlast_control(list(nb_size_prior = list(
    prior = "normal", param = c(0, 1 / 9), initial = 0
  )))
  expect_identical(normal_nb$nb_size_prior$prior, "normal")
  expect_equal(normal_nb$nb_size_prior$param, c(0, 1 / 9), tolerance = 0)
  expect_identical(mgcvST:::.inlast_control(control), control)
})

test_that("default Gaussian variance agrees with a restricted likelihood oracle", {
  skip_on_cran()
  f <- .inlast_fixture(seed = 2401L)
  d <- f$data
  n <- nrow(d)
  y <- 0.6 + 0.4 * d$z + d$offset0 +
    0.35 * sin(2 * pi * d$x) - 0.25 * cos(2 * pi * d$y) +
    rnorm(n, sd = 0.3)
  model <- inlaST.set(
    response ~ z + offset(offset0), d, f$basis, family = gaussian(),
    control = list(fixed_precision = 2.5)
  )
  fit <- inlaST.estimate(
    matrix(y, nrow = 1L, dimnames = list("gaussian_flat", NULL)), model,
    BPPARAM = BiocParallel::SerialParam()
  )
  target <- fit$geometry$target[["global"]]
  X <- fit$geometry$X
  B <- fit$geometry$smooth[[target]]$B
  Q <- fit$geometry$smooth[[target]]$penalties[[1L]]
  G <- B %*% solve(2.5 * Q, t(B))
  residual <- y - model$offset
  # Independently integrate the spatial field and unpenalized fixed effects.
  objective <- function(log_phi) {
    V <- exp(log_phi) * diag(n) + G
    Vi <- solve(V)
    ViX <- Vi %*% X
    XtViX <- crossprod(X, ViX)
    P <- Vi - ViX %*% solve(XtViX, t(ViX))
    as.numeric(determinant(V, logarithm = TRUE)$modulus +
      determinant(XtViX, logarithm = TRUE)$modulus +
      crossprod(residual, P %*% residual))
  }
  exact <- optimize(objective, interval = log(c(0.01, 2)), tol = 1e-10)
  expect_true(fit$diagnostics$converged)
  expect_lt(abs(log(fit$dispersion[[1L]]) - exact$minimum), 2e-3)
})

test_that("NB joint mode satisfies the fixed-hyperparameter penalized score equations", {
  skip_on_cran()
  f <- .inlast_fixture(n = 60L, seed = 808L)
  d <- f$data
  eta <- 0.35 + 0.25 * d$z + d$offset0 + 0.3 * sin(2 * pi * d$x)
  set.seed(809)
  y <- rnbinom(nrow(d), mu = exp(eta), size = 4)
  model <- inlaST.set(
    response ~ z + offset(offset0), d, f$basis,
    family = mgcv::nb(theta = 4)
  )
  tau <- 2
  engine <- mgcvST:::.inlast_fit_feature(
    model$inla_spec, y, offset = model$offset,
    control = list(
      fixed_precision = tau, nb_size = 4
    )
  )

  j <- model$geometry$target[["global"]]
  X <- model$geometry$X
  B <- model$geometry$smooth[[j]]$B
  Q <- model$geometry$smooth[[j]]$penalties[[1L]]
  beta <- engine$fixed_mode
  u <- engine$coefficients$global
  expect_equal(
    engine$eta, as.numeric(model$offset + X %*% beta + B %*% u),
    tolerance = 2e-7
  )

  likelihood_score <- 4 * (y - engine$mu) / (4 + engine$mu)
  fixed_residual <- as.numeric(crossprod(X, likelihood_score))
  random_data_score <- as.numeric(crossprod(B, likelihood_score))
  random_penalty_score <- as.numeric(tau * Q %*% u)
  random_residual <- random_data_score - random_penalty_score
  random_scale <- max(1, abs(random_data_score), abs(random_penalty_score))
  expect_lt(max(abs(fixed_residual)), 2e-4)
  expect_lt(max(abs(random_residual)) / random_scale, 2e-5)
  expect_lt(abs(mean(B %*% u)), 1e-10)
})

test_that("INLA nuisance covariance uses expected curvature without configs", {
  skip_on_cran()
  f <- .inlast_fixture(n = 64L, seed = 918L)
  d <- f$data
  set.seed(919L)
  eta <- -1.1 + 0.2 * d$z + d$offset0 + 0.25 * sin(2 * pi * d$x)
  y <- rnbinom(nrow(d), mu = exp(eta), size = 2)
  model <- inlaST.set(
    response ~ z + offset(offset0), d, f$basis,
    family = mgcv::nb(theta = 2)
  )
  engine <- mgcvST:::.inlast_fit_feature(
    model$inla_spec, y, offset = model$offset,
    control = list(
      fixed_precision = 1.7, nb_size = 2, keep_fit = TRUE
    ), diagnostics = TRUE
  )

  X <- model$geometry$X
  B <- model$geometry$smooth[[1L]]$B
  Q <- model$geometry$smooth[[1L]]$penalties[[1L]]
  T <- cbind(X, B)
  H <- crossprod(T / sqrt(engine$working_variance)) +
    as.matrix(Matrix::bdiag(matrix(0, ncol(X), ncol(X)), 1.7 * Q))
  oracle <- solve(H)[seq_len(ncol(X)), seq_len(ncol(X)), drop = FALSE]
  expect_equal(unname(engine$nuisance_covariance), unname(oracle), tolerance = 1e-9)
  expect_identical(engine$expected_nuisance_covariance, engine$nuisance_covariance)
  expect_null(engine$inla$misc$configs)
})

test_that("INLA null specifications reindex fixed and iid nuisance coefficients", {
  A <- Matrix::Matrix(c(1, 0, 0, 1), 2, 2, sparse = TRUE)
  spec <- list(
    fixed = list(X = matrix(1, 2, 1), names = "(Intercept)"),
    random = list(
      list(name = "global", A = A, Q = Matrix::Diagonal(2), target = TRUE,
           kind = "spde", sp_index = 1L),
      list(name = "batch", A = A, Q = Matrix::Diagonal(2), target = FALSE,
           kind = "nuisance", subtype = "iid", sp_index = 2L)
    ), nuisance_index = c(1L, 4L, 5L)
  )
  null <- mgcvST:::.inlast_null_spec(spec)
  expect_length(null$random, 1L)
  expect_identical(null$random[[1L]]$name, "batch")
  expect_identical(null$random[[1L]]$sp_index, 2L)
  expect_identical(null$nuisance_index, 1:3)
  expect_identical(
    mgcvST:::.inlast_null_control(list(fixed_precision = c(2, 3)), spec)$fixed_precision,
    3
  )
  expect_identical(
    mgcvST:::.inlast_null_control(list(fixed_precision = 2), spec)$fixed_precision,
    2
  )
  fixed_only <- spec
  fixed_only$random <- spec$random[1L]
  fixed_only <- mgcvST:::.inlast_null_spec(fixed_only)
  fixed_only$family <- "gaussian"
  expect_length(mgcvST:::.inlast_validate_spec(fixed_only, c(1, 2), c(0, 0))$random, 0L)
})

test_that("constrained SPDE tau mode uses the m-minus-one normalizer", {
  skip_on_cran()
  f <- .inlast_fixture(n = 80L, seed = 910L)
  d <- f$data
  phi <- 0.16
  set.seed(911L)
  y <- 0.4 + 0.3 * d$z + d$offset0 +
    0.2 * sin(2 * pi * d$x) + rnorm(nrow(d), sd = sqrt(phi))
  model <- inlaST.set(
    response ~ z + offset(offset0), d, f$basis,
    family = gaussian()
  )
  engine <- mgcvST:::.inlast_fit_feature(
    model$inla_spec, y, offset = model$offset,
    control = list(
      gaussian_precision = 1 / phi,
      precision_prior = list(
        prior = "normal", param = c(0, 1 / 9), initial = 0
      )
    )
  )

  target <- model$geometry$target[["global"]]
  X <- model$geometry$X
  B <- model$geometry$smooth[[target]]$B
  Q <- model$geometry$smooth[[target]]$penalties[[1L]]
  design <- cbind(X, B)
  q <- ncol(B)
  expect_identical(q, ncol(model$inla_spec$random[[target]]$A) - 1L)
  residual <- y - model$offset
  base_precision <- crossprod(design) / phi
  rhs <- crossprod(design, residual) / phi
  response_quadratic <- sum(residual^2) / phi
  spatial <- ncol(X) + seq_len(q)

  # Integrate the projected spatial coefficients and flat fixed effects.
  # The constrained prior contributes (m - 1) / 2 * log(tau), while the
  # required log(tau) ~ N(0, 3^2) prior contributes -log(tau)^2 / 18.
  log_posterior <- function(log_tau) {
    precision <- base_precision
    precision[spatial, spatial] <-
      precision[spatial, spatial] + exp(log_tau) * Q
    factor <- chol(precision)
    log_determinant <- 2 * sum(log(diag(factor)))
    solution <- backsolve(
      factor, forwardsolve(t(factor), rhs)
    )
    q / 2 * log_tau - 0.5 * log_determinant -
      0.5 * (response_quadratic - sum(rhs * solution)) -
      log_tau^2 / 18
  }
  exact <- stats::optimize(
    function(log_tau) -log_posterior(log_tau),
    interval = c(-12, 18), tol = 1e-10
  )

  expect_true(is.finite(engine$tau[["global"]]))
  expect_lt(
    abs(log(engine$tau[["global"]]) - exact$minimum),
    2e-3
  )
})
