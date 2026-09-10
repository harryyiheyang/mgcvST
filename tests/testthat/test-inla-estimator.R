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

.inlast_spatial_means <- function(fit, feature = 1L) {
  vapply(fit$score_components, function(component) {
    j <- fit$geometry$target[[component]]
    B <- fit$geometry$smooth[[j]]$B
    coefficient <- fit$smooth_coefficients[[component]][feature, ]
    mean(as.numeric(B %*% coefficient))
  }, numeric(1L))
}

test_that("INLA defaults spatial precision and NB size to flat", {
  control <- mgcvST:::.inlast_control()
  fields <- c("precision_prior", "gaussian_precision_prior", "nb_size_prior")
  for (field in fields) {
    if (field %in% c("precision_prior", "nb_size_prior")) {
      expect_identical(control[[field]]$prior, "flat")
      expect_length(control[[field]]$param, 0L)
    } else {
      expect_identical(control[[field]]$prior, "normal")
      expect_equal(control[[field]]$param, c(0, 1 / 9), tolerance = 0)
    }
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

test_that("inlaST Gaussian mode agrees with a fixed-hyperparameter oracle", {
  skip_on_cran()
  f <- .inlast_fixture()
  d <- f$data
  set.seed(805)
  signal <- 0.35 * sin(2 * pi * d$x) - 0.25 * cos(2 * pi * d$y)
  y <- 0.6 + 0.4 * d$z + d$offset0 + signal + rnorm(nrow(d), sd = 0.3)

  model <- inlaST.set(
    response ~ z + offset(offset0), d, f$basis,
    family = gaussian()
  )
  control <- list(
    fixed_precision = 2.5, gaussian_precision = 1 / 0.09,
    fixed_effect_precision = 0
  )
  engine <- mgcvST:::.inlast_fit_feature(
    model$inla_spec, y, offset = model$offset, control = control
  )
  fit <- inlaST.estimate(
    matrix(y, nrow = 1L, dimnames = list("gaussian_oracle", NULL)),
    model, retain_smooth = TRUE,
    BPPARAM = BiocParallel::SerialParam(),
    control = control
  )

  expect_s3_class(model, "inlaST_model")
  expect_s3_class(fit, "inlaST_fit")
  expect_true(fit$diagnostics$converged)
  expect_true(is.finite(fit$lambda[[1L]]) && fit$lambda[[1L]] > 0)

  target <- fit$geometry$target[["global"]]
  X <- fit$geometry$X
  B <- fit$geometry$smooth[[target]]$B
  Q <- fit$geometry$smooth[[target]]$penalties[[1L]]
  design <- cbind(X, B)
  penalty <- matrix(0, ncol(design), ncol(design))
  q <- ncol(B)
  penalty[ncol(X) + seq_len(q), ncol(X) + seq_len(q)] <- 2.5 * Q
  precision <- crossprod(design) / 0.09 + penalty
  rhs <- crossprod(design, y - model$offset) / 0.09
  oracle <- solve(precision, rhs)
  oracle_covariance <- solve(precision)
  estimated <- c(engine$fixed_mode, engine$coefficients$global)

  expect_equal(engine$lambda[[1L]], 0.09 * 2.5, tolerance = 1e-10)
  expect_equal(unname(estimated), as.numeric(oracle), tolerance = 3e-5)
  expect_equal(
    as.numeric(fit$smooth_coefficients$global[1L, ]),
    tail(as.numeric(oracle), q), tolerance = 3e-5
  )
  expect_equal(
    unname(fit$nuisance_covariance[[1L]]),
    unname(oracle_covariance[seq_len(ncol(X)), seq_len(ncol(X)), drop = FALSE]),
    tolerance = 3e-5
  )
  expect_lt(max(abs(.inlast_spatial_means(fit))), 5e-8)
})

test_that("inlaST enforces observation-mean zero for Poisson and NB fits", {
  skip_on_cran()
  f <- .inlast_fixture()
  d <- f$data
  eta <- 0.4 + 0.3 * d$z + d$offset0 + 0.35 * sin(2 * pi * d$x)
  set.seed(806)
  responses <- list(
    poisson = list(family = poisson(), y = rpois(nrow(d), exp(eta))),
    negative_binomial = list(
      family = mgcv::nb(theta = 4),
      y = rnbinom(nrow(d), mu = exp(eta), size = 4)
    )
  )

  for (case_name in names(responses)) {
    case <- responses[[case_name]]
    if (identical(case_name, "negative_binomial")) {
      local <- spde_basis(
        f$mesh, f$basis$coordinates, kappa = 2,
        project_intercept = TRUE
      )
      model <- inlaST.set(
        response ~ z + offset(offset0), d,
        list(global = f$basis, local = local), family = case$family,
        setting = "global_local"
      )
    } else {
      model <- inlaST.set(
        response ~ z + offset(offset0), d, f$basis,
        family = case$family
      )
    }
    fit <- inlaST.estimate(
      matrix(case$y, nrow = 1L, dimnames = list("feature", NULL)),
      model, retain_smooth = TRUE,
      BPPARAM = BiocParallel::SerialParam(),
      control = list()
    )
    expect_true(fit$diagnostics$converged)
    expect_true(all(is.finite(fit$working_error)))
    expect_true(all(is.finite(fit$working_variance)))
    expect_true(all(fit$working_variance > 0))
    if (identical(case_name, "negative_binomial")) {
      expect_equal(fit$family_parameters[[1L]], 4, tolerance = 1e-10)
    }
    spatial_means <- .inlast_spatial_means(fit)
    expect_identical(names(spatial_means), model$components)
    expect_lt(max(abs(spatial_means)), 5e-8)
  }
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
      fixed_precision = tau, nb_size = 4,
      fixed_effect_precision = 0
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

test_that("INLA nuisance covariance is the constrained conditional posterior block", {
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
      fixed_precision = 1.7, nb_size = 2,
      fixed_effect_precision = 0, keep_fit = TRUE
    ), diagnostics = TRUE
  )

  posterior <- mgcvST:::.inlast_posterior_vp(
    engine$inla, model$inla_spec
  )
  expect_equal(
    engine$nuisance_covariance,
    posterior$nuisance_covariance,
    tolerance = 1e-12
  )
  expect_equal(
    unname(diag(engine$nuisance_covariance)),
    unname(engine$inla$summary.fixed[, "sd"]^2),
    tolerance = 2e-6
  )
  expect_true(is.matrix(engine$expected_nuisance_covariance))
  expect_identical(
    dim(engine$expected_nuisance_covariance),
    dim(engine$nuisance_covariance)
  )

  config <- engine$inla$misc$configs$config[[1L]]
  full <- mgcvST:::.inlast_posterior_covariance_selected(
    engine$inla, seq_len(nrow(config$Q))
  )
  constraint <- engine$inla$misc$configs$constr$A
  expect_lt(max(abs(constraint %*% full$covariance)), 1e-8)

  contents <- engine$inla$misc$configs$contents
  offset <- engine$inla$misc$configs$mnpred
  expected_variance <- rep(NA_real_, nrow(config$Q))
  for (j in seq_along(contents$tag)) {
    tag <- contents$tag[j]
    start <- contents$start[j] - offset
    if (start < 1L || start > length(expected_variance)) next
    index <- start + seq_len(contents$length[j]) - 1L
    if (tag %in% names(engine$inla$summary.random)) {
      expected_variance[index] <- engine$inla$summary.random[[tag]]$sd^2
    } else if (tag %in% rownames(engine$inla$summary.fixed)) {
      expected_variance[index] <- engine$inla$summary.fixed[tag, "sd"]^2
    }
  }
  expect_false(anyNA(expected_variance))
  expect_equal(diag(full$covariance), expected_variance, tolerance = 2e-6)
})

test_that("inlaST compact fits run the existing covariance score path", {
  skip_on_cran()
  f <- .inlast_fixture(n = 54L)
  d <- f$data
  set.seed(807)
  eta <- 0.5 + 0.2 * d$z + d$offset0 + 0.25 * sin(2 * pi * d$x)
  Y <- rbind(
    feature_a = rpois(nrow(d), exp(eta)),
    feature_b = rpois(nrow(d), exp(eta + 0.15 * cos(2 * pi * d$y)))
  )
  model <- inlaST.set(
    response ~ z + offset(offset0), d, f$basis,
    family = poisson()
  )
  fit <- testthat::with_mocked_bindings(
    inlaST.estimate(
      Y, model, BPPARAM = BiocParallel::SerialParam(), control = list(),
      marginal_args = list(method = "liu"), retain_marginal = TRUE
    ),
    gam = function(...) stop("mgcv::gam() was called during INLA estimation"),
    bam = function(...) stop("mgcv::bam() was called during INLA estimation"),
    .package = "mgcv"
  )
  score <- mgcvST.test(
    fit, pairs = matrix(c("feature_a", "feature_b"), nrow = 1L),
    calibration = "liu", BPPARAM = BiocParallel::SerialParam()
  )

  expect_s3_class(score, "mgcvST_test")
  expect_identical(score$discoveries$pairs_tested, 1L)
  expect_true(is.finite(score$results$signed_score))
  expect_true(is.finite(score$results$p_two_sided))
  expect_true(score$results$p_two_sided >= 0 && score$results$p_two_sided <= 1)
  replay <- mgcvST.marginal(
    fit, calibration = "liu", BPPARAM = BiocParallel::SerialParam()
  )
  expect_equal(
    replay$p_value, fit$diagnostics$marginal_p_value,
    tolerance = 1e-10
  )
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
