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
    fixed_precision = 2.5, gaussian_precision = 1 / 0.09
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
    model <- inlaST.set(
      response ~ z + offset(offset0), d, f$basis,
      family = case$family
    )
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
      retain_marginal = TRUE
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

test_that("inlaST scores the null fit before the full fit", {
  skip_if_not_installed("INLA")
  skip_if_not_installed("geometry")
  f <- .inlast_fixture(n = 24L, seed = 1821L)
  model <- inlaST.set(response ~ z + offset(offset0), f$data, f$basis,
                      family = mgcv::nb(theta = 4))
  Y <- rbind(g1 = rep(1, nrow(f$data)), g2 = rep(2, nrow(f$data)))
  calls <- new.env(parent = emptyenv())
  calls$event <- character()
  calls$family <- list()
  fake_fit <- function(spec, y, phase) {
    random_names <- vapply(spec$random, `[[`, character(1L), "name")
    list(
      working_error = rep(if (identical(phase, "null")) y[1L] else 20, length(y)),
      working_variance = rep(if (identical(phase, "null")) y[1L] + 2 else 30, length(y)),
      dispersion = 1, family_parameters = 4, smoothing_parameters = rep(1, 1L),
      nuisance_covariance = diag(ncol(spec$nuisance_design)),
      expected_nuisance_covariance = NULL, converged = TRUE,
      log_marginal_likelihood = 0, fit_seconds = 0,
      constraint_residual = stats::setNames(rep(0, length(random_names)), random_names),
      observation_spatial_mean = stats::setNames(rep(0, length(random_names)), random_names),
      coefficients = list(global = numeric()), estimation = NULL
    )
  }
  task <- function() {
    function(payload, spec, base_offset, control, diagnostics, libpaths) {
      phase <- if (length(spec$random)) "full" else "null"
      calls$event <- c(calls$event, phase)
      calls$family[[length(calls$family) + 1L]] <- if (is.null(payload$poisson)) {
        rep(spec$family, nrow(payload$Y))
      } else ifelse(payload$poisson, "poisson", spec$family)
      lapply(seq_len(nrow(payload$Y)), function(j) {
        fake_fit(spec, payload$Y[j, ], phase)
      })
    }
  }
  marginal <- function(feature_id, score_sparse, nuisance_design, null_state,
                       features, ...) {
    calls$event <- c(calls$event, "score")
    calls$constraint <- score_sparse$constraint
    calls$W <- null_state$working_variance[, features, drop = FALSE]
    calls$statistic <- colSums(1 / calls$W)
    data.frame(
      feature_id = feature_id[features], statistic = calls$statistic,
      p_value = c(.2, .3)[seq_along(features)], method_requested = "liu",
      method_used = "liu", fallback_used = FALSE, fallback_reason = NA_character_,
      davies_ifault = NA_integer_, error_message = NA_character_
    )
  }
  route <- function(Y, X, offset, threshold, active) {
    list(phi = c(1, 2), poisson = c(TRUE, FALSE))
  }
  fit <- testthat::with_mocked_bindings(
    inlaST.estimate(Y, model, BPPARAM = BiocParallel::SerialParam(),
                    retain_marginal = TRUE),
    .inlast_chunk_task = task, .inlast_null_marginal = marginal,
    .mgcvst_prescreen_route = route, .package = "mgcvST"
  )

  expect_identical(calls$event, c("null", "score", "full"))
  expect_identical(calls$family[[1L]], c("poisson", "negative_binomial"))
  expect_identical(calls$family[[2L]], c("poisson", "negative_binomial"))
  expect_identical(calls$constraint, model$inla_spec$random[[1L]]$constraint)
  expect_identical(calls$constraint, fit$score_sparse$constraint)
  expect_true(all(fit$diagnostics$null_converged))
  expect_true(all(fit$diagnostics$converged))
  expect_null(fit$null_score)
  expect_identical(names(fit$marginal_data$null_state),
                   c("working_error", "working_variance", "nuisance_precision"))
  expect_equal(unname(calls$W), matrix(c(rep(3, nrow(f$data)),
    rep(4, nrow(f$data))), nrow(f$data), 2L), tolerance = 1e-12)
  expect_false(identical(calls$W[, 1L], calls$W[, 2L]))
  expect_false(identical(calls$statistic[1L], calls$statistic[2L]))
  expect_false(identical(calls$W, fit$working_variance))
  replay <- mgcvST.marginal(fit, calibration = "liu",
                             BPPARAM = BiocParallel::SerialParam())
  expect_equal(replay$p_value, fit$diagnostics$marginal_p_value, tolerance = 1e-12)
  fit$marginal_data$result <- NULL
  recomputed <- testthat::with_mocked_bindings(
    mgcvST.marginal(fit, calibration = "liu", BPPARAM = BiocParallel::SerialParam()),
    .inlast_null_marginal = marginal, .package = "mgcvST"
  )
  expect_equal(recomputed$p_value, replay$p_value, tolerance = 1e-12)
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

test_that("the direct null score receives precomputed iid nuisance precision", {
  A <- Matrix::Matrix(c(1, 0, 0, 1), 2, 2, sparse = TRUE)
  spec <- list(
    fixed = list(X = matrix(1, 2, 1), names = "(Intercept)"),
    random = list(list(name = "batch", A = A, Q = Matrix::Diagonal(2),
      target = FALSE, kind = "nuisance", subtype = "iid", sp_index = 2L))
  )
  precision <- mgcvST:::.inlast_null_nuisance_precision(
    spec, cbind(NA_real_, c(6, 12)), c(2, 3), 1:2
  )
  expect_equal(precision, rbind(c(0, 0), c(3, 4), c(3, 4)), tolerance = 1e-12)
  seen <- new.env(parent = emptyenv())
  state <- list(
    working_error = matrix(0, 2, 2), working_variance = matrix(1, 2, 2),
    nuisance_precision = precision
  )
  direct <- function(score_sparse, nuisance_design, null_state, features,
                     threads = 1L) {
    seen$precision <- null_state$nuisance_precision[, features, drop = FALSE]
    lapply(features, function(j) list(error = "score fixture"))
  }
  testthat::with_mocked_bindings(
    mgcvST:::.inlast_null_marginal(
      c("g1", "g2"), list(), matrix(0, 2, 3), state, 1:2
    ),
    .inlast_sparse_null_batch = direct, .package = "mgcvST"
  )
  expect_identical(seen$precision, precision)
})

test_that("a null-fit error is reported separately and does not skip the full fit", {
  skip_if_not_installed("INLA")
  skip_if_not_installed("geometry")
  f <- .inlast_fixture(n = 20L, seed = 1822L)
  model <- inlaST.set(response ~ z + offset(offset0), f$data, f$basis,
                      family = poisson())
  Y <- rbind(g1 = rep(1, nrow(f$data)), g2 = rep(2, nrow(f$data)))
  calls <- new.env(parent = emptyenv())
  calls$event <- character()
  fake_fit <- function(spec, y) {
    random_names <- vapply(spec$random, `[[`, character(1L), "name")
    list(
      working_error = rep(1, length(y)), working_variance = rep(2, length(y)),
      dispersion = 1, family_parameters = numeric(), smoothing_parameters = 1,
      nuisance_covariance = diag(ncol(spec$nuisance_design)),
      expected_nuisance_covariance = NULL, converged = TRUE,
      log_marginal_likelihood = 0, fit_seconds = 0,
      constraint_residual = stats::setNames(rep(0, length(random_names)), random_names),
      observation_spatial_mean = stats::setNames(rep(0, length(random_names)), random_names),
      coefficients = list(global = numeric()), estimation = NULL
    )
  }
  task <- function() {
    function(payload, spec, base_offset, control, diagnostics, libpaths) {
      phase <- if (length(spec$random)) "full" else "null"
      calls$event <- c(calls$event, phase)
      lapply(seq_len(nrow(payload$Y)), function(j) {
        if (identical(phase, "null") && identical(j, 1L)) return(simpleError("null failed"))
        fake_fit(spec, payload$Y[j, ])
      })
    }
  }
  marginal <- function(feature_id, score_sparse, nuisance_design, null_state,
                       features, ...) {
    calls$event <- c(calls$event, "score")
    expect_identical(features, 2L)
    data.frame(
      feature_id = feature_id[features], statistic = 1, p_value = .5,
      method_requested = "liu", method_used = "liu", fallback_used = FALSE,
      fallback_reason = NA_character_, davies_ifault = NA_integer_,
      error_message = NA_character_
    )
  }
  fit <- testthat::with_mocked_bindings(
    inlaST.estimate(Y, model, BPPARAM = BiocParallel::SerialParam()),
    .inlast_chunk_task = task, .inlast_null_marginal = marginal, .package = "mgcvST"
  )

  expect_identical(calls$event, c("null", "score", "full"))
  expect_false(fit$diagnostics$null_converged[1L])
  expect_match(fit$diagnostics$null_error_message[1L], "null failed")
  expect_true(all(fit$diagnostics$converged))
  expect_true(is.na(fit$diagnostics$marginal_p_value[1L]))
  expect_equal(fit$diagnostics$marginal_p_value[2L], .5, tolerance = 1e-12)
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
