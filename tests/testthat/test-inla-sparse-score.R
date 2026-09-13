.inlast_sparse_fixture <- function(n = 72L, seed = 1701L) {
  skip_if_not_installed("INLA")
  skip_if_not_installed("geometry")
  set.seed(seed)
  vertices <- as.matrix(expand.grid(
    x = seq(0, 1, length.out = 5L),
    y = seq(0, 1, length.out = 5L)
  ))
  mesh <- list(loc = vertices, graph = list(tv = geometry::delaunayn(vertices)))
  data <- data.frame(
    x = runif(n, 0.02, 0.98), y = runif(n, 0.02, 0.98),
    z = seq(-1, 1, length.out = n), exposure = runif(n, 0.8, 1.3)
  )
  data$offset0 <- log(data$exposure)
  basis <- spde_basis(
    mesh, as.matrix(data[c("x", "y")]), kappa = 1.2,
    project_intercept = TRUE
  )
  eta <- 0.4 + 0.25 * data$z + data$offset0 + 0.2 * sin(2 * pi * data$x)
  Y <- rbind(
    a = eta + rnorm(n, sd = 0.3),
    b = eta + 0.1 * cos(2 * pi * data$y) + rnorm(n, sd = 0.3),
    c = eta - 0.1 * sin(2 * pi * data$y) + rnorm(n, sd = 0.3)
  )
  list(data = data, basis = basis, Y = Y)
}

test_that("single-global sparse score preserves expected-curvature Gram and traces", {
  skip_on_cran()
  f <- .inlast_sparse_fixture()
  model <- inlaST.set(
    response ~ z + offset(offset0), f$data, f$basis, family = gaussian()
  )
  fit <- inlaST.estimate(
    f$Y, model, score_backend = "auto", diagnostics = TRUE,
    BPPARAM = BiocParallel::SerialParam(),
    control = list(fixed_precision = 1.7, gaussian_precision = 1 / 0.09)
  )
  expect_identical(fit$score_backend_requested, "auto")
  expect_identical(fit$score_backend, "sparse")
  expect_true(is.list(fit$score_sparse))
  expect_true(all(vapply(
    mgcvST:::.mgcvst_model_fixed_factors(fit), is.null, logical(1L)
  )))

  states <- mgcvST:::.inlast_sparse_batch(
    fit, seq_len(nrow(f$Y)), threads = 1L, score_only = FALSE
  )
  for (j in seq_len(nrow(f$Y))) {
    singleton <- mgcvST:::.mgcvst_model_sparse_score_state(fit, j)
    score_only <- mgcvST:::.mgcvst_model_sparse_score_state(
      fit, j, score_only = TRUE
    )
    expect_equal(states[[j]]$a, singleton$a, tolerance = 1e-10)
    expect_equal(states[[j]]$M, singleton$M, tolerance = 1e-10)
    expect_equal(score_only$a, singleton$a, tolerance = 1e-10)
    expect_null(score_only$M)
  }
  pairs <- t(combn(rownames(f$Y), 2L))
  sparse_test <- mgcvST.test(
    fit, pairs = pairs, calibration = "liu",
    BPPARAM = BiocParallel::SerialParam()
  )
  local <- matrix(match(pairs, rownames(f$Y)), ncol = 2L)
  M <- lapply(states, `[[`, "M")
  moments <- mgcvST:::mgcvst_pair_trace_powers_cpp(M, local, 4L, 1L)
  for (j in seq_len(nrow(local))) {
    score <- sum(states[[local[j, 1L]]]$a * states[[local[j, 2L]]]$a)
    liu <- mgcvST:::.liu_squared_score_moments(
      abs(score), moments[j, 1L], moments[j, 2L],
      moments[j, 3L], moments[j, 4L]
    )
    expect_equal(sparse_test$results$signed_score[j], score, tolerance = 1e-10)
    expect_equal(sparse_test$results$information[j], moments[j, 1L], tolerance = 1e-9)
    expect_equal(sparse_test$results$p_two_sided[j], liu$p_value, tolerance = 1e-10)
  }

  broken <- fit
  broken$score_sparse$constraint[] <- 0
  expect_error(
    mgcvST:::.mgcvst_model_score_state(broken, 1L),
    "constraint|positive"
  )
  expect_identical(broken$score_backend, "sparse")
})

test_that("sparse INLA downstream rejects SOCK and agrees across OpenMP counts", {
  skip_on_cran()
  f <- .inlast_sparse_fixture(n = 64L, seed = 1711L)
  model <- inlaST.set(
    response ~ z + offset(offset0), f$data, f$basis, family = gaussian()
  )
  fit <- inlaST.estimate(
    f$Y, model, score_backend = "sparse",
    BPPARAM = BiocParallel::SerialParam(),
    control = list(fixed_precision = 2, gaussian_precision = 1 / 0.09)
  )
  pairs <- t(combn(rownames(f$Y), 2L))
  serial <- mgcvST.test(
    fit, pairs = pairs, calibration = "liu",
    BPPARAM = BiocParallel::SerialParam(), chunk_size = 1L
  )
  threaded <- mgcvST.test(
    fit, pairs = pairs, calibration = "liu",
    BPPARAM = BiocParallel::SerialParam(), threads = 2L,
    chunk_size = 1L
  )
  expect_equal(threaded$results$signed_score, serial$results$signed_score,
               tolerance = 1e-10)
  expect_equal(threaded$results$p_two_sided, serial$results$p_two_sided,
               tolerance = 1e-10)
  bp <- BiocParallel::SnowParam(2L, type = "SOCK", progressbar = FALSE)
  expect_error(
    mgcvST.test(fit, pairs = pairs, calibration = "liu", BPPARAM = bp),
    "must be SerialParam"
  )
})

test_that("unsupported sparse structures error or use the dense auto backend", {
  skip_on_cran()
  f <- .inlast_sparse_fixture(n = 56L, seed = 1721L)
  s <- mgcv::s
  model <- inlaST.set(
    response ~ s(z, k = 5) + offset(offset0),
    f$data, f$basis, family = gaussian()
  )
  expect_error(
    inlaST.estimate(
      f$Y[1L, , drop = FALSE], model, score_backend = "sparse",
      BPPARAM = BiocParallel::SerialParam(),
      control = list(fixed_precision = c(2, 2), gaussian_precision = 1 / 0.09)
    ),
    "exactly one random block"
  )
  fit <- inlaST.estimate(
    f$Y[1L, , drop = FALSE], model, score_backend = "auto",
    BPPARAM = BiocParallel::SerialParam(),
    control = list(fixed_precision = c(2, 2), gaussian_precision = 1 / 0.09)
  )
  expect_identical(fit$score_backend_requested, "auto")
  expect_identical(fit$score_backend, "dense")
  expect_null(fit$score_sparse)
  expect_true(is.finite(
    mgcvST:::.mgcvst_model_score_state(fit, 1L)$M[1L, 1L]
  ))
})

test_that("Gaussian sparse units cover zero-X and observation-scale multi-X", {
  skip_on_cran()
  f <- .inlast_sparse_fixture(n = 52L, seed = 1731L)
  zero_signal <- 0.15 * sin(2 * pi * f$data$x) + rnorm(nrow(f$data), sd = 0.3)
  zero_signal <- zero_signal - mean(zero_signal)
  cases <- list(
    zero_X = list(
      formula = response ~ 0 + offset(offset0), scale = "raw",
      fixed_columns = 0L,
      response = matrix(f$data$offset0 + zero_signal, nrow = 1L)
    ),
    observation_multi_X = list(
      formula = response ~ z + x + offset(offset0), scale = "observation",
      fixed_columns = 3L, response = f$Y[1L, , drop = FALSE]
    )
  )
  for (case_name in names(cases)) {
    case <- cases[[case_name]]
    model <- inlaST.set(
      case$formula, f$data, f$basis, family = gaussian(),
      precision_scale = case$scale
    )
    if (identical(case_name, "zero_X")) {
      expect_identical(model$inla_spec$fixed$names, character())
      expect_identical(dim(model$inla_spec$fixed$X), c(nrow(f$data), 0L))
    }
    fit <- inlaST.estimate(
      case$response, model, score_backend = "auto",
      diagnostics = TRUE, BPPARAM = BiocParallel::SerialParam(),
      control = list(fixed_precision = 0.8, gaussian_precision = 5)
    )
    expect_true(fit$diagnostics$converged, info = case_name)
    expect_identical(ncol(fit$geometry$nuisance_design), case$fixed_columns)
    expect_identical(dim(fit$nuisance_covariance[[1L]]),
                     rep(case$fixed_columns, 2L))
    expect_identical(fit$score_backend, "sparse")
    expect_equal(fit$smoothing_parameters[[1L]], 0.8 / 5,
                 tolerance = 1e-12)
    diagnostic <- fit$inla_diagnostics[[1L]]
    expect_equal(diagnostic$tau[[1L]], 0.8, tolerance = 1e-12)
    expect_equal(
      diagnostic$tau_internal[[1L]] * diagnostic$precision_scale[[1L]],
      0.8, tolerance = 1e-12
    )
    expect_lt(abs(diagnostic$constraint_residual_uncorrected[[1L]]), 1e-6)
    expect_lt(abs(diagnostic$constraint_residual[[1L]]), 1e-14)
    expect_lt(abs(diagnostic$observation_spatial_mean[[1L]]), 1e-14)

    dense <- fit
    dense$score_backend <- "dense"
    dense$score_sparse <- dense$.mgcvst_fixed_factors <- NULL
    sparse_state <- mgcvST:::.mgcvst_model_score_state(fit, 1L)
    dense$nuisance_covariance[[1L]] <- sparse_state$expected_vp
    dense_state <- mgcvST:::.mgcvst_model_score_state(dense, 1L)
    expect_equal(sum(sparse_state$a^2), sum(dense_state$a^2), tolerance = 1e-9)
    expect_equal(sum(diag(sparse_state$M)), sum(diag(dense_state$M)), tolerance = 1e-9)
    expect_equal(sum(sparse_state$M^2), sum(dense_state$M^2), tolerance = 1e-9)
  }
})
