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

test_that("single-global sparse score matches dense with native nuisance Vp", {
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

  dense <- fit
  dense$score_backend <- "dense"
  dense$score_sparse <- NULL
  dense$.mgcvst_fixed_factors <- NULL
  for (j in seq_len(nrow(f$Y))) {
    sparse_state <- mgcvST:::.mgcvst_model_score_state(fit, j)
    dense_state <- mgcvST:::.mgcvst_model_score_state(dense, j)
    expect_equal(sparse_state$a, dense_state$a, tolerance = 1e-9)
    expect_equal(sparse_state$M, dense_state$M, tolerance = 1e-9)
    expect_identical(sparse_state$width, dense_state$width)
  }
  pairs <- t(combn(rownames(f$Y), 2L))
  sparse_test <- mgcvST.test(
    fit, pairs = pairs, calibration = "liu",
    BPPARAM = BiocParallel::SerialParam()
  )
  dense_test <- mgcvST.test(
    dense, pairs = pairs, calibration = "liu",
    BPPARAM = BiocParallel::SerialParam()
  )
  numeric <- c("signed_score", "information", "p_two_sided",
               "p_positive", "p_negative")
  for (field in numeric) {
    expect_equal(
      sparse_test$results[[field]], dense_test$results[[field]],
      tolerance = 1e-9
    )
  }

  # Preserve the native posterior block even when it is not exact expected GLS.
  perturbed_sparse <- fit
  perturbed_dense <- dense
  perturbed_sparse$nuisance_covariance <-
    lapply(fit$nuisance_covariance, `*`, 0.997)
  perturbed_dense$nuisance_covariance <- perturbed_sparse$nuisance_covariance
  expect_gt(max(abs(
    perturbed_sparse$nuisance_covariance[[1L]] -
      fit$expected_nuisance_covariance[[1L]]
  )), 0)
  for (j in seq_len(nrow(f$Y))) {
    expect_equal(
      mgcvST:::.mgcvst_model_score_state(perturbed_sparse, j)$a,
      mgcvST:::.mgcvst_model_score_state(perturbed_dense, j)$a,
      tolerance = 1e-9
    )
    expect_equal(
      mgcvST:::.mgcvst_model_score_state(perturbed_sparse, j)$M,
      mgcvST:::.mgcvst_model_score_state(perturbed_dense, j)$M,
      tolerance = 1e-9
    )
  }

  broken <- fit
  broken$score_sparse$constraint[] <- 0
  expect_error(
    mgcvST:::.mgcvst_model_score_state(broken, 1L),
    "constraint|positive"
  )
  expect_identical(broken$score_backend, "sparse")
})

test_that("sparse score workers agree over SOCK", {
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
  bp <- BiocParallel::SnowParam(2L, type = "SOCK", progressbar = FALSE)
  on.exit(BiocParallel::bpstop(bp), add = TRUE)
  socket <- mgcvST.test(
    fit, pairs = pairs, calibration = "liu", BPPARAM = bp,
    chunk_size = 1L
  )
  expect_equal(socket$results$signed_score, serial$results$signed_score,
               tolerance = 1e-10)
  expect_equal(socket$results$p_two_sided, serial$results$p_two_sided,
               tolerance = 1e-10)
  expect_true(all(is.na(socket$results$error_message)))
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
    dense_state <- mgcvST:::.mgcvst_model_score_state(dense, 1L)
    expect_equal(sparse_state$a, dense_state$a, tolerance = 1e-9)
    expect_equal(sparse_state$M, dense_state$M, tolerance = 1e-9)
  }
})
