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

test_that("scan records carry the working response of the INLA null fit", {
  data <- scan_data(9048L, features = 5L)
  spec <- list(fixed = list(X = data$X, names = colnames(data$X)),
               random = list(), family = "negative_binomial")
  expect_true(mgcvST:::.inlast_scan_eligible(spec))
  withr::local_options(mgcvST.scan_null = FALSE)
  expect_false(mgcvST:::.inlast_scan_eligible(spec))
  withr::local_options(mgcvST.scan_null = TRUE)
  with_random <- spec
  with_random$random <- list(list(A = diag(2)))
  expect_false(mgcvST:::.inlast_scan_eligible(with_random))
  extra <- matrix(stats::rnorm(5 * ncol(data$Y), 0, .1), 5)
  fits <- mgcvST:::.inlast_scan_fits(
    seq_len(5L), Matrix::Matrix(data$Y, sparse = TRUE), spec, data$offset,
    extra, control = list(), poisson = c(TRUE, FALSE, FALSE, FALSE, FALSE))
  expect_length(fits, 5L)
  for (j in seq_len(5L)) {
    f <- fits[[j]]
    o <- data$offset + extra[j, ]
    expect_equal(f$eta, as.numeric(data$X %*% f$fixed_mode) + o)
    expect_equal(f$working_error,
                 f$eta + (data$Y[j, ] - f$mu) / f$mu - o)
    if (j == 1L) {
      expect_identical(f$family, "poisson")
      expect_equal(f$working_variance, 1 / f$mu)
    } else {
      expect_identical(f$family, "negative_binomial")
      expect_equal(f$working_variance, 1 / f$mu + 1 / f$family_parameters)
    }
    expect_true(f$converged)
    expect_length(f$random_mode, 0L)
  }
  route <- mgcvST:::.mgcvst_prescreen_route(data$Y, data$X, data$offset, TRUE)
  expect_identical(route$poisson, route$score <= 0)
  fits0 <- mgcvST:::.inlast_scan_fits(
    seq_len(5L), Matrix::Matrix(data$Y, sparse = TRUE), spec, data$offset,
    NULL, control = list())
  expect_identical(vapply(fits0, `[[`, "", "family") == "poisson",
                   route$poisson)
  fixed <- mgcvST:::.inlast_scan_fits(
    2L, Matrix::Matrix(data$Y, sparse = TRUE), spec, data$offset, NULL,
    control = list(nb_size = 3))
  expect_identical(fixed[[1L]]$family_parameters, 3)
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
