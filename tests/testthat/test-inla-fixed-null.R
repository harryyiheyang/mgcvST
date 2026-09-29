test_that("latent modes retain empty random blocks and validate mismatches", {
  tags <- paste0(".inlast_x", 1:6)
  fit <- list(mode = list(x = c(99, 98, seq_len(6L))),
              summary.fixed = data.frame(mean = seq_len(6L), row.names = tags))
  modes <- mgcvST:::.inlast_latent_modes(fit, tags, character(), integer())
  expect_identical(modes$fixed, stats::setNames(as.numeric(1:6), tags))
  expect_length(modes$random, 0L)
  expect_error(mgcvST:::.inlast_latent_modes(
    fit, tags, ".inlast_r", integer()), "block specification is invalid")
})

test_that("fixed-only NB modes preserve covariates and offsets", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  set.seed(9021L)
  d <- data.frame(group = factor(rep(letters[1:6], each = 20L)),
                  offset0 = seq(-.3, .3, length.out = 120L))
  y <- stats::rnbinom(nrow(d), mu = exp(2 + d$offset0 +
    c(0, -.2, .1, .3, -.1, .2)[d$group]), size = 2)
  for (formula in list(~ 1, ~ group)) {
    X <- stats::model.matrix(formula, d)
    spec <- list(fixed = list(X = X, names = colnames(X)),
                  random = list(), family = "negative_binomial")
    fit <- mgcvST:::.inlast_fit_feature(spec, y, offset = d$offset0,
      control = list(nb_size = 2, num_threads = 1L))
    objective <- function(beta) {
      -sum(stats::dnbinom(y, mu = exp(d$offset0 + X %*% beta),
                          size = 2, log = TRUE))
    }
    gradient <- function(beta) {
      mu <- exp(d$offset0 + X %*% beta)
      -as.numeric(crossprod(X, 2 * (y - mu) / (2 + mu)))
    }
    reference <- stats::optim(rep(0, ncol(X)), objective, gradient,
      method = "BFGS", control = list(reltol = 1e-12))
    expect_identical(reference$convergence, 0L)
    expect_true(fit$converged)
    expect_length(fit$random_mode, 0L)
    expect_length(fit$tau, 0L)
    expect_length(fit$coefficients, 0L)
    expect_identical(names(fit$fixed_mode), colnames(X))
    expect_equal(unname(fit$fixed_mode), reference$par, tolerance = 2e-4)
    expect_equal(fit$eta, as.numeric(d$offset0 + X %*% fit$fixed_mode),
      tolerance = 1e-12)
    expect_equal(fit$working_variance, 1 / fit$mu + 1 / 2,
      tolerance = 1e-12)
    expect_equal(fit$working_error,
      fit$eta + (y - fit$mu) / fit$mu - d$offset0, tolerance = 1e-12)
  }
})

test_that("public NB estimates retain fixed-only null scores", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  skip_if_not_installed("geometry")
  set.seed(9022L)
  vertices <- as.matrix(expand.grid(x = seq(0, 1, length.out = 4L),
                                   y = seq(0, 1, length.out = 4L)))
  mesh <- list(loc = vertices, graph = list(tv = geometry::delaunayn(vertices)))
  d <- data.frame(x = runif(120L, .02, .98), y = runif(120L, .02, .98),
    group = factor(rep(letters[1:6], each = 20L)),
    offset0 = seq(-.3, .3, length.out = 120L))
  y <- stats::rnbinom(nrow(d), mu = exp(2 + d$offset0 +
    c(0, -.2, .1, .3, -.1, .2)[d$group]), size = 2)
  extra_offset <- seq(.1, -.1, length.out = nrow(d))
  for (formula in list(response ~ offset(offset0),
                       response ~ group + offset(offset0))) {
    model <- inlaST.set(formula, d, family = mgcv::nb(), mesh = mesh,
      kappa = .7, coordinates = c("x", "y"))
    fit <- inlaST.estimate(matrix(y, nrow = 1L), model,
      offset = extra_offset, BPPARAM = BiocParallel::SerialParam(),
      control = list(fixed_precision = 2, num_threads = 1L), threads = 1L)
    expect_true(fit$diagnostics$converged)
    expect_true(fit$diagnostics$null_converged)
    expect_true(is.na(fit$diagnostics$null_error_message))
    expect_identical(fit$diagnostics$family_used, "negative_binomial")
    expect_true(is.finite(fit$diagnostics$marginal_p_value))
    expect_equal(fit$offset, model$offset + extra_offset, tolerance = 0)
  }
})
