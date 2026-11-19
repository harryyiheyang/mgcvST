spa <- function(q, lambda) spdeST:::mgcvst_spa_quadform_cpp(q, lambda)

test_that("the Stage 1 saddlepoint agrees with an accurate Davies reference", {
  skip_if_not_installed("CompQuadForm")
  set.seed(11)
  lambda <- sort(rexp(300)^2, decreasing = TRUE)
  mu <- sum(lambda); sd <- sqrt(2 * sum(lambda^2))
  for (q in mu + sd * c(-1, 0.5, 3, 6, 10)) {
    d <- CompQuadForm::davies(q, lambda, lim = 1e6, acc = 1e-14)
    expect_identical(d$ifault, 0L)
    expect_lt(abs(log10(spa(q, lambda)) - log10(d$Qq)), 0.05)
  }
})

test_that("the Stage 1 saddlepoint stays finite beyond machine precision", {
  lambda <- c(5, 2, rep(0.5, 50))
  q <- 1200
  sp <- spa(q, lambda)
  expect_true(is.finite(sp) && sp > 0 && sp < 1e-40)
  # exponential tail rate 1/(2 * max(lambda)) (Chen and Lumley 2019, Theorem 1)
  sp2 <- spa(q + 10, lambda)
  expect_lt(abs((log(sp) - log(sp2)) / 10 - 1 / (2 * max(lambda))), 0.01)
})

test_that("the Stage 1 saddlepoint uses its limit at the mean and returns 1 for q <= 0", {
  lambda <- c(5, 2, rep(0.5, 50))
  mu <- sum(lambda)
  rho3 <- 8 * sum(lambda^3) / (2 * sum(lambda^2))^1.5
  centre <- spa(mu, lambda)
  expect_equal(centre, 0.5 - rho3 / (6 * sqrt(2 * pi)), tolerance = 1e-12)
  for (f in c(0.999, 1.001)) expect_lt(abs(spa(mu * f, lambda) - centre), 0.01)
  expect_identical(spa(0, lambda), 1)
  expect_identical(spa(-2, lambda), 1)
  expect_true(is.na(spa(3, numeric(0))))
})

test_that("spde.scan takes factors as 0-1 columns of the null design", {
  skip_if_not_installed("geometry")
  set.seed(17)
  xy <- as.matrix(expand.grid(x = seq(0, 1, length.out = 5), y = seq(0, 1, length.out = 5)))
  mesh <- list(loc = xy, graph = list(tv = geometry::delaunayn(xy)))
  n <- 150L
  data <- data.frame(x = runif(n, .01, .99), y = runif(n, .01, .99),
                     slide = factor(rep(c("a", "b", "c"), length.out = n)),
                     z = rnorm(n), offset0 = 0)
  basis <- spde_basis(mesh, as.matrix(data[, c("x", "y")]), kappa = .7)
  eta <- 1 + c(a = 0, b = .8, c = -.5)[as.character(data$slide)] + .3 * data$z
  data$response <- rnbinom(n, mu = exp(eta), size = 3)
  model <- mgcvST.set(response ~ slide + z + s(x, y, bs = "spde", xt = basis),
                      data = data, family = mgcv::nb())
  expect_identical(ncol(model$null_X), 4L)
  Y <- rbind(g1 = data$response, g2 = rnbinom(n, mu = exp(eta), size = 3))
  withr::local_options(mgcvST.progress = FALSE)
  scan <- spdeST.scan(Y, model, BPPARAM = BiocParallel::SerialParam())
  expect_true(all(is.finite(scan$diagnostics$marginal_p_value)))
  field <- spdeST:::.spde_scan_field(model$G$smooth[[1L]]$basis)
  direct <- spde.scan(Y, model$null_X, field$A, field$Q)
  expect_equal(scan$diagnostics$marginal_p_value, direct$p_value, tolerance = 1e-10)
  expect_error(mgcvST.estimate(Y, mgcvST.set(
                 response ~ s(z, k = 5) + s(x, y, bs = "spde", xt = basis),
                 data = data, family = mgcv::nb()),
                               BPPARAM = BiocParallel::SerialParam()),
               "without smooth terms")
})
