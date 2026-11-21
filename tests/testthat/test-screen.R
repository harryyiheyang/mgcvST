screen_fixture <- function(n = 30L, q = 16L, seed = 5L) {
  set.seed(seed)
  H <- lapply(seq_len(n), function(g) {
    Q <- qr.Q(qr(matrix(rnorm(q * q), q)))
    d <- (1:q)^(-runif(1, 0.3, 2.5)) * 10^runif(1, -2, 2)
    M <- Q %*% (d * t(Q))
    (M + t(M)) / 2
  })
  a <- sapply(H, function(M) {
    e <- eigen(M, symmetric = TRUE)
    drop(e$vectors %*% (sqrt(pmax(e$values, 0)) * rnorm(q)))
  })
  a[, 2L] <- a[, 2L] + 5 * a[, 1L]
  a[, 7L] <- a[, 7L] - 4 * a[, 3L]
  list(H = H, a = a, G = spdeST:::mgcvst_pair_basis_cpp(H, diag(q), 1L),
       pr = t(combn(n, 2L)))
}

test_that("screened pairs have exact p-values above the Gaussian bound", {
  f <- screen_fixture()
  z0 <- 2.6
  full <- spdeST:::mgcvst_pair_spa_cpp(f$H, f$G, f$a, f$pr[, 1L], f$pr[, 2L], 1L, 4L)
  scr <- spdeST:::mgcvst_pair_spa_cpp(f$H, f$G, f$a, f$pr[, 1L], f$pr[, 2L], 1L, 4L,
                                      NULL, z0)
  skipped <- scr$status == 4L
  expect_true(any(skipped) && any(!skipped))
  expect_equal(scr$screened, sum(skipped))
  expect_identical(scr$log_p_two_sided[!skipped], full$log_p_two_sided[!skipped])
  expect_identical(scr$score, full$score)
  # The bound: reported <= exact, and both above the screening tail.
  expect_true(all(scr$log_p_two_sided[skipped] <= full$log_p_two_sided[skipped] + 1e-12))
  expect_true(all(scr$log_p_two_sided[skipped] > log(2 * pnorm(-z0))))
  expect_true(all(full$log_p_two_sided[skipped] > log(2 * pnorm(-z0))))
  expect_error(spdeST:::mgcvst_pair_spa_cpp(f$H, f$G, f$a, f$pr[, 1L], f$pr[, 2L],
                                            1L, 4L, NULL, 1.5), "screen_z")
})

test_that("the Gaussian tail bounds the exact tail from z = 2 on", {
  skip_if_not_installed("CompQuadForm")
  set.seed(9)
  spectra <- list(1, c(1, 1), c(1, 0.3), rep(1, 5), rep(1, 60), (1:80)^-1, (1:80)^-2,
                  (1:80)^-0.5, c(1, rep(0.05, 99)), c(1, rep(0.2, 99)), 0.8^(0:40),
                  rexp(50)^2)
  for (s in spectra) {
    s <- s / sqrt(sum(s^2))
    for (z in c(2, 2.3, 2.7, 3.2, 4)) {
      exact <- 2 * CompQuadForm::imhof(z, c(s / 2, -s / 2), epsabs = 1e-12,
                                       epsrel = 1e-12, limit = 1e5)$Qq
      expect_gte(exact, 2 * pnorm(-z))
    }
  }
})

test_that("the screening value follows the adjustment", {
  by <- spdeST:::.mgcvst_screen_z("BY", 0.05, 3736011)
  expect_equal(by$theta, 0.05 / sum(1 / seq_len(3736011)), tolerance = 1e-10)
  expect_equal(by$z, qnorm(by$theta, lower.tail = FALSE))
  expect_gte(by$z, 2)
  expect_identical(spdeST:::.mgcvst_screen_z("BH", 0.05, 1000)$z, 0)
  expect_gt(spdeST:::.mgcvst_screen_z("BH", 0.01, 1000)$z, 2)
  expect_gt(spdeST:::.mgcvst_screen_z("Sidak", 0.05, 1000)$z, 2)
})

test_that("discoveries and their q-values do not change with the screen", {
  f <- screen_fixture(n = 40L)
  m <- nrow(f$pr)
  full <- spdeST:::mgcvst_pair_spa_cpp(f$H, f$G, f$a, f$pr[, 1L], f$pr[, 2L], 1L, 4L)
  for (adjust in c("BY", "Sidak")) {
    z0 <- spdeST:::.mgcvst_screen_z(adjust, 0.05, m)$z
    expect_gte(z0, 2)
    scr <- spdeST:::mgcvst_pair_spa_cpp(f$H, f$G, f$a, f$pr[, 1L], f$pr[, 2L], 1L,
                                        4L, NULL, z0)
    qf <- spdeST:::.mgcvst_log_adjust(full$log_p_two_sided, adjust)$log_q
    qs <- spdeST:::.mgcvst_log_adjust(scr$log_p_two_sided, adjust)$log_q
    found <- qf <= log(0.05)
    expect_true(any(found))
    expect_identical(qs <= log(0.05), found)
    expect_identical(qs[found], qf[found])
    expect_true(all(scr$status[found] == 0L))
  }
})

test_that("inlaST.test(screen = TRUE) keeps the discoveries of both routes", {
  skip_if_not_installed("geometry")
  set.seed(31)
  n <- 900L
  p <- 24L
  xy <- as.matrix(expand.grid(x = seq(0, 1, length.out = 5), y = seq(0, 1, length.out = 5)))
  mesh <- list(loc = xy, graph = list(tv = geometry::delaunayn(xy)))
  data <- data.frame(x = runif(n, .01, .99), y = runif(n, .01, .99), z = runif(n))
  basis <- spde_basis(mesh, as.matrix(data[, c("x", "y")]), kappa = .7)
  field <- cbind(sin(6 * data$x) + cos(5 * data$y), cos(7 * data$x * data$y))
  load <- matrix(0, p, 2L)
  load[1:5, 1L] <- 1.2
  load[6:9, 2L] <- 1.2
  load[10:11, 1L] <- -1.2
  Y <- t(sapply(seq_len(p), function(g) rpois(n, exp(2 + drop(field %*% load[g, ])))))
  rownames(Y) <- paste0("g", seq_len(p))
  data$response <- Y[1L, ]
  model <- inlaST.set(response ~ z, data, basis, family = stats::poisson())
  sp <- BiocParallel::SerialParam()
  fit <- spdeST:::inlaST.estimate_spatial(inlaST.estimate(Y, model, BPPARAM = sp),
                                          Y, features = "all", BPPARAM = sp)
  full <- inlaST.test(fit, moments = "exact")
  scr <- inlaST.test(fit, moments = "exact", screen = TRUE)
  expect_gte(scr$timing$screen$z, 2)
  skipped <- scr$results$status == 4L
  expect_true(any(skipped))
  expect_identical(scr$timing$screen$screened, as.numeric(sum(skipped)))
  expect_identical(scr$discoveries, full$discoveries)
  found <- full$results$log_q <= log(0.05)
  expect_true(any(found))
  expect_identical(scr$results$log_q <= log(0.05), found)
  expect_identical(scr$results$log_q[found], full$results$log_q[found])
  expect_identical(scr$results[!skipped, c("i", "j", "score", "log_p_two_sided")],
                   full$results[!skipped, c("i", "j", "score", "log_p_two_sided")])
  expect_true(all(scr$results$log_p_two_sided[skipped] <=
                    full$results$log_p_two_sided[skipped] + 1e-12))
  # The PCAlearning route: the same discoveries, bounds below the route's own
  # p-values, and evaluated pairs unchanged.
  full <- inlaST.test(fit, moments = "pcalearning", rank = 5L, k = 10L)
  scr <- inlaST.test(fit, moments = "pcalearning", rank = 5L, k = 10L, screen = TRUE)
  skipped <- scr$results$status == 4L
  expect_true(any(skipped))
  expect_identical(scr$timing$screen$screened, as.numeric(sum(skipped)))
  expect_identical(scr$discoveries, full$discoveries)
  found <- full$results$log_q <= log(0.05)
  expect_true(any(found))
  expect_identical(scr$results$log_q <= log(0.05), found)
  expect_identical(scr$results$log_p_two_sided[!skipped],
                   full$results$log_p_two_sided[!skipped])
  expect_true(all(scr$results$log_p_two_sided[skipped] <=
                    full$results$log_p_two_sided[skipped]))
  expect_true(all(scr$results$log_p_two_sided[skipped] > log(scr$timing$screen$tau)))
})
