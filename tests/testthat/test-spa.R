test_that("the saddlepoint kernel matches an R reference on a full spectrum", {
  set.seed(1)
  s <- sort(abs(rnorm(40L)) * exp(-seq(0, 3, length.out = 40L)), decreasing = TRUE)
  Tm <- .spa_powers(s)
  x <- c(0.05, 0.3, 1, 2, 4, 8, 15, 30) * sqrt(sum(s^2))
  z <- .spa_kernel(x, s, Tm)
  expect_identical(unname(z[, "status"]), rep(0, length(x)))
  expect_identical(unname(z[, "remainder_kind"]), rep(0, length(x)))
  ref <- vapply(x, function(v) .spa_ref(v, s), numeric(1L))
  expect_equal(unname(z[, "log_p_two_sided"]), ref, tolerance = 1e-8)
  # The two one-sided tails partition the probability.
  expect_equal(exp(z[, "log_p_positive"]) + exp(z[, "log_p_negative"]),
               rep(1, length(x)), tolerance = 1e-12)
  expect_true(all(z[, "log_p_positive"] <= -log(2) + 1e-15))
  expect_equal(unname(z[, "log_p_two_sided"]),
               unname(pmin(0, z[, "log_p_positive"] + log(2))), tolerance = 1e-12)
})

test_that("U = 0 gives p = 1 and U < 0 mirrors U > 0", {
  s <- c(1, 0.6, 0.3)
  Tm <- .spa_powers(s)
  zero <- .spa_kernel(0, s, Tm)
  expect_identical(unname(zero[, "status"]), 0)
  expect_identical(unname(zero[, "log_p_two_sided"]), 0)
  expect_equal(unname(zero[, c("log_p_positive", "log_p_negative")]), rep(-log(2), 2L))
  x <- c(0.4, 1.5, 3)
  pos <- .spa_kernel(x, s, Tm)
  neg <- .spa_kernel(-x, s, Tm)
  expect_identical(neg[, "log_p_two_sided"], pos[, "log_p_two_sided"])
  expect_identical(neg[, "log_p_positive"], pos[, "log_p_negative"])
  expect_identical(neg[, "log_p_negative"], pos[, "log_p_positive"])
  # |w| below 1e-3 is the middle of the distribution: two-sided p = 1.
  tiny <- .spa_kernel(1e-6, s, Tm)
  expect_identical(unname(tiny[, "log_p_two_sided"]), 0)
  expect_equal(unname(tiny[, "log_p_positive"]), -log(2))
})

test_that("a tail far beyond the double range is finite in log space", {
  s <- c(1, 0.6, 0.3)
  Tm <- .spa_powers(s)
  x <- c(40, 200, 1500) * sqrt(sum(s^2))
  z <- .spa_kernel(x, s, Tm)
  expect_true(all(z[, "status"] == 0))
  expect_true(all(is.finite(z[, "log_p_two_sided"])))
  expect_lt(z[3L, "log_p_two_sided"], -1000)
  expect_true(all(diff(z[, "log_p_two_sided"]) < 0))
  ref <- vapply(x, function(v) .spa_ref(v, s), numeric(1L))
  expect_equal(unname(z[, "log_p_two_sided"]), ref, tolerance = 1e-8)
})

test_that("each remainder form is used when its moments allow it", {
  s <- c(1, 0.8, 0.5)
  lead <- .spa_powers(s)
  x <- c(0.5, 2, 4)
  check <- function(mu, order, kind) {
    z <- .spa_kernel(x, s, lead + mu, order)
    expect_identical(unique(unname(z[, "remainder_kind"])), kind)
    expect_true(all(z[, "status"] == 0))
    rem <- if (order == 4L) .spa_rem_two(mu, max(s)) else .spa_rem_one(mu, max(s))
    expect_identical(rem$kind, as.integer(kind))
    ref <- vapply(x, function(v) .spa_ref(v, s, rem$m, rem$u, rem$g), numeric(1L))
    expect_equal(unname(z[, "log_p_two_sided"]), ref, tolerance = 1e-8)
    z
  }
  nodes <- list(u = c(0.09, 0.01), m = c(4, 10))
  mu <- vapply(1:4, function(r) sum(nodes$m * nodes$u^r), numeric(1L))
  # Two nodes: four moments, order 4.
  two <- check(mu, 4L, 2)
  # One node: the same moments at order 2 (and at order 4 when mu_4 <= 0).
  one <- check(mu, 2L, 1)
  expect_gt(max(abs(two[, "log_p_two_sided"] - one[, "log_p_two_sided"])), 1e-6)
  check(c(mu[1:3], -1e-4), 4L, 1)
  # Gaussian term: mu_2 = 0 leaves only the variance mu_1.
  check(c(0.3, 0, 0, 0), 4L, 3)
  check(c(0.3, 0, 0, 0), 2L, 3)
  # No remainder: the total sums equal the leading sums, or differ by less
  # than 1e-10 of the total.
  for (eps in c(0, 1e-12)) {
    z <- .spa_kernel(x, s, lead * (1 + eps))
    expect_identical(unique(unname(z[, "remainder_kind"])), 0)
    expect_equal(unname(z[, "log_p_two_sided"]),
                 vapply(x, function(v) .spa_ref(v, s), numeric(1L)), tolerance = 1e-8)
  }
  # Moments above 1e-10 of the total are kept.
  kept <- .spa_kernel(x, s, lead * (1 + 1e-7))
  expect_true(all(kept[, "remainder_kind"] > 0))
})

test_that("only one-node remainders above the leading values use the Gaussian fallback", {
  s <- c(1, 0.5)
  lead <- .spa_powers(s)
  x <- c(5, 20)
  node <- function(u, m) vapply(1:4, function(r) sum(m * u^r), numeric(1L))
  # One node m = 0.5, u = 1.5 > s_1^2 = 1.
  z <- .spa_kernel(x, s, lead + node(1.5, 0.5), 2L)
  expect_identical(unname(z[, "status"]), c(0, 0))
  expect_identical(unname(z[, "remainder_kind"]), c(3, 3))
  ref <- vapply(x, function(v) .spa_ref(v, s, g = 0.75), numeric(1L))
  expect_equal(unname(z[, "log_p_two_sided"]), ref, tolerance = 1e-9)
  expect_identical(attr(z, "nodes_above_leading"), 0)
  # The order-4 fallback obeys the same condition when two nodes cannot be fit.
  z1 <- .spa_kernel(x, s, lead + c(0.75, 1.125, 0, 0), 4L)
  expect_identical(z1, z)
  # Two nodes with u_2 > s_1^2 and m_2 < 1.
  mu <- node(c(0.2, 1.5), c(3, 0.5))
  z2 <- .spa_kernel(x, s, lead + mu, 4L)
  expect_identical(unname(z2[, "status"]), c(0, 0))
  expect_identical(unname(z2[, "remainder_kind"]), c(2, 2))
  ref2 <- vapply(x, function(v) .spa_ref(v, s, m = c(3, 0.5), u = c(0.2, 1.5)), numeric(1L))
  expect_equal(unname(z2[, "log_p_two_sided"]), ref2, tolerance = 1e-9)
  expect_identical(attr(z2, "nodes_above_leading"), 2)
  # No node above the leading values: the diagnostic is zero.
  z3 <- .spa_kernel(x, s, lead + node(c(0.2, 0.05), c(3, 2)), 4L)
  expect_identical(attr(z3, "nodes_above_leading"), 0)
  # Sweep inconsistent one-node moments, including nodes far above the lead.
  xs <- c(0.5, 5, 20, 60)
  for (u in c(1.01, 1.5, 4, 25)) for (m in c(0.02, 0.5, 1, 3)) {
    zz <- .spa_kernel(xs, s, lead + node(u, m), 2L)
    expect_true(all(zz[, "status"] == 0), info = paste(u, m))
    expect_true(all(zz[, "remainder_kind"] == 3), info = paste(u, m))
    refs <- vapply(xs, function(v) .spa_ref(v, s, g = m * u), numeric(1L))
    expect_true(all(is.finite(refs)))
    expect_true(all(diff(zz[, "log_p_two_sided"]) < 0), info = paste(u, m))
    expect_equal(unname(zz[, "log_p_two_sided"]), refs, tolerance = 1e-7,
                 info = paste(u, m))
  }
})

test_that("inconsistent tiny remainder moments agree with Monte Carlo tails", {
  withr::local_seed(9037L)
  s <- 0.63 * exp(-(0:79) / 3.5)
  lead <- .spa_powers(s)
  x <- 0.5
  U <- numeric(200000L)
  for (a in s) U <- U + a * rnorm(length(U)) * rnorm(length(U))
  E <- rnorm(length(U)) * rnorm(length(U))
  for (mu1 in c(1e-4, 1e-6, 1e-8)) {
    # A genuine small remainder has variance mu1 and negligible higher moments.
    mc <- log10(mean(abs(U + sqrt(mu1) * E) > x))
    for (order in c(2L, 4L)) {
      z <- .spa_kernel(x, s, lead + c(mu1, 1e-3, 0, 0), order)
      expect_identical(unname(z[, "status"]), 0)
      expect_identical(unname(z[, "remainder_kind"]), 3)
      expect_identical(attr(z, "nodes_above_leading"), 0)
      lp <- unname(z[, "log_p_two_sided"]) / log(10)
      expect_gt(lp, -10)
      expect_lt(abs(lp - mc), 0.02)
      expect_identical(z, .spa_kernel(x, s, lead + c(mu1, 0, 0, 0), order))
    }
  }
})

test_that("compatible one-node moments retain their saddlepoint calculation", {
  s <- 0.63 * exp(-(0:79) / 3.5)
  mu <- c(1e-2, 1e-3, 0, 0)
  x <- c(-2, -0.5, 0.5, 2)
  z <- .spa_kernel(x, s, .spa_powers(s) + mu, 2L)
  expect_identical(unname(z[, "remainder_kind"]), rep(1, length(x)))
  expect_identical(unname(z[, "status"]), rep(0, length(x)))
  expect_identical(z, .spa_kernel(x, s, .spa_powers(s) + mu, 4L))
  ref <- vapply(abs(x), function(v) .spa_ref(v, s, m = 0.1, u = 0.1), numeric(1L))
  expect_equal(unname(z[, "log_p_two_sided"]), ref, tolerance = 1e-8)
  # Roundoff in subtracting the leading sums must not reject equality.
  z <- .spa_kernel(0.5, s, .spa_powers(s) + c(1, 0.63^2, 0, 0), 2L)
  expect_identical(unname(z[, "remainder_kind"]), 1)
  ref <- .spa_ref(0.5, s, m = 1 / 0.63^2, u = 0.63^2)
  expect_equal(unname(z[, "log_p_two_sided"]), ref, tolerance = 1e-8)
  above <- .spa_kernel(0.5, s, .spa_powers(s) + c(1, 0.63^2 * (1 + 1e-7), 0, 0), 2L)
  expect_identical(unname(above[, "remainder_kind"]), 3)
  # Equality at the leading squared value is still compatible.
  s <- c(0.5, 0.25)
  z <- .spa_kernel(0.5, s, .spa_powers(s) + c(0.125, 0.03125, 0, 0), 2L)
  expect_identical(unname(z[, "remainder_kind"]), 1)
})

test_that("without positive leading values the existing one-node calculation remains available", {
  Tm <- 3 * 0.2^(1:4)
  x <- c(0.5, 2)
  z <- .spa_kernel(x, numeric(0), Tm, 2L)
  expect_identical(unname(z[, "remainder_kind"]), c(1, 1))
  expect_identical(unname(z[, "status"]), c(0, 0))
  expect_identical(z, .spa_kernel(x, c(0, 0), Tm, 2L))
  ref <- vapply(x, function(v) .spa_ref(v, numeric(0), m = 3, u = 0.2), numeric(1L))
  expect_equal(unname(z[, "log_p_two_sided"]), ref, tolerance = 1e-8)
})

test_that("invalid saddlepoint inputs return a status and no p-value", {
  s <- c(1, 0.6)
  Tm <- .spa_powers(s)
  bad <- mgcvST:::mgcvst_spa_cpp(c(NA, 1), matrix(s, ncol = 1L),
                                 matrix(Tm, 4L, 1L), 4L, 1L)
  expect_identical(unname(bad[1L, "status"]), 1)
  expect_true(is.na(bad[1L, "log_p_two_sided"]))
  nonpositive <- mgcvST:::mgcvst_spa_cpp(1, matrix(s, ncol = 1L),
                                         matrix(c(0, Tm[2:4]), 4L, 1L), 4L, 1L)
  expect_identical(unname(nonpositive[1L, "status"]), 1)
  expect_error(mgcvST:::mgcvst_spa_cpp(1, matrix(s, ncol = 1L), matrix(Tm, 4L, 1L),
                                       3L, 1L), "order")
})

.spa_pair_fixture <- function(n = 7L, q = 9L, seed = 11L) {
  withr::local_seed(seed)
  H <- lapply(seq_len(n), function(g) {
    Z <- matrix(rnorm(q * q), q, q) %*% diag(exp(-seq(0, 2.5, length.out = q) * runif(1L, 0.5, 1.5)))
    tcrossprod(Z) * exp(rnorm(1L, 0, 2))
  })
  a <- matrix(rnorm(q * n), q, n)
  idx <- which(upper.tri(matrix(0, n, n)), arr.ind = TRUE)
  idx <- idx[order(idx[, 1L], idx[, 2L]), , drop = FALSE]
  list(H = H, a = a, left = idx[, 1L], right = idx[, 2L], q = q, n = n)
}

test_that("with k = q the exact pair kernel equals the full-spectrum saddlepoint", {
  f <- .spa_pair_fixture()
  G <- mgcvST:::mgcvst_pair_basis_cpp(f$H, diag(f$q), 1L)
  z <- mgcvST:::mgcvst_pair_spa_cpp(f$H, G, f$a, f$left, f$right, 1L, 4L)
  expect_true(all(z$status == 0L))
  expect_identical(unique(z$remainder_kind), 0L)
  U <- colSums(f$a[, f$left, drop = FALSE] * f$a[, f$right, drop = FALSE])
  expect_equal(z$score, U, tolerance = 1e-12)
  ref <- vapply(seq_along(U), function(p) {
    s <- svd(.spa_sqrtm(f$H[[f$left[p]]]) %*% .spa_sqrtm(f$H[[f$right[p]]]))$d
    .spa_ref(abs(U[p]), s[s > 1e-14 * s[1L]])
  }, numeric(1L))
  expect_equal(z$log_p_two_sided, ref, tolerance = 1e-7)
  # The one-sided tails follow the sign of the score.
  expect_equal(pmin(z$log_p_positive, z$log_p_negative), z$log_p_two_sided - log(2),
               tolerance = 1e-12)
  expect_true(all((z$score > 0) == (z$log_p_positive < z$log_p_negative)))
})

test_that("the four-moment remainder of a shared basis matches an R reference", {
  f <- .spa_pair_fixture(n = 6L, q = 12L, seed = 5L)
  S <- mgcvST:::mgcvst_pair_basis_sum_cpp(f$H)
  V <- eigen(S, symmetric = TRUE)$vectors[, 1:3, drop = FALSE]
  G <- mgcvST:::mgcvst_pair_basis_cpp(f$H, V, 1L)
  # G_g = H_g^{1/2} V of the normalized matrix, symmetric square root.
  H1 <- f$H[[1L]]
  expect_equal(G[[1L]], .spa_sqrtm(H1 / max(abs(H1))) %*% V, tolerance = 1e-10)
  for (order in c(4L, 2L)) {
    z <- mgcvST:::mgcvst_pair_spa_cpp(f$H, G, f$a, f$left, f$right, 1L, order)
    expect_true(all(z$status == 0L))
    U <- z$score
    ref <- vapply(seq_along(U), function(p) {
      i <- f$left[p]
      j <- f$right[p]
      Hi <- f$H[[i]] / max(abs(f$H[[i]]))
      Hj <- f$H[[j]] / max(abs(f$H[[j]]))
      units <- sqrt(max(abs(f$H[[i]])) * max(abs(f$H[[j]])))
      s <- svd(.spa_sqrtm(Hi) %*% .spa_sqrtm(Hj))$d
      shat <- svd(crossprod(.spa_sqrtm(Hi) %*% V, .spa_sqrtm(Hj) %*% V))$d
      mu <- .spa_powers(s) - .spa_powers(shat)
      rem <- if (order == 4L) .spa_rem_two(mu, max(shat)) else .spa_rem_one(mu, max(shat))
      .spa_ref(abs(U[p]) / units, shat, rem$m, rem$u, rem$g)
    }, numeric(1L))
    expect_equal(z$log_p_two_sided, ref, tolerance = 1e-6)
    expect_true(all(z$remainder_kind >= 1L))
  }
})

test_that("a zero-matrix gene gets status 1 and no p-value, and does not disturb the others", {
  f <- .spa_pair_fixture(n = 4L)
  f$H[[2L]][] <- 0
  S <- mgcvST:::mgcvst_pair_basis_sum_cpp(f$H)
  expect_identical(attr(S, "used"), 3L)
  V <- eigen(S, symmetric = TRUE)$vectors[, 1:3, drop = FALSE]
  G <- mgcvST:::mgcvst_pair_basis_cpp(f$H, V, 1L)
  expect_true(all(is.na(G[[2L]])))
  z <- mgcvST:::mgcvst_pair_spa_cpp(f$H, G, f$a, f$left, f$right, 1L, 4L)
  bad <- f$left == 2L | f$right == 2L
  expect_true(all(z$status[bad] == 1L))
  expect_true(all(is.na(z$log_p_two_sided[bad])))
  expect_true(all(z$status[!bad] == 0L))
  expect_true(all(is.finite(z$log_p_two_sided[!bad])))
  none <- mgcvST:::mgcvst_pair_spa_cpp(f$H[-2L], G[-2L], f$a[, -2L], c(1L, 1L, 2L),
                                       c(2L, 3L, 3L), 1L, 4L)
  expect_identical(none$log_p_two_sided[1:3], z$log_p_two_sided[!bad])
})

test_that("pair results are bitwise identical for 1 and 4 threads, and whatever the batching", {
  f <- .spa_pair_fixture(n = 14L, q = 10L, seed = 3L)
  S <- mgcvST:::mgcvst_pair_basis_sum_cpp(f$H)
  # The shared-basis sum is accumulated in feature order whatever the batches.
  head <- mgcvST:::mgcvst_pair_basis_sum_cpp(f$H[1:5])
  batched <- mgcvST:::mgcvst_pair_basis_sum_cpp(f$H[6:14], head)
  attr(S, "used") <- attr(batched, "used") <- NULL
  expect_identical(S, batched)
  V <- eigen(S, symmetric = TRUE)$vectors[, 1:4, drop = FALSE]
  G1 <- mgcvST:::mgcvst_pair_basis_cpp(f$H, V, 1L)
  G4 <- mgcvST:::mgcvst_pair_basis_cpp(f$H, V, 4L)
  expect_identical(G1, G4)
  one <- mgcvST:::mgcvst_pair_spa_cpp(f$H, G1, f$a, f$left, f$right, 1L, 4L)
  four <- mgcvST:::mgcvst_pair_spa_cpp(f$H, G4, f$a, f$left, f$right, 4L, 4L)
  expect_identical(one, four)
  # A pair does not depend on the other pairs of its call.
  last <- length(f$left)
  alone <- mgcvST:::mgcvst_pair_spa_cpp(f$H, G1, f$a, f$left[last], f$right[last], 1L, 4L)
  expect_identical(alone$log_p_two_sided, one$log_p_two_sided[last])
})

test_that("pair results are symmetric in the two genes and the sign follows the score", {
  f <- .spa_pair_fixture(n = 3L, q = 8L, seed = 9L)
  G <- mgcvST:::mgcvst_pair_basis_cpp(f$H, diag(f$q), 1L)
  one <- mgcvST:::mgcvst_pair_spa_cpp(f$H, G, f$a, 1L, 2L, 1L, 4L)
  swapped <- mgcvST:::mgcvst_pair_spa_cpp(rev(f$H[1:2]), rev(G[1:2]), f$a[, 2:1], 1L, 2L, 1L, 4L)
  expect_equal(swapped$log_p_two_sided, one$log_p_two_sided, tolerance = 1e-12)
  a <- f$a
  a[, 2L] <- -a[, 2L]
  negative <- mgcvST:::mgcvst_pair_spa_cpp(f$H, G, a, 1L, 2L, 1L, 4L)
  expect_equal(negative$score, -one$score, tolerance = 1e-12)
  expect_equal(negative$log_p_positive, one$log_p_negative, tolerance = 1e-12)
  expect_equal(negative$log_p_two_sided, one$log_p_two_sided, tolerance = 1e-12)
  # Covariance units do not change the p-value.
  scale <- c(1e-120, 1e120)
  H <- list(f$H[[1L]] * scale[1L], f$H[[2L]] * scale[2L])
  as <- sweep(f$a[, 1:2], 2L, sqrt(scale), "*")
  Gs <- mgcvST:::mgcvst_pair_basis_cpp(H, diag(f$q), 1L)
  scaled <- mgcvST:::mgcvst_pair_spa_cpp(H, Gs, as, 1L, 2L, 1L, 4L)
  expect_equal(scaled$log_p_two_sided, one$log_p_two_sided, tolerance = 1e-9)
})
