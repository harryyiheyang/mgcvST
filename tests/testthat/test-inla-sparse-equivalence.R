test_that("raw constrained marginal traces equal the projected TAPS reference", {
  set.seed(2026091301)
  n <- 18L
  q <- 7L
  k <- 2L
  A <- matrix(rnorm(n * q), n, q)
  X <- cbind(1, rnorm(n))
  W <- runif(n, 0.4, 1.8)
  e <- rnorm(n)
  R <- matrix(rnorm(q * q), q, q)
  Q <- crossprod(R) + diag(q)
  g <- runif(q, 0.5, 1.5)

  K <- crossprod(A, W * A)
  L <- crossprod(A, W * X)
  J <- crossprod(X, W * X)
  t <- as.numeric(crossprod(A, W * e))
  xe <- as.numeric(crossprod(X, W * e))
  h0 <- t - as.numeric(L %*% solve(J, xe))
  B0 <- K - L %*% solve(J, t(L))

  ug <- g / sqrt(sum(g^2))
  Q2 <- Q %*% Q
  norm_formula <- sqrt(sum(Q^2) - 2 * as.numeric(crossprod(ug, Q2 %*% ug)) +
    as.numeric(crossprod(ug, Q %*% ug))^2)
  Z <- qr.Q(qr(matrix(g, ncol = 1L)), complete = TRUE)[, -1L, drop = FALSE]
  Qp <- crossprod(Z, Q %*% Z)
  norm_projected <- norm(Qp, "F")
  expect_equal(norm_formula, norm_projected, tolerance = 1e-12)

  Qi <- solve(Q)
  Qig <- as.numeric(Qi %*% g)
  C <- Qi - tcrossprod(Qig) / as.numeric(crossprod(g, Qig))
  tau <- 1 / norm_projected
  stat_raw <- as.numeric(crossprod(h0, C %*% h0)) / tau
  G <- C %*% B0 / tau
  moments_raw <- numeric(4L)
  Gk <- diag(q)
  for (j in seq_len(4L)) {
    Gk <- Gk %*% G
    moments_raw[j] <- sum(diag(Gk))
  }

  Theta <- solve(Qp / norm_projected)
  h_projected <- as.numeric(crossprod(Z, h0))
  B_projected <- crossprod(Z, B0 %*% Z)
  stat_projected <- as.numeric(crossprod(h_projected, Theta %*% h_projected))
  E <- eigen((Theta + t(Theta)) / 2, symmetric = TRUE)
  Theta_sqrt <- E$vectors %*% (sqrt(E$values) * t(E$vectors))
  H <- Theta_sqrt %*% B_projected %*% Theta_sqrt
  lambda <- eigen((H + t(H)) / 2, symmetric = TRUE, only.values = TRUE)$values
  moments_projected <- vapply(seq_len(4L), function(j) sum(lambda^j), numeric(1L))

  expect_equal(C / tau, Z %*% Theta %*% t(Z), tolerance = 1e-10)
  expect_equal(stat_raw, stat_projected, tolerance = 1e-10)
  expect_equal(moments_raw, moments_projected, tolerance = 1e-9)
  expect_equal(
    mgcvST:::.mgcvst_marginal_liu(stat_raw, moments_raw),
    mgcvST:::.mgcvst_marginal_liu(stat_projected, moments_projected),
    tolerance = 1e-10
  )
})

test_that("sparse observation basis matches an independent QR reference", {
  set.seed(2026092201)
  n <- 19L
  q <- 8L
  A0 <- matrix(rnorm(n * q), n, q)
  A0[abs(A0) < 0.55] <- 0
  A <- Matrix::Matrix(A0, sparse = TRUE)
  z <- matrix(rnorm(q * q), q, q)
  Q0 <- crossprod(z) + diag(q)
  Q <- Matrix::Matrix(Q0, sparse = TRUE)
  g <- as.numeric(colMeans(A0))
  fit <- list(score_sparse = list(A = A, Q = Q, constraint = g))
  fit <- mgcvST:::.inlast_sparse_prepare(fit)
  got <- mgcvST:::.inlast_sparse_observation_basis(fit, coverage = 0.995)

  Q2 <- qr.Q(qr(matrix(g, ncol = 1L)), complete = TRUE)[, -1L, drop = FALSE]
  qp <- crossprod(Q2, Q0 %*% Q2)
  ep <- eigen(qp, symmetric = TRUE)
  Bp <- Q2 %*% (ep$vectors %*% (1 / sqrt(ep$values) * t(ep$vectors)))
  ref <- crossprod(Bp, crossprod(A0) %*% Bp)
  er <- eigen(ref, symmetric = TRUE)
  val <- er$values
  vec <- er$vectors
  keep <- cumsum(val) / sum(val)
  r <- which(keep >= 0.995)[1L]
  cref <- Bp %*% vec[, seq_len(r), drop = FALSE]

  expect_equal(got$rank, r)
  expect_equal(got$values, val, tolerance = 2e-10)
  expect_equal(tcrossprod(got$basis), tcrossprod(cref), tolerance = 2e-10)
  expect_equal(got$tail, 1 - keep[r], tolerance = 2e-10)
  expect_equal(crossprod(got$basis, Q0 %*% got$basis), diag(r),
    tolerance = 2e-10)

  all <- mgcvST:::.inlast_sparse_observation_basis(
    fit, coverage = 0.995, full_rank = TRUE
  )
  expect_equal(all$rank, q - 1L)
  expect_equal(crossprod(g, all$basis), matrix(0, 1L, q - 1L),
    tolerance = 2e-12)
  expect_equal(tcrossprod(all$basis), Q2 %*% solve(qp, t(Q2)),
    tolerance = 2e-10)
})

