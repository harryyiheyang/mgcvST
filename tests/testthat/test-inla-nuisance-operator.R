test_that("multiple iid nuisance penalties match dense Woodbury operators", {
  set.seed(2026092001)
  n <- 18L
  q <- 7L
  p1 <- 3L
  p2 <- 4L
  features <- 2L
  A0 <- matrix(rnorm(n * q), n, q)
  A0[abs(A0) < 0.6] <- 0
  A <- Matrix::Matrix(A0, sparse = TRUE)
  A <- methods::as(methods::as(A, "generalMatrix"), "CsparseMatrix")
  R <- matrix(rnorm(q * q), q, q)
  Q0 <- crossprod(R) + diag(q)
  Q <- Matrix::Matrix(Q0, sparse = TRUE)
  Q <- methods::as(methods::as(Q, "generalMatrix"), "CsparseMatrix")
  g <- runif(q, 0.5, 1.5)
  X <- cbind(1, rnorm(n))
  Z1 <- model.matrix(~ factor(rep(seq_len(p1), length.out = n)) - 1)
  Z2 <- model.matrix(~ factor(rep(c(1L, 3L, 2L, 4L), length.out = n)) - 1)
  U <- cbind(X, Z1, Z2)
  E <- matrix(rnorm(n * features), n, features)
  D <- matrix(runif(n * features, 0.5, 1.7), n, features)
  tau <- c(0.8, 1.3)
  tau_b <- matrix(c(1.7, 0.9, 2.1, 1.2), nrow = 2L, byrow = TRUE)
  Sn <- mgcvST:::.inlast_iid_nuisance_precision(
    ncol(X), c(ncol(Z1), ncol(Z2)), tau_b, features
  )

  observed <- mgcvST:::mgcvst_inla_sparse_batch_cpp(
    A, Q, g, U, E, D, tau, threads = 1L,
    nuisance_precision = Sn
  )
  units <- mgcvST:::mgcvst_inla_sparse_units_cpp(
    A, Q, g, U, E, D, tau, threads = 1L,
    nuisance_precision = Sn
  )

  for (f in seq_len(features)) {
    W <- diag(1 / D[, f])
    Qi <- solve(Q0)
    Qig <- as.numeric(Qi %*% g)
    Cq <- Qi - tcrossprod(Qig) / as.numeric(crossprod(g, Qig))
    Vtarget <- diag(D[, f]) + A0 %*% (Cq / tau[f]) %*% t(A0)
    Vti <- solve(Vtarget)
    S <- diag(Sn[, f])
    J <- crossprod(U, Vti %*% U) + S
    Vp <- solve(J)
    P <- Vti - Vti %*% U %*% Vp %*% t(U) %*% Vti
    h <- as.numeric(crossprod(A0, P %*% E[, f]))
    B <- crossprod(A0, P %*% A0)
    G <- qr.Q(qr(matrix(g, ncol = 1L)), complete = TRUE)[, -1L, drop = FALSE]
    Qp <- crossprod(G, Q0 %*% G)
    T <- G %*% backsolve(chol(Qp), diag(q - 1L)) / sqrt(tau[f])
    a <- as.numeric(crossprod(T, h))
    M <- crossprod(T, B %*% T)

    expect_equal(unname(observed[[f]]$expected_vp), unname(Vp), tolerance = 2e-9)
    expect_equal(units[[f]]$nuisance_score,
                 as.numeric(crossprod(U, Vti %*% E[, f])), tolerance = 2e-9)
    expect_equal(observed[[f]]$statistic, sum(a^2), tolerance = 2e-8)
    expect_equal(sum(observed[[f]]$M^2), sum(M^2), tolerance = 2e-7)
    expect_equal(sum(diag(observed[[f]]$M)), sum(diag(M)), tolerance = 2e-8)
  }
})

test_that("iid nuisance marginal null keeps the nuisance covariance", {
  set.seed(2026092002)
  n <- 16L
  q <- 6L
  p <- 4L
  A0 <- matrix(rnorm(n * q), n, q)
  A0[abs(A0) < 0.5] <- 0
  A <- Matrix::Matrix(A0, sparse = TRUE)
  A <- methods::as(methods::as(A, "generalMatrix"), "CsparseMatrix")
  R <- matrix(rnorm(q * q), q, q)
  Q0 <- crossprod(R) + diag(q)
  Q <- Matrix::Matrix(Q0, sparse = TRUE)
  Q <- methods::as(methods::as(Q, "generalMatrix"), "CsparseMatrix")
  g <- runif(q, 0.4, 1.4)
  X <- cbind(1, rnorm(n))
  Z <- model.matrix(~ factor(rep(seq_len(p), length.out = n)) - 1)
  U <- cbind(X, Z)
  e <- matrix(rnorm(n), n, 1L)
  d <- matrix(runif(n, 0.6, 1.5), n, 1L)
  tau <- 1 / 2.3
  tau_b <- 1.4
  Sn <- mgcvST:::.inlast_iid_nuisance_precision(
    ncol(X), ncol(Z), matrix(tau_b), 1L
  )
  observed <- mgcvST:::mgcvst_inla_sparse_batch_cpp(
    A, Q, g, U, e, d, tau, threads = 1L, null_target = TRUE,
    nuisance_precision = Sn
  )[[1L]]

  W <- diag(1 / d[, 1L])
  S <- diag(Sn[, 1L])
  Vp <- solve(crossprod(U, W %*% U) + S)
  P <- W - W %*% U %*% Vp %*% t(U) %*% W
  h <- as.numeric(crossprod(A0, P %*% e[, 1L]))
  B <- crossprod(A0, P %*% A0)
  G <- qr.Q(qr(matrix(g, ncol = 1L)), complete = TRUE)[, -1L, drop = FALSE]
  Qp <- crossprod(G, Q0 %*% G)
  T <- G %*% backsolve(chol(Qp), diag(q - 1L)) / sqrt(tau)
  a <- as.numeric(crossprod(T, h))
  M <- crossprod(T, B %*% T)
  moments <- numeric(4L)
  Mk <- diag(nrow(M))
  for (k in seq_len(4L)) {
    Mk <- Mk %*% M
    moments[k] <- sum(diag(Mk))
  }

  expect_equal(unname(observed$expected_vp), unname(Vp), tolerance = 2e-10)
  expect_equal(observed$statistic, sum(a^2), tolerance = 2e-8)
  expect_equal(observed$moments, moments, tolerance = 2e-7)
})

test_that("zero nuisance precision preserves the old sparse score", {
  set.seed(2026092003)
  n <- 14L
  q <- 6L
  A <- Matrix::rsparsematrix(n, q, 0.7)
  A <- methods::as(methods::as(A, "generalMatrix"), "CsparseMatrix")
  Q <- Matrix::Matrix(crossprod(matrix(rnorm(q * q), q, q)) + diag(q),
                      sparse = TRUE)
  Q <- methods::as(methods::as(Q, "generalMatrix"), "CsparseMatrix")
  g <- runif(q, 0.5, 1.5)
  X <- cbind(1, rnorm(n))
  E <- matrix(rnorm(2L * n), n, 2L)
  D <- matrix(runif(2L * n, 0.7, 1.4), n, 2L)
  tau <- c(0.9, 1.2)
  old <- mgcvST:::mgcvst_inla_sparse_batch_cpp(
    A, Q, g, X, E, D, tau, threads = 1L
  )
  new <- mgcvST:::mgcvst_inla_sparse_batch_cpp(
    A, Q, g, X, E, D, tau, threads = 1L,
    nuisance_precision = matrix(0, ncol(X), ncol(E))
  )
  expect_identical(old, new)
})
