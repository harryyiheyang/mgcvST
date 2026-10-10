test_that("raw constrained marginal spectrum equals the projected TAPS reference", {
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
  # Raw whitened coordinates of the sparse kernel: B'QB = I and the constraint
  # direction u = B'g / |B'g| is projected out of the curvature.
  Bw <- backsolve(chol(Q), diag(q))
  u <- as.numeric(crossprod(Bw, g))
  u <- u / sqrt(sum(u^2))
  Pu <- diag(q) - tcrossprod(u)
  M <- Pu %*% crossprod(Bw, B0 %*% Bw) %*% Pu / tau
  lambda_raw <- eigen((M + t(M)) / 2, symmetric = TRUE, only.values = TRUE)$values
  lambda_raw <- lambda_raw[lambda_raw > 1e-12 * max(lambda_raw)]

  Theta <- solve(Qp / norm_projected)
  h_projected <- as.numeric(crossprod(Z, h0))
  B_projected <- crossprod(Z, B0 %*% Z)
  stat_projected <- as.numeric(crossprod(h_projected, Theta %*% h_projected))
  E <- eigen((Theta + t(Theta)) / 2, symmetric = TRUE)
  Theta_sqrt <- E$vectors %*% (sqrt(E$values) * t(E$vectors))
  H <- Theta_sqrt %*% B_projected %*% Theta_sqrt
  lambda <- eigen((H + t(H)) / 2, symmetric = TRUE, only.values = TRUE)$values

  expect_equal(C, Bw %*% Pu %*% t(Bw), tolerance = 1e-10)
  expect_equal(C / tau, Z %*% Theta %*% t(Z), tolerance = 1e-10)
  expect_equal(stat_raw, stat_projected, tolerance = 1e-10)
  expect_length(lambda_raw, q - 1L)
  expect_equal(lambda_raw, lambda, tolerance = 1e-9)
  raw <- mgcvST:::.mgcvst_marginal_davies(
    list(statistic = stat_raw, lambda = lambda_raw), 1e-8, 1e5)
  projected <- mgcvST:::.mgcvst_marginal_davies(
    list(statistic = stat_projected, lambda = lambda), 1e-8, 1e5)
  expect_identical(raw$method_used, "davies")
  expect_equal(raw, projected, tolerance = 1e-10)
})

test_that("INLA Stage 1 p-values match a dense eigenvalue and Davies reference", {
  set.seed(2026100901)
  n <- 30L
  q <- 8L
  k <- 3L
  A0 <- matrix(rnorm(n * q), n, q)
  A0[abs(A0) < 0.5] <- 0
  A <- methods::as(methods::as(Matrix::Matrix(A0, sparse = TRUE), "generalMatrix"),
                   "CsparseMatrix")
  R <- matrix(rnorm(q * q), q, q)
  Q0 <- crossprod(R) + diag(q)
  Q <- methods::as(methods::as(Matrix::Matrix(Q0, sparse = TRUE), "generalMatrix"),
                   "CsparseMatrix")
  g <- runif(q, 0.5, 1.5)
  X <- cbind(1, rnorm(n))
  E <- matrix(rnorm(n * k), n, k) + A0 %*% matrix(rnorm(q * k), q, k) %*% diag(c(0, 0.15, 0.3))
  V <- matrix(runif(n * k, 0.6, 1.6), n, k)
  geometry <- list(A = A, Q = Q, constraint = g)
  fits <- lapply(seq_len(k), function(j) {
    list(working_error = E[, j], working_variance = V[, j])
  })
  marginal <- function() {
    mgcvST:::.inlast_null_marginal(paste0("g", seq_len(k)), geometry, X, fits,
                                   list(random = list()), rep(1, k),
                                   matrix(1, k, 1L), seq_len(k))
  }
  got <- marginal()
  kernel <- mgcvST:::.inlast_sparse_null_batch(geometry, X, E, V, NULL)

  # Dense reference built from M on the constrained complement of g.
  G <- qr.Q(qr(matrix(g, ncol = 1L)), complete = TRUE)[, -1L, drop = FALSE]
  Qp <- crossprod(G, Q0 %*% G)
  tau <- 1 / norm(Qp, "F")
  T <- G %*% backsolve(chol(Qp), diag(q - 1L)) / sqrt(tau)
  statistic <- numeric(k)
  lambda <- vector("list", k)
  for (j in seq_len(k)) {
    W <- 1 / V[, j]
    Vp <- solve(crossprod(X, W * X))
    P <- diag(W) - (W * X) %*% Vp %*% t(W * X)
    a <- as.numeric(crossprod(T, crossprod(A0, P %*% E[, j])))
    M <- crossprod(T, crossprod(A0, P %*% A0) %*% T)
    ev <- eigen((M + t(M)) / 2, symmetric = TRUE, only.values = TRUE)$values
    lambda[[j]] <- ev[ev > 1e-12 * max(ev)]
    statistic[j] <- sum(a^2)
    reference <- CompQuadForm::davies(statistic[j], lambda[[j]], lim = 1e5, acc = 1e-8)
    expect_identical(reference$ifault, 0L)
    expect_equal(kernel[[j]]$lambda, lambda[[j]], tolerance = 1e-10)
    expect_equal(got$statistic[j], statistic[j], tolerance = 1e-10)
    expect_equal(got$p_value[j], reference$Qq, tolerance = 1e-8)
  }
  expect_true(all(got$p_value > 1e-6 & got$p_value < 1))
  expect_identical(got$method_requested, rep("davies", k))
  expect_identical(got$method_used, rep("davies", k))
  expect_identical(got$fallback_used, rep(FALSE, k))
  expect_identical(got$davies_ifault, rep(0L, k))
  expect_true(all(is.na(got$error_message)))
  expect_null(mgcvST:::mgcvst_inla_sparse_batch_cpp(
    geometry$A, geometry$Q, g, X, E, V, rep(1, k))[[1L]]$lambda)

  local({
    testthat::local_mocked_bindings(davies = function(...) list(Qq = 0, ifault = 0L),
                                    .package = "CompQuadForm")
    failed <- marginal()
    expect_identical(failed$method_used, rep("saddlepoint", k))
    expect_identical(failed$fallback_used, rep(TRUE, k))
    expect_equal(failed$p_value, vapply(seq_len(k), function(j) {
      mgcvST:::.mgcvst_marginal_saddlepoint(statistic[j], lambda[[j]])
    }, numeric(1L)), tolerance = 1e-8)
  })
  local({
    testthat::local_mocked_bindings(davies = function(...) list(Qq = 0.4, ifault = 1L),
                                    .package = "CompQuadForm")
    kept <- marginal()
    expect_identical(kept$p_value, rep(0.4, k))
    expect_identical(kept$method_used, rep("davies", k))
    expect_identical(kept$fallback_used, rep(FALSE, k))
    expect_identical(kept$davies_ifault, rep(1L, k))
  })
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
  got <- mgcvST:::.inlast_sparse_observation_basis(fit)

  Q2 <- qr.Q(qr(matrix(g, ncol = 1L)), complete = TRUE)[, -1L, drop = FALSE]
  qp <- crossprod(Q2, Q0 %*% Q2)
  ep <- eigen(qp, symmetric = TRUE)
  Bp <- Q2 %*% (ep$vectors %*% (1 / sqrt(ep$values) * t(ep$vectors)))
  ref <- crossprod(Bp, crossprod(A0) %*% Bp)
  er <- eigen(ref, symmetric = TRUE)
  val <- er$values
  vec <- er$vectors
  # All q - 1 directions of the constrained field are kept.
  r <- q - 1L
  cref <- Bp %*% vec[, seq_len(r), drop = FALSE]

  expect_equal(got$rank, r)
  expect_identical(got$kind, "full_rank")
  expect_null(got$coverage)
  expect_null(got$tail)
  expect_equal(got$values, val, tolerance = 2e-10)
  expect_equal(tcrossprod(got$basis), tcrossprod(cref), tolerance = 2e-10)
  expect_equal(crossprod(got$basis, Q0 %*% got$basis), diag(r),
    tolerance = 2e-10)
})

