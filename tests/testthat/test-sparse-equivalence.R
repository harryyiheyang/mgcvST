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
  fit <- spdeST:::.spdest_sparse_prepare(fit)
  got <- spdeST:::.spdest_sparse_observation_basis(fit)

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

