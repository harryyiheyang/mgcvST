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

test_that("C++ sparse INLA states equal dense projected score references", {
  set.seed(2026091302)
  n <- 18L
  q <- 7L
  p <- 3L
  A0 <- matrix(rnorm(n * q), n, q)
  A0[abs(A0) < 0.7] <- 0
  A <- Matrix::Matrix(A0, sparse = TRUE)
  X <- cbind(1, rnorm(n))
  R <- matrix(rnorm(q * q), q, q)
  Q <- Matrix::Matrix(crossprod(R) + diag(q), sparse = TRUE)
  g <- runif(q, 0.5, 1.5)
  E <- matrix(rnorm(n * p), n, p)
  D <- matrix(runif(n * p, 0.4, 1.8), n, p)
  tau <- c(0.7, 1.1, 1.6)
  fit <- list(
    score_sparse = list(A = A, Q = Q, constraint = g,
      target = "global", sp_index = 1L),
    dispersion = rep(1, p), smoothing_parameters = matrix(tau, p, 1L),
    working_error = E, working_variance = D,
    geometry = list(nuisance_design = X), feature_id = paste0("g", seq_len(p))
  )

  dense_state <- function(f) {
    W <- 1 / D[, f]
    K <- crossprod(A0, W * A0)
    L <- crossprod(A0, W * X)
    tvec <- as.numeric(crossprod(A0, W * E[, f]))
    xte <- as.numeric(crossprod(X, W * E[, f]))
    H <- tau[f] * as.matrix(Q) + K
    Hi <- solve(H)
    Hig <- as.numeric(Hi %*% g)
    solve_H <- function(B) {
      ans <- Hi %*% B
      ans - tcrossprod(Hig, as.numeric(crossprod(g, ans)) /
        as.numeric(crossprod(g, Hig)))
    }
    SL <- solve_H(L)
    St <- as.numeric(solve_H(tvec))
    U <- L - K %*% SL
    J <- crossprod(X, W * X) - crossprod(L, SL)
    Vp <- solve(J)
    qvec <- xte - as.numeric(crossprod(L, St))
    h <- tvec - as.numeric(K %*% St) - as.numeric(U %*% Vp %*% qvec)
    Z <- qr.Q(qr(matrix(g, ncol = 1L)), complete = TRUE)[, -1L, drop = FALSE]
    Qp <- crossprod(Z, as.matrix(Q) %*% Z)
    T <- Z %*% backsolve(chol(Qp), diag(q - 1L)) / sqrt(tau[f])
    KT <- K %*% T
    SKT <- solve_H(KT)
    M <- crossprod(T, KT) - crossprod(KT, SKT)
    UtT <- crossprod(U, T)
    M <- M - crossprod(UtT, Vp %*% UtT)
    list(a = as.numeric(crossprod(T, h)), M = (M + t(M)) / 2, Vp = Vp)
  }

  ref <- lapply(seq_len(p), dense_state)
  one <- mgcvST:::.inlast_sparse_batch(fit, seq_len(p), threads = 1L)
  two <- mgcvST:::.inlast_sparse_batch(fit, seq_len(p), threads = 2L)
  scores_ref <- vapply(seq_len(p), function(j) sum(ref[[j]]$a^2), numeric(1L))
  scores_one <- vapply(one, `[[`, numeric(1L), "statistic")
  expect_equal(scores_one, scores_ref, tolerance = 2e-9)
  expect_equal(scores_one, vapply(two, `[[`, numeric(1L), "statistic"), tolerance = 2e-12)
  expect_equal(lapply(one, `[[`, "expected_vp"), lapply(ref, `[[`, "Vp"), tolerance = 2e-10)
  expect_equal(lapply(two, `[[`, "a"), lapply(one, `[[`, "a"), tolerance = 2e-12)
  expect_equal(lapply(two, `[[`, "M"), lapply(one, `[[`, "M"), tolerance = 2e-12)
  expect_true(all(is.finite(vapply(one, `[[`, numeric(1L), "constraint_diagnostic"))))
  expect_equal(
    vapply(two, `[[`, numeric(1L), "constraint_diagnostic"),
    vapply(one, `[[`, numeric(1L), "constraint_diagnostic"),
    tolerance = 2e-12
  )
  expect_true(all(vapply(one, function(z) z$normalization == q - 1L, logical(1L))))

  pairs <- rbind(c(1L, 2L), c(1L, 3L), c(2L, 3L))
  for (j in seq_len(nrow(pairs))) {
    i1 <- pairs[j, 1L]
    i2 <- pairs[j, 2L]
    U_ref <- as.numeric(crossprod(ref[[i1]]$a, ref[[i2]]$a))
    U_cpp <- as.numeric(crossprod(one[[i1]]$a, one[[i2]]$a))
    expect_equal(U_cpp, U_ref, tolerance = 2e-9)
    G_ref <- ref[[i1]]$M %*% ref[[i2]]$M
    G_cpp <- one[[i1]]$M %*% one[[i2]]$M
    tr_ref <- tr_cpp <- numeric(4L)
    Gr <- diag(nrow(G_ref))
    Gc <- diag(nrow(G_cpp))
    for (k in seq_len(4L)) {
      Gr <- Gr %*% G_ref
      Gc <- Gc %*% G_cpp
      tr_ref[k] <- sum(diag(Gr))
      tr_cpp[k] <- sum(diag(Gc))
    }
    expect_equal(tr_cpp, tr_ref, tolerance = 5e-8)
    expect_equal(
      mgcvST:::.liu_squared_score_moments(U_cpp, tr_cpp[1L], tr_cpp[2L],
        tr_cpp[3L], tr_cpp[4L])$p_value,
      mgcvST:::.liu_squared_score_moments(U_ref, tr_ref[1L], tr_ref[2L],
        tr_ref[3L], tr_ref[4L])$p_value,
      tolerance = 2e-8
    )
  }

  only <- mgcvST:::.inlast_sparse_batch(
    fit, seq_len(p), threads = 2L, score_only = TRUE
  )
  expect_equal(lapply(only, `[[`, "a"), lapply(one, `[[`, "a"), tolerance = 2e-12)
  expect_true(all(vapply(only, function(z) is.null(z$M), logical(1L))))
  expect_true(all(vapply(only, function(z) is.null(z$moments), logical(1L))))
  prepared <- mgcvST:::.inlast_sparse_prepare(fit)
  again <- mgcvST:::.inlast_sparse_prepare(prepared)
  expect_identical(prepared$score_sparse$cache$prepared,
                   again$score_sparse$cache$prepared)
  units <- mgcvST:::.inlast_sparse_units(prepared, seq_len(p), threads = 2L)
  expect_true(all(vapply(units, function(z) is.null(z$M), logical(1L))))
  expect_true(all(vapply(units, function(z) inherits(z$K, "sparseMatrix") &&
    inherits(z$H_L, "sparseMatrix"), logical(1L))))
  restored <- unserialize(serialize(units, NULL))
  state <- mgcvST:::.inlast_sparse_materialize(prepared, restored, threads = 2L)
  expect_equal(lapply(state, `[[`, "M"), lapply(one, `[[`, "M"), tolerance = 2e-10)
  expect_equal(lapply(state, `[[`, "a"), lapply(one, `[[`, "a"), tolerance = 2e-12)
})

test_that("C++ sparse marginal moments equal the projected TAPS spectrum", {
  set.seed(2026091303)
  n <- 18L
  q <- 7L
  p <- 3L
  A0 <- matrix(rnorm(n * q), n, q)
  A0[abs(A0) < 0.6] <- 0
  A <- Matrix::Matrix(A0, sparse = TRUE)
  X <- cbind(1, rnorm(n))
  R <- matrix(rnorm(q * q), q, q)
  Q0 <- crossprod(R) + diag(q)
  Q <- Matrix::Matrix(Q0, sparse = TRUE)
  g <- runif(q, 0.5, 1.5)
  E <- matrix(rnorm(n * p), n, p)
  D <- matrix(runif(n * p, 0.5, 1.7), n, p)
  fit <- list(
    score_sparse = list(A = A, Q = Q, constraint = g,
      target = "global", sp_index = 1L),
    dispersion = rep(1, p), smoothing_parameters = matrix(1, p, 1L),
    working_error = E, working_variance = D,
    geometry = list(nuisance_design = X), feature_id = paste0("g", seq_len(p))
  )
  null_state <- list(
    working_error = E, working_variance = D, nuisance_precision = NULL
  )
  z1 <- mgcvST:::.inlast_sparse_null_batch(
    fit$score_sparse, fit$geometry$nuisance_design, null_state,
    seq_len(p), threads = 1L
  )
  z2 <- mgcvST:::.inlast_sparse_null_batch(
    fit$score_sparse, fit$geometry$nuisance_design, null_state,
    seq_len(p), threads = 2L
  )
  ug <- g / sqrt(sum(g^2))
  Q2 <- Q0 %*% Q0
  normQp <- sqrt(sum(Q0^2) - 2 * as.numeric(crossprod(ug, Q2 %*% ug)) +
    as.numeric(crossprod(ug, Q0 %*% ug))^2)
  Z <- qr.Q(qr(matrix(g, ncol = 1L)), complete = TRUE)[, -1L, drop = FALSE]
  Theta <- solve(crossprod(Z, Q0 %*% Z) / normQp)

  for (f in seq_len(p)) {
    W <- 1 / D[, f]
    K <- crossprod(A0, W * A0)
    L <- crossprod(A0, W * X)
    J <- crossprod(X, W * X)
    tvec <- as.numeric(crossprod(A0, W * E[, f]))
    xte <- as.numeric(crossprod(X, W * E[, f]))
    h0 <- tvec - as.numeric(L %*% solve(J, xte))
    B0 <- K - L %*% solve(J, t(L))
    hp <- as.numeric(crossprod(Z, h0))
    Bp <- crossprod(Z, B0 %*% Z)
    stat <- as.numeric(crossprod(hp, Theta %*% hp))
    Et <- eigen((Theta + t(Theta)) / 2, symmetric = TRUE)
    Ts <- Et$vectors %*% (sqrt(Et$values) * t(Et$vectors))
    H <- Ts %*% Bp %*% Ts
    lambda <- eigen((H + t(H)) / 2, symmetric = TRUE, only.values = TRUE)$values
    moments <- vapply(seq_len(4L), function(k) sum(lambda^k), numeric(1L))
    expect_null(z1[[f]]$M)
    expect_null(z2[[f]]$M)
    expect_equal(z1[[f]]$statistic, stat, tolerance = 2e-9)
    expect_equal(z1[[f]]$moments, moments, tolerance = 5e-8)
    expect_equal(z2[[f]]$statistic, z1[[f]]$statistic, tolerance = 2e-12)
    expect_equal(z2[[f]]$moments, z1[[f]]$moments, tolerance = 2e-12)
    expect_equal(
      mgcvST:::.mgcvst_marginal_liu(z1[[f]]$statistic, z1[[f]]$moments),
      mgcvST:::.mgcvst_marginal_liu(stat, moments), tolerance = 2e-8
    )
  }
})
