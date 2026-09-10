# Experimental conditioned-kernel score state for one fixed-kappa SPDE target.
#
# This prototype preserves the production native-nuisance-Vp definition while
# avoiding an observation-by-mesh score factor.  Its dense work is restricted
# to the mesh dimension.  It supports one target and fixed-effect nuisance
# columns only; additional random nuisance blocks are deliberately out of scope.

.inlast_conditioned_sparse_solver <- function(H, g) {
  factor <- Matrix::Cholesky(
    Matrix::forceSymmetric(H), LDL = FALSE, super = FALSE
  )
  solve_H <- function(rhs) {
    rhs <- if (is.null(dim(rhs))) matrix(as.numeric(rhs), ncol = 1L) else rhs
    as.matrix(Matrix::solve(factor, rhs))
  }
  Hinv_g <- as.numeric(solve_H(g))
  denominator <- sum(g * Hinv_g)
  if (!is.finite(denominator) || denominator <= 0) {
    stop("The mean constraint has a non-positive H-inverse norm.")
  }
  function(rhs) {
    vector <- is.null(dim(rhs))
    answer <- solve_H(rhs)
    multiplier <- as.numeric(crossprod(g, answer)) / denominator
    answer <- answer - tcrossprod(Hinv_g, multiplier)
    if (vector) as.numeric(answer) else answer
  }
}

inlast_conditioned_sparse_state <- function(
    raw_A, raw_Q, g, projection, working_error, working_variance, tau,
    nuisance_X, nuisance_Vp, tolerance = 1e-9) {
  if (!requireNamespace("Matrix", quietly = TRUE)) stop("Matrix is required.")
  A <- methods::as(raw_A, "CsparseMatrix")
  Q <- Matrix::forceSymmetric(methods::as(raw_Q, "CsparseMatrix"))
  if (any(!is.finite(A@x)) || any(!is.finite(Q@x))) {
    stop("raw_A and raw_Q must be finite sparse matrices.")
  }
  n <- nrow(A)
  m <- ncol(A)
  if (!all(dim(Q) == c(m, m)) ||
      !isTRUE(Matrix::isSymmetric(Q, tol = tolerance))) {
    stop("raw_Q must be symmetric with ncol(raw_A) rows and columns.")
  }
  Z <- as.matrix(projection)
  storage.mode(Z) <- "double"
  if (!all(dim(Z) == c(m, m - 1L)) || any(!is.finite(Z)) ||
      qr(Z, tol = tolerance)$rank != ncol(Z)) {
    stop("projection must be a finite full-rank m by (m - 1) matrix.")
  }
  g <- as.numeric(g)
  expected_g <- as.numeric(Matrix::crossprod(A, rep.int(1 / n, n)))
  if (length(g) != m || any(!is.finite(g)) || !isTRUE(all.equal(
      g, expected_g, tolerance = tolerance, check.attributes = FALSE))) {
    stop("g must equal crossprod(raw_A, 1 / n).")
  }
  if (max(abs(crossprod(g, Z))) >
      tolerance * max(1, max(abs(g)), max(abs(Z)))) {
    stop("projection does not obey the observation mean constraint.")
  }
  e <- as.numeric(working_error)
  D <- as.numeric(working_variance)
  if (length(e) != n || any(!is.finite(e))) {
    stop("working_error must be finite with one value per observation.")
  }
  if (length(D) != n || any(!is.finite(D)) || any(D <= 0)) {
    stop("working_variance must be positive with one value per observation.")
  }
  tau <- as.numeric(tau)
  if (length(tau) != 1L || !is.finite(tau) || tau <= 0) {
    stop("tau must be one positive finite original-Q precision.")
  }
  X <- as.matrix(nuisance_X)
  storage.mode(X) <- "double"
  if (nrow(X) != n || any(!is.finite(X)) ||
      (ncol(X) && qr(X, tol = tolerance)$rank != ncol(X))) {
    stop("nuisance_X must be finite, aligned, and full column rank.")
  }
  Vp <- as.matrix(nuisance_Vp)
  storage.mode(Vp) <- "double"
  if (!all(dim(Vp) == c(ncol(X), ncol(X))) || any(!is.finite(Vp)) ||
      !isTRUE(all.equal(Vp, t(Vp), tolerance = tolerance))) {
    stop("nuisance_Vp must be a finite symmetric nuisance covariance.")
  }
  Vp <- (Vp + t(Vp)) / 2

  W <- 1 / D
  WA <- A * W
  K <- Matrix::forceSymmetric(Matrix::crossprod(A, WA))
  L <- if (ncol(X)) as.matrix(Matrix::crossprod(A, W * X)) else
    matrix(numeric(), m, 0L)
  tvec <- as.numeric(Matrix::crossprod(A, W * e))
  H <- Matrix::forceSymmetric(tau * Q + K)
  solve_S <- .inlast_conditioned_sparse_solver(H, g)

  # Only p + 1 constrained solves are needed for h and the nuisance terms.
  solved_small <- solve_S(cbind(L, tvec))
  SL <- if (ncol(X)) solved_small[, seq_len(ncol(X)), drop = FALSE] else
    matrix(numeric(), m, 0L)
  St <- as.numeric(solved_small[, ncol(X) + 1L])
  U <- if (ncol(X)) L - as.matrix(K %*% SL) else
    matrix(numeric(), m, 0L)
  q <- if (ncol(X)) {
    as.numeric(crossprod(X, W * e) - crossprod(L, St))
  } else numeric()
  h <- tvec - as.numeric(K %*% St)
  if (ncol(X)) h <- h - as.numeric(U %*% Vp %*% q)

  Qprojected <- crossprod(Z, as.matrix(Q %*% Z))
  Qprojected <- (Qprojected + t(Qprojected)) / 2
  R <- chol(Qprojected)
  # Qprojected = R'R, so R^-1 R^-T is its covariance and Z R^-1 is
  # the coefficient-space factor used by production's B R^-1.
  T <- Z %*% backsolve(R, diag(ncol(Z))) / sqrt(tau)
  KT <- as.matrix(K %*% T)
  SKT <- solve_S(KT)

  a <- as.numeric(crossprod(T, h))
  M <- crossprod(T, KT) - crossprod(KT, SKT)
  if (ncol(X)) {
    UtT <- crossprod(U, T)
    M <- M - crossprod(UtT, Vp %*% UtT)
  }
  M <- (M + t(M)) / 2

  XvinvX <- if (ncol(X)) {
    crossprod(X, W * X) - crossprod(L, SL)
  } else matrix(numeric(), 0L, 0L)
  nuisance_identity_error <- if (ncol(X)) {
    max(abs(XvinvX - XvinvX %*% Vp %*% XvinvX))
  } else 0
  nuisance_identity_relative_error <- nuisance_identity_error /
    max(1, max(abs(XvinvX)))
  constraint_error <- max(
    abs(sum(g * St)),
    if (ncol(X)) max(abs(crossprod(g, SL))) else 0,
    max(abs(crossprod(g, SKT)))
  )
  list(
    a = a, M = M, width = ncol(Z),
    audit = list(nuisance_information = XvinvX),
    diagnostics = list(
      constraint_solve_error = constraint_error,
      nuisance_identity_error = nuisance_identity_error,
      nuisance_identity_relative_error = nuisance_identity_relative_error,
      minimum_M_eigenvalue = min(eigen(M, symmetric = TRUE,
                                       only.values = TRUE)$values),
      minimum_M_relative_eigenvalue = min(eigen(
        M, symmetric = TRUE, only.values = TRUE
      )$values) / max(1, max(abs(M))),
      observation_dimension = n,
      mesh_dimension = m,
      nuisance_dimension = ncol(X),
      dense_observation_by_mesh_constructed = FALSE,
      dense_mesh_factor_dimension = dim(T),
      sparse_K_nonzeros = Matrix::nnzero(K)
    )
  )
}

inlast_conditioned_sparse_pair <- function(state1, state2,
                                            method = c("liu", "davies")) {
  method <- match.arg(method)
  if (length(state1$a) != length(state2$a) ||
      !all(dim(state1$M) == dim(state2$M))) {
    stop("The conditioned score states are not aligned.")
  }
  signed_score <- as.numeric(crossprod(state1$a, state2$a))
  calibration <- mgcvST::rkhs_score_calibrate(
    signed_score, state1$M, state2$M, method = method
  )
  c(list(signed_score = signed_score, statistic = signed_score^2), calibration)
}
