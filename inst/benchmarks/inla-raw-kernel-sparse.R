# Experimental raw-kernel score state using sparse mesh operations.
#
# This file is independent of the production score path.  It implements one
# fixed-kappa SPDE component whose fitted null covariance obeys g' u = 0, while
# the tested kernel is the raw covariance A Q^-1 A' / tau.  It never constructs
# an observation-space P, sqrt(P), covariance matrix, projection basis Z, or a
# dense observation-by-mesh score factor.

.inlast_raw_sparse_matrix <- function(x, name) {
  x <- methods::as(x, "CsparseMatrix")
  if (any(!is.finite(x@x))) stop(name, " must be finite.")
  x
}

.inlast_raw_dense_matrix <- function(x, name) {
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  if (length(dim(x)) != 2L || any(!is.finite(x))) {
    stop(name, " must be a finite numeric matrix.")
  }
  x
}

# Return B -> H_c^-1 B, where H_c^-1 is H^-1 restricted to g' u = 0.
.inlast_raw_constrained_solver <- function(H, g) {
  factor <- Matrix::Cholesky(
    Matrix::forceSymmetric(H), perm = FALSE, LDL = FALSE, super = FALSE
  )
  solve_H <- function(B) {
    B <- if (is.null(dim(B))) matrix(as.numeric(B), ncol = 1L) else B
    as.matrix(Matrix::solve(factor, B))
  }
  Hinv_g <- as.numeric(solve_H(g))
  denominator <- sum(g * Hinv_g)
  if (!is.finite(denominator) || denominator <= 0) {
    stop("The constraint has a non-positive H-inverse norm.")
  }
  function(B) {
    was_vector <- is.null(dim(B))
    answer <- solve_H(B)
    multiplier <- as.numeric(crossprod(g, answer)) / denominator
    answer <- answer - tcrossprod(Hinv_g, multiplier)
    if (was_vector) as.numeric(answer) else answer
  }
}

# Construct one feature's raw-kernel score state.
#
# The null covariance is
#   Vc = D + A (Q^-1 - Q^-1 g (g'Q^-1g)^-1 g'Q^-1) A' / tau,
# and P is its exact nuisance-GLS residual precision.  The tested kernel is
# Graw = A Q^-1 A' / tau.  `stabilize_score = TRUE` represents the score basis
# as A0 = A - X (X'X)^-1 X'A.  Since P X = 0, this leaves a and M unchanged in
# exact arithmetic while avoiding cancellation in the raw near-intercept
# direction.  A0 is handled as a sparse matrix plus a rank-p correction; it is
# never materialized.
inlast_raw_kernel_sparse_state <- function(
    raw_A, raw_Q, g, working_error, working_variance, tau,
    nuisance_X = NULL, stabilize_score = TRUE, tolerance = 1e-9) {
  if (!requireNamespace("Matrix", quietly = TRUE)) stop("Matrix is required.")
  A <- .inlast_raw_sparse_matrix(raw_A, "raw_A")
  Q <- Matrix::forceSymmetric(.inlast_raw_sparse_matrix(raw_Q, "raw_Q"))
  n <- nrow(A)
  m <- ncol(A)
  if (!all(dim(Q) == c(m, m)) ||
      !isTRUE(Matrix::isSymmetric(Q, tol = tolerance))) {
    stop("raw_Q must be symmetric with ncol(raw_A) rows and columns.")
  }
  g <- as.numeric(g)
  expected_g <- as.numeric(Matrix::crossprod(A, rep.int(1 / n, n)))
  if (length(g) != m || any(!is.finite(g)) || !isTRUE(all.equal(
      g, expected_g, tolerance = tolerance, check.attributes = FALSE))) {
    stop("g must equal crossprod(raw_A, 1 / n).")
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
    stop("tau must be one positive finite precision.")
  }
  X <- if (is.null(nuisance_X)) {
    matrix(1, n, 1L)
  } else {
    .inlast_raw_dense_matrix(nuisance_X, "nuisance_X")
  }
  if (nrow(X) != n) stop("nuisance_X has the wrong observation dimension.")
  if (ncol(X) && qr(X, tol = tolerance)$rank < ncol(X)) {
    stop("nuisance_X must have full column rank.")
  }
  stabilize_score <- as.logical(stabilize_score)
  if (length(stabilize_score) != 1L || is.na(stabilize_score)) {
    stop("stabilize_score must be TRUE or FALSE.")
  }

  Dinv <- 1 / D
  Dinv_half_A <- Matrix::Diagonal(x = sqrt(Dinv)) %*% A
  S <- Matrix::forceSymmetric(Matrix::crossprod(Dinv_half_A))
  T <- if (ncol(X)) {
    as.matrix(Matrix::crossprod(A, Dinv * X))
  } else {
    matrix(numeric(), m, 0L)
  }
  U <- if (ncol(X)) crossprod(X, Dinv * X) else
    matrix(numeric(), 0L, 0L)
  H <- Matrix::forceSymmetric(tau * Q + S)
  Hc_solve <- .inlast_raw_constrained_solver(H, g)

  HcS <- Hc_solve(S)
  AWA <- as.matrix(S) - as.matrix(S %*% HcS)
  AWA <- (AWA + t(AWA)) / 2
  if (ncol(X)) {
    HcT <- Hc_solve(T)
    AWX <- T - as.matrix(S %*% HcT)
    XWX <- U - crossprod(T, HcT)
    XWX <- (XWX + t(XWX)) / 2
    nuisance_Vp <- solve(XWX, diag(ncol(X)))
  } else {
    HcT <- matrix(numeric(), m, 0L)
    AWX <- matrix(numeric(), m, 0L)
    XWX <- matrix(numeric(), 0L, 0L)
    nuisance_Vp <- matrix(numeric(), 0L, 0L)
  }

  rhs_e <- as.numeric(Matrix::crossprod(A, Dinv * e))
  Hce <- Hc_solve(rhs_e)
  AWe <- rhs_e - as.numeric(S %*% Hce)
  XWe <- if (ncol(X)) {
    as.numeric(crossprod(X, Dinv * e) - crossprod(T, Hce))
  } else numeric()

  # Low-rank representation A0 = A - X E.  This branch is an exact change of
  # score-factor representative because the P below annihilates X.
  E <- if (stabilize_score && ncol(X)) {
    solve(crossprod(X), as.matrix(Matrix::crossprod(X, A)))
  } else {
    matrix(0, ncol(X), m)
  }
  BWB <- AWA
  BWX <- AWX
  BWe <- AWe
  if (ncol(X) && stabilize_score) {
    BWB <- BWB - AWX %*% E - crossprod(E, t(AWX)) +
      crossprod(E, XWX %*% E)
    BWX <- AWX - crossprod(E, XWX)
    BWe <- AWe - as.numeric(crossprod(E, XWe))
  }
  if (ncol(X)) {
    K <- BWB - BWX %*% nuisance_Vp %*% t(BWX)
    b <- BWe - as.numeric(BWX %*% nuisance_Vp %*% XWe)
  } else {
    K <- BWB
    b <- BWe
  }
  K <- (K + t(K)) / 2

  # Q = L L'.  The raw score factor is A0 L^-T / sqrt(tau), hence
  # a = L^-1 A0'Pe/sqrt(tau) and M = L^-1(A0'PA0)L^-T/tau.
  Q_factor <- Matrix::Cholesky(
    Q, perm = FALSE, LDL = FALSE, super = FALSE
  )
  L <- Matrix::expand(Q_factor)$L
  a <- as.numeric(Matrix::solve(L, b)) / sqrt(tau)
  left <- as.matrix(Matrix::solve(L, K))
  M <- as.matrix(Matrix::solve(L, t(left))) / tau
  M <- (M + t(M)) / 2

  # Quantify how close the raw constraint direction A Q^-1 g is to the
  # nuisance span without constructing the raw observation-by-mesh factor.
  qinv_g <- as.numeric(Matrix::solve(Q_factor, g))
  raw_constraint_direction <- as.numeric(A %*% qinv_g)
  raw_direction_residual <- if (ncol(X)) {
    qr.resid(qr(X), raw_constraint_direction)
  } else raw_constraint_direction
  raw_direction_norm <- sqrt(sum(raw_constraint_direction^2))

  nuisance_error <- if (ncol(X)) {
    max(abs(XWX - XWX %*% nuisance_Vp %*% XWX))
  } else 0
  list(
    a = a,
    M = M,
    coefficient_score = b,
    coefficient_information = K,
    nuisance_covariance = nuisance_Vp,
    diagnostics = list(
      nuisance_normal_equation_error = nuisance_error,
      constraint_solve_error = max(
        abs(crossprod(g, HcS)),
        if (ncol(X)) max(abs(crossprod(g, HcT))) else 0,
        abs(sum(g * Hce))
      ),
      minimum_M_eigenvalue = min(eigen(M, symmetric = TRUE,
                                       only.values = TRUE)$values),
      raw_constraint_direction_nuisance_residual_ratio =
        sqrt(sum(raw_direction_residual^2)) /
        max(raw_direction_norm, .Machine$double.eps),
      score_basis_residualized = stabilize_score,
      dense_mesh_dimension = m,
      observation_square_matrix_constructed = FALSE,
      nuisance_dimension = ncol(X)
    )
  )
}

# Pair two aligned states through the package calibration routine.  `moments`
# contains trace((M1 M2)^k), k=1,...,4; no observation-space square root is
# involved.
inlast_raw_kernel_sparse_pair <- function(state1, state2,
                                          method = c("liu", "davies")) {
  method <- match.arg(method)
  if (length(state1$a) != length(state2$a) ||
      !all(dim(state1$M) == dim(state2$M))) {
    stop("The raw score states are not aligned.")
  }
  signed_score <- as.numeric(crossprod(state1$a, state2$a))
  calibration <- rkhs_score_calibrate(
    signed_score, state1$M, state2$M, method = method
  )
  c(list(signed_score = signed_score, statistic = signed_score^2), calibration)
}
