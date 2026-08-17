# Multiply dense matrices through the package's single-thread numerical backend.
.magic_mm <- function(A, B, transA = FALSE, transB = FALSE) {
  CppMatrix::matrixMultiply(
    as.matrix(A), as.matrix(B), transA = transA, transB = transB
  )
}

# Solve a dense linear system through the package numerical backend.
.magic_solve <- function(A, B) {
  CppMatrix::matrixSolve(as.matrix(A), as.matrix(B))
}

# Normalize and validate a finite numeric matrix input.
.as_numeric_matrix <- function(x, name) {
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  if (length(dim(x)) != 2L || any(!is.finite(x))) {
    stop(name, " must be a finite numeric matrix.")
  }
  x
}

#' Construct a compact RKHS score operator
#'
#' Constructs the covariance factor `T = sqrt(field_scale) B Q^{-1/2}` and
#' the compact Woodbury quantities needed to apply the residual precision
#' operator
#' \deqn{P=V^{-1}-V^{-1}X(X^T V^{-1}X)^{-1}X^T V^{-1},}
#' where `V = diag(diag_var) + T T^T`. The function never constructs `V` or
#' `P` as an observation-by-observation matrix.
#'
#' @param B Observation-by-SPDE-coefficient basis matrix.
#' @param Q Symmetric SPDE precision matrix. All features in one RKHS analysis
#'   must use the same fixed `Q` and coefficient ordering.
#' @param diag_var Positive diagonal working variances.
#' @param X Fixed-effect design. `NULL` gives an intercept. Supply an `n`-row,
#'   zero-column matrix for no fixed-effect projection.
#' @param field_scale Positive marginal field covariance scale.
#' @param allow_psd Whether `Q` may be positive semidefinite. Its Moore-Penrose
#'   inverse then defines the retained field covariance, retaining zero-variance
#'   directions in the score coordinates.
#' @return An object of class `rkhs_score_operator`.
#' @export
rkhs_score_operator <- function(B, Q, diag_var, X = NULL,
                                field_scale = 1, allow_psd = FALSE) {
  B <- .as_numeric_matrix(B, "B")
  n <- nrow(B)
  q <- ncol(B)
  if (q < 1L) stop("B must have at least one column.")

  Q <- Matrix::Matrix(Q, sparse = TRUE)
  if (!all(dim(Q) == c(q, q))) {
    stop("Q must be square with ncol(B) rows and columns.")
  }
  if (!isTRUE(Matrix::isSymmetric(Q, tol = 1e-10))) {
    stop("Q must be symmetric.")
  }
  Q <- Matrix::forceSymmetric(Q)
  allow_psd <- as.logical(allow_psd)
  if (length(allow_psd) != 1L || is.na(allow_psd)) {
    stop("allow_psd must be TRUE or FALSE.")
  }

  diag_var <- as.numeric(diag_var)
  if (length(diag_var) != n || any(!is.finite(diag_var)) ||
      any(diag_var <= 0)) {
    stop("diag_var must contain one positive finite value per row of B.")
  }
  field_scale <- as.numeric(field_scale)
  if (length(field_scale) != 1L || !is.finite(field_scale) ||
      field_scale <= 0) {
    stop("field_scale must be one positive finite value.")
  }

  if (is.null(X)) {
    X <- matrix(1, nrow = n, ncol = 1L)
  } else {
    X <- .as_numeric_matrix(X, "X")
    if (nrow(X) != n) stop("X and B must have the same number of rows.")
  }

  if (allow_psd) {
    E <- CppMatrix::matrixEigen(as.matrix(Q))
    d <- as.numeric(E$values)
    V <- as.matrix(E$vectors)
    tol <- sqrt(.Machine$double.eps) * max(1, max(abs(d)))
    if (min(d) < -tol) stop("Q must be positive semidefinite.")
    keep <- d > tol
    if (!any(keep)) stop("Q must have at least one positive eigenvalue.")
    Qinvhalf <- CppMatrix::matrixMultiply(
      sweep(V[, keep, drop = FALSE], 2L, 1 / sqrt(d[keep]), "*"),
      t(V[, keep, drop = FALSE])
    )
    T <- sqrt(field_scale) * .magic_mm(B, Qinvhalf)
  } else {
    d <- as.numeric(CppMatrix::matrixEigen(as.matrix(Q))$values)
    tol <- sqrt(.Machine$double.eps) * max(1, max(abs(d)))
    if (min(d) <= tol) {
      stop("Q must be positive definite unless allow_psd = TRUE.")
    }
    R <- Matrix::chol(Q)
    Rinv <- Matrix::solve(R, Matrix::Diagonal(q))
    T <- sqrt(field_scale) * .magic_mm(B, as.matrix(Rinv))
  }
  .rkhs_score_operator_factor(
    T, diag_var, X, field_scale = field_scale, B = B, Q = Q
  )
}

# Construct an operator when the aligned field factor is already available.
.rkhs_score_operator_factor <- function(T, diag_var, X,
                                        field_scale = 1,
                                        B = NULL, Q = NULL) {
  n <- nrow(T)
  q <- ncol(T)
  Dinv <- 1 / diag_var
  DinvT <- Dinv * T
  M <- diag(q) + .magic_mm(T, DinvT, transA = TRUE)

  if (ncol(X) > 0L) {
    DinvX <- Dinv * X
    VTX <- .magic_mm(T, DinvX, transA = TRUE)
    VinvX <- DinvX - DinvT %*% .magic_solve(M, VTX)
    XVX <- .magic_mm(X, VinvX, transA = TRUE)
    if (qr(XVX)$rank != ncol(X)) {
      stop("X is rank deficient under the marginal covariance metric.")
    }
    XVXinv <- .magic_solve(XVX, diag(ncol(X)))
  } else {
    VinvX <- matrix(numeric(0), nrow = n, ncol = 0L)
    XVXinv <- matrix(numeric(0), nrow = 0L, ncol = 0L)
  }

  structure(
    list(
      n = n,
      q = q,
      B = B,
      Q = Q,
      field_factor = T,
      diag_var = diag_var,
      Dinv = Dinv,
      DinvT = DinvT,
      woodbury = M,
      X = X,
      VinvX = VinvX,
      XtVinvX_inv = XVXinv,
      field_scale = field_scale
    ),
    class = "rkhs_score_operator"
  )
}

#' Apply a compact residual precision operator
#'
#' @param operator An object returned by [rkhs_score_operator()].
#' @param value A vector or matrix with `operator$n` rows.
#' @return A numeric vector when `value` is a vector, otherwise a matrix.
#' @export
rkhs_score_apply_P <- function(operator, value) {
  if (!inherits(operator, "rkhs_score_operator")) {
    stop("operator must be returned by rkhs_score_operator().")
  }
  was_vector <- is.null(dim(value))
  value <- if (was_vector) matrix(as.numeric(value), ncol = 1L) else
    .as_numeric_matrix(value, "value")
  if (nrow(value) != operator$n || any(!is.finite(value))) {
    stop("value must be finite and have operator$n rows.")
  }

  DinvY <- operator$Dinv * value
  TtDinvY <- .magic_mm(operator$field_factor, DinvY, transA = TRUE)
  VinvY <- DinvY - operator$DinvT %*%
    .magic_solve(operator$woodbury, TtDinvY)

  if (ncol(operator$X) > 0L) {
    XtVinvY <- .magic_mm(operator$X, VinvY, transA = TRUE)
    adj <- operator$VinvX %*%
      .magic_mm(operator$XtVinvX_inv, XtVinvY)
    ans <- VinvY - adj
  } else {
    ans <- VinvY
  }

  if (was_vector) as.numeric(ans) else as.matrix(ans)
}
