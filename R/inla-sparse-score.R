# Sparse score backend for one fixed-kappa INLA SPDE target.  The fitted null
# covariance retains the exact observation mean constraint and the nuisance
# adjustment uses INLA's own conditional posterior Vp block.

.inlast_sparse_score_capability <- function(model) {
  spec <- model$inla_spec
  random <- spec$random
  target <- if (is.list(random)) {
    which(vapply(random, function(x) isTRUE(x$target), logical(1L)))
  } else integer()
  reasons <- character()
  if (!is.list(random) || length(random) != 1L) {
    reasons <- c(reasons, "exactly one random block is required")
  }
  if (length(target) != 1L || !identical(target, 1L)) {
    reasons <- c(reasons, "exactly one target SPDE block is required")
  } else if (!identical(random[[target]]$kind, "spde")) {
    reasons <- c(reasons, "the target random block must be an SPDE")
  }
  if (length(model$geometry$target) != 1L ||
      !identical(names(model$geometry$target), "global")) {
    reasons <- c(reasons, "only the single global target is supported")
  }
  list(eligible = !length(reasons), reason = paste(unique(reasons), collapse = "; "))
}

.inlast_sparse_score_geometry <- function(model) {
  capability <- .inlast_sparse_score_capability(model)
  if (!capability$eligible) {
    stop("The sparse score backend is unavailable: ", capability$reason, ".")
  }
  block <- model$inla_spec$random[[1L]]
  A <- methods::as(block$A, "CsparseMatrix")
  Q <- Matrix::forceSymmetric(methods::as(block$Q, "CsparseMatrix"))
  g <- as.numeric(block$constraint)
  Z <- as.matrix(block$projection)
  Qprojected <- crossprod(Z, as.matrix(Q %*% Z))
  Qprojected <- (Qprojected + t(Qprojected)) / 2
  R <- chol(Qprojected)
  Tbase <- Z %*% backsolve(R, diag(ncol(Z)))
  list(
    A = A, Q = Q, constraint = g, projection = Z,
    coefficient_factor = Tbase,
    target = "global", sp_index = as.integer(block$sp_index),
    definition = paste(
      "single fixed-kappa SPDE conditioned on the observation mean;",
      "native INLA nuisance Vp"
    )
  )
}

.mgcvst_model_sparse_constrained_solver <- function(H, g) {
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
    stop("The sparse score constraint has a non-positive H-inverse norm.")
  }
  function(rhs) {
    vector <- is.null(dim(rhs))
    answer <- solve_H(rhs)
    multiplier <- as.numeric(crossprod(g, answer)) / denominator
    answer <- answer - tcrossprod(Hinv_g, multiplier)
    if (vector) as.numeric(answer) else answer
  }
}

.mgcvst_model_sparse_score_state <- function(fit, feature) {
  geometry <- fit$score_sparse
  if (!is.list(geometry) || is.null(geometry$coefficient_factor)) {
    stop("The fit lacks its sparse INLA score geometry.")
  }
  phi <- as.numeric(fit$dispersion[feature])
  sp <- as.numeric(fit$smoothing_parameters[feature, geometry$sp_index])
  if (length(phi) != 1L || length(sp) != 1L ||
      !is.finite(phi) || phi <= 0 || !is.finite(sp) || sp <= 0) {
    stop("The feature has invalid dispersion or smoothing parameters.")
  }
  tau <- sp / phi
  A <- geometry$A
  Q <- geometry$Q
  g <- geometry$constraint
  T <- geometry$coefficient_factor / sqrt(tau)
  e <- as.numeric(fit$working_error[, feature])
  D <- as.numeric(fit$working_variance[, feature])
  if (length(e) != nrow(A) || length(D) != nrow(A) ||
      any(!is.finite(e)) || any(!is.finite(D)) || any(D <= 0)) {
    stop("The feature has an invalid working state.")
  }
  X <- as.matrix(fit$geometry$nuisance_design)
  Vp <- fit$nuisance_covariance[[feature]]
  if (!is.matrix(Vp) || !all(dim(Vp) == ncol(X)) ||
      any(!is.finite(Vp))) {
    stop("The feature has an incompatible conditional nuisance covariance.")
  }
  Vp <- (Vp + t(Vp)) / 2

  W <- 1 / D
  K <- Matrix::forceSymmetric(Matrix::crossprod(A, A * W))
  L <- if (ncol(X)) as.matrix(Matrix::crossprod(A, W * X)) else
    matrix(numeric(), ncol(A), 0L)
  tvec <- as.numeric(Matrix::crossprod(A, W * e))
  solve_S <- .mgcvst_model_sparse_constrained_solver(
    Matrix::forceSymmetric(tau * Q + K), g
  )
  solved_small <- solve_S(cbind(L, tvec))
  SL <- if (ncol(X)) solved_small[, seq_len(ncol(X)), drop = FALSE] else
    matrix(numeric(), ncol(A), 0L)
  St <- as.numeric(solved_small[, ncol(X) + 1L])
  U <- if (ncol(X)) L - as.matrix(K %*% SL) else
    matrix(numeric(), ncol(A), 0L)
  q <- if (ncol(X)) {
    as.numeric(crossprod(X, W * e) - crossprod(L, St))
  } else numeric()
  h <- tvec - as.numeric(K %*% St)
  if (ncol(X)) h <- h - as.numeric(U %*% Vp %*% q)

  KT <- as.matrix(K %*% T)
  SKT <- solve_S(KT)
  a <- as.numeric(crossprod(T, h))
  M <- crossprod(T, KT) - crossprod(KT, SKT)
  if (ncol(X)) {
    UtT <- crossprod(U, T)
    M <- M - crossprod(UtT, Vp %*% UtT)
  }
  M <- (M + t(M)) / 2
  list(
    a = a, M = M,
    width = stats::setNames(ncol(T), geometry$target),
    target = NULL, operator = NULL, backend = "sparse_conditioned_INLA"
  )
}
