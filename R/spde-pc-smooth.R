#' Public principal-component SPDE smooth methods
#'
#' These functions are the S3 constructor and prediction-matrix methods used
#' by `mgcv` for `bs = "spdePC"`. The constructor projects the intercept in
#' coefficient space, eigendecomposes the resulting fixed-kappa covariance
#' `G = Q^{-1}`, and retains the leading directions whose eigenvalues explain
#' `xt$pc_cutoff` of `tr(G)`. The fitted object retains `pc_score_Q`, the
#' corresponding positive-semidefinite precision in the full projected
#' coefficient coordinates, for covariance score tests.
#'
#' The spatial scale is intentionally fixed: supply `sp = c(lambda, kappa)`
#' with a positive `kappa` (normally `sp = c(-1, 0.1)`). The first entry can
#' be `-1` to estimate only the smoothing parameter `lambda` by REML.
#'
#' @param object An `mgcv` principal-component SPDE smooth specification or
#'   fitted smooth.
#' @param data A data-like object containing the smooth coordinates.
#' @param knots Knot information supplied by `mgcv`.
#' @return `smooth.construct.spdePC.smooth.spec()` returns a fitted smooth
#'   specification. `Predict.matrix.spdePC.smooth()` returns the numeric
#'   prediction matrix at `data`.
#' @name spdePC_smooth_methods
#' @method smooth.construct spdePC.smooth.spec
#' @export
smooth.construct.spdePC.smooth.spec <- function(object, data, knots) {
  if (is.null(object$sp)) object$sp <- c(-1, 0.1)
  if (!is.numeric(object$sp) || length(object$sp) != 2L ||
      any(!is.finite(object$sp)) ||
      !(object$sp[1L] == -1 || object$sp[1L] > 0) ||
      object$sp[2L] <= 0) {
    stop(
      "bs = 'spdePC' requires sp = c(lambda, kappa), with lambda equal to -1 or positive and kappa positive and fixed."
    )
  }
  if (is.null(object$xt$mesh)) {
    stop("bs = 'spdePC' requires xt = list(mesh = ..., pc_cutoff = ...).")
  }
  pc.cutoff <- object$xt$pc_cutoff
  if (is.null(pc.cutoff)) pc.cutoff <- 0.999
  pc.cutoff <- as.numeric(pc.cutoff)
  if (length(pc.cutoff) != 1L || !is.finite(pc.cutoff) ||
      pc.cutoff <= 0 || pc.cutoff > 1) {
    stop("xt$pc_cutoff must be one finite number in (0, 1].")
  }

  loc <- cbind(data[[object$term[1L]]], data[[object$term[2L]]])
  x <- .spde_mesh_parts(object$xt$mesh)
  loc <- sweep(loc, 2, x$transform$center, "-") / x$transform$scale
  model <- object$xt$model
  if (is.null(model)) {
    model <- INLA::inla.spde2.matern(x$mesh, alpha = 2, constr = FALSE)
  }
  A <- INLA::inla.spde.make.A(x$mesh, loc = loc)
  if (any(Matrix::rowSums(A) == 0)) {
    stop("Some model locations are outside the supplied SPDE mesh.")
  }
  p <- model$param.inla
  m.raw <- ncol(A)
  project.intercept <- object$xt$project.intercept
  if (is.null(project.intercept)) project.intercept <- TRUE
  if (!is.logical(project.intercept) || length(project.intercept) != 1L ||
      is.na(project.intercept)) {
    stop("xt$project.intercept must be TRUE or FALSE.")
  }

  if (project.intercept) {
    G0 <- as.matrix(Matrix::crossprod(A, matrix(1, nrow(A), 1L))) / nrow(A)
    fitqr <- qr(G0)
    projection.rank <- fitqr$rank
    if (projection.rank < 1L || projection.rank >= m.raw) {
      stop("The intercept projection does not have a valid coefficient-space rank.")
    }
    projection <- qr.Q(fitqr, complete = TRUE)[,
      (projection.rank + 1L):m.raw, drop = FALSE]
  } else {
    projection.rank <- 0L
    projection <- diag(m.raw)
  }

  kappa <- object$sp[2L]
  Q <- kappa^4 * p$M0 + 2 * kappa^2 * p$M1 + p$M2
  Q <- crossprod(projection, Q %*% projection)
  Q <- as.matrix(Matrix::forceSymmetric(Q))
  m <- ncol(Q)
  G <- CppMatrix::matrixSolve(Q, diag(m))
  G <- (G + t(G)) / 2
  E <- CppMatrix::matrixEigen(G)
  evals <- as.numeric(E$values)
  vectors <- as.matrix(E$vectors)
  ord <- order(evals, decreasing = TRUE)
  evals <- evals[ord]
  vectors <- vectors[, ord, drop = FALSE]
  tol <- sqrt(.Machine$double.eps) * max(1, max(abs(evals)))
  if (min(evals) < -tol) {
    stop("The fixed-kappa projected SPDE covariance is not positive semidefinite.")
  }
  evals <- pmax(evals, 0)
  if (!is.finite(sum(evals)) || sum(evals) <= 0 || any(evals <= 0)) {
    stop("The fixed-kappa projected SPDE covariance has invalid eigenvalues.")
  }
  trace.fraction <- cumsum(evals) / sum(evals)
  retained <- which(trace.fraction >= pc.cutoff)[1L]
  pc.vectors <- vectors[, seq_len(retained), drop = FALSE]
  pc.projection <- CppMatrix::matrixMultiply(projection, pc.vectors)
  score.basis <- CppMatrix::matrixMultiply(as.matrix(A), projection)
  pc.score.Q <- CppMatrix::matrixMultiply(
    sweep(pc.vectors, 2L, 1 / evals[seq_len(retained)], "*"),
    t(pc.vectors)
  )
  pc.score.Q <- (pc.score.Q + t(pc.score.Q)) / 2
  object$X <- CppMatrix::matrixMultiply(as.matrix(A), pc.projection)
  object$S <- list(diag(1 / evals[seq_len(retained)]))
  object$sp <- object$sp[1L]
  object$L <- NULL
  object$rank <- as.integer(Matrix::rankMatrix(object$S[[1L]]))
  object$null.space.dim <- 0L
  object$C <- matrix(0, 0L, retained)
  object$side.constrain <- FALSE
  object$no.rescale <- TRUE
  object$df <- object$bs.dim <- retained
  object$mesh <- x$mesh
  object$transform <- x$transform
  object$spde.model <- model
  object$projection <- projection
  object$projection.rank <- projection.rank
  object$score_basis <- score.basis
  object$pc_projection <- pc.projection
  object$pc_eigenvectors <- pc.vectors
  object$pc_score_Q <- pc.score.Q
  object$pc_eigenvalues <- evals
  object$pc_cumulative_trace <- trace.fraction
  object$pc_cutoff <- pc.cutoff
  object$pc_full_dimension <- m
  object$pc_retained_dimension <- retained
  object$project.intercept <- project.intercept
  object$raw.dimension <- m.raw
  object$kappa <- kappa
  object$kappa.estimated <- FALSE
  class(object) <- c("spdePC.smooth", "mgcv.smooth")
  object
}

#' @rdname spdePC_smooth_methods
#' @method Predict.matrix spdePC.smooth
#' @export
Predict.matrix.spdePC.smooth <- function(object, data) {
  loc <- cbind(data[[object$term[1L]]], data[[object$term[2L]]])
  loc <- sweep(loc, 2, object$transform$center, "-") /
    object$transform$scale
  A <- INLA::inla.spde.make.A(object$mesh, loc = loc)
  pc.projection <- object[["pc_projection"]]
  if (!is.matrix(pc.projection) || nrow(pc.projection) != ncol(A)) {
    stop("The stored SPDE principal-component projection is incompatible with the mesh.")
  }
  CppMatrix::matrixMultiply(as.matrix(A), pc.projection)
}
