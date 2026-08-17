# Extract an fmesher mesh and its coordinate transformation.
.spde_mesh_parts <- function(mesh) {
  if (inherits(mesh, "spde_mesh")) {
    return(list(mesh = mesh$mesh, transform = mesh$transform))
  }
  if (inherits(mesh, "fm_mesh_2d")) {
    return(list(mesh = mesh, transform = list(center = c(0, 0), scale = 1)))
  }
  stop("mesh must be a spde_mesh or fm_mesh_2d object.")
}

#' Construct the INLA MatÃ©rn SPDE matrices used by mgcvST
#'
#' @param mesh A `spde_mesh` or `fm_mesh_2d` object.
#' @param alpha SPDE operator order. The current mgcv smooth requires `2`.
#' @param constr Passed to [INLA::inla.spde2.matern()]. Use `FALSE` for the
#'   proper no-null-space smooth.
#' @param ... Additional arguments passed to [INLA::inla.spde2.matern()].
#' @return An INLA SPDE model object.
#' @export
spde_model <- function(mesh, alpha = 2, constr = FALSE, ...) {
  if (!identical(as.numeric(alpha), 2)) {
    stop("The current mgcvST smooth requires alpha = 2.")
  }
  x <- .spde_mesh_parts(mesh)
  INLA::inla.spde2.matern(x$mesh, alpha = alpha, constr = constr, ...)
}

#' Construct a sparse observation-to-mesh projector
#'
#' @param mesh A `spde_mesh` or `fm_mesh_2d` object.
#' @param loc Observation or prediction coordinates.
#' @return A sparse projector matrix with one row per location.
#' @export
spde_project <- function(mesh, loc) {
  x <- .spde_mesh_parts(mesh)
  loc <- .spde_xy(loc)
  loc <- sweep(loc, 2, x$transform$center, "-") / x$transform$scale
  INLA::inla.spde.make.A(x$mesh, loc = loc)
}

#' Construct finite element matrices for an SPDE mesh
#'
#' @param mesh A `spde_mesh` or `fm_mesh_2d` object.
#' @param order Finite element integration order passed to [fmesher::fm_fem()].
#' @return The finite element matrix list returned by `fmesher`.
#' @export
spde_fem <- function(mesh, order = 2L) {
  x <- .spde_mesh_parts(mesh)
  fmesher::fm_fem(x$mesh, order = order)
}

#' Evaluate the proper MatÃ©rn SPDE precision matrix
#'
#' @param model An object returned by `spde_model()`.
#' @param kappa Positive spatial scale parameter.
#' @param tau Positive precision multiplier.
#' @return A sparse symmetric precision matrix.
#' @export
spde_precision <- function(model, kappa, tau = 1) {
  if (length(kappa) != 1L || !is.finite(kappa) || kappa <= 0) {
    stop("kappa must be one positive finite number.")
  }
  if (length(tau) != 1L || !is.finite(tau) || tau <= 0) {
    stop("tau must be one positive finite number.")
  }
  p <- model$param.inla
  if (is.null(p$M0) || is.null(p$M1) || is.null(p$M2)) {
    stop("model does not contain the required M0, M1, and M2 matrices.")
  }
  Matrix::forceSymmetric(
    tau^2 * (kappa^4 * p$M0 + 2 * kappa^2 * p$M1 + p$M2)
  )
}

#' mgcv SPDE smooth methods
#'
#' These functions are the S3 constructor and prediction-matrix methods used
#' by `mgcv` for `bs = "spde"`. They are registered as S3 methods when
#' `mgcvST` is loaded.
#'
#' @param object An `mgcv` SPDE smooth specification or fitted SPDE smooth.
#' @param data A data-like object containing the smooth coordinates.
#' @param knots Knot information supplied by `mgcv`.
#' @return `smooth.construct.spde.smooth.spec()` returns a fitted smooth
#'   specification. `Predict.matrix.spde.smooth()` returns the numeric
#'   prediction matrix at `data`.
#' @name spde_smooth_methods
#' @method smooth.construct spde.smooth.spec
#' @export
smooth.construct.spde.smooth.spec <- function(object, data, knots) {
  if (is.null(object$sp)) {
    object$sp <- c(-1, 0.1)
  }
  if (!is.numeric(object$sp) || length(object$sp) != 2L ||
      any(!is.finite(object$sp)) ||
      (object$sp[1L] != -1 && object$sp[1L] <= 0) || object$sp[2L] <= 0) {
    stop(
      "bs = 'spde' requires sp = c(lambda, kappa), with lambda equal to -1 or positive and kappa positive and fixed."
    )
  }
  if (is.null(object$xt$mesh)) {
    stop("bs = 'spde' requires xt = list(mesh = ...).")
  }

  loc <- cbind(data[[object$term[1]]], data[[object$term[2]]])
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
    G <- as.matrix(Matrix::crossprod(A, matrix(1, nrow(A), 1))) / nrow(A)
    fitqr <- qr(G)
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
  object$X <- as.matrix(A %*% projection)
  m <- ncol(object$X)

  kappa <- object$sp[2L]
  Q <- kappa^4 * p$M0 + 2 * kappa^2 * p$M1 + p$M2
  Q <- crossprod(projection, Q %*% projection)
  object$S <- list(as.matrix(Matrix::forceSymmetric(Q)))
  object$sp <- object$sp[1L]
  object$L <- NULL
  object$kappa <- kappa
  object$kappa.estimated <- FALSE
  object$rank <- as.integer(vapply(object$S, Matrix::rankMatrix, numeric(1)))
  object$null.space.dim <- 0L
  object$C <- matrix(0, 0, m)
  object$side.constrain <- FALSE
  object$no.rescale <- TRUE
  object$df <- object$bs.dim <- m
  object$mesh <- x$mesh
  object$transform <- x$transform
  object$spde.model <- model
  object$projection <- projection
  object$projection.rank <- projection.rank
  object$project.intercept <- project.intercept
  object$raw.dimension <- m.raw
  class(object) <- c("spde.smooth", "mgcv.smooth")
  object
}

#' @rdname spde_smooth_methods
#' @method Predict.matrix spde.smooth
#' @export
Predict.matrix.spde.smooth <- function(object, data) {
  loc <- cbind(data[[object$term[1]]], data[[object$term[2]]])
  loc <- sweep(loc, 2, object$transform$center, "-") /
    object$transform$scale
  A <- INLA::inla.spde.make.A(object$mesh, loc = loc)
  projection <- object[["projection"]]
  if (is.null(projection)) {
    projection <- diag(ncol(A))
  }
  if (!is.matrix(projection) || nrow(projection) != ncol(A)) {
    stop("The stored SPDE projection is incompatible with the mesh.")
  }
  as.matrix(A %*% projection)
}
