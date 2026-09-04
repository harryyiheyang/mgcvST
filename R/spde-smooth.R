#' Construct the finite-element matrices used by mgcvST
#'
#' This helper belongs to the basis-construction path. Fitting a prepared
#' basis does not call it.
#'
#' @param mesh A [spde_mesh()], `fm_mesh_2d`, or compatible mesh list.
#' @param order Finite-element order. Only linear triangles (`2`) are used.
#' @return A list containing `M0`, `M1`, and `M2`.
#' @export
spde_fem <- function(mesh, order = 2L) {
  if (!identical(as.integer(order), 2L)) {
    stop("The current mgcvST basis requires order = 2.")
  }
  .spde_basis_fem(.spde_basis_mesh(mesh))
}

#' Construct a sparse observation-to-mesh projector
#'
#' This helper belongs to the basis-construction path. Fitting a prepared
#' basis does not call it.
#'
#' @param mesh A [spde_mesh()], `fm_mesh_2d`, or compatible mesh list.
#' @param loc Coordinates in the original coordinate system.
#' @return A sparse projector matrix with one row per location.
#' @export
spde_project <- function(mesh, loc) {
  x <- .spde_basis_mesh(mesh)
  z <- .spde_xy(loc)
  z <- sweep(z, 2L, x$transform$center, "-") / x$transform$scale
  .spde_basis_project(x, z)
}

#' Construct the SPDE penalty model
#'
#' @param mesh A [spde_mesh()], `fm_mesh_2d`, or compatible mesh list.
#' @param alpha SPDE operator order. Only `2` is supported.
#' @param constr Retained for compatibility; must be `FALSE`.
#' @param ... Unused.
#' @return A lightweight model containing `M0`, `M1`, and `M2`.
#' @export
spde_model <- function(mesh, alpha = 2, constr = FALSE, ...) {
  if (!identical(as.numeric(alpha), 2)) {
    stop("The current mgcvST basis requires alpha = 2.")
  }
  if (!identical(constr, FALSE)) {
    stop("The current mgcvST basis requires constr = FALSE.")
  }
  list(param.inla = spde_fem(mesh))
}

#' Evaluate the proper Matern SPDE precision matrix
#'
#' @param model An object returned by [spde_model()].
#' @param kappa Positive dimensionless spatial scale.
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
    stop("model does not contain M0, M1, and M2.")
  }
  Matrix::forceSymmetric(
    tau^2 * (kappa^4 * p$M0 + 2 * kappa^2 * p$M1 + p$M2)
  )
}

#' mgcv SPDE smooth methods
#'
#' Both smooths consume a complete object returned by [spde_basis()] through
#' `xt`. The mgcv argument `k` is ignored because basis dimension is fixed at
#' construction time.
#'
#' @param object An `mgcv` SPDE smooth specification or fitted smooth.
#' @param data A data-like object containing the smooth coordinates.
#' @param knots Knot information supplied by `mgcv`.
#' @return The constructed smooth or its prediction matrix.
#' @name spde_smooth_methods
#' @method smooth.construct spde.smooth.spec
#' @export
smooth.construct.spde.smooth.spec <- function(object, data, knots) {
  .spde_basis_warn_k(object, "spde")
  loc <- cbind(data[[object$term[1L]]], data[[object$term[2L]]])
  .spde_basis_validate(object$xt, loc)
  basis <- object$xt
  object$timing <- list2env(list(basis_seconds = 0, basis_calls = 0L,
                                prediction_seconds = 0, prediction_calls = 0L),
                           parent = emptyenv())
  object$X <- .spde_basis_at(basis, loc, timing = object$timing, stage = "basis")

  if (is.null(basis$kappa)) {
    object$S <- basis$penalty
    object$sp <- rep(-1, 2L)
    object$L <- rbind(c(1, 4), c(1, 2), c(1, 0))
    object$kappa.estimated <- TRUE
  } else {
    object$S <- list(basis$Q)
    object$sp <- -1
    object$L <- NULL
    object$kappa.estimated <- FALSE
  }
  object$kappa <- basis$kappa
  object$rank <- as.integer(vapply(object$S, Matrix::rankMatrix, numeric(1)))
  object$null.space.dim <- 0L
  object$C <- matrix(0, 0L, ncol(object$X))
  object$side.constrain <- FALSE
  object$no.rescale <- TRUE
  object$df <- object$bs.dim <- ncol(object$X)
  object$basis <- basis
  object$projection.rank <- basis$projection_rank
  object$project.intercept <- basis$project_intercept
  object$raw.dimension <- basis$raw_dimension
  z <- .spde_basis_component(object)
  object$component <- z$component
  object$score.component <- z$score.component
  class(object) <- c("spde.smooth", "mgcv.smooth")
  object
}

#' @rdname spde_smooth_methods
#' @method Predict.matrix spde.smooth
#' @export
Predict.matrix.spde.smooth <- function(object, data) {
  loc <- cbind(data[[object$term[1L]]], data[[object$term[2L]]])
  .spde_basis_at(object$basis, loc, timing = object$timing)
}
