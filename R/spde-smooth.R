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

  object$S <- list(basis$Q)
  object$sp <- -1
  object$L <- NULL
  kappa <- .spde_kappa_fields(basis)
  object[names(kappa)] <- kappa
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
