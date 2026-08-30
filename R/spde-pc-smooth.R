#' Principal-component SPDE smooth methods
#'
#' The PC dimension is determined only from the cumulative covariance
#' contribution stored in `xt$pc_cutoff`. The mgcv argument `k` is ignored.
#'
#' @param object An `mgcv` principal-component SPDE smooth specification or
#'   fitted smooth.
#' @param data A data-like object containing the smooth coordinates.
#' @param knots Knot information supplied by `mgcv`.
#' @return The constructed smooth or its prediction matrix.
#' @name spdePC_smooth_methods
#' @method smooth.construct spdePC.smooth.spec
#' @export
smooth.construct.spdePC.smooth.spec <- function(object, data, knots) {
  .spde_basis_warn_k(object, "spdePC")
  loc <- cbind(data[[object$term[1L]]], data[[object$term[2L]]])
  .spde_basis_validate(object$xt, loc, pc = TRUE)
  basis <- object$xt
  retained <- which(basis$pc_cumulative >= basis$pc_cutoff)[1L]
  V <- basis$pc_vectors[, seq_len(retained), drop = FALSE]
  evals <- basis$pc_values[seq_len(retained)]

  object$X <- CppMatrix::matrixMultiply(basis$B, V)
  object$S <- list(diag(1 / evals))
  object$sp <- -1
  object$L <- NULL
  object$rank <- retained
  object$null.space.dim <- 0L
  object$C <- matrix(0, 0L, retained)
  object$side.constrain <- FALSE
  object$no.rescale <- TRUE
  object$df <- object$bs.dim <- retained
  object$basis <- basis
  object$score_basis <- basis$B
  object$pc_projection <- V
  object$pc_eigenvectors <- V
  object$pc_score_Q <- (V * rep(1 / evals, each = nrow(V))) %*% t(V)
  object$pc_score_Q <- (object$pc_score_Q + t(object$pc_score_Q)) / 2
  object$pc_eigenvalues <- basis$pc_values
  object$pc_cumulative_trace <- basis$pc_cumulative
  object$pc_cutoff <- basis$pc_cutoff
  object$pc_full_dimension <- ncol(basis$B)
  object$pc_retained_dimension <- retained
  object$projection.rank <- basis$projection_rank
  object$project.intercept <- basis$project_intercept
  object$raw.dimension <- basis$raw_dimension
  object$kappa <- basis$kappa
  object$kappa.estimated <- FALSE
  z <- .spde_basis_component(object)
  object$component <- z$component
  object$score.component <- z$score.component
  class(object) <- c("spdePC.smooth", "mgcv.smooth")
  object
}

#' @rdname spdePC_smooth_methods
#' @method Predict.matrix spdePC.smooth
#' @export
Predict.matrix.spdePC.smooth <- function(object, data) {
  loc <- cbind(data[[object$term[1L]]], data[[object$term[2L]]])
  .spde_basis_validate(object$basis, loc, pc = TRUE)
  CppMatrix::matrixMultiply(object$basis$B, object$pc_projection)
}
