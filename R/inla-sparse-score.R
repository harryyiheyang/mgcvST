# Sparse score backend for one fixed-kappa INLA SPDE target.  The fitted null
# covariance retains the exact observation mean constraint and the nuisance
# adjustment uses the expected working curvature.

.inlast_sparse_score_capability <- function(model) {
  spec <- model$inla_spec
  random <- spec$random
  target <- if (is.list(random)) {
    which(vapply(random, function(x) isTRUE(x$target), logical(1L)))
  } else integer()
  reasons <- character()
  if (!is.list(random) || !length(random)) {
    reasons <- c(reasons, "one SPDE target block is required")
  }
  nuisance <- if (is.list(random) && length(random) > 1L) random[-1L] else list()
  if (length(nuisance) &&
      !all(vapply(nuisance, function(z) {
        identical(z$kind, "nuisance") && identical(z$subtype, "iid")
      }, logical(1L)))) {
    reasons <- c(reasons, "every nuisance random block must be iid")
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
  list(
    A = A, Q = Q, constraint = g,
    cache = new.env(parent = emptyenv()),
    target = "global", sp_index = as.integer(block$sp_index),
    normalization = ncol(Q) - 1L,
    definition = "constrained sparse INLA; expected curvature; penalized nuisance Vp"
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

.inlast_sparse_prepare <- function(fit) {
  geometry <- fit$score_sparse
  if (!is.list(geometry) || is.null(geometry$Q)) {
    stop("The fit lacks its sparse INLA score geometry.")
  }
  if (!is.environment(geometry$cache)) geometry$cache <- new.env(parent = emptyenv())
  cache <- geometry$cache
  valid <- identical(cache$Q, geometry$Q) &&
    identical(cache$constraint, geometry$constraint) &&
    mgcvst_inla_sparse_prepared_valid_cpp(cache$prepared)
  if (!valid) {
    Q <- methods::as(methods::as(geometry$Q, "generalMatrix"), "CsparseMatrix")
    cache$prepared <- mgcvst_inla_sparse_prepare_cpp(Q, as.numeric(geometry$constraint))
    cache$Q <- geometry$Q
    cache$constraint <- geometry$constraint
    cache$general_Q <- Q
    g <- geometry$constraint / sqrt(sum(geometry$constraint^2))
    Qg <- as.numeric(Q %*% g)
    cache$penalty_norm <- sqrt(sum(Q@x^2) - 2 * sum(Qg^2) + sum(g * Qg)^2)
  }
  if (!identical(cache$A, geometry$A)) {
    cache$general_A <- methods::as(methods::as(geometry$A, "generalMatrix"), "CsparseMatrix")
    cache$A <- geometry$A
  }
  fit$score_sparse <- geometry
  fit
}

.inlast_sparse_observation_basis <- function(fit, coverage = 0.995,
                                              full_rank = FALSE) {
  fit <- .inlast_sparse_prepare(fit)
  geometry <- fit$score_sparse
  cache <- geometry$cache
  valid <- identical(cache$observation_basis_A, geometry$A) &&
    identical(cache$observation_basis_Q, geometry$Q) &&
    identical(cache$observation_basis_constraint, geometry$constraint) &&
    identical(cache$observation_basis_coverage, coverage) &&
    identical(cache$observation_basis_full_rank, full_rank)
  if (!valid) {
    cache$observation_basis <- mgcvst_inla_sparse_observation_basis_cpp(
      cache$general_A, as.numeric(geometry$constraint), coverage, full_rank,
      cache$prepared
    )
    cache$observation_basis_A <- geometry$A
    cache$observation_basis_Q <- geometry$Q
    cache$observation_basis_constraint <- geometry$constraint
    cache$observation_basis_coverage <- coverage
    cache$observation_basis_full_rank <- full_rank
  }
  cache$observation_basis
}

.inlast_sparse_nuisance_precision <- function(fit, features) {
  spec <- fit$model$inla_spec
  if (is.null(spec)) spec <- fit$inla_spec
  if (is.null(spec)) return(NULL)
  blocks <- which(vapply(spec$random, function(z) !isTRUE(z$target), logical(1L)))
  if (!length(blocks)) return(NULL)
  width <- vapply(spec$random[blocks], function(z) ncol(z$A), integer(1L))
  sp <- vapply(spec$random[blocks], function(z) as.integer(z$sp_index), integer(1L))
  precision <- t(fit$smoothing_parameters[features, sp, drop = FALSE] /
                   as.numeric(fit$dispersion[features]))
  .inlast_iid_nuisance_precision(ncol(spec$fixed$X), width, precision, length(features))
}

.inlast_null_nuisance_precision <- function(spec, smoothing_parameters,
                                            dispersion, features) {
  blocks <- which(vapply(spec$random, function(z) !isTRUE(z$target), logical(1L)))
  if (!length(blocks)) return(NULL)
  width <- vapply(spec$random[blocks], function(z) ncol(z$A), integer(1L))
  sp <- vapply(spec$random[blocks], function(z) as.integer(z$sp_index), integer(1L))
  precision <- t(smoothing_parameters[features, sp, drop = FALSE] /
                   as.numeric(dispersion[features]))
  .inlast_iid_nuisance_precision(ncol(spec$fixed$X), width, precision, length(features))
}

.inlast_sparse_null_batch <- function(score_sparse, nuisance_design, null_state,
                                      features, threads = 1L) {
  score_fit <- .inlast_sparse_prepare(list(score_sparse = score_sparse))
  geometry <- score_fit$score_sparse
  scale <- geometry$cache$penalty_norm
  Q <- geometry$cache$general_Q
  nuisance_precision <- null_state$nuisance_precision
  if (!is.null(nuisance_precision)) {
    nuisance_precision <- nuisance_precision[, features, drop = FALSE]
  }
  out <- mgcvst_inla_sparse_batch_cpp(
    geometry$cache$general_A, Q,
    as.numeric(geometry$constraint), as.matrix(nuisance_design),
    null_state$working_error[, features, drop = FALSE],
    null_state$working_variance[, features, drop = FALSE],
    rep.int(1 / scale, length(features)), as.integer(threads), FALSE, TRUE,
    32L, geometry$cache$prepared, nuisance_precision = nuisance_precision
  )
  for (j in seq_along(out)) {
    out[[j]]$width <- stats::setNames(length(out[[j]]$a), geometry$target)
    out[[j]]$normalization <- ncol(Q) - 1L
    out[[j]]$backend <- "sparse_conditioned_INLA_OpenMP"
  }
  out
}

.inlast_sparse_batch <- function(fit, features, threads = 1L,
                                 score_only = FALSE) {
  fit <- .inlast_sparse_prepare(fit)
  geometry <- fit$score_sparse
  if (!is.list(geometry) || is.null(geometry$Q)) {
    stop("The fit lacks its sparse INLA score geometry.")
  }
  Q <- geometry$cache$general_Q
  phi <- as.numeric(fit$dispersion[features])
  tau <- as.numeric(fit$smoothing_parameters[features, geometry$sp_index]) / phi
  if (any(!is.finite(tau)) || any(tau <= 0)) {
    stop("The feature has invalid dispersion or smoothing parameters.")
  }
  out <- mgcvst_inla_sparse_batch_cpp(
    geometry$cache$general_A, Q,
    as.numeric(geometry$constraint), as.matrix(fit$geometry$nuisance_design),
    fit$working_error[, features, drop = FALSE],
    fit$working_variance[, features, drop = FALSE], tau,
    as.integer(threads), score_only, FALSE, 32L, geometry$cache$prepared,
    nuisance_precision = .inlast_sparse_nuisance_precision(fit, features)
  )
  for (j in seq_along(out)) {
    out[[j]]$width <- stats::setNames(length(out[[j]]$a), geometry$target)
    out[[j]]$normalization <- ncol(Q) - 1L
    out[[j]]$backend <- "sparse_conditioned_INLA_OpenMP"
  }
  out
}

.mgcvst_model_sparse_score_state <- function(fit, feature, score_only = FALSE) {
  ans <- .inlast_sparse_batch(fit, feature, 1L, score_only)[[1L]]
  if (!is.null(ans$error) && nzchar(ans$error)) stop(ans$error)
  ans
}

.inlast_sparse_units <- function(fit, features, threads = 1L) {
  fit <- .inlast_sparse_prepare(fit)
  geometry <- fit$score_sparse
  tau <- as.numeric(fit$smoothing_parameters[features, geometry$sp_index]) /
    as.numeric(fit$dispersion[features])
  mgcvst_inla_sparse_units_cpp(
    geometry$cache$general_A, geometry$cache$general_Q,
    as.numeric(geometry$constraint), as.matrix(fit$geometry$nuisance_design),
    fit$working_error[, features, drop = FALSE],
    fit$working_variance[, features, drop = FALSE], tau,
    as.integer(threads), geometry$cache$prepared,
    nuisance_precision = .inlast_sparse_nuisance_precision(fit, features)
  )
}

.inlast_sparse_materialize <- function(fit, units, threads = 1L) {
  fit <- .inlast_sparse_prepare(fit)
  geometry <- fit$score_sparse
  mgcvst_inla_sparse_materialize_cpp(
    units, geometry$cache$general_Q, as.numeric(geometry$constraint),
    as.integer(threads), geometry$cache$prepared, 32L
  )
}

.inlast_sparse_materialize_reduced <- function(fit, units, basis,
                                                threads = 1L) {
  fit <- .inlast_sparse_prepare(fit)
  geometry <- fit$score_sparse
  mgcvst_inla_sparse_materialize_reduced_cpp(
    units, geometry$cache$general_Q, as.numeric(geometry$constraint),
    basis$coordinate, basis$basis, as.integer(threads), geometry$cache$prepared
  )
}
