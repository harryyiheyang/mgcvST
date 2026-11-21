# Sparse score backend for one fixed-kappa SPDE target.  The fitted null
# covariance retains the exact observation mean constraint and the nuisance
# adjustment uses the expected working curvature.

.spdest_sparse_score_capability <- function(model) {
  spec <- model$spde_spec
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

.spdest_sparse_score_geometry <- function(model) {
  capability <- .spdest_sparse_score_capability(model)
  if (!capability$eligible) {
    stop("The sparse score backend is unavailable: ", capability$reason, ".")
  }
  .spdest_score_geometry_from_spec(model$spde_spec)
}

# The sparse score geometry is a function of the specification alone, so a
# worker rebuilds it from the spec it already holds.
.spdest_score_geometry_from_spec <- function(spec) {
  block <- spec$random[[1L]]
  A <- methods::as(block$A, "CsparseMatrix")
  Q <- Matrix::forceSymmetric(methods::as(block$Q, "CsparseMatrix"))
  g <- as.numeric(block$constraint)
  list(
    A = A, Q = Q, constraint = g,
    cache = new.env(parent = emptyenv()),
    target = "global", sp_index = as.integer(block$sp_index),
    normalization = ncol(Q) - 1L,
    definition = "constrained sparse SPDE; expected curvature; penalized nuisance Vp"
  )
}

.spdest_sparse_prepare <- function(fit) {
  geometry <- fit$score_sparse
  if (!is.list(geometry) || is.null(geometry$Q)) {
    stop("The fit lacks its sparse score geometry.")
  }
  if (!is.environment(geometry$cache)) geometry$cache <- new.env(parent = emptyenv())
  cache <- geometry$cache
  valid <- identical(cache$Q, geometry$Q) &&
    identical(cache$constraint, geometry$constraint) &&
    mgcvst_sparse_prepared_valid_cpp(cache$prepared)
  if (!valid) {
    Q <- methods::as(methods::as(geometry$Q, "generalMatrix"), "CsparseMatrix")
    cache$prepared <- mgcvst_sparse_prepare_cpp(Q, as.numeric(geometry$constraint))
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
  X <- fit$geometry$nuisance_design
  if (!is.null(X) && !identical(cache$X, X)) {
    cache$general_X <- methods::as(methods::as(methods::as(
      as.matrix(X), "CsparseMatrix"), "generalMatrix"), "dMatrix")
    cache$X <- X
  }
  fit$score_sparse <- geometry
  fit
}

# The observation-kernel basis shared by every pair method and by
# spdeST.wgcna(): all q = m - 1 directions of the constrained field, ordered by
# their observation-kernel eigenvalue. No coverage truncation is applied, so the
# projection is an orthogonal change of coordinates of the m-dimensional
# whitened score and the WGCNA normalization is the same m - 1.
.spdest_sparse_observation_basis <- function(fit) {
  fit <- .spdest_sparse_prepare(fit)
  geometry <- fit$score_sparse
  cache <- geometry$cache
  valid <- identical(cache$observation_basis_A, geometry$A) &&
    identical(cache$observation_basis_Q, geometry$Q) &&
    identical(cache$observation_basis_constraint, geometry$constraint)
  if (!valid) {
    basis <- mgcvst_sparse_observation_basis_cpp(
      cache$general_A, as.numeric(geometry$constraint), cache$prepared
    )
    basis$kind <- "full_rank"
    cache$observation_basis <- basis
    cache$observation_basis_A <- geometry$A
    cache$observation_basis_Q <- geometry$Q
    cache$observation_basis_constraint <- geometry$constraint
  }
  cache$observation_basis
}

# The test and WGCNA must use the observation basis that the estimation
# recorded: the fit stores its kind and rank, and a different basis is refused.
.spdest_check_basis <- function(fit, basis) {
  recorded <- fit$basis_spec
  if (is.null(recorded) || !identical(recorded$kind, basis$kind) ||
      !identical(as.integer(recorded$rank), as.integer(basis$rank))) {
    stop("The observation basis (", basis$kind, ", rank ", basis$rank,
         ") differs from the one recorded by the estimation (",
         if (is.null(recorded)) "none" else paste0(recorded$kind, ", rank ",
                                                    recorded$rank),
         "); re-run spdeST.scan().", call. = FALSE)
  }
  invisible(basis)
}

.spdest_sparse_nuisance_precision <- function(fit, features) {
  spec <- fit$model$spde_spec
  if (is.null(spec)) spec <- fit$spde_spec
  if (is.null(spec)) return(NULL)
  .spdest_nuisance_precision(spec, fit$smoothing_parameters, fit$dispersion, features)
}

.spdest_nuisance_precision <- function(spec, smoothing_parameters,
                                       dispersion, features) {
  blocks <- which(vapply(spec$random, function(z) !isTRUE(z$target), logical(1L)))
  if (!length(blocks)) return(NULL)
  width <- vapply(spec$random[blocks], function(z) ncol(z$A), integer(1L))
  sp <- vapply(spec$random[blocks], function(z) as.integer(z$sp_index), integer(1L))
  precision <- t(smoothing_parameters[features, sp, drop = FALSE] /
                   as.numeric(dispersion[features]))
  .spdest_iid_nuisance_precision(ncol(spec$fixed$X), width, precision, length(features))
}

# Full-space expected-curvature score vectors a_j of the fitted spatial models,
# computed once at estimation from the working models of `fits`.
.spdest_estimate_scores <- function(fits, score_sparse, spec, nuisance_design,
                                    dispersion, smoothing_parameters, features,
                                    threads = 1L, block = 64L) {
  score_fit <- .spdest_sparse_prepare(list(score_sparse = score_sparse))
  geometry <- score_fit$score_sparse
  a <- matrix(NA_real_, ncol(geometry$Q), length(fits))
  error <- rep(NA_character_, length(fits))
  tau <- as.numeric(smoothing_parameters[, geometry$sp_index]) / as.numeric(dispersion)
  bad <- features[!is.finite(tau[features]) | tau[features] <= 0]
  error[bad] <- "The feature has invalid dispersion or smoothing parameters."
  features <- setdiff(features, bad)
  X <- as.matrix(nuisance_design)
  for (rows in split(seq_along(features), ceiling(seq_along(features) / block))) {
    ids <- features[rows]
    out <- mgcvst_sparse_batch_cpp(
      geometry$cache$general_A, geometry$cache$general_Q,
      as.numeric(geometry$constraint), X,
      do.call(cbind, lapply(fits[ids], `[[`, "working_error")),
      do.call(cbind, lapply(fits[ids], `[[`, "working_variance")),
      tau[ids], as.integer(threads), 32L, geometry$cache$prepared,
      nuisance_precision = .spdest_nuisance_precision(
        spec, smoothing_parameters, dispersion, ids
      )
    )
    for (k in seq_along(ids)) {
      if (!is.null(out[[k]]$error)) error[ids[k]] <- out[[k]]$error else
        a[, ids[k]] <- out[[k]]$a
    }
  }
  list(a = a, error = error)
}

.spdest_family_code <- function(family) {
  code <- match(family, c("gaussian", "poisson", "negative_binomial")) - 1L
  if (anyNA(code)) stop("Unknown or missing feature family.")
  code
}

# Saved compact working-model inputs of `features` for the shared
# coefficient -> eta -> mu -> V step.
.spdest_compact_inputs <- function(fit, features, penalty = TRUE) {
  features <- as.integer(features)
  offset <- fit$offset
  n <- nrow(fit$score_sparse$A)
  O <- if (is.null(offset)) matrix(0, n, 1L) else if (is.matrix(offset))
    t(offset[features, , drop = FALSE]) else matrix(as.numeric(offset), n, 1L)
  storage.mode(O) <- "double"
  size <- vapply(fit$family_parameters[features], function(x) {
    if (length(x)) as.numeric(x[1L]) else NA_real_
  }, numeric(1L))
  B <- fit$target_coefficients[, features, drop = FALSE]
  C <- fit$nuisance_coefficients[, features, drop = FALSE]
  a <- fit$score_a[, features, drop = FALSE]
  dimnames(B) <- dimnames(C) <- dimnames(a) <- NULL
  dispersion <- as.numeric(fit$dispersion[features])
  list(
    B = B, C = C, O = O, family = .spdest_family_code(fit$feature_family[features]),
    size = size, dispersion = dispersion,
    tau = as.numeric(fit$smoothing_parameters[features, fit$score_sparse$sp_index]) /
      dispersion,
    a = a,
    nuisance_precision = if (penalty) .spdest_sparse_nuisance_precision(fit, features)
  )
}

# Double reconstruction units (K, U, Vp, H factor) from the shared working
# state, carrying the saved score vector a_j.
.spdest_sparse_units <- function(fit, features, threads = 1L) {
  fit <- .spdest_sparse_prepare(fit)
  geometry <- fit$score_sparse
  z <- .spdest_compact_inputs(fit, features)
  mgcvst_compact_units_cpp(
    geometry$cache$general_A, geometry$cache$general_Q,
    as.numeric(geometry$constraint), geometry$cache$general_X,
    z$B, z$C, z$O, z$family, z$size, z$dispersion, z$tau, z$a,
    as.integer(threads), geometry$cache$prepared,
    nuisance_precision = z$nuisance_precision
  )
}

.spdest_sparse_materialize_reduced <- function(fit, units, basis,
                                                threads = 1L) {
  fit <- .spdest_sparse_prepare(fit)
  geometry <- fit$score_sparse
  mgcvst_sparse_materialize_reduced_cpp(
    units, geometry$cache$general_Q, as.numeric(geometry$constraint),
    basis$coordinate, basis$basis, as.integer(threads), geometry$cache$prepared
  )
}
