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
  X <- fit$geometry$nuisance_design
  if (!is.null(X) && !identical(cache$X, X)) {
    cache$general_X <- methods::as(methods::as(methods::as(
      as.matrix(X), "CsparseMatrix"), "generalMatrix"), "dMatrix")
    cache$X <- X
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
  .inlast_nuisance_precision(spec, fit$smoothing_parameters, fit$dispersion, features)
}

.inlast_nuisance_precision <- function(spec, smoothing_parameters,
                                       dispersion, features) {
  blocks <- which(vapply(spec$random, function(z) !isTRUE(z$target), logical(1L)))
  if (!length(blocks)) return(NULL)
  width <- vapply(spec$random[blocks], function(z) ncol(z$A), integer(1L))
  sp <- vapply(spec$random[blocks], function(z) as.integer(z$sp_index), integer(1L))
  precision <- t(smoothing_parameters[features, sp, drop = FALSE] /
                   as.numeric(dispersion[features]))
  .inlast_iid_nuisance_precision(ncol(spec$fixed$X), width, precision, length(features))
}

# Null-model score vectors and full-space curvature trace moments.
.inlast_sparse_null_batch <- function(score_sparse, nuisance_design, E, V,
                                      nuisance_precision, threads = 1L) {
  score_fit <- .inlast_sparse_prepare(list(score_sparse = score_sparse))
  geometry <- score_fit$score_sparse
  scale <- geometry$cache$penalty_norm
  Q <- geometry$cache$general_Q
  out <- mgcvst_inla_sparse_batch_cpp(
    geometry$cache$general_A, Q,
    as.numeric(geometry$constraint), as.matrix(nuisance_design), E, V,
    rep.int(1 / scale, ncol(E)), as.integer(threads), TRUE,
    32L, geometry$cache$prepared, nuisance_precision = nuisance_precision
  )
  for (j in seq_along(out)) {
    out[[j]]$width <- stats::setNames(length(out[[j]]$a), geometry$target)
    out[[j]]$normalization <- ncol(Q) - 1L
    out[[j]]$backend <- "sparse_conditioned_INLA_OpenMP"
  }
  out
}

# Marginal score test of each converged null fit, in chunks of null working
# models; nothing of the null fits is retained.
.inlast_null_marginal <- function(feature_id, score_sparse, nuisance_design,
                                  null_fits, null_spec, dispersion,
                                  smoothing_parameters, features,
                                  chunk_size = 16L, threads = 1L) {
  ans <- data.frame(feature_id = feature_id[features], statistic = NA_real_, p_value = NA_real_,
    method_requested = "liu", method_used = "liu", fallback_used = FALSE,
    fallback_reason = NA_character_, davies_ifault = NA_integer_,
    error_message = NA_character_, stringsAsFactors = FALSE)
  for (rows in split(seq_along(features), ceiling(seq_along(features) / chunk_size))) {
    ids <- features[rows]
    E <- do.call(cbind, lapply(null_fits[ids], `[[`, "working_error"))
    V <- do.call(cbind, lapply(null_fits[ids], `[[`, "working_variance"))
    z <- .inlast_sparse_null_batch(
      score_sparse, nuisance_design, E, V,
      .inlast_nuisance_precision(null_spec, smoothing_parameters, dispersion, ids),
      threads
    )
    for (j in seq_along(rows)) {
      k <- rows[j]
      if (!is.null(z[[j]]$error) && nzchar(z[[j]]$error)) {
        ans$error_message[k] <- z[[j]]$error
        next
      }
      ans$statistic[k] <- z[[j]]$statistic
      ans$p_value[k] <- .mgcvst_marginal_liu(z[[j]]$statistic, z[[j]]$moments)
      if (!is.finite(ans$p_value[k])) ans$error_message[k] <- "Invalid marginal Liu p-value."
    }
  }
  ans
}

# Full-space expected-curvature score vectors a_j of the fitted spatial models,
# computed once at estimation from the working models of `fits`.
.inlast_estimate_scores <- function(fits, score_sparse, spec, nuisance_design,
                                    dispersion, smoothing_parameters, features,
                                    threads = 1L, block = 64L) {
  score_fit <- .inlast_sparse_prepare(list(score_sparse = score_sparse))
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
    out <- mgcvst_inla_sparse_batch_cpp(
      geometry$cache$general_A, geometry$cache$general_Q,
      as.numeric(geometry$constraint), X,
      do.call(cbind, lapply(fits[ids], `[[`, "working_error")),
      do.call(cbind, lapply(fits[ids], `[[`, "working_variance")),
      tau[ids], as.integer(threads), FALSE, 32L, geometry$cache$prepared,
      nuisance_precision = .inlast_nuisance_precision(
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

.inlast_family_code <- function(family) {
  code <- match(family, c("gaussian", "poisson", "negative_binomial")) - 1L
  if (anyNA(code)) stop("Unknown or missing INLA feature family.")
  code
}

# Saved compact working-model inputs of `features` for the shared
# coefficient -> eta -> mu -> V step.
.inlast_compact_inputs <- function(fit, features, penalty = TRUE) {
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
    B = B, C = C, O = O, family = .inlast_family_code(fit$feature_family[features]),
    size = size, dispersion = dispersion,
    tau = as.numeric(fit$smoothing_parameters[features, fit$score_sparse$sp_index]) /
      dispersion,
    a = a,
    nuisance_precision = if (penalty) .inlast_sparse_nuisance_precision(fit, features)
  )
}

# Shared step: eta, mu and working variance (n x features) recovered from the
# saved coefficients and family parameters.
.inlast_working_state <- function(fit, features, threads = 1L) {
  fit <- .inlast_sparse_prepare(fit)
  cache <- fit$score_sparse$cache
  z <- .inlast_compact_inputs(fit, features, penalty = FALSE)
  mgcvst_inla_working_state_cpp(
    cache$general_A, cache$general_X, z$B, z$C, z$O, z$family, z$size,
    z$dispersion, as.integer(threads)
  )
}

# Double reconstruction units (K, U, Vp, H factor) from the shared working
# state, carrying the saved score vector a_j.
.inlast_sparse_units <- function(fit, features, threads = 1L) {
  fit <- .inlast_sparse_prepare(fit)
  geometry <- fit$score_sparse
  z <- .inlast_compact_inputs(fit, features)
  mgcvst_inla_compact_units_cpp(
    geometry$cache$general_A, geometry$cache$general_Q,
    as.numeric(geometry$constraint), geometry$cache$general_X,
    z$B, z$C, z$O, z$family, z$size, z$dispersion, z$tau, z$a,
    as.integer(threads), geometry$cache$prepared,
    nuisance_precision = z$nuisance_precision
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

# Hash of the compact state that determines pair tests of `features`:
# coefficients, saved scores, family and precision parameters, offsets,
# shared geometry and, when supplied, the common projection basis.
.inlast_compact_signature <- function(fit, features, basis = NULL) {
  spec <- fit$model$inla_spec
  sparse <- fit$score_sparse
  offset <- fit$offset
  if (is.matrix(offset)) offset <- offset[features, , drop = FALSE]
  .mgcvst_pair_input_hash(list(
    version = 1L, estimator = fit$estimator, feature_id = fit$feature_id[features],
    target = fit$target_coefficients[, features, drop = FALSE],
    nuisance = fit$nuisance_coefficients[, features, drop = FALSE],
    score = fit$score_a[, features, drop = FALSE],
    family = fit$feature_family[features],
    family_parameters = fit$family_parameters[features],
    dispersion = fit$dispersion[features],
    smoothing = fit$smoothing_parameters[features, , drop = FALSE],
    offset = offset,
    geometry = list(
      A = sparse$A, Q = sparse$Q, constraint = sparse$constraint,
      sp_index = sparse$sp_index, nuisance_design = fit$geometry$nuisance_design,
      fixed_width = ncol(spec$fixed$X),
      random = lapply(spec$random, function(z) {
        list(target = z$target, kind = z$kind, subtype = z$subtype,
             sp_index = z$sp_index, width = ncol(z$A))
      })
    ),
    basis = if (is.null(basis)) NULL else basis[c("coordinate", "basis", "rank", "coverage")]
  ))
}
