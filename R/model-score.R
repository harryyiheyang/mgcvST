# Return a full-rank column subset without changing its span.
.mgcvst_full_rank_design <- function(X) {
  X <- as.matrix(X)
  if (!ncol(X)) return(X)
  z <- qr(X, tol = 1e-10)
  if (!z$rank) return(matrix(numeric(), nrow(X), 0L))
  X[, z$pivot[seq_len(z$rank)], drop = FALSE]
}

# Factor one positive-definite fixed-kappa SPDE covariance.
.mgcvst_spde_factor <- function(B, Q, scale) {
  sqrt(scale) * .mgcvst_spde_factor_base(B, Q)
}

.mgcvst_spde_factor_base <- function(B, Q) {
  Q <- Matrix::forceSymmetric(Matrix::Matrix(Q, sparse = TRUE))
  R <- Matrix::chol(Q)
  Rinv <- Matrix::solve(R, Matrix::Diagonal(nrow(Q)))
  .magic_mm(B, as.matrix(Rinv))
}

# Per-test geometry only. Cache conditions without moving failures out of
# the original per-pair tryCatch or changing feature-validation precedence.
.mgcvst_model_fixed_factors <- function(fit) {
  lapply(fit$geometry$smooth, function(s) {
    if (s$fixed || is.null(s$score_component) || length(s$penalties) != 1L ||
        length(s$sp_index) != 1L) return(NULL)
    tryCatch(.mgcvst_spde_factor_base(s$B, s$penalties[[1L]]),
             error = function(e) e)
  })
}

# Construct one feature's full nuisance covariance and marked score factors.
.mgcvst_model_operator <- function(fit, feature) {
  geometry <- fit$geometry
  phi <- as.numeric(fit$dispersion[feature])
  sp <- as.numeric(fit$smoothing_parameters[feature, ])
  if (!is.finite(phi) || phi <= 0 || any(!is.finite(sp)) || any(sp <= 0)) {
    stop("The feature has invalid dispersion or smoothing parameters.")
  }
  X <- geometry$X
  factors <- vector("list", length(geometry$smooth))
  target <- vector("list", length(geometry$target))
  names(target) <- names(geometry$target)
  for (j in seq_along(geometry$smooth)) {
    s <- geometry$smooth[[j]]
    if (s$fixed) {
      X <- cbind(X, s$B)
      next
    }
    penalty <- matrix(0, ncol(s$B), ncol(s$B))
    for (k in seq_along(s$penalties)) {
      penalty <- penalty + sp[s$sp_index[k]] * s$penalties[[k]]
    }
    score_name <- s$score_component
    if (!is.null(score_name)) {
      if (length(s$penalties) != 1L || length(s$sp_index) != 1L) {
        stop("A marked SPDE score component must have one penalty.")
      }
      base <- fit$.mgcvst_fixed_factors[[j]]
      if (inherits(base, "condition")) stop(base)
      F <- if (is.null(base)) {
        .mgcvst_spde_factor(s$B, s$penalties[[1L]], phi / sp[s$sp_index])
      } else {
        sqrt(phi / sp[s$sp_index]) * base
      }
      target[[score_name]] <- F
      factors[[j]] <- F
      next
    }
    E <- CppMatrix::matrixEigen((penalty + t(penalty)) / 2)
    d <- as.numeric(E$values)
    V <- as.matrix(E$vectors)
    tol <- sqrt(.Machine$double.eps) * max(1, max(abs(d)))
    if (min(d) < -tol) stop("A fitted nuisance penalty is not positive semidefinite.")
    positive <- d > tol
    if (any(positive)) {
      factors[[j]] <- sqrt(phi) * .magic_mm(
        s$B,
        sweep(V[, positive, drop = FALSE], 2L, 1 / sqrt(d[positive]), "*")
      )
    }
    if (any(!positive)) {
      X <- cbind(X, .magic_mm(s$B, V[, !positive, drop = FALSE]))
    }
  }
  if (any(vapply(target, is.null, logical(1L)))) {
    stop("The feature is missing a marked SPDE score factor.")
  }
  factors <- Filter(Negate(is.null), factors)
  T <- do.call(cbind, factors)
  X <- .mgcvst_full_rank_design(X)
  operator <- .rkhs_score_operator_factor(
    T, fit$working_variance[, feature], X,
    field_scale = 1, B = NULL, Q = NULL
  )
  list(operator = operator, target = target)
}

# Low-rank score state for all marked components of one feature.
.mgcvst_model_score_state <- function(fit, feature) {
  z <- .mgcvst_model_operator(fit, feature)
  F <- do.call(cbind, z$target)
  Pe <- rkhs_score_apply_P(z$operator, fit$working_error[, feature])
  PF <- rkhs_score_apply_P(z$operator, F)
  width <- vapply(z$target, ncol, integer(1L))
  a <- as.numeric(.magic_mm(F, matrix(Pe, ncol = 1L), transA = TRUE))
  M <- .magic_mm(F, PF, transA = TRUE)
  list(
    a = a,
    M = (M + t(M)) / 2,
    width = width,
    target = z$target,
    operator = z$operator
  )
}

# Lazy, successful-state-only cache bounded by the current pair chunk.
.mgcvst_model_cached_state <- function(fit, feature) {
  cache <- fit$.mgcvst_state_cache
  if (is.null(cache)) return(.mgcvst_model_score_state(fit, feature))
  key <- as.character(feature)
  if (!exists(key, envir = cache, inherits = FALSE)) {
    state <- .mgcvst_model_score_state(fit, feature)
    assign(key, state[c("a", "M")], envir = cache)
  }
  get(key, envir = cache, inherits = FALSE)
}

# Single marked-SPDE score with every other smoother in the marginal V.
.mgcvst_model_pair_single <- function(fit, i, j, calibration) {
  s1 <- .mgcvst_model_cached_state(fit, i)
  s2 <- .mgcvst_model_cached_state(fit, j)
  score <- as.numeric(crossprod(s1$a, s2$a))
  cal <- rkhs_score_calibrate(score, s1$M, s2$M, method = calibration)
  list(
    score = score, information = cal$information,
    effective_rank = cal$effective_rank,
    p_two_sided = cal$p_two_sided, p_positive = cal$p_positive,
    p_negative = cal$p_negative
  )
}
