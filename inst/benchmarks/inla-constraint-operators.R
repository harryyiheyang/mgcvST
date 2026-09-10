# Reference score operators for studying the observation mean-zero constraint.
#
# This file is deliberately independent of the production score path.  It
# compares five precisely different tests while keeping the INLA fit itself
# constrained by g'u = 0:
#
#   projected                  constrained null V, constrained score kernel;
#   raw_constrained            same model in raw mesh coordinates;
#   raw_kernel_only            constrained null V, unconstrained score kernel;
#   raw_full_keep_nuisance     unconstrained null V and kernel, but the fitted
#                              constrained nuisance covariance is naively kept;
#   raw_full_recompute_nuisance unconstrained null V and kernel, with the
#                              unpenalized nuisance GLS covariance recomputed.
#
# For Q = R'R, L = R^-1 and h = L'g, the constrained raw-coordinate factor is
#
#   F_c = tau^-1/2 A L {I - h h'/(h'h)}.
#
# Its covariance equals the projected construction exactly:
#
#   A [Q^-1 - Q^-1 g (g'Q^-1 g)^-1 g'Q^-1] A'
#     = A Z (Z'QZ)^-1 Z'A',                       g'Z = 0.
#
# The unconstrained factor F_u = tau^-1/2 A L adds the rank-one direction
#
#   v = tau^-1/2 A Q^-1 g / sqrt(g'Q^-1 g),
#   F_u F_u' = F_c F_c' + v v'.
#
# An intercept removes this difference only when v is in the nuisance-design
# span.  Although barycentric A satisfies A 1 = 1, Q^-1 g is generally not a
# constant mesh vector, so that equivalence must be checked rather than assumed.
# More specifically, Q^-1 g is constant exactly when g is proportional to Q 1.
# For the present Matern FEM precision, M1 1 = M2 1 = 0 and hence
# Q 1 = kappa^4 M0 1.  Uniformly distributed observations can make the empirical
# barycentric weights g approximately proportional to the lumped FEM masses,
# but a finite irregular spatial sample does not make this identity exact.

.inlast_ref_matrix <- function(x, name) {
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  if (length(dim(x)) != 2L || any(!is.finite(x))) {
    stop(name, " must be a finite numeric matrix.")
  }
  x
}

.inlast_ref_solve <- function(A, B) {
  solve(A, B)
}

# Remove the Euclidean nuisance-column component of a score factor.  For the
# two reference constructions whose P operator is defined to satisfy P X = 0,
# this changes neither F' P e nor F' P F in exact arithmetic.  It avoids a
# catastrophic cancellation when the unconstrained raw factor contains a very
# large near-intercept direction.  QR is used because nuisance_X can contain
# several columns with quite different scales.
.inlast_ref_residualize_score <- function(F, X) {
  F <- .inlast_ref_matrix(F, "score_factor")
  X <- .inlast_ref_matrix(X, "nuisance_X")
  if (!ncol(X)) return(F)
  F - qr.fitted(qr(X), F)
}

.inlast_ref_vsolve <- function(field_factor, working_variance) {
  F <- .inlast_ref_matrix(field_factor, "field_factor")
  D <- as.numeric(working_variance)
  if (length(D) != nrow(F) || any(!is.finite(D)) || any(D <= 0)) {
    stop("working_variance must be positive with one value per observation.")
  }
  Dinv <- 1 / D
  DinvF <- Dinv * F
  woodbury <- diag(ncol(F)) + crossprod(F, DinvF)
  function(value) {
    was_vector <- is.null(dim(value))
    value <- if (was_vector) matrix(as.numeric(value), ncol = 1L) else
      .inlast_ref_matrix(value, "value")
    if (nrow(value) != nrow(F)) stop("value has the wrong observation dimension.")
    DinvY <- Dinv * value
    out <- DinvY - DinvF %*%
      .inlast_ref_solve(woodbury, crossprod(F, DinvY))
    if (was_vector) as.numeric(out) else as.matrix(out)
  }
}

.inlast_ref_state <- function(null_factor, score_factor, working_error,
                              working_variance, nuisance_X, nuisance_Vp,
                              nuisance = c("keep", "recompute")) {
  nuisance <- match.arg(nuisance)
  F0 <- .inlast_ref_matrix(null_factor, "null_factor")
  Fs <- .inlast_ref_matrix(score_factor, "score_factor")
  if (nrow(F0) != nrow(Fs)) stop("The null and score factors must align by observation.")
  e <- as.numeric(working_error)
  if (length(e) != nrow(F0) || any(!is.finite(e))) {
    stop("working_error must be finite with one value per observation.")
  }
  X <- if (is.null(nuisance_X)) {
    matrix(1, nrow(F0), 1L)
  } else {
    .inlast_ref_matrix(nuisance_X, "nuisance_X")
  }
  if (nrow(X) != nrow(F0)) stop("nuisance_X has the wrong observation dimension.")

  Vsolve <- .inlast_ref_vsolve(F0, working_variance)
  if (!ncol(X)) {
    Vp <- matrix(numeric(), 0L, 0L)
    WX <- matrix(numeric(), nrow(F0), 0L)
  } else {
    WX <- Vsolve(X)
    Vp <- if (nuisance == "recompute") {
      .inlast_ref_solve(crossprod(X, WX), diag(ncol(X)))
    } else {
      .inlast_ref_matrix(nuisance_Vp, "nuisance_Vp")
    }
    if (!all(dim(Vp) == ncol(X))) {
      stop("nuisance_Vp must be square with ncol(nuisance_X) rows.")
    }
  }
  Papply <- function(value) {
    was_vector <- is.null(dim(value))
    value <- if (was_vector) matrix(as.numeric(value), ncol = 1L) else
      .inlast_ref_matrix(value, "value")
    out <- Vsolve(value)
    if (ncol(X)) out <- out - WX %*% Vp %*% crossprod(X, out)
    if (was_vector) as.numeric(out) else as.matrix(out)
  }
  Pe <- Papply(e)
  PF <- Papply(Fs)
  a <- as.numeric(crossprod(Fs, Pe))
  M <- crossprod(Fs, PF)
  M <- (M + t(M)) / 2
  list(
    a = a, M = M, null_factor = F0, score_factor = Fs,
    nuisance_covariance = Vp,
    nuisance_annihilation_error = if (ncol(X)) max(abs(Papply(X))) else 0,
    apply_P = Papply
  )
}

# Construct all reference states for one fitted feature.
#
# `tau` is the INLA generic0 precision multiplying raw_Q.  Because mgcvST uses
# lambda = phi * tau, the score field scale phi/lambda is exactly 1/tau for
# Gaussian, Poisson and negative-binomial working models.
inlast_constraint_reference_states <- function(
    raw_A, raw_Q, g, Z, working_error, working_variance, tau,
    nuisance_X = NULL, nuisance_Vp = NULL, tolerance = 1e-9) {
  if (!requireNamespace("Matrix", quietly = TRUE)) stop("Matrix is required.")
  A_sparse <- methods::as(raw_A, "CsparseMatrix")
  Q_sparse <- Matrix::forceSymmetric(methods::as(raw_Q, "CsparseMatrix"))
  A <- .inlast_ref_matrix(A_sparse, "raw_A")
  Q <- .inlast_ref_matrix(Q_sparse, "raw_Q")
  n <- nrow(A)
  m <- ncol(A)
  if (!all(dim(Q) == c(m, m)) || !isTRUE(isSymmetric(Q, tol = tolerance))) {
    stop("raw_Q must be symmetric with ncol(raw_A) rows and columns.")
  }
  tau <- as.numeric(tau)
  if (length(tau) != 1L || !is.finite(tau) || tau <= 0) {
    stop("tau must be one positive finite precision.")
  }
  g <- as.numeric(g)
  expected_g <- as.numeric(crossprod(A, rep.int(1 / n, n)))
  if (length(g) != m || !isTRUE(all.equal(
      g, expected_g, tolerance = tolerance, check.attributes = FALSE))) {
    stop("g must equal crossprod(raw_A, 1 / n).")
  }
  Z <- .inlast_ref_matrix(Z, "Z")
  if (!all(dim(Z) == c(m, m - 1L)) ||
      max(abs(crossprod(Z) - diag(m - 1L))) > tolerance ||
      max(abs(crossprod(g, Z))) > tolerance) {
    stop("Z must be an orthonormal basis for null(t(g)).")
  }
  R <- chol(Q)
  L <- backsolve(R, diag(m))
  h <- as.numeric(crossprod(L, g))
  h2 <- sum(h^2)
  if (!is.finite(h2) || h2 <= 0) stop("The constraint has zero Q-inverse norm.")
  innovation_projection <- diag(m) - tcrossprod(h) / h2
  scale <- 1 / sqrt(tau)
  raw_unconstrained <- scale * A %*% L
  raw_constrained <- raw_unconstrained %*% innovation_projection

  Q_projected <- crossprod(Z, Q %*% Z)
  R_projected <- chol(Q_projected)
  projected <- scale * (A %*% Z) %*%
    backsolve(R_projected, diag(m - 1L))
  extra_direction <- as.numeric(
    scale * A %*% L %*% (h / sqrt(h2))
  )

  X <- if (is.null(nuisance_X)) matrix(1, n, 1L) else
    .inlast_ref_matrix(nuisance_X, "nuisance_X")
  raw_score_residualized <- .inlast_ref_residualize_score(
    raw_unconstrained, X
  )
  if (is.null(nuisance_Vp)) {
    centered_vsolve <- .inlast_ref_vsolve(raw_constrained, working_variance)
    nuisance_Vp <- if (ncol(X)) {
      .inlast_ref_solve(crossprod(X, centered_vsolve(X)), diag(ncol(X)))
    } else matrix(numeric(), 0L, 0L)
  }
  centered_args <- list(
    working_error = working_error, working_variance = working_variance,
    nuisance_X = X, nuisance_Vp = nuisance_Vp, nuisance = "keep"
  )
  raw_args <- centered_args

  states <- list(
    projected = do.call(
      .inlast_ref_state,
      c(list(null_factor = projected, score_factor = projected), centered_args)
    ),
    raw_constrained = do.call(
      .inlast_ref_state,
      c(list(null_factor = raw_constrained,
             score_factor = raw_constrained), centered_args)
    ),
    raw_kernel_only = do.call(
      .inlast_ref_state,
      c(list(null_factor = raw_constrained,
             score_factor = raw_score_residualized), centered_args)
    ),
    raw_full_keep_nuisance = do.call(
      .inlast_ref_state,
      c(list(null_factor = raw_unconstrained,
             score_factor = raw_unconstrained), raw_args)
    ),
    raw_full_recompute_nuisance = .inlast_ref_state(
      null_factor = raw_unconstrained,
      score_factor = raw_score_residualized,
      working_error = working_error,
      working_variance = working_variance,
      nuisance_X = X,
      nuisance_Vp = nuisance_Vp,
      nuisance = "recompute"
    )
  )

  kernel_projected <- tcrossprod(projected)
  kernel_constrained <- tcrossprod(raw_constrained)
  kernel_unconstrained <- tcrossprod(raw_unconstrained)
  nuisance_fit <- if (ncol(X)) qr.fitted(qr(X), extra_direction) else
    numeric(n)
  extra_norm <- sqrt(sum(extra_direction^2))
  list(
    states = states,
    factors = list(
      projected = projected, raw_constrained = raw_constrained,
      raw_unconstrained = raw_unconstrained,
      raw_unconstrained_score_residualized = raw_score_residualized,
      extra_direction = extra_direction
    ),
    diagnostics = list(
      constrained_kernel_relative_error =
        max(abs(kernel_projected - kernel_constrained)) /
        max(1, max(abs(kernel_projected))),
      rank_one_identity_relative_error =
        max(abs(kernel_unconstrained - kernel_constrained -
                  tcrossprod(extra_direction))) /
        max(1, max(abs(kernel_unconstrained))),
      extra_direction_nuisance_residual_ratio = if (extra_norm > 0) {
        sqrt(sum((extra_direction - nuisance_fit)^2)) / extra_norm
      } else 0,
      qinv_g_node_constant_ratio = {
        qinv_g <- as.numeric(L %*% h)
        sqrt(sum((qinv_g - mean(qinv_g))^2)) /
          max(sqrt(sum(qinv_g^2)), .Machine$double.eps)
      },
      raw_score_nuisance_projection_relative_norm =
        sqrt(sum((raw_unconstrained - raw_score_residualized)^2)) /
        max(sqrt(sum(raw_unconstrained^2)), .Machine$double.eps),
      row_sum_error = max(abs(rowSums(A) - 1))
    ),
    assumptions = list(
      fitted_field_constraint = "crossprod(raw_A, 1 / n)' u = 0",
      positive_parameter_prior = "log(parameter) ~ N(0, 3^2)",
      recomputed_nuisance = paste(
        "raw_full_recompute_nuisance treats every nuisance_X column as",
        "unpenalized; the intended type-I simulation uses an intercept only"
      ),
      score_factor_stabilization = paste(
        "raw_kernel_only and raw_full_recompute_nuisance residualize the raw",
        "score factor against nuisance_X by QR; this preserves a and M because",
        "their defining P operators satisfy P nuisance_X = 0"
      )
    )
  )
}

# Combine two feature states using the production calibration routine.  The
# two states must come from the same named construction and aligned raw mesh.
inlast_constraint_reference_pair <- function(state1, state2,
                                              method = c("liu", "davies")) {
  method <- match.arg(method)
  if (length(state1$a) != length(state2$a) ||
      !all(dim(state1$M) == dim(state2$M))) {
    stop("The feature states do not use aligned innovation coordinates.")
  }
  score <- as.numeric(crossprod(state1$a, state2$a))
  calibration <- rkhs_score_calibrate(
    score, state1$M, state2$M, method = method
  )
  c(list(signed_score = score, statistic = score^2), calibration)
}

# Fast analytic checks used before a simulation run.  A constrained raw score
# and the projected score can have rotated innovation coordinates, so compare
# invariant pair scores and calibration matrices rather than their `a` vectors.
inlast_constraint_reference_equivalence <- function(reference1, reference2,
                                                     tolerance = 1e-8) {
  p1 <- inlast_constraint_reference_pair(
    reference1$states$projected, reference2$states$projected, method = "liu"
  )
  c1 <- inlast_constraint_reference_pair(
    reference1$states$raw_constrained,
    reference2$states$raw_constrained, method = "liu"
  )
  errors <- c(
    feature1_kernel =
      reference1$diagnostics$constrained_kernel_relative_error,
    feature2_kernel =
      reference2$diagnostics$constrained_kernel_relative_error,
    signed_score = abs(p1$signed_score - c1$signed_score) /
      max(1, abs(p1$signed_score)),
    information = abs(p1$information - c1$information) /
      max(1, abs(p1$information))
  )
  list(
    passed = all(is.finite(errors)) && max(errors) <= tolerance,
    relative_errors = errors,
    projected = p1,
    raw_constrained = c1
  )
}
