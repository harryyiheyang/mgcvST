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
  Q <- Matrix::forceSymmetric(Matrix::Matrix(Q, sparse = TRUE))
  R <- Matrix::chol(Q)
  Rinv <- Matrix::solve(R, Matrix::Diagonal(nrow(Q)))
  sqrt(scale) * .magic_mm(B, as.matrix(Rinv))
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
      F <- .mgcvst_spde_factor(
        s$B, s$penalties[[1L]], phi / sp[s$sp_index]
      )
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

# Calibrate a general bilinear kernel x' C y using an SVD coordinate bridge.
.mgcvst_bilinear_calibrate <- function(a1, a2, M1, M2, C,
                                       calibration, tol) {
  z <- svd(C)
  keep <- z$d > tol * max(1, max(z$d))
  if (!any(keep)) {
    return(list(
      p_value = NA_real_, p_two_sided = NA_real_, p_positive = NA_real_,
      p_negative = NA_real_, method = calibration, information = 0,
      effective_rank = 0, status = "degenerate_information",
      davies_ifault = NA_integer_, U = 0
    ))
  }
  d <- sqrt(z$d[keep])
  U <- z$u[, keep, drop = FALSE]
  V <- z$v[, keep, drop = FALSE]
  a1s <- d * as.numeric(crossprod(U, a1))
  a2s <- d * as.numeric(crossprod(V, a2))
  H1 <- sweep(crossprod(U, M1 %*% U), 1L, d, "*")
  H1 <- sweep(H1, 2L, d, "*")
  H2 <- sweep(crossprod(V, M2 %*% V), 1L, d, "*")
  H2 <- sweep(H2, 2L, d, "*")
  score <- as.numeric(crossprod(a1s, a2s))
  out <- rkhs_score_calibrate(
    score, H1, H2, method = calibration, tol = tol
  )
  out$U <- score
  out
}

# Low-rank score state for all marked components of one feature.
.mgcvst_model_score_state <- function(fit, feature) {
  z <- .mgcvst_model_operator(fit, feature)
  F <- do.call(cbind, z$target)
  Pe <- rkhs_score_apply_P(z$operator, fit$working_error[, feature])
  PF <- rkhs_score_apply_P(z$operator, F)
  width <- vapply(z$target, ncol, integer(1L))
  list(
    a = as.numeric(.magic_mm(F, matrix(Pe, ncol = 1L), transA = TRUE)),
    M = (.magic_mm(F, PF, transA = TRUE) +
           t(.magic_mm(F, PF, transA = TRUE))) / 2,
    width = width,
    target = z$target,
    operator = z$operator
  )
}

# Single marked-SPDE score with every other smoother in the marginal V.
.mgcvst_model_pair_single <- function(fit, i, j, calibration,
                                      tol) {
  s1 <- .mgcvst_model_score_state(fit, i)
  s2 <- .mgcvst_model_score_state(fit, j)
  q <- s1$width[["global"]]
  C <- diag(q)
  cal <- .mgcvst_bilinear_calibrate(
    s1$a, s2$a, s1$M, s2$M, C, calibration, tol
  )
  list(
    score = cal$U, global_score = cal$U, local_score = NA_real_,
    information = cal$information, global_information = cal$information,
    local_information = NA_real_, cross_information = NA_real_,
    efficient_information = cal$information,
    schur_coefficient = NA_real_, schur_fraction = 0,
    effective_rank = cal$effective_rank, p_value = cal$p_two_sided,
    p_two_sided = cal$p_two_sided, p_positive = cal$p_positive,
    p_negative = cal$p_negative,
    calibration = cal$method, status = cal$status,
    davies_ifault = cal$davies_ifault
  )
}

# Lin-style global score projected off the local SPDE score direction.
.mgcvst_model_pair_global_local <- function(fit, i, j, calibration,
                                            tol) {
  s1 <- .mgcvst_model_score_state(fit, i)
  s2 <- .mgcvst_model_score_state(fit, j)
  qg <- s1$width[["global"]]
  ql <- s1$width[["local"]]
  q <- qg + ql
  Cg <- matrix(0, q, q)
  Cl <- matrix(0, q, q)
  Cg[seq_len(qg), seq_len(qg)] <- diag(qg)
  il <- qg + seq_len(ql)
  Cl[il, il] <- diag(ql)
  score <- function(C) as.numeric(crossprod(s1$a, C %*% s2$a))
  info <- function(Ca, Cb) as.numeric(sum(diag(
    s1$M %*% Ca %*% s2$M %*% t(Cb)
  )))
  Ug <- score(Cg)
  Ul <- score(Cl)
  Igg <- info(Cg, Cg)
  Ill <- info(Cl, Cl)
  Igl <- info(Cg, Cl)
  if (!is.finite(Ill) || Ill <= tol) {
    stop("The local SPDE score has no positive nuisance information.")
  }
  coefficient <- Igl / Ill
  Ceff <- Cg - coefficient * Cl
  Ieff <- Igg - Igl^2 / Ill
  if (!is.finite(Ieff) || Ieff <= tol * max(1, Igg)) {
    stop("The Lin-Schur projection leaves no positive global information.")
  }
  cal <- .mgcvst_bilinear_calibrate(
    s1$a, s2$a, s1$M, s2$M, Ceff, calibration, tol
  )
  if (abs(cal$information - Ieff) >
      1e-7 * max(1, abs(cal$information), abs(Ieff))) {
    stop("The calibrated kernel and Schur information do not agree.")
  }
  list(
    score = cal$U, global_score = Ug, local_score = Ul,
    information = Ieff, global_information = Igg,
    local_information = Ill, cross_information = Igl,
    efficient_information = Ieff,
    schur_coefficient = coefficient,
    schur_fraction = max(0, min(1, Igl^2 / (Igg * Ill))),
    effective_rank = cal$effective_rank, p_value = cal$p_two_sided,
    p_two_sided = cal$p_two_sided, p_positive = cal$p_positive,
    p_negative = cal$p_negative,
    calibration = cal$method, status = cal$status,
    davies_ifault = cal$davies_ifault
  )
}
