# Resolve and validate the observation-by-innovation score factor.
.score_factor <- function(operator, score_factor) {
  if (is.null(score_factor)) score_factor <- operator$field_factor
  score_factor <- .as_numeric_matrix(score_factor, "score_factor")
  if (nrow(score_factor) != operator$n) {
    stop("score_factor must have operator$n rows.")
  }
  score_factor
}

#' Compute a feature's low-rank score summary
#'
#' @param error Working-response error vector.
#' @param operator A compact score operator.
#' @param score_factor Optional observation-by-innovation factor defining the
#'   tested cross-covariance direction. It defaults to the full RKHS factor in
#'   `operator`.
#' @return A list containing `a`, `H`, and the score factor.
#' @export
rkhs_score_summary <- function(error, operator, score_factor = NULL) {
  error <- as.numeric(error)
  if (length(error) != operator$n || any(!is.finite(error))) {
    stop("error must contain one finite value per observation.")
  }
  F <- .score_factor(operator, score_factor)
  Pe <- rkhs_score_apply_P(operator, error)
  PF <- rkhs_score_apply_P(operator, F)
  a <- as.numeric(.magic_mm(F, matrix(Pe, ncol = 1L), transA = TRUE))
  H <- .magic_mm(F, PF, transA = TRUE)
  H <- (H + t(H)) / 2
  list(a = a, H = H, factor = F)
}

#' Compute low-rank Fisher information
#'
#' @param H1,H2 Aligned innovation-space covariance summaries.
#' @return The scalar information `tr(H1 H2)`.
#' @export
rkhs_score_information <- function(H1, H2) {
  H1 <- .as_numeric_matrix(H1, "H1")
  H2 <- .as_numeric_matrix(H2, "H2")
  if (!all(dim(H1) == dim(H2)) || nrow(H1) != ncol(H1)) {
    stop("H1 and H2 must be square matrices with identical dimensions.")
  }
  as.numeric(sum(H1 * t(H2)))
}

# Return a numerical factor for a positive-semidefinite matrix.
.psd_factor <- function(H, tol) {
  H <- (H + t(H)) / 2
  E <- CppMatrix::matrixEigen(H)
  d <- as.numeric(E$values)
  cutoff <- tol * max(1, max(abs(d)))
  if (min(d) < -cutoff) stop("H is not positive semidefinite.")
  keep <- d > cutoff
  if (!any(keep)) return(matrix(numeric(0), nrow(H), 0L))
  sweep(as.matrix(E$vectors[, keep, drop = FALSE]), 2L,
        sqrt(pmax(d[keep], 0)), "*")
}

#' Singular values governing the Gaussian null score distribution
#'
#' @param H1,H2 Aligned innovation-space score covariance summaries.
#' @param tol Relative numerical-rank tolerance.
#' @return The positive singular values in decreasing order.
#' @export
rkhs_score_singular_values <- function(H1, H2, tol = 1e-10) {
  H1 <- .as_numeric_matrix(H1, "H1")
  H2 <- .as_numeric_matrix(H2, "H2")
  if (!all(dim(H1) == dim(H2)) || nrow(H1) != ncol(H1)) {
    stop("H1 and H2 must be square matrices with identical dimensions.")
  }
  tol <- as.numeric(tol)
  if (length(tol) != 1L || !is.finite(tol) || tol <= 0) {
    stop("tol must be one positive finite value.")
  }
  F1 <- .psd_factor(H1, tol)
  F2 <- .psd_factor(H2, tol)
  if (ncol(F1) == 0L || ncol(F2) == 0L) return(numeric(0))
  A <- .magic_mm(F1, F2, transA = TRUE)
  s <- as.numeric(CppMatrix::matrixSVD(A)$d)
  if (!length(s)) return(numeric(0))
  s[s > tol * max(s)]
}

#' Calibrate a signed bilinear Gaussian score
#'
#' Under the Gaussian null, `U` has the distribution
#' `sum(s * Z * W)`, equivalently a signed quadratic form with weights
#' `c(s / 2, -s / 2)`.
#'
#' @param U Observed bilinear score.
#' @param H1,H2 Score covariance summaries.
#' @param method Either `"liu"` (the default fast calibration) or `"davies"`.
#' @param tol Numerical-rank tolerance.
#' @param davies_accuracy Accuracy passed to [CompQuadForm::davies()].
#' @param davies_limit Integration limit passed to [CompQuadForm::davies()].
#' @return Calibration diagnostics and simultaneous two-sided, positive, and
#'   negative p-values.
#' @export
rkhs_score_calibrate <- function(U, H1, H2,
                                 method = c("liu", "davies"),
                                 tol = 1e-10,
                                 davies_accuracy = 1e-7,
                                 davies_limit = 10000L) {
  method <- match.arg(method)
  davies_accuracy <- as.numeric(davies_accuracy)
  davies_limit <- as.integer(davies_limit)
  if (length(davies_accuracy) != 1L || !is.finite(davies_accuracy) ||
      davies_accuracy <= 0) {
    stop("davies_accuracy must be one positive finite value.")
  }
  if (length(davies_limit) != 1L || is.na(davies_limit) || davies_limit < 1L) {
    stop("davies_limit must be one positive integer.")
  }
  U <- as.numeric(U)
  if (length(U) != 1L || !is.finite(U)) stop("U must be finite.")
  moments <- .rkhs_score_moments(H1, H2)
  information <- moments[1L]
  if (!is.finite(information) || information <= tol) {
    return(list(
      p_value = NA_real_, p_two_sided = NA_real_, p_positive = NA_real_,
      p_negative = NA_real_, method = method, information = information,
      effective_rank = 0, singular_values = numeric(0), moments = moments,
      status = "degenerate_information", davies_ifault = NA_integer_
    ))
  }
  effective_rank <- information^2 / moments[2L]
  s <- numeric(0)

  if (method == "liu") {
    liu <- .liu_squared_score_moments(
      abs(U), moments[1L], moments[2L], moments[3L], moments[4L]
    )
    p.two <- liu$p_value
    ifault <- NA_integer_
  } else {
    if (!requireNamespace("CompQuadForm", quietly = TRUE)) {
      stop("method = 'davies' requires the optional CompQuadForm package.")
    }
    s <- rkhs_score_singular_values(H1, H2, tol = tol)
    if (abs(U) <= tol * sqrt(information)) {
      p.two <- 1
      ifault <- 0L
      liu <- NULL
    } else {
      weights <- c(s / 2, -s / 2)
      fit <- CompQuadForm::davies(
        abs(U), lambda = weights, lim = davies_limit,
        acc = davies_accuracy
      )
      ifault <- as.integer(fit$ifault)
      valid <- ifault == 0L && is.finite(fit$Qq) &&
        fit$Qq >= 0 && fit$Qq <= 0.5 + tol
      if (valid) {
        p.two <- min(1, max(0, 2 * as.numeric(fit$Qq)))
        liu <- NULL
      } else {
        return(list(
          p_value = NA_real_, p_two_sided = NA_real_,
          p_positive = NA_real_, p_negative = NA_real_, method = "davies",
          information = information, effective_rank = effective_rank,
          singular_values = s, moments = moments, status = "davies_failed",
          davies_ifault = ifault, liu_parameters = NULL
        ))
      }
    }
  }

  p.positive <- if (U >= 0) p.two / 2 else 1 - p.two / 2
  p.negative <- if (U <= 0) p.two / 2 else 1 - p.two / 2

  list(
    p_value = p.two,
    p_two_sided = p.two,
    p_positive = p.positive,
    p_negative = p.negative,
    method = method,
    information = information,
    effective_rank = effective_rank,
    singular_values = s,
    moments = moments,
    status = "ok",
    davies_ifault = ifault,
    liu_parameters = liu
  )
}

# Compute trace((H1 H2)^k), k = 1, ..., 4, without a spectral decomposition.
.rkhs_score_moments <- function(H1, H2) {
  H1 <- .as_numeric_matrix(H1, "H1")
  H2 <- .as_numeric_matrix(H2, "H2")
  if (!all(dim(H1) == dim(H2)) || nrow(H1) != ncol(H1)) {
    stop("H1 and H2 must be square matrices with identical dimensions.")
  }
  as.numeric(mgcvst_pair_trace_powers_cpp(
    list(H1, H2), matrix(c(1L, 2L), nrow = 1L), maxPower = 4L,
    threads = 1L
  ))
}

# Match the first four moments of the squared score to a noncentral chi-square.
.liu_squared_score <- function(U, s) {
  .liu_squared_score_moments(
    U, sum(s^2), sum(s^4), sum(s^6), sum(s^8)
  )
}

# Match squared-score moments using trace powers of H1 H2.
.liu_squared_score_moments <- function(U, A, B, C, D) {
  c1 <- A
  c2 <- A^2 + 3 * B
  c3 <- A^3 + 9 * A * B + 15 * C
  c4 <- A^4 + 18 * A^2 * B + 60 * A * C + 24 * B^2 + 105 * D
  s1 <- c3 / c2^(3 / 2)
  s2 <- c4 / c2^2
  tstar <- (U^2 - c1) / sqrt(2 * c2)

  noncentral <- s1^2 > s2
  a <- 1 / s1
  delta <- rep(0, length(s1))
  df <- 1 / s1^2
  if (any(noncentral)) {
    a[noncentral] <- 1 / (
      s1[noncentral] - sqrt(s1[noncentral]^2 - s2[noncentral])
    )
    delta[noncentral] <- s1[noncentral] * a[noncentral]^3 -
      a[noncentral]^2
    df[noncentral] <- a[noncentral]^2 - 2 * delta[noncentral]
  }

  muX <- df + delta
  sigmaX <- sqrt(2) * a
  x <- tstar * sigmaX + muX
  p_value <- stats::pchisq(x, df = df, ncp = delta, lower.tail = FALSE)
  list(
    p_value = pmin(1, pmax(0, p_value)),
    c1 = c1, c2 = c2, c3 = c3, c4 = c4,
    skewness_scale = s1, kurtosis_scale = s2,
    scale = a, df = df, ncp = delta, transformed = x
  )
}

#' Low-rank RKHS covariance score test
#'
#' @param error1,error2 Working-response error vectors for two features.
#' @param operator1,operator2 Compact marginal score operators.
#' @param score_factor1,score_factor2 Optional aligned factors defining
#'   `C12 = score_factor1 %*% t(score_factor2)`.
#' @param method Calibration method passed to [rkhs_score_calibrate()].
#' @param tol Numerical tolerance.
#' @return An object of class `rkhs_covariance_score`. `statistic` and
#'   `quadratic_statistic` are the primary quadratic statistic `U^2`;
#'   `signed_score` (and the backward-compatible `score`) retain `U`.
#'   `normal_statistic = U / sqrt(I)` is populated only when the normal
#'   approximation is explicitly used.
#' @export
rkhs_covariance_score <- function(error1, error2, operator1, operator2,
                                  score_factor1 = NULL,
                                  score_factor2 = NULL,
                                  method = c("liu", "davies"),
                                  tol = 1e-10) {
  method <- match.arg(method)
  S1 <- rkhs_score_summary(error1, operator1, score_factor1)
  S2 <- rkhs_score_summary(error2, operator2, score_factor2)
  if (length(S1$a) != length(S2$a)) {
    stop("The two score factors must use aligned innovation coordinates.")
  }
  U <- as.numeric(crossprod(S1$a, S2$a))
  cal <- rkhs_score_calibrate(
    U, S1$H, S2$H, method = method, tol = tol
  )
  structure(
    c(list(
      score = U,
      signed_score = U,
      statistic = U^2,
      quadratic_statistic = U^2,
      normal_statistic = NA_real_,
      summary1 = S1,
      summary2 = S2
    ),
      cal),
    class = "rkhs_covariance_score"
  )
}

# Direct observation-space reference retained for development verification.
.rkhs_covariance_score_direct <- function(error1, error2, operator1, operator2,
                                          score_factor1 = NULL,
                                          score_factor2 = NULL) {
  F1 <- .score_factor(operator1, score_factor1)
  F2 <- .score_factor(operator2, score_factor2)
  if (ncol(F1) != ncol(F2)) {
    stop("The two score factors must use aligned innovation coordinates.")
  }
  error1 <- as.numeric(error1)
  error2 <- as.numeric(error2)
  if (length(error1) != operator1$n || length(error2) != operator2$n) {
    stop("Each error vector must match its operator dimension.")
  }
  C12 <- .magic_mm(F1, F2, transB = TRUE)
  P1e1 <- rkhs_score_apply_P(operator1, error1)
  P2e2 <- rkhs_score_apply_P(operator2, error2)
  U <- as.numeric(crossprod(P1e1, C12 %*% P2e2))
  P1C12 <- rkhs_score_apply_P(operator1, C12)
  P2C21 <- rkhs_score_apply_P(operator2, t(C12))
  information <- as.numeric(sum(P1C12 * t(P2C21)))
  list(score = U, information = information)
}

#' Print an RKHS covariance score result
#'
#' @param x An `rkhs_covariance_score` object.
#' @param ... Unused.
#' @return `x`, invisibly.
#' @export
print.rkhs_covariance_score <- function(x, ...) {
  cat("Quadratic-form RKHS covariance score test\n")
  cat("  signed score:", format(x$signed_score), "\n")
  cat("  quadratic statistic:", format(x$quadratic_statistic), "\n")
  cat("  information:", format(x$information), "\n")
  cat("  calibration:", x$method, "\n")
  cat("  two-sided p-value:", format.pval(x$p_two_sided), "\n")
  cat("  positive p-value:", format.pval(x$p_positive), "\n")
  cat("  negative p-value:", format.pval(x$p_negative), "\n")
  invisible(x)
}
