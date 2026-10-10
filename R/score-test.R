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

# Return a numerical factor for a positive-semidefinite score matrix.
.psd_factor <- function(H) {
  H <- (H + t(H)) / 2
  E <- CppMatrix::matrixEigen(H)
  d <- as.numeric(E$values)
  if (min(d) < -1e-10) stop("H is not positive semidefinite.")
  d[d < 0] <- 0
  keep <- d > 0
  if (!any(keep)) return(matrix(numeric(0), nrow(H), 0L))
  sweep(as.matrix(E$vectors[, keep, drop = FALSE]), 2L,
        sqrt(pmax(d[keep], 0)), "*")
}

#' Singular values governing the Gaussian null score distribution
#'
#' @param H1,H2 Aligned innovation-space score covariance summaries.
#' @return The positive singular values in decreasing order.
#' @export
rkhs_score_singular_values <- function(H1, H2) {
  H1 <- .as_numeric_matrix(H1, "H1")
  H2 <- .as_numeric_matrix(H2, "H2")
  if (!all(dim(H1) == dim(H2)) || nrow(H1) != ncol(H1)) {
    stop("H1 and H2 must be square matrices with identical dimensions.")
  }
  .rkhs_score_spectrum(H1, H2)
}

#' Calibrate a signed bilinear Gaussian score
#'
#' Under the Gaussian null, `U` has the distribution
#' `sum(s * Z * W)`, equivalently a signed quadratic form with weights
#' `c(s / 2, -s / 2)`. The calibration is the Lugannani-Rice saddlepoint
#' approximation of the exact distribution on the full singular spectrum `s`
#' of the pair, computed in log space. It is the single-pair reference of the
#' saddlepoint calibration that [mgcvST.test()] and [inlaST.test()] apply to
#' a shared basis of `k` leading singular values: with `k` equal to the
#' dimension, the pair tests reproduce this value.
#'
#' @param U Observed bilinear score.
#' @param H1,H2 Score covariance summaries.
#' @return Simultaneous two-sided, positive, and negative p-values and their
#'   natural logarithms (`log_p_two_sided`, `log_p_positive`,
#'   `log_p_negative`), with the information, effective rank, trace moments
#'   `tr((H1 H2)^j)`, `j = 1, ..., 4`, and the saddlepoint record `spa`
#'   (number of singular values, `remainder_kind` and `status`).
#' @export
rkhs_score_calibrate <- function(U, H1, H2) {
  U <- as.numeric(U)
  if (length(U) != 1L || !is.finite(U)) stop("U must be finite.")
  H1 <- .as_numeric_matrix(H1, "H1")
  H2 <- .as_numeric_matrix(H2, "H2")
  if (!all(dim(H1) == dim(H2)) || nrow(H1) != ncol(H1) || nrow(H1) < 1L) {
    stop("H1 and H2 must be non-empty square matrices with identical dimensions.")
  }
  scale1 <- max(abs(H1))
  scale2 <- max(abs(H2))
  units <- sqrt(scale1) * sqrt(scale2)
  normalized <- if (scale1 > 0 && scale2 > 0) {
    .rkhs_score_moments(H1 / scale1, H2 / scale2)
  } else rep(0, 4L)
  moments <- normalized
  # Restore units in steps, avoiding premature overflow/underflow of units^k.
  for (j in 1:4) for (k in seq_len(2L * j)) moments[j] <- moments[j] * units
  information <- moments[1L]
  invalid <- list(
    p_two_sided = NA_real_, p_positive = NA_real_, p_negative = NA_real_,
    log_p_two_sided = NA_real_, log_p_positive = NA_real_,
    log_p_negative = NA_real_, information = information, effective_rank = 0,
    moments = moments, spa = NULL
  )
  if (!all(is.finite(normalized)) || any(normalized <= 0)) return(invalid)
  s <- .rkhs_score_spectrum(H1 / scale1, H2 / scale2)
  if (!length(s)) return(invalid)
  spectrum_sums <- vapply(1:4, function(r) sum(s^(2 * r)), numeric(1L))
  z <- mgcvst_spa_cpp(U / units, matrix(s, ncol = 1L),
                      matrix(spectrum_sums, 4L, 1L), 4L, 1L)
  if (z[1L, "status"] != 0) return(invalid)
  list(
    p_two_sided = unname(exp(z[1L, "log_p_two_sided"])),
    p_positive = unname(exp(z[1L, "log_p_positive"])),
    p_negative = unname(exp(z[1L, "log_p_negative"])),
    log_p_two_sided = unname(z[1L, "log_p_two_sided"]),
    log_p_positive = unname(z[1L, "log_p_positive"]),
    log_p_negative = unname(z[1L, "log_p_negative"]),
    information = information,
    effective_rank = normalized[1L]^2 / normalized[2L],
    moments = moments,
    spa = list(k = length(s), remainder_kind = unname(z[1L, "remainder_kind"]),
               status = unname(z[1L, "status"]))
  )
}

# Full singular spectrum of the pair (H1, H2): the singular values of
# F1' F2 for PSD factors F1 F1' = H1 and F2 F2' = H2.
.rkhs_score_spectrum <- function(H1, H2) {
  F1 <- .psd_factor(H1)
  F2 <- .psd_factor(H2)
  if (ncol(F1) == 0L || ncol(F2) == 0L) return(numeric(0))
  s <- as.numeric(CppMatrix::matrixSVD(.magic_mm(F1, F2, transA = TRUE))$d)
  s[s > 0]
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

#' Low-rank RKHS covariance score test
#'
#' @param error1,error2 Working-response error vectors for two features.
#' @param operator1,operator2 Compact marginal score operators.
#' @param score_factor1,score_factor2 Optional aligned factors defining
#'   `C12 = score_factor1 %*% t(score_factor2)`.
#' @return An object of class `rkhs_covariance_score`. `signed_score` retains
#'   `U`, and `statistic` is the primary quadratic statistic `U^2`.
#' @export
rkhs_covariance_score <- function(error1, error2, operator1, operator2,
                                  score_factor1 = NULL,
                                  score_factor2 = NULL) {
  S1 <- rkhs_score_summary(error1, operator1, score_factor1)
  S2 <- rkhs_score_summary(error2, operator2, score_factor2)
  if (length(S1$a) != length(S2$a)) {
    stop("The two score factors must use aligned innovation coordinates.")
  }
  U <- as.numeric(crossprod(S1$a, S2$a))
  cal <- rkhs_score_calibrate(U, S1$H, S2$H)
  structure(
    c(list(
      signed_score = U,
      statistic = U^2,
      calibration = "saddlepoint",
      summary1 = S1,
      summary2 = S2
    ),
      cal),
    class = "rkhs_covariance_score"
  )
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
  cat("  quadratic statistic:", format(x$statistic), "\n")
  cat("  information:", format(x$information), "\n")
  cat("  calibration:", x$calibration, "\n")
  cat("  two-sided p-value:", format.pval(x$p_two_sided), "\n")
  cat("  positive p-value:", format.pval(x$p_positive), "\n")
  cat("  negative p-value:", format.pval(x$p_negative), "\n")
  invisible(x)
}
