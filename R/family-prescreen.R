# Per-feature count-family screen shared by the INLA and mgcv paths.
#
#' Poisson prescreen for count families
#'
#' @name mgcvst-family-prescreen
#' @keywords internal
#'
#' @section The default: the overdispersion score:
#' For one count feature we fit a plain Poisson GLM with the model's offset and
#' its parametric (covariate) design only -- the spatial term is deliberately
#' left out -- and compute the overdispersion score of the NB2 family
#' (variance `mu + alpha mu^2`) at its Poisson end point,
#'
#'     score = sum((y - mu)^2 - y) / 2,
#'
#' the derivative of the NB2 profile log-likelihood in `alpha` at `alpha = 0`
#' (Cameron and Trivedi). A score that is not positive is the KKT condition of
#' that end point: the negative-binomial maximum-likelihood fit is the Poisson
#' fit, and the feature is routed to the Poisson family. Otherwise the size is
#' finite and the feature keeps the negative binomial. The fits are batched in
#' C++ (src/nb_null.cpp). The Pearson dispersion
#' `phi = sum((y - mu)^2 / mu) / (n - p)` of the same fit is reported as well.
#'
#' @section A phi threshold:
#' `poisson_screen_phi = <number>` restores the earlier screen: a feature is
#' routed when its Pearson dispersion is at most that number. `0` disables the
#' routing; `NULL` or `"score"` is the default above.
#'
#' @section Which family a routed gene actually gets:
#' The routing decision is shared, the routing target is not. In the **mgcv**
#' path a routed gene is fitted with `stats::quasipoisson(link = "log")`, which
#' `mgcv::bam()` supports: the Poisson and quasipoisson point estimates are the
#' same. The routed mgcv null score fit and full spatial fit each estimate their
#' quasipoisson scale and smoothing parameters. In the **INLA** path a routed
#' gene keeps plain Poisson, because INLA has no quasi-likelihood families.
NULL

# The package default for control$poisson_screen_phi: the exact
# overdispersion-score rule.
.mgcvst_prescreen_default <- "score"

# Normalise the knob. Absent/NULL means "use the default"; "score" is the exact
# rule; 0 (the documented "off" setting) disables the screen; otherwise one
# finite positive phi threshold.
.mgcvst_prescreen_threshold <- function(value) {
  if (is.null(value)) value <- .mgcvst_prescreen_default
  if (identical(value, "score")) return("score")
  if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
      !is.finite(value) || value < 0) {
    stop("control$poisson_screen_phi must be \"score\" or one finite ",
         "non-negative number (0 disables the Poisson prescreen; NULL ",
         "restores the default \"score\").")
  }
  value <- as.numeric(value)
  if (value <= 0) return(NULL)
  value
}

# Pearson dispersion of an offset + covariate-only Poisson GLM, one value per
# row of Y. `X` is the dense parametric design (no spatial columns); `offset` is
# NULL, one observation-length vector, or a matrix matching Y. A feature whose
# GLM cannot be fitted returns NA, which never routes (conservative).
.mgcvst_prescreen_dispersion <- function(Y, X, offset = NULL) {
  if (!is.matrix(Y)) Y <- as.matrix(Y)
  n <- ncol(Y)
  X <- if (is.null(X)) matrix(0, n, 0L) else as.matrix(X)
  storage.mode(X) <- "double"
  if (nrow(X) != n) stop("The prescreen design and Y disagree in length.")
  if (!is.null(offset) && !is.matrix(offset) && length(offset) != n) {
    stop("The prescreen offset and Y disagree in length.")
  }
  family <- stats::poisson()
  glm_control <- stats::glm.control(epsilon = 1e-8, maxit = 25L)
  vapply(seq_len(nrow(Y)), function(j) {
    y <- as.numeric(Y[j, ])
    off <- if (is.null(offset)) numeric(n) else
      if (is.matrix(offset)) as.numeric(offset[j, ]) else as.numeric(offset)
    if (!ncol(X)) {
      mu <- exp(off)
      p <- 0L
    } else {
      fit <- tryCatch(
        suppressWarnings(stats::glm.fit(
          x = X, y = y, offset = off, family = family, control = glm_control
        )),
        error = function(e) NULL
      )
      if (is.null(fit)) return(NA_real_)
      mu <- as.numeric(fit$fitted.values)
      p <- as.integer(fit$rank)
    }
    if (length(mu) != n || any(!is.finite(mu)) || any(mu <= 0) || n - p <= 0L) {
      return(NA_real_)
    }
    sum((y - mu)^2 / mu) / (n - p)
  }, numeric(1L))
}

# Overdispersion score and Pearson dispersion of the offset + covariate-only
# Poisson fit of every row of Y, batched in C++. A feature whose fit fails
# returns NA, which never routes (conservative).
.mgcvst_prescreen_score <- function(Y, X, offset = NULL, threads = 1L) {
  if (!methods::is(Y, "Matrix")) Y <- as.matrix(Y)
  n <- ncol(Y)
  p <- nrow(Y)
  X <- if (is.null(X)) matrix(0, n, 0L) else as.matrix(X)
  storage.mode(X) <- "double"
  if (nrow(X) != n) stop("The prescreen design and Y disagree in length.")
  if (!is.null(offset) && !is.matrix(offset) && length(offset) != n) {
    stop("The prescreen offset and Y disagree in length.")
  }
  O <- if (is.null(offset)) matrix(0, n, 1L) else if (is.matrix(offset))
    t(offset) else matrix(as.numeric(offset), n, 1L)
  if (!ncol(X)) {
    # No coefficients: the Poisson mean is the offset itself.
    Yd <- as.matrix(Y)
    M <- exp(if (ncol(O) == 1L) matrix(O[, 1L], p, n, byrow = TRUE) else t(O))
    R <- Yd - M
    return(list(score = rowSums(R^2 - Yd) / 2, phi = rowSums(R^2 / M) / n))
  }
  Yt <- methods::as(methods::as(methods::as(
    Matrix::t(Y), "CsparseMatrix"), "generalMatrix"), "dMatrix")
  z <- mgcvst_nb_null_batch_cpp(Yt, X, O, rep(TRUE, p),
                                threads = as.integer(threads))
  list(score = z$score0, phi = z$phi)
}

# Routing decision for a count model: TRUE means "fit this feature with the
# Poisson family instead of the negative binomial".
.mgcvst_prescreen_route <- function(Y, X, offset, threshold, active,
                                    threads = 1L) {
  p <- nrow(as.matrix(Y))
  out <- list(phi = rep(NA_real_, p), poisson = rep(FALSE, p))
  if (!isTRUE(active) || is.null(threshold)) return(out)
  if (identical(threshold, "score")) {
    z <- .mgcvst_prescreen_score(Y, X, offset, threads)
    out$phi <- z$phi
    out$score <- z$score
    out$poisson <- is.finite(z$score) & z$score <= 0
    return(out)
  }
  out$phi <- .mgcvst_prescreen_dispersion(Y, X, offset)
  out$poisson <- is.finite(out$phi) & out$phi <= threshold
  out
}

# Parametric (non-smooth) columns of a frozen mgcv design. Every smooth block,
# spatial or not, is excluded: the screen is deliberately covariate-only.
.mgcvst_prescreen_design <- function(G) {
  X <- G$X
  if (is.null(X)) return(NULL)
  smooth_columns <- unlist(lapply(G$smooth, function(s) {
    seq.int(s$first.para, s$last.para)
  }), use.names = FALSE)
  keep <- setdiff(seq_len(ncol(X)), smooth_columns)
  as.matrix(X[, keep, drop = FALSE])
}
