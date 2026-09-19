# Per-feature count-family screen shared by the INLA and mgcv paths.
#
#' Poisson prescreen for count families
#'
#' @name mgcvst-family-prescreen
#' @keywords internal
#'
#' @section What phi is:
#' For one count feature we fit a plain Poisson GLM with the model's offset and
#' its parametric (covariate) design only -- the spatial term is deliberately
#' left out -- and report the Pearson dispersion
#'
#'     phi = sum((y - mu)^2 / mu) / (n - p),
#'
#' with `mu` the fitted Poisson mean and `p` the rank of the covariate design.
#' The calculation reuses the parametric design and is performed before the
#' feature fits are distributed.
#'
#' @section How the screen is used:
#' The screen fits an offset-and-covariate-only Poisson GLM and uses its Pearson
#' dispersion as a routing statistic. The spatial term is omitted so that the
#' calculation is inexpensive and identical across the two estimation paths.
#' This statistic is a screening criterion rather than the scale estimate of the
#' subsequent spatial fit. The threshold is calibrated for the package workflow;
#' its default is `1.1`.
#' Genes above the threshold keep the negative binomial family and are fitted
#' exactly as before.
#'
#' @section Which family a routed gene actually gets:
#' The routing decision is shared, the routing target is not. In the **mgcv**
#' path a routed gene is fitted with `stats::quasipoisson(link = "log")`, which
#' `mgcv::gam()` supports: the Poisson and quasipoisson point estimates are the
#' same. The routed mgcv fit estimates its quasipoisson scale and smoothing
#' parameters from the full spatial model. The screening `phi` is retained in
#' diagnostics, whereas the fitted scale enters the working variance and the
#' marginal score test. In the **INLA** path a routed gene keeps plain Poisson,
#' because INLA has no quasi-likelihood families.
#'
#' The Poisson and quasipoisson GLM point estimates agree, so the Poisson fit
#' supplies the routing statistic. The screen selects a family; it does not
#' provide an initial value or fixed scale for the subsequent spatial fit.
NULL

# The package default for control$poisson_screen_phi.
.mgcvst_prescreen_default <- 1.1

# Normalise the knob. Absent/NULL means "use the default"; 0 (the documented
# "off" setting) disables the screen; otherwise one finite positive threshold.
.mgcvst_prescreen_threshold <- function(value) {
  if (is.null(value)) value <- .mgcvst_prescreen_default
  if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
      !is.finite(value) || value < 0) {
    stop("control$poisson_screen_phi must be one finite non-negative number ",
         "(0 disables the Poisson prescreen; NULL restores the default).")
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
  Y <- as.matrix(Y)
  storage.mode(Y) <- "double"
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
    y <- Y[j, ]
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

# Routing decision for a count model: TRUE means "fit this feature with the
# Poisson family instead of the negative binomial".
.mgcvst_prescreen_route <- function(Y, X, offset, threshold, active) {
  p <- nrow(as.matrix(Y))
  out <- list(phi = rep(NA_real_, p), poisson = rep(FALSE, p))
  if (!isTRUE(active) || is.null(threshold)) return(out)
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
