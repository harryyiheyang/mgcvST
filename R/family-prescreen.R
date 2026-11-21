# Per-feature count-family screen.
#
#' Poisson routing for count families
#'
#' @name mgcvst-family-prescreen
#' @keywords internal
#'
#' @section The rule: the overdispersion score:
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
#' There is no user threshold: the earlier `poisson_screen_phi` control, a
#' cutoff on the Pearson dispersion, was removed because the score rule is
#' exact.
#'
#' @section Which family a routed gene actually gets:
#' A routed gene is fitted with the Poisson family (dispersion 1).
NULL

# The removed control is refused by name, so an old script fails loudly
# instead of silently routing differently.
.mgcvst_refuse_screen_control <- function(control) {
  if (!is.null(names(control)) && "poisson_screen_phi" %in% names(control)) {
    stop("control$poisson_screen_phi was removed in spdeST 0.0.1.9047: a ",
         "negative-binomial feature is fitted as Poisson exactly when its ",
         "overdispersion score sum((y - mu)^2 - y) / 2 is not positive.",
         call. = FALSE)
  }
  invisible(control)
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
# Poisson family instead of the negative binomial", exactly when the NB2
# maximum-likelihood fit of the covariate-only model is the Poisson fit.
.mgcvst_prescreen_route <- function(Y, X, offset, active, threads = 1L) {
  p <- nrow(as.matrix(Y))
  out <- list(phi = rep(NA_real_, p), score = rep(NA_real_, p),
              poisson = rep(FALSE, p))
  if (!isTRUE(active)) return(out)
  z <- .mgcvst_prescreen_score(Y, X, offset, threads)
  out$phi <- z$phi
  out$score <- z$score
  out$poisson <- is.finite(z$score) & z$score <= 0
  out
}

