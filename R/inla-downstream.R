#' Test cross-feature spatial covariance from an INLA fit
#'
#' This is the INLA-named entry point for [mgcvST.test()].  It delegates to
#' the score engine registered on the `inlaST_fit` object, including its sparse
#' score state when available.  No model is refitted and no score algorithm is
#' duplicated here.
#' For sparse INLA fits, Liu uses a constrained observation-kernel basis that
#' retains at least 0.995 of its eigenvalue sum. This is a pairwise-test
#' approximation; the fitted sparse field is unchanged. Basis and pair-stage
#' timings are stored in `timing$inla_projection`.
#'
#' @param fitinlaST An object returned by [inlaST.estimate()].
#' @inheritParams mgcvST.test
#' @param calibration Only `"liu"` is supported for the existing INLA pair path.
#' @param chunk_size For conditional pairs, the number of genes materialized
#'   together (1 to 64). The Liu path groups tested pairs.
#' @param BPPARAM Compatibility argument; only `SerialParam()` is accepted.
#' @param threads Positive number of OpenMP threads for sparse INLA feature
#'   preparation and reduced materialization. `NULL` uses one thread.
#' @return The `mgcvST_test` object returned by [mgcvST.test()].
#' @export
inlaST.test <- function(
    fitinlaST, q.value = 0.05, FDR = TRUE, method = "BH",
    BPPARAM = BiocParallel::SerialParam(), ...,
    pairs = NULL, highlight = NULL,
    calibration = "liu",
    chunk_size = NULL, threads = NULL, verbose = FALSE,
    pairwise_method = c("liu", "conditional"),
    checkpoint_dir = NULL, resume = TRUE) {
  pairwise_method <- match.arg(pairwise_method)
  if (pairwise_method == "conditional" && missing(method)) method <- "BY"
  mgcvST.test(
    fitmgcvST = fitinlaST, q.value = q.value, FDR = FDR, method = method,
    BPPARAM = BPPARAM, ..., pairs = pairs, highlight = highlight,
    calibration = calibration, chunk_size = chunk_size,
    threads = threads, verbose = verbose,
    pairwise_method = pairwise_method,
    checkpoint_dir = checkpoint_dir, resume = resume
  )
}

#' Re-evaluate retained marginal tests from an INLA fit
#'
#' This is the INLA-named entry point for [mgcvST.marginal()].  It delegates
#' directly to the existing frozen-state replay implementation.  The original
#' INLA models are not refitted.
#'
#' @param fitinlaST An object returned by [inlaST.estimate()] with
#'   `retain_marginal = TRUE`.
#' @inheritParams mgcvST.marginal
#' @param calibration Only `"liu"` is supported for INLA.
#' @param BPPARAM Compatibility argument; only `SerialParam()` is accepted.
#' @return The data frame returned by [mgcvST.marginal()].
#' @export
inlaST.marginal <- function(
    fitinlaST, features = NULL,
    calibration = "liu", fallback = c("none", "liu"),
    BPPARAM = BiocParallel::SerialParam(), chunk_size = 100L, threads = 1L,
    null.tol = 1e-10, max_eps = 1e-8, max_iter = 1e5) {
  mgcvST.marginal(
    fitmgcvST = fitinlaST, features = features,
    calibration = calibration, fallback = fallback, BPPARAM = BPPARAM,
    chunk_size = chunk_size, threads = threads, null.tol = null.tol,
    max_eps = max_eps, max_iter = max_iter
  )
}
