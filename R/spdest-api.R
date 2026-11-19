# spdeST user interface: one pipeline for the mgcv and INLA implementations.
#
#   model <- mgcvST.set(...)  or  inlaST.set(...)   # every fixed setting
#   scan  <- spdeST.scan(Y, model, ...)              # Stage 1, both engines
#   fit   <- spdeST.estimate(scan, Y, ...)           # Stage 2, chosen features
#   res   <- spdeST.test(fit, ...)                   # pair test
#
# The model object is the pointer: its class selects the implementation, by
# S3 dispatch of the internal generics below to R/spdest-mgcv.R and
# R/spdest-inla.R. This file holds no engine-specific code.

.spdest_flag <- function(value, name) {
  if (!is.logical(value) || length(value) != 1L || is.na(value)) {
    stop(name, " must be TRUE or FALSE.")
  }
  value
}

# Only the arguments the caller supplied, so that an implementation keeps its
# own defaults for the others.
.spdest_supplied <- function(...) {
  args <- list(...)
  args[!vapply(args, is.null, logical(1L))]
}

# spde.scan() controls supplied by name; the others keep its defaults.
.spdest_scan_control <- function(scan_control) {
  if (!is.list(scan_control)) stop("scan_control must be a named list.")
  if (!length(scan_control)) return(list())
  allowed <- c("tol", "maxit_outer", "maxit_inner", "alpha_max", "alpha_min",
               "step_bound", "step_fallback", "maxit_halving", "eta_max", "rank_tol",
               "spectrum_tol", "spa_centre_tol", "spa_root_tol",
               "spa_root_maxit", "block_size")
  bad <- setdiff(names(scan_control), allowed)
  if (is.null(names(scan_control)) || any(!nzchar(names(scan_control))) || length(bad)) {
    stop("scan_control takes named spde.scan() controls: ",
         paste(allowed, collapse = ", "), ".")
  }
  scan_control
}

.spdest_scan <- function(model, ...) UseMethod(".spdest_scan")
.spdest_scan.default <- function(model, ...) {
  stop("model must be returned by mgcvST.set() or inlaST.set().")
}
.spdest_estimate <- function(scan, ...) UseMethod(".spdest_estimate")
.spdest_test <- function(fit, ...) UseMethod(".spdest_test")

#' Stage 1 spatial score scan
#'
#' Tests every feature for a spatial SPDE variance component against its null
#' model. The model from [mgcvST.set()] or [inlaST.set()] carries the
#' implementation and every fixed setting; Stage 1 is the same for both: GLM
#' nulls by maximum likelihood (NB2, or Poisson when the overdispersion score
#' `sum((y - mu)^2 - y) / 2` is not positive) and the score test of the raw
#' sparse SPDE field with the null projection `P0`, calibrated by the
#' saddlepoint approximation ([spde.scan()]). The null is a GLM: covariates
#' enter as linear terms or factors (0-1 columns of the model matrix); a null
#' with smooth terms or random effects is refused.
#'
#' @param Y Feature-by-observation response matrix.
#' @param model Result of [mgcvST.set()] or [inlaST.set()].
#' @param offset `NULL`, a length-n vector, or a feature-by-observation matrix
#'   added to the model offset. It is kept for [spdeST.estimate()].
#' @param feature_id Feature identifiers.
#' @param adjust,q.value Multiple-testing adjustment of the Stage 1 p-values
#'   (`"BY"`, `"BH"`, `"Sidak"` or `"none"`) and the q-value that defines the
#'   discoveries of [spdeST.estimate()].
#' @param BPPARAM BiocParallel backend; `NULL` uses the implementation default
#'   (`bpparam()` for mgcv, serial for INLA).
#' @param chunk_size Features per task; `NULL` chooses it from `BPPARAM`.
#' @param threads OpenMP threads per worker of the compiled INLA kernels; the
#'   mgcv implementation runs one thread per worker.
#' @param checkpoint_dir,resume Chunk checkpoints and their reuse.
#' @param scan_control Named list of [spde.scan()] controls, for example
#'   `list(tol = 1e-8, maxit_outer = 10)`.
#' @return An `spdeST_scan` object: the Stage 1 fit with
#'   `diagnostics$marginal_p_value`, `diagnostics$marginal_q_value`,
#'   `diagnostics$family_used` and `stage1`, carrying the model for
#'   [spdeST.estimate()].
#' @export
spdeST.scan <- function(Y, model, offset = NULL, feature_id = rownames(Y),
                        adjust = c("BY", "BH", "Sidak", "none"), q.value = 0.05,
                        BPPARAM = NULL, chunk_size = NULL, threads = 1L,
                        checkpoint_dir = NULL, resume = TRUE,
                        scan_control = list()) {
  fit <- .spdest_scan(model, Y = Y, offset = offset, feature_id = feature_id,
                      adjust = match.arg(adjust), q.value = q.value,
                      BPPARAM = BPPARAM, chunk_size = chunk_size,
                      threads = threads, checkpoint_dir = checkpoint_dir,
                      resume = resume, scan_control = scan_control)
  class(fit) <- c("spdeST_scan", class(fit))
  fit
}

#' Stage 2 spatial estimation of the selected features
#'
#' Fits the spatial model of the features chosen by `features` with the
#' implementation of the model behind `scan` (`mgcv::bam(method = "fREML",
#' discrete = TRUE)` or INLA), with the settings fixed in [mgcvST.set()] or
#' [inlaST.set()] and the offset given to [spdeST.scan()]. Stage 1 is not
#' recomputed. The fit keeps the score vectors and the per-feature materials
#' that [spdeST.test()] needs.
#'
#' With an [inlaST.set()] model, separate R processes can give slightly
#' different Stage 2 results for the same input even with one thread: on a
#' feature whose marginal likelihood is flat near its mode, the INLA
#' optimiser can stop at one of two nearby points (in one MISO run, a
#' relative difference of about 1.5e-3 in the smoothing parameter and 9e-4 in
#' the log marginal likelihood, moving log p-values of tests that include the
#' feature by about 0.09). Other features agree to about 1e-5. Repeated fits
#' in one process are identical. Tightening `control.inla$tolerance` does not
#' remove it; treat it as the precision of the INLA fit.
#'
#' @param scan Result of [spdeST.scan()], or of an earlier `spdeST.estimate()`
#'   to add features.
#' @param Y The response matrix given to [spdeST.scan()].
#' @param features `"discoveries"` (Stage 1 q-value at most `q.value`),
#'   `"all"`, `"none"`, feature IDs, one-based indices, or one logical per
#'   feature.
#' @param adjust,q.value `NULL` keeps those of [spdeST.scan()].
#' @inheritParams spdeST.scan
#' @return An `spdeST_fit` object for [spdeST.test()].
#' @export
spdeST.estimate <- function(scan, Y, features = "discoveries", adjust = NULL,
                            q.value = NULL, BPPARAM = NULL, chunk_size = NULL,
                            threads = 1L, checkpoint_dir = NULL, resume = TRUE) {
  if (!inherits(scan, "spdeST_scan")) {
    stop("scan must be returned by spdeST.scan().")
  }
  fit <- .spdest_estimate(scan, Y = Y, features = features, adjust = adjust,
                          q.value = q.value, BPPARAM = BPPARAM,
                          chunk_size = chunk_size, threads = threads,
                          checkpoint_dir = checkpoint_dir, resume = resume)
  class(fit) <- unique(c("spdeST_fit", "spdeST_scan", class(fit)))
  fit
}

#' Pairwise spatial covariance test
#'
#' Tests the spatial covariance of feature pairs estimated by
#' [spdeST.estimate()], with the implementation of their model. The result
#' keeps the Haseman-Elston materials: `t1` per pair and, in `features`,
#' `trace_M` and `trace_M2` per feature.
#'
#' @param fit Result of [spdeST.estimate()].
#' @param pairs,q.value,adjust,threads,chunk_size,checkpoint_dir,resume,verbose,moments,rank,n_per_cell,seed,k,screen
#'   Pair test settings, as described in [mgcvST.test()].
#' @return An `mgcvST_test` object.
#' @export
spdeST.test <- function(fit, pairs = NULL, q.value = 0.05,
                        adjust = c("BY", "BH", "Sidak", "none"),
                        threads = NULL, chunk_size = NULL, checkpoint_dir = NULL,
                        resume = TRUE, verbose = FALSE, moments,
                        rank = .mgcvst_pca_defaults$rank,
                        n_per_cell = .mgcvst_pca_defaults$n_per_cell,
                        seed = .mgcvst_pca_defaults$seed, k = NULL,
                        screen = FALSE) {
  if (!inherits(fit, "spdeST_fit")) {
    stop("fit must be returned by spdeST.estimate().")
  }
  args <- list(fit, pairs = pairs, q.value = q.value, adjust = match.arg(adjust),
               threads = threads, chunk_size = chunk_size,
               checkpoint_dir = checkpoint_dir, resume = resume,
               verbose = verbose, rank = rank, n_per_cell = n_per_cell,
               seed = seed, k = k, screen = screen)
  if (!missing(moments)) args$moments <- moments
  do.call(.spdest_test, args)
}
