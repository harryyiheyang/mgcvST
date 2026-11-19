.mgcvst_thread_vars <- c(
  "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
  "BLIS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "RCPP_PARALLEL_NUM_THREADS"
)

# Set process-level numerical libraries and optional data.table to one thread.
# Returns (invisibly) the settings found, for .mgcvst_thread_restore(). A
# worker process keeps the limit for its life; a user-facing function restores
# the caller's settings when it returns:
#   old <- .mgcvst_thread_limit(); on.exit(.mgcvst_thread_restore(old), add = TRUE)
.mgcvst_thread_limit <- function() {
  old <- list(
    pid = Sys.getpid(),
    env = Sys.getenv(.mgcvst_thread_vars, unset = NA_character_, names = TRUE),
    blas = tryCatch(RhpcBLASctl::blas_get_num_procs(), error = function(e) NA_integer_),
    dt = if (requireNamespace("data.table", quietly = TRUE)) {
      data.table::getDTthreads()
    } else NA_integer_
  )
  do.call(Sys.setenv, as.list(stats::setNames(
    rep("1", length(.mgcvst_thread_vars)), .mgcvst_thread_vars)))
  RhpcBLASctl::blas_set_num_threads(1L)
  if (!is.na(old$dt)) data.table::setDTthreads(1L)
  invisible(old)
}

# Put back the thread settings that .mgcvst_thread_limit() found, in the
# process that recorded them (never in a forked worker).
.mgcvst_thread_restore <- function(old) {
  if (!is.list(old) || !identical(old$pid, Sys.getpid())) return(invisible(NULL))
  unset <- is.na(old$env)
  if (any(unset)) Sys.unsetenv(names(old$env)[unset])
  if (any(!unset)) do.call(Sys.setenv, as.list(old$env[!unset]))
  if (length(old$blas) == 1L && is.finite(old$blas) && old$blas >= 1L) {
    RhpcBLASctl::blas_set_num_threads(as.integer(old$blas))
  }
  if (length(old$dt) == 1L && !is.na(old$dt)) data.table::setDTthreads(as.integer(old$dt))
  invisible(NULL)
}

# Load the package and limit the numerical libraries to one thread in a worker.
.mgcvst_worker_initialize <- function() {
  loadNamespace("spdeST")
  .mgcvst_thread_limit()
  invisible(NULL)
}

# Recover stable training-row identifiers or use their stored positions.
.mgcvst_row_id <- function(fit, n) {
  id <- NULL
  if (!is.null(fit$model)) id <- rownames(fit$model)
  if (is.null(id) || length(id) != n) id <- names(fit$residuals)
  if (is.null(id) || length(id) != n) id <- as.character(seq_len(n))
  as.character(id)
}

# Identify features with a complete compact working model.
.mgcvst_feature_available <- function(fit) {
  n <- length(fit$feature_id)
  if (identical(fit$estimator, "INLA") && identical(fit$score_backend, "sparse")) {
    if (length(fit$dispersion) != n || length(fit$lambda) != n ||
        ncol(fit$score_a) != n) {
      stop("The compact fit dimensions are incompatible with feature_id.")
    }
    return(is.finite(fit$dispersion) & fit$dispersion > 0 &
      is.finite(fit$lambda) & fit$lambda > 0 &
      colSums(!is.finite(fit$score_a)) == 0L)
  }
  if (length(fit$dispersion) != n || length(fit$lambda) != n ||
      ncol(fit$working_error) != n || ncol(fit$working_variance) != n) {
    stop("The compact fit dimensions are incompatible with feature_id.")
  }
  is.finite(fit$dispersion) & fit$dispersion > 0 &
    is.finite(fit$lambda) & fit$lambda > 0 &
    colSums(!is.finite(fit$working_error)) == 0L &
    colSums(!is.finite(fit$working_variance)) == 0L
}

# Convert a condition to serializable class, message, and call fields.
.mgcvst_condition <- function(e) {
  call <- conditionCall(e)
  list(
    class = paste(class(e), collapse = "/"),
    message = conditionMessage(e),
    call = if (is.null(call)) "" else paste(deparse(call), collapse = " ")
  )
}

# The Stage 1 score test of a null GAM, calibrated by Davies (the saddlepoint
# when Davies fails).
.mgcvst_marginal_score <- function(fit, setup) {
  score <- .mgcvst_null_score_test(fit, setup)
  p_value <- as.numeric(score$smooth.pvalue)
  if (length(p_value) != 1L || !is.finite(p_value) ||
      p_value < 0 || p_value > 1) {
    stop("The corrected marginal spatial score returned an invalid p-value.")
  }
  used <- score$method
  if (!is.character(used) || length(used) != 1L) used <- NA_character_
  list(p_value = p_value, requested_method = "davies", method = used,
       fallback = if (is.na(used)) NA else !identical(used, "davies"))
}

# Clone the minimal package function closure needed on remote workers: the
# model fit chunk and everything it calls.
.mgcvst_worker_bundle <- function() {
  names <- c(
    ".mgcvst_thread_limit", ".mgcvst_thread_vars", ".mgcvst_worker_initialize", ".mgcvst_row_id",
    ".mgcvst_condition", ".mgcvst_null_chunk", ".mgcvst_spatial_chunk",
    ".mgcvst_model_fit_one", ".mgcvst_chunk_save",
    ".mgcvst_marginal_score", ".mgcvst_marginal_working",
    ".mgcvst_marginal_matrixsqrt", ".mgcvst_marginal_saddlepoint",
    ".mgcvst_marginal_davies", ".mgcvst_null_score_setup",
    ".mgcvst_null_score_spectrum", ".mgcvst_null_score_test",
    ".mgcvst_fit_null", ".working_family_id",
    ".gam_training_lpmatrix", ".mgcvst_expand_penalty",
    ".mgcvst_model_sp", ".mgcvst_cached_model_geometry",
    ".mgcvst_nuisance_state", "rkhs_extract_working_model"
  )
  source_env <- environment(.mgcvst_worker_bundle)
  bundle <- new.env(parent = baseenv())
  for (name in names) {
    value <- get(name, envir = source_env, inherits = TRUE)
    if (is.function(value)) environment(value) <- bundle
    assign(name, value, envir = bundle)
  }
  bundle
}

#' Stage 1 of the mgcv implementation
#'
#' The method behind [spdeST.scan()] for a model from [mgcvST.set()]: the null
#' fit and Stage 1 score test of every feature. When the null has no smooth
#' term, the family is `nb()` with an estimated size and the tested smooth is
#' an spde smooth, Stage 1 is [spde.scan()] on the raw sparse field of the
#' smooth's basis; otherwise every feature gets its null GAM
#' (`mgcv::bam(method = "fREML", discrete = TRUE)`, or `mgcv::gam(method =
#' "REML")` for a null with one parametric coefficient and no smooth) and the
#' score test of the target block, calibrated by Davies with the saddlepoint
#' when Davies fails. `fit$stage1$null_engine` records `"spde.scan"` or
#' `"gam"`; `options(mgcvST.scan_null = FALSE)` keeps the GAM nulls. With a
#' negative-binomial family a feature whose Poisson overdispersion score
#' `sum((y - mu)^2 - y) / 2` is not positive is fitted as Poisson. No spatial
#' model is fitted here: [mgcvST.estimate_spatial()] fits Stage 2.
#'
#' With `checkpoint_dir`, every completed chunk is saved by the worker that
#' computed it, and a repeated call resumes from the saved chunks. Chunks are
#' keyed by their step, features, responses and family routing; pass
#' `chunk_size` explicitly for a resumable run. A checkpoint directory written
#' for another model, offset or control is refused. Every worker sets the
#' BLAS/OpenMP thread controls to one, and `control` is forced to
#' `nthreads = 1` and `ncv.threads = 1`.
#'
#' @param Y Numeric feature-by-observation matrix.
#' @param G A model returned by [mgcvST.set()].
#' @param feature_id Unique feature identifiers. Defaults to `rownames(Y)` or
#'   sequential identifiers.
#' @param BPPARAM A `BiocParallelParam`; defaults to the registered `bpparam()`.
#' @param chunk_size Positive integer number of features per task, or `NULL`
#'   (at most one chunk per worker; at most 50 features with a checkpoint).
#' @param retain_smooth Retain the score-component coefficients of the spatial
#'   fits of Stage 2.
#' @param control An `mgcv::gam.control()` object.
#' @param ... Additional arguments passed to `mgcv::bam()`.
#' @param diagnostics Compute `summary.gam()` diagnostics in Stage 2.
#' @param offset Optional additional log-scale offset: an observation-length
#'   vector or a feature-by-observation matrix in the order of `Y`. A matrix
#'   is kept on the fit for Stage 2.
#' @param adjust Multiple-testing adjustment of the Stage 1 p-values:
#'   `"BY"` (the default), `"BH"`, `"Sidak"` or `"none"`.
#' @param q.value Stage 1 discovery threshold in `(0, 1]`.
#' @param checkpoint_dir Optional directory for the resumable chunk checkpoints.
#' @param resume Reuse the completed chunks of a compatible checkpoint.
#' @return A compact `mgcvST_model_fit` (also `mgcvST_fit`) with the Stage 1
#'   p-values and q-values in `diagnostics`, `y_digest` of each response row
#'   and the estimation context that [mgcvST.estimate_spatial()] reuses.
#' @keywords internal
mgcvST.estimate <- function(
    Y, G, feature_id = rownames(Y),
    BPPARAM = BiocParallel::bpparam(), chunk_size = NULL,
    retain_smooth = FALSE,
    control = mgcv::gam.control(nthreads = 1L), ...,
    diagnostics = FALSE, offset = NULL,
    adjust = c("BY", "BH", "Sidak", "none"),
    q.value = 0.05, checkpoint_dir = NULL, resume = TRUE) {
  call <- match.call()
  adjust <- match.arg(adjust)
  q.value <- .mgcvst_check_q_value(q.value)
  # One thread in this process before any worker is forked, so that a forked
  # worker never inherits a threaded BLAS; the caller's settings return on exit.
  old_threads <- .mgcvst_thread_limit()
  on.exit(.mgcvst_thread_restore(old_threads), add = TRUE)
  .spdest_flag(diagnostics, "diagnostics")
  .spdest_flag(retain_smooth, "retain_smooth")
  if (!inherits(G, "mgcvST_model")) {
    stop("G must be a model prepared by mgcvST.set().")
  }
  .mgcvst_estimate_model(
    Y = Y, model = G, feature_id = feature_id, BPPARAM = BPPARAM,
    chunk_size = chunk_size, retain_smooth = retain_smooth, control = control,
    gam_args = list(...), call = call,
    diagnostics = diagnostics, offset = offset,
    adjust = adjust, q.value = q.value,
    checkpoint_dir = checkpoint_dir, resume = resume
  )
}

#' Print compact feature-fit diagnostics
#'
#' @param x An `mgcvST_fit` object.
#' @param ... Unused.
#' @return `x`, invisibly.
#' @export
print.mgcvST_fit <- function(x, ...) {
  cat("Compact mgcvST feature fit\n")
  cat("  features:", length(x$feature_id), "\n")
  cat("  fitted:", sum(.mgcvst_feature_available(x)), "\n")
  if (!is.null(x$diagnostics$spatial_selected)) {
    cat("  spatial models:", sum(x$diagnostics$spatial_fitted), "of",
        length(x$feature_id), "features (",
        sum(x$diagnostics$spatial_selected), "selected )\n")
  }
  cat("  backend:", x$timing$backend, "with", x$timing$workers, "worker(s)\n")
  if (!is.null(x$kappa_unit)) {
    cat("  kappa (unit scale, fixed):", format(x$kappa_unit), "\n")
    cat("  unit length L:", format(x$unit_length), "\n")
  }
  cat("  elapsed seconds:", format(x$timing$elapsed), "\n")
  cat("  object size:", format(utils::object.size(x), units = "auto"), "\n")
  invisible(x)
}

# Normalize feature-ID or feature-index pairs to integer indices i < j.
.mgcvst_pair_index <- function(pairs, feature_id) {
  pairs <- as.matrix(pairs)
  if (length(dim(pairs)) != 2L || ncol(pairs) != 2L || !nrow(pairs)) {
    stop("pairs must be a non-empty two-column matrix or data frame.")
  }
  if (is.numeric(pairs)) {
    if (any(!is.finite(pairs)) || any(pairs != floor(pairs)) ||
        any(pairs < 1L) || any(pairs > length(feature_id))) {
      stop("Numeric pairs must contain valid integer feature indices.")
    }
    index <- matrix(as.integer(pairs), ncol = 2L)
  } else {
    pairs <- matrix(as.character(pairs), ncol = 2L)
    index <- matrix(match(pairs, feature_id), ncol = 2L)
    if (anyNA(index)) {
      missing <- unique(pairs[is.na(index)])
      stop("pairs contains unknown feature IDs: ", paste(missing, collapse = ", "))
    }
  }
  if (any(index[, 1L] == index[, 2L])) {
    stop("Each requested pair must contain two different features.")
  }
  index <- cbind(
    pmin(index[, 1L], index[, 2L]),
    pmax(index[, 1L], index[, 2L])
  )
  key <- (index[, 1L] - 1) * length(feature_id) + index[, 2L]
  if (anyDuplicated(key)) {
    stop("pairs contains duplicated tests, including reversed duplicates.")
  }
  colnames(index) <- c("feature1", "feature2")
  index
}

