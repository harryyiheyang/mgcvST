# Set process-level numerical libraries and optional data.table to one thread.
.mgcvst_thread_limit <- function() {
  Sys.setenv(
    OMP_NUM_THREADS = "1",
    OPENBLAS_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1",
    BLIS_NUM_THREADS = "1",
    VECLIB_MAXIMUM_THREADS = "1",
    RCPP_PARALLEL_NUM_THREADS = "1"
  )
  RhpcBLASctl::blas_set_num_threads(1L)
  if (requireNamespace("data.table", quietly = TRUE)) {
    data.table::setDTthreads(1L)
  }
  invisible(NULL)
}

# Expand ordered source-file and package-directory inputs to normalized R files.
.mgcvst_source_files <- function(source_files) {
  if (is.null(source_files)) return(character())
  source_files <- as.character(source_files)
  if (anyNA(source_files) || any(!nzchar(source_files))) {
    stop("source_files must contain non-empty paths.")
  }

  out <- character()
  for (path in source_files) {
    if (dir.exists(path)) {
      rdir <- file.path(path, "R")
      if (dir.exists(rdir)) path <- rdir
      files <- list.files(
        path, pattern = "\\.[Rr]$", full.names = TRUE
      )
      if (!length(files)) {
        stop("No R source files were found under: ", path)
      }
      out <- c(out, sort(files))
    } else {
      if (!file.exists(path)) stop("Source file does not exist: ", path)
      out <- c(out, path)
    }
  }
  unique(normalizePath(out, winslash = "/", mustWork = TRUE))
}

# Build a stable worker-initialization key from source state and initializer code.
.mgcvst_init_key <- function(source_files, worker_init) {
  info <- if (length(source_files)) file.info(source_files) else NULL
  source_key <- if (length(source_files)) paste(
    source_files, info$size, as.numeric(info$mtime), collapse = "|"
  ) else "no-source-files"
  init_key <- if (is.function(worker_init)) {
    paste(deparse(body(worker_init)), collapse = "")
  } else {
    "no-worker-init"
  }
  paste("mgcvST", source_key, init_key, sep = "::")
}

# Apply thread limits and source custom methods once for each worker/key pair.
.mgcvst_worker_initialize <- function(source_files, worker_init, init_key) {
  loadNamespace("mgcvST")
  .mgcvst_thread_limit()
  state <- getOption("mgcvST.worker_initialization")
  if (!is.list(state) || !identical(state$pid, Sys.getpid()) ||
      !is.environment(state$guard)) {
    state <- list(pid = Sys.getpid(), guard = new.env(parent = emptyenv()))
    options(mgcvST.worker_initialization = state)
  }
  guard <- state$guard
  if (!exists(init_key, envir = guard, inherits = FALSE)) {
    for (path in source_files) sys.source(path, envir = .GlobalEnv)
    if (is.function(worker_init)) worker_init()
    assign(init_key, TRUE, envir = guard)
  }
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

# Run the requested corrected marginal score calibration.
.mgcvst_marginal_score <- function(fit, marginal_test, marginal_args,
                                   test_component = 1L, setup = NULL) {
  cacheable <- is.null(marginal_test) && is.null(marginal_args$lpmatrix)
  if (is.null(marginal_test) && !is.null(setup)) {
    score <- do.call(
      .mgcvst_null_score_test,
      utils::modifyList(
        list(null_fit = fit, setup = setup),
        marginal_args
      )
    )
  } else {
    if (is.null(marginal_test)) marginal_test <- taps_score_test
    if (!is.function(marginal_test)) {
      stop("marginal_test must be NULL or a function.")
    }
    args <- utils::modifyList(
      list(
        fit = fit, test.component = test_component, n_threads = 1L
      ),
      marginal_args
    )
    score <- do.call(marginal_test, args)
  }
  p_value <- as.numeric(score$smooth.pvalue)
  if (length(p_value) != 1L || !is.finite(p_value) ||
      p_value < 0 || p_value > 1) {
    stop("The corrected marginal spatial score returned an invalid p-value.")
  }
  requested <- marginal_args$method
  if (is.null(requested) && is.null(marginal_test) && !is.null(setup)) {
    requested <- "davies"
  }
  if (is.null(requested)) requested <- formals(marginal_test)[["method"]]
  if (!is.character(requested) || length(requested) != 1L) requested <- NA_character_
  used <- score$method
  if (!is.character(used) || length(used) != 1L) used <- NA_character_
  fallback <- if (is.na(requested) || is.na(used)) NA else
    identical(requested, "davies") && !identical(used, "davies")
  cache <- if (cacheable) attr(score, "marginal_spectrum", exact = TRUE) else NULL
  list(p_value = p_value, requested_method = requested,
       method = used, fallback = fallback, cache = cache)
}

# Clone the minimal package function closure needed on remote workers: the
# model fit chunk and everything it calls.
.mgcvst_worker_bundle <- function() {
  names <- c(
    ".mgcvst_thread_limit", ".mgcvst_worker_initialize", ".mgcvst_row_id",
    ".mgcvst_condition", ".mgcvst_null_chunk", ".mgcvst_spatial_chunk",
    ".mgcvst_model_fit_one", ".mgcvst_chunk_save",
    ".mgcvst_marginal_score", ".mgcvst_marginal_geometry",
    ".mgcvst_marginal_spectrum", ".mgcvst_marginal_working",
    ".mgcvst_marginal_matrixsqrt", ".mgcvst_marginal_saddlepoint",
    ".mgcvst_marginal_davies", ".mgcvst_null_score_setup",
    ".mgcvst_null_score_spectrum", ".mgcvst_null_score_test",
    ".mgcvst_fit_null", ".working_family_id", "taps_score_test",
    ".gam_training_lpmatrix", ".mgcvst_expand_penalty",
    ".mgcvst_model_geometry", ".mgcvst_geometry_signature", ".mgcvst_model_sp",
    ".mgcvst_cached_model_geometry", ".mgcvst_training_design",
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

#' Estimate compact covariance working summaries
#'
#' Fits each row of `Y` from a reusable shared design and retains the fixed
#' numerical summaries needed by [mgcvST.test()]. `G` is a model prepared by
#' [mgcvST.set()] or [model.set()], or a reusable
#' `mgcv::gam(..., fit = FALSE)` setup, which is converted as by
#' `mgcvST.set(G = G)`; both enter the same estimation and testing path.
#'
#' Estimation has two steps. Step 1 fits the null model with the spatial score
#' smooth removed for every feature; the marginal screening uses that null
#' PIRLS state, the prepared `G$X`, and the target penalty, and the Stage 1
#' null-first p-values are adjusted by `adjust` into the Stage 1 q-values
#' `diagnostics$marginal_q_value`. No null fit keeps a working vector. Step 2
#' fits the spatial model of the features chosen by `spatial` only. Features
#' without a spatial model are marked `spatial_fitted = FALSE` in the
#' diagnostics, have missing working errors, variances, `dispersion` and
#' `lambda`, and are unavailable to [mgcvST.test()] and [mgcvST.wgcna()];
#' `pairs = NULL` in [mgcvST.test()] means all pairs among the features with a
#' spatial model. [mgcvST.estimate_spatial()] adds spatial models for more
#' features after step 1.
#'
#' With `checkpoint_dir`, every completed chunk of either step is saved by the
#' worker that computed it, and a repeated call with the same arguments resumes
#' from the saved chunks. Chunks are keyed by their step, features, responses
#' and family routing, so pass `chunk_size` explicitly for a resumable run:
#' without it the chunks depend on the number of workers (and step 2 on the
#' number of selected features), and a changed chunking recomputes rather than
#' reuses. A chunk whose result holds a failed feature is computed again on
#' resume. A checkpoint directory written for another model, offset or control
#' is refused. Without `chunk_size`, a run with a checkpoint directory uses
#' chunks of at most 50 features.
#'
#' Wood diagnostics are optional (`diagnostics = TRUE`).
#' The marginal test is the package-local `taps_score_test()`, using the TAPS
#' arithmetic included in mgcvST. No external mgcv.taps installation or sourced
#' score function is required. Full `gam` objects are never
#' retained. Genes are processed in chunks through `BiocParallel`; the default
#' uses the registered backend. On Windows, use a persistent
#' `SnowParam(type = "SOCK")`. [mgcvST.test()] runs in the calling process with
#' OpenMP `threads` and takes no `BiocParallel` backend.
#' A fitting, compaction, or marginal-test error is recorded for that feature;
#' the other features continue. A marginal-test error leaves the already
#' constructed compact working model intact and only its marginal p-value
#' unavailable.
#' The built-in marginal test is calibrated by Davies. When Davies errors,
#' returns a missing or non-finite p-value, or returns a value outside
#' (0, 1], the saddlepoint approximation (Kuonen 1999) is used instead.
#' Marginal diagnostics retain the requested method, the method actually used
#' (`"davies"` or `"saddlepoint"`), and a Davies-to-saddlepoint fallback flag.
#' Custom callbacks that omit method metadata leave the corresponding
#' diagnostics as `NA`.
#'
#' `source_files` supports source-first custom smooths. Each SOCK worker
#' sources the ordered files once, before it evaluates a chunk. Multicore
#' children inherit the parent state, while the same guard remains harmless.
#' Every worker sets common BLAS/OpenMP thread controls to one, and `control`
#' is forced to `nthreads = 1` and `ncv.threads = 1`.
#'
#' Keep `SerialParam()` for reproducible small jobs. On Windows, the standard
#' parallel backend is a persistent `SnowParam(type = "SOCK")`. On a
#' single-node Linux Slurm allocation, prefer `MulticoreParam()` when forking
#' is available; do not force SOCK workers on Linux. Persistent SOCK workers
#' are useful when fork is unavailable and their startup cost is amortized
#' across a large scan.
#'
#' @param Y Numeric feature-by-observation matrix.
#' @param G A model returned by [mgcvST.set()] or [model.set()], or a reusable
#'   setup returned by `mgcv::gam(..., fit = FALSE)` with one marked SPDE
#'   smooth (Gaussian, negative-binomial, Poisson or quasi-Poisson family). Its
#'   response is replaced by each row of `Y`.
#' @param feature_id Unique feature identifiers. Defaults to `rownames(Y)` or
#'   sequential identifiers.
#' @param BPPARAM A `BiocParallelParam`; defaults to the registered `bpparam()`.
#' @param chunk_size Positive integer number of genes per task, or `NULL`. The
#'   default creates at most one chunk per worker, limiting repeated
#'   serialization on SOCK workers. A chunk payload is built when a worker takes
#'   the chunk, so at most about one chunk of `Y` per worker is copied.
#' @param source_files Ordered R files, or directories containing R files, to
#'   source once per worker. Use this for sourced custom S3 smooth methods.
#' @param worker_init Optional zero-argument initialization function run once
#'   per worker after `source_files`.
#' @param marginal_test Optional function implementing the corrected marginal
#'   spatial score interface. `NULL` uses mgcvST's package-local implementation.
#'   A custom function is used only when explicitly supplied.
#' @param offset Optional additional log/link-scale offset for a model prepared
#'   by [mgcvST.set()] (or a `gam(fit = FALSE)` setup passed as `G`): a shared
#'   observation-length vector or a feature-by-observation matrix in exactly
#'   the same order as `Y`. Added to the shared formula/setup offset. A matrix
#'   is kept on the fit (in the total offset and for a later step 2).
#'   Covariates and smooths are fixed by `mgcvST.set()` and cannot vary by
#'   gene.
#' @param diagnostics Logical; compute `summary.gam()`/Wood diagnostics.
#'   FALSE (default) leaves Wood fields NA without calling summary.
#'   Basic convergence and fitting diagnostics are still retained.
#' @param retain_marginal Logical; retain minimal null-fit inputs for a later
#'   `mgcvST.marginal()` recalibration call. FALSE by default. Marginal testing
#'   during estimation is always performed; full gam objects are never retained.
#' @param marginal_args Named list of additional marginal-score arguments.
#'   `fit`, `test.component`, and `n_threads` are controlled by mgcvST.
#'   The built-in test uses Davies with the saddlepoint fallback described in
#'   Details, typically needed in the extreme upper tail; its `method` accepts
#'   only `"davies"`, and `max_eps` and `max_iter` set the Davies accuracy and
#'   integration limit. No additional marginal call is needed.
#' @param retain_smooth Logical; retain the feature-by-coefficient smooth
#'   coefficient matrix of the score component. This opt-in representation
#'   supports prediction and other downstream uses without retaining full
#'   `gam` objects.
#' @param method Retained for API compatibility. Null and full fits use
#'   `mgcv::bam(method = "fREML", discrete = TRUE)`, except a null model with
#'   one parametric coefficient and no smooths uses `mgcv::gam(method = "REML")`
#'   to avoid BAM's one-column QR dimension error.
#' @param control An `mgcv::gam.control()` object. Internal thread counts are
#'   always forced to one. It may additionally carry `poisson_screen_phi`
#'   (default `1.01`), the Poisson prescreen threshold: with a negative-binomial
#'   family, each feature first gets an offset-and-covariate-only Poisson GLM,
#'   and a feature whose Pearson dispersion
#'   `phi = sum((y - mu)^2 / mu) / (n - p)` is at most the threshold is fitted
#'   with `stats::quasipoisson(link = "log")` instead. The Poisson and
#'   quasipoisson point estimates agree. The null score fit and the full
#'   spatial fit both estimate their own scale and smoothing parameters. The
#'   screening phi is retained in diagnostics and is not passed as a fixed
#'   `fit$sig2`. In the INLA path, routed features use plain Poisson. Set the
#'   threshold to `0` to disable the screen (`NULL` restores the default);
#'   other families ignore it. The entry is removed before `control` reaches
#'   `mgcv::bam()`. Per-feature `prescreen_phi` and `family_used` are reported
#'   in the diagnostics, and `family_used` reads `"quasipoisson"` for a routed
#'   gene.
#' @param ... Additional arguments passed to `mgcv::bam()`.
#' @param spatial Features that receive a spatial model in step 2:
#'   `"discoveries"` (the default; the features with Stage 1 q-value at most
#'   `q.value`), `"all"`, `"none"`, a vector of feature IDs or one-based
#'   indices, or a logical vector with one value per feature. A user-given
#'   vector allows, for example, a Stage 1 adjustment within a modality or a
#'   family of features that the user performs outside this function.
#' @param adjust Multiple-testing adjustment of the Stage 1 p-values:
#'   `"BY"` (the default), `"BH"`, `"Sidak"` or `"none"`.
#' @param q.value Stage 1 discovery threshold in `(0, 1]`.
#' @param checkpoint_dir Optional directory for the resumable chunk checkpoints
#'   of both steps.
#' @param resume Reuse the completed chunks of a compatible checkpoint.
#' @return A compact object of class `mgcvST_model_fit` (also `mgcvST_fit`)
#'   containing marginal score p-values, feature IDs, working errors and
#'   variances, separate per-feature `dispersion` and `lambda`, shared score
#'   geometry, optimized fit criterion and its original mgcv name, exact
#'   residual degrees of freedom, timing, and convergence diagnostics. One
#'   shared `geometry$nuisance_design` and one small conditional nuisance
#'   block per feature in `nuisance_covariance` are retained; full GAM and
#'   `Vp` objects are discarded. When `retain_smooth = TRUE`, it also contains
#'   `smooth_coefficients`. Score methods derive the field scale as
#'   `dispersion / lambda`. The diagnostics table holds `marginal_q_value`,
#'   `spatial_selected` and `spatial_fitted`, the effective degrees of freedom
#'   of the spatial smooth `edf_spatial` and the flag `spatial_degenerate`
#'   (`edf_spatial` below 3: the gene is fitted with essentially no spatial
#'   field, and [mgcvST.test()] gives every pair with such a gene p = 1 and
#'   status 4), and `y_digest` records a digest of each response row, which
#'   [mgcvST.estimate_spatial()] checks.
#' @seealso [mgcvST.estimate_spatial()] to add spatial models after step 1.
#' @export
mgcvST.estimate <- function(
    Y, G, feature_id = rownames(Y),
    BPPARAM = BiocParallel::bpparam(), chunk_size = NULL,
    source_files = NULL, worker_init = NULL,
    marginal_test = NULL, marginal_args = list(), method = "REML",
    retain_smooth = FALSE,
    control = mgcv::gam.control(nthreads = 1L), ...,
    diagnostics = FALSE, retain_marginal = FALSE, offset = NULL,
    spatial = "discoveries", adjust = c("BY", "BH", "Sidak", "none"),
    q.value = 0.05, checkpoint_dir = NULL, resume = TRUE) {
  if ("marginal" %in% names(list(...))) {
    stop("mgcvST.estimate() always runs the marginal score test; remove marginal.")
  }
  call <- match.call()
  adjust <- match.arg(adjust)
  q.value <- .mgcvst_check_q_value(q.value)
  for (name in c("diagnostics", "retain_marginal")) {
    value <- get(name)
    if (!is.logical(value) || length(value) != 1L || is.na(value)) {
      stop(name, " must be TRUE or FALSE.")
    }
  }
  if (!is.null(marginal_test) && !is.function(marginal_test)) {
    stop("marginal_test must be NULL or a function.")
  }
  if (!is.list(marginal_args) ||
      (length(marginal_args) && (is.null(names(marginal_args)) ||
                               any(!nzchar(names(marginal_args)))))) {
    stop("marginal_args must be a named list.")
  }
  if (any(names(marginal_args) %in% c("fit", "test.component", "n_threads"))) {
    stop("Do not supply fit, test.component or n_threads through marginal_args.")
  }
  if (is.null(marginal_test) && !is.null(marginal_args[["method"]]) &&
      !identical(marginal_args[["method"]], "davies")) {
    stop("marginal_args$method = ", paste(deparse(marginal_args[["method"]]), collapse = ""),
         " is not available: the marginal score test is calibrated by Davies, ",
         "with a saddlepoint approximation when Davies fails. ",
         "Remove method from marginal_args.", call. = FALSE)
  }
  if (!is.logical(retain_smooth) || length(retain_smooth) != 1L ||
      is.na(retain_smooth)) {
    stop("retain_smooth must be TRUE or FALSE.")
  }
  if (!is.function(worker_init) && !is.null(worker_init)) {
    stop("worker_init must be NULL or a zero-argument function.")
  }
  if (is.function(worker_init) && length(formals(worker_init))) {
    stop("worker_init must be a zero-argument function.")
  }
  if (!inherits(G, "mgcvST_model")) {
    if (!is.list(G) || is.null(G$y) || is.null(G$family) || is.null(G$smooth)) {
      stop("G must be a model prepared by mgcvST.set() or model.set(), or a ",
           "reusable setup returned by mgcv::gam(..., fit = FALSE).")
    }
    G <- .mgcvst_set_prepare(G = G, .allow_poisson = TRUE)
  }
  .mgcvst_estimate_model(
    Y = Y, model = G, feature_id = feature_id, BPPARAM = BPPARAM,
    chunk_size = chunk_size, source_files = source_files,
    worker_init = worker_init, marginal_test = marginal_test,
    marginal_args = marginal_args, method = method,
    retain_smooth = retain_smooth, control = control,
    gam_args = list(...), call = call,
    diagnostics = diagnostics, retain_marginal = retain_marginal, offset = offset,
    spatial = spatial, adjust = adjust, q.value = q.value,
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

