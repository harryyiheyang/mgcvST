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
  guard_name <- ".mgcvST_worker_initialization"
  if (!exists(guard_name, envir = .GlobalEnv, inherits = FALSE)) {
    assign(guard_name, new.env(parent = emptyenv()), envir = .GlobalEnv)
  }
  guard <- get(guard_name, envir = .GlobalEnv, inherits = FALSE)
  if (!exists(init_key, envir = guard, inherits = FALSE)) {
    for (path in source_files) sys.source(path, envir = .GlobalEnv)
    if (is.function(worker_init)) worker_init()
    assign(init_key, TRUE, envir = guard)
  }
  invisible(NULL)
}

# Locate the sole SPDE smooth supported by the high-throughput API.
.mgcvst_spde_index <- function(fit) {
  if (!inherits(fit, "gam")) stop("fit must inherit from class 'gam'.")
  idx <- which(vapply(
    fit$smooth, inherits, logical(1), what = "spde.smooth"
  ))
  if (!length(idx)) stop("The fit does not contain an SPDE smooth.")
  if (length(idx) != 1L) stop("The fit must contain exactly one SPDE smooth.")
  if (length(fit$smooth) != 1L) {
    stop(
      "The fit contains additional nuisance smooths; their covariance is not ",
      "implemented by the validated score operator."
    )
  }
  idx
}

# Locate the sole SPDE or SPDE-PC smooth supported by the pair-score API.
.mgcvst_score_smooth_index <- function(fit) {
  if (!inherits(fit, "gam")) stop("fit must inherit from class 'gam'.")
  idx <- which(vapply(
    fit$smooth,
    function(s) inherits(s, "spde.smooth") || inherits(s, "spdePC.smooth"),
    logical(1)
  ))
  if (!length(idx)) stop("The fit does not contain an SPDE or SPDE-PC smooth.")
  if (length(idx) != 1L) stop("The fit must contain exactly one SPDE or SPDE-PC smooth.")
  if (length(fit$smooth) != 1L) {
    stop(
      "The fit contains additional nuisance smooths; their covariance is not ",
      "implemented by the validated score operator."
    )
  }
  idx
}

# Recover stable training-row identifiers, with a positional fallback.
.mgcvst_row_id <- function(fit, n) {
  id <- NULL
  if (!is.null(fit$model)) id <- rownames(fit$model)
  if (is.null(id) || length(id) != n) id <- names(fit$residuals)
  if (is.null(id) || length(id) != n) id <- as.character(seq_len(n))
  as.character(id)
}

# Stop unless two compact fits share row, basis, penalty, and fixed design.
.mgcvst_geometry_equal <- function(x, y, tol) {
  if (!identical(x$row_id, y$row_id)) {
    stop("Feature fits do not share the same row identifiers and ordering.")
  }
  if (!all(dim(x$B) == dim(y$B))) {
    stop("Feature fits have different SPDE basis dimensions.")
  }
  B_scale <- max(1, max(abs(x$B)), max(abs(y$B)))
  if (max(abs(x$B - y$B)) > tol * B_scale) {
    stop("Feature fits do not share the same SPDE basis and row ordering.")
  }
  Qx <- as.matrix(x$Q)
  Qy <- as.matrix(y$Q)
  if (!all(dim(Qx) == dim(Qy))) {
    stop("Feature fits have different unscaled Q dimensions.")
  }
  Q_scale <- max(1, max(abs(Qx)), max(abs(Qy)))
  if (max(abs(Qx - Qy)) > tol * Q_scale) {
    stop("Feature fits do not share the same unscaled SPDE precision Q.")
  }
  if (!all(dim(x$X) == dim(y$X))) {
    stop("Feature fits have different parametric design dimensions.")
  }
  X_difference <- if (length(x$X)) max(abs(x$X - y$X)) else 0
  X_scale <- if (length(x$X)) {
    max(1, max(abs(x$X)), max(abs(y$X)))
  } else {
    1
  }
  if (X_difference > tol * X_scale) {
    stop("Feature fits do not share the same parametric design and row ordering.")
  }
  invisible(TRUE)
}

# Reduce one gam fit to the working summaries and shared geometry used downstream.
.mgcvst_compact_fit <- function(fit) {
  W <- rkhs_extract_working_model(fit)
  L <- .gam_training_lpmatrix(fit)
  smooth_index <- .mgcvst_spde_index(fit)
  S <- .gam_single_smooth(fit, smooth_index, L)
  smooth_cols <- unique(unlist(lapply(
    fit$smooth, function(s) seq.int(s$first.para, s$last.para)
  )))
  X <- L[, setdiff(seq_len(ncol(L)), smooth_cols), drop = FALSE]
  geometry <- list(
    B = S$B,
    Q = S$Q,
    X = X,
    row_id = .mgcvst_row_id(fit, nrow(L)),
    smooth_label = S$label,
    smooth_columns = S$columns
  )
  fit_summary <- summary(fit)
  smooth_table <- fit_summary$s.table[smooth_index, , drop = FALSE]
  wood_p_value <- as.numeric(smooth_table[1L, "p-value"])
  smooth_edf <- as.numeric(smooth_table[1L, "edf"])
  statistic_column <- ncol(smooth_table) - 1L
  wood_statistic <- as.numeric(smooth_table[1L, statistic_column])
  parameters <- W$family_parameters
  if (is.null(parameters)) parameters <- numeric()
  list(
    working_error = W$working_error,
    working_variance = W$working_variance,
    field_scale = W$dispersion / S$sp,
    smoothing_parameter = S$sp,
    dispersion = W$dispersion,
    family = W$family,
    family_parameters = paste(format(parameters, digits = 16), collapse = ","),
    smooth_label = S$label,
    wood_p_value = wood_p_value,
    wood_statistic = wood_statistic,
    smooth_edf = smooth_edf,
    geometry = geometry
  )
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

# Run the corrected marginal score with an explicit Davies-to-Liu fallback.
.mgcvst_marginal_score <- function(fit, marginal_test, marginal_args) {
  if (is.null(marginal_test)) {
    marginal_test <- get0(
      "taps_score_test", envir = .GlobalEnv, mode = "function",
      inherits = TRUE
    )
    if (is.null(marginal_test) &&
        requireNamespace("mgcv.taps", quietly = TRUE)) {
      marginal_test <- getExportedValue("mgcv.taps", "taps_score_test")
    }
    if (is.null(marginal_test)) {
      stop(
        "The corrected marginal spatial score requires taps_score_test(). ",
        "Install mgcv.taps, source it through source_files, or supply ",
        "marginal_test explicitly."
      )
    }
  }
  if (!is.function(marginal_test)) {
    stop("marginal_test must be NULL or a function.")
  }
  args <- utils::modifyList(
    list(
      fit = fit, test.component = 1L, method = "davies",
      max_eps = 1e-8, max_iter = 100000L, n_threads = 1L
    ),
    marginal_args
  )
  score <- do.call(marginal_test, args)
  p_value <- as.numeric(score$smooth.pvalue)
  used <- "davies"
  if (length(p_value) != 1L || !is.finite(p_value) ||
      p_value <= 0 || p_value > 1) {
    args$method <- "liu"
    args$n_threads <- 1L
    score <- do.call(marginal_test, args)
    p_value <- as.numeric(score$smooth.pvalue)
    used <- "liu_fallback"
  }
  if (length(p_value) != 1L || !is.finite(p_value) ||
      p_value < 0 || p_value > 1) {
    stop("The corrected marginal spatial score returned an invalid p-value.")
  }
  list(p_value = p_value, method = used)
}

# Clone the minimal package function closure needed on remote workers.
.mgcvst_worker_bundle <- function() {
  names <- c(
    ".mgcvst_thread_limit", ".mgcvst_worker_initialize",
    ".mgcvst_spde_index", ".mgcvst_row_id", ".mgcvst_geometry_equal",
    ".mgcvst_compact_fit", ".mgcvst_condition", ".mgcvst_fit_chunk",
    ".mgcvst_test_chunk", ".mgcvst_marginal_score", ".working_family_id",
    ".gam_training_lpmatrix", ".gam_single_smooth",
    "rkhs_extract_working_model", ".magic_mm", ".magic_solve",
    ".as_numeric_matrix", ".rkhs_score_operator_factor",
    "rkhs_score_operator", "rkhs_score_apply_P",
    ".score_factor", "rkhs_score_summary", "rkhs_score_information",
    ".psd_factor", "rkhs_score_singular_values", ".rkhs_score_moments",
    "rkhs_score_calibrate", ".liu_squared_score",
    ".liu_squared_score_moments", "rkhs_covariance_score"
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

# Fit and compact one feature chunk while preserving every per-feature failure.
.mgcvst_fit_chunk <- function(payload, G0, family_raw, method, control,
                              gam_args, source_files, worker_init, init_key,
                              geometry_tol, marginal_test, marginal_args) {
  .mgcvst_worker_initialize(source_files, worker_init, init_key)
  ids <- payload$index
  Y <- payload$Y
  feature_id <- payload$feature_id
  n <- ncol(Y)
  k <- nrow(Y)
  E <- matrix(NA_real_, nrow = n, ncol = k)
  V <- matrix(NA_real_, nrow = n, ncol = k)
  field_scale <- smoothing_parameter <- dispersion <- rep(NA_real_, k)
  marginal_p_value <- wood_p_value <- wood_statistic <- smooth_edf <-
    rep(NA_real_, k)
  success <- converged <- smoothing_converged <- rep(FALSE, k)
  smoothing_converged[] <- NA
  family <- family_parameters <- smooth_label <- marginal_method <-
    rep(NA_character_, k)
  irls_iterations <- rep(NA_integer_, k)
  outer_convergence <- rep(NA_character_, k)
  fit_seconds <- summary_seconds <- marginal_seconds <- rep(0, k)
  error_class <- error_message <- error_call <- rep(NA_character_, k)
  geometry <- NULL
  response_index <- attr(G0$terms, "response")
  if (length(response_index) != 1L || response_index < 1L ||
      is.null(G0$mf) || response_index > ncol(G0$mf)) {
    stop("G does not contain a reusable model-frame response column.")
  }

  G1 <- G0
  for (j in seq_len(k)) {
    response <- as.numeric(Y[j, ])
    G1$y <- response
    G1$mf[[response_index]] <- response
    G1$family <- unserialize(family_raw)

    t0 <- proc.time()[["elapsed"]]
    fit_result <- tryCatch(
      list(
        value = do.call(
          mgcv::gam,
          c(list(G = G1, method = method, control = control), gam_args)
        ),
        error = NULL
      ),
      error = function(e) list(value = NULL, error = e)
    )
    fit_seconds[j] <- proc.time()[["elapsed"]] - t0
    if (!is.null(fit_result$error)) {
      err <- .mgcvst_condition(fit_result$error)
      error_class[j] <- err$class
      error_message[j] <- err$message
      error_call[j] <- err$call
      next
    }

    fit <- fit_result$value
    converged[j] <- isTRUE(fit$converged)
    if (!is.null(fit$mgcv.conv)) {
      smoothing_converged[j] <- isTRUE(fit$mgcv.conv)
    }
    if (!is.null(fit$iter) && length(fit$iter)) {
      irls_iterations[j] <- as.integer(fit$iter[[1L]])
    }
    if (!is.null(fit$outer.info$conv)) {
      outer_convergence[j] <- paste(fit$outer.info$conv, collapse = "; ")
    }

    t0 <- proc.time()[["elapsed"]]
    compact_result <- tryCatch(
      list(value = .mgcvst_compact_fit(fit), error = NULL),
      error = function(e) list(value = NULL, error = e)
    )
    summary_seconds[j] <- proc.time()[["elapsed"]] - t0
    if (!is.null(compact_result$error)) {
      err <- .mgcvst_condition(compact_result$error)
      error_class[j] <- err$class
      error_message[j] <- err$message
      error_call[j] <- err$call
      next
    }

    compact <- compact_result$value
    geometry_result <- tryCatch(
      {
        if (is.null(geometry)) {
          geometry <- compact$geometry
        } else {
          .mgcvst_geometry_equal(geometry, compact$geometry, geometry_tol)
        }
        NULL
      },
      error = function(e) e
    )
    if (!is.null(geometry_result)) {
      err <- .mgcvst_condition(geometry_result)
      error_class[j] <- err$class
      error_message[j] <- err$message
      error_call[j] <- err$call
      next
    }

    smoothing_parameter[j] <- compact$smoothing_parameter
    dispersion[j] <- compact$dispersion
    family[j] <- compact$family
    family_parameters[j] <- compact$family_parameters
    smooth_label[j] <- compact$smooth_label
    wood_p_value[j] <- compact$wood_p_value
    wood_statistic[j] <- compact$wood_statistic
    smooth_edf[j] <- compact$smooth_edf

    t0 <- proc.time()[["elapsed"]]
    marginal_result <- tryCatch(
      list(
        value = .mgcvst_marginal_score(
          fit, marginal_test = marginal_test,
          marginal_args = marginal_args
        ),
        error = NULL
      ),
      error = function(e) list(value = NULL, error = e)
    )
    marginal_seconds[j] <- proc.time()[["elapsed"]] - t0
    if (!is.null(marginal_result$error)) {
      err <- .mgcvst_condition(marginal_result$error)
      error_class[j] <- err$class
      error_message[j] <- err$message
      error_call[j] <- err$call
      next
    }

    success[j] <- TRUE
    E[, j] <- compact$working_error
    V[, j] <- compact$working_variance
    field_scale[j] <- compact$field_scale
    marginal_p_value[j] <- marginal_result$value$p_value
    marginal_method[j] <- marginal_result$value$method
  }

  diagnostics <- data.frame(
    index = ids,
    feature_id = feature_id,
    success = success,
    family = family,
    smooth_label = smooth_label,
    smoothing_parameter = smoothing_parameter,
    dispersion = dispersion,
    field_scale = field_scale,
    marginal_p_value = marginal_p_value,
    marginal_method = marginal_method,
    wood_p_value = wood_p_value,
    wood_statistic = wood_statistic,
    smooth_edf = smooth_edf,
    family_parameters = family_parameters,
    converged = converged,
    smoothing_converged = smoothing_converged,
    irls_iterations = irls_iterations,
    outer_convergence = outer_convergence,
    fit_seconds = fit_seconds,
    summary_seconds = summary_seconds,
    marginal_seconds = marginal_seconds,
    total_seconds = fit_seconds + summary_seconds + marginal_seconds,
    pid = Sys.getpid(),
    error_class = error_class,
    error_message = error_message,
    error_call = error_call,
    stringsAsFactors = FALSE
  )
  list(
    index = ids,
    working_error = E,
    working_variance = V,
    field_scale = field_scale,
    diagnostics = diagnostics,
    geometry = geometry
  )
}

#' Estimate compact marginal-screen and covariance working summaries
#'
#' Fits each row of `Y` with a reusable `mgcv::gam(fit = FALSE)` setup and
#' computes the corrected marginal spatial score used for feature screening,
#' and retains only the fixed numerical summaries needed by [mgcvST.test()].
#' The default marginal test is `mgcv.taps::taps_score_test()`, or a sourced
#' `taps_score_test()` found on the worker. Full `gam` objects are never
#' retained. Genes are processed in chunks through `BiocParallel`; the default
#' is serial. On Windows, use a persistent `SnowParam(type = "SOCK")` object
#' and pass the same object to estimation and testing.
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
#' @param G A reusable setup returned by `mgcv::gam(..., fit = FALSE)`. Its
#'   response is replaced by each row of `Y`.
#' @param feature_id Unique feature identifiers. Defaults to `rownames(Y)` or
#'   sequential identifiers.
#' @param BPPARAM A `BiocParallelParam`; defaults to `SerialParam()`.
#' @param chunk_size Positive number of genes per task. The default creates at
#'   most one chunk per worker, limiting repeated serialization on SOCK
#'   workers.
#' @param source_files Ordered R files, or directories containing R files, to
#'   source once per worker. Use this for sourced custom S3 smooth methods.
#' @param worker_init Optional zero-argument initialization function run once
#'   per worker after `source_files`.
#' @param marginal_test Optional function implementing the corrected marginal
#'   spatial score interface of `mgcv.taps::taps_score_test()`. `NULL` locates
#'   a sourced `taps_score_test()` or the installed `mgcv.taps` export on each
#'   worker.
#' @param marginal_args Named list of additional marginal-score arguments.
#'   `fit`, `test.component`, `method`, and `n_threads` are controlled by
#'   mgcvST; Davies is primary, Liu is the explicit numerical fallback, and
#'   score threads remain one.
#' @param method Fitting method passed to `mgcv::gam()`. `"REML"` is the
#'   default.
#' @param control An `mgcv::gam.control()` object. Internal thread counts are
#'   always forced to one.
#' @param geometry_tol Positive relative tolerance used to verify shared row,
#'   basis, fixed-design, and unscaled precision geometry.
#' @param ... Additional arguments passed to `mgcv::gam(G = G, ...)`.
#' @return A compact object of class `mgcvST_fit` containing marginal score
#'   p-values, Wood audit values, feature IDs, working errors and variances,
#'   field scales, shared geometry, timing and convergence diagnostics, and an
#'   explicit failures table.
#' @export
mgcvST.estimate <- function(
    Y, G, feature_id = rownames(Y),
    BPPARAM = BiocParallel::SerialParam(), chunk_size = NULL,
    source_files = NULL, worker_init = NULL,
    marginal_test = NULL, marginal_args = list(), method = "REML",
    control = mgcv::gam.control(nthreads = 1L), geometry_tol = 1e-9, ...) {
  call <- match.call()
  Y <- as.matrix(Y)
  storage.mode(Y) <- "double"
  if (length(dim(Y)) != 2L || !nrow(Y) || !ncol(Y) || any(!is.finite(Y))) {
    stop("Y must be a non-empty finite numeric feature-by-observation matrix.")
  }
  if (!is.list(G) || is.null(G$y) || is.null(G$family) || is.null(G$smooth)) {
    stop("G must be a reusable setup returned by mgcv::gam(..., fit = FALSE).")
  }
  if (ncol(Y) != length(G$y)) {
    stop("ncol(Y) must equal the number of observations in G.")
  }
  if (is.null(feature_id)) feature_id <- as.character(seq_len(nrow(Y)))
  feature_id <- as.character(feature_id)
  if (length(feature_id) != nrow(Y) || anyNA(feature_id) ||
      any(!nzchar(feature_id)) || anyDuplicated(feature_id)) {
    stop("feature_id must contain one unique non-empty identifier per row of Y.")
  }
  if (!inherits(BPPARAM, "BiocParallelParam")) {
    stop("BPPARAM must inherit from 'BiocParallelParam'.")
  }
  if (!is.function(worker_init) && !is.null(worker_init)) {
    stop("worker_init must be NULL or a zero-argument function.")
  }
  if (is.function(worker_init) && length(formals(worker_init))) {
    stop("worker_init must be a zero-argument function.")
  }
  if (!is.null(marginal_test) && !is.function(marginal_test)) {
    stop("marginal_test must be NULL or a function.")
  }
  if (!is.list(marginal_args) ||
      (length(marginal_args) &&
       (is.null(names(marginal_args)) || any(!nzchar(names(marginal_args)))))) {
    stop("marginal_args must be a named list.")
  }
  marginal_forbidden <- intersect(
    names(marginal_args),
    c("fit", "test.component", "method", "n_threads")
  )
  if (length(marginal_forbidden)) {
    stop(
      "Do not supply these arguments through marginal_args: ",
      paste(marginal_forbidden, collapse = ", ")
    )
  }
  geometry_tol <- as.numeric(geometry_tol)
  if (length(geometry_tol) != 1L || !is.finite(geometry_tol) ||
      geometry_tol <= 0) {
    stop("geometry_tol must be one positive finite value.")
  }
  if (!is.list(control)) stop("control must be returned by mgcv::gam.control().")
  control$nthreads <- 1L
  control$ncv.threads <- 1L
  gam_args <- list(...)
  forbidden <- intersect(names(gam_args), c("G", "family", "method", "control"))
  if (length(forbidden)) {
    stop("Do not supply these arguments through ...: ", paste(forbidden, collapse = ", "))
  }

  source_files <- .mgcvst_source_files(source_files)
  init_key <- .mgcvst_init_key(source_files, worker_init)
  workers <- max(1L, min(nrow(Y), BiocParallel::bpworkers(BPPARAM)))
  if (is.null(chunk_size)) chunk_size <- ceiling(nrow(Y) / workers)
  chunk_size <- as.integer(chunk_size)
  if (length(chunk_size) != 1L || is.na(chunk_size) || chunk_size < 1L) {
    stop("chunk_size must be one positive integer.")
  }
  chunk_id <- ceiling(seq_len(nrow(Y)) / chunk_size)
  ids <- split(seq_len(nrow(Y)), chunk_id)
  payload <- lapply(ids, function(i) list(
    index = i,
    feature_id = feature_id[i],
    Y = Y[i, , drop = FALSE]
  ))
  family_raw <- serialize(G$family, NULL)
  worker_bundle <- .mgcvst_worker_bundle()
  fit_chunk <- get(".mgcvst_fit_chunk", envir = worker_bundle,
                   inherits = FALSE)

  t0 <- proc.time()[["elapsed"]]
  chunks <- BiocParallel::bplapply(
    payload, fit_chunk,
    G0 = G, family_raw = family_raw, method = method, control = control,
    gam_args = gam_args, source_files = source_files,
    worker_init = worker_init, init_key = init_key,
    geometry_tol = geometry_tol, marginal_test = marginal_test,
    marginal_args = marginal_args, BPPARAM = BPPARAM
  )
  elapsed <- proc.time()[["elapsed"]] - t0

  n <- ncol(Y)
  p <- nrow(Y)
  E <- matrix(NA_real_, nrow = n, ncol = p,
              dimnames = list(NULL, feature_id))
  V <- matrix(NA_real_, nrow = n, ncol = p,
              dimnames = list(NULL, feature_id))
  field_scale <- stats::setNames(rep(NA_real_, p), feature_id)
  diagnostics <- vector("list", length(chunks))
  geometry <- NULL
  for (j in seq_along(chunks)) {
    z <- chunks[[j]]
    i <- z$index
    E[, i] <- z$working_error
    V[, i] <- z$working_variance
    field_scale[i] <- z$field_scale
    diagnostics[[j]] <- z$diagnostics
    if (!is.null(z$geometry)) {
      if (is.null(geometry)) {
        geometry <- z$geometry
      } else {
        .mgcvst_geometry_equal(geometry, z$geometry, geometry_tol)
      }
    }
  }
  diagnostics <- do.call(rbind, diagnostics)
  diagnostics <- diagnostics[order(diagnostics$index), , drop = FALSE]
  rownames(diagnostics) <- NULL
  failures <- diagnostics[!diagnostics$success, c(
    "index", "feature_id", "error_class", "error_message", "error_call"
  ), drop = FALSE]
  rownames(failures) <- NULL

  structure(
    list(
      feature_id = feature_id,
      working_error = E,
      working_variance = V,
      field_scale = field_scale,
      geometry = geometry,
      diagnostics = diagnostics,
      failures = failures,
      timing = list(
        elapsed = elapsed,
        workers = workers,
        chunks = length(chunks),
        chunk_size = chunk_size,
        backend = class(BPPARAM)[1L]
      ),
      geometry_tol = geometry_tol,
      source_files = source_files,
      call = call
    ),
    class = c("mgcvST_fit", "mgcvST")
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
  cat("  successful:", sum(x$diagnostics$success), "\n")
  cat("  failures:", nrow(x$failures), "\n")
  cat("  backend:", x$timing$backend, "with", x$timing$workers, "worker(s)\n")
  cat("  elapsed seconds:", format(x$timing$elapsed), "\n")
  cat("  object size:", format(utils::object.size(x), units = "auto"), "\n")
  invisible(x)
}

#' Direct covariance score test for two SPDE mgcv fits
#'
#' Finds the sole SPDE smooth in each ordinary `mgcv::gam` fit, validates exact
#' row identifiers and relative agreement of the training basis and unscaled
#' precision matrices, and delegates to [rkhs_covariance_score_irls()]. No
#' score mathematics is reimplemented here.
#'
#' The maximum absolute difference in `B` or unscaled `Q` must not exceed
#' `precision_tol * max(1, max(abs(M1)), max(abs(M2)))`. Thus changed
#' dimensions, coefficient order, row order, or fixed-kappa precision stop
#' with an error rather than being coerced.
#'
#' @param fit1,fit2 Converged `mgcv::gam` fits with exactly one smooth, which
#'   must inherit from `spde.smooth` or `spdePC.smooth`.
#' @param method Calibration method passed to
#'   [rkhs_covariance_score_irls()].
#' @param normal_min_eff_rank Effective-rank gate for the normal approximation.
#' @param precision_tol Positive relative tolerance for basis and unscaled-Q
#'   comparison. The default is `1e-9`.
#' @return An `rkhs_covariance_score` object with conditional IRLS metadata.
#' @export
score_test <- function(
    fit1, fit2,
    method = c("liu", "davies", "normal", "auto", "signed"),
    normal_min_eff_rank = 30, precision_tol = 1e-9) {
  method <- match.arg(method)
  precision_tol <- as.numeric(precision_tol)
  if (length(precision_tol) != 1L || !is.finite(precision_tol) ||
      precision_tol <= 0) {
    stop("precision_tol must be one positive finite value.")
  }
  i1 <- .mgcvst_score_smooth_index(fit1)
  i2 <- .mgcvst_score_smooth_index(fit2)
  n1 <- length(fit1$linear.predictors)
  n2 <- length(fit2$linear.predictors)
  if (n1 != n2) stop("The two fits use different numbers of rows.")
  if (!identical(.mgcvst_row_id(fit1, n1), .mgcvst_row_id(fit2, n2))) {
    stop("The two fits do not share identical row identifiers and ordering.")
  }
  rkhs_covariance_score_irls(
    fit1, fit2, smooth_index1 = i1, smooth_index2 = i2,
    method = method, normal_min_eff_rank = normal_min_eff_rank,
    geometry_tol = precision_tol
  )
}

# Normalize explicit feature-ID or feature-index pairs to integer indices.
.mgcvst_pair_index <- function(pairs, feature_id) {
  if (is.null(pairs)) {
    stop(
      "pairs is required and must explicitly identify tests; mgcvST.test() ",
      "never expands all feature combinations."
    )
  }
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

# Construct each feature-level score summary once for the Liu pair engine.
.mgcvst_liu_summaries <- function(fitmgcvST, used, verbose) {
  geometry <- fitmgcvST$geometry
  Q <- Matrix::forceSymmetric(Matrix::Matrix(geometry$Q, sparse = TRUE))
  R <- Matrix::chol(Q)
  Rinv <- Matrix::solve(R, Matrix::Diagonal(nrow(Q)))
  T0 <- .magic_mm(geometry$B, as.matrix(Rinv))

  a <- H <- vector("list", length(used))
  success <- rep(FALSE, length(used))
  error <- rep(NA_character_, length(used))
  t0 <- proc.time()[["elapsed"]]
  for (j in seq_along(used)) {
    i <- used[j]
    z <- tryCatch(
      {
        T <- sqrt(fitmgcvST$field_scale[i]) * T0
        op <- .rkhs_score_operator_factor(
          T, fitmgcvST$working_variance[, i], geometry$X,
          field_scale = fitmgcvST$field_scale[i]
        )
        S <- rkhs_score_summary(fitmgcvST$working_error[, i], op)
        list(a = S$a, H = S$H, error = NULL)
      },
      error = function(e) list(a = NULL, H = NULL, error = e)
    )
    if (is.null(z$error)) {
      a[[j]] <- z$a
      H[[j]] <- z$H
      success[j] <- TRUE
    } else {
      error[j] <- .mgcvst_condition(z$error)$message
    }
    if (verbose && (j %% 100L == 0L || j == length(used))) {
      message("Constructed Liu summaries for ", j, " of ", length(used),
              " features.")
    }
  }
  list(
    used = used, a = a, H = H, success = success, error = error,
    elapsed = proc.time()[["elapsed"]] - t0
  )
}

# Evaluate Liu-calibrated pairs from feature summaries in bounded C++ blocks.
.mgcvst_liu_pairs <- function(index, pair_index, feature_id, summaries, threads,
                              chunk_size, tol, verbose) {
  active <- which(summaries$success)
  active_feature <- summaries$used[active]
  local <- matrix(match(index, active_feature), ncol = 2L)
  summary_failed <- !stats::complete.cases(local)

  k <- nrow(index)
  score <- information <- effective_rank <- p_value <- rep(NA_real_, k)
  status <- rep("pending", k)
  error <- rep(NA_character_, k)
  if (any(summary_failed)) {
    status[summary_failed] <- "score_summary_failed"
    for (j in which(summary_failed)) {
      failed <- index[j, ][!(index[j, ] %in% active_feature)]
      failed_local <- match(failed, summaries$used)
      error[j] <- paste(
        paste0(feature_id[failed], ": ",
               summaries$error[failed_local]),
        collapse = " | "
      )
    }
  }

  rows <- which(!summary_failed)
  if (!length(rows)) {
    return(list(
      result = data.frame(
        pair_index = pair_index, score = score, information = information,
        effective_rank = effective_rank, p_value = p_value,
        status = status, error = error, stringsAsFactors = FALSE
      ),
      elapsed = 0
    ))
  }

  avec <- do.call(cbind, summaries$a[active])
  H <- summaries$H[active]
  t0 <- proc.time()[["elapsed"]]
  dense <- ncol(avec)^2 <= 4 * length(rows)
  if (dense) {
    G <- CppMatrix::matrixMultiply(avec, avec, transA = TRUE)
    score[rows] <- G[cbind(local[rows, 1L], local[rows, 2L])]
    rm(G)
  }

  starts <- seq.int(1L, length(rows), by = chunk_size)
  for (b in seq_along(starts)) {
    z <- rows[starts[b]:min(length(rows), starts[b] + chunk_size - 1L)]
    pair <- local[z, , drop = FALSE]
    if (!dense) {
      score[z] <- colSums(
        avec[, pair[, 1L], drop = FALSE] *
          avec[, pair[, 2L], drop = FALSE]
      )
    }
    moments <- mgcvst_pair_trace_powers_cpp(
      H, pair, maxPower = 4L, threads = threads
    )
    good <- apply(is.finite(moments), 1L, all) &
      moments[, 1L] > tol & moments[, 2L] > 0 &
      moments[, 3L] > 0 & moments[, 4L] > 0
    if (any(good)) {
      zi <- z[good]
      M <- moments[good, , drop = FALSE]
      liu <- .liu_squared_score_moments(
        abs(score[zi]), M[, 1L], M[, 2L], M[, 3L], M[, 4L]
      )
      information[zi] <- M[, 1L]
      effective_rank[zi] <- M[, 1L]^2 / M[, 2L]
      p_value[zi] <- liu$p_value
      valid <- is.finite(p_value[zi]) & p_value[zi] >= 0 & p_value[zi] <= 1
      status[zi[valid]] <- "ok"
      status[zi[!valid]] <- "score_failed"
      error[zi[!valid]] <- "Liu calibration returned an invalid p-value."
    }
    if (any(!good)) {
      status[z[!good]] <- "score_failed"
      error[z[!good]] <- "Liu trace moments were non-finite or non-positive."
    }
    if (verbose && (b %% 10L == 0L || b == length(starts))) {
      message("Evaluated Liu block ", b, " of ", length(starts), ".")
    }
  }

  list(
    result = data.frame(
      pair_index = pair_index, score = score, information = information,
      effective_rank = effective_rank, p_value = p_value,
      status = status, error = error, stringsAsFactors = FALSE
    ),
    elapsed = proc.time()[["elapsed"]] - t0
  )
}

# Evaluate one pair chunk from compact feature summaries.
.mgcvst_test_chunk <- function(payload, geometry, calibration,
                               normal_min_eff_rank, tol) {
  .mgcvst_thread_limit()
  pairs <- payload$pairs
  k <- nrow(pairs)
  used <- sort(unique(as.vector(pairs)))
  summaries <- vector("list", length(payload$feature_id))
  for (i in used) {
    op <- rkhs_score_operator(
      geometry$B, geometry$Q, payload$working_variance[, i],
      geometry$X, payload$field_scale[i]
    )
    summaries[[i]] <- rkhs_score_summary(
      payload$working_error[, i], op
    )
  }
  out <- vector("list", k)
  for (j in seq_len(k)) {
    i1 <- pairs[j, 1L]
    i2 <- pairs[j, 2L]
    test_result <- tryCatch(
      {
        S1 <- summaries[[i1]]
        S2 <- summaries[[i2]]
        U <- as.numeric(crossprod(S1$a, S2$a))
        cal <- rkhs_score_calibrate(
          U, S1$H, S2$H, method = calibration,
          normal_min_eff_rank = normal_min_eff_rank, tol = tol
        )
        value <- c(list(
          score = U,
          signed_score = U,
          statistic = U^2,
          quadratic_statistic = U^2,
          normal_statistic = if (identical(cal$method, "normal") &&
                                 cal$information > tol) {
            U / sqrt(cal$information)
          } else {
            NA_real_
          }
        ), cal)
        list(value = value, error = NULL)
      },
      error = function(e) list(value = NULL, error = e)
    )
    if (!is.null(test_result$error)) {
      err <- .mgcvst_condition(test_result$error)
      out[[j]] <- data.frame(
        pair_index = payload$pair_index[j],
        feature1 = payload$feature_id[i1],
        feature2 = payload$feature_id[i2],
        score = NA_real_, signed_score = NA_real_,
        statistic = NA_real_, quadratic_statistic = NA_real_,
        normal_statistic = NA_real_, information = NA_real_,
        effective_rank = NA_real_, p_value = NA_real_,
        calibration = NA_character_, status = "score_failed",
        davies_ifault = NA_integer_, error = err$message,
        stringsAsFactors = FALSE
      )
      next
    }
    value <- test_result$value
    out[[j]] <- data.frame(
      pair_index = payload$pair_index[j],
      feature1 = payload$feature_id[i1],
      feature2 = payload$feature_id[i2],
      score = value$score,
      signed_score = value$signed_score,
      statistic = value$statistic,
      quadratic_statistic = value$quadratic_statistic,
      normal_statistic = value$normal_statistic,
      information = value$information,
      effective_rank = value$effective_rank,
      p_value = value$p_value,
      calibration = value$method,
      status = value$status,
      davies_ifault = if (is.null(value$davies_ifault)) NA_integer_ else
        value$davies_ifault,
      error = NA_character_,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, out)
}

#' Covariance score tests for an explicit pair universe
#'
#' Tests every pair in `pairs`, plus any additional `highlight` pairs, using
#' the fixed compact summaries in `fitmgcvST`. The primary statistic is the
#' squared signed score `U^2`. The default Liu calibration matches the first
#' four moments of `U^2` directly from trace powers of `H1 H2`, without a
#' per-pair spectral decomposition. Davies remains available when explicitly
#' requested for a smaller pair universe. No GAM is refitted.
#'
#' `q.value`, `FDR`, and `method` apply to the resulting pair-covariance
#' p-values. Automatic discoveries are the pairs passing that threshold.
#' Highlighted pairs are always retained in the returned set regardless of
#' significance, while any fit or score failure remains explicit. The
#' function does not generate all pairwise combinations: callers must define
#' the tested universe, which prevents accidental allocation of
#' `choose(10000, 2)` rows.
#'
#' Under Liu calibration, each feature summary is constructed once. Dense pair
#' universes use `CppMatrix::matrixMultiply(..., transA = TRUE)` for all signed
#' scores, followed by bounded, multithreaded internal C++ blocks for the Liu trace
#' moments. This avoids copying the full collection of `H` matrices to Snow
#' workers. Other calibrations retain BiocParallel pair chunks.
#'
#' @param fitmgcvST A compact object returned by [mgcvST.estimate()].
#' @param q.value Pair-covariance discovery threshold in `(0, 1]`.
#' @param FDR Logical; if `TRUE`, adjust pair p-values before applying
#'   `q.value`.
#' @param method Method passed to `stats::p.adjust()` when `FDR = TRUE`.
#' @param BPPARAM A `BiocParallelParam`; defaults to `SerialParam()`.
#' @param ... Reserved for future score options; currently unused arguments
#'   are rejected.
#' @param pairs Required two-column tested-pair universe using feature IDs or
#'   one-based indices.
#' @param highlight Optional two-column pairs that must be tested and retained.
#'   They may overlap `pairs`; non-overlapping highlights are appended.
#' @param calibration Score calibration passed to [rkhs_score_calibrate()].
#' @param normal_min_eff_rank Effective-rank gate for normal calibration.
#' @param tol Positive numerical score tolerance.
#' @param chunk_size Positive number of tested pairs per task or C++ block.
#'   Liu calibration defaults to 10,000 pairs per interruptible block.
#' @param threads Positive number of OpenMP threads for the Liu C++ kernel.
#'   `NULL` uses `bpworkers(BPPARAM)` without launching Snow workers.
#' @param verbose Logical; report Liu summary and block progress.
#' @return A compact `mgcvST_test` object containing pair results, adjusted
#'   p-values, discovery/highlight/retention flags, threshold metadata,
#'   explicit failures, and timing metadata.
#' @export
mgcvST.test <- function(
    fitmgcvST, q.value = 0.05, FDR = TRUE, method = "BH",
    BPPARAM = BiocParallel::SerialParam(), ...,
    pairs = NULL, highlight = NULL,
    calibration = c("liu", "davies", "normal", "auto", "signed"),
    normal_min_eff_rank = 30, tol = 1e-10, chunk_size = NULL,
    threads = NULL, verbose = FALSE) {
  if (!inherits(fitmgcvST, "mgcvST_fit")) {
    stop("fitmgcvST must be returned by mgcvST.estimate().")
  }
  if (is.null(fitmgcvST$geometry)) {
    stop("fitmgcvST has no successful feature geometry to test.")
  }
  q.value <- as.numeric(q.value)
  if (length(q.value) != 1L || !is.finite(q.value) ||
      q.value <= 0 || q.value > 1) {
    stop("q.value must be one finite value in (0, 1].")
  }
  if (!is.logical(FDR) || length(FDR) != 1L || is.na(FDR)) {
    stop("FDR must be TRUE or FALSE.")
  }
  if (!is.character(method) || length(method) != 1L || is.na(method) ||
      !(method %in% stats::p.adjust.methods)) {
    stop("method must be one of stats::p.adjust.methods.")
  }
  if (!inherits(BPPARAM, "BiocParallelParam")) {
    stop("BPPARAM must inherit from 'BiocParallelParam'.")
  }
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("verbose must be TRUE or FALSE.")
  }
  unused <- list(...)
  if (length(unused)) {
    stop("Unused arguments in ...: ", paste(names(unused), collapse = ", "))
  }
  calibration <- match.arg(calibration)
  normal_min_eff_rank <- as.numeric(normal_min_eff_rank)
  if (length(normal_min_eff_rank) != 1L ||
      !is.finite(normal_min_eff_rank) || normal_min_eff_rank <= 0) {
    stop("normal_min_eff_rank must be one positive finite value.")
  }
  tol <- as.numeric(tol)
  if (length(tol) != 1L || !is.finite(tol) || tol <= 0) {
    stop("tol must be one positive finite value.")
  }
  if (is.null(threads)) threads <- BiocParallel::bpworkers(BPPARAM)
  threads <- as.integer(threads)
  if (length(threads) != 1L || is.na(threads) || threads < 1L) {
    stop("threads must be one positive integer.")
  }
  .mgcvst_thread_limit()

  index <- .mgcvst_pair_index(pairs, fitmgcvST$feature_id)
  highlight_index <- matrix(integer(), nrow = 0L, ncol = 2L)
  if (!is.null(highlight)) {
    highlight_index <- .mgcvst_pair_index(
      highlight, fitmgcvST$feature_id
    )
  }
  n_feature <- length(fitmgcvST$feature_id)
  key <- (index[, 1L] - 1) * n_feature + index[, 2L]
  highlight_key <- if (nrow(highlight_index)) {
    (highlight_index[, 1L] - 1) * n_feature + highlight_index[, 2L]
  } else {
    character()
  }
  extra <- which(!(highlight_key %in% key))
  if (length(extra)) {
    index <- rbind(index, highlight_index[extra, , drop = FALSE])
    key <- c(key, highlight_key[extra])
  }
  highlighted <- key %in% highlight_key

  success <- as.logical(fitmgcvST$diagnostics$success)
  i1 <- index[, 1L]
  i2 <- index[, 2L]
  fit_failed <- !success[i1] | !success[i2]
  status <- ifelse(fit_failed, "feature_fit_failed", "pending")
  error <- rep(NA_character_, nrow(index))
  for (j in which(fit_failed)) {
    failed <- c(i1[j], i2[j])[!success[c(i1[j], i2[j])]]
    error[j] <- paste(
      paste0(
        fitmgcvST$feature_id[failed], ": ",
        fitmgcvST$diagnostics$error_message[failed]
      ),
      collapse = " | "
    )
  }

  result <- data.frame(
    pair_index = seq_len(nrow(index)),
    feature1 = fitmgcvST$feature_id[i1],
    feature2 = fitmgcvST$feature_id[i2],
    score = NA_real_,
    signed_score = NA_real_,
    statistic = NA_real_,
    quadratic_statistic = NA_real_,
    normal_statistic = NA_real_,
    information = NA_real_,
    effective_rank = NA_real_,
    p_value = NA_real_,
    p_positive = NA_real_,
    p_negative = NA_real_,
    p_adjusted = NA_real_,
    p_positive_adjusted = NA_real_,
    p_negative_adjusted = NA_real_,
    discovered = FALSE,
    discovered_positive = FALSE,
    discovered_negative = FALSE,
    highlighted = highlighted,
    retained = highlighted,
    calibration = NA_character_,
    status = status,
    davies_ifault = NA_integer_,
    error = error,
    stringsAsFactors = FALSE
  )

  tested_rows <- which(!fit_failed)
  workers <- if (length(tested_rows) && calibration != "liu") {
    max(1L, min(length(tested_rows), BiocParallel::bpworkers(BPPARAM)))
  } else if (length(tested_rows)) {
    threads
  } else {
    0L
  }
  if (is.null(chunk_size)) {
    chunk_size <- if (calibration == "liu") {
      10000L
    } else if (workers > 0L) {
      ceiling(length(tested_rows) / workers)
    } else {
      1L
    }
  }
  chunk_size <- as.integer(chunk_size)
  if (length(chunk_size) != 1L || is.na(chunk_size) || chunk_size < 1L) {
    stop("chunk_size must be one positive integer.")
  }

  elapsed <- summary_elapsed <- 0
  chunks <- list()
  if (length(tested_rows)) {
    if (calibration == "liu") {
      used <- sort(unique(as.vector(index[tested_rows, , drop = FALSE])))
      summaries <- .mgcvst_liu_summaries(fitmgcvST, used, verbose)
      summary_elapsed <- summaries$elapsed
      evaluated <- .mgcvst_liu_pairs(
        index[tested_rows, , drop = FALSE], tested_rows,
        fitmgcvST$feature_id, summaries, threads, chunk_size, tol, verbose
      )
      elapsed <- summary_elapsed + evaluated$elapsed
      evaluated <- evaluated$result
      target <- evaluated$pair_index
      result$score[target] <- evaluated$score
      result$signed_score[target] <- evaluated$score
      result$statistic[target] <- evaluated$score^2
      result$quadratic_statistic[target] <- evaluated$score^2
      result$information[target] <- evaluated$information
      result$effective_rank[target] <- evaluated$effective_rank
      result$p_value[target] <- evaluated$p_value
      result$calibration[target[evaluated$status == "ok"]] <- "liu"
      result$status[target] <- evaluated$status
      result$error[target] <- evaluated$error
      chunks <- seq_len(ceiling(length(tested_rows) / chunk_size))
    } else {
      chunks <- split(
        tested_rows, ceiling(seq_along(tested_rows) / chunk_size)
      )
      payload <- lapply(chunks, function(rows) {
        used <- unique(as.vector(t(index[rows, , drop = FALSE])))
        local <- match(index[rows, , drop = FALSE], used)
        local <- matrix(local, ncol = 2L)
        list(
          pair_index = rows,
          pairs = local,
          feature_id = fitmgcvST$feature_id[used],
          working_error = fitmgcvST$working_error[, used, drop = FALSE],
          working_variance = fitmgcvST$working_variance[, used, drop = FALSE],
          field_scale = fitmgcvST$field_scale[used]
        )
      })
      worker_bundle <- .mgcvst_worker_bundle()
      test_chunk <- get(".mgcvst_test_chunk", envir = worker_bundle,
                        inherits = FALSE)
      t0 <- proc.time()[["elapsed"]]
      evaluated <- BiocParallel::bplapply(
        payload, test_chunk,
        geometry = fitmgcvST$geometry, calibration = calibration,
        normal_min_eff_rank = normal_min_eff_rank, tol = tol,
        BPPARAM = BPPARAM
      )
      elapsed <- proc.time()[["elapsed"]] - t0
      evaluated <- do.call(rbind, evaluated)
      core <- c(
        "score", "signed_score", "statistic", "quadratic_statistic",
        "normal_statistic", "information", "effective_rank", "p_value",
        "calibration", "status", "davies_ifault", "error"
      )
      target <- match(evaluated$pair_index, result$pair_index)
      result[target, core] <- evaluated[, core, drop = FALSE]
    }
  }

  valid <- result$status == "ok" & is.finite(result$p_value) &
    result$p_value >= 0 & result$p_value <= 1
  result$p_positive[valid] <- ifelse(
    result$signed_score[valid] >= 0,
    result$p_value[valid] / 2,
    1 - result$p_value[valid] / 2
  )
  result$p_negative[valid] <- ifelse(
    result$signed_score[valid] <= 0,
    result$p_value[valid] / 2,
    1 - result$p_value[valid] / 2
  )
  if (FDR) {
    result$p_adjusted[valid] <- stats::p.adjust(
      result$p_value[valid], method = method
    )
    result$p_positive_adjusted[valid] <- stats::p.adjust(
      result$p_positive[valid], method = method
    )
    result$p_negative_adjusted[valid] <- stats::p.adjust(
      result$p_negative[valid], method = method
    )
  } else {
    result$p_adjusted[valid] <- result$p_value[valid]
    result$p_positive_adjusted[valid] <- result$p_positive[valid]
    result$p_negative_adjusted[valid] <- result$p_negative[valid]
  }
  result$discovered <- valid & result$p_adjusted <= q.value
  result$discovered_positive <- valid &
    result$p_positive_adjusted <= q.value
  result$discovered_negative <- valid &
    result$p_negative_adjusted <= q.value
  result$retained <- result$highlighted | result$discovered
  raw_threshold <- if (any(result$discovered)) {
    max(result$p_value[result$discovered])
  } else {
    NA_real_
  }

  structure(
    list(
      results = result,
      threshold = list(
        q_value = q.value,
        FDR = FDR,
        adjustment_method = if (FDR) method else "none",
        raw_p_threshold = raw_threshold
      ),
      discoveries = list(
        pairs_requested = nrow(index),
        pairs_fit_failed = sum(fit_failed),
        pairs_tested = length(tested_rows),
        pairs_with_p_value = sum(valid),
        pairs_failed = sum(!fit_failed & !valid),
        pairs_discovered = sum(result$discovered),
        pairs_discovered_positive = sum(result$discovered_positive),
        pairs_discovered_negative = sum(result$discovered_negative),
        pairs_highlighted = sum(result$highlighted),
        pairs_retained = sum(result$retained)
      ),
      pair_contract = paste0(
        "explicit_tested_pair_universe_with_FDR_discoveries_",
        "union_force_retained_highlights"
      ),
      timing = list(
        elapsed = elapsed,
        summary_elapsed = summary_elapsed,
        pair_elapsed = elapsed - summary_elapsed,
        workers = workers,
        chunks = length(chunks),
        backend = if (calibration == "liu") {
          "mgcvSTOpenMP"
        } else {
          class(BPPARAM)[1L]
        }
      ),
      calibration = calibration,
      call = match.call()
    ),
    class = "mgcvST_test"
  )
}

#' Print covariance-test diagnostics
#'
#' @param x An `mgcvST_test` object.
#' @param ... Unused.
#' @return `x`, invisibly.
#' @export
print.mgcvST_test <- function(x, ...) {
  cat("mgcvST quadratic-form covariance tests\n")
  cat("  tested pair universe:", x$discoveries$pairs_requested, "\n")
  cat("  pairs with p-value:", x$discoveries$pairs_with_p_value, "\n")
  cat("  FDR discoveries:", x$discoveries$pairs_discovered, "\n")
  cat("  highlighted pairs:", x$discoveries$pairs_highlighted, "\n")
  cat("  retained union:", x$discoveries$pairs_retained, "\n")
  cat("  decision threshold:", format(x$threshold$q_value), "\n")
  invisible(x)
}
