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
  if (length(fit$dispersion) != n || length(fit$lambda) != n ||
      ncol(fit$score_a) != n) {
    stop("The compact fit dimensions are incompatible with feature_id.")
  }
  is.finite(fit$dispersion) & fit$dispersion > 0 &
    is.finite(fit$lambda) & fit$lambda > 0 &
    colSums(!is.finite(fit$score_a)) == 0L
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

