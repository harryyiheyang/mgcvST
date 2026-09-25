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

# Locate the sole SPDE or SPDE-PC smooth supported by the high-throughput API.
.mgcvst_spde_index <- function(fit) {
  if (!inherits(fit, "gam")) stop("fit must inherit from class 'gam'.")
  idx <- which(vapply(
    fit$smooth,
    function(s) inherits(s, "spde.smooth") || inherits(s, "spdePC.smooth"),
    logical(1)
  ))
  if (!length(idx)) stop("The fit does not contain an SPDE or SPDE-PC smooth.")
  if (length(idx) != 1L) {
    stop("The fit must contain exactly one SPDE or SPDE-PC smooth.")
  }
  if (length(fit$smooth) != 1L) {
    stop(
      "The fit contains additional nuisance smooths; their covariance is not ",
      "implemented by the validated score operator."
    )
  }
  idx
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

# Derive per-feature field scales from primary fit parameters.
.mgcvst_field_scale <- function(fit) {
  n <- length(fit$feature_id)
  if (is.null(fit$dispersion) || is.null(fit$lambda)) {
    stop("The fit must contain dispersion and lambda.")
  }
  available <- .mgcvst_feature_available(fit)
  dispersion <- as.numeric(fit$dispersion)
  lambda <- as.numeric(fit$lambda)
  if (length(dispersion) != n || length(lambda) != n) {
    stop("dispersion and lambda must contain one value per feature.")
  }
  scale <- dispersion / lambda
  if (any(!is.finite(scale[available])) || any(scale[available] <= 0)) {
    stop("Available features require positive finite dispersion/lambda scales.")
  }
  stats::setNames(scale, fit$feature_id)
}

# Reduce one gam fit to the working summaries and shared geometry used downstream.
.mgcvst_compact_fit <- function(fit, retain_smooth = FALSE, diagnostics = TRUE) {
  W <- rkhs_extract_working_model(fit)
  L <- fit$.taps_score_X
  if (is.null(L)) L <- .gam_training_lpmatrix(fit)
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
    smooth_columns = S$columns,
    score_precision_psd = S$score_precision_psd
  )
  wood_p_value <- smooth_edf <- wood_statistic <- NA_real_
  if (diagnostics) {
    fit_summary <- summary(fit)
    smooth_table <- fit_summary$s.table[smooth_index, , drop = FALSE]
    wood_p_value <- as.numeric(smooth_table[1L, "p-value"])
    smooth_edf <- as.numeric(smooth_table[1L, "edf"])
    statistic_column <- ncol(smooth_table) - 1L
    wood_statistic <- as.numeric(smooth_table[1L, statistic_column])
  }
  residual_df <- as.numeric(fit$df.residual)
  if (length(residual_df) != 1L) residual_df <- NA_real_
  parameters <- W$family_parameters
  if (is.null(parameters)) parameters <- numeric()
  ans <- list(
    working_error = W$working_error,
    working_variance = W$working_variance,
    smoothing_parameter = S$sp,
    dispersion = W$dispersion,
    family = W$family,
    family_parameters = as.numeric(parameters),
    family_parameters_diagnostic = paste(
      format(parameters, digits = 16), collapse = ","
    ),
    smooth_label = S$label,
    wood_p_value = wood_p_value,
    wood_statistic = wood_statistic,
    smooth_edf = smooth_edf,
    residual_df = residual_df,
    geometry = geometry
  )
  if (retain_smooth) {
    s <- fit$smooth[[smooth_index]]
    fit_basis <- L[, S$columns, drop = FALSE]
    fit_penalty <- as.matrix(s$S[[1L]])
    coefficients <- as.numeric(stats::coef(fit)[S$columns])
    if (!all(dim(fit_penalty) == c(ncol(fit_basis), ncol(fit_basis))) ||
        length(coefficients) != ncol(fit_basis) ||
        any(!is.finite(c(fit_basis, fit_penalty, coefficients)))) {
      stop("The retained reduced smooth fit is not finite and dimensionally aligned.")
    }
    metadata <- list(
      smooth_class = class(s)[1L],
      smooth_label = S$label,
      coefficient_columns = as.integer(S$columns),
      reduced_dimension = ncol(fit_basis),
      raw_dimension = s[["raw.dimension"]],
      projected_dimension = s[["pc_full_dimension"]],
      pc_retained_dimension = s[["pc_retained_dimension"]],
      pc_cutoff = s[["pc_cutoff"]],
      kappa = s[["kappa"]]
    )
    ans$smooth_coefficients <- coefficients
    ans$fit_geometry <- list(
      B = fit_basis,
      P = fit_penalty,
      row_id = .mgcvst_row_id(fit, nrow(L)),
      metadata = metadata
    )
  }
  ans
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
    identical(requested, "davies") && identical(used, "liu")
  cache <- if (cacheable) attr(score, "marginal_spectrum", exact = TRUE) else NULL
  list(p_value = p_value, requested_method = requested,
       method = used, fallback = fallback, cache = cache)
}

# Clone the minimal package function closure needed on remote workers.
.mgcvst_worker_bundle <- function() {
  names <- c(
    ".mgcvst_thread_limit", ".mgcvst_worker_initialize",
    ".mgcvst_spde_index", ".mgcvst_row_id", ".mgcvst_compact_fit",
    ".mgcvst_condition", ".mgcvst_fit_chunk",
    ".mgcvst_capture_marginal", ".mgcvst_marginal_geometry",
    ".mgcvst_null_score_setup", ".mgcvst_null_score_spectrum",
    ".mgcvst_null_score_test",
    ".mgcvst_test_chunk", ".mgcvst_marginal_score", ".working_family_id",
    "taps_score_test", ".mgcvst_marginal_spectrum", ".mgcvst_marginal_working",
    ".mgcvst_marginal_matrixsqrt", ".mgcvst_marginal_moments",
    ".mgcvst_marginal_liu", ".mgcvst_marginal_davies",
    ".gam_training_lpmatrix", ".gam_single_smooth",
    ".mgcvst_expand_penalty", ".mgcvst_model_geometry",
    ".mgcvst_geometry_signature", ".mgcvst_model_sp",
    ".mgcvst_cached_model_geometry", ".mgcvst_training_design",
    ".mgcvst_nuisance_state",
    ".mgcvst_model_fit_one",
    ".mgcvst_model_fit_chunk",
    ".mgcvst_spde_factor", ".mgcvst_spde_factor_base",
    ".mgcvst_model_fixed_factors",
    ".mgcvst_model_sparse_constrained_solver",
    ".mgcvst_pack_symmetric", ".mgcvst_unpack_symmetric",
    ".mgcvst_pack_score_state", ".mgcvst_unpack_score_state",
    ".mgcvst_dense_temp_dir", ".mgcvst_dense_cleanup",
    ".mgcvst_dense_pair_groups",
    ".mgcvst_dense_pair_chunk", ".mgcvst_model_state_shard",
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
  worker_trace_powers <- function(matrixList, pairs, maxPower = 4L,
                                  threads = 1L) {
    loadNamespace("mgcvST")
    .Call(
      "_mgcvST_mgcvst_pair_trace_powers_cpp",
      matrixList, pairs, maxPower, threads, PACKAGE = "mgcvST"
    )
  }
  environment(worker_trace_powers) <- bundle
  assign(
    "mgcvst_pair_trace_powers_cpp", worker_trace_powers,
    envir = bundle
  )
  bundle
}

# Fit and compact one feature chunk.
.mgcvst_fit_chunk <- function(payload, G0, family_raw, method, control,
                              gam_args, source_files, worker_init, init_key,
                              marginal_test, marginal_args,
                              retain_smooth,
                              diagnostics = TRUE, retain_marginal = FALSE,
                              routed_family_raw = NULL) {
  .mgcvst_worker_initialize(source_files, worker_init, init_key)
  ids <- payload$index
  Y <- payload$Y
  feature_id <- payload$feature_id
  n <- ncol(Y)
  k <- nrow(Y)
  E <- matrix(NA_real_, nrow = n, ncol = k)
  V <- matrix(NA_real_, nrow = n, ncol = k)
  smoothing_parameter <- dispersion <- rep(NA_real_, k)
  marginal_p_value <- wood_p_value <- wood_statistic <- smooth_edf <-
    residual_df <- criterion <- rep(NA_real_, k)
  converged <- smoothing_converged <- rep(FALSE, k)
  smoothing_converged[] <- NA
  family <- family_parameters_diagnostic <- smooth_label <- criterion_name <-
    rep(NA_character_, k)
  family_parameters <- vector("list", k)
  irls_iterations <- rep(NA_integer_, k)
  outer_convergence <- rep(NA_character_, k)
  fit_seconds <- summary_seconds <- marginal_seconds <- rep(0, k)
  error_class <- error_message <- error_call <- rep(NA_character_, k)
  geometry <- NULL
  fit_geometry <- NULL
  marginal_geometry <- NULL
  marginal_requested_method <- marginal_method <- rep(NA_character_, k)
  marginal_fallback <- rep(NA, k)
  design_cache <- new.env(parent = emptyenv())
  L0 <- attr(G0, "training_X")
  if (!is.null(L0)) {
    pseudo <- G0
    pseudo$model <- G0$mf
    pseudo$coefficients <- stats::setNames(numeric(ncol(L0)), G0$term.names)
    pseudo$linear.predictors <- numeric(nrow(L0))
    class(pseudo) <- c("gam", "glm", "lm")
    design_cache$frozen <- TRUE
    design_cache$L <- L0
    design_cache$geometry <- .mgcvst_model_geometry(pseudo, L0)
  }
  marginal_state <- if (retain_marginal) vector("list", k) else NULL
  coefficient_count <- length(seq.int(
    G0$smooth[[1L]]$first.para, G0$smooth[[1L]]$last.para
  ))
  smooth_coefficients <- if (retain_smooth) {
    matrix(NA_real_, nrow = k, ncol = coefficient_count)
  } else {
    NULL
  }
  response_index <- attr(G0$terms, "response")
  if (length(response_index) != 1L || response_index < 1L ||
      is.null(G0$mf) || response_index > ncol(G0$mf)) {
    stop("G does not contain a reusable model-frame response column.")
  }

  # Poisson prescreen routing (R/family-prescreen.R). Absent flags keep the
  # model family for every feature, so the unscreened path is unchanged.
  poisson_route <- payload$poisson
  if (is.null(poisson_route)) poisson_route <- rep(FALSE, k)
  prescreen_phi <- payload$prescreen_phi
  if (is.null(prescreen_phi)) prescreen_phi <- rep(NA_real_, k)
  family_used <- rep(NA_character_, k)

  target_index <- which(vapply(
    G0$smooth, function(s) identical(s$score.component, "global"), logical(1L)
  ))
  null_setup <- .mgcvst_null_score_setup(G0, target_index, attr(G0, "null_spec"))
  for (j in seq_len(k)) {
    response <- as.numeric(Y[j, ])
    family_serialized <- if (isTRUE(poisson_route[j])) routed_family_raw else family_raw
    family_used[j] <- .working_family_id(unserialize(family_serialized)$family)

    null_data <- null_setup$spec$data
    null_data[[null_setup$spec$response]] <- response
    null_fit <- NULL
    t0 <- proc.time()[["elapsed"]]
    marginal_result <- tryCatch(
      {
        null_fit <- do.call(mgcv::bam, c(list(
          formula = null_setup$spec$formula, data = null_data,
          family = unserialize(family_serialized), offset = null_setup$spec$offset,
          method = "fREML", discrete = TRUE, nthreads = 1L,
          control = control), gam_args))
        .mgcvst_marginal_score(
          null_fit, marginal_test = marginal_test,
          marginal_args = marginal_args,
          test_component = null_setup$target_index, setup = null_setup
        )
      },
      error = function(e) e
    )
    marginal_seconds[j] <- proc.time()[["elapsed"]] - t0
    if (inherits(marginal_result, "condition")) {
      err <- .mgcvst_condition(marginal_result)
      error_class[j] <- err$class
      error_message[j] <- err$message
      error_call[j] <- err$call
    } else {
      marginal_p_value[j] <- marginal_result$p_value
      marginal_requested_method[j] <- marginal_result$requested_method
      marginal_method[j] <- marginal_result$method
      marginal_fallback[j] <- marginal_result$fallback
    }
    if (retain_marginal && !inherits(marginal_result, "condition")) {
      marginal_state[[j]] <- list(marginal_cache = marginal_result$cache)
    }

    t0 <- proc.time()[["elapsed"]]
    full_data <- attr(G0, "full_spec")$data
    full_data[[attr(G0, "full_spec")$response]] <- response
    fit_result <- tryCatch(
      list(
        value = do.call(
          mgcv::bam,
           c(list(formula = attr(G0, "full_spec")$formula, data = full_data,
                  family = unserialize(family_serialized), method = "fREML",
                  discrete = TRUE, nthreads = 1L, control = control), gam_args)
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
    if (length(fit$gcv.ubre) == 1L && is.finite(fit$gcv.ubre)) {
      criterion[j] <- as.numeric(fit$gcv.ubre)
      nm <- names(fit$gcv.ubre)
      if (length(nm) == 1L && nzchar(nm)) criterion_name[j] <- nm
    }
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
      list(
        value = {
          fit$.taps_score_X <- .mgcvst_training_design(fit, design_cache)
          if (is.null(design_cache$signature) &&
              !length(source_files) && is.null(worker_init)) {
            design_cache$signature <- .mgcvst_geometry_signature(fit)
            design_cache$L <- fit$.taps_score_X
          }
          .mgcvst_compact_fit(fit, retain_smooth = retain_smooth,
                             diagnostics = diagnostics)
        },
        error = NULL
      ),
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
    if (is.null(geometry)) geometry <- compact$geometry
    if (retain_smooth) {
      if (is.null(fit_geometry)) fit_geometry <- compact$fit_geometry
      smooth_coefficients[j, ] <- compact$smooth_coefficients
    }

    smoothing_parameter[j] <- compact$smoothing_parameter
    dispersion[j] <- compact$dispersion
    family[j] <- compact$family
    family_parameters[[j]] <- compact$family_parameters
    family_parameters_diagnostic[j] <- compact$family_parameters_diagnostic
    smooth_label[j] <- compact$smooth_label
    wood_p_value[j] <- compact$wood_p_value
    wood_statistic[j] <- compact$wood_statistic
    smooth_edf[j] <- compact$smooth_edf
    residual_df[j] <- compact$residual_df

    E[, j] <- compact$working_error
    V[, j] <- compact$working_variance

  }

  diagnostics <- data.frame(
    index = ids,
    feature_id = feature_id,
    family = family,
    smooth_label = smooth_label,
    smoothing_parameter = smoothing_parameter,
    dispersion = dispersion,
    marginal_p_value = marginal_p_value,
    marginal_requested_method = marginal_requested_method,
    marginal_method = marginal_method,
    marginal_fallback = marginal_fallback,
    wood_p_value = wood_p_value,
    wood_statistic = wood_statistic,
    smooth_edf = smooth_edf,
    residual_df = residual_df,
    criterion = criterion,
    criterion_name = criterion_name,
    family_parameters = family_parameters_diagnostic,
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
    prescreen_phi = prescreen_phi,
    family_used = family_used,
    stringsAsFactors = FALSE
  )
  list(
    index = ids,
    working_error = E,
    working_variance = V,
    lambda = smoothing_parameter,
    dispersion = dispersion,
    family_parameters = family_parameters,
    diagnostics = diagnostics,
    geometry = geometry,
    smooth_coefficients = smooth_coefficients,
    fit_geometry = fit_geometry,
    marginal_geometry = marginal_geometry,
    marginal_state = marginal_state
  )
}

#' Estimate compact covariance working summaries
#'
#' Fits each row of `Y` from a reusable `mgcv::gam(fit = FALSE)` setup and
#' retains the fixed numerical summaries needed by [mgcvST.test()]. Each
#' feature first fits the null model with the spatial score smooth removed;
#' marginal screening uses that null PIRLS state, the prepared `G$X`, and the
#' target penalty before the full spatial fit. Wood diagnostics are optional
#' (`diagnostics = TRUE`).
#' The marginal test is the package-local `taps_score_test()`, using the TAPS
#' arithmetic included in mgcvST. No external mgcv.taps installation or sourced
#' score function is required. Full `gam` objects are never
#' retained. Genes are processed in chunks through `BiocParallel`; the default
#' uses the registered backend. On Windows, use a persistent
#' `SnowParam(type = "SOCK")` and pass it to estimation and testing.
#' A fitting, compaction, or marginal-test error is recorded for that feature;
#' the other features continue. A marginal-test error leaves the already
#' constructed compact working model intact and only its marginal p-value
#' unavailable.
#' Marginal diagnostics retain the requested method, actual method, and a
#' Davies-to-Liu fallback flag. Custom callbacks that omit method metadata
#' leave the corresponding diagnostics as `NA`.
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
#' @param BPPARAM A `BiocParallelParam`; defaults to the registered `bpparam()`.
#' @param chunk_size Positive number of genes per task. The default creates at
#'   most one chunk per worker, limiting repeated serialization on SOCK
#'   workers.
#' @param source_files Ordered R files, or directories containing R files, to
#'   source once per worker. Use this for sourced custom S3 smooth methods.
#' @param worker_init Optional zero-argument initialization function run once
#'   per worker after `source_files`.
#' @param marginal_test Optional function implementing the corrected marginal
#'   spatial score interface. `NULL` uses mgcvST's package-local implementation.
#'   A custom function is used only when explicitly supplied.
#' @param offset Optional additional log/link-scale offset for an
#'   [mgcvST.set()] model: a shared observation-length vector or a
#'   feature-by-observation matrix in exactly the same order as `Y`.
#'   Added to the shared formula/setup offset. Covariates and smooths are fixed
#'   by `mgcvST.set()` and cannot vary by gene.
#' @param diagnostics Logical; compute `summary.gam()`/Wood diagnostics.
#'   FALSE (default) leaves Wood fields NA without calling summary.
#'   Basic convergence and fitting diagnostics are still retained.
#' @param retain_marginal Logical; retain minimal null-fit inputs for a later
#'   `mgcvST.marginal()` recalibration call. FALSE by default. Marginal testing
#'   during estimation is always performed; full gam objects are never retained.
#' @param marginal_args Named list of additional marginal-score arguments.
#'   `fit`, `test.component`, and `n_threads` are controlled by mgcvST.
#'   Use `list(method = "liu")` for direct Liu, or the default Davies with
#'   Liu fallback within mgcvST. No additional marginal call is needed.
#' @param retain_smooth Logical; retain the feature-by-coefficient smooth
#'   coefficient matrix and one shared reduced fit basis and unscaled penalty.
#'   This opt-in representation supports prediction and other downstream uses
#'   without retaining full `gam` objects.
#' @param method Retained for API compatibility. Null and full fits use
#'   `mgcv::bam(method = "fREML", discrete = TRUE)`.
#' @param control An `mgcv::gam.control()` object. Internal thread counts are
#'   always forced to one. It may additionally carry `poisson_screen_phi`
#'   (default `1.1`), the Poisson prescreen threshold: with a negative-binomial
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
#' @return A compact object of class `mgcvST_fit` containing marginal score
#'   p-values, Wood audit values, feature IDs, working errors and variances,
#'   separate per-feature `dispersion` and `lambda`, shared score geometry,
#'   optimized fit criterion and its original mgcv name, exact residual degrees
#'   of freedom, timing, and convergence diagnostics. A validated `model.set()`
#'   fit retains one shared `geometry$nuisance_design` and one small conditional
#'   nuisance block per feature in `nuisance_covariance`; full GAM and `Vp`
#'   objects are discarded. When
#'   `retain_smooth = TRUE`, it also contains `smooth_coefficients`,
#'   `fit_basis`, `fit_penalty`, `row_id`, and `basis_metadata`. New objects do
#'   not store a redundant `field_scale`; score methods derive it as
#'   `dispersion / lambda`.
#' @export
mgcvST.estimate <- function(
    Y, G, feature_id = rownames(Y),
    BPPARAM = BiocParallel::bpparam(), chunk_size = NULL,
    source_files = NULL, worker_init = NULL,
    marginal_test = NULL, marginal_args = list(), method = "REML",
    retain_smooth = FALSE,
    control = mgcv::gam.control(nthreads = 1L), ...,
    diagnostics = FALSE, retain_marginal = FALSE, offset = NULL) {
  if ("marginal" %in% names(list(...))) {
    stop("mgcvST.estimate() always runs the marginal score test; remove marginal.")
  }
  call <- match.call()
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
  if (inherits(G, "mgcvST_model")) {
    return(.mgcvst_estimate_model(
      Y = Y, model = G, feature_id = feature_id, BPPARAM = BPPARAM,
      chunk_size = chunk_size, source_files = source_files,
      worker_init = worker_init, marginal_test = marginal_test,
      marginal_args = marginal_args, method = method,
      retain_smooth = retain_smooth, control = control,
      gam_args = list(...), call = call,
      diagnostics = diagnostics, retain_marginal = retain_marginal, offset = offset
    ))
  }
  if (!is.null(offset)) stop("Additional offsets require a model prepared by mgcvST.set().")
  Y <- as.matrix(Y)
  storage.mode(Y) <- "double"
  if (length(dim(Y)) != 2L || !nrow(Y) || !ncol(Y) || any(!is.finite(Y))) {
    stop("Y must be a non-empty finite numeric feature-by-observation matrix.")
  }
  if (!is.list(G) || is.null(G$y) || is.null(G$family) || is.null(G$smooth)) {
    stop("G must be a reusable setup returned by mgcv::gam(..., fit = FALSE).")
  }
  G <- .mgcvst_mark_global(G)
  spec <- .mgcvst_external_spec(G)
  G <- mgcv::gam(spec$formula, data = spec$data, family = G$family,
                 fit = FALSE, na.action = stats::na.fail)
  G <- .mgcvst_mark_global(G)
  target_index <- which(vapply(G$smooth,
    function(s) identical(s$score.component, "global"), logical(1L)))
  null_formula <- .mgcvst_null_formula(spec$formula, target_index)
  null_G <- mgcv::gam(null_formula, data = spec$data, family = G$family,
                      fit = FALSE, na.action = stats::na.fail)
  attr(G, "training_X") <- as.matrix(G$X)
  attr(G, "null_spec") <- list(formula = null_formula, data = spec$data,
    response = names(G$mf)[attr(G$terms, "response")], X0 = as.matrix(null_G$X))
  attr(G, "full_spec") <- list(formula = spec$formula, data = spec$data,
    response = names(G$mf)[attr(G$terms, "response")])
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
  if (!is.logical(retain_smooth) || length(retain_smooth) != 1L ||
      is.na(retain_smooth)) {
    stop("retain_smooth must be TRUE or FALSE.")
  }
  if (!is.list(marginal_args) ||
      (length(marginal_args) &&
       (is.null(names(marginal_args)) || any(!nzchar(names(marginal_args)))))) {
    stop("marginal_args must be a named list.")
  }
  marginal_forbidden <- intersect(
    names(marginal_args),
    c("fit", "test.component", "n_threads")
  )
  if (length(marginal_forbidden)) {
    stop(
      "Do not supply these arguments through marginal_args: ",
      paste(marginal_forbidden, collapse = ", ")
    )
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
  # Poisson prescreen (R/family-prescreen.R). The knob lives in control and is
  # removed before control reaches mgcv::bam(), which rejects unknown entries.
  screen_threshold <- .mgcvst_prescreen_threshold(
    control[["poisson_screen_phi", exact = TRUE]]
  )
  control[["poisson_screen_phi"]] <- NULL
  prescreen <- .mgcvst_prescreen_route(
    Y, X = .mgcvst_prescreen_design(G), offset = G$offset,
    threshold = screen_threshold,
    active = isTRUE(tryCatch(
      identical(.working_family_id(G$family$family), "negative_binomial"),
      error = function(e) FALSE
    ))
  )
  chunk_id <- ceiling(seq_len(nrow(Y)) / chunk_size)
  ids <- split(seq_len(nrow(Y)), chunk_id)
  payload <- lapply(ids, function(i) list(
    index = i,
    feature_id = feature_id[i],
    Y = Y[i, , drop = FALSE],
    poisson = prescreen$poisson[i],
    prescreen_phi = prescreen$phi[i]
  ))
  attr(G, "training_X") <- as.matrix(G$X)
  family_raw <- serialize(G$family, NULL)
  # A prescreen-routed gene is fitted in the mgcv path with quasipoisson.
  # mgcv estimates the routed scale from the full spatial model. The INLA path
  # uses plain poisson because INLA has no quasi families.
  routed_family_raw <- serialize(stats::quasipoisson(link = "log"), NULL)
  worker_bundle <- .mgcvst_worker_bundle()
  fit_chunk <- get(".mgcvst_fit_chunk", envir = worker_bundle,
                   inherits = FALSE)

  t0 <- proc.time()[["elapsed"]]
  chunks <- BiocParallel::bplapply(
    payload, fit_chunk,
    G0 = G, family_raw = family_raw, method = method, control = control,
    gam_args = gam_args, source_files = source_files,
    worker_init = worker_init, init_key = init_key,
    marginal_test = marginal_test, marginal_args = marginal_args,
    retain_smooth = retain_smooth,
    diagnostics = diagnostics,
    retain_marginal = retain_marginal,
    routed_family_raw = routed_family_raw,
    BPPARAM = BPPARAM
  )
  elapsed <- proc.time()[["elapsed"]] - t0

  n <- ncol(Y)
  p <- nrow(Y)
  E <- matrix(NA_real_, nrow = n, ncol = p,
              dimnames = list(NULL, feature_id))
  V <- matrix(NA_real_, nrow = n, ncol = p,
              dimnames = list(NULL, feature_id))
  lambda <- dispersion <- stats::setNames(rep(NA_real_, p), feature_id)
  family_parameters <- stats::setNames(vector("list", p), feature_id)
  smooth_coefficients <- NULL
  fit_geometry <- NULL
  diagnostics <- vector("list", length(chunks))
  geometry <- NULL
  for (j in seq_along(chunks)) {
    z <- chunks[[j]]
    i <- z$index
    E[, i] <- z$working_error
    V[, i] <- z$working_variance
    lambda[i] <- z$lambda
    dispersion[i] <- z$dispersion
    family_parameters[i] <- z$family_parameters
    diagnostics[[j]] <- z$diagnostics
    if (!is.null(z$geometry)) {
      if (is.null(geometry)) geometry <- z$geometry
    }
    if (retain_smooth && !is.null(z$fit_geometry)) {
      if (is.null(fit_geometry)) {
        fit_geometry <- z$fit_geometry
        smooth_coefficients <- matrix(
          NA_real_, nrow = p, ncol = ncol(fit_geometry$B),
          dimnames = list(feature_id, NULL)
        )
      }
      smooth_coefficients[i, ] <- z$smooth_coefficients
    }
  }
  diagnostics <- do.call(rbind, diagnostics)
  diagnostics <- diagnostics[order(diagnostics$index), , drop = FALSE]
  rownames(diagnostics) <- NULL
  ans <- list(
      feature_id = feature_id,
      working_error = E,
      working_variance = V,
      dispersion = dispersion,
      lambda = lambda,
      family_parameters = family_parameters,
      geometry = geometry,
      row_id = if (is.null(geometry)) NULL else geometry$row_id,
      diagnostics = diagnostics,
      timing = list(
        elapsed = elapsed,
        workers = workers,
        chunks = length(chunks),
        chunk_size = chunk_size,
        backend = class(BPPARAM)[1L]
      ),
      source_files = source_files,
      retain_smooth = retain_smooth,
      test_engine = "spde",
      call = call
    )
  if (retain_marginal) {
    ans$marginal_data <- .mgcvst_collect_marginal(chunks, p, feature_id)
  }
  if (retain_smooth) {
    if (is.null(fit_geometry)) {
      ans$smooth_coefficients <- matrix(
        NA_real_, nrow = p, ncol = 0L,
        dimnames = list(feature_id, NULL)
      )
      ans$fit_basis <- NULL
      ans$fit_penalty <- NULL
      ans$basis_metadata <- NULL
    } else {
      ans$smooth_coefficients <- smooth_coefficients
      ans$fit_basis <- fit_geometry$B
      ans$fit_penalty <- fit_geometry$P
      ans$row_id <- fit_geometry$row_id
      ans$basis_metadata <- fit_geometry$metadata
    }
  }
  structure(ans, class = c("mgcvST_fit", "mgcvST"))
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
  cat("  backend:", x$timing$backend, "with", x$timing$workers, "worker(s)\n")
  cat("  elapsed seconds:", format(x$timing$elapsed), "\n")
  cat("  object size:", format(utils::object.size(x), units = "auto"), "\n")
  invisible(x)
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

# Construct the shared dense score factor once for every requested feature.
.mgcvst_legacy_shared_score_factor <- function(geometry) {
  Q <- as.matrix(geometry$Q)
  allow_psd <- isTRUE(geometry$score_precision_psd)
  if (allow_psd) {
    E <- CppMatrix::matrixEigen((Q + t(Q)) / 2)
    d <- as.numeric(E$values)
    tol <- sqrt(.Machine$double.eps) * max(1, max(abs(d)))
    if (min(d) < -tol) stop("The shared score precision is not positive semidefinite.")
    keep <- d > tol
    if (!any(keep)) stop("The shared score precision has no positive eigenvalues.")
    Qinvhalf <- CppMatrix::matrixMultiply(
      sweep(as.matrix(E$vectors[, keep, drop = FALSE]),
            2L, 1 / sqrt(d[keep]), "*"),
      t(as.matrix(E$vectors[, keep, drop = FALSE]))
    )
    T0 <- .magic_mm(geometry$B, Qinvhalf)
  } else {
    Q <- Matrix::forceSymmetric(Matrix::Matrix(Q, sparse = TRUE))
    R <- Matrix::chol(Q)
    Rinv <- Matrix::solve(R, Matrix::Diagonal(nrow(Q)))
    T0 <- .magic_mm(geometry$B, as.matrix(Rinv))
  }
  T0
}

# Construct each feature-level score summary once for the Liu pair engine.
.mgcvst_liu_summaries <- function(fitmgcvST, used, verbose, threads = 1L) {
  t0 <- proc.time()[["elapsed"]]
  geometry <- fitmgcvST$geometry
  field_scale <- .mgcvst_field_scale(fitmgcvST)
  T0 <- .mgcvst_legacy_shared_score_factor(geometry)

  a <- H <- vector("list", length(used))
  has_summary <- rep(FALSE, length(used))
  error_message <- rep(NA_character_, length(used))
  for (first in seq.int(1L, length(used), by = 32L)) {
    rows <- first:min(length(used), first + 31L)
    ids <- used[rows]
    units <- mgcvst_dense_score_batch_cpp(
      T0, fitmgcvST$working_variance[, ids, drop = FALSE],
      fitmgcvST$working_error[, ids, drop = FALSE], field_scale[ids],
      geometry$X, list(), threads
    )
    for (k in seq_along(rows)) {
      j <- rows[k]
      z <- units[[k]]
      if (!is.null(z$error)) {
        error_message[j] <- z$error
      } else {
        a[[j]] <- z$a
        H[[j]] <- z$H
        has_summary[j] <- TRUE
      }
    }
    if (verbose) {
      message("Constructed Liu summaries for ", max(rows), " of ", length(used),
              " features.")
    }
  }
  list(
    used = used, a = a, H = H, has_summary = has_summary,
    error_message = error_message,
    elapsed = proc.time()[["elapsed"]] - t0
  )
}

# Construct and write packed legacy score states for one feature-first chunk.
.mgcvst_legacy_score_unit_chunk <- function(payload, T0, X_fixed, threads = 1L) {
  .mgcvst_thread_limit()
  paths <- character(length(payload$feature))
  units <- mgcvst_dense_score_batch_cpp(
    T0, payload$working_variance, payload$working_error, payload$field_scale,
    X_fixed, list(), threads
  )
  for (j in seq_along(payload$feature)) {
    z <- units[[j]]
    state <- if (is.null(z$error)) {
      .mgcvst_pack_score_state(list(
        a = z$a, M = z$H, width = length(z$a)
      ))
    } else list(error = z$error)
    key <- as.character(payload$feature[j])
    path <- file.path(payload$directory, paste0("feature-", key, ".rds"))
    saveRDS(state, path, compress = FALSE)
    paths[j] <- path
  }
  stats::setNames(paths, as.character(payload$feature))
}

# Evaluate Liu-calibrated pairs from feature summaries with the fused
# double-precision C++ kernel (mgcvst_pair_liu_cpp): score, trace moments,
# and the log-space Liu tail in one call per chunk.
.mgcvst_liu_pairs <- function(index, pair_index, feature_id, summaries, threads,
                              chunk_size, verbose) {
  active <- which(summaries$has_summary)
  active_feature <- summaries$used[active]
  local <- matrix(match(index, active_feature), ncol = 2L)
  missing_summary <- !stats::complete.cases(local)

  k <- nrow(index)
  score <- information <- effective_rank <- p_value <- rep(NA_real_, k)
  log_p_two_sided <- log_p_positive <- log_p_negative <- rep(NA_real_, k)
  error_message <- rep(NA_character_, k)
  for (j in which(missing_summary)) {
    missing_feature <- index[j, ][!(index[j, ] %in% active_feature)]
    missing_local <- match(missing_feature, summaries$used)
    error_message[j] <- paste(
      paste0(feature_id[missing_feature], ": ",
             summaries$error_message[missing_local]),
      collapse = " | "
    )
  }
  rows <- which(!missing_summary)
  if (!length(rows)) {
    return(list(
      result = data.frame(
        pair_index = pair_index, score = score, information = information,
        effective_rank = effective_rank, p_value = p_value,
        log_p_two_sided = log_p_two_sided, log_p_positive = log_p_positive,
        log_p_negative = log_p_negative,
        error_message = error_message, stringsAsFactors = FALSE
      ),
      elapsed = 0
    ))
  }

  avec <- do.call(cbind, summaries$a[active])
  H <- summaries$H[active]
  t0 <- proc.time()[["elapsed"]]

  starts <- seq.int(1L, length(rows), by = chunk_size)
  for (b in seq_along(starts)) {
    z <- rows[starts[b]:min(length(rows), starts[b] + chunk_size - 1L)]
    pair <- local[z, , drop = FALSE]
    perm <- order(pair[, 1L])
    zp <- z[perm]
    left <- pair[perm, 1L]
    right <- pair[perm, 2L]
    res <- mgcvst_pair_liu_cpp(H, avec, left, right, threads)

    score[zp] <- res$score
    information[zp] <- res$information
    effective_rank[zp] <- res$effective_rank
    log_p_two_sided[zp] <- res$log_p_two_sided
    log_p_positive[zp] <- res$log_p_positive
    log_p_negative[zp] <- res$log_p_negative
    p_value[zp] <- exp(res$log_p_two_sided)

    bad_moments <- res$status == 1L
    bad_p <- res$status == 2L
    if (any(bad_moments)) {
      error_message[zp[bad_moments]] <-
        "Liu trace moments were non-finite or non-positive."
    }
    if (any(bad_p)) {
      error_message[zp[bad_p]] <-
        "Liu calibration returned an invalid p-value."
    }
    if (verbose && (b %% 10L == 0L || b == length(starts))) {
      message("Evaluated Liu block ", b, " of ", length(starts), ".")
    }
  }

  list(
    result = data.frame(
      pair_index = pair_index, score = score, information = information,
      effective_rank = effective_rank, p_value = p_value,
      log_p_two_sided = log_p_two_sided, log_p_positive = log_p_positive,
      log_p_negative = log_p_negative,
      error_message = error_message,
      stringsAsFactors = FALSE
    ),
    elapsed = proc.time()[["elapsed"]] - t0
  )
}

# Evaluate one pair chunk from compact feature summaries.
.mgcvst_test_chunk <- function(payload, geometry, calibration) {
  .mgcvst_thread_limit()
  pairs <- payload$pairs
  k <- nrow(pairs)
  used <- sort(unique(as.vector(pairs)))
  summaries <- vector("list", length(payload$feature_id))
  summary_error <- rep(NA_character_, length(payload$feature_id))
  for (i in used) {
    z <- tryCatch(
      {
        op <- rkhs_score_operator(
          geometry$B, geometry$Q, payload$working_variance[, i],
          geometry$X, payload$field_scale[i],
          allow_psd = isTRUE(geometry$score_precision_psd)
        )
        rkhs_score_summary(payload$working_error[, i], op)
      },
      error = function(e) e
    )
    if (inherits(z, "condition")) {
      summary_error[i] <- conditionMessage(z)
    } else {
      summaries[[i]] <- z
    }
  }
  signed_score <- statistic <- information <- effective_rank <-
    p_two_sided <- p_positive <- p_negative <- rep(NA_real_, k)
  error_message <- rep(NA_character_, k)
  for (j in seq_len(k)) {
    i1 <- pairs[j, 1L]
    i2 <- pairs[j, 2L]
    missing <- c(i1, i2)[vapply(summaries[c(i1, i2)], is.null, logical(1L))]
    if (length(missing)) {
      error_message[j] <- paste(
        paste0(payload$feature_id[missing], ": ", summary_error[missing]),
        collapse = " | "
      )
      next
    }
    S1 <- summaries[[i1]]
    S2 <- summaries[[i2]]
    z <- tryCatch(
      {
        U <- as.numeric(crossprod(S1$a, S2$a))
        cal <- rkhs_score_calibrate(U, S1$H, S2$H, method = calibration)
        list(U = U, calibration = cal)
      },
      error = function(e) e
    )
    if (inherits(z, "condition")) {
      error_message[j] <- conditionMessage(z)
      next
    }
    signed_score[j] <- z$U
    statistic[j] <- z$U^2
    information[j] <- z$calibration$information
    effective_rank[j] <- z$calibration$effective_rank
    p_two_sided[j] <- z$calibration$p_two_sided
    p_positive[j] <- z$calibration$p_positive
    p_negative[j] <- z$calibration$p_negative
  }
  data.frame(
    pair_index = payload$pair_index,
    feature1 = payload$feature_id[pairs[, 1L]],
    feature2 = payload$feature_id[pairs[, 2L]],
    signed_score = signed_score, statistic = statistic,
    information = information, effective_rank = effective_rank,
    p_two_sided = p_two_sided, p_positive = p_positive,
    p_negative = p_negative, error_message = error_message,
    stringsAsFactors = FALSE
  )
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
#' significance. When one feature or pair cannot be evaluated, its numerical
#' fields remain `NA` and its concrete error message is retained without
#' interrupting the other pairs. The
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
#' @param chunk_size Positive number of tested pairs per task or C++ block.
#'   Liu calibration defaults to 10,000 pairs per interruptible block.
#' @param threads Positive number of OpenMP threads for feature score preparation
#'   and the Liu pair kernel. `NULL` uses `bpworkers(BPPARAM)`. Model pair tasks
#'   and Davies calibration use `BPPARAM` after preparation has completed.
#' @param verbose Logical; report Liu summary and block progress.
#' @param cache_bytes Internal byte ceiling for resident Liu score states.
#'   `NULL` adapts to available system and job memory.
#' @param checkpoint_dir Optional Liu checkpoint directory; `NULL` uses
#'   temporary storage removed on exit.
#' @param resume Reuse compatible completed checkpoint entries.
#' @return A compact `mgcvST_test` object containing pair results, adjusted
#'   p-values, discovery/highlight/retention flags, threshold metadata, and
#'   timing metadata.
#' @keywords internal
.mgcvst_test_spde <- function(
    fitmgcvST, q.value = 0.05, FDR = TRUE, method = "BH",
    BPPARAM = BiocParallel::SerialParam(), ...,
    pairs = NULL, highlight = NULL,
    calibration = c("liu", "davies"),
    chunk_size = NULL,
    threads = NULL, verbose = FALSE, cache_bytes = NULL,
    checkpoint_dir = NULL, resume = TRUE) {
  if (!inherits(fitmgcvST, "mgcvST_fit")) {
    stop("fitmgcvST must be returned by mgcvST.estimate().")
  }
  if (is.null(fitmgcvST$geometry)) {
    stop("fitmgcvST has no feature geometry to test.")
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
  if (calibration != "liu" && (!is.null(cache_bytes) ||
      !is.null(checkpoint_dir) || !isTRUE(resume))) {
    stop("cache_bytes, checkpoint_dir and resume currently require calibration = 'liu'.")
  }
  if (calibration == "davies" &&
      !requireNamespace("CompQuadForm", quietly = TRUE)) {
    stop("calibration = 'davies' requires the optional CompQuadForm package.")
  }
  if (is.null(threads)) threads <- BiocParallel::bpworkers(BPPARAM)
  threads <- as.integer(threads)
  if (length(threads) != 1L || is.na(threads) || threads < 1L) {
    stop("threads must be one positive integer.")
  }
  .mgcvst_thread_limit()
  field_scale <- .mgcvst_field_scale(fitmgcvST)

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

  available <- .mgcvst_feature_available(fitmgcvST)
  i1 <- index[, 1L]
  i2 <- index[, 2L]
  pair_available <- available[i1] & available[i2]

  result <- data.frame(
    pair_index = seq_len(nrow(index)),
    feature1 = fitmgcvST$feature_id[i1],
    feature2 = fitmgcvST$feature_id[i2],
    signed_score = NA_real_,
    statistic = NA_real_,
    information = NA_real_,
    effective_rank = NA_real_,
    p_two_sided = NA_real_,
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
    error_message = NA_character_,
    stringsAsFactors = FALSE
  )
  unavailable_rows <- which(!pair_available)
  for (j in unavailable_rows) {
    missing_feature <- c(i1[j], i2[j])[!available[c(i1[j], i2[j])]]
    result$error_message[j] <- paste(
      paste0(
        fitmgcvST$feature_id[missing_feature], ": ",
        fitmgcvST$diagnostics$error_message[missing_feature]
      ),
      collapse = " | "
    )
  }

  tested_rows <- which(pair_available)
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
  pipeline <- NULL
  if (length(tested_rows)) {
    if (calibration == "liu") {
      evaluated <- .mgcvst_pair_pipeline(
        fitmgcvST,
        index[tested_rows, , drop = FALSE], tested_rows,
        threads, chunk_size, verbose, cache_bytes = cache_bytes,
        checkpoint_dir = checkpoint_dir, resume = resume
      )
      pipeline <- evaluated$metadata
      summary_elapsed <- pipeline$preparation_elapsed
      elapsed <- summary_elapsed + evaluated$elapsed
      evaluated <- evaluated$result
      target <- evaluated$pair_index
      result$signed_score[target] <- evaluated$score
      result$statistic[target] <- evaluated$score^2
      result$information[target] <- evaluated$information
      result$effective_rank[target] <- evaluated$effective_rank
      result$p_two_sided[target] <- evaluated$p_value
      result$error_message[target] <- evaluated$error_message
      chunks <- seq_len(pipeline$chunks)
    } else {
      used <- sort(unique(as.vector(index[tested_rows, , drop = FALSE])))
      T0 <- .mgcvst_legacy_shared_score_factor(fitmgcvST$geometry)
      state_dir <- .mgcvst_dense_temp_dir()
      on.exit(.mgcvst_dense_cleanup(state_dir), add = TRUE)
      feature_chunk_size <- 32L
      feature_groups <- split(
        used, ceiling(seq_along(used) / feature_chunk_size)
      )
      feature_payload <- lapply(feature_groups, function(features) list(
        feature = features,
        working_error = fitmgcvST$working_error[, features, drop = FALSE],
        working_variance = fitmgcvST$working_variance[, features, drop = FALSE],
        field_scale = field_scale[features], directory = state_dir
      ))
      t0 <- proc.time()[["elapsed"]]
      shard_paths <- lapply(
        feature_payload, .mgcvst_legacy_score_unit_chunk,
        T0 = T0, X_fixed = fitmgcvST$geometry$X, threads = threads
      )
      summary_elapsed <- proc.time()[["elapsed"]] - t0
      shard_paths <- unlist(shard_paths, use.names = FALSE)
      names(shard_paths) <- as.character(used)
      packed_bytes <- 8 * ncol(T0) * (ncol(T0) + 1) / 2
      unpacked_bytes <- 8 * ncol(T0)^2
      max_features <- max(1L, floor(
        (256 * 1024^2) / (packed_bytes + unpacked_bytes)
      ))
      chunks <- .mgcvst_dense_pair_groups(
        tested_rows, index, chunk_size, max_features
      )
      payload <- lapply(chunks, function(rows) {
        local_used <- unique(as.vector(t(index[rows, , drop = FALSE])))
        list(
          rows = rows,
          pairs = index[rows, , drop = FALSE],
          shards = shard_paths[as.character(local_used)]
        )
      })
      t0 <- proc.time()[["elapsed"]]
      evaluated <- BiocParallel::bplapply(
        payload, .mgcvst_dense_pair_chunk,
        calibration = calibration, BPPARAM = BPPARAM
      )
      elapsed <- summary_elapsed + proc.time()[["elapsed"]] - t0
      evaluated <- do.call(rbind, unlist(evaluated, recursive = FALSE))
      core <- c(
        "signed_score", "information", "effective_rank",
        "p_two_sided", "p_positive", "p_negative", "error_message"
      )
      target <- match(evaluated$pair_index, result$pair_index)
      result[target, core] <- evaluated[, core, drop = FALSE]
      result$statistic[target] <- evaluated$signed_score^2
    }
  }

  valid <- is.finite(result$p_two_sided) &
    result$p_two_sided >= 0 & result$p_two_sided <= 1
  result$p_positive[valid] <- ifelse(
    result$signed_score[valid] >= 0,
    result$p_two_sided[valid] / 2,
    1 - result$p_two_sided[valid] / 2
  )
  result$p_negative[valid] <- ifelse(
    result$signed_score[valid] <= 0,
    result$p_two_sided[valid] / 2,
    1 - result$p_two_sided[valid] / 2
  )
  if (FDR) {
    result$p_adjusted[valid] <- stats::p.adjust(
      result$p_two_sided[valid], method = method
    )
    result$p_positive_adjusted[valid] <- stats::p.adjust(
      result$p_positive[valid], method = method
    )
    result$p_negative_adjusted[valid] <- stats::p.adjust(
      result$p_negative[valid], method = method
    )
  } else {
    result$p_adjusted[valid] <- result$p_two_sided[valid]
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
    max(result$p_two_sided[result$discovered])
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
        pairs_tested = length(tested_rows),
        pairs_with_p_value = sum(valid),
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
        pair_pipeline = pipeline,
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
