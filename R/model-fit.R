# Expand a penalty to the coefficient dimension used by one smooth.
.mgcvst_expand_penalty <- function(S, n_coef, label) {
  S <- as.matrix(S)
  if (nrow(S) != ncol(S)) stop("Penalty for '", label, "' must be square.")
  if (nrow(S) == n_coef) return(S)
  if (!nrow(S) || n_coef %% nrow(S)) {
    stop("Penalty for '", label, "' is incompatible with its basis.")
  }
  kronecker(diag(n_coef %/% nrow(S)), S)
}

# Extract shared smooth geometry and fitted smoothing parameters.
.mgcvst_model_geometry <- function(fit) {
  L <- .gam_training_lpmatrix(fit)
  smooth_columns <- lapply(
    fit$smooth, function(s) seq.int(s$first.para, s$last.para)
  )
  all_smooth <- unique(unlist(smooth_columns, use.names = FALSE))
  parametric <- setdiff(seq_len(ncol(L)), all_smooth)
  smooth <- vector("list", length(fit$smooth))
  sp_value <- fit$sp
  if (!is.null(fit$full.sp) && length(fit$full.sp) == length(fit$sp)) {
    sp_value <- fit$full.sp
  }
  for (j in seq_along(fit$smooth)) {
    s <- fit$smooth[[j]]
    columns <- smooth_columns[[j]]
    fixed <- isTRUE(s$fixed) || is.null(s$S) || !length(s$S) ||
      (!is.null(s$rank) && isTRUE(all(s$rank == 0)))
    sp_index <- if (fixed) integer() else seq.int(s$first.sp, s$last.sp)
    if (!fixed && length(sp_index) != length(s$S)) {
      stop("Smooth '", s$label, "' has inconsistent penalty indexes.")
    }
    penalties <- if (fixed) list() else lapply(
      s$S, .mgcvst_expand_penalty,
      n_coef = length(columns), label = s$label
    )
    smooth[[j]] <- list(
      label = s$label,
      B = as.matrix(L[, columns, drop = FALSE]),
      penalties = penalties,
      sp_index = sp_index,
      fixed = fixed,
      score_component = s$score.component,
      columns = columns
    )
  }
  score_component <- vapply(
    smooth,
    function(s) if (is.null(s$score_component)) "" else s$score_component,
    character(1L)
  )
  marked <- which(nzchar(score_component))
  if (!length(marked) || sum(score_component == "global") != 1L ||
      sum(score_component == "local") > 1L) {
    stop("The model must mark one global and at most one local SPDE component.")
  }
  target <- stats::setNames(marked, score_component[marked])
  expected <- if ("local" %in% names(target)) c("global", "local") else "global"
  target <- target[expected]
  for (name in names(target)) {
    z <- smooth[[target[[name]]]]
    if (z$fixed || length(z$penalties) != 1L) {
      stop("A score SPDE component must have one fitted full-rank penalty.")
    }
  }
  offset <- fit$offset
  if (is.null(offset)) offset <- numeric(nrow(L))
  list(
    X = as.matrix(L[, parametric, drop = FALSE]),
    smooth = smooth,
    target = target,
    score_components = names(target),
    offset = as.numeric(offset),
    row_id = .mgcvst_row_id(fit, nrow(L)),
    sp = as.numeric(sp_value)
  )
}

# Check that two fitted genes use the same model geometry.
.mgcvst_model_geometry_equal <- function(x, y, tol) {
  if (!identical(x$row_id, y$row_id) ||
      !identical(x$score_components, y$score_components) ||
      !identical(x$target, y$target)) {
    stop("Feature fits do not share model rows or score components.")
  }
  if (!all(dim(x$X) == dim(y$X)) ||
      max(abs(x$X - y$X), 0) > tol * max(1, abs(x$X), abs(y$X))) {
    stop("Feature fits do not share the same linear design.")
  }
  if (length(x$offset) != length(y$offset) ||
      max(abs(x$offset - y$offset), 0) >
      tol * max(1, abs(x$offset), abs(y$offset))) {
    stop("Feature fits do not share the same offset.")
  }
  if (length(x$smooth) != length(y$smooth)) {
    stop("Feature fits do not share the same smoother collection.")
  }
  for (j in seq_along(x$smooth)) {
    a <- x$smooth[[j]]
    b <- y$smooth[[j]]
    if (!identical(a$label, b$label) || !identical(a$fixed, b$fixed) ||
        !identical(a$score_component, b$score_component) ||
        !identical(a$sp_index, b$sp_index) ||
        !all(dim(a$B) == dim(b$B)) ||
        max(abs(a$B - b$B), 0) > tol * max(1, abs(a$B), abs(b$B)) ||
        length(a$penalties) != length(b$penalties)) {
      stop("Feature fits do not share smoother geometry.")
    }
    for (k in seq_along(a$penalties)) {
      if (!all(dim(a$penalties[[k]]) == dim(b$penalties[[k]])) ||
          max(abs(a$penalties[[k]] - b$penalties[[k]]), 0) >
          tol * max(1, abs(a$penalties[[k]]), abs(b$penalties[[k]]))) {
        stop("Feature fits do not share smoother penalties.")
      }
    }
  }
  invisible(TRUE)
}

# Fit and reduce one feature under a model.set() setup.
.mgcvst_model_fit_one <- function(response, G0, family_raw, method, control,
                                  gam_args, marginal_test, marginal_args,
                                  retain_smooth) {
  G <- G0
  response_index <- attr(G$terms, "response")
  if (length(response_index) != 1L || response_index < 1L ||
      is.null(G$mf) || response_index > ncol(G$mf)) {
    stop("The reusable model does not contain a response bridge.")
  }
  G$y <- as.numeric(response)
  G$mf[[response_index]] <- as.numeric(response)
  G$family <- unserialize(family_raw)
  t0 <- proc.time()[["elapsed"]]
  fit <- do.call(
    mgcv::gam,
    c(list(G = G, method = method, control = control), gam_args)
  )
  fit_seconds <- proc.time()[["elapsed"]] - t0
  if (!isTRUE(fit$converged)) stop("The feature GAM did not converge.")
  W <- rkhs_extract_working_model(fit)
  geometry <- .mgcvst_model_geometry(fit)
  target_index <- unname(geometry$target[["global"]])
  marginal <- .mgcvst_marginal_score(
    fit, marginal_test, marginal_args, test_component = target_index
  )
  fit_summary <- summary(fit)
  criterion <- if (length(fit$gcv.ubre) == 1L) as.numeric(fit$gcv.ubre) else NA_real_
  criterion_name <- if (length(fit$gcv.ubre) == 1L) names(fit$gcv.ubre) else NA_character_
  coefficients <- NULL
  if (retain_smooth) {
    coefficients <- lapply(
      geometry$target,
      function(j) as.numeric(stats::coef(fit)[geometry$smooth[[j]]$columns])
    )
  }
  list(
    working_error = W$working_error,
    working_variance = W$working_variance,
    dispersion = W$dispersion,
    family = W$family,
    family_parameters = if (is.null(W$family_parameters)) numeric() else
      as.numeric(W$family_parameters),
    geometry = geometry,
    sp = geometry$sp,
    coefficients = coefficients,
    marginal_p_value = marginal$p_value,
    marginal_method = marginal$method,
    residual_df = as.numeric(fit$df.residual),
    criterion = criterion,
    criterion_name = criterion_name,
    converged = TRUE,
    outer_convergence = paste(fit$outer.info$conv, collapse = "; "),
    fit_seconds = fit_seconds,
    smooth_table = fit_summary$s.table
  )
}

# Parallel chunk for model.set() fits.
.mgcvst_model_fit_chunk <- function(payload, G0, family_raw, method, control,
                                    gam_args, source_files, worker_init,
                                    init_key, marginal_test, marginal_args,
                                    retain_smooth) {
  .mgcvst_worker_initialize(source_files, worker_init, init_key)
  out <- vector("list", length(payload$index))
  for (j in seq_along(payload$index)) {
    fit <- tryCatch(
      .mgcvst_model_fit_one(
        payload$Y[j, ], G0, family_raw, method, control, gam_args,
        marginal_test, marginal_args, retain_smooth
      ),
      error = function(e) e
    )
    if (inherits(fit, "condition")) {
      out[[j]] <- list(
        success = FALSE, error = .mgcvst_condition(fit),
        index = payload$index[j], feature_id = payload$feature_id[j]
      )
    } else {
      fit$success <- TRUE
      fit$index <- payload$index[j]
      fit$feature_id <- payload$feature_id[j]
      out[[j]] <- fit
    }
  }
  out
}

# Fit all features for a model.set() object.
.mgcvst_estimate_model <- function(
    Y, model, feature_id, BPPARAM, chunk_size, source_files, worker_init,
    marginal_test, marginal_args, method, retain_smooth, control,
    geometry_tol, gam_args, call) {
  Y <- as.matrix(Y)
  storage.mode(Y) <- "double"
  if (length(dim(Y)) != 2L || !nrow(Y) || !ncol(Y) || any(!is.finite(Y))) {
    stop("Y must be a non-empty finite numeric feature-by-observation matrix.")
  }
  if (ncol(Y) != length(model$G$y)) {
    stop("ncol(Y) must equal the number of observations in model.")
  }
  if (is.null(feature_id)) feature_id <- rownames(Y)
  if (is.null(feature_id)) feature_id <- as.character(seq_len(nrow(Y)))
  feature_id <- as.character(feature_id)
  if (length(feature_id) != nrow(Y) || anyNA(feature_id) ||
      any(!nzchar(feature_id)) || anyDuplicated(feature_id)) {
    stop("feature_id must contain one unique non-empty identifier per feature.")
  }
  if (!inherits(BPPARAM, "BiocParallelParam")) {
    stop("BPPARAM must inherit from 'BiocParallelParam'.")
  }
  if (!is.list(control)) stop("control must be returned by mgcv::gam.control().")
  control$nthreads <- 1L
  control$ncv.threads <- 1L
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
  ids <- split(seq_len(nrow(Y)), ceiling(seq_len(nrow(Y)) / chunk_size))
  payload <- lapply(ids, function(i) list(
    index = i, feature_id = feature_id[i], Y = Y[i, , drop = FALSE]
  ))
  family_raw <- serialize(model$G$family, NULL)
  worker_bundle <- .mgcvst_worker_bundle()
  fit_chunk <- get(".mgcvst_model_fit_chunk", envir = worker_bundle,
                   inherits = FALSE)
  t0 <- proc.time()[["elapsed"]]
  chunks <- BiocParallel::bplapply(
    payload, fit_chunk,
    G0 = model$G, family_raw = family_raw, method = method,
    control = control, gam_args = gam_args, source_files = source_files,
    worker_init = worker_init, init_key = init_key,
    marginal_test = marginal_test, marginal_args = marginal_args,
    retain_smooth = retain_smooth, BPPARAM = BPPARAM
  )
  elapsed <- proc.time()[["elapsed"]] - t0
  fits <- unlist(chunks, recursive = FALSE)
  fits <- fits[order(vapply(fits, `[[`, integer(1L), "index"))]
  success <- vapply(fits, `[[`, logical(1L), "success")
  good <- which(success)
  geometry <- NULL
  if (length(good)) {
    geometry <- fits[[good[1L]]]$geometry
    for (j in good[-1L]) {
      .mgcvst_model_geometry_equal(geometry, fits[[j]]$geometry, geometry_tol)
    }
  }
  n <- ncol(Y)
  p <- nrow(Y)
  E <- V <- matrix(NA_real_, n, p, dimnames = list(NULL, feature_id))
  dispersion <- stats::setNames(rep(NA_real_, p), feature_id)
  family_parameters <- stats::setNames(vector("list", p), feature_id)
  n_sp <- if (is.null(geometry)) length(model$G$sp) else length(geometry$sp)
  smoothing_parameters <- matrix(
    NA_real_, p, n_sp,
    dimnames = list(feature_id, names(model$G$sp))
  )
  diagnostics <- data.frame(
    index = seq_len(p), feature_id = feature_id, success = success,
    converged = FALSE, marginal_p_value = NA_real_,
    marginal_method = NA_character_, residual_df = NA_real_,
    criterion = NA_real_, criterion_name = NA_character_,
    fit_seconds = NA_real_, outer_convergence = NA_character_,
    error_class = NA_character_, error_message = NA_character_,
    error_call = NA_character_, stringsAsFactors = FALSE
  )
  coefficient <- if (retain_smooth && !is.null(geometry)) {
    lapply(geometry$score_components, function(name) {
      q <- ncol(geometry$smooth[[geometry$target[[name]]]]$B)
      matrix(NA_real_, p, q, dimnames = list(feature_id, NULL))
    })
  } else {
    NULL
  }
  if (!is.null(coefficient)) names(coefficient) <- geometry$score_components
  for (j in seq_len(p)) {
    z <- fits[[j]]
    if (!z$success) {
      diagnostics$error_class[j] <- z$error$class
      diagnostics$error_message[j] <- z$error$message
      diagnostics$error_call[j] <- z$error$call
      next
    }
    E[, j] <- z$working_error
    V[, j] <- z$working_variance
    dispersion[j] <- z$dispersion
    family_parameters[[j]] <- z$family_parameters
    smoothing_parameters[j, ] <- z$sp
    diagnostics$converged[j] <- z$converged
    diagnostics$marginal_p_value[j] <- z$marginal_p_value
    diagnostics$marginal_method[j] <- z$marginal_method
    diagnostics$residual_df[j] <- z$residual_df
    diagnostics$criterion[j] <- z$criterion
    diagnostics$criterion_name[j] <- z$criterion_name
    diagnostics$fit_seconds[j] <- z$fit_seconds
    diagnostics$outer_convergence[j] <- z$outer_convergence
    if (!is.null(coefficient)) {
      for (name in names(coefficient)) coefficient[[name]][j, ] <- z$coefficients[[name]]
    }
  }
  failures <- diagnostics[!success, c(
    "index", "feature_id", "error_class", "error_message", "error_call"
  ), drop = FALSE]
  target_lambda <- if (!is.null(geometry)) {
    vapply(geometry$target, function(j) geometry$smooth[[j]]$sp_index, integer(1L))
  } else {
    integer()
  }
  component_lambda <- if (length(target_lambda)) {
    smoothing_parameters[, target_lambda, drop = FALSE]
  } else {
    matrix(NA_real_, p, 0L)
  }
  colnames(component_lambda) <- names(target_lambda)
  structure(
    list(
      feature_id = feature_id,
      working_error = E,
      working_variance = V,
      dispersion = dispersion,
      lambda = if ("global" %in% colnames(component_lambda))
        component_lambda[, "global"] else rep(NA_real_, p),
      component_lambda = component_lambda,
      smoothing_parameters = smoothing_parameters,
      family_parameters = family_parameters,
      geometry = geometry,
      row_id = if (is.null(geometry)) NULL else geometry$row_id,
      offset = if (is.null(geometry)) model$offset else geometry$offset,
      linear_design = if (is.null(geometry)) NULL else geometry$X,
      score_components = if (is.null(geometry)) model$components else geometry$score_components,
      model_setting = model$setting,
      model = model,
      diagnostics = diagnostics,
      failures = failures,
      timing = list(
        elapsed = elapsed, workers = workers, chunks = length(chunks),
        chunk_size = chunk_size, backend = class(BPPARAM)[1L]
      ),
      geometry_tol = geometry_tol,
      smooth_coefficients = coefficient,
      retain_smooth = retain_smooth,
      test_engine = if (length(model$components) == 2L) "global_local" else "single_model",
      call = call
    ),
    class = c("mgcvST_model_fit", "mgcvST_fit", "mgcvST")
  )
}
