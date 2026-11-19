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
.mgcvst_model_geometry <- function(fit, L = NULL) {
  if (is.null(L)) L <- .gam_training_lpmatrix(fit)
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
  if (length(marked) != 1L || sum(score_component == "global") != 1L) {
    stop("The model must mark exactly one global SPDE score component. The ",
         "second 'local' geographic process was removed from mgcvST.")
  }
  target <- stats::setNames(marked, "global")
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

# Fit and reduce one feature under an mgcvST.set() model.
.mgcvst_model_fit_one <- function(response, G0, family_raw, control,
                                  gam_args, retain_smooth, diagnostics = TRUE,
                                  geometry_cache, offset = NULL,
                                  poisson = FALSE, routed_family_raw = NULL) {
  full_spec <- attr(G0, "full_spec")
  response_index <- attr(G0$terms, "response")
  if (length(response_index) != 1L || response_index < 1L ||
       is.null(G0$mf) || response_index > ncol(G0$mf)) {
    stop("The reusable model does not contain a response bridge.")
  }
  full_data <- full_spec$data
  full_data[[full_spec$response]] <- as.numeric(response)
  # Poisson prescreen routing (R/family-prescreen.R); FALSE keeps the model
  # family. TRUE selects the Poisson routing family.
  family <- mgcv::fix.family.ls(
    unserialize(if (isTRUE(poisson)) routed_family_raw else family_raw)
  )
  t0 <- proc.time()[["elapsed"]]
  fit <- do.call(
    mgcv::bam,
    c(list(formula = full_spec$formula, data = full_data, family = family,
           offset = offset, method = "fREML", discrete = TRUE, nthreads = 1L,
           control = control), gam_args)
  )
  fit_seconds <- proc.time()[["elapsed"]] - t0
  W <- rkhs_extract_working_model(fit)
  geometry <- .mgcvst_cached_model_geometry(fit, geometry_cache)
  nuisance <- .mgcvst_nuisance_state(fit, geometry, geometry_cache)
  if (!is.null(nuisance$error)) {
    stop("nuisance covariance unavailable: ", nuisance$error)
  }
  geometry$nuisance_columns <- nuisance$columns
  geometry$nuisance_design <- nuisance$design
  geometry$nuisance_projection <- "conditional_Vp_block"
  fit_summary <- if (diagnostics) summary(fit) else NULL
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
    gam = fit,
    working_error = W$working_error,
    working_variance = W$working_variance,
    dispersion = W$dispersion,
    family = W$family,
    family_parameters = if (is.null(W$family_parameters)) numeric() else
      as.numeric(W$family_parameters),
    geometry = geometry,
    nuisance_covariance = nuisance$covariance,
    sp = geometry$sp,
    coefficients = coefficients,
    residual_df = as.numeric(fit$df.residual),
    criterion = criterion,
    criterion_name = criterion_name,
    family_used = W$family,
    converged = isTRUE(fit$converged),
    outer_convergence = paste(fit$outer.info$conv, collapse = "; "),
    fit_seconds = fit_seconds,
    smooth_table = fit_summary$s.table
  )
}

# Step 1 worker for mgcvST.set() fits: the null fit and the corrected marginal
# score test of every feature of a chunk. Only the marginal result is
# returned: no working
# vector of any null fit leaves the worker. The chunk is saved by the worker
# when the payload names a checkpoint file.
.mgcvst_null_chunk <- function(payload, G0, family_raw, control, gam_args,
                               routed_family_raw = NULL, scan = FALSE) {
  .mgcvst_worker_initialize()
  out <- vector("list", length(payload$index))
  target_index <- which(vapply(
    G0$smooth, function(s) identical(s$score.component, "global"), logical(1L)
  ))
  null_setup <- .mgcvst_null_score_setup(G0, target_index, attr(G0, "null_spec"))
  # spde.scan (R/spde-scan.R): the Stage 1 of the chunk, as inlaST.estimate().
  if (isTRUE(scan)) {
    out <- get(".mgcvst_scan_chunk", envir = asNamespace("spdeST"))(
      payload, G0, target_index, null_setup$X0)
    .mgcvst_chunk_save(payload$chunk_file, payload$chunk_key, out)
    return(out)
  }
  for (j in seq_along(payload$index)) {
    response <- as.numeric(payload$Y[j, ])
    feature_offset <- if (is.matrix(payload$offset)) payload$offset[j, ] else payload$offset
    null_data <- null_setup$spec$data
    null_data[[null_setup$spec$response]] <- response
    marginal_result <- tryCatch(
      {
        null_fit <- .mgcvst_fit_null(
          null_setup, null_data,
          unserialize(if (isTRUE(payload$poisson[j])) routed_family_raw else family_raw),
          feature_offset, control, gam_args
        )
        .mgcvst_marginal_score(null_fit, null_setup)
      },
      error = function(e) e
    )
    failed <- inherits(marginal_result, "condition")
    out[[j]] <- list(
      marginal_p_value = if (failed) NA_real_ else marginal_result$p_value,
      marginal_requested_method = if (failed) NA_character_ else
        marginal_result$requested_method,
      marginal_method = if (failed) NA_character_ else marginal_result$method,
      marginal_fallback = if (failed) NA else marginal_result$fallback,
      marginal_error = if (failed) .mgcvst_condition(marginal_result) else NULL
    )
  }
  .mgcvst_chunk_save(payload$chunk_file, payload$chunk_key, out)
  out
}

# Step 2 worker for mgcvST.set() fits: the spatial fit of every feature of a
# chunk, reduced to the compact working model. Chunk attributes carry the shared
# prediction geometry.
.mgcvst_spatial_chunk <- function(payload, G0, family_raw, control,
                                  gam_args, retain_smooth,
                                  diagnostics = TRUE, geometry_seed,
                                  routed_family_raw = NULL) {
  .mgcvst_worker_initialize()
  out <- vector("list", length(payload$index))
  geometry_cache <- list2env(geometry_seed, parent = emptyenv())
  shared_geometry <- NULL
  for (j in seq_along(payload$index)) {
    feature_offset <- if (is.matrix(payload$offset)) payload$offset[j, ] else payload$offset
    fit <- tryCatch(
      .mgcvst_model_fit_one(
        payload$Y[j, ], G0, family_raw, control, gam_args,
        retain_smooth, diagnostics = diagnostics, geometry_cache = geometry_cache,
        offset = feature_offset,
        poisson = isTRUE(payload$poisson[j]),
        routed_family_raw = routed_family_raw
      ),
      error = function(e) e
    )
    if (inherits(fit, "condition")) {
      out[[j]] <- list(
        error = .mgcvst_condition(fit),
        index = payload$index[j], feature_id = payload$feature_id[j]
      )
    } else {
      fit$gam <- NULL
      if (is.null(shared_geometry)) shared_geometry <- fit$geometry
      fit$geometry <- NULL
      fit$index <- payload$index[j]
      fit$feature_id <- payload$feature_id[j]
      out[[j]] <- fit
    }
  }
  attr(out, "model_geometry") <- shared_geometry
  .mgcvst_chunk_save(payload$chunk_file, payload$chunk_key, out)
  out
}

# Checkpoint signature of an mgcv estimation: the design, penalties, family,
# offsets, controls and options that the chunk results depend on.
.mgcvst_model_signature <- function(model, control, gam_args,
                                    offset, feature_id, retain_smooth,
                                    diagnostics) {
  # The family enters by name, link and fixed parameters: serialized closures
  # differ between sessions (byte code, environments) and would make a
  # checkpoint look foreign.
  family <- model$G$family
  family <- list(
    family = family$family, link = family$link,
    theta = tryCatch(as.numeric(family$getTheta()), error = function(e) NULL)
  )
  smooth <- lapply(model$G$smooth, function(s) {
    s[intersect(c("S", "sp", "first.para", "last.para", "label", "fixed"), names(s))]
  })
  .mgcvst_pair_input_hash(list(
    format = .mgcvst_fit_format, X = model$G$X, offset = model$G$offset,
    smooth = smooth, null_X = model$null_X, family = family,
    formulas = c(paste(deparse(model$full_formula), collapse = ""),
                 paste(deparse(model$null_formula), collapse = "")),
    control = control, gam_args = gam_args,
    extra_offset = offset, feature_id = feature_id,
    retain_smooth = retain_smooth, diagnostics = diagnostics
  ))
}

# Stage 1 of an mgcvST.set() model: the null fit and Stage 1 p-value of every
# feature. The spatial fits are added by .mgcvst_apply_spatial().
.mgcvst_estimate_model <- function(
    Y, model, feature_id, BPPARAM, chunk_size, retain_smooth, control,
    gam_args, call, diagnostics = TRUE, offset = NULL,
    adjust = "BY", q.value = 0.05, checkpoint_dir = NULL, resume = TRUE) {
  if (!is.matrix(Y)) Y <- as.matrix(Y)
  if (length(dim(Y)) != 2L || !nrow(Y) || !ncol(Y)) {
    stop("Y must be a non-empty finite numeric feature-by-observation matrix.")
  }
  if (ncol(Y) != length(model$G$y)) {
    stop("ncol(Y) must equal the number of observations in model.")
  }
  # Checked in blocks of rows: no temporary of the size of Y, no copy of a
  # double or integer matrix.
  Y <- .mgcvst_check_response_matrix(Y)
  chunk_size <- .mgcvst_check_chunk_size(chunk_size)
  if (!is.null(offset)) {
    if (!is.numeric(offset) || any(!is.finite(offset)) ||
        (is.null(dim(offset)) && length(offset) != ncol(Y)) ||
        (!is.null(dim(offset)) && (!is.matrix(offset) || !identical(dim(offset), dim(Y))))) {
      stop("offset must be finite numeric: an observation-length vector or a matrix matching Y.")
    }
    if (is.matrix(offset) &&
        ((!is.null(rownames(offset)) && !identical(rownames(offset), rownames(Y))) ||
         (!is.null(colnames(offset)) && !identical(colnames(offset), colnames(Y))))) {
      stop("Named offset rows and columns must match Y exactly.")
    }
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
  p <- nrow(Y)
  n <- ncol(Y)
  .mgcvst_refuse_screen_control(control)
  control$nthreads <- 1L
  control$ncv.threads <- 1L
  forbidden <- intersect(names(gam_args),
    c("G", "family", "method", "control", "formula", "data", "weights",
      "subset", "na.action", "knots", "paraPen", "H"))
  if (length(forbidden)) {
    stop("Do not supply these arguments through ...: ", paste(forbidden, collapse = ", "))
  }
  family_raw <- serialize(model$G$family, NULL)
  # Routing target: Poisson (dispersion 1), as in inlaST.estimate().
  routed_family_raw <- serialize(stats::poisson(link = "log"), NULL)
  base_offset <- if (is.null(model$G$offset)) numeric(n) else
    as.numeric(model$G$offset)
  prescreen <- .mgcvst_prescreen_route(
    Y, X = .mgcvst_prescreen_design(model$G),
    offset = if (is.null(offset)) base_offset else if (is.matrix(offset))
      sweep(offset, 2L, base_offset, "+") else base_offset + offset,
    active = isTRUE(tryCatch(
      identical(.working_family_id(model$G$family$family), "negative_binomial"),
      error = function(e) FALSE
    ))
  )
  t0 <- proc.time()[["elapsed"]]
  y_digest <- .mgcvst_row_digests(Y)
  signature <- .mgcvst_model_signature(
    model, control, gam_args, offset, feature_id, retain_smooth, diagnostics
  )
  # spde.scan: batched maximum-likelihood nulls when the null has no smooth.
  scan <- .mgcvst_scan_eligible(
    model$G, which(vapply(model$G$smooth, function(s)
      identical(s$score.component, "global"), logical(1L))),
    model$null_formula, model$G$family, gam_args)
  if (scan) signature <- .mgcvst_pair_input_hash(list(signature, null = "spde.scan"))
  store <- .mgcvst_chunk_store(checkpoint_dir, "mgcv", signature, resume)

  # Step 1: the null fit and Stage 1 p-value of every feature.
  worker_bundle <- .mgcvst_worker_bundle()
  null_chunk <- get(".mgcvst_null_chunk", envir = worker_bundle, inherits = FALSE)
  groups <- .mgcvst_feature_chunks(seq_len(p), chunk_size, BPPARAM,
                                   checkpoint = !is.null(store))
  # A payload is built only when its chunk is computed and a worker is free.
  make_payload <- function(i) list(
    index = i, feature_id = feature_id[i], Y = .mgcvst_double_rows(Y, i),
    offset = if (is.matrix(offset)) offset[i, , drop = FALSE] else offset,
    poisson = prescreen$poisson[i]
  )
  null_t0 <- proc.time()[["elapsed"]]
  step1 <- .mgcvst_run_chunks(
    groups, make_payload, "null", store, y_digest, prescreen$poisson, BPPARAM,
    null_chunk,
    G0 = structure(model$G,
      null_spec = list(formula = model$null_formula, data = model$null_data,
        response = model$null_response, X0 = model$null_X)),
    family_raw = family_raw, control = control, gam_args = gam_args,
    routed_family_raw = routed_family_raw,
    scan = scan
  )
  null_elapsed <- proc.time()[["elapsed"]] - null_t0
  nulls <- unlist(step1$results, recursive = FALSE)
  marginal_p <- vapply(nulls, `[[`, numeric(1L), "marginal_p_value")
  marginal_q <- .mgcvst_stage1_q(marginal_p, adjust)

  family_id <- tryCatch(.working_family_id(model$G$family$family),
                        error = function(e) NA_character_)
  table <- data.frame(
    index = seq_len(p), feature_id = feature_id,
    converged = FALSE, marginal_p_value = marginal_p,
    marginal_q_value = marginal_q,
    marginal_requested_method = vapply(nulls, `[[`, character(1L),
                                       "marginal_requested_method"),
    marginal_method = vapply(nulls, `[[`, character(1L), "marginal_method"),
    marginal_fallback = vapply(nulls, `[[`, logical(1L), "marginal_fallback"),
    residual_df = NA_real_,
    criterion = NA_real_, criterion_name = NA_character_,
    fit_seconds = NA_real_, outer_convergence = NA_character_,
    error_class = NA_character_, error_message = NA_character_,
    error_call = NA_character_,
    spatial_selected = FALSE, spatial_fitted = FALSE,
    prescreen_phi = prescreen$phi,
    family_used = ifelse(prescreen$poisson, "poisson", family_id),
    stringsAsFactors = FALSE
  )
  for (j in seq_len(p)) {
    if (is.null(nulls[[j]]$marginal_error)) next
    table$error_class[j] <- nulls[[j]]$marginal_error$class
    table$error_message[j] <- nulls[[j]]$marginal_error$message
    table$error_call[j] <- nulls[[j]]$marginal_error$call
  }
  rm(nulls)

  total_offset <- if (is.null(offset)) model$offset else if (is.matrix(offset))
    sweep(offset, 2L, model$offset, "+") else model$offset + offset
  ans <- structure(
    list(
      feature_id = feature_id,
      working_error = matrix(NA_real_, n, p, dimnames = list(NULL, feature_id)),
      working_variance = matrix(NA_real_, n, p, dimnames = list(NULL, feature_id)),
      dispersion = stats::setNames(rep(NA_real_, p), feature_id),
      lambda = stats::setNames(rep(NA_real_, p), feature_id),
      component_lambda = matrix(NA_real_, p, 0L),
      smoothing_parameters = matrix(
        NA_real_, p, length(model$G$sp),
        dimnames = list(feature_id, names(model$G$sp))
      ),
      family_parameters = stats::setNames(vector("list", p), feature_id),
      geometry = NULL,
      nuisance_covariance = stats::setNames(vector("list", p), feature_id),
      row_id = NULL,
      offset = total_offset,
      linear_design = NULL,
      score_components = model$components,
      model_setting = model$setting,
      model = model,
      diagnostics = table,
      timing = list(
        elapsed = NA_real_, null_fit_elapsed = null_elapsed, fit_elapsed = 0,
        chunks = length(groups), chunk_size = chunk_size,
        backend = class(BPPARAM)[1L],
        workers = max(1L, min(length(groups), BiocParallel::bpworkers(BPPARAM))),
        resumed_null_chunks = step1$resumed, null_payloads = step1$built,
        spatial_chunks = 0L, resumed_spatial_chunks = 0L, spatial_payloads = 0L
      ),
      smooth_coefficients = NULL,
      retain_smooth = retain_smooth,
      test_engine = "single_model",
      y_digest = y_digest,
      n_observation = n,
      format = .mgcvst_fit_format,
      signature = signature,
      stage1 = list(adjust = adjust, q.value = q.value,
                    null_engine = if (scan) "spde.scan" else "gam"),
      estimation_context = list(
        family_raw = family_raw, routed_family_raw = routed_family_raw,
        control = control, gam_args = gam_args,
        diagnostics = diagnostics, offset = offset,
        poisson = prescreen$poisson
      ),
      call = call
    ),
    class = c("mgcvST_model_fit", "mgcvST_fit", "spdeST")
  )
  kappa <- .spde_kappa_fields(model)
  ans[names(kappa)] <- kappa
  ans$timing$setup <- model$timing
  rm(step1)
  ans$timing$elapsed <- proc.time()[["elapsed"]] - t0
  ans
}

# The geometry-dependent fields of a fit, set when the first spatial fit
# establishes the shared geometry.
.mgcvst_fit_set_geometry <- function(fit, geometry) {
  model <- fit$model
  p <- length(fit$feature_id)
  if (ncol(fit$smoothing_parameters) != length(geometry$sp)) {
    fit$smoothing_parameters <- matrix(NA_real_, p, length(geometry$sp),
                                       dimnames = list(fit$feature_id, NULL))
  }
  fit$geometry <- geometry
  fit$row_id <- geometry$row_id
  fit$linear_design <- geometry$X
  fit$score_components <- geometry$score_components
  fit$geometry$offset <- model$offset
  if (fit$retain_smooth) {
    coefficient <- lapply(geometry$score_components, function(name) {
      q <- ncol(geometry$smooth[[geometry$target[[name]]]]$B)
      matrix(NA_real_, p, q, dimnames = list(fit$feature_id, NULL))
    })
    names(coefficient) <- geometry$score_components
    fit$smooth_coefficients <- coefficient
  }
  fit
}

# Step 2: fit the spatial model of the features `index` that have none yet and
# write the compact working models into `fit`.
.mgcvst_apply_spatial <- function(fit, Y, index, BPPARAM, chunk_size, store) {
  context <- fit$estimation_context
  model <- fit$model
  index <- setdiff(index, which(fit$diagnostics$spatial_fitted))
  if (!length(index)) return(fit)
  offset <- context$offset
  worker_bundle <- .mgcvst_worker_bundle()
  fit_chunk <- get(".mgcvst_spatial_chunk", envir = worker_bundle, inherits = FALSE)
  fit_args <- list(
    G0 = structure(model$G,
      null_spec = list(formula = model$null_formula, data = model$null_data,
        response = model$null_response, X0 = model$null_X),
      full_spec = list(formula = model$full_formula, data = model$full_data,
        response = model$null_response)),
    family_raw = context$family_raw,
    control = context$control, gam_args = context$gam_args,
    retain_smooth = fit$retain_smooth,
    diagnostics = context$diagnostics,
    routed_family_raw = context$routed_family_raw
  )
  payload_for <- function(i) list(
    index = i, feature_id = fit$feature_id[i], Y = .mgcvst_double_rows(Y, i),
    offset = if (is.matrix(offset)) offset[i, , drop = FALSE] else offset,
    poisson = context$poisson[i]
  )
  run <- function(chunk_groups, param, ...) {
    do.call(.mgcvst_run_chunks, c(
      list(chunk_groups, payload_for, "spatial", store, fit$y_digest,
           context$poisson, param, fit_chunk),
      fit_args, list(...)))
  }
  fit_t0 <- proc.time()[["elapsed"]]
  # The shared prediction geometry of the model seeds every worker.
  seed <- list(frozen = TRUE, L = model$L, geometry = model$geometry)
  groups <- .mgcvst_feature_chunks(index, chunk_size, BPPARAM,
                                   checkpoint = !is.null(store))
  tail <- run(groups, BPPARAM, geometry_seed = seed)
  chunks <- tail$results
  fit_elapsed <- proc.time()[["elapsed"]] - fit_t0

  if (is.null(fit$geometry)) {
    for (chunk in chunks) {
      geometry <- attr(chunk, "model_geometry")
      if (!is.null(geometry)) {
        fit <- .mgcvst_fit_set_geometry(fit, geometry)
        break
      }
    }
  }
  fits <- unlist(chunks, recursive = FALSE)
  rm(chunks)
  # The n-by-p working matrices are filled in place from local variables.
  E <- fit$working_error
  V <- fit$working_variance
  fit["working_error"] <- list(NULL)
  fit["working_variance"] <- list(NULL)
  table <- fit$diagnostics
  for (z in fits) {
    j <- z$index
    if (!is.null(z$error)) {
      table$error_class[j] <- z$error$class
      table$error_message[j] <- z$error$message
      table$error_call[j] <- z$error$call
      next
    }
    E[, j] <- z$working_error
    V[, j] <- z$working_variance
    fit$dispersion[j] <- z$dispersion
    fit$family_parameters[[j]] <- z$family_parameters
    fit$smoothing_parameters[j, ] <- z$sp
    fit$nuisance_covariance[[j]] <- z$nuisance_covariance
    table$converged[j] <- z$converged
    table$residual_df[j] <- z$residual_df
    table$criterion[j] <- z$criterion
    table$criterion_name[j] <- z$criterion_name
    table$fit_seconds[j] <- z$fit_seconds
    table$outer_convergence[j] <- z$outer_convergence
    table$family_used[j] <- z$family_used
    table$spatial_fitted[j] <- TRUE
    if (!is.null(fit$smooth_coefficients)) {
      for (name in names(fit$smooth_coefficients)) {
        fit$smooth_coefficients[[name]][j, ] <- z$coefficients[[name]]
      }
    }
  }
  fit$working_error <- E
  fit$working_variance <- V
  fit$diagnostics <- table
  geometry <- fit$geometry
  target_lambda <- if (!is.null(geometry)) {
    vapply(geometry$target, function(j) geometry$smooth[[j]]$sp_index, integer(1L))
  } else integer()
  component_lambda <- if (length(target_lambda)) {
    fit$smoothing_parameters[, target_lambda, drop = FALSE]
  } else {
    matrix(NA_real_, length(fit$feature_id), 0L)
  }
  colnames(component_lambda) <- names(target_lambda)
  fit$component_lambda <- component_lambda
  fit$lambda <- if ("global" %in% colnames(component_lambda))
    component_lambda[, "global"] else rep(NA_real_, length(fit$feature_id))
  fit$timing$fit_elapsed <- fit$timing$fit_elapsed + fit_elapsed
  fit$timing$spatial_chunks <- fit$timing$spatial_chunks + tail$chunks
  fit$timing$resumed_spatial_chunks <- fit$timing$resumed_spatial_chunks +
    tail$resumed
  fit$timing$spatial_payloads <- fit$timing$spatial_payloads + tail$built
  fit
}

#' Add spatial models after step 1 of mgcv estimation
#'
#' Fits the spatial model of further features of a fit returned by
#' [mgcvST.estimate()] and returns the updated fit; the supplied fit is not
#' changed. Features that already have a spatial model are skipped. The
#' responses `Y` must be those of step 1, which is checked against the digests
#' stored in the fit. The model, controls, offset, family routing and retention
#' options of the original call are reused.
#'
#' @param fitmgcvST A fit returned by [mgcvST.estimate()] (version 0.0.1.9032
#'   or later).
#' @param Y The feature-by-observation matrix given to [mgcvST.estimate()].
#' @param features The features to add: `"discoveries"` (the default; Stage 1
#'   q-value at most `q.value` under `adjust`, computed from the stored Stage 1
#'   p-values), `"all"`, a vector of feature IDs or one-based indices, or a
#'   logical vector with one value per feature.
#' @param adjust,q.value Adjustment and threshold used when
#'   `features = "discoveries"`. The defaults are those of the original call.
#' @param BPPARAM A `BiocParallelParam`; defaults to the registered `bpparam()`.
#' @param chunk_size Positive number of features per task; see
#'   [mgcvST.estimate()].
#' @param checkpoint_dir,resume Resumable chunk checkpoints as in
#'   [mgcvST.estimate()]. The directory of the original call can be reused:
#'   its step 2 chunks are keyed by their features and responses. Pass the
#'   same `chunk_size` as for the original call, or an explicit one: without it
#'   the chunks depend on the number of workers.
#' @return The updated `mgcvST_model_fit`.
#' @keywords internal
mgcvST.estimate_spatial <- function(
    fitmgcvST, Y, features = "discoveries",
    adjust = c("BY", "BH", "Sidak", "none"), q.value = 0.05,
    BPPARAM = BiocParallel::bpparam(), chunk_size = NULL,
    checkpoint_dir = NULL, resume = TRUE) {
  fit <- fitmgcvST
  if (!inherits(fit, "mgcvST_model_fit") || identical(fit$estimator, "INLA")) {
    stop("fitmgcvST must be returned by mgcvST.estimate(); use ",
         "inlaST.estimate_spatial() for an inlaST.estimate() fit.")
  }
  .mgcvst_check_fit_format(fit)
  if (is.null(fit$estimation_context)) {
    stop("The fit was estimated before the two-step estimator and cannot be ",
         "extended; re-run mgcvST.estimate().")
  }
  adjust <- if (missing(adjust)) fit$stage1$adjust else match.arg(adjust)
  q.value <- if (missing(q.value)) fit$stage1$q.value else
    .mgcvst_check_q_value(q.value)
  if (!inherits(BPPARAM, "BiocParallelParam")) {
    stop("BPPARAM must inherit from 'BiocParallelParam'.")
  }
  chunk_size <- .mgcvst_check_chunk_size(chunk_size)
  if (!is.matrix(Y)) Y <- as.matrix(Y)
  marginal_q <- .mgcvst_stage1_q(fit$diagnostics$marginal_p_value, adjust)
  index <- .mgcvst_select_spatial(features, fit$feature_id, marginal_q, q.value)
  index <- setdiff(index, which(fit$diagnostics$spatial_fitted))
  .mgcvst_check_responses(fit, Y, index)
  store <- .mgcvst_chunk_store(checkpoint_dir, "mgcv", fit$signature, resume)
  t0 <- proc.time()[["elapsed"]]
  fit$diagnostics$spatial_selected[index] <- TRUE
  fit <- .mgcvst_apply_spatial(fit, Y, index, BPPARAM, chunk_size, store)
  fit$timing$elapsed <- fit$timing$elapsed + proc.time()[["elapsed"]] - t0
  fit
}

# A null of mgcvST.estimate() is Stage 1 of spde.scan() when it has no smooth,
# the family is nb() with an estimated size and log link, and the tested smooth
# is an spde smooth; options(mgcvST.scan_null = FALSE) keeps the GAM nulls.
.mgcvst_scan_eligible <- function(G, target_index, null_formula, family,
                                  gam_args) {
  if (!isTRUE(getOption("mgcvST.scan_null", TRUE))) return(FALSE)
  if (length(gam_args)) return(FALSE)
  if (length(target_index) != 1L ||
      !inherits(G$smooth[[target_index]], "spde.smooth")) return(FALSE)
  smooth <- tryCatch(length(mgcv::interpret.gam(null_formula)$smooth.spec),
                     error = function(e) NA_integer_)
  if (!identical(smooth, 0L)) return(FALSE)
  id <- tryCatch(.working_family_id(family$family), error = function(e) NA)
  isTRUE(identical(id, "negative_binomial")) && isTRUE(family$n.theta == 1L) &&
    identical(family$link, "log")
}

# Stage 1 rows of one payload by spde.scan() on the raw field of the spde
# smooth: the fields that .mgcvst_null_chunk reports per feature.
.mgcvst_scan_chunk <- function(payload, G0, target_index, X0, max_eps = 1e-8,
                               max_iter = 1e5) {
  field <- .spde_scan_field(G0$smooth[[target_index]]$basis)
  base <- if (is.null(G0$offset)) numeric(ncol(payload$Y)) else as.numeric(G0$offset)
  offset <- if (is.null(payload$offset)) base else if (is.matrix(payload$offset))
    sweep(payload$offset, 2L, base, "+") else base + payload$offset
  scan <- spde.scan(payload$Y, X0, field$A, field$Q, offset = offset,
                    poisson = payload$poisson, feature_id = payload$feature_id,
                    max_eps = max_eps, max_iter = max_iter)
  lapply(seq_len(nrow(scan)), function(j) {
    failed <- is.na(scan$p_value[j])
    list(marginal_p_value = scan$p_value[j],
         marginal_requested_method = if (failed) NA_character_ else "davies",
         marginal_method = scan$method[j], marginal_fallback = scan$fallback[j],
         marginal_error = if (failed) list(class = "simpleError",
           message = paste("spde.scan:", scan$error[j]), call = "") else NULL)
  })
}
