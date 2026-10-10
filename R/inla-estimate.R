# Two-step sparse INLA estimation.
#
# Step 1 fits the null model of every feature and computes its Stage 1
# null-first p-value inside the worker; step 2 fits the spatial model of the
# selected features only and computes their score vector a_j and the mean
# mu_bar of the fitted mean inside the worker as well. A worker therefore
# returns compact per-feature results, and the manager never holds an
# observation-length working vector of any feature.

# Cluster workers start with their own default library stack: a package
# installed into a non-default library (a scratch or site library) is invisible
# to them, and -- worse -- a DIFFERENT build of the same package sitting in the
# default library would be picked up silently. Sending the manager's .libPaths()
# with every task fixes both. The wrapper lives in baseenv() so that a worker
# which cannot yet load mgcvST is still able to deserialize it; only after the
# library stack is set does it reach into the namespace for the real worker
# function, named by `fun_name`. Every serialized argument (sparse Matrix
# blocks, plain lists and vectors) is likewise free of mgcvST classes.
.inlast_chunk_task <- function() {
  task <- function(payload, fun_name, args, libpaths) {
    if (length(libpaths)) .libPaths(unique(c(libpaths, .libPaths())))
    # A worker can arrive with mgcvST already loaded from another library:
    # deserialising exported globals (e.g. testthat's topLevelEnvironment
    # option under SnowParam's exportglobals) loads the namespace by name
    # BEFORE this body runs, resolving against the worker's default paths.
    # A stale build then shadows the manager's; reload from the right path.
    expected <- find.package("mgcvST", lib.loc = .libPaths(), quiet = TRUE)
    if ("mgcvST" %in% loadedNamespaces() && length(expected)) {
      current <- getNamespaceInfo(asNamespace("mgcvST"), "path")
      same <- identical(normalizePath(current, winslash = "/", mustWork = FALSE),
                        normalizePath(expected[[1L]], winslash = "/", mustWork = FALSE))
      if (!same) try(unloadNamespace("mgcvST"), silent = TRUE)
    }
    fun <- get(fun_name, envir = asNamespace("mgcvST"))
    do.call(fun, c(list(payload), args))
  }
  environment(task) <- baseenv()
  task
}

# The part of a fit that is not an observation-length vector.
.inlast_compact <- function(z) {
  z[c("working_error", "working_variance", "eta", "mu", "inla")] <- NULL
  z
}

# Step 1 worker: null fits and Stage 1 statistics of one chunk of features, in
# sub-blocks of `block` features so that a worker holds at most `block`
# features' working vectors. Returns, per feature, the compact null state and
# the marginal row, or the error record.
.inlast_null_chunk <- function(payload, spec, null_spec, base_offset,
                               null_control, threads = 1L, block = 16L) {
  geometry <- .inlast_score_geometry_from_spec(spec)
  nuisance_design <- as.matrix(spec$nuisance_design)
  n_sp <- spec$geometry_sp_length
  k <- length(payload$index)
  out <- vector("list", k)
  for (first in seq.int(1L, k, by = block)) {
    rows <- first:min(k, first + block - 1L)
    fits <- .inlast_fit_chunk(
      rows, payload$Y, null_spec, base_offset, payload$extra_offset,
      null_control, FALSE, payload$poisson
    )
    failed <- vapply(fits, inherits, logical(1L), what = "condition")
    ok <- !failed & vapply(fits, function(z) isTRUE(z$converged), logical(1L))
    dispersion <- rep(NA_real_, length(rows))
    smoothing <- matrix(NA_real_, length(rows), n_sp)
    for (kk in which(!failed)) {
      dispersion[kk] <- fits[[kk]]$dispersion
      sp <- fits[[kk]]$smoothing_parameters
      if (length(sp)) smoothing[kk, seq_along(sp)] <- sp
    }
    marginal <- NULL
    row_of <- integer(length(rows))
    if (any(ok)) {
      marginal <- .inlast_null_marginal(
        payload$feature_id[rows], geometry, nuisance_design, fits, null_spec,
        dispersion, smoothing, which(ok), chunk_size = length(rows),
        threads = threads
      )
      row_of[which(ok)] <- seq_len(sum(ok))
    }
    for (kk in seq_along(rows)) {
      z <- fits[[kk]]
      out[[rows[kk]]] <- if (failed[kk]) {
        list(error = .mgcvst_condition(z))
      } else {
        list(null = .inlast_compact(z),
             marginal = if (row_of[kk]) as.list(marginal[row_of[kk], ]) else NULL)
      }
    }
    rm(fits)
  }
  .mgcvst_chunk_save(payload$chunk_file, payload$chunk_key, out)
  out
}

# What the manager keeps of one spatial fit: the compact estimates, never an
# observation-length vector.
.inlast_spatial_record <- function(z, spec, scores_a, score_error, mu_bar,
                                   retain_smooth, diagnostics, estimation) {
  record <- list(
    dispersion = z$dispersion, family_parameters = z$family_parameters,
    smoothing_parameters = z$smoothing_parameters,
    target = as.numeric(z$random_mode[[1L]]),
    nuisance = if (ncol(spec$nuisance_design)) .inlast_nuisance_mode(z, spec) else
      numeric(),
    converged = isTRUE(z$converged),
    log_marginal_likelihood = z$log_marginal_likelihood,
    fit_seconds = z$fit_seconds, spatial_fallback = z$spatial_fallback,
    constraint_residual = z$constraint_residual,
    observation_spatial_mean = z$observation_spatial_mean,
    score_a = scores_a, score_error = score_error, mu_bar = mu_bar
  )
  if (estimation) record$estimation <- z$estimation
  if (retain_smooth) record$coefficients <- z$coefficients
  if (diagnostics) {
    record$diagnostics <- z[c(
      "tau", "tau_internal", "precision_scale", "lambda", "mode_status",
      "mode_status_text", "constraint_residual", "constraint_residual_uncorrected",
      "observation_spatial_mean", "estimation", "spatial_fallback"
    )]
  }
  record
}

# Step 2 worker: spatial fits of one chunk of features, their score vectors
# a_j and the mean mu_bar of the fitted mean. A native INLA crash falls back
# on the compact null state of the feature (payload$null).
.inlast_spatial_chunk <- function(payload, spec, base_offset, control,
                                  diagnostics, retain_smooth, threads = 1L,
                                  block = 16L) {
  geometry <- .inlast_score_geometry_from_spec(spec)
  nuisance_design <- as.matrix(spec$nuisance_design)
  n_sp <- spec$geometry_sp_length
  k <- length(payload$index)
  out <- vector("list", k)
  have_estimation <- FALSE
  for (first in seq.int(1L, k, by = block)) {
    rows <- first:min(k, first + block - 1L)
    fits <- .inlast_fit_chunk(
      rows, payload$Y, spec, base_offset, payload$extra_offset, control,
      diagnostics, payload$poisson
    )
    for (kk in seq_along(rows)) {
      if (!inherits(fits[[kk]], "inlaCrashError")) next
      j <- rows[kk]
      feature_offset <- base_offset
      if (!is.null(payload$extra_offset)) feature_offset <- feature_offset +
        if (is.matrix(payload$extra_offset)) payload$extra_offset[j, ] else
          payload$extra_offset
      fits[[kk]] <- .inlast_spatial_fallback(
        fits[[kk]], payload$null[[j]], spec, payload$Y[j, ], feature_offset
      )
    }
    failed <- vapply(fits, inherits, logical(1L), what = "condition")
    converged <- !failed & vapply(fits, function(z) isTRUE(z$converged), logical(1L))
    dispersion <- rep(NA_real_, length(rows))
    smoothing <- matrix(NA_real_, length(rows), n_sp)
    for (kk in which(!failed)) {
      dispersion[kk] <- fits[[kk]]$dispersion
      smoothing[kk, ] <- fits[[kk]]$smoothing_parameters
    }
    scores <- .inlast_estimate_scores(
      fits, geometry, spec, nuisance_design, dispersion, smoothing,
      which(converged), threads = threads
    )
    for (kk in seq_along(rows)) {
      z <- fits[[kk]]
      out[[rows[kk]]] <- if (failed[kk]) {
        list(error = .mgcvst_condition(z))
      } else {
        record <- .inlast_spatial_record(
          z, spec, scores$a[, kk], scores$error[kk], mean(z$mu),
          retain_smooth, diagnostics, estimation = !have_estimation
        )
        have_estimation <- TRUE
        record
      }
    }
    rm(fits, scores)
  }
  .mgcvst_chunk_save(payload$chunk_file, payload$chunk_key, out)
  out
}

# Checkpoint signature of an INLA estimation: the model, offset and controls
# the chunk results depend on. The chunking and the family routing of a feature
# are not part of it: a chunk file is keyed by its step, its features, their
# responses and their routing, so a changed chunk_size or a feature that is
# routed differently recomputes its chunk, and never mixes.
.inlast_estimation_signature <- function(model, control, offset, feature_id,
                                         retain_smooth, diagnostics) {
  .mgcvst_pair_input_hash(list(
    format = .mgcvst_fit_format, spec = model$inla_spec, offset = model$offset,
    extra_offset = offset, control = control, feature_id = feature_id,
    retain_smooth = retain_smooth, diagnostics = diagnostics
  ))
}

.inlast_check_threads <- function(threads) {
  if (length(threads) != 1L || !is.numeric(threads) || !is.finite(threads) ||
      threads < 1 || threads > .Machine$integer.max ||
      threads != as.integer(threads)) stop("threads must be a positive integer.")
  as.integer(threads)
}

#' Estimate mgcvST working models with sparse INLA
#'
#' Fits one latent Gaussian model per feature with INLA and returns the compact
#' working-model contract consumed by [inlaST.test()]. No GAM fit is used.
#'
#' Estimation has two steps. Step 1 fits only the null model (no spatial
#' field) of every feature and computes the Stage 1 null-first p-value of its
#' marginal score test inside the worker; the p-values are then adjusted by
#' `adjust` into the Stage 1 q-values `diagnostics$marginal_q_value`. Step 2
#' fits the spatial model of the features chosen by `spatial` only, and
#' computes their score vector and the mean of the fitted mean in the worker.
#' No worker returns an observation-length vector: the manager holds the
#' compact per-feature estimates only, so its memory grows with the number of
#' features and the mesh size, not with the number of observations. Features
#' without a spatial model are marked `spatial_fitted = FALSE` in the
#' diagnostics table and are unavailable to [inlaST.test()] and
#' [inlaST.wgcna()]; `pairs = NULL` in [inlaST.test()] means all pairs among
#' the features with a spatial model. Use [inlaST.estimate_spatial()] to add
#' spatial models for more features after step 1.
#'
#' With `checkpoint_dir`, every completed chunk of either step is saved by the
#' worker that computed it, and a repeated call with the same arguments
#' resumes from the saved chunks. Chunks are keyed by their step, features,
#' responses and family routing, so pass `chunk_size` explicitly for a
#' resumable run: without it the chunks depend on the number of workers (and
#' step 2 on the number of selected features), and a changed chunking
#' recomputes rather than reuses. A chunk whose result holds a failed feature
#' is computed again on resume. A checkpoint directory written for another
#' model, offset, control or version is refused. Without `chunk_size`, a run
#' with a checkpoint directory uses chunks of at most 50 features.
#'
#' @param Y Numeric feature-by-observation matrix.
#' @param model An object returned by [inlaST.set()].
#' @param feature_id Unique feature identifiers.
#' @param BPPARAM A `BiocParallelParam` distributing feature chunks over
#'   workers. Each worker fits its own features with INLA using
#'   `control$num_threads` (default `1L`). Setting it to `NULL` lets INLA choose
#'   its own thread count, which can oversubscribe BiocParallel workers.
#'   `SerialParam()` keeps everything in one process.
#' @param chunk_size Positive integer number of features per task, or `NULL`.
#'   The default creates one chunk per worker in each step; a smaller value
#'   makes the checkpoint finer. A chunk payload is built when a worker takes
#'   the chunk, so at most about one chunk of `Y` per worker is copied.
#' @param offset Optional shared observation offset or matrix matching `Y`. A
#'   feature-by-observation matrix is kept on the fit (`extra_offset`) for a
#'   later step 2, in addition to the total offset.
#' @param control Named INLA engine overrides for controls saved by
#'   [inlaST.set()]. Omitted entries inherit model settings. Each supplied
#'   prior list replaces the whole prior; numerical `control.inla` entries
#'   merge by name. Explicit `NULL` resets optional fixed parameters, except
#'   NB size fixed by the model family. The supported approximation is
#'   `int_strategy = "eb"`, `latent_strategy = "gaussian"`. `num_threads`
#'   defaults to `1L`; set it explicitly to raise INLA's thread count per worker.
#'   Fixed effects always use explicit zero precision in INLA.
#'   Optional positive `fixed_precision`,
#'   `gaussian_precision`, and `nb_size` fix latent precision multipliers,
#'   inverse Gaussian residual variance, and NB size, respectively.
#'   `fixed_precision` values always refer to the original FEM multiplier,
#'   including when the model uses observation-scale precision priors. Prior
#'   lists `precision_prior`, `gaussian_precision_prior` and `nb_size_prior`
#'   contain `prior`, `param`, and logarithmic `initial` values (default zero).
#'   Spatial precision and NB size default to `prior = "flat"` with no
#'   parameters on the log scale. Gaussian observation precision defaults to
#'   `prior = "flat", param = numeric(), initial = 0`, a flat prior on the
#'   internal log precision `-log(variance)`. On the variance scale this has
#'   density proportional to `1/variance`; it is a Gaussian dispersion prior,
#'   distinct from the negative-binomial size prior. Custom normal, registered
#'   scalar INLA priors, and INLA expression/table priors remain supported.
#'   Normal parameters are mean and precision. A flat log-hyperparameter
#'   objective corresponds to density proportional to `1/parameter` on its
#'   positive scale.
#'   `control.inla` accepts supported numerical tuning, with Gaussian latent
#'   strategy and EB integration enforced. `poisson_screen_phi` (default `1.01`)
#'   is the Poisson prescreen threshold: with a negative-binomial family, each
#'   feature first gets an offset-and-covariate-only Poisson GLM, and a feature
#'   whose Pearson dispersion `phi = sum((y - mu)^2 / mu) / (n - p)` is at most
#'   the threshold is fitted with the Poisson family instead. Set it to `0` to
#'   disable the screen (`NULL` restores the default); other families ignore it.
#'   The per-feature
#'   `phi` and the family actually used are reported in the diagnostics as
#'   `prescreen_phi` and `family_used`. Unknown controls are rejected.
#' @param retain_smooth Retain estimated score-component coefficients.
#' @param diagnostics Retain the per-feature INLA optimizer/hyperparameter
#'   diagnostics in `inla_diagnostics`.
#' @param threads OpenMP threads per worker for the sparse marginal and score
#'   computations of the workers.
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
#'   of both steps; see Details.
#' @param resume Reuse the completed chunks of a compatible checkpoint.
#' @details Hyperparameter priors remain part of INLA's empirical-Bayes
#' estimates; these are not mgcv REML estimates. `mgcv::nb(theta = value)`
#' fixes NB size, and a conflicting `control$nb_size` is rejected. The
#' working model uses conditional latent estimates and expected Fisher
#' variances with INLA `config = FALSE` and variational-Bayes correction
#' disabled. The small coefficient covariance
#' is reconstructed by sparse precision solves without an observation-level inverse.
#' Cross-feature iid effects are treated as independent in pairwise calibration.
#' Sparse downstream scores rebuild the nuisance covariance from
#' the same expected Fisher matrix;
#' this does not establish finite-sample calibration after hyperparameter
#' estimation. The marginal score test of each null fit
#' (`diagnostics$marginal_p_value`) is calibrated by Davies on the positive
#' eigenvalues of the full-space expected curvature. When Davies errors,
#' returns a missing or non-finite p-value, or returns a value outside
#' (0, 1], the saddlepoint approximation (Kuonen 1999) is used instead;
#' `marginal_method` reports `"davies"` or `"saddlepoint"` and
#' `marginal_fallback` flags the fallback. Pair tests in [inlaST.test()] are
#' calibrated by Liu moment matching. Every spatial component's constraint
#' residual and observed spatial mean are retained in the result.
#' The score uses the SPDE covariance conditioned on observation mean zero.
#' The sparse score uses a matching expected-curvature nuisance adjustment.
#' The INLA path is sparse-only: the sparse kernel is the sole score and
#' marginal implementation, with no dense score or custom marginal callback.
#' If INLA's native program crashes during the spatial fit (`inlaCrashError`)
#' after a successful null fit, that null fit supplies the fixed and iid
#' nuisance effects, family parameters and working state. The spatial effect
#' is set to zero and its original-FEM precision is assigned `1e8`;
#' `lambda = dispersion * 1e8`. No additional fit is attempted. The null score
#' p-value is retained and pair scores use the usual working-model formulas.
#' `diagnostics$spatial_fallback` identifies these usable fallback states;
#' `outer_convergence` is `"null_zero_spatial_fallback"` and the original
#' error remains in `error_class`, `error_message` and `error_call`. Their
#' spatial likelihood criterion is missing. Input errors and unsuccessful
#' null fits do not use this fallback. A design the sparse capability gate rejects makes
#' [inlaST.set()] error at setup, and a model that somehow reaches this
#' function without that geometry errors here.
#' Flat hyperpriors need not yield proper hyperparameter posteriors. They are
#' supported only as empirical-Bayes optimization objectives in the current
#' single-configuration engine. A returned finite precision or zero optimizer
#' status does not establish that the maximum is interior; zero spatial
#' variance and the Poisson limit of the NB model require boundary checks.
#' @return An `inlaST_fit` that is also an `mgcvST_model_fit`. In addition to
#'   the estimates, it holds `mu_bar` (the mean of each fitted mean),
#'   `null_state` (the compact null estimates used by a later step 2),
#'   `y_digest` (a digest of each response row, which a later step 2 checks)
#'   and the diagnostics columns `marginal_q_value`, `spatial_selected` and
#'   `spatial_fitted`.
#' @seealso [inlaST.estimate_spatial()] to add spatial models after step 1.
#' @export
inlaST.estimate <- function(
    Y, model, feature_id = rownames(Y),
    BPPARAM = BiocParallel::SerialParam(), chunk_size = NULL,
    offset = NULL, control = list(), retain_smooth = FALSE,
    diagnostics = FALSE, threads = 1L, spatial = "discoveries",
    adjust = c("BY", "BH", "Sidak", "none"), q.value = 0.05,
    checkpoint_dir = NULL, resume = TRUE) {
  if (!requireNamespace("INLA", quietly = TRUE)) {
    stop("inlaST.estimate() requires the INLA package.")
  }
  adjust <- match.arg(adjust)
  q.value <- .mgcvst_check_q_value(q.value)
  threads <- .inlast_check_threads(threads)
  checked <- .inlast_validate_estimate(
    Y, model, feature_id, BPPARAM, chunk_size, offset, control,
    retain_smooth, diagnostics
  )
  # Sparse-only: the geometry builder is the capability gate. It stops with the
  # reason when the model cannot be scored; there is no dense alternative.
  score_sparse <- .inlast_sparse_score_geometry(model)
  # Validate once in the parent so a misspelled or invalid control cannot turn
  # every feature into an otherwise opaque per-feature failure.
  control <- .inlast_control(.inlast_merge_control(model$inla_control, control))
  control <- .inlast_family_control(model, control)
  Y <- checked$Y
  feature_id <- checked$feature_id
  chunk_size <- checked$chunk_size
  p <- nrow(Y)
  # Validate the spatial selection before any fit is paid for.
  if (!(is.character(spatial) && length(spatial) == 1L &&
        spatial %in% c("discoveries", "all", "none"))) {
    .mgcvst_select_spatial(spatial, feature_id, rep(NA_real_, p), q.value)
  }
  # Cheap per-feature overdispersion prescreen (R/family-prescreen.R): a
  # near-Poisson feature drifts to a flat NB size likelihood, so fit it with the
  # Poisson family instead. The family is fixed inside the spec, so a routed
  # feature is fitted against a Poisson twin of the same spec
  # (.inlast_poisson_spec()); the twin is rebuilt once per chunk in the worker
  # rather than shipped, because it shares every sparse block with the original.
  screen_threshold <- .mgcvst_prescreen_threshold(
    control[["poisson_screen_phi", exact = TRUE]]
  )
  prescreen <- .mgcvst_prescreen_route(
    Y, X = model$inla_spec$fixed$X,
    offset = if (is.null(offset)) model$offset else
      if (is.matrix(offset)) sweep(offset, 2L, model$offset, "+") else
        model$offset + offset,
    threshold = screen_threshold,
    active = identical(model$inla_spec$family, "negative_binomial")
  )
  # Validate the twin once in the parent so a spec the switch cannot handle
  # fails here rather than inside every affected worker.
  if (any(prescreen$poisson)) .inlast_poisson_spec(model$inla_spec)
  family_used <- rep(model$inla_spec$family, p)
  family_used[prescreen$poisson] <- "poisson"
  null_spec <- .inlast_null_spec(model$inla_spec)
  null_control <- .inlast_null_control(control, model$inla_spec)
  if (any(prescreen$poisson)) .inlast_poisson_spec(null_spec)

  t0 <- proc.time()[["elapsed"]]
  y_digest <- .mgcvst_row_digests(Y)
  signature <- .inlast_estimation_signature(
    model, control, offset, feature_id, retain_smooth, diagnostics
  )
  store <- .mgcvst_chunk_store(checkpoint_dir, "inla", signature, resume)
  groups <- .mgcvst_feature_chunks(seq_len(p), chunk_size, BPPARAM,
                                   checkpoint = !is.null(store))
  # A payload is built only when its chunk is computed and a worker is free.
  make_payload <- function(index) {
    extra_offset <- if (is.null(offset) || !is.matrix(offset)) offset else
      offset[index, , drop = FALSE]
    list(index = index, feature_id = feature_id[index],
         Y = .mgcvst_double_rows(Y, index), extra_offset = extra_offset,
         poisson = if (any(prescreen$poisson)) prescreen$poisson[index] else NULL)
  }
  null_t0 <- proc.time()[["elapsed"]]
  step1 <- .mgcvst_run_chunks(
    groups, make_payload, "null", store, y_digest, prescreen$poisson, BPPARAM,
    .inlast_chunk_task(),
    fun_name = ".inlast_null_chunk",
    args = list(spec = model$inla_spec, null_spec = null_spec,
                base_offset = model$offset, null_control = null_control,
                threads = threads),
    libpaths = .libPaths()
  )
  null_elapsed <- proc.time()[["elapsed"]] - null_t0
  nulls <- unlist(step1$results, recursive = FALSE)

  n_sp <- model$inla_spec$geometry_sp_length
  m <- ncol(score_sparse$A)
  px <- ncol(model$geometry$nuisance_design)
  null_converged <- rep(FALSE, p)
  null_fit_seconds <- rep(NA_real_, p)
  null_error_class <- null_error_message <- null_error_call <- rep(NA_character_, p)
  null_state <- stats::setNames(vector("list", p), feature_id)
  marginal <- data.frame(
    p_value = rep(NA_real_, p), method_requested = NA_character_,
    method_used = NA_character_, fallback_used = NA, error_message = NA_character_,
    stringsAsFactors = FALSE
  )
  for (j in seq_len(p)) {
    z <- nulls[[j]]
    if (!is.null(z$error)) {
      null_error_class[j] <- sub("/.*$", "", z$error$class)
      null_error_message[j] <- z$error$message
      null_error_call[j] <- z$error$call
      next
    }
    null_state[j] <- list(z$null)
    null_converged[j] <- isTRUE(z$null$converged)
    null_fit_seconds[j] <- z$null$fit_seconds
    if (!is.null(z$marginal)) {
      marginal$p_value[j] <- z$marginal$p_value
      marginal$method_requested[j] <- z$marginal$method_requested
      marginal$method_used[j] <- z$marginal$method_used
      marginal$fallback_used[j] <- z$marginal$fallback_used
      marginal$error_message[j] <- z$marginal$error_message
    }
  }
  rm(nulls)
  marginal_q <- .mgcvst_stage1_q(marginal$p_value, adjust)

  diagnostics_table <- data.frame(
    index = seq_len(p), feature_id = feature_id, converged = FALSE,
    marginal_p_value = marginal$p_value, marginal_q_value = marginal_q,
    marginal_requested_method = marginal$method_requested,
    marginal_method = marginal$method_used,
    marginal_fallback = marginal$fallback_used,
    residual_df = NA_real_, criterion = NA_real_,
    criterion_name = "INLA log marginal likelihood", fit_seconds = NA_real_,
    outer_convergence = NA_character_, error_class = NA_character_,
    error_message = NA_character_, error_call = NA_character_,
    spatial_selected = FALSE, spatial_fitted = FALSE,
    spatial_fallback = FALSE, spatial_fallback_method = NA_character_,
    spatial_precision_assigned = NA_real_,
    null_converged = null_converged, null_fit_seconds = null_fit_seconds,
    null_error_class = null_error_class, null_error_message = null_error_message,
    null_error_call = null_error_call,
    score_error_class = NA_character_, score_error_message = NA_character_,
    score_error_call = NA_character_,
    prescreen_phi = prescreen$phi, family_used = family_used,
    stringsAsFactors = FALSE
  )
  failed_score <- !is.na(marginal$error_message) & nzchar(marginal$error_message)
  if (any(failed_score)) {
    diagnostics_table$score_error_class[failed_score] <- "sparse_score"
    diagnostics_table$score_error_message[failed_score] <- marginal$error_message[failed_score]
  }
  coefficient <- if (retain_smooth) {
    lapply(model$geometry$target, function(j) {
      matrix(NA_real_, p, .inlast_target_width(model, j),
             dimnames = list(feature_id, NULL))
    })
  } else NULL
  if (!is.null(coefficient)) names(coefficient) <- names(model$geometry$target)
  component_lambda <- matrix(NA_real_, p, length(model$geometry$target),
    dimnames = list(feature_id, names(model$geometry$target)))
  constraint_residual <- observation_spatial_mean <- matrix(
    NA_real_, p, length(model$inla_spec$random),
    dimnames = list(feature_id, vapply(model$inla_spec$random, `[[`, character(1L), "name"))
  )
  total_offset <- if (is.null(offset)) model$offset else if (is.matrix(offset))
    sweep(offset, 2L, model$offset, "+") else model$offset + offset
  ans <- structure(list(
    feature_id = feature_id,
    dispersion = stats::setNames(rep(NA_real_, p), feature_id),
    lambda = component_lambda[, "global"],
    component_lambda = component_lambda,
    smoothing_parameters = matrix(NA_real_, p, n_sp,
      dimnames = list(feature_id, model$inla_spec$sp_names)),
    family_parameters = stats::setNames(vector("list", p), feature_id),
    feature_family = family_used,
    target_coefficients = matrix(NA_real_, m, p, dimnames = list(NULL, feature_id)),
    nuisance_coefficients = matrix(NA_real_, px, p, dimnames = list(NULL, feature_id)),
    score_a = matrix(NA_real_, m, p, dimnames = list(NULL, feature_id)),
    mu_bar = stats::setNames(rep(NA_real_, p), feature_id),
    null_state = null_state,
    y_digest = y_digest,
    n_observation = ncol(Y),
    basis_spec = list(kind = "full_rank", rank = m - 1L),
    format = .mgcvst_fit_format,
    signature = signature,
    control = control,
    extra_offset = offset,
    routed_poisson = prescreen$poisson,
    geometry = model$geometry,
    row_id = model$geometry$row_id,
    offset = total_offset,
    linear_design = model$geometry$X,
    score_components = model$geometry$score_components,
    model_setting = model$setting,
    model = model,
    diagnostics = diagnostics_table,
    timing = list(elapsed = NA_real_, null_fit_elapsed = null_elapsed,
                  fit_elapsed = 0, chunks = length(groups),
                  chunk_size = chunk_size, backend = class(BPPARAM)[1L],
                  workers = max(1L, min(length(groups),
                                        BiocParallel::bpworkers(BPPARAM))),
                  resumed_null_chunks = step1$resumed,
                  null_payloads = step1$built,
                  spatial_chunks = 0L, resumed_spatial_chunks = 0L,
                  spatial_payloads = 0L),
    smooth_coefficients = coefficient,
    retain_smooth = retain_smooth,
    test_engine = "single_model",
    score_backend = "sparse",
    score_sparse = score_sparse,
    estimator = "INLA",
    estimation = NULL,
    constraint_residual = constraint_residual,
    observation_spatial_mean = observation_spatial_mean,
    inla_diagnostics = if (diagnostics) stats::setNames(vector("list", p), feature_id) else NULL,
    mean_constraint = "observation",
    mean_constraint_active = TRUE,
    stage1 = list(adjust = adjust, q.value = q.value),
    call = match.call()
  ), class = c("inlaST_fit", "mgcvST_model_fit", "mgcvST_fit", "mgcvST"))
  kappa <- .spde_kappa_fields(model)
  ans[names(kappa)] <- kappa
  rm(step1)

  index <- .mgcvst_select_spatial(spatial, feature_id, marginal_q, q.value)
  ans$diagnostics$spatial_selected[index] <- TRUE
  ans <- .inlast_apply_spatial(ans, Y, index, BPPARAM, chunk_size, threads, store)
  ans$timing$elapsed <- proc.time()[["elapsed"]] - t0
  ans
}

#' Add spatial models after step 1 of INLA estimation
#'
#' Fits the spatial model of further features of a fit returned by
#' [inlaST.estimate()] and returns the updated fit; the supplied fit is not
#' changed. Features that already have a spatial model are skipped. The
#' responses `Y` must be those of step 1, which is checked against the digests
#' stored in the fit. The controls, offset, family routing and retention
#' options of the original call are reused.
#'
#' @param fitinlaST A fit returned by [inlaST.estimate()] (version 0.0.1.9032 or
#'   later).
#' @param Y The feature-by-observation matrix given to [inlaST.estimate()].
#' @param features The features to add: `"discoveries"` (the default; Stage 1
#'   q-value at most `q.value` under `adjust`, computed from the stored Stage 1
#'   p-values), `"all"`, a vector of feature IDs or one-based indices, or a
#'   logical vector with one value per feature.
#' @param adjust,q.value Adjustment and threshold used when
#'   `features = "discoveries"`. The defaults are those of the original call.
#' @inheritParams inlaST.estimate
#' @param checkpoint_dir,resume Resumable chunk checkpoints as in
#'   [inlaST.estimate()]. The directory of the original call can be reused:
#'   its step 2 chunks are keyed by their features and responses. Pass the
#'   same `chunk_size` as for the original call, or an explicit one: without it
#'   the chunks depend on the number of workers.
#' @return The updated `inlaST_fit`.
#' @export
inlaST.estimate_spatial <- function(
    fitinlaST, Y, features = "discoveries",
    adjust = c("BY", "BH", "Sidak", "none"), q.value = 0.05,
    BPPARAM = BiocParallel::SerialParam(), chunk_size = NULL,
    checkpoint_dir = NULL, resume = TRUE, threads = 1L) {
  if (!requireNamespace("INLA", quietly = TRUE)) {
    stop("inlaST.estimate_spatial() requires the INLA package.")
  }
  fit <- fitinlaST
  if (!inherits(fit, "inlaST_fit")) {
    stop("fitinlaST must be returned by inlaST.estimate().")
  }
  .mgcvst_check_fit_format(fit)
  adjust <- if (missing(adjust)) fit$stage1$adjust else match.arg(adjust)
  q.value <- if (missing(q.value)) fit$stage1$q.value else
    .mgcvst_check_q_value(q.value)
  threads <- .inlast_check_threads(threads)
  if (!inherits(BPPARAM, "BiocParallelParam")) {
    stop("BPPARAM must inherit from 'BiocParallelParam'.")
  }
  chunk_size <- .mgcvst_check_chunk_size(chunk_size)
  if (!is.matrix(Y)) Y <- as.matrix(Y)
  marginal_q <- .mgcvst_stage1_q(fit$diagnostics$marginal_p_value, adjust)
  index <- .mgcvst_select_spatial(features, fit$feature_id, marginal_q, q.value)
  index <- setdiff(index, which(fit$diagnostics$spatial_fitted))
  .mgcvst_check_responses(fit, Y, index)
  store <- .mgcvst_chunk_store(checkpoint_dir, "inla", fit$signature, resume)
  t0 <- proc.time()[["elapsed"]]
  fit$diagnostics$spatial_selected[index] <- TRUE
  fit <- .inlast_apply_spatial(fit, Y, index, BPPARAM, chunk_size, threads, store)
  fit$timing$elapsed <- fit$timing$elapsed + proc.time()[["elapsed"]] - t0
  fit
}

# Step 2: fit the spatial model of the features `index` that have none yet and
# write the compact results into `fit`.
.inlast_apply_spatial <- function(fit, Y, index, BPPARAM, chunk_size, threads,
                                  store) {
  model <- fit$model
  spec <- model$inla_spec
  index <- setdiff(index, which(fit$diagnostics$spatial_fitted))
  if (!length(index)) return(fit)
  offset <- fit$extra_offset
  poisson <- fit$routed_poisson
  groups <- .mgcvst_feature_chunks(index, chunk_size, BPPARAM,
                                   checkpoint = !is.null(store))
  make_payload <- function(i) {
    extra_offset <- if (is.null(offset) || !is.matrix(offset)) offset else
      offset[i, , drop = FALSE]
    list(index = i, feature_id = fit$feature_id[i], Y = .mgcvst_double_rows(Y, i),
         extra_offset = extra_offset, null = fit$null_state[i],
         poisson = if (any(poisson)) poisson[i] else NULL)
  }
  fit_t0 <- proc.time()[["elapsed"]]
  step2 <- .mgcvst_run_chunks(
    groups, make_payload, "spatial", store, fit$y_digest, poisson, BPPARAM,
    .inlast_chunk_task(),
    fun_name = ".inlast_spatial_chunk",
    args = list(spec = spec, base_offset = model$offset, control = fit$control,
                diagnostics = !is.null(fit$inla_diagnostics),
                retain_smooth = fit$retain_smooth, threads = threads),
    libpaths = .libPaths()
  )
  fit_elapsed <- proc.time()[["elapsed"]] - fit_t0
  fits <- unlist(step2$results, recursive = FALSE)
  index <- unlist(groups, use.names = FALSE)

  diag <- fit$diagnostics
  for (k in seq_along(index)) {
    j <- index[k]
    z <- fits[[k]]
    if (!is.null(z$error)) {
      diag$error_class[j] <- sub("/.*$", "", z$error$class)
      diag$error_message[j] <- z$error$message
      diag$error_call[j] <- z$error$call
      next
    }
    fit$dispersion[j] <- z$dispersion
    fit$family_parameters[[j]] <- z$family_parameters
    fit$smoothing_parameters[j, ] <- z$smoothing_parameters
    fit$target_coefficients[, j] <- z$target
    if (length(z$nuisance)) fit$nuisance_coefficients[, j] <- z$nuisance
    fit$score_a[, j] <- z$score_a
    fit$mu_bar[j] <- z$mu_bar
    diag$converged[j] <- z$converged
    diag$criterion[j] <- z$log_marginal_likelihood
    diag$fit_seconds[j] <- z$fit_seconds
    diag$outer_convergence[j] <- if (z$converged) "converged" else "failed"
    diag$spatial_fitted[j] <- TRUE
    if (!is.null(z$spatial_fallback)) {
      fallback <- z$spatial_fallback
      diag$spatial_fallback[j] <- TRUE
      diag$spatial_fallback_method[j] <- fallback$method
      diag$spatial_precision_assigned[j] <- fallback$spatial_precision
      diag$outer_convergence[j] <- "null_zero_spatial_fallback"
      diag$error_class[j] <- fallback$error_class[1L]
      diag$error_message[j] <- fallback$error_message
      diag$error_call[j] <- fallback$error_call
    }
    if (!is.na(z$score_error) && nzchar(z$score_error)) {
      diag$score_error_class[j] <- "sparse_score"
      diag$score_error_message[j] <- if (is.na(diag$score_error_message[j])) {
        z$score_error
      } else paste(diag$score_error_message[j], z$score_error, sep = " | ")
    }
    fit$constraint_residual[j, ] <- z$constraint_residual
    fit$observation_spatial_mean[j, ] <- z$observation_spatial_mean
    if (!is.null(fit$inla_diagnostics)) fit$inla_diagnostics[j] <- list(z$diagnostics)
    if (!is.null(fit$smooth_coefficients)) {
      for (name in names(fit$smooth_coefficients)) {
        fit$smooth_coefficients[[name]][j, ] <- z$coefficients[[name]]
      }
    }
    if (is.null(fit$estimation) && !is.null(z$estimation)) {
      fit$estimation <- z$estimation
      fit$estimation$control <- fit$control
    }
  }
  fit$diagnostics <- diag
  target_sp <- vapply(model$geometry$target, function(j) {
    model$geometry$smooth[[j]]$sp_index
  }, integer(1L))
  fit$component_lambda <- fit$smoothing_parameters[, target_sp, drop = FALSE]
  colnames(fit$component_lambda) <- names(target_sp)
  fit$lambda <- fit$component_lambda[, "global"]
  fit$timing$fit_elapsed <- fit$timing$fit_elapsed + fit_elapsed
  fit$timing$spatial_chunks <- fit$timing$spatial_chunks + step2$chunks
  fit$timing$resumed_spatial_chunks <- fit$timing$resumed_spatial_chunks +
    step2$resumed
  fit$timing$spatial_payloads <- fit$timing$spatial_payloads + step2$built
  fit
}
