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
  task <- function(payload, fun_name, args, libpaths, manager = NULL) {
    if (length(libpaths)) .libPaths(unique(c(libpaths, .libPaths())))
    # A worker can arrive with mgcvST already loaded from another library:
    # deserialising exported globals (e.g. testthat's topLevelEnvironment
    # option under SnowParam's exportglobals) loads the namespace by name
    # BEFORE this body runs, resolving against the worker's default paths.
    # A stale build then shadows the manager's; reload from the right path.
    expected <- find.package("spdeST", lib.loc = .libPaths(), quiet = TRUE)
    # The manager process itself (a serial backend) keeps its namespace: it is
    # the build that made this call, and under devtools::load_all() or
    # testthat::test_local() unloading it would drop the caller's own objects.
    if (!identical(manager, Sys.getpid()) &&
        "spdeST" %in% loadedNamespaces() && length(expected)) {
      current <- getNamespaceInfo(asNamespace("spdeST"), "path")
      same <- identical(normalizePath(current, winslash = "/", mustWork = FALSE),
                        normalizePath(expected[[1L]], winslash = "/", mustWork = FALSE))
      if (!same) try(unloadNamespace("spdeST"), silent = TRUE)
    }
    fun <- get(fun_name, envir = asNamespace("spdeST"))
    # A worker process keeps one BLAS/OpenMP thread for its life, so that the
    # LAPACK calls of the sparse kernels never oversubscribe the node.
    if (!identical(manager, Sys.getpid())) {
      get(".mgcvst_thread_limit", envir = asNamespace("spdeST"))()
    }
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

# Step 1 worker: Stage 1 of one chunk of features by spde.scan()
# (.inlast_scan_chunk). Returns, per feature, the compact null state and the
# marginal row, or the error record.
.inlast_null_chunk <- function(payload, spec, null_spec, base_offset,
                               null_control, threads = 1L, scan_args = list()) {
  out <- .inlast_scan_chunk(payload, spec, null_spec, base_offset,
                            null_control, threads, scan_args)
  .mgcvst_chunk_save(payload$chunk_file, payload$chunk_key, out)
  out
}

# Stage 1 is spde.scan(): the null must have no random block, a count family
# and a fixed-effects design.
.inlast_scan_check <- function(null_spec) {
  family <- null_spec$family
  if (length(null_spec$random) ||
      !(length(family) == 1L && family %in% c("negative_binomial", "poisson")) ||
      is.null(null_spec$fixed$X) || ncol(as.matrix(null_spec$fixed$X)) < 1L) {
    stop("Stage 1 (spde.scan) needs a GLM null: a negative binomial or Poisson ",
         "family, fixed effects only (linear terms or factors), no random block.")
  }
  invisible(TRUE)
}

# Stage 1 of one chunk by spde.scan() on the raw field of the target block,
# with the compact null state that step 2 reuses after a native INLA crash.
.inlast_scan_chunk <- function(payload, spec, null_spec, base_offset,
                               null_control, threads = 1L, scan_args = list()) {
  t0 <- proc.time()[["elapsed"]]
  block <- spec$random[[1L]]
  X <- as.matrix(null_spec$fixed$X)
  k <- length(payload$index)
  offset <- if (is.null(payload$extra_offset)) base_offset else
    if (is.matrix(payload$extra_offset))
      sweep(payload$extra_offset, 2L, base_offset, "+") else
      base_offset + as.numeric(payload$extra_offset)
  poisson <- if (identical(null_spec$family, "poisson")) rep(TRUE, k) else
    if (is.null(payload$poisson)) NULL else as.logical(payload$poisson)
  scan <- do.call(spde.scan, c(list(
    payload$Y, X, block$A, block$Q, offset = offset, poisson = poisson,
    size = null_control[["nb_size", exact = TRUE]],
    feature_id = payload$feature_id, threads = threads), scan_args))
  seconds <- (proc.time()[["elapsed"]] - t0) / max(1L, k)
  beta <- attr(scan, "coefficients")
  fixed_names <- colnames(X)
  if (is.null(fixed_names)) fixed_names <- paste0("x", seq_len(ncol(X)))
  n_sp <- length(null_spec[["sp", exact = TRUE]])
  empty <- stats::setNames(numeric(), character())
  lapply(seq_len(k), function(j) {
    if (is.na(scan$family[j])) {
      return(list(error = list(class = "simpleError", call = "",
                               message = paste("spde.scan null fit:", scan$error[j]))))
    }
    b <- stats::setNames(beta[, j], fixed_names)
    nb <- identical(scan$family[j], "negative_binomial")
    null <- list(
      dispersion = 1, family = scan$family[j],
      family_parameters = if (nb) scan$size[j] else numeric(),
      tau = empty, tau_internal = empty, precision_scale = empty, lambda = empty,
      smoothing_parameters = rep(NA_real_, n_sp), coefficients = list(),
      fixed_mode = b, fixed_mean = b, random_mode = list(), random_mean = list(),
      constraint_residual_uncorrected = empty, constraint_residual = empty,
      observation_spatial_mean = empty, converged = isTRUE(scan$converged[j]),
      log_marginal_likelihood = NA_real_, log_likelihood = scan$log_likelihood[j],
      fit_seconds = seconds,
      estimation = list(engine = "spde.scan",
                        overdispersion_score = scan$overdispersion_score[j])
    )
    marginal <- if (isTRUE(scan$converged[j])) list(
      feature_id = scan$feature_id[j], statistic = scan$statistic[j],
      p_value = scan$p_value[j],
      error_message = if (is.na(scan$p_value[j])) scan$error[j] else NA_character_
    ) else NULL
    list(null = null, marginal = marginal)
  })
}

# What the manager keeps of one spatial fit: the compact estimates, never an
# observation-length vector.
.inlast_spatial_record <- function(z, spec, scores_a, score_error, mu_bar,
                                   retain_smooth, diagnostics, estimation) {
  record <- list(
    dispersion = z$dispersion, family = z$family,
    family_parameters = z$family_parameters,
    smoothing_parameters = z$smoothing_parameters,
    target = as.numeric(z$random_mode[[1L]]),
    nuisance = if (ncol(spec$nuisance_design)) .inlast_nuisance_mode(z, spec) else
      numeric(),
    converged = isTRUE(z$converged),
    log_marginal_likelihood = z$log_marginal_likelihood,
    fit_seconds = z$fit_seconds, null_size = z$null_size, gradient = z$gradient,
    spatial_fallback = z$spatial_fallback,
    constraint_residual = z$constraint_residual,
    observation_spatial_mean = z$observation_spatial_mean,
    score_a = scores_a, score_error = score_error, mu_bar = mu_bar
  )
  if (estimation) record$estimation <- z$estimation
  if (retain_smooth) record$coefficients <- z$coefficients
  if (diagnostics) {
    record$diagnostics <- z[c(
      "tau", "tau_internal", "precision_scale", "lambda", "constraint_residual", "constraint_residual_uncorrected",
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
      diagnostics, payload$poisson, threads, payload$null
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
#' This function is Stage 1, the method behind [spdeST.scan()]: it fits only
#' the null model (no spatial field) of every feature and computes the Stage 1
#' null-first p-value of its marginal score test inside the worker; the
#' p-values are then adjusted by `adjust` into the Stage 1 q-values
#' `diagnostics$marginal_q_value`. Stage 2, [inlaST.estimate_spatial()],
#' fits the spatial model of the chosen features and computes their score vector and the mean of the fitted mean in the worker.
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
#'   Without an `initial` value in `nb_size_prior`, the NB size of a feature
#'   starts at the size of its Stage 1 null model: from `log(size) = 0` the
#'   search can leave a weakly overdispersed feature for an unbounded size.
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
#'   strategy and EB integration enforced. With a negative-binomial family each
#'   feature first gets its offset-and-covariate-only Poisson GLM, and it is
#'   fitted with the Poisson family when its overdispersion score
#'   `sum((y - mu)^2 - y) / 2` is not positive, i.e. exactly when its NB2
#'   maximum-likelihood fit is the Poisson fit (R/family-prescreen.R; there is
#'   no threshold to set). Stage 1 is [spde.scan()] on the raw sparse field
#'   of the target block, the same function [mgcvST.estimate()] calls; the
#'   null must have fixed effects only. The
#'   per-feature Pearson dispersion of that Poisson fit and the family actually
#'   used are reported in the diagnostics as
#'   `prescreen_phi` and `family_used`. Unknown controls are rejected.
#'   `engine` selects the Stage 2 fitter. `"spde"`, the default, is
#'   [spde.estimate()]: count families, fixed effects and one SPDE field, no
#'   prior and no fixed precision or size; `spde` is then a named list of its
#'   controls, for example `list(tol = 1e-8)`, and the INLA package is not
#'   needed. `"inla"` fits with INLA and takes the other controls of this
#'   list. `options(spdeST.engine = )` sets the default of a session.
#' @param retain_smooth Retain estimated score-component coefficients.
#' @param diagnostics Retain the per-feature INLA optimizer/hyperparameter
#'   diagnostics in `inla_diagnostics`.
#' @param threads OpenMP threads per worker for the sparse marginal and score
#'   computations of the workers.
#' @param adjust Multiple-testing adjustment of the Stage 1 p-values:
#'   `"BY"` (the default), `"BH"`, `"Sidak"` or `"none"`.
#' @param q.value Stage 1 discovery threshold in `(0, 1]`.
#' @param checkpoint_dir Optional directory for the resumable chunk checkpoints
#'   of both steps; see Details.
#' @param resume Reuse the completed chunks of a compatible checkpoint.
#' @param scan_control Named list of [spde.scan()] controls.
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
#' estimation. Stage 1 is [spde.scan()]: GLM nulls and the marginal score
#' test (`diagnostics$marginal_p_value`) calibrated by the saddlepoint
#' approximation (Kuonen 1999). Pair tests in [inlaST.test()] are
#' calibrated by a saddlepoint approximation. Every spatial component's constraint
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
#' `diagnostics$spatial_fallback` identifies these usable fallback states and
#' the original error remains in `error_class`, `error_message` and `error_call`. Their
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
#' @keywords internal
inlaST.estimate <- function(
    Y, model, feature_id = rownames(Y),
    BPPARAM = BiocParallel::SerialParam(), chunk_size = NULL,
    offset = NULL, control = list(), retain_smooth = FALSE,
    diagnostics = FALSE, threads = 1L,
    adjust = c("BY", "BH", "Sidak", "none"), q.value = 0.05,
    checkpoint_dir = NULL, resume = TRUE, scan_control = list()) {
  scan_args <- .spdest_scan_control(scan_control)
  adjust <- match.arg(adjust)
  q.value <- .mgcvst_check_q_value(q.value)
  threads <- .inlast_check_threads(threads)
  # One BLAS thread in this process before any worker is forked, so that a
  # forked worker never inherits a threaded BLAS; restored on exit.
  old_threads <- .mgcvst_thread_limit()
  on.exit(.mgcvst_thread_restore(old_threads), add = TRUE)
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
  if (!identical(control$engine, "spde") && !requireNamespace("INLA", quietly = TRUE)) {
    stop("inlaST.estimate() requires the INLA package unless control$engine = 'spde'.")
  }
  control <- .inlast_family_control(model, control)
  if (identical(control$engine, "spde")) .spde_estimate_check(model$inla_spec, control)
  Y <- checked$Y
  feature_id <- checked$feature_id
  chunk_size <- checked$chunk_size
  p <- nrow(Y)
  # Exact Poisson routing (R/family-prescreen.R): a feature whose NB2
  # maximum-likelihood fit is the Poisson fit (overdispersion score <= 0) is
  # fitted with the Poisson family. The family is fixed inside the spec, so a routed
  # feature is fitted against a Poisson twin of the same spec
  # (.inlast_poisson_spec()); the twin is rebuilt once per chunk in the worker
  # rather than shipped, because it shares every sparse block with the original.
  screen_offset <- if (is.null(offset)) model$offset else
    if (is.matrix(offset)) sweep(offset, 2L, model$offset, "+") else
      model$offset + offset
  prescreen <- .mgcvst_prescreen_route(
    Y, X = model$inla_spec$fixed$X, offset = screen_offset,
    active = identical(model$inla_spec$family, "negative_binomial"),
    threads = threads
  )
  null_spec <- .inlast_null_spec(model$inla_spec)
  .inlast_scan_check(null_spec)
  # Validate the twin once in the parent so a spec the switch cannot handle
  # fails here rather than inside every affected worker.
  if (any(prescreen$poisson)) .inlast_poisson_spec(model$inla_spec)
  family_used <- rep(model$inla_spec$family, p)
  family_used[prescreen$poisson] <- "poisson"
  null_control <- .inlast_null_control(control, model$inla_spec)
  if (any(prescreen$poisson)) .inlast_poisson_spec(null_spec)

  t0 <- proc.time()[["elapsed"]]
  y_digest <- .mgcvst_row_digests(Y)
  signature <- .inlast_estimation_signature(
    model, control, offset, feature_id, retain_smooth, diagnostics
  )
  signature <- .mgcvst_pair_input_hash(list(signature, scan_args = scan_args))
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
                threads = threads, scan_args = scan_args),
    libpaths = .libPaths(), manager = Sys.getpid()
  )
  null_elapsed <- proc.time()[["elapsed"]] - null_t0
  nulls <- unlist(step1$results, recursive = FALSE)

  n_sp <- model$inla_spec$geometry_sp_length
  m <- ncol(score_sparse$A)
  px <- ncol(model$geometry$nuisance_design)
  null_error_class <- null_error_message <- null_error_call <- rep(NA_character_, p)
  null_state <- stats::setNames(vector("list", p), feature_id)
  marginal <- data.frame(
    p_value = rep(NA_real_, p), error_message = NA_character_,
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
    if (!is.null(z$marginal)) {
      marginal$p_value[j] <- z$marginal$p_value
      marginal$error_message[j] <- z$marginal$error_message
    }
  }
  rm(nulls)
  marginal_q <- .mgcvst_stage1_q(marginal$p_value, adjust)

  diagnostics_table <- data.frame(
    index = seq_len(p), feature_id = feature_id, converged = FALSE,
    marginal_p_value = marginal$p_value, marginal_q_value = marginal_q,
    residual_df = NA_real_, criterion = NA_real_,
    criterion_name = if (identical(control$engine, "spde"))
      "spde.estimate Laplace marginal likelihood" else "INLA log marginal likelihood",
    fit_seconds = NA_real_,
    null_size = vapply(null_state, function(z) {
      if (is.null(z)) NA_real_ else if (length(z$family_parameters)) z$family_parameters[1L] else Inf
    }, numeric(1L)),
    spatial_size = NA_real_,
    gradient_log_tau = NA_real_, gradient_log_alpha = NA_real_,
    error_class = NA_character_,
    error_message = NA_character_, error_call = NA_character_,
    spatial_selected = FALSE, spatial_fitted = FALSE,
    spatial_fallback = FALSE, spatial_fallback_method = NA_character_,
    spatial_precision_assigned = NA_real_,
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
    stage1 = list(adjust = adjust, q.value = q.value,
                  null_engine = "spde.scan"),
    call = match.call()
  ), class = c("inlaST_fit", "mgcvST_model_fit", "mgcvST_fit", "spdeST"))
  kappa <- .spde_kappa_fields(model)
  ans[names(kappa)] <- kappa
  rm(step1)
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
#' @keywords internal
inlaST.estimate_spatial <- function(
    fitinlaST, Y, features = "discoveries",
    adjust = c("BY", "BH", "Sidak", "none"), q.value = 0.05,
    BPPARAM = BiocParallel::SerialParam(), chunk_size = NULL,
    checkpoint_dir = NULL, resume = TRUE, threads = 1L) {
  fit <- fitinlaST
  if (!inherits(fit, "inlaST_fit")) {
    stop("fitinlaST must be returned by inlaST.estimate().")
  }
  if (!identical(fit$control$engine, "spde") &&
      !requireNamespace("INLA", quietly = TRUE)) {
    stop("inlaST.estimate_spatial() requires the INLA package.")
  }
  .mgcvst_check_fit_format(fit)
  adjust <- if (missing(adjust)) fit$stage1$adjust else match.arg(adjust)
  q.value <- if (missing(q.value)) fit$stage1$q.value else
    .mgcvst_check_q_value(q.value)
  threads <- .inlast_check_threads(threads)
  # One BLAS thread in this process before any worker is forked, so that a
  # forked worker never inherits a threaded BLAS; restored on exit.
  old_threads <- .mgcvst_thread_limit()
  on.exit(.mgcvst_thread_restore(old_threads), add = TRUE)
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
  # spde.estimate() fits the features of a chunk in parallel threads, so a
  # checkpointed chunk holds several features per thread.
  groups <- .mgcvst_feature_chunks(
    index, chunk_size, BPPARAM, checkpoint = !is.null(store),
    cap = if (identical(fit$control$engine, "spde")) max(50L, 4L * threads) else 50L)
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
    libpaths = .libPaths(), manager = Sys.getpid()
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
    # The family of the spatial fit: spde.estimate() returns the Poisson fit
    # of a feature whose criterion does not rise along 1 / size.
    fit$feature_family[j] <- diag$family_used[j] <- z$family
    diag$spatial_size[j] <- if (length(z$family_parameters)) z$family_parameters[1L] else Inf
    if (length(z$gradient) == 2L) {
      diag$gradient_log_tau[j] <- z$gradient[[1L]]
      diag$gradient_log_alpha[j] <- z$gradient[[2L]]
    }
    fit$smoothing_parameters[j, ] <- z$smoothing_parameters
    fit$target_coefficients[, j] <- z$target
    if (length(z$nuisance)) fit$nuisance_coefficients[, j] <- z$nuisance
    fit$score_a[, j] <- z$score_a
    fit$mu_bar[j] <- z$mu_bar
    diag$converged[j] <- z$converged
    diag$criterion[j] <- z$log_marginal_likelihood
    diag$fit_seconds[j] <- z$fit_seconds
    diag$spatial_fitted[j] <- TRUE
    if (!is.null(z$spatial_fallback)) {
      fallback <- z$spatial_fallback
      diag$spatial_fallback[j] <- TRUE
      diag$spatial_fallback_method[j] <- fallback$method
      diag$spatial_precision_assigned[j] <- fallback$spatial_precision
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
