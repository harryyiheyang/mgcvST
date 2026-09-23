# Resolve the score engine recorded by the fitted object.
.mgcvst_test_engine <- function(fit) {
  engine <- fit$test_engine
  if (is.null(engine)) {
    stop("fitmgcvST does not record a score-test engine.")
  }
  engine <- as.character(engine)
  registry <- c(
    spde = ".mgcvst_test_spde",
    single_model = ".mgcvst_test_model_single"
  )
  target <- unname(registry[engine])
  if (length(target) != 1L || is.na(target)) {
    stop("fitmgcvST contains an unknown score-test engine: ", engine)
  }
  get(target, envir = environment(.mgcvst_test_engine), inherits = TRUE)
}

# Construct each requested model score state once and write one packed shard.
.mgcvst_model_state_shard <- function(features, fit, paths, threads = 1L,
                                      native = NULL) {
  .mgcvst_thread_limit()
  if (!is.null(native)) {
    for (first in seq.int(1L, length(features), by = 32L)) {
      rows <- first:min(length(features), first + 31L)
      ids <- features[rows]
      phi <- fit$dispersion[ids]
      sp <- fit$smoothing_parameters[ids, , drop = FALSE]
      bad <- !is.finite(phi) | phi <= 0 |
        rowSums(!is.finite(sp) | sp <= 0) > 0L
      units <- mgcvst_dense_score_batch_cpp(
        native$T0, fit$working_variance[, ids, drop = FALSE],
        fit$working_error[, ids, drop = FALSE],
        fit$dispersion[ids] / fit$smoothing_parameters[ids, native$sp_index],
        native$X, fit$nuisance_covariance[ids], threads
      )
      for (k in seq_along(ids)) {
        z <- units[[k]]
        if (bad[k]) z <- list(error =
          "The feature has invalid dispersion or smoothing parameters.")
        unit <- if (is.null(z$error)) {
          .mgcvst_pack_score_state(list(a = z$a, M = z$H, width = native$width))
        } else list(error = z$error)
        saveRDS(unit, paths[rows[k]])
      }
    }
    return(features)
  }
  for (k in seq_along(features)) {
    z <- tryCatch(.mgcvst_model_score_state(fit, features[k]),
                  error = function(e) e)
    if (inherits(z, "condition")) {
      unit <- list(error = conditionMessage(z))
    } else {
      unit <- .mgcvst_pack_score_state(z)
    }
    saveRDS(unit, paths[k])
  }
  features
}

# The current mgcv model has one marked SPDE and conditional nuisance covariance.
.mgcvst_model_dense_preparation <- function(fit, features) {
  geometry <- fit$geometry
  if (identical(fit$score_backend, "sparse") || length(geometry$target) != 1L ||
      is.null(geometry$nuisance_design) ||
      length(fit$nuisance_covariance) < max(features) ||
      any(vapply(fit$nuisance_covariance[features], is.null, logical(1L)))) {
    return(NULL)
  }
  j <- unname(geometry$target[[1L]])
  s <- geometry$smooth[[j]]
  T0 <- fit$.mgcvst_fixed_factors[[j]]
  if (s$fixed || length(s$sp_index) != 1L || !is.matrix(T0)) return(NULL)
  list(T0 = T0, X = geometry$nuisance_design, sp_index = s$sp_index,
       width = stats::setNames(ncol(T0), names(geometry$target)))
}

# Shared orchestration for model.set() score engines.
.mgcvst_test_model <- function(
    fitmgcvST, pair_function, test_definition,
    q.value = 0.05, FDR = TRUE, method = "BH",
    BPPARAM = BiocParallel::SerialParam(), ...,
    pairs = NULL, highlight = NULL,
    calibration = c("liu", "davies"),
    chunk_size = NULL,
    threads = NULL, verbose = FALSE, cache_bytes = NULL,
    checkpoint_dir = NULL, resume = TRUE, approximate = FALSE,
    n_ref = 100L, ref_method = c("random", "score", "hyper"),
    ref_seed = 1L, ref_tol = 1e-6,
    diagnostic_pairs = 0L) {
  if (!inherits(fitmgcvST, "mgcvST_model_fit")) {
    stop("The model score engine requires a fit from mgcvST.estimate(Y, model).")
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
  # INLA fits use the sparse score kernel, so the INLA predicate selects the
  # corresponding downstream execution path.
  inla_fit <- .mgcvst_inla_downstream(fitmgcvST)
  if (inla_fit && calibration != "liu") {
    stop("INLA downstream tests support calibration = 'liu' only.")
  }
  if (inla_fit) {
    .mgcvst_inla_serial_backend(BPPARAM)
    .mgcvst_inla_require_sparse(fitmgcvST)
  }
  if (calibration == "davies" &&
      !requireNamespace("CompQuadForm", quietly = TRUE)) {
    stop("calibration = 'davies' requires the optional CompQuadForm package.")
  }
  if (is.null(threads)) {
    threads <- if (inla_fit) 1L else BiocParallel::bpworkers(BPPARAM)
  }
  threads <- as.integer(threads)
  if (length(threads) != 1L || is.na(threads) || threads < 1L) {
    stop("threads must be one positive integer.")
  }
  .mgcvst_thread_limit()
  if (calibration != "liu" && (!is.null(cache_bytes) ||
      !is.null(checkpoint_dir) || !isTRUE(resume))) {
    stop("cache_bytes, checkpoint_dir and resume currently require calibration = 'liu'.")
  }

  index <- .mgcvst_pair_index(pairs, fitmgcvST$feature_id)
  highlight_index <- matrix(integer(), nrow = 0L, ncol = 2L)
  if (!is.null(highlight)) {
    highlight_index <- .mgcvst_pair_index(highlight, fitmgcvST$feature_id)
  }
  n_feature <- length(fitmgcvST$feature_id)
  key <- (index[, 1L] - 1L) * n_feature + index[, 2L]
  highlight_key <- if (nrow(highlight_index)) {
    (highlight_index[, 1L] - 1L) * n_feature + highlight_index[, 2L]
  } else {
    numeric()
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
    statistic = NA_real_, information = NA_real_,
    effective_rank = NA_real_, p_two_sided = NA_real_,
    p_positive = NA_real_, p_negative = NA_real_,
    p_adjusted = NA_real_, p_positive_adjusted = NA_real_,
    p_negative_adjusted = NA_real_, discovered = FALSE,
    discovered_positive = FALSE, discovered_negative = FALSE,
    highlighted = highlighted, retained = highlighted,
    error_message = NA_character_,
    stringsAsFactors = FALSE
  )
  unavailable_rows <- which(!pair_available)
  for (k in unavailable_rows) {
    missing_feature <- c(i1[k], i2[k])[!available[c(i1[k], i2[k])]]
    result$error_message[k] <- paste(
      paste0(
        fitmgcvST$feature_id[missing_feature], ": ",
        fitmgcvST$diagnostics$error_message[missing_feature]
      ),
      collapse = " | "
    )
  }

  tested_rows <- which(pair_available)
  workers <- if (length(tested_rows)) {
    max(1L, min(length(tested_rows), BiocParallel::bpworkers(BPPARAM)))
  } else {
    0L
  }
  inla_projection <- NULL
  inla_basis_elapsed <- 0
  inla_test_started <- NULL
  if (inla_fit && length(tested_rows)) {
    inla_test_started <- proc.time()[["elapsed"]]
    fitmgcvST <- .inlast_sparse_prepare(fitmgcvST)
    t_basis <- proc.time()[["elapsed"]]
    inla_projection <- .inlast_sparse_observation_basis(fitmgcvST)
    inla_basis_elapsed <- proc.time()[["elapsed"]] - t_basis
  }
  if (is.null(chunk_size)) {
    chunk_size <- if (approximate) 10000L else if (inla_fit) .mgcvst_inla_pair_chunk_size(
      fitmgcvST, basis = inla_projection
    ) else if (calibration == "liu") 10000L else if (workers > 0L)
      ceiling(length(tested_rows) / workers) else 1L
  }
  chunk_size <- as.integer(chunk_size)
  if (length(chunk_size) != 1L || is.na(chunk_size) || chunk_size < 1L) {
    stop("chunk_size must be one positive integer.")
  }
  chunks <- list()
  chunk_count <- 0L
  elapsed <- summary_elapsed <- 0
  native_preparation <- FALSE
  pipeline <- NULL
  if (length(tested_rows)) {
    chunks <- if (inla_fit) {
      NULL
    } else split(tested_rows, ceiling(seq_along(tested_rows) / chunk_size))
    if (inla_fit) chunk_count <- ceiling(length(tested_rows) /
      min(chunk_size, 128L))
    if (inla_fit) {
      t0 <- proc.time()[["elapsed"]]
      evaluated <- .mgcvst_inla_test_pairs(
        fitmgcvST, index[tested_rows, , drop = FALSE], tested_rows,
        threads, chunk_size, verbose, basis = inla_projection,
        cache_bytes = cache_bytes, checkpoint_dir = checkpoint_dir,
        resume = resume, approximate = approximate, n_ref = n_ref,
        ref_method = ref_method, ref_seed = ref_seed, ref_tol = ref_tol,
        diagnostic_pairs = diagnostic_pairs
      )
      elapsed <- proc.time()[["elapsed"]] - t0
      inla_projection <- attr(evaluated$result, "inla_pairwise")
      inla_projection$basis_elapsed <- inla_basis_elapsed
      inla_projection$test_wall_elapsed <- proc.time()[["elapsed"]] -
        inla_test_started
      chunk_count <- inla_projection$chunks
      summary_elapsed <- inla_projection$preparation_elapsed
      target <- evaluated$result$pair_index
      columns <- intersect(names(evaluated$result), names(result))
      result[target, columns] <- evaluated$result[, columns, drop = FALSE]
      result$statistic[target] <- result$signed_score[target]^2
      evaluated <- list()
    } else if (calibration == "liu") {
      evaluate <- if (approximate) .mgcvst_pair_approximate else .mgcvst_pair_pipeline
      args <- list(
        fitmgcvST, index[tested_rows, , drop = FALSE], tested_rows,
        threads, chunk_size, verbose, cache_bytes = cache_bytes,
        checkpoint_dir = checkpoint_dir, resume = resume
      )
      if (approximate) args <- c(args, list(n_ref = n_ref, ref_method = ref_method,
        ref_seed = ref_seed, ref_tol = ref_tol,
        diagnostic_pairs = diagnostic_pairs))
      evaluated <- do.call(evaluate, args)
      pipeline <- evaluated$metadata
      native_preparation <- pipeline$preparation_backend %in%
        c("sparse", "model_native", "legacy_native")
      summary_elapsed <- pipeline$preparation_elapsed
      elapsed <- summary_elapsed + evaluated$elapsed
      chunk_count <- pipeline$chunks
      target <- evaluated$result$pair_index
      result$signed_score[target] <- evaluated$result$score
      result$statistic[target] <- evaluated$result$score^2
      result$information[target] <- evaluated$result$information
      result$effective_rank[target] <- evaluated$result$effective_rank
      result$p_two_sided[target] <- evaluated$result$p_value
      result$error_message[target] <- evaluated$result$error_message
      evaluated <- list()
    } else {
    chunks <- .mgcvst_dense_pair_groups(tested_rows, index, chunk_size)
    used <- sort(unique(as.vector(index[tested_rows, , drop = FALSE])))
    feature_workers <- max(1L, min(length(used), BiocParallel::bpworkers(BPPARAM)))
    feature_groups <- split(used, ceiling(seq_along(used) /
      ceiling(length(used) / feature_workers)))
    cache_dir <- .mgcvst_dense_temp_dir()
    on.exit(.mgcvst_dense_cleanup(cache_dir), add = TRUE)
    worker_bundle <- .mgcvst_worker_bundle()
    state_shard <- get(".mgcvst_model_state_shard", envir = worker_bundle,
                       inherits = FALSE)
    test_chunk <- get(".mgcvst_dense_pair_chunk", envir = worker_bundle,
                      inherits = FALSE)
    t0 <- proc.time()[["elapsed"]]
    test_fit <- fitmgcvST
    test_fit$.mgcvst_fixed_factors <- .mgcvst_model_fixed_factors(test_fit)
    native <- .mgcvst_model_dense_preparation(test_fit, used)
    native_preparation <- !is.null(native)
    shard_paths <- file.path(cache_dir, paste0("state-", used, ".rds"))
    names(shard_paths) <- as.character(used)
    if (native_preparation) {
      .mgcvst_model_state_shard(used, test_fit, shard_paths, threads, native)
    } else {
    BiocParallel::bplapply(
      seq_along(feature_groups), function(k, groups, paths, fit, worker_fun) {
        feature <- groups[[k]]
        worker_fun(feature, fit, paths[as.character(feature)])
      }, groups = feature_groups, paths = shard_paths, fit = test_fit,
      worker_fun = state_shard, BPPARAM = BPPARAM
    )
    }
    summary_elapsed <- proc.time()[["elapsed"]] - t0
    payload <- lapply(chunks, function(rows) {
      pair <- index[rows, , drop = FALSE]
      feature <- sort(unique(as.vector(pair)))
      list(
        rows = rows, pairs = pair,
        shards = shard_paths[as.character(feature)]
      )
    })
    evaluated <- BiocParallel::bplapply(
      payload, test_chunk, calibration = calibration, BPPARAM = BPPARAM
    )
    elapsed <- proc.time()[["elapsed"]] - t0
    evaluated <- unlist(evaluated, recursive = FALSE)
    }
    for (z in evaluated) {
      target <- z$pair_index
      names <- intersect(names(z), names(result))
      result[target, names] <- z[1L, names, drop = FALSE]
      result$statistic[target] <- result$signed_score[target]^2
    }
  }

  valid <- is.finite(result$p_two_sided) &
    result$p_two_sided >= 0 & result$p_two_sided <= 1
  result$p_positive[valid] <- ifelse(result$signed_score[valid] >= 0,
    result$p_two_sided[valid] / 2, 1 - result$p_two_sided[valid] / 2)
  result$p_negative[valid] <- ifelse(result$signed_score[valid] <= 0,
    result$p_two_sided[valid] / 2, 1 - result$p_two_sided[valid] / 2)
  if (FDR) {
    result$p_adjusted[valid] <- stats::p.adjust(result$p_two_sided[valid], method)
    result$p_positive_adjusted[valid] <- stats::p.adjust(
      result$p_positive[valid], method
    )
    result$p_negative_adjusted[valid] <- stats::p.adjust(
      result$p_negative[valid], method
    )
  } else {
    result$p_adjusted[valid] <- result$p_two_sided[valid]
    result$p_positive_adjusted[valid] <- result$p_positive[valid]
    result$p_negative_adjusted[valid] <- result$p_negative[valid]
  }
  result$discovered <- valid & result$p_adjusted <= q.value
  result$discovered_positive <- valid & result$p_positive_adjusted <= q.value
  result$discovered_negative <- valid & result$p_negative_adjusted <= q.value
  result$retained <- result$highlighted | result$discovered
  raw_threshold <- if (any(result$discovered)) {
    max(result$p_two_sided[result$discovered])
  } else {
    NA_real_
  }
  timing <- list(
    elapsed = elapsed, summary_elapsed = summary_elapsed,
    pair_elapsed = elapsed - summary_elapsed,
    workers = if (calibration == "liu") threads else workers,
    chunks = if (calibration == "liu") chunk_count else length(chunks),
    backend = if (calibration == "liu") "C++ OpenMP" else class(BPPARAM)[1L],
    preparation_backend = if (inla_fit || native_preparation)
      "C++ OpenMP" else class(BPPARAM)[1L],
    preparation_threads = if (inla_fit || native_preparation) threads else workers
  )
  if (inla_fit) timing$inla_projection <- inla_projection
  if (!is.null(pipeline)) timing$pair_pipeline <- pipeline
  structure(
    list(
      results = result,
      threshold = list(
        q_value = q.value, FDR = FDR,
        adjustment_method = if (FDR) method else "none",
        raw_p_threshold = raw_threshold
      ),
      discoveries = list(
        pairs_requested = nrow(index),
        pairs_tested = length(tested_rows), pairs_with_p_value = sum(valid),
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
      test_definition = test_definition,
      timing = timing,
      calibration = calibration,
      call = match.call()
    ),
    class = "mgcvST_test"
  )
}

# Single marked-SPDE entry point.
.mgcvst_test_model_single <- function(...) {
  .mgcvst_test_model(
    ..., pair_function = .mgcvst_model_pair_single,
    test_definition = "single_global_cross_gene_covariance_at_independence"
  )
}

#' Test cross-feature spatial covariance
#'
#' Dispatches a compact fit to its registered score engine. Standard one-SPDE
#' fits use the SPDE score path. One-component fits constructed from
#' [model.set()] use the model score path.
#'
#' @inheritParams .mgcvst_test_spde
#' @param pairwise_method `"liu"` keeps the existing pair test. For sparse INLA
#'   fits, `"conditional"` uses both conditional-normal directions, combines
#'   their p-values by the equal-weight Cauchy rule, and applies BY across the
#'   tested pair family. With `pairs = NULL`, it tests every available gene pair.
#' @param conditional_precision Precision used for conditional variance
#'   multiplication: `"double"` or `"float32"`. Scores, p-values, and BY
#'   adjustment remain in double precision.
#' @param cache_bytes Optional byte ceiling for resident Liu score states.
#'   The default adapts to available system and job memory, with space reserved
#'   for native working buffers. This is a cache budget, not a process limit.
#' @param checkpoint_dir Optional checkpoint directory. Conditional testing
#'   saves per-gene variance rows; Liu testing saves reusable score states and
#'   pair batches. With `NULL`, temporary storage is removed on exit.
#' @param resume Reuse compatible completed checkpoint entries.
#' @param approximate Use real-gene landmark CUR trace reconstruction with
#'   Liu calibration. Landmark-to-feature traces use float32 matrix products
#'   with double accumulation and storage; landmark block `W` uses double.
#'   Approximate p-values remain approximate throughout the test and are used
#'   for multiple-testing adjustment. The default `FALSE` computes exact traces
#'   in the existing common coordinates.
#' @param n_ref Maximum number of real genes used as landmarks in approximate
#'   mode.
#' @param ref_method Landmark selection: uniform random sampling, k-means on
#'   unnormalized score vectors (`"score"`), or k-means on standardized fitted
#'   covariance variance scales (`"hyper"`).
#' @param ref_seed Non-negative integer seed for reference selection; the
#'   caller's random-number state is restored.
#' @param ref_tol Relative eigenvalue cutoff for each normalized reference trace
#'   matrix `W`. Both positive and negative retained eigenvalues are inverted.
#' @param diagnostic_pairs Number of pairs with two non-landmark endpoints
#'   sampled for optional exact approximation diagnostics. Defaults to `0`,
#'   so no exact diagnostics run. Positive values add diagnostics without
#'   changing any returned pair p-value or adjusted p-value.
#' @details Let `S_ij = a_i' a_j` and `v_(i|j) = a_j' M_i a_j`.
#'   Under independent Gaussian null scores, `S_ij | a_j` is normal with
#'   variance `v_(i|j)`, so each directional two-sided normal p-value is exactly
#'   uniform. The directions are combined using
#'   `T = (tan((0.5-p_(i|j))*pi) + tan((0.5-p_(j|i))*pi))/2` and the standard
#'   Cauchy upper tail. Because `T` cannot exceed its larger component,
#'   `p_ij >= min(p_(i|j), p_(j|i))` and the null rejection probability is at
#'   most `2*alpha` under any dependence.
#'
#'   For `X = 1/p_(i|j)` and `Y = 1/p_(j|i)`, both directional tails satisfy
#'   `Pr(X > x) = Pr(Y > x) = 1/x`; in the far tail, `1/p_ij` approaches
#'   `(X+Y)/2`. If both directions become extreme together, their normal
#'   z-scores obey `z_2 = R*z_1`, where `R = sqrt(v_(i|j)/v_(j|i))`.
#'   Under a continuous, nondegenerate distribution of `R`, unequal extremes
#'   occur together only when `abs(R-1)` is of order `1/log(1/alpha)`.
#'   Thus the combined tail approaches the nominal tail as `alpha` tends to
#'   zero under this condition.
#'   Conditional output includes log-scale p-values so tails below the
#'   floating-point range remain available for BY calculations and reporting.
#' @export
mgcvST.test <- function(
    fitmgcvST, q.value = 0.05, FDR = TRUE, method = "BH",
    BPPARAM = BiocParallel::SerialParam(), ...,
    pairs = NULL, highlight = NULL,
    calibration = c("liu", "davies"),
    chunk_size = NULL,
    threads = NULL, verbose = FALSE, cache_bytes = NULL,
    checkpoint_dir = NULL, resume = TRUE, approximate = FALSE,
    n_ref = 100L, ref_method = c("random", "score", "hyper"),
    ref_seed = 1L, ref_tol = 1e-6,
    diagnostic_pairs = 0L,
    pairwise_method = c("liu", "conditional"),
    conditional_precision = c("double", "float32")) {
  pairwise_method <- match.arg(pairwise_method)
  conditional_precision <- match.arg(conditional_precision)
  if (!is.logical(approximate) || length(approximate) != 1L || is.na(approximate)) {
    stop("approximate must be TRUE or FALSE.")
  }
  if (pairwise_method == "conditional") {
    if (missing(method)) method <- "BY"
    if (!identical(method, "BY")) {
      stop("pairwise_method = 'conditional' requires method = 'BY'.")
    }
    if (!isTRUE(FDR)) {
      stop("pairwise_method = 'conditional' requires FDR = TRUE.")
    }
    if (!is.null(highlight)) {
      stop("highlight is unavailable for conditional pairwise results.")
    }
    if (!identical(calibration, c("liu", "davies")) &&
        !identical(calibration, "liu")) {
      stop("Conditional pairs do not use a non-Liu calibration argument.")
    }
    if (approximate || !is.null(cache_bytes)) {
      stop("approximate and cache_bytes require pairwise_method = 'liu'.")
    }
    if (length(list(...))) stop("Unused arguments in ... for conditional pairs.")
    .mgcvst_inla_serial_backend(BPPARAM)
    return(.mgcvst_conditional_test(
      fitmgcvST, pairs, q.value, threads, chunk_size,
      checkpoint_dir, resume, conditional_precision, match.call()
    ))
  }
  if (!identical(conditional_precision, "double")) {
    stop("conditional_precision requires pairwise_method = 'conditional'.")
  }
  calibration <- match.arg(calibration)
  if (approximate && calibration != "liu") {
    stop("approximate = TRUE requires calibration = 'liu'.")
  }
  if (approximate) ref_method <- match.arg(ref_method)
  engine <- .mgcvst_test_engine(fitmgcvST)
  args <- list(
    fitmgcvST = fitmgcvST, q.value = q.value, FDR = FDR, method = method,
    BPPARAM = BPPARAM, ..., pairs = pairs, highlight = highlight,
    calibration = calibration, chunk_size = chunk_size,
    threads = threads, verbose = verbose, cache_bytes = cache_bytes,
    checkpoint_dir = checkpoint_dir, resume = resume, approximate = approximate,
    n_ref = n_ref, ref_method = ref_method, ref_seed = ref_seed, ref_tol = ref_tol,
    diagnostic_pairs = diagnostic_pairs
  )
  do.call(engine, args)
}
