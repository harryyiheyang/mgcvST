# Bound serialization buffers while hashing every input value and its attributes.
.mgcvst_pair_input_hash <- function(x) {
  if (is.list(x) && !isS4(x)) {
    parts <- vapply(x, .mgcvst_pair_input_hash, character(1L))
    return(digest::digest(list(type = typeof(x), attributes = attributes(x),
                               parts = parts), algo = "md5"))
  }
  if (is.atomic(x) && length(x) > 1048576L) {
    starts <- seq.int(1, length(x), by = 1048576L)
    parts <- vapply(starts, function(first) {
      last <- min(length(x), first + 1048575L)
      digest::digest(x[seq.int(first, last)], algo = "md5")
    }, character(1L))
    return(digest::digest(list(type = typeof(x), length = length(x),
                               attributes = attributes(x), parts = parts),
                          algo = "md5"))
  }
  digest::digest(x, algo = "md5")
}

# Hash only inputs that determine feature scores and their common coordinates.
.mgcvst_pair_signature <- function(fit, basis = NULL) {
  geometry <- fit$geometry
  smooth <- if (is.list(geometry$smooth)) lapply(geometry$smooth, function(z) {
    z[c("B", "penalties", "sp_index", "fixed", "score_component")]
  }) else NULL
  sparse_backend <- identical(fit$score_backend, "sparse")
  sparse <- if (sparse_backend) {
    z <- fit$score_sparse
    list(A = z$A, Q = z$Q, constraint = z$constraint,
         sp_index = z$sp_index,
         nuisance_design = geometry$nuisance_design,
         basis = basis[c("coordinate", "basis", "rank", "coverage")],
         target = fit$target_coefficients, nuisance = fit$nuisance_coefficients,
         score = fit$score_a, family = fit$feature_family,
         family_parameters = fit$family_parameters, dispersion = fit$dispersion,
         smoothing_parameters = fit$smoothing_parameters)
  } else NULL
  spec <- if (is.list(fit$model)) fit$model$inla_spec else NULL
  if (is.null(spec)) spec <- fit$inla_spec
  random <- if (is.list(spec$random)) lapply(spec$random, function(z) {
    list(target = z$target, kind = z$kind, subtype = z$subtype,
         sp_index = z$sp_index,
         width = if (is.null(z$A)) NULL else ncol(z$A))
  }) else NULL
  inputs <- list(
    pipeline_version = 3L, estimator = fit$estimator,
    test_engine = fit$test_engine, score_backend = fit$score_backend,
    feature_id = fit$feature_id,
    working_error = if (sparse_backend) NULL else fit$working_error,
    working_variance = if (sparse_backend) NULL else fit$working_variance,
    dispersion = if (sparse_backend) NULL else fit$dispersion,
    lambda = fit$lambda,
    smoothing_parameters = if (sparse_backend) NULL else fit$smoothing_parameters,
    nuisance_covariance = if (sparse_backend) NULL else fit$nuisance_covariance,
    geometry = list(B = geometry$B, Q = geometry$Q, X = geometry$X,
                    score_precision_psd = geometry$score_precision_psd,
                    nuisance_design = geometry$nuisance_design,
                    target = geometry$target, smooth = smooth),
    sparse = sparse, random = random,
    inla_fixed_width = if (is.null(spec$fixed$X)) NULL else ncol(spec$fixed$X)
  )
  list(version = 3L, md5 = .mgcvst_pair_input_hash(inputs))
}

# Build one bounded feature batch with the existing full-precision score
# kernels. Callers only ever pass mode "model_native" or "legacy_native";
# the former per-feature R-loop fallback ("model_fallback", built on
# .mgcvst_model_score_state()) is deleted: a model.set() fit without a
# usable native dense preparation is a hard error in .mgcvst_pair_pipeline()
# below, never a fallback here.
.mgcvst_pair_build_batch <- function(fit, ids, threads, mode, native = NULL,
                                     T0 = NULL, field_scale = NULL) {
  if (identical(mode, "model_native")) {
    phi <- fit$dispersion[ids]
    sp <- fit$smoothing_parameters[ids, , drop = FALSE]
    bad <- !is.finite(phi) | phi <= 0 |
      rowSums(!is.finite(sp) | sp <= 0) > 0L
    z <- mgcvst_dense_score_batch_cpp(
      native$T0, fit$working_variance[, ids, drop = FALSE],
      fit$working_error[, ids, drop = FALSE],
      phi / sp[, native$sp_index], native$X,
      fit$nuisance_covariance[ids], threads
    )
    return(lapply(seq_along(ids), function(k) {
      if (bad[k]) return(list(error =
        "The feature has invalid dispersion or smoothing parameters."))
      if (!is.null(z[[k]]$error)) return(list(error = z[[k]]$error))
      list(a = z[[k]]$a, M = z[[k]]$H, width = native$width)
    }))
  }
  z <- mgcvst_dense_score_batch_cpp(
    T0, fit$working_variance[, ids, drop = FALSE],
    fit$working_error[, ids, drop = FALSE], field_scale[ids],
    fit$geometry$X, list(), threads
  )
  lapply(z, function(x) {
    if (!is.null(x$error)) list(error = x$error) else
      list(a = x$a, M = x$H, width = length(x$a))
  })
}

# Evaluate existing Liu calibration from resumable feature-first score states.
.mgcvst_pair_pipeline <- function(fit, index, pair_index, threads, chunk_size,
                                  verbose, cache_bytes = NULL,
                                  checkpoint_dir = NULL, resume = TRUE,
                                  state_store = NULL) {
  if (!is.matrix(index) || ncol(index) != 2L || !nrow(index) ||
      anyNA(index) || any(index != floor(index)) ||
      any(index < 1L) || any(index > length(fit$feature_id))) {
    stop("index must contain valid two-column feature indices.")
  }
  if (length(pair_index) != nrow(index) || anyNA(pair_index)) {
    stop("pair_index must identify every requested pair.")
  }
  if (!is.numeric(threads) || length(threads) != 1L ||
      !is.finite(threads) || threads < 1L || threads != floor(threads)) {
    stop("threads must be one positive integer.")
  }
  if (!is.numeric(chunk_size) || length(chunk_size) != 1L ||
      !is.finite(chunk_size) || chunk_size < 1L ||
      chunk_size != floor(chunk_size)) {
    stop("chunk_size must be one positive integer.")
  }
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("verbose must be TRUE or FALSE.")
  }
  if (!is.null(cache_bytes) && (!is.numeric(cache_bytes) ||
      length(cache_bytes) != 1L || !is.finite(cache_bytes) ||
      cache_bytes < 0)) {
    stop("cache_bytes must be NULL or a non-negative byte count.")
  }
  if (identical(fit$score_backend, "sparse")) {
    stop("The sparse INLA exact score_liu path no longer uses ",
         ".mgcvst_pair_pipeline(); it is served by the fp16 pipeline.")
  }
  used <- sort(unique(as.vector(index)))
  if (identical(fit$test_engine, "spde")) {
    width <- ncol(fit$geometry$B)
    mode <- "legacy_native"
  } else {
    target <- fit$geometry$target
    width <- if (length(target)) ncol(fit$geometry$smooth[[unname(target[[1L]])]]$B) else 1L
    mode <- "model_native"
  }
  if (!is.finite(width) || width < 1L) stop("The score coordinate width is invalid.")
  if (is.null(state_store)) {
    signature <- if (is.null(checkpoint_dir)) {
      list(version = 2L, temporary = tempfile("mgcvst-pair-run-"))
    } else .mgcvst_pair_signature(fit)
    store <- .mgcvst_store_open(checkpoint_dir, signature, fit$feature_id,
                                storage = "double", resume = resume)
    if (isTRUE(store$temporary)) on.exit(.mgcvst_store_cleanup(store), add = TRUE)
  } else {
    if (!inherits(state_store, "mgcvst_score_store") ||
        !identical(state_store$storage, "double") ||
        !identical(state_store$feature_ids, fit$feature_id)) {
      stop("An existing state_store must contain aligned double score states.")
    }
    store <- state_store
    signature <- store$signature
    if (!all(vapply(used, function(id) .mgcvst_store_has(store, id), logical(1L)))) {
      stop("The existing state_store is incomplete; exact rechecks cannot rebuild states.")
    }
  }
  pair_path <- .mgcvst_pair_checkpoint(store, index, pair_index)

  state_estimate <- 8 * (width^2 + width) + 2048
  automatic_cache <- is.null(cache_bytes)
  probe <- if (automatic_cache) .mgcvst_memory_probe() else NULL
  reserve <- 3 * 8 * width^2 * min(threads, chunk_size) +
    2 * state_estimate + 256 * nrow(index) + 64 * 1024^2
  if (automatic_cache) {
    available <- probe$available
    cache_bytes <- if (is.finite(available)) {
      max(0, 0.7 * available - reserve)
    } else 512 * 1024^2
  }

  existing <- vapply(used, function(id) .mgcvst_store_has(store, id), logical(1L))
  resume_count <- sum(existing)
  missing <- used[!existing]
  builds <- 0L
  preparation_started <- proc.time()[["elapsed"]]
  if (length(missing)) {
    native <- NULL
    T0 <- field_scale <- NULL
    if (identical(mode, "legacy_native")) {
      T0 <- .mgcvst_legacy_shared_score_factor(fit$geometry)
      field_scale <- .mgcvst_field_scale(fit)
    } else {
      fit$.mgcvst_fixed_factors <- .mgcvst_model_fixed_factors(fit)
      native <- .mgcvst_model_dense_preparation(fit, missing)
      if (is.null(native)) {
        stop("Model score states require the conditional nuisance covariance; ",
             "re-estimate with the current mgcvST.estimate().")
      }
    }
    q <- width
    p <- if (!is.null(native)) ncol(native$X) else if (mode == "legacy_native")
      ncol(fit$geometry$X) else ncol(fit$geometry$nuisance_design)
    if (is.null(p)) p <- 0L
    n <- nrow(fit$working_variance)
    feature_work <- 8 * (8 * n * q + 4 * n * p + 8 * q^2 + 4 * q * p + 4 * p^2)
    first <- 1L
    while (first <= length(missing)) {
      probe <- .mgcvst_memory_probe()
      headroom <- if (!is.null(probe) && is.finite(probe$available)) {
        0.3 * probe$available
      } else max(cache_bytes, 512 * 1024^2)
      if (headroom < feature_work) {
        stop("Insufficient available memory for one score-state preparation. ",
             "Increase the job memory allocation.")
      }
      batch_size <- as.integer(max(1L, min(32L, floor(headroom / feature_work))))
      ids <- missing[first:min(length(missing), first + batch_size - 1L)]
      states <- .mgcvst_pair_build_batch(
        fit, ids, threads, mode, native, T0, field_scale
      )
      if (length(states) != length(ids)) {
        stop("The score-state backend returned the wrong feature count.")
      }
      for (k in seq_along(ids)) .mgcvst_store_write(store, ids[k], states[[k]])
      builds <- builds + length(ids)
      rm(states)
      if (verbose) message("Stored score states for ", builds, " of ",
                           length(missing), " remaining features.")
      first <- first + length(ids)
    }
    rm(native, T0, field_scale)
    fit$.mgcvst_fixed_factors <- NULL
  }
  preparation_elapsed <- proc.time()[["elapsed"]] - preparation_started

  ord <- .mgcvst_pair_order(index, used,
    max(2L, floor(cache_bytes / state_estimate)), pair_path)
  index <- index[ord, , drop = FALSE]
  pair_index <- pair_index[ord]

  n_pairs <- nrow(index)
  r_score <- rep(NA_real_, n_pairs)
  r_information <- rep(NA_real_, n_pairs)
  r_effective_rank <- rep(NA_real_, n_pairs)
  r_p_value <- rep(NA_real_, n_pairs)
  r_log_p_two_sided <- rep(NA_real_, n_pairs)
  r_log_p_positive <- rep(NA_real_, n_pairs)
  r_log_p_negative <- rep(NA_real_, n_pairs)
  r_error_message <- rep(NA_character_, n_pairs)
  cache <- new.env(parent = emptyenv())
  cache$state <- list()
  cache$bytes <- 0
  cache$last <- numeric()
  cache$clock <- 0
  cache$hits <- 0L
  cache$misses <- 0L
  cache$evictions <- 0L
  elapsed <- 0
  chunks <- 0L
  resumed_pairs <- 0L
  first <- 1L
  last_probe <- proc.time()[["elapsed"]]
  while (first <= nrow(index)) {
    saved <- .mgcvst_pair_checkpoint_read(pair_path, first, pair_index)
    if (!is.null(saved)) {
      idx <- seq.int(first, saved$last)
      r_score[idx] <- saved$result$score
      r_information[idx] <- saved$result$information
      r_effective_rank[idx] <- saved$result$effective_rank
      r_p_value[idx] <- saved$result$p_value
      r_log_p_two_sided[idx] <- saved$result$log_p_two_sided
      r_log_p_positive[idx] <- saved$result$log_p_positive
      r_log_p_negative[idx] <- saved$result$log_p_negative
      r_error_message[idx] <- saved$result$error_message
      resumed_pairs <- resumed_pairs + nrow(saved$result)
      chunks <- chunks + 1L
      first <- saved$last + 1L
      rm(saved)
      next
    }
    if (automatic_cache &&
        proc.time()[["elapsed"]] - last_probe >= 2) {
      probe <- .mgcvst_memory_probe()
      if (is.finite(probe$available)) {
        cache_bytes <- max(0, 0.7 * (probe$available + cache$bytes) - reserve)
      }
      last_probe <- proc.time()[["elapsed"]]
    }
    max_features <- max(2L, min(length(used),
      floor(cache_bytes / state_estimate)))
    while (cache$bytes > cache_bytes && length(cache$last)) {
      drop <- names(cache$last)[which.min(cache$last)]
      cache$bytes <- cache$bytes - as.numeric(object.size(cache$state[[drop]]))
      cache$state[[drop]] <- NULL
      cache$last <- cache$last[names(cache$last) != drop]
      cache$evictions <- cache$evictions + 1L
    }
    cap <- as.integer(min(nrow(index), first + as.double(chunk_size) - 1))
    endpoints <- as.vector(t(index[seq.int(first, cap), , drop = FALSE]))
    first_seen <- which(!duplicated(endpoints))
    last <- if (length(first_seen) > max_features) {
      first + (first_seen[max_features + 1L] - 1L) %/% 2L - 1L
    } else cap
    active <- sort(unique(as.vector(index[seq.int(first, last), , drop = FALSE])))
    rows <- first:last
    states <- vector("list", length(active))
    for (k in seq_along(active)) {
      key <- as.character(active[k])
      if (!is.null(cache$state[[key]])) {
        cache$hits <- cache$hits + 1L
        states[[k]] <- cache$state[[key]]
      } else {
        cache$misses <- cache$misses + 1L
        state <- .mgcvst_store_read(store, active[k])
        size <- as.numeric(object.size(state))
        state_estimate <- max(state_estimate, size)
        while (cache$bytes + size > cache_bytes && length(cache$last)) {
          candidate <- setdiff(names(cache$last), as.character(active))
          if (!length(candidate)) break
          drop <- candidate[which.min(cache$last[candidate])]
          cache$bytes <- cache$bytes - as.numeric(object.size(cache$state[[drop]]))
          cache$state[[drop]] <- NULL
          cache$last <- cache$last[names(cache$last) != drop]
          cache$evictions <- cache$evictions + 1L
        }
        if (cache$bytes + size <= cache_bytes) {
          cache$state[[key]] <- state
          cache$bytes <- cache$bytes + size
        }
        states[[k]] <- state
      }
      if (!is.null(cache$state[[key]])) {
        cache$clock <- cache$clock + 1
        cache$last[key] <- cache$clock
      }
    }
    failed <- vapply(states, function(z) !is.null(z$error), logical(1L))
    summaries <- list(
      used = active, a = lapply(states, `[[`, "a"),
      H = lapply(states, `[[`, "M"), has_summary = !failed,
      error_message = vapply(states, function(z) {
        if (is.null(z$error)) NA_character_ else z$error
      }, character(1L)), elapsed = 0
    )
    z <- .mgcvst_liu_pairs(index[rows, , drop = FALSE], pair_index[rows],
                           fit$feature_id, summaries, threads, length(rows), FALSE)
    r_score[rows] <- z$result$score
    r_information[rows] <- z$result$information
    r_effective_rank[rows] <- z$result$effective_rank
    r_p_value[rows] <- z$result$p_value
    r_log_p_two_sided[rows] <- z$result$log_p_two_sided
    r_log_p_positive[rows] <- z$result$log_p_positive
    r_log_p_negative[rows] <- z$result$log_p_negative
    r_error_message[rows] <- z$result$error_message
    .mgcvst_pair_checkpoint_write(pair_path, first, last, z$result)
    elapsed <- elapsed + z$elapsed
    chunks <- chunks + 1L
    if (verbose && (chunks %% 10L == 0L || last == nrow(index))) {
      message("Evaluated Liu pair group ", chunks, ".")
    }
    rm(states, summaries, z)
    if (exists("state", inherits = FALSE)) rm(state)
    first <- last + 1L
  }
  back <- order(ord)
  result <- data.frame(
    pair_index = pair_index[back], score = r_score[back],
    information = r_information[back], effective_rank = r_effective_rank[back],
    p_value = r_p_value[back], log_p_two_sided = r_log_p_two_sided[back],
    log_p_positive = r_log_p_positive[back], log_p_negative = r_log_p_negative[back],
    error_message = r_error_message[back], stringsAsFactors = FALSE
  )
  rownames(result) <- NULL
  list(result = result, elapsed = elapsed,
       metadata = list(path = if (isTRUE(store$temporary)) NULL else store$path,
                       builds = builds, resume_count = resume_count,
                       resumed_pairs = resumed_pairs, pair_path = pair_path,
                       preparation_backend = mode,
                       pair_schedule = "resident_gene_blocks",
                       chunks = chunks, cache_hits = cache$hits,
                       cache_misses = cache$misses,
                       cache_evictions = cache$evictions,
                       cache_bytes = cache_bytes,
                       resident_bytes = cache$bytes,
                       preparation_elapsed = preparation_elapsed,
                       signature = signature, storage = "double"))
}
