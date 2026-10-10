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
         basis = basis[c("coordinate", "basis", "rank", "kind")],
         basis_spec = fit$basis_spec,
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

# Build one bounded feature batch with the full-precision dense score kernel.
# A model.set() fit without a usable conditional nuisance covariance is a hard
# error at the call site of .mgcvst_pair_pipeline(), never a per-feature R-loop
# fallback.
.mgcvst_pair_build_batch <- function(fit, ids, threads, native) {
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
  lapply(seq_along(ids), function(k) {
    if (bad[k]) return(list(error =
      "The feature has invalid dispersion or smoothing parameters."))
    if (!is.null(z[[k]]$error)) return(list(error = z[[k]]$error))
    list(a = z[[k]]$a, M = z[[k]]$H, width = native$width)
  })
}

# Scores and Liu log p-values of the pairs `index` (global feature indices)
# from the score states of the features `active`, one state per feature. A
# pair with a feature whose state failed is returned with status 3.
.mgcvst_liu_pairs <- function(index, active, states, threads) {
  out <- .mgcvst_pairs_frame(index[, 1L], index[, 2L],
                             status = .mgcvst_pair_status[["feature"]])
  good <- vapply(states, function(z) is.null(z$error), logical(1L))
  local <- matrix(match(index, active), ncol = 2L)
  rows <- which(good[local[, 1L]] & good[local[, 2L]])
  if (!length(rows)) return(out)
  keep <- which(good)
  avec <- do.call(cbind, lapply(states[keep], `[[`, "a"))
  H <- lapply(states[keep], `[[`, "M")
  pair <- matrix(match(local[rows, , drop = FALSE], keep), ncol = 2L)
  perm <- order(pair[, 1L])
  res <- mgcvst_pair_liu_cpp(H, avec, pair[perm, 1L], pair[perm, 2L], threads)
  at <- rows[perm]
  out$score[at] <- res$score
  out$log_p_two_sided[at] <- res$log_p_two_sided
  out$log_p_positive[at] <- res$log_p_positive
  out$log_p_negative[at] <- res$log_p_negative
  out$status[at] <- res$status
  # A pair that is not evaluated carries no p-value, as in the PCAlearning
  # route; a non-finite log p-value (-Inf) would otherwise enter the
  # adjustment as p = 0 and count as a discovery.
  invalid <- which(out$status != .mgcvst_pair_status[["ok"]])
  if (length(invalid)) {
    out$log_p_two_sided[invalid] <- NA_real_
    out$log_p_positive[invalid] <- NA_real_
    out$log_p_negative[invalid] <- NA_real_
  }
  out
}

# Evaluate the Liu-calibrated pairs of a model.set() fit from resumable
# feature-first score states and stream them to raw Parquet shards. `index`
# is NULL for every pair of the available features, or a two-column matrix
# of available feature indices i < j. The resident score states adapt to the
# memory available to the process.
.mgcvst_pair_pipeline <- function(fit, index, threads, chunk_size, verbose,
                                  checkpoint_dir = NULL, resume = TRUE) {
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
  if (identical(fit$score_backend, "sparse")) {
    stop("Sparse INLA fits are tested by inlaST.test().")
  }
  available <- .mgcvst_feature_available(fit)
  all_pairs <- is.null(index)
  if (all_pairs) {
    used <- which(available)
    if (length(used) < 2L) stop("At least two available features are required.")
  } else {
    if (!is.matrix(index) || ncol(index) != 2L || !nrow(index) ||
        anyNA(index) || any(index != floor(index)) ||
        any(index < 1L) || any(index > length(fit$feature_id)) ||
        any(index[, 1L] >= index[, 2L]) || !all(available[index])) {
      stop("index must contain two-column indices i < j of available features.")
    }
    used <- sort(unique(as.vector(index)))
  }
  target <- fit$geometry$target
  width <- if (length(target)) {
    ncol(fit$geometry$smooth[[unname(target[[1L]])]]$B)
  } else 1L
  if (!is.finite(width) || width < 1L) stop("The score coordinate width is invalid.")

  contract <- .mgcvst_contract("exact")
  signature <- if (is.null(checkpoint_dir)) {
    list(version = 2L, temporary = tempfile("mgcvst-pair-run-"))
  } else .mgcvst_pair_signature(fit)
  store <- .mgcvst_store_open(checkpoint_dir, signature, fit$feature_id,
                              resume = resume)
  if (isTRUE(store$temporary)) on.exit(.mgcvst_store_cleanup(store), add = TRUE)
  root <- if (is.null(checkpoint_dir)) {
    tempfile("mgcvst-pairs-")
  } else store$path
  # Pair directories of another algorithm contract are refused before any
  # work; the directory of this run is opened once its score states exist.
  .mgcvst_pairs_refuse_stale(root, contract)
  universe <- if (all_pairs) {
    list(all = TRUE, used = used, n_feature = length(fit$feature_id))
  } else list(index = index)

  n_pairs_hint <- if (all_pairs) length(used) * (length(used) - 1) / 2 else
    nrow(index)
  state_estimate <- 8 * (width^2 + width) + 2048
  probe <- .mgcvst_memory_probe()
  reserve <- 4 * 8 * width^2 * min(threads, chunk_size) +
    2 * state_estimate + 256 * min(n_pairs_hint, chunk_size) + 64 * 1024^2
  cache_bytes <- if (is.finite(probe$available)) {
    max(0, 0.7 * probe$available - reserve)
  } else 512 * 1024^2

  existing <- vapply(used, function(id) .mgcvst_store_has(store, id), logical(1L))
  resume_count <- sum(existing)
  missing <- used[!existing]
  builds <- 0L
  preparation_started <- proc.time()[["elapsed"]]
  if (length(missing)) {
    fit$.mgcvst_fixed_factors <- .mgcvst_model_fixed_factors(fit)
    native <- .mgcvst_model_dense_preparation(fit, missing)
    if (is.null(native)) {
      stop("Model score states require the conditional nuisance covariance; ",
           "re-estimate with the current mgcvST.estimate().")
    }
    q <- width
    p <- ncol(native$X)
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
      states <- .mgcvst_pair_build_batch(fit, ids, threads, native)
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
    rm(native)
    fit$.mgcvst_fixed_factors <- NULL
  }
  preparation_elapsed <- proc.time()[["elapsed"]] - preparation_started
  pair_dir <- .mgcvst_pairs_open(root, universe, contract, resume)

  # Resident score states: least-recently-used eviction under `cache_bytes`.
  cache <- new.env(parent = emptyenv())
  cache$state <- list()
  cache$bytes <- 0
  cache$last <- numeric()
  cache$clock <- 0
  cache$hits <- 0L
  cache$misses <- 0L
  cache$evictions <- 0L
  failed <- character()
  drop_oldest <- function(protect = character()) {
    candidate <- setdiff(names(cache$last), protect)
    if (!length(candidate)) return(FALSE)
    drop <- candidate[which.min(cache$last[candidate])]
    cache$bytes <- cache$bytes - as.numeric(object.size(cache$state[[drop]]))
    cache$state[[drop]] <- NULL
    cache$last <- cache$last[names(cache$last) != drop]
    cache$evictions <- cache$evictions + 1L
    TRUE
  }
  fetch <- function(active) {
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
        state_estimate <<- max(state_estimate, size)
        while (cache$bytes + size > cache_bytes && length(cache$last)) {
          if (!drop_oldest(as.character(active))) break
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
      if (!is.null(states[[k]]$error)) failed[key] <<- states[[k]]$error
    }
    states
  }

  shard_files <- character()
  shard_rows <- integer()
  elapsed <- 0
  chunks <- 0L
  resumed_pairs <- 0
  evaluate <- function(window, active, id) {
    t0 <- proc.time()[["elapsed"]]
    states <- fetch(active)
    frame <- .mgcvst_liu_pairs(window, active, states, threads)
    .mgcvst_write_parquet(frame, .mgcvst_shard_file(pair_dir, id))
    elapsed <<- elapsed + proc.time()[["elapsed"]] - t0
    nrow(frame)
  }

  if (all_pairs) {
    capacity <- max(2L, min(length(used), floor(cache_bytes / state_estimate)))
    schedule <- .mgcvst_all_pair_schedule(pair_dir, length(used), capacity,
                                          chunk_size)
    blocks <- schedule$blocks
    for (b in seq_len(nrow(blocks))) {
      pos <- .mgcvst_block_pairs(blocks[b, ])
      window <- cbind(used[pos[, 1L]], used[pos[, 2L]])
      if (.mgcvst_shard_complete(pair_dir, b, window[, 1L], window[, 2L])) {
        resumed_pairs <- resumed_pairs + nrow(window)
      } else {
        evaluate(window, sort(unique(as.vector(window))), b)
      }
      shard_files <- c(shard_files, .mgcvst_shard_file(pair_dir, b))
      shard_rows <- c(shard_rows, nrow(window))
      chunks <- chunks + 1L
      if (verbose && (chunks %% 10L == 0L || b == nrow(blocks))) {
        message("Evaluated Liu pair block ", b, " of ", nrow(blocks), ".")
      }
    }
    n_pairs <- sum(shard_rows)
    pair_schedule <- "resident_gene_tiles"
  } else {
    ord <- .mgcvst_pair_order(index, used,
      max(2L, floor(cache_bytes / state_estimate)), pair_dir)
    index <- index[ord, , drop = FALSE]
    n_pairs <- nrow(index)
    first <- 1L
    last_probe <- proc.time()[["elapsed"]]
    while (first <= n_pairs) {
      saved <- .mgcvst_shard_window(pair_dir, first, index)
      if (!is.na(saved)) {
        resumed_pairs <- resumed_pairs + saved
        shard_files <- c(shard_files, .mgcvst_shard_file(pair_dir, first))
        shard_rows <- c(shard_rows, saved)
        chunks <- chunks + 1L
        first <- first + saved
        next
      }
      if (proc.time()[["elapsed"]] - last_probe >= 2) {
        probe <- .mgcvst_memory_probe()
        if (is.finite(probe$available)) {
          cache_bytes <- max(0, 0.7 * (probe$available + cache$bytes) - reserve)
        }
        last_probe <- proc.time()[["elapsed"]]
      }
      max_features <- max(2L, min(length(used),
        floor(cache_bytes / state_estimate)))
      while (cache$bytes > cache_bytes && drop_oldest()) NULL
      cap <- as.integer(min(n_pairs, first + as.double(chunk_size) - 1))
      endpoints <- as.vector(t(index[seq.int(first, cap), , drop = FALSE]))
      first_seen <- which(!duplicated(endpoints))
      last <- if (length(first_seen) > max_features) {
        first + (first_seen[max_features + 1L] - 1L) %/% 2L - 1L
      } else cap
      window <- index[seq.int(first, last), , drop = FALSE]
      evaluate(window, sort(unique(as.vector(window))), first)
      shard_files <- c(shard_files, .mgcvst_shard_file(pair_dir, first))
      shard_rows <- c(shard_rows, nrow(window))
      chunks <- chunks + 1L
      if (verbose && (chunks %% 10L == 0L || last == n_pairs)) {
        message("Evaluated Liu pair group ", chunks, ".")
      }
      first <- last + 1L
    }
    pair_schedule <- "resident_gene_blocks"
  }

  list(
    pair_dir = pair_dir, shards = shard_files, rows = shard_rows,
    n_pairs = n_pairs, temporary = is.null(checkpoint_dir),
    failed = data.frame(
      feature_id = fit$feature_id[as.integer(names(failed))],
      error = unname(failed), stringsAsFactors = FALSE
    ),
    elapsed = elapsed,
    metadata = list(
      path = if (isTRUE(store$temporary)) NULL else store$path,
      pair_dir = pair_dir, builds = builds, resume_count = resume_count,
      resumed_pairs = resumed_pairs, preparation_backend = "model_native",
      pair_schedule = pair_schedule, chunks = chunks,
      cache_hits = cache$hits, cache_misses = cache$misses,
      cache_evictions = cache$evictions, cache_bytes = cache_bytes,
      resident_bytes = cache$bytes, preparation_elapsed = preparation_elapsed,
      signature = signature, contract = contract
    )
  )
}
