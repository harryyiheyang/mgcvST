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

# Defaults of the exact-moment route: the number of leading singular values
# kept from the shared basis.
.mgcvst_exact_defaults <- list(k = 20L)

# Scores and saddlepoint log p-values of the pairs `index` (global feature
# indices) from the score states of the features `active`, one state per
# feature, and their pair bases G (indexed by global feature). `degenerate`
# flags (by global feature) the genes with a degenerate spatial fit; they have
# no state here.
#   * a pair with a feature whose state failed has status 3 and no p-value;
#   * otherwise a pair with a degenerate gene has status 4 and p = 1, and no
#     score;
#   * a pair the kernel could not evaluate (status 1 or 2) keeps its status and
#     has p = 1: two-sided and both one-sided log p are 0.
# The attribute `nodes_above_leading` counts the pairs whose remainder has a
# node above the largest leading value.
.mgcvst_spa_pairs <- function(index, active, states, G, threads, order = 4L,
                              degenerate = NULL) {
  out <- .mgcvst_pairs_frame(index[, 1L], index[, 2L],
                             status = .mgcvst_pair_status[["feature"]])
  good <- vapply(states, function(z) is.null(z$error), logical(1L))
  local <- matrix(match(index, active), ncol = 2L)
  deg <- if (is.null(degenerate)) rep(FALSE, nrow(index)) else
    degenerate[index[, 1L]] | degenerate[index[, 2L]]
  # A gene that has no state because it is degenerate is not a failure.
  state_failed <- (!is.na(local[, 1L]) & !good[local[, 1L]]) |
    (!is.na(local[, 2L]) & !good[local[, 2L]])
  rows <- which(!deg & !is.na(local[, 1L]) & !is.na(local[, 2L]))
  rows <- rows[good[local[rows, 1L]] & good[local[rows, 2L]]]
  attr(out, "nodes_above_leading") <- 0
  set_degenerate <- function() {
    at <- which(deg & !state_failed)
    out$status[at] <<- .mgcvst_pair_status[["degenerate"]]
    out$log_p_two_sided[at] <<- 0
    out$log_p_positive[at] <<- 0
    out$log_p_negative[at] <<- 0
    out
  }
  if (!length(rows)) return(set_degenerate())
  keep <- which(good)
  avec <- do.call(cbind, lapply(states[keep], `[[`, "a"))
  H <- lapply(states[keep], `[[`, "M")
  pair <- matrix(match(local[rows, , drop = FALSE], keep), ncol = 2L)
  perm <- order(pair[, 1L])
  res <- mgcvst_pair_spa_cpp(H, G[active[keep]], avec, pair[perm, 1L],
                             pair[perm, 2L], threads, order)
  at <- rows[perm]
  out$score[at] <- res$score
  out$log_p_two_sided[at] <- res$log_p_two_sided
  out$log_p_positive[at] <- res$log_p_positive
  out$log_p_negative[at] <- res$log_p_negative
  out$remainder_kind[at] <- res$remainder_kind
  out$status[at] <- res$status
  if (!is.null(res$nodes_above_leading)) {
    attr(out, "nodes_above_leading") <- res$nodes_above_leading
  }
  # A pair the kernel could not evaluate has p = 1 and keeps its status code,
  # so that a non-finite log p-value never enters the adjustment and the pair
  # stays in its family.
  invalid <- which(out$status %in% c(.mgcvst_pair_status[["moments"]],
                                     .mgcvst_pair_status[["p_value"]]))
  if (length(invalid)) {
    out$log_p_two_sided[invalid] <- 0
    out$log_p_positive[invalid] <- 0
    out$log_p_negative[invalid] <- 0
  }
  set_degenerate()
}

# How the dense score state (a, M) of a feature is produced: from the working
# model of an mgcv fit, or by reconstructing the reduced curvature of a sparse
# INLA fit in its observation basis. `missing` are the features to build.
.mgcvst_state_builder <- function(fit, missing, threads, basis = NULL) {
  if (identical(fit$score_backend, "sparse")) {
    fit <- .inlast_sparse_prepare(fit)
    r <- basis$rank
    m <- ncol(fit$score_sparse$Q)
    return(list(
      work = 8 * (6 * r^2 + 4 * m * r),
      build = function(ids) {
        units <- .inlast_sparse_units(fit, ids, threads = threads)
        bad <- vapply(units, function(z) !is.null(z$error), logical(1L))
        states <- vector("list", length(ids))
        if (any(!bad)) {
          z <- .inlast_sparse_materialize_reduced(fit, units[!bad], basis,
                                                  threads = threads)
          states[!bad] <- lapply(z, function(s) {
            if (!is.null(s$error)) list(error = as.character(s$error)) else
              list(a = as.numeric(s$a), M = s$M, width = r)
          })
        }
        for (j in which(bad)) states[[j]] <- list(error = as.character(units[[j]]$error))
        states
      },
      finish = function() invisible(NULL)
    ))
  }
  fit$.mgcvst_fixed_factors <- .mgcvst_model_fixed_factors(fit)
  native <- .mgcvst_model_dense_preparation(fit, missing)
  if (is.null(native)) {
    stop("Model score states require the conditional nuisance covariance; ",
         "re-estimate with the current mgcvST.estimate().")
  }
  q <- .mgcvst_state_width(fit)
  p <- ncol(native$X)
  if (is.null(p)) p <- 0L
  n <- nrow(fit$working_variance)
  list(
    work = 8 * (8 * n * q + 4 * n * p + 8 * q^2 + 4 * q * p + 4 * p^2),
    build = function(ids) .mgcvst_pair_build_batch(fit, ids, threads, native),
    finish = function() invisible(NULL)
  )
}

# Width q of the score coordinates of the states of a fit.
.mgcvst_state_width <- function(fit, basis = NULL) {
  if (identical(fit$score_backend, "sparse")) return(as.integer(basis$rank))
  target <- fit$geometry$target
  width <- if (length(target)) {
    ncol(fit$geometry$smooth[[unname(target[[1L]])]]$B)
  } else 1L
  if (!is.finite(width) || width < 1L) stop("The score coordinate width is invalid.")
  as.integer(width)
}

# Shared basis V of the exact route: the k leading eigenvectors of the sum of
# the normalized state matrices over `used`, added in feature order whatever the
# batching (so the basis, and its sha, are reproducible). `sum` is the running
# sum when the states were accumulated while they were built, or NULL, in which
# case the stored states are read in one pass. A saved basis for the same
# feature set is reused.
.mgcvst_exact_basis <- function(store, used, width, k, path, sum, batch_size,
                                verbose) {
  key <- digest::digest(list(used = used, k = k, width = width), algo = "sha256")
  file <- if (is.null(path)) NULL else file.path(path, "basis.rds")
  if (!is.null(file) && file.exists(file)) {
    saved <- tryCatch(readRDS(file), error = function(e) NULL)
    if (is.list(saved) && identical(saved$key, key) &&
        identical(saved$sha, digest::digest(saved$V, algo = "sha256"))) {
      if (verbose) message("Reused the shared basis (k = ", k, ").")
      return(saved)
    }
  }
  S <- sum
  if (is.null(S)) {
    S <- matrix(0, width, width)
    for (first in seq.int(1L, length(used), by = batch_size)) {
      ids <- used[first:min(length(used), first + batch_size - 1L)]
      states <- lapply(ids, function(id) .mgcvst_store_read(store, id))
      good <- vapply(states, function(z) is.null(z$error), logical(1L))
      if (any(good)) {
        S <- mgcvst_pair_basis_sum_cpp(lapply(states[good], `[[`, "M"), S)
      }
    }
  }
  if (!all(is.finite(S)) || !any(S != 0)) {
    stop("No usable score state to build the shared basis.")
  }
  eg <- eigen(S, symmetric = TRUE)
  V <- eg$vectors[, seq_len(k), drop = FALSE]
  basis <- list(V = V, k = k, q = width, values = eg$values[seq_len(k)],
                n_features = length(used), key = key,
                sha = digest::digest(V, algo = "sha256"), version = 1L)
  if (!is.null(file)) {
    tmp <- tempfile("basis-", tmpdir = path, fileext = ".tmp")
    on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
    saveRDS(basis, tmp, compress = FALSE)
    if (file.exists(file)) unlink(file)
    if (!file.rename(tmp, file)) stop("Could not commit the shared basis.")
  }
  if (verbose) message("Built the shared basis (k = ", k, ") from ", length(used),
                       " features.")
  basis
}

# Evaluate the saddlepoint-calibrated pairs of a fit from resumable
# feature-first score states and stream them to raw Parquet shards. `index`
# is NULL for every pair of the available features, or a two-column matrix
# of available feature indices i < j. The resident score states adapt to the
# memory available to the process. The states are those of the working model
# for an mgcv fit and the reduced curvature in `basis` for a sparse INLA fit.
.mgcvst_pair_pipeline <- function(fit, index, threads, chunk_size, verbose,
                                  checkpoint_dir = NULL, resume = TRUE,
                                  k = .mgcvst_exact_defaults$k, basis = NULL,
                                  route = NULL) {
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
  sparse <- identical(fit$score_backend, "sparse")
  if (sparse && is.null(basis)) stop("A sparse INLA fit needs its observation basis.")
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
  width <- .mgcvst_state_width(fit, basis)
  k <- as.integer(min(k, width))
  if (!is.finite(k) || k < 1L) stop("k must be a positive integer.")
  # Genes with a degenerate spatial fit have no state and no pair basis: their
  # pairs are written with p = 1 and status 4.
  degenerate <- .mgcvst_degenerate_features(fit)
  used_kernel <- used[!degenerate[used]]

  contract_early <- .mgcvst_contract("exact")
  signature <- if (is.null(checkpoint_dir)) {
    list(version = 2L, temporary = tempfile("mgcvst-pair-run-"))
  } else .mgcvst_pair_signature(fit, basis)
  store <- .mgcvst_store_open(checkpoint_dir, signature, fit$feature_id,
                              resume = resume)
  if (isTRUE(store$temporary)) on.exit(.mgcvst_store_cleanup(store), add = TRUE)
  if (!isTRUE(store$temporary)) .mgcvst_route_save(store$path, route)
  root <- if (is.null(checkpoint_dir)) {
    tempfile("mgcvst-pairs-")
  } else store$path
  # Pair directories of another algorithm contract are refused before any
  # work; the directory of this run is opened once the shared basis exists.
  .mgcvst_pairs_refuse_stale(root, contract_early)
  universe <- if (all_pairs) {
    list(all = TRUE, used = used, n_feature = length(fit$feature_id),
         degenerate = which(degenerate[used]))
  } else list(index = index, degenerate = which(degenerate[used]))

  n_pairs_hint <- if (all_pairs) length(used) * (length(used) - 1) / 2 else
    nrow(index)
  state_estimate <- 8 * (width^2 + width) + 2048
  basis_bytes <- 8 * width * k * length(used_kernel)
  probe <- .mgcvst_memory_probe()
  if (is.finite(probe$available) && basis_bytes > 0.3 * probe$available) {
    stop(sprintf(paste0("The pair bases of %d features need about %.1f GB, ",
                        "above 30%% of the %.1f GB available; use ",
                        "moments = \"pcalearning\", fewer features or a smaller k."),
                 length(used_kernel), basis_bytes / 1024^3, probe$available / 1024^3))
  }
  reserve <- 4 * 8 * width^2 * min(threads, chunk_size) +
    2 * state_estimate + 256 * min(n_pairs_hint, chunk_size) + 64 * 1024^2 +
    basis_bytes
  cache_bytes <- if (is.finite(probe$available)) {
    max(0, 0.7 * probe$available - reserve)
  } else 512 * 1024^2

  existing <- vapply(used_kernel, function(id) .mgcvst_store_has(store, id),
                     logical(1L))
  resume_count <- sum(existing)
  missing <- used_kernel[!existing]
  builds <- 0L
  preparation_started <- proc.time()[["elapsed"]]
  # When every state is built in this run, in feature order, the sum behind the
  # shared basis is accumulated as they are built.
  running <- if (identical(missing, used_kernel)) matrix(0, width, width) else NULL
  if (length(missing)) {
    builder <- .mgcvst_state_builder(fit, missing, threads, basis)
    first <- 1L
    while (first <= length(missing)) {
      probe <- .mgcvst_memory_probe()
      headroom <- if (!is.null(probe) && is.finite(probe$available)) {
        0.3 * probe$available
      } else max(cache_bytes, 512 * 1024^2)
      if (headroom < builder$work) {
        stop("Insufficient available memory for one score-state preparation. ",
             "Increase the job memory allocation.")
      }
      batch_size <- as.integer(max(1L, min(32L, floor(headroom / builder$work))))
      ids <- missing[first:min(length(missing), first + batch_size - 1L)]
      states <- builder$build(ids)
      if (length(states) != length(ids)) {
        stop("The score-state backend returned the wrong feature count.")
      }
      for (j in seq_along(ids)) .mgcvst_store_write(store, ids[j], states[[j]])
      if (!is.null(running)) {
        good <- vapply(states, function(z) is.null(z$error), logical(1L))
        if (any(good)) {
          running <- mgcvst_pair_basis_sum_cpp(lapply(states[good], `[[`, "M"),
                                               running)
        }
      }
      builds <- builds + length(ids)
      rm(states)
      if (verbose) message("Stored score states for ", builds, " of ",
                           length(missing), " remaining features.")
      first <- first + length(ids)
    }
    builder$finish()
    rm(builder)
  }
  read_batch <- 32L
  G <- vector("list", length(fit$feature_id))
  if (length(used_kernel)) {
    shared <- .mgcvst_exact_basis(store, used_kernel, width, k,
      path = if (isTRUE(store$temporary)) NULL else store$path,
      sum = running, batch_size = read_batch, verbose = verbose)
    # Pair bases G_g = H_g^(1/2) V of every used feature, resident.
    for (first in seq.int(1L, length(used_kernel), by = read_batch)) {
      ids <- used_kernel[first:min(length(used_kernel), first + read_batch - 1L)]
      states <- lapply(ids, function(id) .mgcvst_store_read(store, id))
      good <- vapply(states, function(z) is.null(z$error), logical(1L))
      if (any(good)) {
        G[ids[good]] <- mgcvst_pair_basis_cpp(lapply(states[good], `[[`, "M"),
                                              shared$V, threads)
      }
    }
  } else {
    # Every gene of the pairs is degenerate: no basis, every pair has p = 1.
    shared <- list(V = NULL, sha = digest::digest(list(no_kernel_gene = TRUE),
                                                 algo = "sha256"))
  }
  rm(running)
  preparation_elapsed <- proc.time()[["elapsed"]] - preparation_started
  contract <- .mgcvst_contract("exact", k = k, basis_sha = shared$sha)
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
    for (j in seq_along(active)) {
      key <- as.character(active[j])
      if (!is.null(cache$state[[key]])) {
        cache$hits <- cache$hits + 1L
        states[[j]] <- cache$state[[key]]
      } else {
        cache$misses <- cache$misses + 1L
        state <- .mgcvst_store_read(store, active[j])
        size <- as.numeric(object.size(state))
        state_estimate <<- max(state_estimate, size)
        while (cache$bytes + size > cache_bytes && length(cache$last)) {
          if (!drop_oldest(as.character(active))) break
        }
        if (cache$bytes + size <= cache_bytes) {
          cache$state[[key]] <- state
          cache$bytes <- cache$bytes + size
        }
        states[[j]] <- state
      }
      if (!is.null(cache$state[[key]])) {
        cache$clock <- cache$clock + 1
        cache$last[key] <- cache$clock
      }
      if (!is.null(states[[j]]$error)) failed[key] <<- states[[j]]$error
    }
    states
  }

  shard_files <- character()
  shard_rows <- integer()
  elapsed <- 0
  chunks <- 0L
  resumed_pairs <- 0
  nodes_above <- 0
  evaluate <- function(window, active, id) {
    t0 <- proc.time()[["elapsed"]]
    active <- setdiff(active, which(degenerate))
    states <- fetch(active)
    frame <- .mgcvst_spa_pairs(window, active, states, G, threads, order = 4L,
                               degenerate = degenerate)
    nodes_above <<- nodes_above + as.numeric(attr(frame, "nodes_above_leading"))
    attr(frame, "nodes_above_leading") <- NULL
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
        message("Evaluated pair block ", b, " of ", nrow(blocks), ".")
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
        message("Evaluated pair group ", chunks, ".")
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
      resumed_pairs = resumed_pairs,
      preparation_backend = if (sparse) "sparse_reduced" else "model_native",
      pair_schedule = pair_schedule, chunks = chunks, k = k,
      basis_sha = shared$sha, q = width, nodes_above_leading = nodes_above,
      degenerate_genes = sum(degenerate[used]),
      cache_hits = cache$hits, cache_misses = cache$misses,
      cache_evictions = cache$evictions, cache_bytes = cache_bytes,
      resident_bytes = cache$bytes, preparation_elapsed = preparation_elapsed,
      signature = signature, contract = contract
    )
  )
}
