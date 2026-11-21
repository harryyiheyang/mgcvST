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

# Defaults of the exact-moment route: the number of leading singular values
# kept from the shared basis.
.mgcvst_exact_defaults <- list(k = 20L)

# Scores and saddlepoint log p-values of the pairs `index` (global feature
# indices) from the score states of the features `active`, one state per
# feature, and their pair bases G (indexed by global feature). A pair with a
# feature whose state failed is returned with status 3 and no p-value; a pair
# the kernel could not evaluate (status 1 or 2) keeps its status and has p = 1:
# two-sided and both one-sided log p are 0. The attribute `nodes_above_leading`
# counts the pairs whose remainder has a node above the largest leading value.
.mgcvst_spa_pairs <- function(index, active, states, G, threads, order = 4L,
                              screen_z = 0) {
  out <- .mgcvst_pairs_frame(index[, 1L], index[, 2L],
                             status = .mgcvst_pair_status[["feature"]])
  good <- vapply(states, function(z) is.null(z$error), logical(1L))
  local <- matrix(match(index, active), ncol = 2L)
  rows <- which(good[local[, 1L]] & good[local[, 2L]])
  attr(out, "nodes_above_leading") <- 0
  attr(out, "screened") <- 0
  if (!length(rows)) return(out)
  keep <- which(good)
  avec <- do.call(cbind, lapply(states[keep], `[[`, "a"))
  H <- lapply(states[keep], `[[`, "M")
  pair <- matrix(match(local[rows, , drop = FALSE], keep), ncol = 2L)
  perm <- order(pair[, 1L])
  res <- if (screen_z > 0) {
    mgcvst_pair_spa_cpp(H, G[active[keep]], avec, pair[perm, 1L],
                        pair[perm, 2L], threads, order, NULL, screen_z)
  } else {
    mgcvst_pair_spa_cpp(H, G[active[keep]], avec, pair[perm, 1L],
                        pair[perm, 2L], threads, order)
  }
  if (!is.null(res$screened)) attr(out, "screened") <- res$screened
  at <- rows[perm]
  out$score[at] <- res$score
  out$t1[at] <- res$t1
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
  out
}

# The score state (a, M) of a feature: the reduced curvature of its fit in the
# observation basis. `missing` are the features to build.
.mgcvst_state_builder <- function(fit, missing, threads, basis) {
  fit <- .spdest_sparse_prepare(fit)
  r <- basis$rank
  m <- ncol(fit$score_sparse$Q)
  list(
    work = 8 * (6 * r^2 + 4 * m * r),
    build = function(ids) {
      units <- .spdest_sparse_units(fit, ids, threads = threads)
      bad <- vapply(units, function(z) !is.null(z$error), logical(1L))
      states <- vector("list", length(ids))
      if (any(!bad)) {
        z <- .spdest_sparse_materialize_reduced(fit, units[!bad], basis,
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
  )
}

# Shared basis V of the exact route: the k leading eigenvectors of the sum
# `S` of the normalized state matrices, added in feature order whatever the
# batching (so the basis, and its sha, are reproducible).
.mgcvst_exact_basis <- function(S, n_features, k, verbose) {
  if (!all(is.finite(S)) || !any(S != 0)) {
    stop("No usable score state to build the shared basis.")
  }
  eg <- eigen(S, symmetric = TRUE)
  V <- eg$vectors[, seq_len(k), drop = FALSE]
  if (verbose) message("Built the shared basis (k = ", k, ") from ", n_features,
                       " features.")
  list(V = V, k = k, q = nrow(S), values = eg$values[seq_len(k)],
       n_features = n_features, sha = digest::digest(V, algo = "sha256"))
}

# Evaluate the saddlepoint-calibrated pairs of a fit from feature-first score
# states, kept in a temporary directory, and stream them to raw Parquet shards. `index`
# is NULL for every pair of the available features, or a two-column matrix
# of available feature indices i < j. The resident score states adapt to the
# memory available to the process. The states are the reduced curvature of
# the fit in its observation `basis`.
.mgcvst_pair_pipeline <- function(fit, index, threads, chunk_size, verbose,
                                  k = .mgcvst_exact_defaults$k, basis,
                                  screen_z = 0) {
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
  width <- as.integer(basis$rank)
  k <- as.integer(min(k, width))
  if (!is.finite(k) || k < 1L) stop("k must be a positive integer.")

  store <- .mgcvst_store_open(fit$feature_id)
  on.exit(.mgcvst_store_cleanup(store), add = TRUE)

  n_pairs_hint <- if (all_pairs) length(used) * (length(used) - 1) / 2 else
    nrow(index)
  state_estimate <- 8 * (width^2 + width) + 2048
  basis_bytes <- 8 * width * k * length(used)
  probe <- .mgcvst_memory_probe()
  if (is.finite(probe$available) && basis_bytes > 0.3 * probe$available) {
    stop(sprintf(paste0("The pair bases of %d features need about %.1f GB, ",
                        "above 30%% of the %.1f GB available; use ",
                        "moments = \"pcalearning\", fewer features or a smaller k."),
                 length(used), basis_bytes / 1024^3, probe$available / 1024^3))
  }
  reserve <- 4 * 8 * width^2 * min(threads, chunk_size) +
    2 * state_estimate + 256 * min(n_pairs_hint, chunk_size) + 64 * 1024^2 +
    basis_bytes
  cache_bytes <- if (is.finite(probe$available)) {
    max(0, 0.7 * probe$available - reserve)
  } else 512 * 1024^2

  missing <- used
  builds <- 0L
  preparation_started <- proc.time()[["elapsed"]]
  # The sum behind the shared basis is accumulated as the states are built, in
  # feature order.
  running <- matrix(0, width, width)
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
      good <- vapply(states, function(z) is.null(z$error), logical(1L))
      if (any(good)) {
        running <- mgcvst_pair_basis_sum_cpp(lapply(states[good], `[[`, "M"),
                                             running)
      }
      builds <- builds + length(ids)
      rm(states)
      if (verbose) message("Stored score states for ", builds, " of ",
                           length(missing), " features.")
      first <- first + length(ids)
    }
    builder$finish()
    rm(builder)
  }
  read_batch <- 32L
  shared <- .mgcvst_exact_basis(running, length(used), k, verbose)
  rm(running)
  # Pair bases G_g = H_g^(1/2) V of every used feature, resident, and the
  # traces tr(M_g) and tr(M_g^2) of each state, kept with the test.
  G <- vector("list", length(fit$feature_id))
  trace_M <- trace_M2 <- rep(NA_real_, length(used))
  for (first in seq.int(1L, length(used), by = read_batch)) {
    at <- first:min(length(used), first + read_batch - 1L)
    ids <- used[at]
    states <- lapply(ids, function(id) .mgcvst_store_read(store, id))
    good <- vapply(states, function(z) is.null(z$error), logical(1L))
    if (any(good)) {
      G[ids[good]] <- mgcvst_pair_basis_cpp(lapply(states[good], `[[`, "M"),
                                            shared$V, threads)
      for (g in which(good)) {
        M <- states[[g]]$M
        trace_M[at[g]] <- sum(diag(M))
        trace_M2[at[g]] <- sum(M * t(M))
      }
    }
  }
  preparation_elapsed <- proc.time()[["elapsed"]] - preparation_started
  contract <- .mgcvst_contract("exact", k = k, basis_sha = shared$sha,
                               screen_z = screen_z)
  pair_dir <- .mgcvst_pairs_open()
  finished <- FALSE
  on.exit(if (!finished) unlink(pair_dir, recursive = TRUE), add = TRUE)

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
  nodes_above <- 0
  screened <- 0
  evaluate <- function(window, active, id) {
    t0 <- proc.time()[["elapsed"]]
    states <- fetch(active)
    frame <- if (screen_z > 0) {
      .mgcvst_spa_pairs(window, active, states, G, threads, order = 4L,
                        screen_z = screen_z)
    } else .mgcvst_spa_pairs(window, active, states, G, threads, order = 4L)
    nodes_above <<- nodes_above + as.numeric(attr(frame, "nodes_above_leading"))
    if (!is.null(attr(frame, "screened"))) {
      screened <<- screened + as.numeric(attr(frame, "screened"))
    }
    attr(frame, "nodes_above_leading") <- NULL
    attr(frame, "screened") <- NULL
    .mgcvst_write_parquet(frame, .mgcvst_shard_file(pair_dir, id))
    elapsed <<- elapsed + proc.time()[["elapsed"]] - t0
    nrow(frame)
  }

  if (all_pairs) {
    progress <- .mgcvst_progress("pairs", length(used) * (length(used) - 1) / 2,
                                 "pairs")
    capacity <- max(2L, min(length(used), floor(cache_bytes / state_estimate)))
    schedule <- .mgcvst_all_pair_schedule(length(used), capacity, chunk_size)
    blocks <- schedule$blocks
    for (b in seq_len(nrow(blocks))) {
      pos <- .mgcvst_block_pairs(blocks[b, ])
      window <- cbind(used[pos[, 1L]], used[pos[, 2L]])
      evaluate(window, sort(unique(as.vector(window))), b)
      progress$add(nrow(window))
      shard_files <- c(shard_files, .mgcvst_shard_file(pair_dir, b))
      shard_rows <- c(shard_rows, nrow(window))
      chunks <- chunks + 1L
      if (verbose && (chunks %% 10L == 0L || b == nrow(blocks))) {
        message("Evaluated pair block ", b, " of ", nrow(blocks), ".")
      }
    }
    progress$done()
    n_pairs <- sum(shard_rows)
    pair_schedule <- "resident_gene_tiles"
  } else {
    ord <- .mgcvst_pair_order(index, used,
      max(2L, floor(cache_bytes / state_estimate)))
    index <- index[ord, , drop = FALSE]
    n_pairs <- nrow(index)
    progress <- .mgcvst_progress("pairs", n_pairs, "pairs")
    first <- 1L
    last_probe <- proc.time()[["elapsed"]]
    while (first <= n_pairs) {
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
      progress$add(nrow(window))
      shard_files <- c(shard_files, .mgcvst_shard_file(pair_dir, first))
      shard_rows <- c(shard_rows, nrow(window))
      chunks <- chunks + 1L
      if (verbose && (chunks %% 10L == 0L || last == n_pairs)) {
        message("Evaluated pair group ", chunks, ".")
      }
      first <- last + 1L
    }
    progress$done()
    pair_schedule <- "resident_gene_blocks"
  }

  finished <- TRUE
  list(
    pair_dir = pair_dir, shards = shard_files, rows = shard_rows,
    n_pairs = n_pairs,
    failed = data.frame(
      feature_id = fit$feature_id[as.integer(names(failed))],
      error = unname(failed), stringsAsFactors = FALSE
    ),
    elapsed = elapsed,
    metadata = list(
      builds = builds,
      preparation_backend = "sparse_reduced",
      pair_schedule = pair_schedule, chunks = chunks, k = k,
      basis_sha = shared$sha, q = width, nodes_above_leading = nodes_above,
      screen_z = screen_z, screened = screened,
      cache_hits = cache$hits, cache_misses = cache$misses,
      cache_evictions = cache$evictions, cache_bytes = cache_bytes,
      resident_bytes = cache$bytes, preparation_elapsed = preparation_elapsed,
      contract = contract
    ),
    features = data.frame(feature_index = used, trace_M = trace_M,
                          trace_M2 = trace_M2)
  )
}
