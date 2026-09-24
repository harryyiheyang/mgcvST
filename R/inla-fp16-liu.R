# Exact score_liu for the sparse INLA backend, backed by the fp16 pair kernel
# (src/inla_fp16.cpp). Per-gene reduced curvature is built once in double and
# packed to fp16; pairs are scored by fp32 GEMM with double traces and Liu
# tail. Both an explicit (normalized, bounds-checked) pair block and the full
# gene-pair universe (pairs = NULL, generated and scored in left-gene blocks,
# never materialized as one pair matrix or character data frame) share the
# same compact result columns: integer i, j; double score, mlog10p.

# Signature covering everything that determines the fp16 gene states of `used`.
.mgcvst_inla_fp16_signature <- function(fit, used, basis) {
  digest::digest(list(version = 1L, format = "fp16_liu",
                      compact = .inlast_compact_signature(fit, used, basis),
                      r = basis$rank), algo = "sha256")
}

# Per-gene batch size that keeps one build call's added fp16 bytes within a
# fraction of available memory; mgcvst_fp16_build_cpp enforces the exact safe
# line and stops with an explicit error if memory is actually insufficient.
.mgcvst_inla_fp16_batch_size <- function(r) {
  per_gene <- r * (r + 1) + 8 * r
  probe <- .mgcvst_memory_probe()
  available <- if (is.finite(probe$available)) probe$available else 8 * 1024^3
  as.integer(max(1L, min(2048L, floor(0.3 * available / max(1, per_gene)))))
}

# Build or resume the fp16 gene states of `used` (global feature indices) into
# one cache, checkpointed as raw fp16 shards under checkpoint_dir/states. With
# checkpoint_dir = NULL, a temporary directory is used and removed on release.
.mgcvst_inla_fp16_states <- function(fit, used, basis, threads, checkpoint_dir,
                                     resume, verbose) {
  n <- length(used)
  r <- basis$rank
  cache <- mgcvst_fp16_cache_cpp(n, r)
  released <- FALSE
  release <- function() {
    if (!released) {
      mgcvst_fp16_cache_release_cpp(cache)
      released <<- TRUE
    }
  }
  sig <- .mgcvst_inla_fp16_signature(fit, used, basis)
  batch_size <- .mgcvst_inla_fp16_batch_size(r)

  temporary <- is.null(checkpoint_dir)
  root <- if (temporary) tempfile("mgcvst-fp16-") else checkpoint_dir
  if (!dir.exists(root) && !dir.create(root, recursive = TRUE)) {
    stop("Could not create the fp16 checkpoint directory.")
  }
  root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  state_dir <- file.path(root, "states")
  if (!dir.exists(state_dir) && !dir.create(state_dir, recursive = TRUE)) {
    stop("Could not create the fp16 gene-state checkpoint directory.")
  }
  manifest_path <- file.path(root, "fp16-states-manifest.rds")
  manifest <- list(version = 1L, signature = sig,
                   feature_id = fit$feature_id[used], r = r,
                   batch_size = batch_size)
  if (file.exists(manifest_path)) {
    if (!resume) {
      stop("An fp16 checkpoint already exists in ", root, "; use a new directory.")
    }
    if (!identical(readRDS(manifest_path), manifest)) {
      stop("The fp16 checkpoint in ", root, " was written for a different fit, ",
           "basis, gene universe, or batch size.")
    }
  } else {
    tmp <- tempfile("fp16-states-manifest-", tmpdir = root, fileext = ".tmp")
    saveRDS(manifest, tmp, compress = FALSE)
    if (!file.rename(tmp, manifest_path)) {
      stop("Could not commit the fp16 gene-state checkpoint manifest.")
    }
  }

  built <- 0L
  resumed <- 0L
  starts <- if (n) seq.int(1L, n, by = batch_size) else integer()
  for (b in seq_along(starts)) {
    first <- starts[b]
    last <- min(n, first + batch_size - 1L)
    slots <- first:last
    ids <- used[slots]
    shard <- file.path(state_dir, sprintf("state-%06d.bin", b))
    loaded <- FALSE
    if (file.exists(shard)) {
      loaded <- tryCatch({
        mgcvst_fp16_read_cpp(cache, slots, ids, fit$feature_id[ids], sig, shard, NA_real_)
        TRUE
      }, error = function(e) {
        if (verbose) message("Rebuilding stale fp16 shard ", basename(shard), ": ",
                             conditionMessage(e))
        unlink(shard)
        FALSE
      })
    }
    if (loaded) {
      resumed <- resumed + length(slots)
      next
    }
    z <- .inlast_compact_inputs(fit, ids)
    geometry <- fit$score_sparse
    mgcvst_fp16_build_cpp(
      cache, slots, geometry$cache$general_A, geometry$cache$general_Q,
      as.numeric(geometry$constraint), geometry$cache$general_X,
      z$B, z$C, z$O, z$family, z$size, z$dispersion, z$tau, z$a,
      basis$coordinate, basis$basis, NA_real_, threads,
      nuisance_precision = z$nuisance_precision
    )
    mgcvst_fp16_write_cpp(cache, slots, ids, fit$feature_id[ids], sig, shard)
    built <- built + length(slots)
    if (verbose) {
      message("Built fp16 score states for ", last, " of ", n, " genes.")
    }
  }
  info <- mgcvst_fp16_cache_info_cpp(cache)
  failed <- which(info$state == 2L)
  list(
    cache = cache, release = release, used = used, n = n, r = r,
    built = built, resumed = resumed, bytes = info$bytes,
    failed = data.frame(feature_id = fit$feature_id[used[failed]],
                        error = info$error[failed], stringsAsFactors = FALSE),
    root = root, temporary = temporary
  )
}

# Normalize an explicit `pairs` argument (feature IDs or 1-based indices) into
# a deduplicated, bounds-checked block of global feature-index pairs with
# i < j.
.mgcvst_inla_fp16_normalize_pairs <- function(fit, pairs) {
  index <- .mgcvst_pair_index(pairs, fit$feature_id)
  lo <- pmin(index[, 1L], index[, 2L])
  hi <- pmax(index[, 1L], index[, 2L])
  unique(cbind(i = lo, j = hi))
}

# Exact fp16 Liu score/mlog10p for a deduplicated explicit block of global
# feature-index pairs (`pairs_u`, columns i < j), scored against the fp16
# gene states of `state`; failed genes carry through as NA mlog10p, exactly
# as in the streaming path.
.mgcvst_inla_fp16_pairs_result <- function(state, pairs_u, threads) {
  local <- cbind(match(pairs_u[, "i"], state$used), match(pairs_u[, "j"], state$used))
  ord <- order(local[, 1L])
  out <- mgcvst_fp16_pairs_cpp(state$cache, state$used, 1L, state$n,
                               left = local[ord, 1L], right = local[ord, 2L],
                               threads = threads)
  score <- numeric(nrow(pairs_u))
  mlog10p <- numeric(nrow(pairs_u))
  score[ord] <- out$score
  mlog10p[ord] <- out$mlog10p
  data.frame(i = pairs_u[, "i"], j = pairs_u[, "j"], score = score,
             mlog10p = mlog10p)
}

# Read or write the one parquet shard caching an explicit pair block, in the
# same shard format as .mgcvst_inla_fp16_stream(); only used when a
# checkpoint_dir is supplied for explicit pairs (there is no parquet or
# streaming requirement otherwise).
.mgcvst_inla_fp16_explicit_shard <- function(state, pairs_u, threads, resume, verbose) {
  if (!requireNamespace("arrow", quietly = TRUE)) {
    stop("Checkpointing explicit score_liu pairs requires the arrow package.")
  }
  pairs_dir <- file.path(state$root, "pairs")
  if (!dir.exists(pairs_dir) && !dir.create(pairs_dir, recursive = TRUE)) {
    stop("Could not create the fp16 pair checkpoint directory.")
  }
  manifest_path <- file.path(pairs_dir, "manifest.rds")
  manifest <- list(version = 1L, explicit = TRUE,
                   signature = digest::digest(pairs_u, algo = "sha256"))
  shard <- file.path(pairs_dir, "pairs-explicit.parquet")
  if (file.exists(manifest_path)) {
    if (!identical(readRDS(manifest_path), manifest)) {
      stop("The fp16 pair checkpoint in ", pairs_dir, " was written for a ",
           "different explicit pair block; use a new checkpoint_dir.")
    }
    if (resume && file.exists(shard)) {
      if (verbose) message("Reusing the checkpointed explicit pair shard.")
      return(list(result = as.data.frame(arrow::read_parquet(shard)), shard = shard))
    }
  } else {
    tmp <- tempfile("pairs-manifest-", tmpdir = pairs_dir, fileext = ".tmp")
    saveRDS(manifest, tmp, compress = FALSE)
    if (!file.rename(tmp, manifest_path)) {
      stop("Could not commit the fp16 pair checkpoint manifest.")
    }
  }
  result <- .mgcvst_inla_fp16_pairs_result(state, pairs_u, threads)
  tmp <- tempfile("pairs-", tmpdir = pairs_dir, fileext = ".parquet.tmp")
  arrow::write_parquet(result, tmp)
  if (!file.rename(tmp, shard)) {
    stop("Could not commit the fp16 pair shard ", basename(shard), ".")
  }
  list(result = result, shard = shard)
}

# Deterministic left-gene block boundaries so that resumed runs regenerate the
# same shard schedule: each block's total right-partner count stays at or
# below `pairs_per_block`, except a block of one left gene may exceed it.
.mgcvst_fp16_pair_blocks <- function(n, pairs_per_block) {
  blocks <- list()
  left <- 1L
  while (left < n) {
    last <- left
    total <- n - left
    while (last + 1L < n) {
      extra <- n - (last + 1L)
      if (total + extra > pairs_per_block) break
      last <- last + 1L
      total <- total + extra
    }
    blocks[[length(blocks) + 1L]] <- c(left, last, total)
    left <- last + 1L
  }
  blocks
}

# BH step-up directly on mlog10p = -log10(p), avoiding underflow for very
# small p; equivalent to stats::p.adjust(10^-mlog10p, "BH") on -log10 scale.
.mgcvst_bh_mlog10p <- function(mlog10p) {
  m <- length(mlog10p)
  o <- order(mlog10p)
  rank_from_small_p <- rev(seq_len(m))
  raw <- mlog10p[o] - log10(m) + log10(rank_from_small_p)
  out <- numeric(m)
  out[o] <- pmax(cummax(raw), 0)
  out
}

# BH (or another stats::p.adjust method) adjustment of mlog10p already held in
# memory; shared by the explicit (in-memory) and streamed (shard-backed) BH
# helpers below.
.mgcvst_inla_fp16_bh_adjust <- function(mlog10p, q.value, method) {
  valid <- is.finite(mlog10p)
  adjusted <- rep(NA_real_, length(mlog10p))
  discovered <- rep(FALSE, length(mlog10p))
  if (any(valid)) {
    adjusted[valid] <- if (identical(method, "BH")) {
      .mgcvst_bh_mlog10p(mlog10p[valid])
    } else {
      -log10(stats::p.adjust(10^(-mlog10p[valid]), method))
    }
    discovered[valid] <- adjusted[valid] >= -log10(q.value)
  }
  list(adjusted_mlog10p = adjusted, discovered = discovered)
}

# BH results for an explicit pair block already materialized in memory.
.mgcvst_inla_fp16_bh_from_vectors <- function(i, j, mlog10p, total_pairs,
                                              q.value, FDR, method) {
  if (!FDR || !total_pairs) {
    return(list(computed = FALSE,
               reason = if (!FDR) "FDR = FALSE" else "no tested pairs"))
  }
  adj <- .mgcvst_inla_fp16_bh_adjust(mlog10p, q.value, method)
  list(computed = TRUE, i = i, j = j, mlog10p = mlog10p,
       adjusted_mlog10p = adj$adjusted_mlog10p, discovered = adj$discovered)
}

# BH results over every tested pair, read back from the pair shards; only
# numeric columns are ever materialized. Skipped when the mlog10p/i/j columns
# would not comfortably fit in available memory.
.mgcvst_inla_fp16_bh_from_shards <- function(shard_paths, total_pairs, q.value,
                                             FDR, method) {
  if (!FDR || !total_pairs) {
    return(list(computed = FALSE,
               reason = if (!FDR) "FDR = FALSE" else "no tested pairs"))
  }
  probe <- .mgcvst_memory_probe()
  if (is.finite(probe$available) && 24 * total_pairs > 0.4 * probe$available) {
    return(list(computed = FALSE,
               reason = "the pairwise mlog10p column would exceed the safe memory line"))
  }
  i <- integer(total_pairs)
  j <- integer(total_pairs)
  mlog10p <- numeric(total_pairs)
  at <- 0L
  for (path in shard_paths) {
    z <- arrow::read_parquet(path, col_select = c("i", "j", "mlog10p"))
    k <- nrow(z)
    rows <- at + seq_len(k)
    i[rows] <- z$i
    j[rows] <- z$j
    mlog10p[rows] <- z$mlog10p
    at <- at + k
  }
  adj <- .mgcvst_inla_fp16_bh_adjust(mlog10p, q.value, method)
  list(computed = TRUE, i = i, j = j, mlog10p = mlog10p,
       adjusted_mlog10p = adj$adjusted_mlog10p, discovered = adj$discovered)
}

# Streaming exact fp16 Liu pairs over every pair of `state$used`: each
# left-gene block is generated and scored on the fly and written as one
# parquet shard (i, j, score, mlog10p); no pair matrix or per-pair strings are
# ever held in memory. Completed shards are skipped on resume.
.mgcvst_inla_fp16_stream <- function(state, threads, chunk_size, verbose) {
  if (!requireNamespace("arrow", quietly = TRUE)) {
    stop("Streaming score_liu pairs (pairs = NULL) requires the arrow package.")
  }
  n <- state$n
  pairs_dir <- file.path(state$root, "pairs")
  if (!dir.exists(pairs_dir) && !dir.create(pairs_dir, recursive = TRUE)) {
    stop("Could not create the fp16 pair checkpoint directory.")
  }
  manifest_path <- file.path(pairs_dir, "manifest.rds")
  manifest <- list(version = 1L, n = n, chunk_size = chunk_size)
  if (file.exists(manifest_path)) {
    if (!identical(readRDS(manifest_path), manifest)) {
      stop("The fp16 pair checkpoint in ", pairs_dir, " was written for a ",
           "different gene universe or chunk_size; use a new checkpoint_dir.")
    }
  } else {
    tmp <- tempfile("pairs-manifest-", tmpdir = pairs_dir, fileext = ".tmp")
    saveRDS(manifest, tmp, compress = FALSE)
    if (!file.rename(tmp, manifest_path)) {
      stop("Could not commit the fp16 pair checkpoint manifest.")
    }
  }
  blocks <- .mgcvst_fp16_pair_blocks(n, chunk_size)
  shard_paths <- character(length(blocks))
  total_pairs <- 0
  resumed_pairs <- 0
  built_blocks <- 0L
  for (b in seq_along(blocks)) {
    left <- blocks[[b]][1L]
    last <- blocks[[b]][2L]
    count <- blocks[[b]][3L]
    shard <- file.path(pairs_dir, sprintf("pairs-%06d.parquet", b))
    shard_paths[b] <- shard
    total_pairs <- total_pairs + count
    if (file.exists(shard)) {
      resumed_pairs <- resumed_pairs + count
      next
    }
    out <- mgcvst_fp16_pairs_cpp(state$cache, state$used, left, last, threads = threads)
    df <- data.frame(i = out$i, j = out$j, score = out$score, mlog10p = out$mlog10p)
    tmp <- tempfile("pairs-", tmpdir = pairs_dir, fileext = ".parquet.tmp")
    arrow::write_parquet(df, tmp)
    if (!file.rename(tmp, shard)) stop("Could not commit the fp16 pair shard ", basename(shard), ".")
    built_blocks <- built_blocks + 1L
    if (verbose) {
      message("Wrote fp16 pair shard ", b, " of ", length(blocks), ".")
    }
  }
  list(shard_paths = shard_paths, total_pairs = total_pairs,
       resumed_pairs = resumed_pairs, built_blocks = built_blocks,
       blocks = length(blocks), pairs_dir = pairs_dir)
}

# Exact fp16 score_liu pairs for a sparse INLA fit, in the compact result
# format shared by `pairs = NULL` (every gene pair, streamed to Parquet
# shards) and an explicit pair block (materialized in memory, and also
# checkpointed to one Parquet shard when `checkpoint_dir` is supplied). Called
# by inlaST.test() and by mgcvST.test()'s dispatch on INLA fits; internal.
#
# Returns a list: `result` (a data frame with integer `i`, `j` and double
# `score`, `mlog10p`; only set for an explicit pair block), `shards` (Parquet
# shard paths with the same four columns; only set when pairs were streamed
# or an explicit block was checkpointed), `feature_id` (the fit's full
# feature-identifier lookup vector; `i`/`j` are indices into it directly),
# `failed` (a `feature_id`/`error` table of genes whose fp16 state could not
# be built; pairs touching them have `mlog10p = NA` and the run continues),
# `bh` (BH results, or `$computed = FALSE` and `$reason`), `n_genes`,
# `total_pairs`, `checkpoint_dir`, and timing.
.mgcvst_inla_fp16_run <- function(fit, pairs = NULL, checkpoint_dir = NULL,
                                  resume = TRUE, threads = 1L,
                                  chunk_size = 4000000L, verbose = FALSE,
                                  q.value = 0.05, FDR = TRUE, method = "BH") {
  if (!inherits(fit, "mgcvST_model_fit") || !identical(fit$estimator, "INLA")) {
    stop("The exact fp16 score_liu path requires a fit from inlaST.estimate().")
  }
  .mgcvst_inla_require_sparse(fit)
  if (!is.numeric(threads) || length(threads) != 1L || !is.finite(threads) ||
      threads < 1 || threads != floor(threads)) {
    stop("threads must be one positive integer.")
  }
  threads <- as.integer(threads)
  if (!is.numeric(chunk_size) || length(chunk_size) != 1L ||
      !is.finite(chunk_size) || chunk_size < 1) {
    stop("chunk_size must be one positive number of pairs.")
  }
  chunk_size <- as.integer(min(chunk_size, .Machine$integer.max))
  if (!is.logical(resume) || length(resume) != 1L || is.na(resume)) {
    stop("resume must be TRUE or FALSE.")
  }
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("verbose must be TRUE or FALSE.")
  }
  q.value <- as.numeric(q.value)
  if (length(q.value) != 1L || !is.finite(q.value) || q.value <= 0 || q.value > 1) {
    stop("q.value must be one finite value in (0, 1].")
  }
  if (!is.logical(FDR) || length(FDR) != 1L || is.na(FDR)) {
    stop("FDR must be TRUE or FALSE.")
  }
  if (!is.character(method) || length(method) != 1L || is.na(method) ||
      !(method %in% stats::p.adjust.methods)) {
    stop("method must be one of stats::p.adjust.methods.")
  }
  .mgcvst_thread_limit()
  RhpcBLASctl::blas_set_num_threads(1L)

  fit <- .inlast_sparse_prepare(fit)
  basis <- .inlast_sparse_observation_basis(fit, coverage = 0.995)
  explicit <- !is.null(pairs)
  pairs_u <- NULL
  if (explicit) {
    pairs_u <- .mgcvst_inla_fp16_normalize_pairs(fit, pairs)
    used <- sort(unique(c(pairs_u[, "i"], pairs_u[, "j"])))
  } else {
    used <- which(.mgcvst_feature_available(fit))
    if (length(used) < 2L) stop("At least two available INLA genes are required.")
  }

  t0 <- proc.time()[["elapsed"]]
  state <- .mgcvst_inla_fp16_states(fit, used, basis, threads, checkpoint_dir,
                                    resume, verbose)
  on.exit(state$release(), add = TRUE)
  state_elapsed <- proc.time()[["elapsed"]] - t0

  t1 <- proc.time()[["elapsed"]]
  result <- NULL
  shards <- NULL
  pair_blocks <- pair_blocks_built <- pair_blocks_resumed <- resumed_pairs <- NULL
  if (explicit) {
    if (is.null(checkpoint_dir)) {
      result <- .mgcvst_inla_fp16_pairs_result(state, pairs_u, threads)
    } else {
      shard <- .mgcvst_inla_fp16_explicit_shard(state, pairs_u, threads, resume, verbose)
      result <- shard$result
      shards <- shard$shard
    }
    total_pairs <- nrow(result)
    pair_elapsed <- proc.time()[["elapsed"]] - t1
    bh <- .mgcvst_inla_fp16_bh_from_vectors(result$i, result$j, result$mlog10p,
                                            total_pairs, q.value, FDR, method)
  } else {
    stream <- .mgcvst_inla_fp16_stream(state, threads, chunk_size, verbose)
    shards <- stream$shard_paths
    total_pairs <- stream$total_pairs
    pair_elapsed <- proc.time()[["elapsed"]] - t1
    bh <- .mgcvst_inla_fp16_bh_from_shards(shards, total_pairs, q.value, FDR, method)
    pair_blocks <- stream$blocks
    pair_blocks_built <- stream$built_blocks
    pair_blocks_resumed <- stream$blocks - stream$built_blocks
    resumed_pairs <- stream$resumed_pairs
  }

  list(
    result = result, shards = shards, feature_id = fit$feature_id,
    n_genes = length(used), total_pairs = total_pairs,
    failed = state$failed, bh = bh,
    threshold = list(q_value = q.value, FDR = FDR,
                     adjustment_method = if (FDR) method else "none"),
    checkpoint_dir = if (state$temporary) NULL else state$root,
    timing = list(
      gene_states_elapsed = state_elapsed, pair_elapsed = pair_elapsed,
      elapsed = proc.time()[["elapsed"]] - t0, threads = threads,
      gene_states_built = state$built, gene_states_resumed = state$resumed,
      pair_blocks = pair_blocks, pair_blocks_built = pair_blocks_built,
      pair_blocks_resumed = pair_blocks_resumed, resumed_pairs = resumed_pairs
    ),
    call = match.call()
  )
}
