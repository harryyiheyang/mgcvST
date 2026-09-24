# Exact score_liu for the sparse INLA backend, backed by the fp16 pair kernel
# (src/inla_fp16.cpp). Per-gene reduced curvature is built once in double and
# packed to fp16; pairs are scored by fp32 GEMM with double traces and Liu
# tail. With pairs = NULL, gene pairs are generated and scored in left-gene
# blocks and never materialized as one pair matrix or character data frame.

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

# Error text for pairs touching a gene whose fp16 state failed to build.
.mgcvst_inla_fp16_pair_error <- function(feature_id1, feature_id2, failed) {
  m1 <- match(feature_id1, failed$feature_id)
  m2 <- match(feature_id2, failed$feature_id)
  msg1 <- ifelse(!is.na(m1), paste0(feature_id1, ": ", failed$error[m1]), NA_character_)
  msg2 <- ifelse(!is.na(m2), paste0(feature_id2, ": ", failed$error[m2]), NA_character_)
  ifelse(!is.na(msg1) & !is.na(msg2), paste(msg1, msg2, sep = " | "),
        ifelse(!is.na(msg1), msg1, msg2))
}

# Exact fp16 Liu pairs for an explicit `index` (two-column feature indices),
# in the .mgcvst_pair_pipeline() result contract: information and
# effective_rank are not computed on this path.
.mgcvst_inla_fp16_pairs_explicit <- function(fit, state, index, pair_index, threads) {
  local <- matrix(match(index, state$used), ncol = 2L)
  lo <- pmin(local[, 1L], local[, 2L])
  hi <- pmax(local[, 1L], local[, 2L])
  ord <- order(lo)
  out <- mgcvst_fp16_pairs_cpp(state$cache, state$used, 1L, state$n,
                               left = lo[ord], right = hi[ord], threads = threads)
  n <- nrow(index)
  score <- numeric(n)
  mlog10p <- numeric(n)
  score[ord] <- out$score
  mlog10p[ord] <- out$mlog10p
  error_message <- rep(NA_character_, n)
  bad <- is.na(mlog10p)
  if (any(bad)) {
    msg <- .mgcvst_inla_fp16_pair_error(
      fit$feature_id[index[bad, 1L]], fit$feature_id[index[bad, 2L]], state$failed
    )
    error_message[bad] <- ifelse(!is.na(msg), msg,
      "The fp16 Liu pair calibration is unavailable for this pair.")
  }
  data.frame(
    pair_index = pair_index, score = score, information = NA_real_,
    effective_rank = NA_real_, p_value = 10^(-mlog10p),
    error_message = error_message, stringsAsFactors = FALSE
  )
}

# Entry point used by .mgcvst_inla_test_pairs() for liu_approximation = "exact"
# with an explicit pair list; gene states are built once and released here.
.mgcvst_inla_fp16_test_explicit <- function(fit, index, pair_index, threads,
                                            verbose, basis, checkpoint_dir, resume) {
  used <- sort(unique(as.vector(index)))
  t0 <- proc.time()[["elapsed"]]
  state <- .mgcvst_inla_fp16_states(fit, used, basis, threads, checkpoint_dir,
                                    resume, verbose)
  on.exit(state$release(), add = TRUE)
  preparation_elapsed <- proc.time()[["elapsed"]] - t0
  t1 <- proc.time()[["elapsed"]]
  result <- .mgcvst_inla_fp16_pairs_explicit(fit, state, index, pair_index, threads)
  elapsed <- proc.time()[["elapsed"]] - t1
  list(
    result = result, elapsed = elapsed,
    metadata = list(
      path = if (state$temporary) NULL else state$root,
      builds = state$built, resume_count = state$resumed,
      resumed_pairs = 0L, pair_path = NULL,
      preparation_backend = "fp16_sparse_openmp",
      pair_schedule = "fp16_explicit_pairs",
      chunks = 1L, cache_hits = NA_integer_, cache_misses = NA_integer_,
      cache_evictions = NA_integer_, cache_bytes = state$bytes,
      resident_bytes = state$bytes, preparation_elapsed = preparation_elapsed,
      signature = NULL, storage = "fp16"
    )
  )
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

# BH (or another stats::p.adjust method) over every tested pair, read back
# from the pair shards; only numeric columns are ever materialized. Skipped
# when the mlog10p/i/j columns would not comfortably fit in available memory.
.mgcvst_inla_fp16_bh <- function(shard_paths, total_pairs, q.value, FDR, method) {
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
  valid <- is.finite(mlog10p)
  adjusted <- rep(NA_real_, total_pairs)
  discovered <- rep(FALSE, total_pairs)
  if (any(valid)) {
    adjusted[valid] <- if (identical(method, "BH")) {
      .mgcvst_bh_mlog10p(mlog10p[valid])
    } else {
      -log10(stats::p.adjust(10^(-mlog10p[valid]), method))
    }
    discovered[valid] <- adjusted[valid] >= -log10(q.value)
  }
  list(computed = TRUE, i = i, j = j, mlog10p = mlog10p,
       adjusted_mlog10p = adjusted, discovered = discovered)
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

#' Streaming exact fp16 score_liu pairs for a sparse INLA fit
#'
#' Every gene pair among the available genes of `fit` is tested with the exact
#' Liu calibration, using the fp16 pair kernel of [inlaST.test()]. Unlike
#' [inlaST.test()] with an explicit `pairs` argument, gene pairs are generated
#' and scored in left-gene blocks and are never assembled as one pair matrix
#' or per-pair character table; each block is written as one Parquet shard
#' (integer `i`, `j`; double `score`, `mlog10p = -log10(p)`).
#'
#' @param fit An object returned by [inlaST.estimate()].
#' @param checkpoint_dir Directory for fp16 gene-state and pair shards.
#'   `NULL` uses a temporary directory removed on exit.
#' @param resume Reuse compatible completed gene-state and pair shards.
#' @param threads OpenMP threads for building gene states and scoring pairs.
#' @param chunk_size Approximate number of pairs per Parquet shard; block
#'   boundaries are a deterministic function of `chunk_size` and the number of
#'   available genes, so resumed runs regenerate the same shard schedule.
#' @param verbose Report progress.
#' @param q.value,FDR,method BH adjustment (or another
#'   `stats::p.adjust.methods` entry) of `mlog10p` over every tested pair,
#'   read back from the pair shards; skipped when that would not comfortably
#'   fit in available memory (see `$bh$reason`).
#' @return A list: `shards` (Parquet shard paths, `i`/`j`/`score`/`mlog10p`),
#'   `feature_id` (the fit's full feature-identifier lookup vector; shard `i`
#'   and `j` are indices into it directly), `failed` (a `feature_id`/`error`
#'   table of genes whose fp16 state could not be built; pairs touching them
#'   have `mlog10p = NA` and the run continues), `bh` (BH results, or
#'   `$computed = FALSE` and `$reason`), `n_genes`, `total_pairs`, and timing.
#' @export
inlaST.fp16Liu <- function(fit, checkpoint_dir = NULL, resume = TRUE,
                           threads = 1L, chunk_size = 4000000L, verbose = FALSE,
                           q.value = 0.05, FDR = TRUE, method = "BH") {
  if (!inherits(fit, "mgcvST_model_fit") || !identical(fit$estimator, "INLA")) {
    stop("inlaST.fp16Liu() requires a fit from inlaST.estimate().")
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
  used <- which(.mgcvst_feature_available(fit))
  if (length(used) < 2L) stop("At least two available INLA genes are required.")

  t0 <- proc.time()[["elapsed"]]
  state <- .mgcvst_inla_fp16_states(fit, used, basis, threads, checkpoint_dir,
                                    resume, verbose)
  on.exit(state$release(), add = TRUE)
  state_elapsed <- proc.time()[["elapsed"]] - t0

  t1 <- proc.time()[["elapsed"]]
  stream <- .mgcvst_inla_fp16_stream(state, threads, chunk_size, verbose)
  pair_elapsed <- proc.time()[["elapsed"]] - t1

  bh <- .mgcvst_inla_fp16_bh(stream$shard_paths, stream$total_pairs, q.value, FDR, method)

  list(
    shards = stream$shard_paths, feature_id = fit$feature_id,
    n_genes = length(used), total_pairs = stream$total_pairs,
    failed = state$failed, bh = bh,
    threshold = list(q_value = q.value, FDR = FDR,
                     adjustment_method = if (FDR) method else "none"),
    checkpoint_dir = if (state$temporary) NULL else state$root,
    timing = list(
      gene_states_elapsed = state_elapsed, pair_elapsed = pair_elapsed,
      elapsed = proc.time()[["elapsed"]] - t0, threads = threads,
      gene_states_built = state$built, gene_states_resumed = state$resumed,
      pair_blocks = stream$blocks, pair_blocks_built = stream$built_blocks,
      pair_blocks_resumed = stream$blocks - stream$built_blocks,
      resumed_pairs = stream$resumed_pairs
    ),
    call = match.call()
  )
}
