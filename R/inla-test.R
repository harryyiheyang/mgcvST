# Sparse INLA downstream evaluation uses one OpenMP layer in the manager.
.mgcvst_inla_downstream <- function(fit) {
  identical(fit$estimator, "INLA")
}

.mgcvst_inla_require_sparse <- function(fit) {
  if (!identical(fit$score_backend, "sparse")) {
    stop("INLA downstream tests require the single-global sparse score ",
         "geometry built by inlaST.set()/inlaST.estimate(); the dense INLA ",
         "score no longer exists.")
  }
  invisible(NULL)
}

.mgcvst_inla_pair_chunk_size <- function(fit, memory_bytes = 512 * 1024^2,
                                         basis = NULL) {
  if (is.null(basis)) basis <- .inlast_sparse_observation_basis(fit)
  r <- basis$rank
  as.integer(max(1L, min(128L, floor(memory_bytes / (2 * 8 * r^2)))))
}

.mgcvst_inla_serial_backend <- function(BPPARAM) {
  if (!inherits(BPPARAM, "SerialParam")) {
    stop(
      "Sparse INLA downstream tests use C++ OpenMP; BPPARAM must be SerialParam()."
    )
  }
  invisible(NULL)
}

.mgcvst_inla_test_pairs <- function(fit, index, pair_index, threads,
                                    chunk_size, verbose, coverage = 0.995,
                                    full_rank = FALSE, basis = NULL,
                                    cache_bytes = 512 * 1024^2) {
  fit <- .inlast_sparse_prepare(fit)
  if (is.null(basis)) basis <- .inlast_sparse_observation_basis(
    fit, coverage = coverage, full_rank = full_rank
  )
  state_estimate <- 8 * basis$rank^2
  starts <- seq.int(1L, nrow(index), by = min(chunk_size, 128L))
  result <- vector("list", length(starts))
  elapsed <- 0
  cache_limit <- cache_bytes
  cache <- new.env(parent = emptyenv())
  cache$state <- list()
  cache$bytes <- 0
  cache$last <- numeric()
  cache$clock <- 0
  cache$hits <- 0L
  cache$misses <- 0L
  cache$evictions <- 0L
  cache$unit_builds <- 0L
  cache$materializations <- 0L
  cache$transient_unit_bytes <- 0
  cache$unit_build_elapsed <- 0
  cache$reduced_materialize_elapsed <- 0
  for (b in seq_along(starts)) {
    rows <- starts[b]:min(nrow(index), starts[b] + min(chunk_size, 128L) - 1L)
    block_result <- vector("list", length(rows))
    for (pair_pos in seq_along(rows)) {
    pair_row <- rows[pair_pos]
    block_index <- index[pair_row, , drop = FALSE]
    block_used <- sort(unique(as.vector(block_index)))
    block_states <- vector("list", length(block_used))
    missing <- integer()
    for (feature_pos in seq_along(block_used)) {
      key <- as.character(block_used[feature_pos])
      if (!is.null(cache$state[[key]])) {
        cache$clock <- cache$clock + 1
        cache$last[key] <- cache$clock
        cache$hits <- cache$hits + 1L
        block_states[[feature_pos]] <- cache$state[[key]]
      } else {
        cache$misses <- cache$misses + 1L
        missing <- c(missing, feature_pos)
      }
    }
    if (length(missing)) {
      for (position in missing) {
        active <- as.character(block_used)
        while (cache$bytes + state_estimate > cache_limit) {
          drop <- setdiff(names(cache$last), active)
          if (!length(drop)) break
          drop <- drop[which.min(cache$last[drop])]
          cache$bytes <- cache$bytes - as.numeric(object.size(cache$state[[drop]]))
          cache$state[[drop]] <- NULL
          cache$last <- cache$last[names(cache$last) != drop]
          cache$evictions <- cache$evictions + 1L
        }
        t_unit <- proc.time()[["elapsed"]]
        unit <- .inlast_sparse_units(
          fit, block_used[position], threads = threads
        )[[1L]]
        cache$unit_build_elapsed <- cache$unit_build_elapsed +
          proc.time()[["elapsed"]] - t_unit
        cache$unit_builds <- cache$unit_builds + 1L
        cache$transient_unit_bytes <- max(
          cache$transient_unit_bytes, as.numeric(object.size(unit))
        )
        if (!is.null(unit$error) && nzchar(unit$error)) {
          block_states[[position]] <- list(error = unit$error)
          rm(unit)
          next
        }
        t_materialize <- proc.time()[["elapsed"]]
        state <- .inlast_sparse_materialize_reduced(
          fit, list(unit), basis, threads = threads
        )[[1L]]
        cache$reduced_materialize_elapsed <- cache$reduced_materialize_elapsed +
          proc.time()[["elapsed"]] - t_materialize
        rm(unit)
        cache$materializations <- cache$materializations + 1L
        key <- as.character(block_used[position])
        state_bytes <- as.numeric(object.size(state))
        while (cache$bytes + state_bytes > cache_limit) {
          drop <- setdiff(names(cache$last), active)
          if (!length(drop)) break
          drop <- drop[which.min(cache$last[drop])]
          cache$bytes <- cache$bytes - as.numeric(object.size(cache$state[[drop]]))
          cache$state[[drop]] <- NULL
          cache$last <- cache$last[names(cache$last) != drop]
          cache$evictions <- cache$evictions + 1L
        }
        if (state_bytes <= cache_limit && cache$bytes + state_bytes <= cache_limit) {
          cache$state[[key]] <- state
          cache$bytes <- cache$bytes + state_bytes
          cache$clock <- cache$clock + 1
          cache$last[key] <- cache$clock
        }
        block_states[[position]] <- state
      }
    }
    failed <- vapply(block_states, function(z) {
      !is.null(z$error) && length(z$error) == 1L && !is.na(z$error) && nzchar(z$error)
    }, logical(1L))
    summaries <- list(
      used = block_used,
      a = lapply(block_states, `[[`, "a"),
      H = lapply(block_states, `[[`, "M"),
      has_summary = !failed,
      error_message = vapply(block_states, function(z) {
        if (is.null(z$error) || !length(z$error)) NA_character_ else z$error
      }, character(1L)),
      elapsed = 0
    )
    z <- .mgcvst_liu_pairs(
      block_index, pair_index[pair_row], fit$feature_id, summaries,
      threads, 1L, FALSE
    )
    block_result[[pair_pos]] <- z$result
    elapsed <- elapsed + z$elapsed
    }
    result[[b]] <- do.call(rbind, block_result)
    if (verbose && (b %% 10L == 0L || b == length(starts))) {
      message("Evaluated sparse INLA Liu block ", b, " of ", length(starts), ".")
    }
  }
  out <- do.call(rbind, result)
  names(out)[names(out) == "score"] <- "signed_score"
  names(out)[names(out) == "p_value"] <- "p_two_sided"
  valid <- is.finite(out$p_two_sided) & out$p_two_sided >= 0 &
    out$p_two_sided <= 1
  out$p_positive <- out$p_negative <- NA_real_
  out$p_positive[valid] <- ifelse(
    out$signed_score[valid] >= 0,
    out$p_two_sided[valid] / 2,
    1 - out$p_two_sided[valid] / 2
  )
  out$p_negative[valid] <- ifelse(
    out$signed_score[valid] <= 0,
    out$p_two_sided[valid] / 2,
    1 - out$p_two_sided[valid] / 2
  )
  attr(out, "inla_pairwise") <- list(
    q = ncol(fit$score_sparse$Q), r = basis$rank,
    target_coverage = basis$coverage, kept_coverage = basis$kept,
    tail = basis$tail,
    basis = "constrained_observation_kernel_A_Qg_inverse_At",
    cache_limit_bytes = cache_limit, cache_bytes = cache$bytes,
    cache_hits = cache$hits, cache_misses = cache$misses,
    cache_evictions = cache$evictions, unit_builds = cache$unit_builds,
    materializations = cache$materializations,
    transient_unit_bytes = cache$transient_unit_bytes,
    unit_build_elapsed = cache$unit_build_elapsed,
    reduced_materialize_elapsed = cache$reduced_materialize_elapsed,
    liu_elapsed = elapsed, unit_cache = "none",
    pair_schedule = "one_pair_microblocks_no_pair_level_OpenMP_batch"
  )
  list(result = out, elapsed = elapsed)
}

.mgcvst_inla_wgcna_scores <- function(fit, used, threads, verbose) {
  group <- "global"
  fit <- .inlast_sparse_prepare(fit)
  blocks <- split(seq_along(used), ceiling(seq_along(used) / 32L))
  states <- vector("list", length(used))
  for (rows in blocks) {
    states[rows] <- .inlast_sparse_batch(
      fit, features = used[rows], threads = threads, score_only = TRUE
    )
  }
  if (length(states) != length(used)) {
    stop("The sparse INLA batch returned an incompatible feature count.")
  }
  failed <- vapply(states, function(z) {
    !is.null(z$error) && length(z$error) == 1L && !is.na(z$error) && nzchar(z$error)
  }, logical(1L))
  if (any(failed)) {
    stop("Sparse INLA score construction failed: ", paste(
      paste0(fit$feature_id[used[failed]], ": ",
             vapply(states[failed], `[[`, character(1L), "error")),
      collapse = " | "
    ))
  }
  coordinate_width <- length(states[[1L]]$a)
  normalization <- as.integer(states[[1L]]$normalization)
  width <- stats::setNames(coordinate_width, "global")
  A <- do.call(cbind, lapply(states, `[[`, "a"))
  colnames(A) <- fit$feature_id[used]
  if (nrow(A) != coordinate_width || any(!is.finite(A))) {
    stop("The sparse INLA batch returned invalid score coordinates.")
  }
  if (verbose) {
    message("Constructed scores for ", length(used), " features with C++ OpenMP.")
  }
  list(
    A = A, group = group, width = width,
    normalization = normalization,
    feature_id = fit$feature_id[used]
  )
}
