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

.mgcvst_inla_pair_chunk_size <- function(fit, memory_bytes = 512 * 1024^2) {
  q <- ncol(fit$score_sparse$Q)
  if (length(q) != 1L || !is.finite(q) || q < 1L) {
    stop("The sparse INLA score geometry has an invalid dimension.")
  }
  as.integer(max(1L, min(128L, floor(memory_bytes / (2 * 8 * q^2)))))
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
                                    chunk_size, verbose) {
  fit <- .inlast_sparse_prepare(fit)
  used <- sort(unique(as.vector(index)))
  feature_batch_size <- 32L
  feature_blocks <- split(
    seq_along(used), ceiling(seq_along(used) / feature_batch_size)
  )
  units <- vector("list", length(used))
  unit_bytes <- 0
  memory_limit <- 512 * 1024^2
  unit_dir <- NULL
  for (rows in feature_blocks) {
    batch <- .inlast_sparse_units(
      fit, features = used[rows], threads = threads
    )
    if (length(batch) != length(rows)) {
      stop("The sparse INLA unit batch returned an incompatible feature count.")
    }
    batch_bytes <- sum(vapply(batch, function(x) as.numeric(object.size(x)),
                              numeric(1L)))
    if (is.null(unit_dir) && unit_bytes + batch_bytes <= memory_limit) {
      units[rows] <- batch
      unit_bytes <- unit_bytes + batch_bytes
    } else {
      if (is.null(unit_dir)) {
        unit_dir <- tempfile("mgcvst-inla-units-")
        dir.create(unit_dir)
        unit_dir <- normalizePath(unit_dir, winslash = "/", mustWork = TRUE)
        temp_root <- normalizePath(tempdir(), winslash = "/", mustWork = TRUE)
        if (!startsWith(unit_dir, paste0(temp_root, "/"))) {
          stop("The INLA temporary unit cache is outside the R temporary directory.")
        }
        on.exit(unlink(unit_dir, recursive = TRUE, force = TRUE), add = TRUE)
        resident <- which(!vapply(units, is.null, logical(1L)))
        for (j in resident) {
          saveRDS(units[[j]], file.path(unit_dir, paste0(j, ".rds")),
                  compress = FALSE)
          units[j] <- list(NULL)
        }
      }
      for (j in seq_along(rows)) {
        saveRDS(batch[[j]], file.path(unit_dir, paste0(rows[j], ".rds")),
                compress = FALSE)
      }
    }
  }
  load_units <- if (is.null(unit_dir)) {
    function(position) units[position]
  } else {
    function(position) lapply(position, function(j) {
      readRDS(file.path(unit_dir, paste0(j, ".rds")))
    })
  }
  starts <- seq.int(1L, nrow(index), by = chunk_size)
  result <- vector("list", length(starts))
  elapsed <- 0
  for (b in seq_along(starts)) {
    rows <- starts[b]:min(nrow(index), starts[b] + chunk_size - 1L)
    block_index <- index[rows, , drop = FALSE]
    block_used <- sort(unique(as.vector(block_index)))
    block_position <- match(block_used, used)
    block_units <- load_units(block_position)
    unit_failed <- vapply(block_units, function(z) {
      !is.null(z$error) && length(z$error) == 1L && !is.na(z$error) && nzchar(z$error)
    }, logical(1L))
    block_states <- vector("list", length(block_units))
    good_units <- which(!unit_failed)
    if (length(good_units)) {
      block_states[good_units] <- .inlast_sparse_materialize(
        fit, block_units[good_units], threads = threads
      )
    }
    if (any(unit_failed)) {
      block_states[unit_failed] <- lapply(block_units[unit_failed], function(z) {
        list(error = z$error)
      })
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
      block_index, pair_index[rows], fit$feature_id, summaries,
      threads, length(rows), FALSE
    )
    result[[b]] <- z$result
    elapsed <- elapsed + z$elapsed
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
