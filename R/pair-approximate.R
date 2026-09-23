.mgcvst_approx_summary_open <- function(path, signature, feature_ids,
                                        resume) {
  if (!dir.exists(path) && !dir.create(path, recursive = TRUE)) {
    stop("Could not create the approximate score summary directory.")
  }
  manifest <- file.path(path, "manifest.rds")
  expected <- list(format = 3L, signature = signature, feature_ids = feature_ids)
  if (file.exists(manifest)) {
    if (!resume) stop("The approximate score summary already exists.")
    saved <- tryCatch(readRDS(manifest), error = function(e)
      stop("The approximate score summary manifest is unreadable: ",
           conditionMessage(e)))
    if (!identical(saved, expected)) {
      stop("The approximate score summary signature or feature IDs do not match.")
    }
  } else {
    if (length(list.files(path, all.files = TRUE, no.. = TRUE))) {
      stop("The approximate score summary lacks a manifest and is not empty.")
    }
    temp <- tempfile("manifest-", tmpdir = path, fileext = ".tmp")
    on.exit(if (file.exists(temp)) unlink(temp), add = TRUE)
    saveRDS(expected, temp, compress = FALSE)
    if (!file.rename(temp, manifest)) {
      stop("Could not commit the approximate score summary manifest.")
    }
  }
  list(path = normalizePath(path, winslash = "/", mustWork = TRUE),
       signature = signature, feature_ids = feature_ids)
}

.mgcvst_approx_summary_file <- function(store, feature) {
  file.path(store$path, sprintf("summary-%010d.rds", feature))
}

.mgcvst_approx_summary_native_file <- function(store, feature) {
  file.path(store$path, sprintf("summary-%010d.bin", feature))
}

.mgcvst_approx_summary_read <- function(store, feature, n_ref) {
  path <- .mgcvst_approx_summary_file(store, feature)
  native_path <- .mgcvst_approx_summary_native_file(store, feature)
  if (file.exists(path)) {
    z <- tryCatch(readRDS(path), error = function(e)
      stop("The approximate score summary is unreadable: ", basename(path),
           ": ", conditionMessage(e)))
  } else if (file.exists(native_path)) {
    z <- tryCatch(mgcvst_trace_read_cpp(
      native_path, digest::digest(store$signature, algo = "sha256"),
      store$feature_ids[feature], n_ref
    ), error = function(e) {
      stop("The native approximate score summary is unreadable: ",
           basename(native_path), ": ", conditionMessage(e))
    })
    if (!is.list(z)) {
      stop("The native approximate score summary returned an invalid record: ",
           basename(native_path), ".")
    }
    z$feature_id <- store$feature_ids[feature]
    z$signature <- store$signature
  } else {
    return(NULL)
  }
  summary_name <- if (file.exists(path)) basename(path) else basename(native_path)
  if (!is.list(z) || !identical(z$feature_id, store$feature_ids[feature]) ||
      !identical(z$signature, store$signature)) {
    stop("The approximate score summary has incompatible metadata: ",
         summary_name, ".")
  }
  if (!is.null(z$error)) {
    if (!is.character(z$error) || length(z$error) != 1L ||
        is.na(z$error) || !nzchar(z$error)) {
      stop("The approximate score summary has an invalid error record.")
    }
  } else if (!is.numeric(z$a) || !is.matrix(z$cross) ||
             !identical(dim(z$cross), c(n_ref, 4L)) ||
             (!is.null(z$self) && length(z$self) != 4L) || any(!is.finite(z$a)) ||
             any(!is.finite(z$cross)) || any(!is.finite(z$self))) {
    stop("The approximate score summary has invalid numeric data: ",
         summary_name, ".")
  }
  z
}

.mgcvst_approx_summary_write <- function(store, feature, z) {
  path <- .mgcvst_approx_summary_file(store, feature)
  native_path <- .mgcvst_approx_summary_native_file(store, feature)
  if (file.exists(path) || file.exists(native_path)) {
    stop("The approximate score summary already exists.")
  }
  z$feature_id <- store$feature_ids[feature]
  z$signature <- store$signature
  temp <- tempfile(paste0(basename(path), "-"), tmpdir = store$path,
                   fileext = ".tmp")
  on.exit(if (file.exists(temp)) unlink(temp), add = TRUE)
  saveRDS(z, temp, compress = FALSE)
  if (file.exists(path) || !file.rename(temp, path)) {
    stop("Could not commit the approximate score summary: ",
         basename(path), ".")
  }
  invisible(path)
}

# Retain signs of the reference trace matrix; higher-order trace kernels
# need not be positive semidefinite even for valid feature curvatures.
.mgcvst_approx_inverse <- function(W, tolerance) {
  W <- (W + t(W)) / 2
  E <- CppMatrix::matrixEigen(W)
  scale <- max(abs(E$values))
  if (!is.finite(scale) || scale <= 0) {
    stop("The reference trace matrix has no usable eigenvalues.")
  }
  keep <- abs(E$values) > tolerance * scale
  if (!any(keep)) stop("The reference trace matrix is numerically singular.")
  V <- E$vectors[, keep, drop = FALSE]
  inverse <- CppMatrix::matrixMultiply(
    V * rep(1 / E$values[keep], each = nrow(V)), V, transB = TRUE
  )
  list(inverse = inverse, rank = sum(keep),
       negative = sum(E$values[keep] < 0),
       dropped = sum(!keep))
}

.mgcvst_pair_approximate <- function(fit, index, pair_index, threads,
                                     chunk_size, verbose, basis = NULL,
                                     n_ref = 100L,
                                     ref_method = c("random", "score", "hyper"),
                                     ref_seed = 1L, cache_bytes = NULL,
                                     checkpoint_dir = NULL, resume = TRUE,
                                     ref_tol = 1e-6,
                                     diagnostic_pairs = 0L) {
  ref_method <- match.arg(ref_method)
  if (!is.matrix(index) || ncol(index) != 2L || !nrow(index) ||
      anyNA(index) || any(index != floor(index)) ||
      any(index < 1L) || any(index > length(fit$feature_id))) {
    stop("index must contain valid two-column feature indices.")
  }
  if (length(pair_index) != nrow(index) || anyNA(pair_index)) {
    stop("pair_index must identify every requested pair.")
  }
  for (z in list(threads, chunk_size, n_ref)) {
    if (!is.numeric(z) || length(z) != 1L || !is.finite(z) ||
        z < 1L || z != floor(z)) {
      stop("threads, chunk_size, and n_ref must be positive integers.")
    }
  }
  for (z in list(ref_seed, diagnostic_pairs)) {
    if (!is.numeric(z) || length(z) != 1L || !is.finite(z) ||
        z < 0 || z > .Machine$integer.max || z != floor(z)) {
      stop("ref_seed and diagnostic_pairs must be non-negative integers.")
    }
  }
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose) ||
      !is.logical(resume) || length(resume) != 1L || is.na(resume)) {
    stop("verbose and resume must be TRUE or FALSE.")
  }
  if (!is.numeric(ref_tol) || length(ref_tol) != 1L ||
      !is.finite(ref_tol) || ref_tol < 0 || ref_tol >= 1) {
    stop("ref_tol must be one finite value in [0, 1).")
  }
  if (!is.null(cache_bytes) && (!is.numeric(cache_bytes) ||
      length(cache_bytes) != 1L || !is.finite(cache_bytes) || cache_bytes <= 0)) {
    stop("cache_bytes must be NULL or one positive byte count.")
  }
  .mgcvst_thread_limit()
  started <- proc.time()[["elapsed"]]
  used <- sort(unique(as.vector(index)))
  sparse <- identical(fit$score_backend, "sparse")
  native <- T0 <- field_scale <- NULL
  if (sparse) {
    fit <- .inlast_sparse_prepare(fit)
    if (is.null(basis)) basis <- .inlast_sparse_observation_basis(fit)
    width <- basis$rank
    mode <- "sparse"
  } else if (identical(fit$test_engine, "spde")) {
    T0 <- .mgcvst_legacy_shared_score_factor(fit$geometry)
    field_scale <- .mgcvst_field_scale(fit)
    width <- ncol(T0)
    mode <- "legacy_native"
  } else {
    fit$.mgcvst_fixed_factors <- .mgcvst_model_fixed_factors(fit)
    native <- .mgcvst_model_dense_preparation(fit, used)
    mode <- if (is.null(native)) "model_fallback" else "model_native"
    width <- if (!is.null(native)) ncol(native$T0) else sum(vapply(
      fit$geometry$smooth[unname(fit$geometry$target)],
      function(z) ncol(z$B), integer(1L)))
  }
  signature <- if (is.null(checkpoint_dir)) {
    list(version = 2L, temporary = tempfile("mgcvst-cur-run-"))
  } else .mgcvst_pair_signature(fit, basis)
  state_store <- .mgcvst_store_open(
    if (is.null(checkpoint_dir)) NULL else file.path(checkpoint_dir, "states"),
    signature, fit$feature_id, storage = "double", resume = resume,
    encoding = "native"
  )
  if (isTRUE(state_store$temporary)) {
    on.exit(.mgcvst_store_cleanup(state_store), add = TRUE)
  }
  n_ref <- min(as.integer(n_ref), length(used))
  state_bytes <- 8 * (width^2 + width) + 2048
  q <- if (sparse) ncol(fit$score_sparse$Q) else width
  p <- if (sparse) ncol(fit$geometry$nuisance_design) else if (
    !is.null(native)) ncol(native$X) else ncol(fit$geometry$X)
  if (is.null(p)) p <- 0L
  n <- nrow(fit$working_variance)
  feature_work <- if (sparse) {
    8 * (8 * q^2 + 4 * q * p + 4 * p^2 + 4 * n + 4 * width^2)
  } else 8 * (8 * n * q + 4 * n * p + 8 * q^2 + 4 * q * p + 4 * p^2)
  reserve <- 4 * width^2 * n_ref + 24 * width^2 * threads +
    128 * nrow(index) + 64 * 1024^2
  automatic_cache <- is.null(cache_bytes)
  probe <- .mgcvst_memory_probe()
  if (automatic_cache) cache_bytes <- if (is.finite(probe$available)) {
    max(0, 0.7 * probe$available - reserve)
  } else 512 * 1024^2
  if (n_ref * state_bytes > cache_bytes) {
    stop("The reference curvature matrices exceed cache_bytes; use fewer references or more memory.")
  }
  last_probe <- proc.time()[["elapsed"]]
  batch_limit <- function(resident, remaining) {
    now <- proc.time()[["elapsed"]]
    if (now - last_probe >= 2) {
      probe <<- .mgcvst_memory_probe()
      last_probe <<- now
    }
    free <- cache_bytes - resident
    if (is.finite(probe$available)) free <- min(free, 0.7 * probe$available - reserve)
    if (!is.finite(free) || free < feature_work) {
      stop("Insufficient memory for one approximate score-state batch after ",
           "the resident references and working buffers.")
    }
    as.integer(min(max(32L, 2L * threads), remaining, floor(free / feature_work)))
  }

  units <- NULL
  unit_builds <- builds <- resume_count <- 0L
  selection_started <- proc.time()[["elapsed"]]
  scores <- NULL
  if (ref_method == "score" && n_ref < length(used)) {
    units <- .mgcvst_unit_store_open(file.path(state_store$path, "units"),
      signature, fit$feature_id, resume = resume,
      kind = if (sparse) "sparse" else "dense")
    scores <- matrix(NA_real_, width, length(fit$feature_id),
                      dimnames = list(NULL, fit$feature_id))
    first <- 1L
    while (first <= length(used)) {
      count <- batch_limit(as.numeric(object.size(scores)), length(used) - first + 1L)
      ids <- used[first:min(length(used), first + count - 1L)]
      has_state <- vapply(ids, function(id) .mgcvst_store_has(state_store, id), logical(1L))
      missing <- !has_state & !vapply(ids, function(id)
        .mgcvst_unit_store_has(units, id), logical(1L))
      if (any(missing)) {
        todo <- ids[missing]
        if (sparse) {
          z <- .inlast_sparse_units(fit, todo, threads = threads)
        } else {
          z <- lapply(todo, function(id) {
            if (mode == "legacy_native") {
              F <- sqrt(field_scale[id]) * T0
              op <- .rkhs_score_operator_factor(F, fit$working_variance[, id],
                                                fit$geometry$X)
              out <- list(operator = op, target = list(F))
            } else out <- .mgcvst_model_operator(fit, id)
            F <- do.call(cbind, out$target)
            Pe <- .mgcvst_model_apply_P(out$operator, fit$working_error[, id])
            out$a <- as.numeric(.magic_mm(F, matrix(Pe, ncol = 1L), transA = TRUE))
            out
          })
        }
        for (k in seq_along(todo)) .mgcvst_unit_store_write(units, todo[k], z[[k]])
        unit_builds <- unit_builds + length(todo)
        rm(z)
      }
      for (id in ids) {
        if (has_state[match(id, ids)]) {
          z <- .mgcvst_store_read(state_store, id)
          if (!is.null(z$error)) stop("Score landmark preparation failed for ",
                                      fit$feature_id[id], ": ", z$error)
          scores[, id] <- z$a
          next
        }
        z <- .mgcvst_unit_store_read(units, id)
        if (!is.null(z$error)) stop("Score landmark preparation failed for ",
                                    fit$feature_id[id], ": ", z$error)
        scores[, id] <- if (sparse) as.numeric(.magic_mm(
          basis$coordinate, matrix(z$a, ncol = 1L), transA = TRUE)) else z$a
      }
      first <- first + length(ids)
    }
  }
  references <- .mgcvst_landmark_select(fit, used, n_ref, ref_method,
                                       as.integer(ref_seed), scores = scores)
  selection_elapsed <- proc.time()[["elapsed"]] - selection_started
  rm(scores)
  summary_signature <- list(fit = signature, references = references,
    reference_method = ref_method, seed = as.integer(ref_seed), version = 2L)
  summary_path <- file.path(state_store$path, paste0("cur-",
    digest::digest(summary_signature, algo = "md5")))
  summary_store <- .mgcvst_approx_summary_open(summary_path, summary_signature,
                                               fit$feature_id, resume)
  build_batch <- function(ids, load = TRUE) {
    missing <- !vapply(ids, function(id) .mgcvst_store_has(state_store, id), logical(1L))
    resume_count <<- resume_count + sum(!missing)
    if (any(missing)) {
      todo <- ids[missing]
      if (is.null(units)) {
        z <- .mgcvst_pair_build_batch(fit, todo, threads, basis, mode,
                                      native, T0, field_scale)
      } else {
        u <- lapply(todo, function(id) .mgcvst_unit_store_read(units, id))
        names(u) <- unname(fit$feature_id[todo])
        if (sparse) {
          z <- .mgcvst_pair_build_batch(fit, todo, threads, basis, mode,
                                        sparse_units = u)
        } else {
          z <- lapply(u, function(z) {
            F <- do.call(cbind, z$target)
            PF <- .mgcvst_model_apply_P(z$operator, F)
            M <- .magic_mm(F, PF, transA = TRUE)
            list(a = z$a, M = (M + t(M)) / 2,
                 width = vapply(z$target, ncol, integer(1L)))
          })
        }
      }
      for (k in seq_along(todo)) .mgcvst_store_write(state_store, todo[k], z[[k]])
      builds <<- builds + length(todo)
      rm(z)
      if (!is.null(units)) rm(u)
    }
    if (load) lapply(ids, function(id) .mgcvst_store_read(state_store, id)) else invisible(NULL)
  }

  ref_state <- vector("list", n_ref)
  first <- 1L
  while (first <= n_ref) {
    resident <- sum(vapply(ref_state, function(z) as.numeric(object.size(z)), numeric(1L)))
    count <- batch_limit(resident, n_ref - first + 1L)
    pos <- first:min(n_ref, first + count - 1L)
    ref_state[pos] <- build_batch(references[pos])
    first <- max(pos) + 1L
  }
  failed <- vapply(ref_state, function(z) !is.null(z$error), logical(1L))
  if (any(failed)) stop("Reference score construction failed: ",
                         paste(fit$feature_id[references[failed]], collapse = ", "))
  ref_M <- lapply(ref_state, `[[`, "M")
  saved_ref <- lapply(references, function(id)
    .mgcvst_approx_summary_read(summary_store, id, n_ref))
  W <- array(NA_real_, c(n_ref, n_ref, 4L))
  if (all(vapply(saved_ref, function(z) !is.null(z), logical(1L)))) {
    for (k in seq_len(n_ref)) W[k, , ] <- saved_ref[[k]]$cross
  } else {
    pairs <- which(upper.tri(matrix(FALSE, n_ref, n_ref), diag = TRUE), arr.ind = TRUE)
    values <- mgcvst_pair_trace_powers_cpp(ref_M, pairs, maxPower = 4L, threads = threads)
    for (k in seq_len(nrow(pairs))) {
      i <- pairs[k, 1L]; j <- pairs[k, 2L]
      W[i, j, ] <- W[j, i, ] <- values[k, ]
    }
  }
  if (any(!is.finite(W))) stop("Reference trace moments are non-finite.")
  summaries <- vector("list", length(fit$feature_id))
  for (k in seq_along(references)) {
    id <- references[k]
    summaries[[id]] <- list(a = ref_state[[k]]$a, self = W[k, k, ],
                             cross = matrix(W[k, , ], n_ref, 4L))
    if (is.null(saved_ref[[k]])) .mgcvst_approx_summary_write(summary_store, id, summaries[[id]])
  }
  rm(ref_state, saved_ref)
  others <- setdiff(used, references)
  ref_bytes <- sum(vapply(ref_M, function(z) as.numeric(object.size(z)), numeric(1L)))
  complete <- logical(length(others))
  for (k in seq_along(others)) {
    id <- others[k]
    z <- .mgcvst_approx_summary_read(summary_store, id, n_ref)
    complete[k] <- !is.null(z)
    if (complete[k]) {
      if (!.mgcvst_store_has(state_store, id)) {
        stop("A completed landmark summary is missing its double state shard.")
      }
      summaries[[id]] <- z
    }
  }
  resume_count <- resume_count + sum(complete)
  todo <- others[!complete]
  first <- 1L
  while (first <= length(todo)) {
    count <- batch_limit(ref_bytes + as.numeric(object.size(summaries)),
                          length(todo) - first + 1L)
    ids <- todo[first:min(length(todo), first + count - 1L)]
    build_batch(ids, load = FALSE)
    first <- first + length(ids)
  }
  stream <- NULL
  if (length(todo)) {
    if (verbose) message("Submitting ", length(todo), " genes to the native landmark queue.")
    stream <- mgcvst_landmark_stream_cpp(
      vapply(todo, function(id) .mgcvst_store_file(state_store, id), character(1L)),
      unname(fit$feature_id[todo]),
      digest::digest(state_store$signature, algo = "sha256"), ref_M,
      vapply(todo, function(id) .mgcvst_approx_summary_native_file(summary_store, id),
             character(1L)),
      digest::digest(summary_store$signature, algo = "sha256"), threads = threads
    )
    for (id in todo) {
      summaries[[id]] <- .mgcvst_approx_summary_read(summary_store, id, n_ref)
      if (is.null(summaries[[id]])) stop("A native landmark task did not produce its summary.")
    }
  }
  preparation_elapsed <- proc.time()[["elapsed"]] - started
  rm(ref_M, native, T0, field_scale)
  fit$.mgcvst_fixed_factors <- NULL

  failed <- vapply(summaries[used], function(z) !is.null(z$error), logical(1L))
  valid_ids <- used[!failed]
  if (!length(valid_ids)) stop("No valid gene score states are available.")
  C <- array(NA_real_, c(length(valid_ids), n_ref, 4L))
  for (k in seq_along(valid_ids)) C[k, , ] <- summaries[[valid_ids[k]]]$cross
  left <- right <- vector("list", 4L)
  ranks <- negatives <- dropped <- integer(4L)
  for (k in seq_len(4L)) {
    Wk <- matrix(W[, , k], n_ref, n_ref)
    d <- sqrt(diag(Wk))
    if (any(!is.finite(d)) || any(d <= 0)) stop("A reference self trace is invalid.")
    normalized <- sweep(matrix(C[, , k], length(valid_ids), n_ref), 2L, d, "/")
    inverse <- .mgcvst_approx_inverse(Wk / (d %o% d), ref_tol)
    left[[k]] <- t(.magic_mm(normalized, inverse$inverse))
    right[[k]] <- t(normalized)
    ranks[k] <- inverse$rank; negatives[k] <- inverse$negative; dropped[k] <- inverse$dropped
  }
  location <- integer(length(fit$feature_id))
  location[valid_ids] <- seq_along(valid_ids)
  A <- do.call(cbind, lapply(summaries[valid_ids], `[[`, "a"))
  scales <- matrix(1, length(valid_ids), 4L)
  diagnostic_seed <- as.integer((as.double(ref_seed) + 104729) %% .Machine$integer.max)
  hold <- .mgcvst_landmark_holdout(index, references, as.integer(diagnostic_pairs),
                                  diagnostic_seed)
  hold_moments <- matrix(NA_real_, length(hold), 4L)
  hold_logp <- rep(NA_real_, length(hold))
  result <- data.frame(pair_index = pair_index, score = rep(NA_real_, nrow(index)),
    information = NA_real_, effective_rank = NA_real_, p_value = NA_real_,
    error_message = NA_character_, stringsAsFactors = FALSE)
  reconstruction_started <- proc.time()[["elapsed"]]
  pair_path <- .mgcvst_pair_checkpoint(state_store, index, pair_index,
    calibration = list(version = 1L, summary = summary_signature, ref_tol = ref_tol))
  chunks <- resumed_pairs <- 0L
  first <- 1L
  while (first <= nrow(index)) {
    saved <- .mgcvst_pair_checkpoint_read(pair_path, first, pair_index)
    if (!is.null(saved)) {
      result[seq.int(first, saved$last), ] <- saved$result
      resumed_pairs <- resumed_pairs + nrow(saved$result)
      chunks <- chunks + 1L
      first <- saved$last + 1L
      next
    }
    rows <- first:min(nrow(index), first + chunk_size - 1L)
    mapped <- matrix(location[index[rows, ]], ncol = 2L)
    good <- mapped[, 1L] > 0L & mapped[, 2L] > 0L
    for (k in which(!good)) {
      bad <- index[rows[k], ][location[index[rows[k], ]] == 0L]
      result$error_message[rows[k]] <- paste(vapply(summaries[bad], `[[`,
        character(1L), "error"), collapse = " | ")
    }
    if (any(good)) {
      selected <- rows[good]
      values <- mgcvst_pair_lowrank_cpp(A, left, right, scales,
                                        mapped[good, , drop = FALSE], threads)
      result$score[selected] <- values[, 1L]
      m <- values[, 2:5, drop = FALSE]
      okay <- apply(is.finite(m) & m > 0, 1L, all) & m[, 1L] > 1e-10
      result$error_message[selected[!okay]] <- "Approximate Liu trace moments are invalid."
      if (any(okay)) {
        z <- selected[okay]
        mm <- m[okay, , drop = FALSE]
        liu <- .liu_squared_score_moments(abs(values[okay, 1L]), mm[, 1L], mm[, 2L],
                                          mm[, 3L], mm[, 4L])
        result$information[z] <- mm[, 1L]
        result$effective_rank[z] <- mm[, 1L]^2 / mm[, 2L]
        result$p_value[z] <- liu$p_value
        invalid <- !is.finite(liu$p_value) | liu$p_value < 0 | liu$p_value > 1
        result$error_message[z[invalid]] <- "Approximate Liu calibration returned an invalid p-value."
      }
    }
    .mgcvst_pair_checkpoint_write(pair_path, first, max(rows), result[rows, , drop = FALSE])
    chunks <- chunks + 1L
    first <- max(rows) + 1L
  }
  if (length(hold)) {
    mapped <- matrix(location[index[hold, ]], ncol = 2L)
    good <- mapped[, 1L] > 0L & mapped[, 2L] > 0L
    if (any(good)) {
      values <- mgcvst_pair_lowrank_cpp(A, left, right, scales,
        mapped[good, , drop = FALSE], threads)
      hold_moments[good, ] <- values[, 2:5, drop = FALSE]
      okay <- apply(is.finite(hold_moments) & hold_moments > 0, 1L, all) &
        hold_moments[, 1L] > 1e-10
      if (any(okay)) {
        m <- hold_moments[okay, , drop = FALSE]
        liu <- .liu_squared_score_moments(abs(result$score[hold[okay]]),
          m[, 1L], m[, 2L], m[, 3L], m[, 4L])
        hold_logp[okay] <- .mgcvst_liu_logp(liu)
      }
    }
  }
  reconstruction_elapsed <- proc.time()[["elapsed"]] - reconstruction_started
  diagnostics <- NULL
  diagnostic_elapsed <- 0
  if (length(hold)) {
    diagnostic_index <- index[hold, , drop = FALSE]
    diagnostic_file <- file.path(summary_store$path, paste0("diagnostics-log1-",
      digest::digest(diagnostic_index, algo = "md5"), ".rds"))
    if (file.exists(diagnostic_file)) {
      saved <- readRDS(diagnostic_file)
      if (!identical(saved$version, 1L) || !identical(saved$index, diagnostic_index) ||
          !identical(saved$signature, summary_signature)) {
        stop("The saved landmark diagnostic pairs or signature do not match.")
      }
      exact <- saved$exact
    } else {
      exact <- .mgcvst_landmark_exact(state_store, diagnostic_index,
        threads, chunk_size, cache_bytes)
      diagnostic_elapsed <- exact$elapsed
      temp <- tempfile("diagnostics-", tmpdir = summary_store$path, fileext = ".tmp")
      saveRDS(list(version = 1L, index = diagnostic_index, signature = summary_signature,
                   exact = exact), temp, compress = FALSE)
      if (!file.rename(temp, diagnostic_file)) stop("Could not commit landmark diagnostics.")
    }
    diagnostics <- list(pair_index = pair_index[hold], exact_builds = exact$builds,
      error = .mgcvst_landmark_error_table(hold_moments, exact$moments,
                                           hold_logp, exact$log10_p),
      exact_errors = exact$error_message)
  }
  list(result = result,
    elapsed = reconstruction_elapsed + diagnostic_elapsed,
    metadata = list(path = checkpoint_dir, references = references,
      reference_method = ref_method, reference_seed = ref_seed, reference_count = n_ref,
      reference_tol = ref_tol, reference_rank = ranks, reference_negative = negatives,
      reference_dropped = dropped, builds = builds, unit_builds = unit_builds,
      resume_count = resume_count,
      preparation_backend = if (!is.null(units) && !sparse) "dense_score_units_serial" else mode,
      gene_submission = "native_file_queue", landmark_queue = stream, approximate = TRUE,
      chunks = chunks, resumed_pairs = resumed_pairs,
      cache_bytes = cache_bytes, preparation_elapsed = preparation_elapsed,
      selection_elapsed = selection_elapsed, reconstruction_elapsed = reconstruction_elapsed,
      diagnostic_elapsed = diagnostic_elapsed,
      total_elapsed = proc.time()[["elapsed"]] - started,
      diagnostic_seed = diagnostic_seed, diagnostics = diagnostics,
      storage = "double", landmark_precision = "float32",
      calibration = "Liu landmark trace CUR"),
    internal = list(score = A, left = left, right = right, scale = scales,
      cross = C, reference_trace = W, valid_ids = valid_ids))
}
