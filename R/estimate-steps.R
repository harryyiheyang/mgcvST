# Two-step estimation shared by mgcvST.estimate() and inlaST.estimate().
#
# Step 1 fits only the null model of every feature and returns the Stage 1
# null-first p-value; step 2 fits the spatial model of the selected features
# only. This file holds what both branches share: the Stage 1 q-values, the
# feature selection, the content digests that tie a later step to the Y of
# step 1, and resumable chunk checkpoints written by the workers.

# Format of the fit objects of the two-step estimators. Fits without it were
# estimated before 0.0.1.9032.
.mgcvst_fit_format <- 2L

# A sparse INLA fit estimated before 0.0.1.9032 lacks the stored mean mu_bar and
# the two-step bookkeeping that the tests and the later steps read, so it is
# refused rather than reused. An mgcv fit of the earlier format holds every
# quantity the mgcv tests use and is accepted.
.mgcvst_check_fit_format <- function(fit) {
  format <- fit$format
  if (!is.null(format) && (!is.numeric(format) || length(format) != 1L ||
                           format > .mgcvst_fit_format)) {
    stop("The fit was written by a newer version of mgcvST; update the package.",
         call. = FALSE)
  }
  if (identical(fit$estimator, "INLA") &&
      (is.null(format) || format < .mgcvst_fit_format || !is.numeric(fit$mu_bar))) {
    stop("This inlaST fit was estimated before mgcvST 0.0.1.9032 and lacks the ",
         "stored mean mu_bar that the pair test reads; re-run inlaST.estimate().",
         call. = FALSE)
  }
  invisible(fit)
}

# Stage 1 q-values on the natural scale; NA where the p-value is missing.
# The adjustment is the log-space kernel of the pair tests, so a p-value below
# the double range keeps its rank.
.mgcvst_stage1_q <- function(p, adjust) {
  p <- as.numeric(p)
  lp <- rep(NA_real_, length(p))
  ok <- !is.na(p) & p >= 0 & p <= 1
  lp[ok] <- log(p[ok])
  exp(.mgcvst_log_adjust(lp, adjust)$log_q)
}

.mgcvst_check_q_value <- function(q.value) {
  q.value <- as.numeric(q.value)
  if (length(q.value) != 1L || !is.finite(q.value) ||
      q.value <= 0 || q.value > 1) {
    stop("q.value must be one finite value in (0, 1].")
  }
  q.value
}

# Indices of the features that receive a spatial model. `spatial` is
# "discoveries" (Stage 1 q <= q.value), "all", "none", a vector of feature
# IDs or one-based indices, or a logical vector with one value per feature.
.mgcvst_select_spatial <- function(spatial, feature_id, q, q.value) {
  p <- length(feature_id)
  if (is.null(spatial)) stop("spatial must not be NULL; use \"none\" for no spatial fits.")
  if (is.character(spatial) && length(spatial) == 1L &&
      spatial %in% c("discoveries", "all", "none")) {
    index <- switch(spatial,
      discoveries = which(!is.na(q) & q <= q.value),
      all = seq_len(p), none = integer())
  } else if (is.logical(spatial)) {
    if (length(spatial) != p || anyNA(spatial)) {
      stop("A logical spatial selection needs one non-missing value per feature.")
    }
    index <- which(spatial)
  } else if (is.character(spatial)) {
    index <- match(spatial, feature_id)
    if (anyNA(index)) {
      stop("spatial contains unknown feature IDs: ",
           paste(utils::head(spatial[is.na(index)], 5L), collapse = ", "), ".")
    }
  } else if (is.numeric(spatial)) {
    if (anyNA(spatial) || any(spatial != floor(spatial)) ||
        any(spatial < 1) || any(spatial > p)) {
      stop("Numeric spatial entries must be valid one-based feature indices.")
    }
    index <- as.integer(spatial)
  } else {
    stop("spatial must be \"discoveries\", \"all\", \"none\", feature IDs, ",
         "feature indices or a logical vector.")
  }
  sort(unique(as.integer(index)))
}

# One digest per feature row of Y, so that a later step can verify that it is
# given the responses that step 1 used.
.mgcvst_row_digests <- function(Y) {
  vapply(seq_len(nrow(Y)), function(j) digest::digest(Y[j, ], algo = "md5"),
         character(1L))
}

.mgcvst_check_responses <- function(fit, Y, index) {
  if (!is.matrix(Y) || nrow(Y) != length(fit$feature_id) ||
      ncol(Y) != fit$n_observation) {
    stop("Y must be the feature-by-observation matrix used by step 1.")
  }
  if (is.null(fit$y_digest)) {
    stop("The fit does not record the responses of step 1; re-estimate it.")
  }
  now <- vapply(index, function(j) digest::digest(Y[j, ], algo = "md5"),
                character(1L))
  bad <- index[now != fit$y_digest[index]]
  if (length(bad)) {
    stop("Y differs from the responses of step 1 for ",
         paste(utils::head(fit$feature_id[bad], 5L), collapse = ", "),
         if (length(bad) > 5L) paste0(" and ", length(bad) - 5L, " more"), ".")
  }
  invisible(NULL)
}

# ---- resumable chunk checkpoints --------------------------------------------

# Open (or resume) the checkpoint directory of an estimation. The manifest
# records everything that determines the chunk results; a directory written for
# another model, control or chunking is refused, never resumed.
.mgcvst_chunk_store <- function(path, kind, signature, resume = TRUE) {
  if (is.null(path)) return(NULL)
  if (!is.character(path) || length(path) != 1L || is.na(path) || !nzchar(path)) {
    stop("checkpoint_dir must be NULL or one directory name.")
  }
  if (!is.logical(resume) || length(resume) != 1L || is.na(resume)) {
    stop("resume must be TRUE or FALSE.")
  }
  if (!dir.exists(path) && !dir.create(path, recursive = TRUE)) {
    stop("Could not create the estimation checkpoint directory.")
  }
  path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  manifest <- file.path(path, "estimation-manifest.rds")
  record <- list(format = .mgcvst_fit_format, kind = kind, signature = signature)
  if (file.exists(manifest)) {
    if (!resume) {
      stop("The estimation checkpoint ", path,
           " already exists; use resume = TRUE or a new checkpoint_dir.")
    }
    old <- tryCatch(readRDS(manifest), error = function(e) NULL)
    if (!identical(old$format, .mgcvst_fit_format) || !identical(old$kind, kind)) {
      stop("The estimation checkpoint ", path, " was written by another ",
           "estimator or by a version before 0.0.1.9032; use a new checkpoint_dir.")
    }
    if (!identical(old$signature, signature)) {
      stop("The estimation checkpoint ", path, " was written for a different ",
           "model, offset or controls; use a new checkpoint_dir.")
    }
  } else {
    if (length(list.files(path, all.files = TRUE, no.. = TRUE))) {
      stop("The estimation checkpoint directory ", path,
           " has no manifest and is not empty.")
    }
    tmp <- tempfile("manifest-", tmpdir = path, fileext = ".tmp")
    saveRDS(record, tmp, compress = FALSE)
    if (!file.rename(tmp, manifest)) stop("Could not commit the checkpoint manifest.")
  }
  list(path = path, kind = kind)
}

# Key of one chunk: the step, its features and their responses.
.mgcvst_chunk_key <- function(step, index, y_digest) {
  substr(digest::digest(list(step, as.integer(index), y_digest[index]),
                        algo = "sha256"), 1L, 24L)
}

.mgcvst_chunk_file <- function(store, step, key) {
  if (is.null(store)) return(NA_character_)
  file.path(store$path, sprintf("%s-%s.rds", step, key))
}

# Called by the worker when a chunk is complete: the result becomes a
# checkpoint before it is returned to the manager.
.mgcvst_chunk_save <- function(file, key, result) {
  if (is.null(file) || is.na(file)) return(invisible(NULL))
  tmp <- tempfile("chunk-", tmpdir = dirname(file), fileext = ".tmp")
  on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
  saveRDS(list(key = key, result = result), tmp, compress = FALSE)
  if (file.exists(file)) unlink(file)
  if (!file.rename(tmp, file)) stop("Could not commit the estimation chunk ", basename(file), ".")
  invisible(NULL)
}

# Read the completed chunks of a step and run the missing ones. `FUN` is a
# worker function receiving a payload (with chunk_file and chunk_key) and
# `...`; a result is a list with one entry per feature of its chunk.
.mgcvst_run_chunks <- function(payloads, step, store, y_digest, BPPARAM, FUN, ...) {
  n <- length(payloads)
  keys <- vapply(payloads, function(z) .mgcvst_chunk_key(step, z$index, y_digest),
                 character(1L))
  results <- vector("list", n)
  todo <- integer()
  for (k in seq_len(n)) {
    file <- .mgcvst_chunk_file(store, step, keys[k])
    if (!is.na(file) && file.exists(file)) {
      z <- tryCatch(readRDS(file), error = function(e) NULL)
      if (!is.list(z) || !identical(z$key, keys[k]) ||
          length(z$result) != length(payloads[[k]]$index)) {
        stop("The estimation checkpoint chunk ", basename(file),
             " is damaged; delete it to recompute the chunk.")
      }
      results[[k]] <- z$result
    } else {
      payloads[[k]]$chunk_file <- if (is.na(file)) NULL else file
      payloads[[k]]$chunk_key <- keys[k]
      todo <- c(todo, k)
    }
  }
  if (length(todo)) {
    results[todo] <- BiocParallel::bplapply(payloads[todo], FUN, ..., BPPARAM = BPPARAM)
  }
  list(results = results, resumed = n - length(todo), chunks = n)
}

# Contiguous feature chunks of at most `chunk_size` features. The default is
# one chunk per worker; with a checkpoint directory it is capped at 50 features,
# so that a lost chunk costs little.
.mgcvst_feature_chunks <- function(index, chunk_size, BPPARAM, checkpoint = FALSE) {
  if (!length(index)) return(list())
  if (is.null(chunk_size)) {
    workers <- max(1L, min(length(index), BiocParallel::bpworkers(BPPARAM)))
    chunk_size <- ceiling(length(index) / workers)
    if (checkpoint) chunk_size <- min(chunk_size, 50L)
  }
  split(index, ceiling(seq_along(index) / chunk_size))
}
