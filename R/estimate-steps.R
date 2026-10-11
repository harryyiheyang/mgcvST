# Two-step estimation shared by mgcvST.estimate() and inlaST.estimate().
#
# Step 1 fits only the null model of every feature and returns the Stage 1
# null-first p-value; step 2 fits the spatial model of the selected features
# only. This file holds what both branches share: the Stage 1 q-values, the
# feature selection, the content digests that tie a later step to the Y of
# step 1, and resumable chunk checkpoints written by the workers.

# Format of the fit objects of the two-step estimators. Fits without it were
# estimated before 0.0.1.9032; format 3 (0.0.1.9034) adds the effective degrees
# of freedom of the spatial smooth to the diagnostics.
.mgcvst_fit_format <- 3L

# A fit whose spatial smooth has fewer effective degrees of freedom than this
# is degenerate: the gene is fitted with essentially no spatial field (the
# smooth is penalized to the boundary, or to less than a plane), and every pair
# with such a gene gets p = 1 and status 4. One constant for both estimators.
#
# Chosen from the copula SVG simulation (11 datasets, 11,000 genes, 50 x 50
# spots, q = 132): the edf of the genes with field variance v = 0.30 / 0.10 /
# 0.03 has medians 82 / 56 / 26 and minima 49 / 28 / 4e-5, and the edf of the
# genes without a field has median 3e-4 and maximum 11. At 3, 0% of the strong
# and medium genes and 0.76% of the weak genes (v = 0.03) are flagged, against
# 96% of the genes without a field. The nine pairs of genes without a field
# among the 9033 BY discoveries of these datasets (p between 1e-5 and 1e-7) each
# have a gene below 3.
.mgcvst_edf_min <- 3

# Flag of the degenerate spatial fits; missing where the effective degrees of
# freedom are missing.
.mgcvst_degenerate <- function(edf) as.numeric(edf) < .mgcvst_edf_min

# Degenerate spatial fits as a logical vector over the features of a fit; a
# fit without the flag has none.
.mgcvst_degenerate_features <- function(fit) {
  d <- fit$diagnostics$spatial_degenerate
  if (is.null(d)) rep(FALSE, length(fit$feature_id)) else d %in% TRUE
}

# The pair test reads the effective degrees of freedom of the spatial smooth.
.mgcvst_check_edf <- function(fit, estimator = "mgcvST") {
  if (is.null(fit$format) || fit$format < 3L ||
      is.null(fit$diagnostics$spatial_degenerate)) {
    stop("This fit was estimated before mgcvST 0.0.1.9034 and lacks the effective ",
         "degrees of freedom of the spatial smooth (diagnostics$edf_spatial) that ",
         "the pair test reads; re-run ", estimator, ".estimate().", call. = FALSE)
  }
  invisible(fit)
}

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
      (is.null(format) || format < 2L || !is.numeric(fit$mu_bar))) {
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

# The digest of one feature row of Y. It is taken over the plain double values,
# so that neither the names nor the storage type of Y matter.
.mgcvst_row_digest <- function(Y, j) {
  digest::digest(as.numeric(Y[j, ]), algo = "md5")
}

# One digest per feature row of Y, so that a later step can verify that it is
# given the responses that step 1 used.
.mgcvst_row_digests <- function(Y) {
  vapply(seq_len(nrow(Y)), function(j) .mgcvst_row_digest(Y, j), character(1L))
}

# The one validation of chunk_size shared by the estimators and the add-later
# functions: NULL, or one positive integer (not truncated).
.mgcvst_check_chunk_size <- function(chunk_size) {
  if (is.null(chunk_size)) return(NULL)
  if (!is.numeric(chunk_size) || length(chunk_size) != 1L ||
      !is.finite(chunk_size) || chunk_size < 1 || chunk_size != floor(chunk_size) ||
      chunk_size > .Machine$integer.max) {
    stop("chunk_size must be one positive integer.")
  }
  as.integer(chunk_size)
}

# Validate the response matrix in blocks of rows: finite and, for counts,
# non-negative and integer-valued. Only block-sized temporaries are made, and
# the matrix is returned as given (an integer matrix stays integer; its blocks
# are converted to double when a payload is built).
.mgcvst_check_response_matrix <- function(Y, counts = FALSE, block_elements = 4e6) {
  if (!is.matrix(Y)) Y <- as.matrix(Y)
  if (length(dim(Y)) != 2L || !nrow(Y) || !ncol(Y) ||
      !(is.numeric(Y) || is.logical(Y))) {
    stop("Y must be a non-empty finite numeric feature-by-observation matrix.")
  }
  step <- max(1L, as.integer(block_elements %/% ncol(Y)))
  for (first in seq.int(1L, nrow(Y), by = step)) {
    block <- Y[first:min(nrow(Y), first + step - 1L), , drop = FALSE]
    if (anyNA(block) || (is.double(block) && any(!is.finite(block)))) {
      stop("Y must be a non-empty finite numeric feature-by-observation matrix.")
    }
    if (counts && (any(block < 0) ||
                   (is.double(block) && any(block != round(block))))) {
      stop("Count responses must be non-negative integers.")
    }
  }
  Y
}

# Rows `index` of Y as a double matrix (a copy of the block only).
.mgcvst_double_rows <- function(Y, index) {
  block <- Y[index, , drop = FALSE]
  if (!is.double(block)) storage.mode(block) <- "double"
  block
}

.mgcvst_check_responses <- function(fit, Y, index) {
  if (!is.matrix(Y) || nrow(Y) != length(fit$feature_id) ||
      ncol(Y) != fit$n_observation) {
    stop("Y must be the feature-by-observation matrix used by step 1.")
  }
  if (is.null(fit$y_digest)) {
    stop("The fit does not record the responses of step 1; re-estimate it.")
  }
  now <- vapply(index, function(j) .mgcvst_row_digest(Y, j), character(1L))
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
           "estimator or by a version before 0.0.1.9034 (its chunks lack the ",
           "effective degrees of freedom of the spatial smooth); use a new ",
           "checkpoint_dir.")
    }
    if (!identical(old$signature, signature)) {
      stop("The estimation checkpoint ", path, " was written for a different ",
           "model, offset or controls; use a new checkpoint_dir.")
    }
    # Files left by an interrupted write are never part of a checkpoint.
    unlink(list.files(path, "[.]tmp$", full.names = TRUE))
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

# Key of one chunk: the step, its features, their responses and, where a
# feature can be routed to another family (the Poisson prescreen), its routing.
.mgcvst_chunk_key <- function(step, index, y_digest, route = NULL) {
  substr(digest::digest(list(step, as.integer(index), y_digest[index],
                             if (is.null(route)) NULL else route[index]),
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

# TRUE when a chunk result holds a failed feature (a fit error, or in the mgcv
# null step a failed null fit); a resumed run computes such a chunk again.
.mgcvst_chunk_has_error <- function(result) {
  any(vapply(result, function(z) {
    is.list(z) && (!is.null(z[["error"]]) || !is.null(z[["marginal_error"]]))
  }, logical(1L)))
}

# Worker-side wrapper that returns the result of a chunk together with the key
# of its payload. It lives in baseenv(), like the worker functions it calls, so
# that a worker needs nothing from mgcvST to deserialize it.
.mgcvst_keyed_task <- function() {
  task <- function(payload, .worker, ...) {
    list(key = payload$chunk_key, result = .worker(payload, ...))
  }
  environment(task) <- baseenv()
  task
}

# Place the keyed results of bpiterate at the positions of the dispatched
# chunks. `keys` are the keys of the chunks that were dispatched. A missing,
# duplicated or unknown key is an error: a result must never reach the wrong
# chunk.
.mgcvst_place_by_key <- function(returned, keys) {
  got <- vapply(returned, function(z) {
    if (is.list(z) && !inherits(z, "condition") && is.character(z[["key"]]) &&
        length(z[["key"]]) == 1L) z[["key"]] else NA_character_
  }, character(1L))
  if (anyNA(got)) {
    stop("A chunk worker returned a result without the key of its chunk.")
  }
  if (anyDuplicated(got)) {
    stop("Two chunk results carry the same key: ",
         paste(unique(got[duplicated(got)]), collapse = ", "), ".")
  }
  unknown <- setdiff(got, keys)
  if (length(unknown)) {
    stop("A chunk result carries an unknown key: ", unknown[[1L]], ".")
  }
  missing <- setdiff(keys, got)
  if (length(missing)) {
    stop("No result returned for the chunk with key ", missing[[1L]], ".")
  }
  lapply(match(keys, got), function(i) returned[[i]][["result"]])
}

# Read the completed chunks of a step and run the missing ones. `groups` are
# the feature indices of the chunks and `make_payload(index)` builds the payload
# of one chunk; it is called only for a chunk that has to be computed, when a
# worker is free (BiocParallel::bpiterate), so that no more than about one
# payload per worker is alive and a resumed chunk is never copied. `FUN` is a
# worker function receiving a payload (with chunk_file and chunk_key) and `...`;
# a result is a list with one entry per feature of its chunk. A chunk that
# holds a failed feature is computed again. Every dispatched result returns
# with the key of its chunk and is placed by that key, so the outcome does not
# depend on whether bpiterate returns results in iteration or completion order.
.mgcvst_run_chunks <- function(groups, make_payload, step, store, y_digest,
                               route, BPPARAM, FUN, ...) {
  n <- length(groups)
  # Unnamed keys: a chunk is found by its key whatever its position in the list.
  keys <- unname(vapply(groups, function(index) {
    .mgcvst_chunk_key(step, index, y_digest, route)
  }, character(1L)))
  results <- vector("list", n)
  todo <- integer()
  for (k in seq_len(n)) {
    file <- .mgcvst_chunk_file(store, step, keys[k])
    if (!is.na(file) && file.exists(file)) {
      z <- tryCatch(readRDS(file), error = function(e) NULL)
      if (!is.list(z) || !identical(unname(z$key), keys[k]) ||
          length(z$result) != length(groups[[k]])) {
        stop("The estimation checkpoint chunk ", basename(file),
             " is damaged; delete it to recompute the chunk.")
      }
      if (!.mgcvst_chunk_has_error(z$result)) {
        results[[k]] <- z$result
        next
      }
    }
    todo <- c(todo, k)
  }
  built <- 0L
  if (length(todo)) {
    position <- 0L
    iterate <- function() {
      if (position >= length(todo)) return(NULL)
      position <<- position + 1L
      k <- todo[position]
      payload <- make_payload(groups[[k]])
      file <- .mgcvst_chunk_file(store, step, keys[k])
      payload$chunk_file <- if (is.na(file)) NULL else file
      payload$chunk_key <- keys[k]
      built <<- built + 1L
      payload
    }
    if (anyDuplicated(keys[todo])) {
      stop("The chunks of a step must hold distinct features.")
    }
    returned <- BiocParallel::bpiterate(iterate, .mgcvst_keyed_task(), ...,
                                        .worker = FUN, BPPARAM = BPPARAM)
    results[todo] <- .mgcvst_place_by_key(returned, keys[todo])
  }
  list(results = results, resumed = n - length(todo), chunks = n, built = built)
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
