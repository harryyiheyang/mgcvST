# The two-step estimation.
#
# Step 1 fits only the null model of every feature and returns the Stage 1
# null-first p-value; step 2 fits the spatial model of the selected features
# only. This file holds what both branches share: the Stage 1 q-values, the
# feature selection, the content digests that tie a later step to the Y of
# step 1, and the keyed dispatch of chunks to the workers.

# Format of the fit objects of the two-step estimators. Fits without it were
# estimated before 0.0.1.9032.
.mgcvst_fit_format <- 2L

# A fit estimated before 0.0.1.9032 lacks the stored mean mu_bar and the
# two-step bookkeeping that the tests and the later steps read, so it is
# refused rather than reused.
.mgcvst_check_fit_format <- function(fit) {
  format <- fit$format
  if (!is.null(format) && (!is.numeric(format) || length(format) != 1L ||
                           format > .mgcvst_fit_format)) {
    stop("The fit was written by a newer version of spdeST; update the package.",
         call. = FALSE)
  }
  if (is.null(format) || format < .mgcvst_fit_format || !is.numeric(fit$mu_bar)) {
    stop("This fit was estimated before 0.0.1.9032 and lacks the stored mean ",
         "mu_bar that the pair test reads; re-run spdeST.estimate().",
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

# ---- keyed chunks ------------------------------------------------------------

# Key of one chunk: the step, its features, their responses and, where a
# feature can be routed to another family (the Poisson prescreen), its routing.
.mgcvst_chunk_key <- function(step, index, y_digest, route = NULL) {
  substr(digest::digest(list(step, as.integer(index), y_digest[index],
                             if (is.null(route)) NULL else route[index]),
                        algo = "sha256"), 1L, 24L)
}

# Worker-side wrapper that returns the result of a chunk together with the key
# of its payload. It lives in baseenv(), like the worker functions it calls, so
# that a worker needs nothing from spdeST to deserialize it.
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

# Run the chunks of a step. `groups` are the feature indices of the chunks and
# `make_payload(index)` builds the payload of one chunk; it is called when a
# worker is free (BiocParallel::bpiterate), so that no more than about one
# payload per worker is alive. `FUN` is a worker function receiving a payload
# (with chunk_key) and `...`; a result is a list with one entry per feature of
# its chunk. Every result returns with the key of its chunk and is placed by
# that key, so the outcome does not depend on whether bpiterate returns
# results in iteration or completion order. Nothing is written to disk.
.mgcvst_run_chunks <- function(groups, make_payload, step, y_digest,
                               route, BPPARAM, FUN, ...) {
  n <- length(groups)
  keys <- unname(vapply(groups, function(index) {
    .mgcvst_chunk_key(step, index, y_digest, route)
  }, character(1L)))
  results <- vector("list", n)
  built <- 0L
  progress <- .mgcvst_progress(step, sum(lengths(groups)), "features")
  if (n) {
    position <- 0L
    iterate <- function() {
      if (position >= n) return(NULL)
      position <<- position + 1L
      payload <- make_payload(groups[[position]])
      payload$chunk_key <- keys[position]
      built <<- built + 1L
      payload
    }
    if (anyDuplicated(keys)) {
      stop("The chunks of a step must hold distinct features.")
    }
    # The manager collects each result as it arrives and reports progress.
    collect <- function(acc, x) {
      acc[[length(acc) + 1L]] <- x
      progress$add(if (is.list(x) && !inherits(x, "condition"))
        length(x[["result"]]) else 0L)
      acc
    }
    returned <- BiocParallel::bpiterate(iterate, .mgcvst_keyed_task(), ...,
                                        .worker = FUN, REDUCE = collect,
                                        init = list(), BPPARAM = BPPARAM)
    results <- .mgcvst_place_by_key(returned, keys)
  }
  progress$done()
  list(results = results, chunks = n, built = built)
}

# Contiguous feature chunks of at most `chunk_size` features. The default is
# one chunk per worker.
.mgcvst_feature_chunks <- function(index, chunk_size, BPPARAM) {
  if (!length(index)) return(list())
  if (is.null(chunk_size)) {
    workers <- max(1L, min(length(index), BiocParallel::bpworkers(BPPARAM)))
    chunk_size <- ceiling(length(index) / workers)
  }
  split(index, ceiling(seq_along(index) / chunk_size))
}
