# Pack one symmetric score matrix without retaining a duplicated triangle.
.mgcvst_pack_symmetric <- function(M) {
  M <- as.matrix(M)
  if (nrow(M) != ncol(M)) stop("A score-state matrix must be square.")
  list(n = nrow(M), upper = M[upper.tri(M, diag = TRUE)])
}

.mgcvst_unpack_symmetric <- function(x) {
  n <- as.integer(x$n)
  M <- matrix(0, n, n)
  keep <- upper.tri(M, diag = TRUE)
  if (length(x$upper) != sum(keep)) stop("Packed score state has an invalid size.")
  M[keep] <- x$upper
  M[lower.tri(M)] <- t(M)[lower.tri(M)]
  M
}

.mgcvst_pack_score_state <- function(state) {
  list(a = as.numeric(state$a), M = .mgcvst_pack_symmetric(state$M),
       width = state$width)
}

.mgcvst_unpack_score_state <- function(unit) {
  list(a = unit$a, M = .mgcvst_unpack_symmetric(unit$M), width = unit$width)
}

# Create and remove only a verified child of the R session temporary directory.
.mgcvst_dense_temp_dir <- function() {
  path <- tempfile("mgcvst-dense-cache-", tmpdir = tempdir())
  if (!dir.create(path)) stop("Could not create the dense score cache directory.")
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

.mgcvst_dense_cleanup <- function(path) {
  root <- paste0(normalizePath(tempdir(), winslash = "/", mustWork = TRUE), "/")
  target <- paste0(normalizePath(path, winslash = "/", mustWork = TRUE), "/")
  if (!startsWith(target, root) || identical(target, root)) {
    stop("Refusing to remove an unverified dense score cache directory.")
  }
  unlink(sub("/$", "", target), recursive = TRUE, force = TRUE)
}

# Bound both pair count and simultaneously restored dense feature states.
.mgcvst_dense_pair_groups <- function(rows, index, chunk_size,
                                      max_features = 32L) {
  out <- list()
  current <- integer()
  used <- integer()
  for (row in rows) {
    next_used <- union(used, index[row, ])
    if (length(current) &&
        (length(current) >= chunk_size || length(next_used) > max_features)) {
      out[[length(out) + 1L]] <- current
      current <- integer()
      used <- integer()
    }
    current <- c(current, row)
    used <- union(used, index[row, ])
  }
  if (length(current)) out[[length(out) + 1L]] <- current
  out
}

# Read and unpack each feature needed by one pair task exactly once.
.mgcvst_dense_pair_chunk <- function(payload, calibration) {
  .mgcvst_thread_limit()
  states <- list()
  for (key in names(payload$shards)) {
    unit <- readRDS(payload$shards[[key]])
    states[[key]] <- if (is.null(unit$error)) .mgcvst_unpack_score_state(unit) else unit
  }
  out <- vector("list", length(payload$rows))
  for (k in seq_along(payload$rows)) {
    keys <- as.character(payload$pairs[k, ])
    u1 <- states[[keys[1L]]]
    u2 <- states[[keys[2L]]]
    if (is.null(u1) || is.null(u2)) stop("A pair task is missing a feature score state.")
    if (!is.null(u1$error) || !is.null(u2$error)) {
      error <- c(u1$error, u2$error)
      out[[k]] <- data.frame(
        pair_index = payload$rows[k], signed_score = NA_real_,
        information = NA_real_, effective_rank = NA_real_,
        p_two_sided = NA_real_, p_positive = NA_real_, p_negative = NA_real_,
        error_message = paste(error[!is.na(error)], collapse = " | "),
        stringsAsFactors = FALSE
      )
      next
    }
    z <- tryCatch({
      score <- as.numeric(crossprod(u1$a, u2$a))
      cal <- rkhs_score_calibrate(score, u1$M, u2$M, method = calibration)
      list(score = score, cal = cal)
    }, error = function(e) e)
    if (inherits(z, "condition")) {
      out[[k]] <- data.frame(
        pair_index = payload$rows[k], signed_score = NA_real_,
        information = NA_real_, effective_rank = NA_real_,
        p_two_sided = NA_real_, p_positive = NA_real_, p_negative = NA_real_,
        error_message = conditionMessage(z), stringsAsFactors = FALSE
      )
      next
    }
    out[[k]] <- data.frame(
      pair_index = payload$rows[k], signed_score = z$score,
      information = z$cal$information, effective_rank = z$cal$effective_rank,
      p_two_sided = z$cal$p_two_sided, p_positive = z$cal$p_positive,
      p_negative = z$cal$p_negative, error_message = NA_character_,
      stringsAsFactors = FALSE
    )
  }
  out
}
