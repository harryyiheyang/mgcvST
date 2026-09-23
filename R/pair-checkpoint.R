# Completed pair batches are separate from reusable per-feature score states.
.mgcvst_pair_checkpoint <- function(store, index, pair_index, calibration = NULL) {
  if (isTRUE(store$temporary)) return(NULL)
  inputs <- list(version = 1L, index = index, pair_index = pair_index)
  if (!is.null(calibration)) inputs$calibration <- calibration
  signature <- digest::digest(inputs, algo = "sha256")
  path <- file.path(store$path, paste0("pairs-", signature))
  if (!dir.exists(path) && !dir.create(path)) {
    stop("Could not create the pair checkpoint directory.")
  }
  path
}

.mgcvst_pair_checkpoint_read <- function(path, first, pair_index) {
  if (is.null(path)) return(NULL)
  file <- file.path(path, sprintf("block-%010d.rds", first))
  if (!file.exists(file)) return(NULL)
  z <- readRDS(file)
  columns <- c("pair_index", "score", "information", "effective_rank",
               "p_value", "error_message")
  if (!is.list(z) || !identical(z$first, first) ||
      length(z$last) != 1L || !is.finite(z$last) ||
      z$last < first || z$last > length(pair_index) ||
      z$last != floor(z$last) || !is.data.frame(z$result) ||
      !identical(names(z$result), columns) ||
      !identical(z$result$pair_index, unname(pair_index[seq.int(first, z$last)])) ||
      !identical(z$checksum, digest::digest(z$result, algo = "sha256"))) {
    stop("The pair checkpoint is damaged or incompatible: ", basename(file), ".")
  }
  z
}

.mgcvst_pair_checkpoint_write <- function(path, first, last, result) {
  if (is.null(path)) return(invisible(NULL))
  file <- file.path(path, sprintf("block-%010d.rds", first))
  if (file.exists(file)) stop("The pair checkpoint already exists.")
  tmp <- tempfile("pair-", tmpdir = path, fileext = ".tmp")
  on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
  saveRDS(list(first = first, last = last, result = result,
              checksum = digest::digest(result, algo = "sha256")),
          tmp, compress = FALSE)
  if (!file.rename(tmp, file)) stop("Could not commit the pair checkpoint.")
  invisible(NULL)
}

# Reuse two resident gene blocks before loading the next pair of blocks.
.mgcvst_pair_order <- function(index, used, capacity, path = NULL) {
  file <- if (!is.null(path)) file.path(path, "schedule.rds") else NULL
  if (!is.null(file) && file.exists(file)) {
    z <- readRDS(file)
    if (!is.list(z) || length(z$order) != nrow(index) ||
        !identical(z$checksum, digest::digest(z$order, algo = "sha256")) ||
        !identical(sort(z$order), seq_len(nrow(index)))) {
      stop("The stored pair schedule is damaged or incompatible.")
    }
    return(z$order)
  }
  width <- max(1L, floor(capacity / 2))
  block1 <- (match(index[, 1L], used) - 1L) %/% width
  block2 <- (match(index[, 2L], used) - 1L) %/% width
  ord <- if (capacity >= length(used)) seq_len(nrow(index)) else
    order(pmin(block1, block2), pmax(block1, block2), method = "radix")
  if (!is.null(file)) {
    tmp <- tempfile("schedule-", tmpdir = path, fileext = ".tmp")
    on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
    saveRDS(list(order = ord, checksum = digest::digest(ord, algo = "sha256")),
            tmp, compress = FALSE)
    if (!file.rename(tmp, file)) stop("Could not commit the pair schedule.")
  }
  ord
}
