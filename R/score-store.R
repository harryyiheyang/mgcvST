# One score state per feature, in files of a temporary directory, so that the
# pair stage holds only the states its memory allows. The directory is removed
# by .mgcvst_store_cleanup() when the pair stage ends.
.mgcvst_store_open <- function(feature_ids) {
  if (!is.character(feature_ids) || !length(feature_ids) ||
      anyNA(feature_ids) || any(!nzchar(feature_ids)) ||
      anyDuplicated(feature_ids)) {
    stop("feature_ids must be distinct, non-empty character identifiers.")
  }
  structure(list(path = .mgcvst_dense_temp_dir(), feature_ids = feature_ids,
                 storage = "double"),
            class = "mgcvst_score_store")
}

.mgcvst_store_index <- function(store, feature) {
  if (!inherits(store, "mgcvst_score_store")) {
    stop("store must be returned by .mgcvst_store_open().")
  }
  if (is.character(feature) && length(feature) == 1L && !is.na(feature)) {
    index <- match(feature, store$feature_ids)
  } else if (is.numeric(feature) && length(feature) == 1L &&
             is.finite(feature) && feature >= 1 &&
             feature <= length(store$feature_ids) &&
             feature == floor(feature)) {
    index <- as.integer(feature)
  } else {
    stop("feature must be one feature ID or one integer index.")
  }
  if (is.na(index) || index < 1L || index > length(store$feature_ids)) {
    stop("feature is not in this score-state store.")
  }
  index
}

.mgcvst_store_file <- function(store, feature) {
  index <- .mgcvst_store_index(store, feature)
  file.path(store$path, sprintf("feature-%010d.rds", index))
}

.mgcvst_store_read <- function(store, feature) {
  path <- .mgcvst_store_file(store, feature)
  index <- .mgcvst_store_index(store, feature)
  if (!file.exists(path)) {
    stop("The requested score-state shard has not been written: ", store$feature_ids[index], ".")
  }
  unit <- tryCatch(readRDS(path), error = function(e) {
    stop("The score-state shard is unreadable: ", basename(path), ": ",
         conditionMessage(e))
  })
  if (!is.list(unit) || !identical(unit$format, 1L) ||
      !identical(unit$storage, store$storage) ||
      !identical(unit$feature_id,
                 store$feature_ids[.mgcvst_store_index(store, feature)])) {
    stop("The score-state shard has incompatible metadata: ", basename(path), ".")
  }
  if (!is.null(unit$error)) {
    if (!is.character(unit$error) || length(unit$error) != 1L ||
        is.na(unit$error) || !nzchar(unit$error)) {
      stop("The score-state shard has an invalid error record: ", basename(path), ".")
    }
    return(list(error = unit$error))
  }
  state <- tryCatch(.mgcvst_unpack_score_state(unit$state), error = function(e) {
    stop("The score-state shard is damaged: ", basename(path), ": ",
         conditionMessage(e))
  })
  if (!is.numeric(state$a) || length(state$a) != nrow(state$M) ||
      any(!is.finite(state$a)) || any(!is.finite(state$M))) {
    stop("The score-state shard has invalid numeric data: ", basename(path), ".")
  }
  state
}

.mgcvst_store_write <- function(store, feature, state) {
  path <- .mgcvst_store_file(store, feature)
  if (file.exists(path)) {
    stop("The score-state shard already exists: ", basename(path), ".")
  }
  if (!is.list(state)) stop("state must be a score-state list.")
  if (!is.null(state$error)) {
    if (!is.character(state$error) || length(state$error) != 1L ||
        is.na(state$error) || !nzchar(state$error)) {
      stop("state$error must be one non-empty error message.")
    }
    payload <- NULL
    error <- state$error
  } else {
    a <- as.numeric(state$a)
    M <- as.matrix(state$M)
    if (!is.numeric(state$a) || !is.numeric(state$M) ||
        !is.matrix(M) || nrow(M) < 1L || nrow(M) != ncol(M) ||
        length(a) != nrow(M) || any(!is.finite(a)) ||
        any(!is.finite(M)) ||
        !isTRUE(all.equal(M, t(M), tolerance = 1e-8,
                          check.attributes = FALSE))) {
      stop("state must contain a finite score vector and aligned symmetric matrix.")
    }
    payload <- .mgcvst_pack_score_state(state)
    error <- NULL
  }
  unit <- list(format = 1L, storage = store$storage,
               feature_id = store$feature_ids[.mgcvst_store_index(store, feature)],
               error = error, state = payload)
  temp <- tempfile(paste0(basename(path), "-"), tmpdir = store$path,
                   fileext = ".tmp")
  on.exit(if (file.exists(temp)) unlink(temp), add = TRUE)
  saveRDS(unit, temp, compress = FALSE)
  if (file.exists(path) || !file.rename(temp, path)) {
    stop("Could not commit the score-state shard: ", basename(path), ".")
  }
  invisible(path)
}

.mgcvst_store_cleanup <- function(store) {
  if (!inherits(store, "mgcvst_score_store")) {
    stop("store must be returned by .mgcvst_store_open().")
  }
  .mgcvst_dense_cleanup(store$path)
  invisible(TRUE)
}
