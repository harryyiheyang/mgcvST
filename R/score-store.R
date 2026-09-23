# Store one score state per feature so completed states survive interrupted tests.
.mgcvst_store_open <- function(path = NULL, signature, feature_ids,
                               storage = c("double", "float32"), resume = TRUE,
                               encoding = c("rds", "native")) {
  storage <- match.arg(storage)
  encoding <- match.arg(encoding)
  if (encoding == "native" && storage != "double") {
    stop("Native score-state encoding requires storage = 'double'.")
  }
  if (!is.list(signature)) stop("signature must be a list of stable fit identifiers.")
  if (!is.character(feature_ids) || !length(feature_ids) ||
      anyNA(feature_ids) || any(!nzchar(feature_ids)) ||
      anyDuplicated(feature_ids)) {
    stop("feature_ids must be distinct, non-empty character identifiers.")
  }
  if (!is.logical(resume) || length(resume) != 1L || is.na(resume)) {
    stop("resume must be TRUE or FALSE.")
  }

  temporary <- is.null(path)
  if (temporary) {
    path <- .mgcvst_dense_temp_dir()
    committed <- FALSE
    on.exit(if (!committed && dir.exists(path))
      .mgcvst_dense_cleanup(path), add = TRUE)
  } else {
    if (!is.character(path) || length(path) != 1L || is.na(path) ||
        !nzchar(path)) stop("path must be one non-empty directory name.")
    if (file.exists(path) && !dir.exists(path)) {
      stop("The score-state store path exists and is not a directory.")
    }
    if (!dir.exists(path) && !dir.create(path, recursive = TRUE)) {
      stop("Could not create the score-state store directory.")
    }
    path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  }

  manifest_path <- file.path(path, "manifest.rds")
  requested <- if (encoding == "rds") {
    list(format = 1L, signature = signature,
         feature_ids = feature_ids, storage = storage)
  } else {
    list(format = 2L, encoding = "native", signature = signature,
         feature_ids = feature_ids, storage = storage)
  }
  if (file.exists(manifest_path)) {
    if (!resume) stop("The score-state store already exists; use resume = TRUE.")
    existing <- tryCatch(readRDS(manifest_path), error = function(e) {
      stop("The score-state store manifest is unreadable: ", conditionMessage(e))
    })
    if (!identical(existing, requested)) {
      stop("The score-state store signature, feature IDs, or storage do not match.")
    }
  } else {
    contents <- list.files(path, all.files = TRUE, no.. = TRUE)
    if (length(contents)) {
      stop("The score-state store directory lacks a manifest and is not empty.")
    }
    temp <- tempfile("manifest-", tmpdir = path, fileext = ".tmp")
    on.exit(if (file.exists(temp)) unlink(temp), add = TRUE)
    saveRDS(requested, temp, compress = FALSE)
    if (!file.rename(temp, manifest_path)) {
      stop("Could not commit the score-state store manifest.")
    }
  }
  if (temporary) committed <- TRUE
  structure(list(path = path, feature_ids = feature_ids, storage = storage,
                 signature = signature, temporary = temporary,
                 encoding = encoding),
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
  extension <- if (identical(store$encoding, "native")) ".bin" else ".rds"
  file.path(store$path, paste0(sprintf("feature-%010d", index), extension))
}

.mgcvst_store_has <- function(store, feature) {
  file.exists(.mgcvst_store_file(store, feature))
}

.mgcvst_store_read <- function(store, feature) {
  path <- .mgcvst_store_file(store, feature)
  index <- .mgcvst_store_index(store, feature)
  if (!file.exists(path)) {
    stop("The requested score-state shard has not been written: ", store$feature_ids[index], ".")
  }
  if (identical(store$encoding, "native")) {
    state <- tryCatch(mgcvst_state_read_cpp(
      path, digest::digest(store$signature, algo = "sha256"),
      store$feature_ids[index]
    ), error = function(e) {
      stop("The native score-state shard is unreadable: ", basename(path),
           ": ", conditionMessage(e))
    })
    if (!is.list(state)) {
      stop("The native score-state shard returned an invalid record: ", basename(path), ".")
    }
    if (!is.null(state$error) && length(state$error) == 1L &&
        !is.na(state$error) && nzchar(state$error)) {
      return(list(error = state$error))
    }
    if (!is.numeric(state$a) || !is.matrix(state$M) ||
        length(state$a) != nrow(state$M) || nrow(state$M) != ncol(state$M)) {
      stop("The native score-state shard has invalid numeric data: ", basename(path), ".")
    }
    return(list(a = as.numeric(state$a), M = state$M, width = state$width))
  }
  unit <- tryCatch(readRDS(path), error = function(e) {
    stop("The score-state shard is unreadable: ", basename(path), ": ",
         conditionMessage(e))
  })
  if (!is.list(unit) || !identical(unit$format, 1L) ||
      !identical(unit$storage, store$storage) ||
      !identical(unit$signature, store$signature) ||
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
  state <- tryCatch({
    if (identical(store$storage, "float32")) {
      n <- unit$state$M$n
      raw_upper <- unit$state$M$upper_raw
      if (!is.integer(n) || length(n) != 1L || is.na(n) || n < 1L ||
          !is.raw(raw_upper) ||
          length(raw_upper) != 4 * (as.double(n) * (n + 1) / 2)) {
        stop("invalid float32 matrix payload")
      }
      unit$state$M$upper <- readBin(raw_upper, what = "numeric",
                                     n = n * (n + 1) / 2, size = 4L,
                                     endian = "little")
      unit$state$M$upper_raw <- NULL
    }
    .mgcvst_unpack_score_state(unit$state)
  }, error = function(e) {
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
  if (identical(store$encoding, "native")) {
    if (!is.null(state$error)) {
      if (!is.character(state$error) || length(state$error) != 1L ||
          is.na(state$error) || !nzchar(state$error)) {
        stop("state$error must be one non-empty error message.")
      }
      a <- numeric()
      M <- matrix(numeric(), 0L, 0L)
      width <- integer()
      error <- state$error
    } else {
      a <- as.numeric(state$a)
      M <- as.matrix(state$M)
      if (!is.numeric(state$a) || !is.numeric(state$M) || !is.matrix(M) ||
          nrow(M) < 1L || nrow(M) != ncol(M) || length(a) != nrow(M)) {
        stop("state must contain a finite score vector and aligned symmetric matrix.")
      }
      width_names <- names(state$width)
      width <- as.integer(state$width)
      if (!is.null(width_names)) names(width) <- width_names
      error <- ""
    }
    mgcvst_state_write_cpp(
      path, digest::digest(store$signature, algo = "sha256"),
      store$feature_ids[.mgcvst_store_index(store, feature)],
      a, M, width, error
    )
    if (!file.exists(path)) {
      stop("The native score-state writer did not create the shard: ", basename(path), ".")
    }
    return(invisible(path))
  }
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
    if (identical(store$storage, "float32")) {
      payload$M$upper_raw <- writeBin(payload$M$upper, raw(), size = 4L,
                                      endian = "little")
      payload$M$upper <- NULL
    }
    error <- NULL
  }
  unit <- list(format = 1L, storage = store$storage,
               signature = store$signature,
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
  if (!isTRUE(store$temporary)) return(invisible(FALSE))
  .mgcvst_dense_cleanup(store$path)
  invisible(TRUE)
}
