# Store raw sparse reconstruction units separately from materialized score states.
.mgcvst_unit_store_open <- function(path, signature, feature_ids,
                                    resume = TRUE,
                                    kind = c("sparse", "dense")) {
  kind <- match.arg(kind)
  if (!is.character(path) || length(path) != 1L || is.na(path) ||
      !nzchar(path)) {
    stop("path must be one non-empty unit-store directory.")
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
  if (file.exists(path) && !dir.exists(path)) {
    stop("The unit-store path exists and is not a directory.")
  }
  if (!dir.exists(path) && !dir.create(path, recursive = TRUE)) {
    stop("Could not create the sparse unit-store directory.")
  }
  path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  manifest_path <- file.path(path, "manifest.rds")
  manifest <- list(schema = "mgcvst_unit_store", version = 1L,
                   kind = kind, signature = signature,
                   feature_ids = feature_ids)
  if (file.exists(manifest_path)) {
    if (!resume) stop("The sparse unit store already exists; use resume = TRUE.")
    saved <- tryCatch(readRDS(manifest_path), error = function(e) {
      stop("The sparse unit-store manifest is unreadable: ",
           conditionMessage(e))
    })
    if (!identical(saved, manifest)) {
      stop("The sparse unit-store signature, version, or feature IDs do not match.")
    }
  } else {
    contents <- list.files(path, all.files = TRUE, no.. = TRUE)
    if (length(contents)) {
      stop("The sparse unit-store directory lacks a manifest and is not empty.")
    }
    temp <- tempfile("manifest-", tmpdir = path, fileext = ".tmp")
    on.exit(if (file.exists(temp)) unlink(temp), add = TRUE)
    saveRDS(manifest, temp, compress = FALSE)
    if (file.exists(manifest_path) || !file.rename(temp, manifest_path)) {
      stop("Could not commit the sparse unit-store manifest.")
    }
  }
  structure(list(path = path, signature = signature,
                 feature_ids = feature_ids, kind = kind, version = 1L),
            class = "mgcvst_unit_store")
}

.mgcvst_unit_store_index <- function(store, feature) {
  if (!inherits(store, "mgcvst_unit_store")) {
    stop("store must be returned by .mgcvst_unit_store_open().")
  }
  if (is.character(feature) && length(feature) == 1L && !is.na(feature)) {
    index <- match(feature, store$feature_ids)
  } else if (is.numeric(feature) && length(feature) == 1L &&
             is.finite(feature) && feature >= 1 &&
             feature <= length(store$feature_ids) && feature == floor(feature)) {
    index <- as.integer(feature)
  } else {
    stop("feature must be one feature ID or one integer index.")
  }
  if (is.na(index) || index < 1L || index > length(store$feature_ids)) {
    stop("feature is not in this unit store.")
  }
  index
}

.mgcvst_unit_store_file <- function(store, feature) {
  index <- .mgcvst_unit_store_index(store, feature)
  file.path(store$path, sprintf("unit-%010d.rds", index))
}

.mgcvst_unit_store_has <- function(store, feature) {
  file.exists(.mgcvst_unit_store_file(store, feature))
}

.mgcvst_unit_store_validate <- function(unit) {
  if (!is.list(unit)) stop("unit payload must be a list.")
  if (!is.null(unit$error)) {
    if (!is.character(unit$error) || length(unit$error) != 1L ||
        is.na(unit$error) || !nzchar(unit$error)) {
      stop("unit$error must be one non-empty error message.")
    }
    return(invisible(TRUE))
  }
  if (!length(unit)) stop("unit payload must not be empty.")
  invisible(TRUE)
}

.mgcvst_unit_store_write <- function(store, feature, unit) {
  path <- .mgcvst_unit_store_file(store, feature)
  if (file.exists(path)) stop("The sparse unit shard already exists: ", basename(path), ".")
  .mgcvst_unit_store_validate(unit)
  index <- .mgcvst_unit_store_index(store, feature)
  payload <- list(schema = "mgcvst_reconstruction_unit", version = 1L,
                  kind = store$kind,
                  signature = store$signature,
                  feature_id = store$feature_ids[index], unit = unit)
  temp <- tempfile(paste0(basename(path), "-"), tmpdir = store$path,
                   fileext = ".tmp")
  on.exit(if (file.exists(temp)) unlink(temp), add = TRUE)
  saveRDS(payload, temp, compress = FALSE)
  if (file.exists(path) || !file.rename(temp, path)) {
    stop("Could not commit the sparse unit shard: ", basename(path), ".")
  }
  invisible(path)
}

.mgcvst_unit_store_read <- function(store, feature) {
  path <- .mgcvst_unit_store_file(store, feature)
  if (!file.exists(path)) {
    stop("The requested sparse unit shard has not been written: ",
         store$feature_ids[.mgcvst_unit_store_index(store, feature)], ".")
  }
  saved <- tryCatch(readRDS(path), error = function(e) {
    stop("The sparse unit shard is unreadable: ", basename(path), ": ",
         conditionMessage(e))
  })
  index <- .mgcvst_unit_store_index(store, feature)
  if (!is.list(saved) || !identical(saved$schema,
      "mgcvst_reconstruction_unit") ||
      !identical(saved$version, 1L) ||
      !identical(saved$kind, store$kind) ||
      !identical(saved$signature, store$signature) ||
      !identical(saved$feature_id, store$feature_ids[index])) {
    stop("The sparse unit shard has incompatible metadata: ", basename(path), ".")
  }
  tryCatch(.mgcvst_unit_store_validate(saved$unit), error = function(e) {
    stop("The sparse unit shard is damaged: ", basename(path), ": ",
         conditionMessage(e))
  })
  saved$unit
}
