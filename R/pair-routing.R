# The route of a pair test is chosen by the user: `moments = "exact"` or
# `moments = "pcalearning"`. A checkpoint directory records the route it was
# written with, so that a resumed call with another route stops with a clear
# message instead of mixing results of two calibrations.

# The route record of a checkpoint directory, or NULL when it holds none.
.mgcvst_route_stored <- function(checkpoint_dir) {
  file <- file.path(checkpoint_dir, "route.rds")
  if (!file.exists(file)) return(NULL)
  z <- tryCatch(readRDS(file), error = function(e) NULL)
  if (is.list(z) && is.character(z$moments) && length(z$moments) == 1L &&
      z$moments %in% c("exact", "pcalearning")) z else NULL
}

# Stop when a resumed run asks for another route than the directory records.
.mgcvst_route_check <- function(checkpoint_dir, resume, moments) {
  if (is.null(checkpoint_dir) || !isTRUE(resume)) return(invisible(NULL))
  stored <- .mgcvst_route_stored(checkpoint_dir)
  if (!is.null(stored) && !identical(stored$moments, moments)) {
    stop("The checkpoint directory ", checkpoint_dir, " was written with moments = \"",
         stored$moments, "\", and this call asks for moments = \"", moments,
         "\". Resume with moments = \"", stored$moments, "\" or use a new ",
         "checkpoint_dir.", call. = FALSE)
  }
  invisible(NULL)
}

.mgcvst_route_save <- function(checkpoint_dir, route) {
  if (is.null(route)) return(invisible(NULL))
  file <- file.path(checkpoint_dir, "route.rds")
  if (file.exists(file)) return(invisible(NULL))
  tmp <- tempfile("route-", tmpdir = checkpoint_dir, fileext = ".tmp")
  on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
  saveRDS(route[c("moments", "k", "q")], tmp, compress = FALSE)
  if (!file.rename(tmp, file) && !file.exists(file)) {
    stop("Could not commit the route record.")
  }
  invisible(NULL)
}
