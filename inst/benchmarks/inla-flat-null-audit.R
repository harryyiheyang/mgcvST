#!/usr/bin/env Rscript

# Read-only integrity audit for a running or completed flat-prior Monte Carlo.
# It takes a fixed file snapshot at startup and never fits or rewrites results.

root <- Sys.getenv("MGCVST_FLAT_AUDIT_ROOT",
                   "artifacts/flat-prior-investigation/null-500")
completed <- file.path(root, "completed")
stopifnot(dir.exists(completed))
files <- sort(list.files(completed, pattern = "^rep-[0-9]+[.]rds$",
                         full.names = TRUE))
stopifnot(length(files) > 0L)

objects <- lapply(files, function(path) {
  tryCatch(readRDS(path), error = function(e) structure(
    list(path = path, message = conditionMessage(e)), class = "read_error"
  ))
})
read_ok <- !vapply(objects, inherits, logical(1L), "read_error")
read_errors <- if (all(read_ok)) data.frame() else do.call(rbind, lapply(
  objects[!read_ok], function(x) data.frame(file = basename(x$path),
                                             error = x$message)
))
objects <- objects[read_ok]
files_ok <- files[read_ok]

valid_object <- vapply(objects, function(x) {
  is.list(x) && is.data.frame(x$rows) && is.data.frame(x$metadata)
}, logical(1L))
if (!all(valid_object)) stop("Readable result file lacks rows/metadata data frames: ",
                             paste(basename(files_ok[!valid_object]),
                                   collapse = ", "))
rows <- do.call(rbind, lapply(objects, `[[`, "rows"))
metadata <- do.call(rbind, lapply(objects, `[[`, "metadata"))

file_replicate <- as.integer(sub("^rep-([0-9]+)[.]rds$", "\\1",
                                 basename(files_ok)))
object_replicate <- vapply(objects, function(x) {
  ids <- unique(c(x$rows$replicate, x$metadata$replicate))
  if (length(ids) == 1L) as.integer(ids) else NA_integer_
}, integer(1L))

expected_variants <- c("flat_spatial", "flat_both")
expected_kernels <- c("conditioned", "raw_centered")
completeness <- do.call(rbind, Map(function(x, file_id, object_id, file) {
  data.frame(
    file = basename(file), file_replicate = file_id,
    object_replicate = object_id,
    id_match = identical(file_id, object_id),
    row_count = nrow(x$rows), metadata_count = nrow(x$metadata),
    row_keys_complete = setequal(
      paste(x$rows$variant, x$rows$kernel),
      as.vector(outer(expected_variants, expected_kernels, paste))
    ),
    metadata_keys_complete = setequal(
      paste(x$metadata$variant, x$metadata$feature),
      as.vector(outer(expected_variants, 1:2, paste))
    )
  )
}, objects, file_replicate, object_replicate, files_ok))

qfun <- function(x, probability) {
  if (!length(x) || !any(is.finite(x))) return(NA_real_)
  unname(quantile(x[is.finite(x)], probability, names = FALSE))
}
variant_audit <- do.call(rbind, lapply(expected_variants, function(variant) {
  m <- metadata[metadata$variant == variant, , drop = FALSE]
  z <- rows[rows$variant == variant, , drop = FALSE]
  data.frame(
    variant = variant, snapshot_replicates = length(objects),
    metadata_rows = nrow(m), score_rows = nrow(z),
    converged_false = sum(!m$converged, na.rm = TRUE),
    nonzero_mode_status = sum(m$mode_status != 0L, na.rm = TRUE),
    nonfinite_tau = sum(!is.finite(m$tau)),
    tau_min = qfun(m$tau, 0), tau_median = qfun(m$tau, .5),
    tau_p99 = qfun(m$tau, .99), tau_max = qfun(m$tau, 1),
    tau_below_1e_minus8 = sum(m$tau < 1e-8, na.rm = TRUE),
    tau_above_1e8 = sum(m$tau > 1e8, na.rm = TRUE),
    nonfinite_nb_size = sum(!is.finite(m$nb_size)),
    nb_size_min = qfun(m$nb_size, 0), nb_size_median = qfun(m$nb_size, .5),
    nb_size_p95 = qfun(m$nb_size, .95), nb_size_p99 = qfun(m$nb_size, .99),
    nb_size_max = qfun(m$nb_size, 1),
    nb_size_above_1e6 = sum(m$nb_size > 1e6, na.rm = TRUE),
    nb_size_above_1e12 = sum(m$nb_size > 1e12, na.rm = TRUE),
    mean_constraint_above_1e_minus8 = sum(m$mean_error > 1e-8,
                                          na.rm = TRUE),
    invalid_scores = sum(!is.finite(z$p_value)),
    davies_fallbacks = sum(z$fallback, na.rm = TRUE),
    information_min = qfun(z$information, 0),
    information_p01 = qfun(z$information, .01),
    information_median = qfun(z$information, .5),
    maximum_P1_max = qfun(z$maximum_P1, 1)
  )
}))

invalid_reasons <- rows[!is.finite(rows$p_value),
                        c("replicate", "variant", "kernel", "information",
                          "fallback", "error"), drop = FALSE]
if (nrow(invalid_reasons)) {
  invalid_reason_counts <- as.data.frame(table(
    variant = invalid_reasons$variant,
    kernel = invalid_reasons$kernel,
    error = invalid_reasons$error, useNA = "ifany"
  ), stringsAsFactors = FALSE)
  invalid_reason_counts <- invalid_reason_counts[
    invalid_reason_counts$Freq > 0L, , drop = FALSE
  ]
} else {
  invalid_reason_counts <- data.frame()
}

audit <- file.path(root, "audit")
dir.create(audit, showWarnings = FALSE)
tag <- sprintf("snapshot-%04d", length(files))
write.csv(variant_audit, file.path(audit, paste0(tag, "-variants.csv")),
          row.names = FALSE)
write.csv(completeness, file.path(audit, paste0(tag, "-completeness.csv")),
          row.names = FALSE)
write.csv(invalid_reasons, file.path(audit, paste0(tag, "-invalid.csv")),
          row.names = FALSE)
write.csv(invalid_reason_counts,
          file.path(audit, paste0(tag, "-invalid-reason-counts.csv")),
          row.names = FALSE)
write.csv(read_errors, file.path(audit, paste0(tag, "-read-errors.csv")),
          row.names = FALSE)

cat("Snapshot files:", length(files), "readable:", length(objects), "\n")
cat("Complete file schemas:",
    sum(completeness$id_match & completeness$row_keys_complete &
          completeness$metadata_keys_complete), "of", nrow(completeness), "\n")
print(variant_audit, row.names = FALSE)
if (nrow(invalid_reason_counts)) print(invalid_reason_counts, row.names = FALSE)
