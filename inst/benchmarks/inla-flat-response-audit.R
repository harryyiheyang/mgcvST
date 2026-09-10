#!/usr/bin/env Rscript

# Read-only response and invalid-result audit. No estimator is called.

bam_root <- Sys.getenv(
  "MGCVST_FLAT_BAM_NULL_ROOT",
  "artifacts/flat-prior-investigation/bam-null-500"
)
inla_root <- Sys.getenv(
  "MGCVST_FLAT_INLA_NULL_ROOT",
  "artifacts/flat-prior-investigation/null-500"
)
bam <- read.csv(file.path(bam_root, "replicates.csv"),
                stringsAsFactors = FALSE)
payload <- readRDS(file.path(inla_root, "payload.rds"))
stopifnot(payload$seed_base == 161000L, payload$mean_count == .3,
          payload$rho == 0)

RNGkind("L'Ecuyer-CMRG")
regenerate <- function(i, payload) {
  set.seed(payload$seed_base + i)
  innovation <- matrix(rnorm(ncol(payload$truth) * 2L), ncol = 2L)
  if (payload$rho != 0) {
    innovation[, 2L] <- payload$rho * innovation[, 1L] +
      sqrt(1 - payload$rho^2) * innovation[, 2L]
  }
  signal <- payload$truth %*% innovation
  beta <- log(payload$mean_count) - log(mean(exp(
    payload$offset + 0.5 * rowSums(payload$truth^2)
  )))
  y <- lapply(1:2, function(j) rnbinom(
    nrow(payload$X),
    mu = exp(beta + .25 * (j - 1L) + payload$offset + signal[, j]),
    size = 2
  ))
  y
}

ids <- sort(unique(bam$replicate))
response_rows <- lapply(ids, function(i) {
  y <- regenerate(i, payload)
  observed <- bam[bam$replicate == i, , drop = FALSE][1L, ]
  data.frame(
    replicate = i,
    regenerated_mean_y1 = mean(y[[1L]]),
    regenerated_mean_y2 = mean(y[[2L]]),
    regenerated_checksum1 = sum(seq_along(y[[1L]]) * y[[1L]]),
    regenerated_checksum2 = sum(seq_along(y[[2L]]) * y[[2L]]),
    bam_mean_y1 = observed$mean_y1, bam_mean_y2 = observed$mean_y2,
    bam_checksum1 = observed$y_checksum1,
    bam_checksum2 = observed$y_checksum2,
    exact_bam_match =
      mean(y[[1L]]) == observed$mean_y1 &&
      mean(y[[2L]]) == observed$mean_y2 &&
      sum(seq_along(y[[1L]]) * y[[1L]]) == observed$y_checksum1 &&
      sum(seq_along(y[[2L]]) * y[[2L]]) == observed$y_checksum2
  )
})
response_audit <- do.call(rbind, response_rows)

# The INLA worker saves response vectors for its first two predeclared
# replicates. Verify those vectors exactly, rather than only through checksums.
cache_audit <- lapply(1:2, function(i) {
  path <- file.path(inla_root, "cache",
                    sprintf("rep-%04d-flat_spatial.rds", i))
  if (!file.exists(path)) return(data.frame(
    replicate = i, cache_exists = FALSE, exact_inla_match = NA
  ))
  stored <- readRDS(path)$y
  generated <- regenerate(i, payload)
  data.frame(
    replicate = i, cache_exists = TRUE,
    exact_inla_match = identical(stored, generated)
  )
})
cache_audit <- do.call(rbind, cache_audit)

invalid <- bam[!is.finite(bam$p_value), , drop = FALSE]
invalid$inferred_reason <- ifelse(
  is.finite(invalid$signed_score) & !invalid$fallback &
    is.na(invalid$error),
  "calibrator early return: information non-finite or <=1e-10",
  "other"
)

all_attempt <- do.call(rbind, lapply(split(bam, bam$variant), function(z) {
  valid <- is.finite(z$p_value)
  rejected <- sum(z$p_value[valid] < .05)
  data.frame(
    variant = z$variant[1L], attempted = nrow(z), valid = sum(valid),
    invalid = sum(!valid), rejected = rejected,
    valid_only_rate = rejected / sum(valid),
    all_attempt_lower = rejected / nrow(z),
    all_attempt_upper = (rejected + sum(!valid)) / nrow(z)
  )
}))

audit <- file.path(bam_root, "audit")
dir.create(audit, showWarnings = FALSE)
write.csv(response_audit, file.path(audit, "response-match.csv"),
          row.names = FALSE)
write.csv(cache_audit, file.path(audit, "inla-cache-match.csv"),
          row.names = FALSE)
write.csv(invalid, file.path(audit, "invalid-results.csv"), row.names = FALSE)
write.csv(all_attempt, file.path(audit, "all-attempt-bounds.csv"),
          row.names = FALSE)

cat("BAM regenerated response matches:", sum(response_audit$exact_bam_match),
    "of", nrow(response_audit), "\n")
cat("INLA exact cached response matches:",
    sum(cache_audit$exact_inla_match, na.rm = TRUE), "of",
    sum(cache_audit$cache_exists), "\n")
print(invalid[, c("replicate", "variant", "signed_score", "sp1", "sp2",
                  "inferred_reason")], row.names = FALSE)
print(all_attempt, row.names = FALSE)
