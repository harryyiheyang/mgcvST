#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) {
  stop("Usage: conditional-hpc-run.R FIT_RDS GENE_IDS_FILE G CHECKPOINT_DIR TIMING_CSV")
}
fit_file <- args[1L]
ids_file <- args[2L]
G <- as.integer(args[3L])
checkpoint_dir <- args[4L]
output_file <- args[5L]
if (is.na(G) || G < 2L) stop("G must be an integer of at least 2.")
if (!file.exists(fit_file) || !file.exists(ids_file)) stop("Fit or gene ID file is missing.")
if (dir.exists(checkpoint_dir)) stop("Checkpoint directory must be fresh for a scale run.")

library(mgcvST)
t0 <- proc.time()[["elapsed"]]
fit <- readRDS(fit_file)
ids <- readLines(ids_file, warn = FALSE)
if (length(ids) < G || anyDuplicated(ids) || anyNA(match(ids[seq_len(G)], fit$feature_id))) {
  stop("The gene ID file must contain at least G distinct IDs present in the fit.")
}
ids <- ids[seq_len(G)]
pairs <- t(utils::combn(ids, 2L))
t1 <- proc.time()[["elapsed"]]
result <- inlaST.test(
  fit, pairs = pairs, pairwise_method = "conditional_cauchy",
  conditional_precision = "float32", method = "BY", threads = 40L,
  checkpoint_dir = checkpoint_dir, resume = FALSE
)
t2 <- proc.time()[["elapsed"]]
m <- G * (G - 1) / 2
if (nrow(result$results) != m ||
    !all(c("feature1", "feature2", "signed_score", "p_two_sided", "log_p_two_sided") %in%
         names(result$results))) {
  stop("Conditional result has an unexpected shape or column set.")
}
timing <- result$timing
stages <- c(timing$basis_elapsed, timing$score_unit_elapsed,
            timing$reduced_materialize_elapsed, timing$variance_elapsed)
names(stages) <- c("basis_seconds", "score_seconds", "materialize_seconds",
                   "variance_seconds")
row <- data.frame(
  G = G, pairs = m, rank = timing$score_rank,
  precision = timing$conditional_precision, threads = timing$threads,
  load_and_pairs_seconds = t1 - t0,
  basis_seconds = stages[1L], score_seconds = stages[2L],
  materialize_seconds = stages[3L], variance_seconds = stages[4L],
  remaining_test_seconds = (t2 - t1) - sum(stages),
  test_seconds = t2 - t1, total_seconds = t2 - t0
)
write.csv(row, output_file, row.names = FALSE)
print(row)
