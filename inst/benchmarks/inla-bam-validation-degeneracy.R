Sys.setenv(LC_ALL = "C")

repo <- normalizePath(file.path(getwd()), winslash = "/", mustWork = TRUE)
lib <- normalizePath(Sys.getenv("MGCVST_VALIDATION_LIBRARY",
  "artifacts/inla-bam-validation/library"), mustWork = TRUE)
.libPaths(c(lib, .libPaths()))
library(mgcvST)

summary_dir <- file.path(repo, "artifacts/inla-bam-validation/summary")
pair_file <- file.path(summary_dir, "inference-paired-results.csv")
D <- read.csv(pair_file, stringsAsFactors = FALSE)

bad <- D[!D$valid_bam, ]
key <- paste(bad$case, bad$replicate, bad$feature1, bad$feature2, sep = "|")
stopifnot(nrow(D) == 600L)
stopifnot(nrow(bad) == 34L)
stopifnot(length(unique(key)) == 17L)
stopifnot(all(table(bad$calibration) == 17L))
stopifnot(all(is.finite(bad$information_bam)))
stopifnot(all(bad$information_bam > 0 & bad$information_bam <= 1e-10))
stopifnot(all(bad$effective_rank_bam == 0))

case_counts <- aggregate(
  key ~ case,
  data = data.frame(case = bad$case[!duplicated(key)], key = unique(key)),
  FUN = length
)
names(case_counts)[2L] <- "invalid_unique_pairs"

task_dir <- file.path(
  repo,
  "artifacts/inla-bam-validation/inference/task-002-gaussian_spatial_null-r02"
)
fit <- readRDS(file.path(task_dir, "bam-compact-fit.rds"))
s1 <- mgcvST:::.mgcvst_model_score_state(fit, 1L)
s2 <- mgcvST:::.mgcvst_model_score_state(fit, 3L)
U <- sum(s1$a * s2$a)
moments <- mgcvST:::.rkhs_score_moments(s1$M, s2$M)
s <- mgcvST:::rkhs_score_singular_values(s1$M, s2$M)
liu <- mgcvST:::.liu_squared_score_moments(
  abs(U), moments[1L], moments[2L], moments[3L], moments[4L]
)
davies <- CompQuadForm::davies(abs(U), lambda = c(s / 2, -s / 2))

example <- data.frame(
  case = "gaussian_spatial_null",
  replicate = 2L,
  feature1 = "feature1",
  feature2 = "feature3",
  score = U,
  information = moments[1L],
  effective_rank_without_guard = moments[1L]^2 / moments[2L],
  positive_singular_values = length(s),
  maximum_singular_value = max(s),
  minimum_positive_singular_value = min(s),
  standardized_score = U / sqrt(moments[1L]),
  liu_p_without_guard = liu$p_value,
  davies_Qq_without_guard = davies$Qq,
  davies_ifault_without_guard = davies$ifault
)

overview <- data.frame(
  quantity = c(
    "all_method_rows", "rows_per_calibration", "unique_dataset_pairs",
    "invalid_bam_rows_both_methods", "invalid_bam_rows_per_method",
    "both_valid_rows_both_methods", "both_valid_rows_per_method",
    "minimum_invalid_information", "maximum_invalid_information"
  ),
  value = c(
    nrow(D), nrow(D) / 2, length(unique(paste(D$case, D$replicate,
      D$feature1, D$feature2, sep = "|"))), nrow(bad),
    unname(table(bad$calibration)[[1L]]), sum(D$both),
    sum(D$both & D$calibration == "liu"), min(bad$information_bam),
    max(bad$information_bam)
  )
)

write.csv(overview, file.path(summary_dir, "bam-pair-degeneracy-overview.csv"),
          row.names = FALSE)
write.csv(case_counts, file.path(summary_dir, "bam-pair-degeneracy-cases.csv"),
          row.names = FALSE)
write.csv(example, file.path(summary_dir, "bam-pair-degeneracy-example.csv"),
          row.names = FALSE)
