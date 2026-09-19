# Compact, reviewable records; retain full fits and spectra in artifacts.
root <- "artifacts/inla-stress-calibration"
src <- file.path(root, "summary")
out <- "inst/validation/inla-stress"
N <- read.csv(file.path(src, "null-completion.csv"))
S <- read.csv(file.path(src, "stress-completion.csv"))
stopifnot(nrow(N) == 8L, all(N$completed == 500L), all(N$final),
  all(S$summary_status == "final"), S$recorded[S$lane == "all"] == 70L)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
files <- c("null-completion.csv", "null-hyperparameter-summary.csv",
  "null-negative-eigenvalues.csv", "null-p-values.csv", "null-rejection-summary.csv",
  "null-spectrum-summary.csv", "null-state-matrix-aggregate.csv",
  "stress-completion.csv", "stress-fit-summary.csv", "stress-independent-agreement.csv",
  "stress-joint-repeatability.csv", "stress-supervisor-groups.csv", "stress-supervisor-jobs.csv")
stopifnot(all(file.copy(file.path(src, files), file.path(out, files), overwrite = TRUE)))
H <- read.csv(file.path(src, "null-state-matrix-summary.csv"))
H <- H[, c("case", "replicate", "state", "minimum_eigenvalue", "maximum_eigenvalue",
  "relative_minimum_eigenvalue", "negative_below_existing_cutoff", "negative_mass_fraction",
  "Vp_ratio", "expected_spectrum_available", "expected_spectrum_unavailable_reason",
  "symmetry_max_abs", "reconstruction_max_abs")]
stopifnot(nrow(H) == 4000L)
write.csv(H, file.path(out, "null-state-diagnostics.csv"), row.names = FALSE)

K <- list()
P <- list()
for (case in c("independent-0.3", "independent-3", "joint-0.3", "joint-3")) {
  folder <- file.path(root, "null-sensitivity", case)
  K[[case]] <- read.csv(file.path(folder, "calibration-summary.csv"))
  stopifnot(all(K[[case]]$analysis_status == "complete"))
  R <- read.csv(file.path(folder, "replicate-results.csv"))
  R <- R[R$route != "native", c("case", "replicate", "calibration", "route",
    "p_value", "statistic", "information", "unavailable_reason")]
  stopifnot(nrow(R) == 3000L, !anyDuplicated(R[, 1:4]))
  R$id <- paste(R$route, R$calibration, sep = "_")
  R$calibration <- R$route <- NULL
  P[[case]] <- reshape(R, direction = "wide", idvar = c("case", "replicate"), timevar = "id")
}
write.csv(do.call(rbind, K), file.path(out, "posthoc-calibration-summary.csv"), row.names = FALSE)
write.csv(do.call(rbind, P), file.path(out, "posthoc-paired-results.csv"), row.names = FALSE)

pub <- "artifacts/inla-public-parallel-stress"
dest <- file.path(out, "public-parallel")
dir.create(dest, recursive = TRUE, showWarnings = FALSE)
files <- c("comparison-summary.csv", "feature-comparison.csv", "features.csv",
  "pair-comparison.csv", "provenance.csv", "session-info.txt", "timing.csv")
stopifnot(all(file.copy(file.path(pub, files), file.path(dest, files), overwrite = TRUE)))
for (config in c("serial", "snow2", "snow4")) {
  for (file in c("fit-diagnostics.csv", "hyper-diagnostics.csv", "workers.csv")) {
    stopifnot(file.copy(file.path(pub, config, file),
      file.path(dest, paste(config, file, sep = "-")), overwrite = TRUE))
  }
}

files <- c("input-source-sha256.csv", "session.txt")
stopifnot(all(file.copy(file.path(root, files), file.path(out, files), overwrite = TRUE)))
C <- read.csv(file.path(root, "verification", "checks.csv"))
stopifnot(all(C$passed[C$check != "joint_completed_status"]))
stopifnot(file.copy(file.path(root, "verification", "checks-summary.csv"),
  file.path(out, "verification-summary.csv"), overwrite = TRUE))
write.csv(C[C$check %in% c("paired_Y_replicates_compared", "paired_Y_genes_1_2_identical",
  "pair_is_one_replicate", "joint_completed_status") | !C$passed, ],
  file.path(out, "verification-special-checks.csv"), row.names = FALSE)

M <- read.csv(file.path(root, "full-process-samples.csv"))
if (nrow(M)) {
  T <- aggregate(private_bytes ~ utc, M, sum)
  peak <- T$utc[which.max(T$private_bytes)]
  write.csv(M[M$utc == peak, ], file.path(out, "full-process-peak.csv"), row.names = FALSE)
  H <- list()
  for (pid in unique(M$pid)) {
    D <- M[M$pid == pid, ]
    tt <- as.POSIXct(sub("[+]00:00$", "", D$utc), format = "%Y-%m-%dT%H:%M:%OS", tz = "UTC")
    dt <- as.numeric(difftime(max(tt), min(tt), units = "secs"))
    cpu <- max(D$cpu_seconds) - min(D$cpu_seconds)
    H[[as.character(pid)]] <- data.frame(pid = pid, name = D$name[1L], samples = nrow(D),
      sampled_wall_seconds = dt, sampled_cpu_seconds = cpu,
      average_cpu_cores = if (dt > 0) cpu / dt else NA_real_,
      maximum_threads = max(D$threads), peak_private_bytes = max(D$private_bytes))
  }
  write.csv(do.call(rbind, H), file.path(out, "full-process-summary.csv"), row.names = FALSE)
}
