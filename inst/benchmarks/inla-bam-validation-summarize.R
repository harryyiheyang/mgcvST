library(data.table)
setDTthreads(1L)
out <- "artifacts/inla-bam-validation"
dest <- file.path(out, "summary")
dir.create(dest, recursive = TRUE, showWarnings = FALSE)
M <- fread(file.path(out, "inference/manifest.csv"))
P <- J <- T <- vector("list", nrow(M))
for (i in seq_len(nrow(M))) {
  z <- M[i]
  folder <- sprintf("task-%03d-%s-r%02d", z$task, z$case, z$replicate)
  lane <- if (z$nuisance == "random") "nuisance-correction" else "inference"
  src <- file.path(out, lane, folder)
  p <- fread(file.path(src, "pair-tests.csv"))
  m <- fread(file.path(src, "marginal-tests.csv"))
  t <- fread(file.path(src, "timings.csv"))
  for (d in list(p, m, t)) {
    d[, `:=`(case = z$case, replicate = z$replicate, seed = z$seed,
              rho = z$rho, source_lane = lane)]
  }
  m[, p := fifelse(backend == "bam", smooth.pvalue, p_value)]
  P[[i]] <- p
  J[[i]] <- m
  T[[i]] <- t
}
P <- rbindlist(P, fill = TRUE)
J <- rbindlist(J, fill = TRUE)
T <- rbindlist(T, fill = TRUE)
P[, valid := is.finite(p_two_sided) & p_two_sided >= 0 & p_two_sided <= 1]
J[, valid := is.finite(p) & p >= 0 & p <= 1]
PS <- P[, .(independent_datasets = uniqueN(replicate), attempted_pairs = .N,
  valid_pairs = sum(valid), invalid_pairs = sum(!valid),
  raw_rejections = sum(p_two_sided[valid] < .05),
  adjusted_rejections = sum(p_adjusted[valid] <= .05),
  median_p = if (any(valid)) median(p_two_sided[valid]) else NA_real_),
  by = .(case, backend, calibration)]
JS <- J[, .(independent_datasets = uniqueN(replicate), attempted_features = .N,
  valid_features = sum(valid), invalid_features = sum(!valid),
  rejected_features = sum(p[valid] < .05),
  median_p = if (any(valid)) median(p[valid]) else NA_real_),
  by = .(case, backend, calibration)]
PJ <- merge(P[backend == "bam"], P[backend == "inla"],
  by = c("case", "replicate", "seed", "feature1", "feature2", "calibration"),
  suffixes = c("_bam", "_inla"), sort = FALSE)
PJ[, both := valid_bam & valid_inla]
PC <- PJ[, .(independent_datasets = uniqueN(replicate), attempted_pairs = .N,
  both_valid = sum(both),
  median_abs_p_difference = if (any(both))
    median(abs(p_two_sided_bam[both] - p_two_sided_inla[both])) else NA_real_,
  max_abs_p_difference = if (any(both))
    max(abs(p_two_sided_bam[both] - p_two_sided_inla[both])) else NA_real_,
  raw_decision_agreement = if (any(both))
    mean((p_two_sided_bam[both] < .05) == (p_two_sided_inla[both] < .05)) else NA_real_,
  adjusted_decision_agreement = if (any(both))
    mean((p_adjusted_bam[both] <= .05) == (p_adjusted_inla[both] <= .05)) else NA_real_,
  score_sign_agreement = if (any(both))
    mean(sign(signed_score_bam[both]) == sign(signed_score_inla[both])) else NA_real_),
  by = .(case, calibration)]
MJ <- merge(J[backend == "bam"], J[backend == "inla"],
  by = c("case", "replicate", "seed", "feature_id", "calibration"),
  suffixes = c("_bam", "_inla"), sort = FALSE)
MJ[, both := valid_bam & valid_inla]
JC <- MJ[, .(independent_datasets = uniqueN(replicate), attempted_features = .N,
  both_valid = sum(both),
  median_abs_p_difference = if (any(both)) median(abs(p_bam[both] - p_inla[both])) else NA_real_,
  max_abs_p_difference = if (any(both)) max(abs(p_bam[both] - p_inla[both])) else NA_real_,
  decision_agreement = if (any(both)) mean((p_bam[both] < .05) == (p_inla[both] < .05)) else NA_real_),
  by = .(case, calibration)]
fwrite(PS, file.path(dest, "inference-pair-counts.csv"))
fwrite(JS, file.path(dest, "inference-marginal-counts.csv"))
fwrite(PC, file.path(dest, "inference-pair-agreement.csv"))
fwrite(JC, file.path(dest, "inference-marginal-agreement.csv"))
fwrite(PJ, file.path(dest, "inference-paired-results.csv"))
fwrite(MJ, file.path(dest, "inference-paired-marginals.csv"))
fwrite(T[, .(independent_datasets = uniqueN(replicate),
  median_fit_seconds = median(fit_seconds),
  min_fit_seconds = min(fit_seconds), max_fit_seconds = max(fit_seconds)),
  by = .(case, backend)], file.path(dest, "inference-timing.csv"))

S <- fread(file.path(out, "scaling/manifest.csv"))
L <- vector("list", nrow(S))
for (i in seq_len(nrow(S))) {
  src <- file.path(out, "scaling", sprintf("n%d_r%02d", S$n[i], S$replicate[i]))
  if (!file.exists(file.path(src, "result.rds"))) stop("Incomplete scaling case: ", src)
  L[[i]] <- fread(file.path(src, "metrics.csv"))
}
L <- rbindlist(L)
LS <- L[, .(independent_datasets = uniqueN(replicate), paired_feature_comparisons = .N,
  median_backend_field_correlation = median(field_correlation),
  min_backend_field_correlation = min(field_correlation),
  median_field_rmse_bam = median(field_rmse_bam),
  median_field_rmse_inla = median(field_rmse_inla),
  median_inla_estimator_seconds = median(inla_estimator_seconds),
  median_bam_fit_seconds = median(bam_fit_seconds),
  median_bam_estimator_seconds = median(bam_fit_seconds + bam_marginal_seconds + bam_compact_seconds),
  median_bam_over_inla_estimator_time =
    median((bam_fit_seconds + bam_marginal_seconds + bam_compact_seconds) / inla_estimator_seconds),
  median_inla_total_seconds = median(inla_setup_seconds + inla_estimator_seconds + inla_pair_seconds),
  median_bam_total_seconds = median(bam_setup_seconds + bam_fit_seconds +
    bam_marginal_seconds + bam_compact_seconds + bam_pair_seconds),
  median_bam_over_inla_total_time = median((bam_setup_seconds + bam_fit_seconds +
    bam_marginal_seconds + bam_compact_seconds + bam_pair_seconds) /
    (inla_setup_seconds + inla_estimator_seconds + inla_pair_seconds)),
  median_inla_pair_seconds = median(inla_pair_seconds),
  median_bam_pair_seconds = median(bam_pair_seconds),
  median_abs_pair_p_difference = median(abs(bam_pair_p - inla_pair_p)),
  max_abs_pair_p_difference = max(abs(bam_pair_p - inla_pair_p))), by = .(n, q)]
fwrite(LS, file.path(dest, "scaling-summary.csv"))
fwrite(L, file.path(dest, "scaling-replicates.csv"))
