# Inspect completed flat-prior fits without repeating the optimization.
MAGIC <- readRDS("artifacts/datasets/MAGIC/MAGIC.rds")
out <- "artifacts/magic-slide-exploration/fits"
cfg <- readRDS(file.path(out, "configuration.rds"))
priors <- c(cfg$spatial_prior$prior, cfg$nb_prior$prior, cfg$slide_prior$prec$prior,
  cfg$ou_prior$prec$prior, cfg$ou_prior$phi$prior)
stopifnot(length(priors) == 5L, all(priors == "flat"))
R1 <- list()
V <- matrix(NA_real_, nrow(MAGIC$slices), ncol(MAGIC$expression),
  dimnames = list(MAGIC$slices$slice_id, colnames(MAGIC$expression)))
j <- 0L
for (gene in colnames(MAGIC$expression)) {
  f0 <- readRDS(file.path(out, paste0(gene, "-spatial.rds")))
  for (model in c("spatial", "iid", "ou")) {
    j <- j + 1L
    f <- readRDS(file.path(out, paste0(gene, "-", model, ".rds")))
    stopifnot(identical(f$point_id, MAGIC$covariates$point_id),
      f$metrics$mode_status == 0, f$metrics$warnings == 0,
      max(abs(f$mu - MAGIC$covariates$exposure * exp(f$eta))) < 1e-10)
    R1[[j]] <- data.frame(gene = gene, model = model,
      spatial_correlation_to_baseline = cor(f$spatial, f0$spatial),
      spatial_rms_difference = sqrt(mean((f$spatial - f0$spatial)^2)),
      fitted_rate_correlation_to_baseline = cor(exp(f$eta), exp(f0$eta)))
    if (model == "ou") V[, gene] <- f$slide_effect
  }
}
write.csv(do.call(rbind, R1), file.path(out, "component-comparison.csv"), row.names = FALSE)
write.csv(cor(V), file.path(out, "estimated-slide-component-correlations.csv"))
old <- readRDS("artifacts/inla3d-transfer/generic-flat/fits/fit-1.rds")
new <- readRDS(file.path(out, "Snap25-spatial.rds"))
write.csv(data.frame(n = nrow(MAGIC$covariates), genes = ncol(MAGIC$expression),
  metadata_fields = ncol(MAGIC$covariates), slices = nrow(MAGIC$slices),
  all_priors_flat = all(priors == "flat"), fits_validated = j,
  snap25_max_abs_eta_difference_from_saved_baseline = max(abs(new$eta - old$eta)),
  total_count_fields_differ = sum(MAGIC$covariates$total_counts != MAGIC$covariates$total_umi),
  section_order_reproduces_z = all(MAGIC$slices$z ==
    10 * (MAGIC$slices$section_seq_id - MAGIC$slices$section_seq_id[1L]))),
  file.path(out, "validation.csv"), row.names = FALSE)
print(read.csv(file.path(out, "validation.csv")))
print(do.call(rbind, R1))
