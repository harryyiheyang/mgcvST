Sys.setenv(
  LC_ALL = "C",
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1"
)

repo <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
lib <- normalizePath(Sys.getenv("MGCVST_VALIDATION_LIBRARY",
  "artifacts/inla-bam-validation/library"), mustWork = TRUE)
.libPaths(c(lib, .libPaths()))
library(mgcvST)
library(WGCNA)
library(mclust)
WGCNA::allowWGCNAThreads(nThreads = 4L)

source_root <- file.path(repo, "artifacts/inla-bam-validation/components")
out <- file.path(repo, "artifacts/inla-bam-validation/real-network-audit")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

datasets <- c("celltype", "unadjusted")
excluded <- c(
  celltype = "ENSDARG00000098458",
  unadjusted = "ENSDARG00000058105"
)
fit_rows <- list()
comparison_rows <- list()
fit_position <- 0L
comparison_position <- 0L

for (dataset in datasets) {
  path <- file.path(source_root, dataset)
  bam <- readRDS(file.path(path, "bam-compact-fit.rds"))
  inla <- readRDS(file.path(path, "inla-fit.rds"))
  stopifnot(identical(as.character(bam$feature_id), as.character(inla$feature_id)))
  stopifnot(excluded[[dataset]] %in% bam$feature_id)
  ids_full <- as.character(bam$feature_id)
  ids_sensitivity <- ids_full[ids_full != excluded[[dataset]]]

  fits <- list(bam = bam, inla = inla)
  for (scope in c("full", "exclude_nonconverged_inla")) {
    ids <- if (scope == "full") ids_full else ids_sensitivity
    results <- vector("list", 2L)
    names(results) <- c("bam", "inla")

    for (backend in c("bam", "inla")) {
      started <- proc.time()[["elapsed"]]
      network <- mgcvST.wgcna(
        fits[[backend]], indices = ids, group = "global", verbose = TRUE
      )
      elapsed <- proc.time()[["elapsed"]] - started
      results[[backend]] <- network
      saveRDS(
        network,
        file.path(out, paste(dataset, backend, scope, "wgcna.rds", sep = "-")),
        compress = FALSE
      )
      write.csv(
        network$modules,
        file.path(out, paste(dataset, backend, scope, "modules.csv", sep = "-")),
        row.names = FALSE
      )

      labels <- network$networks$selected$labels
      fit_position <- fit_position + 1L
      fit_rows[[fit_position]] <- data.frame(
        dataset = dataset,
        backend = backend,
        scope = scope,
        feature_count = length(ids),
        excluded_feature = if (scope == "full") NA_character_ else excluded[[dataset]],
        network_status = network$networks$selected$status,
        score_coordinate_count = network$networks$selected$q,
        assigned_feature_count = sum(labels != 0L),
        grey_feature_count = sum(labels == 0L),
        non_grey_module_count = length(unique(labels[labels != 0L])),
        score_seconds = network$timing$score_seconds,
        network_seconds = network$timing$network_seconds,
        reported_total_seconds = network$timing$total_seconds,
        measured_elapsed_seconds = elapsed,
        package_version = as.character(packageVersion("mgcvST")),
        stringsAsFactors = FALSE
      )
    }

    bam_network <- results$bam$networks$selected
    inla_network <- results$inla$networks$selected
    stopifnot(identical(bam_network$feature_id, inla_network$feature_id))
    R_difference <- bam_network$correlation - inla_network$correlation
    covariance_difference <- bam_network$covariance - inla_network$covariance
    comparison_position <- comparison_position + 1L
    comparison_rows[[comparison_position]] <- data.frame(
      dataset = dataset,
      scope = scope,
      feature_count = length(ids),
      excluded_feature = if (scope == "full") NA_character_ else excluded[[dataset]],
      module_ARI = mclust::adjustedRandIndex(
        bam_network$labels, inla_network$labels
      ),
      R_difference_frobenius = sqrt(sum(R_difference^2)),
      R_difference_mean_absolute = mean(abs(R_difference)),
      R_difference_max_absolute = max(abs(R_difference)),
      covariance_difference_frobenius = sqrt(sum(covariance_difference^2)),
      covariance_difference_mean_absolute = mean(abs(covariance_difference)),
      covariance_difference_max_absolute = max(abs(covariance_difference)),
      bam_non_grey_modules = length(unique(
        bam_network$labels[bam_network$labels != 0L]
      )),
      inla_non_grey_modules = length(unique(
        inla_network$labels[inla_network$labels != 0L]
      )),
      stringsAsFactors = FALSE
    )
    rm(results, bam_network, inla_network, R_difference, covariance_difference)
    gc()
  }
  rm(bam, inla, fits)
  gc()
}

fit_summary <- do.call(rbind, fit_rows)
comparison_summary <- do.call(rbind, comparison_rows)
write.csv(fit_summary, file.path(out, "network-fit-summary.csv"), row.names = FALSE)
write.csv(
  comparison_summary,
  file.path(out, "network-backend-comparison.csv"),
  row.names = FALSE
)

sensitivity_rows <- list()
position <- 0L
for (dataset in datasets) {
  for (backend in c("bam", "inla")) {
    full <- readRDS(file.path(out, paste(dataset, backend, "full-wgcna.rds", sep = "-")))
    sensitivity <- readRDS(file.path(
      out, paste(dataset, backend, "exclude_nonconverged_inla-wgcna.rds", sep = "-")
    ))
    full_network <- full$networks$selected
    sensitivity_network <- sensitivity$networks$selected
    keep <- match(sensitivity_network$feature_id, full_network$feature_id)
    R_difference <- full_network$correlation[keep, keep] -
      sensitivity_network$correlation
    position <- position + 1L
    sensitivity_rows[[position]] <- data.frame(
      dataset = dataset,
      backend = backend,
      retained_feature_count = length(keep),
      excluded_feature = excluded[[dataset]],
      module_ARI_full_vs_sensitivity = mclust::adjustedRandIndex(
        full_network$labels[keep], sensitivity_network$labels
      ),
      R_difference_mean_absolute = mean(abs(R_difference)),
      R_difference_max_absolute = max(abs(R_difference)),
      stringsAsFactors = FALSE
    )
    rm(full, sensitivity, full_network, sensitivity_network, R_difference)
    gc()
  }
}
write.csv(
  do.call(rbind, sensitivity_rows),
  file.path(out, "network-sensitivity-stability.csv"),
  row.names = FALSE
)
