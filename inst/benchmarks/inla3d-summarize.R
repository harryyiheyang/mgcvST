out <- "artifacts/inla3d/eb"
dest <- "inst/validation/inla3d"
dir.create(dest, recursive = TRUE, showWarnings = FALSE)
R1 <- read.csv(file.path(out, "fits.csv"))
M <- read.csv("artifacts/inla3d/mesh/manifest.csv")
input <- readRDS(file.path(out, "input.rds"))
if (nrow(R1) != 80L || anyDuplicated(R1[c("scenario", "replicate", "mesh")]) ||
    any(table(R1$scenario, R1$mesh) != 10L)) stop("Expected 80 completed, unique fits.")
if (!identical(dim(input$Y), c(30000L, 2L, 10L))) stop("Unexpected input dimensions.")
S <- list()
D <- list()
P <- list()
T <- list()
j <- k <- 0L
for (s in c("broad", "focal")) {
  for (m in M$mesh) {
    z <- R1[R1$scenario == s & R1$mesh == m, ]
    j <- j + 1L
    S[[j]] <- data.frame(scenario = s, mesh = m, nodes = z$nodes[1L],
      replicates = nrow(z), converged = sum(z$mode_status == 0),
      median_fit_seconds = median(z$fit_seconds), min_fit_seconds = min(z$fit_seconds),
      max_fit_seconds = max(z$fit_seconds), median_correlation = median(z$correlation),
      median_rmse = median(z$rmse), median_roi_rmse = median(z$roi_rmse),
      median_outside_rmse = median(z$outside_rmse),
      min_nb_size = min(z$nb_size_mode), max_nb_size = max(z$nb_size_mode),
      max_constraint_error = max(abs(z$constraint_error)))
  }
  for (r in seq_len(10L)) {
    est <- list()
    for (m in M$mesh) {
      f <- readRDS(file.path(out, "fits", paste0(s, "-", r, "-", m, ".rds")))
      stopifnot(identical(f$flat$prior, "flat"), !f$flat$fixed,
        identical(f$integration, "eb"), identical(f$latent, "gaussian"),
        identical(names(f$mode$theta), c(
          "log size for the nbinomial observations (1/overdispersion)",
          "Log precision for field")), length(f$field) == 30000L)
      est[[m]] <- f$field
      k <- k + 1L
      T[[k]] <- data.frame(scenario = s, replicate = r, mesh = m,
        pre_seconds = unname(f$cpu["Pre"]), running_seconds = unname(f$cpu["Running"]),
        post_seconds = unname(f$cpu["Post"]), total_seconds = unname(f$cpu["Total"]))
      g <- input$truth[, s]
      roi <- input$roi
      D[[k]] <- data.frame(scenario = s, replicate = r, mesh = m,
        truth_contrast = mean(g[roi]) - mean(g[!roi]),
        estimated_contrast = mean(f$field[roi]) - mean(f$field[!roi]))
    }
    pairs <- list(c("uniform-1500", "adaptive-1500"),
                  c("uniform-2800", "adaptive-2800"),
                  c("uniform-1500", "uniform-2800"),
                  c("adaptive-1500", "adaptive-2800"))
    for (pair in pairs) {
      a <- R1[R1$scenario == s & R1$replicate == r & R1$mesh == pair[1L], ]
      b <- R1[R1$scenario == s & R1$replicate == r & R1$mesh == pair[2L], ]
      P[[length(P) + 1L]] <- data.frame(scenario = s, replicate = r,
        from = pair[1L], to = pair[2L],
        field_correlation = cor(est[[pair[1L]]], est[[pair[2L]]]),
        field_rms_change = sqrt(mean((est[[pair[1L]]] - est[[pair[2L]]])^2)),
        rmse_change_percent = 100 * (b$rmse / a$rmse - 1),
        roi_rmse_change_percent = 100 * (b$roi_rmse / a$roi_rmse - 1),
        outside_rmse_change_percent = 100 * (b$outside_rmse / a$outside_rmse - 1))
    }
  }
}
D <- do.call(rbind, D)
D$contrast_error <- D$estimated_contrast - D$truth_contrast
write.csv(R1, file.path(dest, "fits.csv"), row.names = FALSE)
write.csv(do.call(rbind, S), file.path(dest, "summary.csv"), row.names = FALSE)
write.csv(do.call(rbind, P), file.path(dest, "paired-comparisons.csv"), row.names = FALSE)
write.csv(D, file.path(dest, "region-contrasts.csv"), row.names = FALSE)
write.csv(do.call(rbind, T), file.path(dest, "native-timings.csv"), row.names = FALSE)
file.copy("artifacts/inla3d/mesh/manifest.csv", file.path(dest, "mesh-manifest.csv"), overwrite = TRUE)
file.copy(file.path(out, "geometry.csv"), file.path(dest, "geometry.csv"), overwrite = TRUE)
dir.create(file.path(dest, "meshes"), showWarnings = FALSE)
mesh_files <- list.files("artifacts/inla3d/mesh", "-(vertices|tetrahedra)\\.csv$", full.names = TRUE)
if (length(mesh_files) != 8L) stop("Expected the eight fixed mesh input files.")
copied <- file.copy(mesh_files, file.path(dest, "meshes"), overwrite = TRUE)
if (!all(copied)) stop("Could not publish the fixed mesh inputs.")
if (file.exists(file.path(out, "memory.csv"))) {
  mem <- read.csv(file.path(out, "memory.csv"))
  mm <- aggregate(cbind(rss_MiB, private_MiB) ~ task, mem, max)
  write.csv(mm, file.path(dest, "sampled-memory.csv"), row.names = FALSE)
} else {
  message("No optional process memory samples supplied; fit summaries remain complete.")
}
files <- c("inst/benchmarks/inla3d-mesh.py", "inst/benchmarks/inla3d-validation.R",
           "inst/benchmarks/inla3d-geometry-check.R", "inst/benchmarks/inla3d-summarize.R",
           file.path(out, "input.rds"), list.files("artifacts/inla3d/mesh", "csv$", full.names = TRUE))
write.csv(data.frame(file = files, md5 = unname(tools::md5sum(files))),
          file.path(dest, "input-source-hashes.csv"), row.names = FALSE)
print(do.call(rbind, S))
