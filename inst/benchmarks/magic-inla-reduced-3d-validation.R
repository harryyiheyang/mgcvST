#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1", BLIS_NUM_THREADS = "1")
if (!requireNamespace("pkgload", quietly = TRUE)) stop("pkgload is required.")
pkgload::load_all(".", export_all = FALSE, quiet = TRUE)
library(BiocParallel)

out <- Sys.getenv("MGCVST_MAGIC_REDUCED_OUTPUT", "artifacts/magic-inla-reduced-3d-validation")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
magic.file <- "artifacts/datasets/MAGIC/MAGIC.rds"
if (!file.exists(magic.file)) stop("MAGIC.rds is unavailable.")
M0 <- readRDS(magic.file)
d <- M0$covariates
id <- c("Snap25", "Foxp1", "Tfap2b")
if (!all(id %in% colnames(M0$expression))) stop("The three registered MAGIC genes are unavailable.")
Y <- t(M0$expression[, id, drop = FALSE])
mesh <- fmesher::fm_mesh_3d(loc = M0$meshes$native3d$loc, tv = M0$meshes$native3d$tv)
kappa <- M0$meshes$native3d$contract$kappa_fixed
flat <- list(prior = "flat", param = numeric(), initial = 0)
S <- mgcvST::inlaST.set(response ~ offset(log(exposure)), d, mesh = mesh,
  kappa = kappa, coordinates = c("x_mm", "y_mm", "z_mm"), family = mgcv::nb())
if (!identical(S$mean_constraint, "observation")) stop("The native observation-mean constraint is absent.")
ctl <- list(precision_prior = flat, nb_size_prior = flat, num_threads = 1L, keep_fit = FALSE)

t0 <- proc.time()[["elapsed"]]
checkpoint.file <- file.path(out, "fit-after-estimate.rds")
if (file.exists(checkpoint.file)) {
  F <- readRDS(checkpoint.file)
  fit.seconds <- NA_real_
  fit_source <- "existing checkpoint"
} else {
  F <- mgcvST::inlaST.estimate(Y, S, feature_id = id, control = ctl, diagnostics = TRUE,
    retain_smooth = TRUE, BPPARAM = SerialParam(), threads = 1L)
  fit.seconds <- proc.time()[["elapsed"]] - t0
  fit_source <- "new current API fit"
  if (as.numeric(object.size(F)) > 800 * 1024^2) stop("The recoverable fit exceeds the 0.8 GiB pre-write budget.")
  saveRDS(F, checkpoint.file, compress = FALSE)
}
if (!all(F$diagnostics$converged)) stop("At least one MAGIC fit did not converge.")
if (file.info(checkpoint.file)$size > 1000000000) stop("The recoverable fit exceeded the disk budget.")

P <- t(combn(seq_along(id), 2L))
t0 <- proc.time()[["elapsed"]]
public <- mgcvST::inlaST.test(F, pairs = t(combn(id, 2L)), calibration = "liu",
  BPPARAM = SerialParam(), threads = 1L)
public.seconds <- proc.time()[["elapsed"]] - t0

F <- mgcvST:::.inlast_sparse_prepare(F)
t0 <- proc.time()[["elapsed"]]
B <- mgcvST:::.inlast_sparse_observation_basis(F)
basis.cached.seconds <- proc.time()[["elapsed"]] - t0
t0 <- proc.time()[["elapsed"]]
U0 <- mgcvST:::.inlast_sparse_units(F, seq_len(nrow(Y)), threads = 1L)
direct <- mgcvST:::.inlast_sparse_materialize_reduced(F, U0, B, threads = 1L)
direct.seconds <- proc.time()[["elapsed"]] - t0
t0 <- proc.time()[["elapsed"]]
full.r <- mgcvST:::.mgcvst_inla_test_pairs(F, P, seq_len(nrow(P)), threads = 1L,
  chunk_size = 1L, verbose = FALSE, full_rank = TRUE)
full.r.seconds <- proc.time()[["elapsed"]] - t0
t0 <- proc.time()[["elapsed"]]
Bfull <- mgcvST:::.inlast_sparse_observation_basis(F, full_rank = TRUE)
direct.full <- mgcvST:::.inlast_sparse_materialize_reduced(F, U0, Bfull, threads = 1L)
full.direct.seconds <- proc.time()[["elapsed"]] - t0
t0 <- proc.time()[["elapsed"]]
old <- mgcvST:::.inlast_sparse_materialize(F, U0, threads = 1L)
old.materialize.seconds <- proc.time()[["elapsed"]] - t0
t0 <- proc.time()[["elapsed"]]
z <- mgcvST:::.inlast_sparse_batch(F, seq_len(nrow(Y)), threads = 1L)
old.batch.seconds <- proc.time()[["elapsed"]] - t0

bad <- vapply(c(U0, direct, direct.full, old, z), function(x) !is.null(x$error) && nzchar(x$error), logical(1L))
if (any(bad)) stop("A sparse score stage failed.")
M <- lapply(z, `[[`, "M")
Md <- lapply(direct, `[[`, "M")
Mfull <- lapply(direct.full, `[[`, "M")
Mr <- lapply(old, `[[`, "M")
a <- lapply(z, `[[`, "a")
t0 <- proc.time()[["elapsed"]]; T0 <- mgcvST:::mgcvst_pair_trace_powers_cpp(M, P, 4L, 1L); trace.full.seconds <- proc.time()[["elapsed"]] - t0
t0 <- proc.time()[["elapsed"]]; Told <- mgcvST:::mgcvst_pair_trace_powers_cpp(Mr, P, 4L, 1L); trace.old.materialize.seconds <- proc.time()[["elapsed"]] - t0
t0 <- proc.time()[["elapsed"]]; Td <- mgcvST:::mgcvst_pair_trace_powers_cpp(Md, P, 4L, 1L); trace.reduced.seconds <- proc.time()[["elapsed"]] - t0
t0 <- proc.time()[["elapsed"]]; Tfull <- mgcvST:::mgcvst_pair_trace_powers_cpp(Mfull, P, 4L, 1L); trace.full.r.seconds <- proc.time()[["elapsed"]] - t0

R1 <- list()
for (j in seq_len(nrow(P))) {
  i1 <- P[j, 1L]; i2 <- P[j, 2L]
  u.full <- sum(a[[i1]] * a[[i2]])
  u.red <- sum(direct[[i1]]$a * direct[[i2]]$a)
  p.full <- mgcvST:::.liu_squared_score_moments(u.full, T0[j, 1L], T0[j, 2L], T0[j, 3L], T0[j, 4L])$p_value
  p.red <- mgcvST:::.liu_squared_score_moments(u.red, Td[j, 1L], Td[j, 2L], Td[j, 3L], Td[j, 4L])$p_value
  R1[[j]] <- data.frame(feature1 = id[i1], feature2 = id[i2], U_full = u.full, U_full_r = full.r$result$signed_score[j], U_reduced = u.red, U_public = public$results$signed_score[j],
    max_full_materialize_trace_difference = max(abs(Told[j, ] - T0[j, ])), max_full_r_trace_difference = max(abs(Tfull[j, ] - T0[j, ])),
    max_reduced_projection_difference = max(abs(Md[[i1]] - crossprod(B$coordinate, M[[i1]] %*% B$coordinate)), abs(Md[[i2]] - crossprod(B$coordinate, M[[i2]] %*% B$coordinate))),
    c1_full = T0[j, 1L], c2_full = T0[j, 2L], c3_full = T0[j, 3L], c4_full = T0[j, 4L],
    c1_reduced = Td[j, 1L], c2_reduced = Td[j, 2L], c3_reduced = Td[j, 3L], c4_reduced = Td[j, 4L],
    liu_full = p.full, liu_full_r = full.r$result$p_two_sided[j], liu_reduced = p.red, liu_public = public$results$p_two_sided[j],
    delta_log10_p = log10(p.red) - log10(p.full), delta_log10_public_p = log10(public$results$p_two_sided[j]) - log10(p.red))
}
write.csv(do.call(rbind, R1), file.path(out, "pairwise-comparison.csv"), row.names = FALSE)
write.csv(data.frame(observations = nrow(d), q = ncol(M[[1L]]), r = B$rank,
  target_coverage = B$coverage, kept_coverage = B$kept, tail = B$tail, fit_source = fit_source,
  fit_seconds = fit.seconds, basis_cached_seconds = basis.cached.seconds, direct_unit_and_M_seconds = direct.seconds,
  full_r_private_seconds = full.r.seconds, full_r_direct_seconds = full.direct.seconds, old_materialize_seconds = old.materialize.seconds,
  old_batch_seconds = old.batch.seconds, trace_full_seconds_3_pairs = trace.full.seconds,
  trace_old_materialize_seconds_3_pairs = trace.old.materialize.seconds,
  trace_reduced_seconds_3_pairs = trace.reduced.seconds, trace_full_r_seconds_3_pairs = trace.full.r.seconds,
  public_test_wall_seconds = public.seconds,
  public_basis_seconds = public$timing$inla_projection$basis_elapsed,
  public_unit_build_seconds = public$timing$inla_projection$unit_build_elapsed,
  public_reduced_materialize_seconds = public$timing$inla_projection$reduced_materialize_elapsed,
  public_liu_seconds = public$timing$inla_projection$liu_elapsed,
  public_test_wall_recorded_seconds = public$timing$inla_projection$test_wall_elapsed,
  public_pair_schedule = public$timing$inla_projection$pair_schedule),
  file.path(out, "timing.csv"), row.names = FALSE)
g <- gc()
p <- as.numeric(system2("powershell", c("-NoProfile", "-Command", paste0("(Get-Process -Id ", Sys.getpid(), ").PeakWorkingSet64")), stdout = TRUE))
files <- list.files(out, recursive = TRUE, full.names = TRUE)
write.csv(data.frame(output_bytes = sum(file.info(files)$size), r_heap_gc_max_mib_approx = sum(g[, "max used"] * c(56, 8)) / 1024^2,
  windows_peak_working_set_bytes = p, disk_budget_bytes = 1000000000,
  note = "F is recoverable; its external-pointer sparse cache is rebuilt by .inlast_sparse_prepare after readRDS."),
  file.path(out, "resource-use.csv"), row.names = FALSE)
capture.output(sessionInfo(), file = file.path(out, "session-info.txt"))
