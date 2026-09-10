#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
Sys.setenv(LC_ALL = "C", OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1", BLIS_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1", RCPP_PARALLEL_NUM_THREADS = "1")

out <- Sys.getenv("MGCVST_COMPONENT_OUTPUT",
  "artifacts/inla-bam-validation/components")
lib <- normalizePath(Sys.getenv("MGCVST_VALIDATION_LIBRARY",
  "artifacts/inla-bam-validation/library"), mustWork = TRUE)
.libPaths(c(lib, .libPaths()))
library(BiocParallel)
library(mgcvST)

flat <- list(prior = "flat", param = numeric(), initial = 0)
ctl <- list(precision_prior = flat, nb_size_prior = flat,
  num_threads = 1L, keep_fit = FALSE)
datasets <- strsplit(Sys.getenv("MGCVST_COMPONENT_DATASETS",
  "unadjusted,celltype"), ",", fixed = TRUE)[[1L]]

for (dataset in datasets) {
  d.out <- file.path(out, dataset)
  input <- readRDS(file.path(d.out, "input.rds"))
  model <- readRDS(file.path(d.out, "inla-model.rds"))
  original <- readRDS(file.path(d.out, "inla-fit.rds"))
  status <- vapply(original$inla_diagnostics, function(x) x$mode_status,
    integer(1L))
  jj <- which(status != 0L)
  if (length(jj) != 1L) stop("Expected one nonzero mode status for ", dataset, ".")
  t0 <- proc.time()[["elapsed"]]
  repeat.fit <- inlaST.estimate(input$Y[jj, , drop = FALSE], model,
    feature_id = input$feature_id[jj], retain_marginal = TRUE,
    marginal_args = list(method = "liu"), diagnostics = TRUE,
    BPPARAM = SerialParam(), chunk_size = 1L, control = ctl,
    score_backend = "sparse")
  seconds <- proc.time()[["elapsed"]] - t0
  result <- list(dataset = dataset, feature_id = input$feature_id[jj],
    original_diagnostics = original$diagnostics[jj, , drop = FALSE],
    original_inla_diagnostics = original$inla_diagnostics[jj],
    repeat_fit = repeat.fit, repeat_seconds = seconds, control = ctl)
  saveRDS(result, file.path(d.out, "inla-mode2-repeat.rds"), compress = FALSE)

  om <- mgcvST.marginal(original, features = input$feature_id[jj],
    calibration = "liu", BPPARAM = SerialParam())
  rm <- mgcvST.marginal(repeat.fit, calibration = "liu",
    BPPARAM = SerialParam())
  z <- data.frame(dataset = dataset, feature_id = input$feature_id[jj],
    original_converged = original$diagnostics$converged[jj],
    repeat_converged = repeat.fit$diagnostics$converged[1L],
    original_mode_status = status[jj],
    repeat_mode_status = repeat.fit$inla_diagnostics[[1L]]$mode_status,
    original_lambda = original$lambda[jj], repeat_lambda = repeat.fit$lambda[1L],
    original_dispersion = original$dispersion[jj],
    repeat_dispersion = repeat.fit$dispersion[1L],
    original_nb_size = original$family_parameters[[jj]],
    repeat_nb_size = repeat.fit$family_parameters[[1L]],
    max_abs_working_error = max(abs(original$working_error[, jj] -
      repeat.fit$working_error[, 1L])),
    max_abs_working_variance = max(abs(original$working_variance[, jj] -
      repeat.fit$working_variance[, 1L])),
    max_abs_constraint_residual = abs(original$constraint_residual[jj] -
      repeat.fit$constraint_residual[1L]),
    original_marginal_p = om$p_value[1L], repeat_marginal_p = rm$p_value[1L],
    seconds = seconds)
  write.csv(z, file.path(d.out, "inla-mode2-repeat.csv"), row.names = FALSE)
}
