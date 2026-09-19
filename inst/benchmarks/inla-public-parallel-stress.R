#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
Sys.setenv(LC_ALL = "C", OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1", BLIS_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1", RCPP_PARALLEL_NUM_THREADS = "1")

library(BiocParallel)
library(mgcvST)

if (as.character(packageVersion("mgcvST")) != "0.0.1.9006") {
  stop("The public parallel stress check requires mgcvST 0.0.1.9006.")
}
Sys.setenv(R_LIBS_USER = paste(.libPaths(), collapse = .Platform$path.sep))

source.dir <- Sys.getenv("MGCVST_COMPONENT_SOURCE",
  "artifacts/inla-bam-validation/components/unadjusted")
out <- Sys.getenv("MGCVST_PUBLIC_PARALLEL_OUTPUT",
  "artifacts/inla-public-parallel-stress")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
input <- readRDS(file.path(source.dir, "input.rds"))
S <- readRDS(file.path(source.dir, "inla-model.rds"))
ids <- input$feature_id[1:10]
Y <- input$Y[1:10, , drop = FALSE]
if (!identical(rownames(Y), ids)) rownames(Y) <- ids
if (nrow(Y) != 10L || ncol(Y) != 2125L || ncol(input$basis$B) != 298L) {
  stop("The frozen Visium-B component input dimensions changed.")
}

flat <- list(prior = "flat", param = numeric(), initial = 0)
ctl <- list(precision_prior = flat, nb_size_prior = flat,
  num_threads = 1L, keep_fit = FALSE)
source.provenance <- read.csv(file.path(dirname(source.dir),
  "package-provenance.csv"))
write.csv(data.frame(package = "mgcvST",
  version = as.character(packageVersion("mgcvST")),
  library = find.package("mgcvST"), source_input = normalizePath(source.dir),
  source_model_package_version = source.provenance$version[1L],
  source_model_package_commit = source.provenance$package_source_commit[1L],
  genes = nrow(Y), observations = ncol(Y), coefficients = ncol(input$basis$B),
  precision_prior = "flat", nb_size_prior = "flat", inla_threads = 1L),
  file.path(out, "provenance.csv"), row.names = FALSE)
write.csv(data.frame(index = seq_along(ids), feature_id = ids),
  file.path(out, "features.csv"), row.names = FALSE)

configs <- data.frame(config = c("serial", "snow2", "snow4"),
  workers = c(1L, 2L, 4L), stringsAsFactors = FALSE)
fits <- list()
tests <- list()
timing <- list()
pairs <- t(combn(ids, 2L))
colnames(pairs) <- c("feature1", "feature2")

for (i in seq_len(nrow(configs))) {
  config <- configs$config[i]
  workers <- configs$workers[i]
  d.out <- file.path(out, config)
  dir.create(d.out, recursive = TRUE, showWarnings = FALSE)
  fit.file <- file.path(d.out, "fit.rds")
  test.file <- file.path(d.out, "pair-tests.rds")
  timing.file <- file.path(d.out, "timing.rds")
  if (config == "serial") {
    BP <- SerialParam(stop.on.error = FALSE)
    worker <- data.frame(worker = 1L, pid = Sys.getpid(),
      version = as.character(packageVersion("mgcvST")),
      library = find.package("mgcvST"))
  } else {
    BP <- SnowParam(workers = workers, type = "SOCK", tasks = 0L,
      stop.on.error = FALSE, progressbar = FALSE)
    BP <- bpstart(BP)
    worker <- do.call(rbind, bplapply(seq_len(workers), function(j) {
      data.frame(worker = j, pid = Sys.getpid(),
        version = as.character(packageVersion("mgcvST")),
        library = find.package("mgcvST"))
    }, BPPARAM = BP))
  }
  write.csv(worker, file.path(d.out, "workers.csv"), row.names = FALSE)
  if (nrow(worker) != workers || length(unique(worker$pid)) != workers ||
      any(worker$version != "0.0.1.9006") ||
      length(unique(worker$library)) != 1L) {
    stop("Worker library validation failed for ", config, ".")
  }

  if (file.exists(fit.file)) {
    fit <- readRDS(fit.file)
    fit.seconds <- if (file.exists(timing.file)) readRDS(timing.file)$fit_seconds else
      NA_real_
  } else {
    t0 <- proc.time()[["elapsed"]]
    fit <- inlaST.estimate(Y, S, feature_id = ids, control = ctl,
      diagnostics = TRUE, retain_smooth = TRUE, retain_marginal = TRUE,
      BPPARAM = BP, chunk_size = ceiling(length(ids) / workers))
    fit.seconds <- proc.time()[["elapsed"]] - t0
    saveRDS(fit, fit.file, compress = FALSE)
  }
  write.csv(fit$diagnostics, file.path(d.out, "fit-diagnostics.csv"),
    row.names = FALSE)
  hyper <- data.frame(feature_id = fit$feature_id,
    mode_status = vapply(fit$inla_diagnostics, function(x) x$mode_status,
      integer(1L)),
    warning_count = vapply(fit$inla_diagnostics, function(x)
      length(x$estimation$hyper_mode_diagnostics$warnings), integer(1L)),
    lambda = fit$lambda,
    nb_size = vapply(fit$family_parameters, function(x) as.numeric(x[1L]),
      numeric(1L)))
  write.csv(hyper, file.path(d.out, "hyper-diagnostics.csv"), row.names = FALSE)
  if (file.exists(test.file)) {
    test <- readRDS(test.file)
    pair.seconds <- if (file.exists(timing.file)) readRDS(timing.file)$pair_seconds else
      NA_real_
  } else {
    t0 <- proc.time()[["elapsed"]]
    test <- mgcvST.test(fit, pairs = pairs, calibration = "liu",
      BPPARAM = BP, chunk_size = ceiling(nrow(pairs) / workers))
    pair.seconds <- proc.time()[["elapsed"]] - t0
    saveRDS(test, test.file, compress = FALSE)
  }
  write.csv(test$results, file.path(d.out, "pair-tests.csv"), row.names = FALSE)
  timing[[i]] <- data.frame(config = config, workers = workers,
    feature_fits = length(ids), pairs = nrow(pairs),
    fit_seconds = fit.seconds, pair_seconds = pair.seconds)
  saveRDS(timing[[i]], timing.file)
  fits[[config]] <- fit
  tests[[config]] <- test
  if (config != "serial") BP <- bpstop(BP)
}
write.csv(do.call(rbind, timing), file.path(out, "timing.csv"), row.names = FALSE)

ref <- fits$serial
feature.comparison <- list()
for (config in c("snow2", "snow4")) {
  fit <- fits[[config]]
  if (!identical(ref$feature_id, fit$feature_id)) {
    stop("Feature order differs for ", config, ".")
  }
  for (j in seq_along(ids)) {
    ref.size <- as.numeric(ref$family_parameters[[j]][1L])
    fit.size <- as.numeric(fit$family_parameters[[j]][1L])
    ref.u <- ref$smooth_coefficients$global[j, ]
    fit.u <- fit$smooth_coefficients$global[j, ]
    feature.comparison[[length(feature.comparison) + 1L]] <- data.frame(
      config = config, feature_id = ids[j],
      converged_serial = ref$diagnostics$converged[j],
      converged_parallel = fit$diagnostics$converged[j],
      max_abs_u_difference = max(abs(fit.u - ref.u)),
      max_abs_working_error_difference = max(abs(
        fit$working_error[, j] - ref$working_error[, j])),
      max_abs_working_variance_difference = max(abs(
        fit$working_variance[, j] - ref$working_variance[, j])),
      max_abs_nuisance_covariance_difference = max(abs(
        fit$nuisance_covariance[[j]] - ref$nuisance_covariance[[j]])),
      lambda_difference = fit$lambda[j] - ref$lambda[j],
      dispersion_difference = fit$dispersion[j] - ref$dispersion[j],
      nb_size_difference = fit.size - ref.size,
      observation_spatial_mean_difference = fit$observation_spatial_mean[j] -
        ref$observation_spatial_mean[j],
      marginal_p_difference = fit$diagnostics$marginal_p_value[j] -
        ref$diagnostics$marginal_p_value[j])
  }
}
write.csv(do.call(rbind, feature.comparison), file.path(out,
  "feature-comparison.csv"), row.names = FALSE)

pair.comparison <- list()
ref.pair <- tests$serial$results
for (config in c("snow2", "snow4")) {
  z <- tests[[config]]$results
  if (!identical(ref.pair[, c("feature1", "feature2")],
      z[, c("feature1", "feature2")])) {
    stop("Pair order differs for ", config, ".")
  }
  pair.comparison[[config]] <- data.frame(config = config,
    feature1 = z$feature1, feature2 = z$feature2,
    signed_score_difference = z$signed_score - ref.pair$signed_score,
    information_difference = z$information - ref.pair$information,
    p_two_sided_difference = z$p_two_sided - ref.pair$p_two_sided,
    p_positive_difference = z$p_positive - ref.pair$p_positive,
    p_negative_difference = z$p_negative - ref.pair$p_negative,
    available_serial = is.finite(ref.pair$p_two_sided),
    available_parallel = is.finite(z$p_two_sided))
}
write.csv(do.call(rbind, pair.comparison), file.path(out,
  "pair-comparison.csv"), row.names = FALSE)
feature.comparison <- do.call(rbind, feature.comparison)
pair.comparison <- do.call(rbind, pair.comparison)
summary <- list()
for (config in c("snow2", "snow4")) {
  F <- feature.comparison[feature.comparison$config == config, , drop = FALSE]
  P <- pair.comparison[pair.comparison$config == config, , drop = FALSE]
  summary[[config]] <- data.frame(config = config, features = nrow(F),
    pairs = nrow(P), converged_parallel = sum(F$converged_parallel),
    pair_availability_mismatches = sum(P$available_serial != P$available_parallel),
    max_abs_u_difference = max(abs(F$max_abs_u_difference)),
    max_abs_working_error_difference = max(abs(F$max_abs_working_error_difference)),
    max_abs_working_variance_difference = max(abs(F$max_abs_working_variance_difference)),
    max_abs_lambda_difference = max(abs(F$lambda_difference)),
    max_abs_pair_score_difference = max(abs(P$signed_score_difference)),
    max_abs_pair_information_difference = max(abs(P$information_difference)),
    max_abs_pair_p_difference = max(abs(P$p_two_sided_difference)))
}
write.csv(do.call(rbind, summary), file.path(out, "comparison-summary.csv"),
  row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out, "session-info.txt"))
