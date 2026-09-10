# Run from an installed mgcvST package with INLA available.
# Optionally set MGCVST_INLA_OUTPUT and MGCVST_INLA_DATASETS (comma separated).
suppressPackageStartupMessages(library(mgcvST))
out <- Sys.getenv("MGCVST_INLA_OUTPUT", "inla-estimator-results")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
datasets <- strsplit(Sys.getenv("MGCVST_INLA_DATASETS", "MISO_E13,Visium_B"),
                     ",", fixed = TRUE)[[1L]]
stopifnot(all(datasets %in% c("MISO_E13", "Visium_B")))
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1")
results <- list()
for (dataset in datasets) {
  e <- new.env()
  utils::data(list = dataset, package = "mgcvST", envir = e)
  slice <- e[[dataset]]
  d <- slice$covariates
  Y <- t(as.matrix(slice$expression))
  t0 <- proc.time()[["elapsed"]]
  basis <- spde_basis(slice$meshes$spde, as.matrix(d[c("x", "y")]),
                      kappa = 0.1, project_intercept = TRUE)
  model <- inlaST.set(
    response ~ offset(offset0) + s(x, y, bs = "spde", xt = basis),
    d, family = mgcv::nb(),
    control = list(control.inla = list(tolerance = 5e-4))
  )
  setup_seconds <- proc.time()[["elapsed"]] - t0
  fit <- inlaST.estimate(
    Y, model, retain_smooth = TRUE, retain_marginal = TRUE,
    marginal_args = list(method = "liu"),
    control = list(control.inla = list(tolerance = 1e-4)),
    BPPARAM = BiocParallel::SerialParam()
  )
  if (!all(fit$diagnostics$converged)) {
    print(fit$diagnostics)
    stop("An INLA feature did not converge for ", dataset)
  }
  pairs <- t(utils::combn(rownames(Y), 2L))
  tested <- inlaST.test(fit, pairs = pairs, calibration = "liu", threads = 1L)
  recalibrated <- inlaST.marginal(
    fit, calibration = "liu", BPPARAM = BiocParallel::SerialParam()
  )
  stopifnot(all(is.finite(tested$results$p_two_sided)),
            all(is.finite(fit$diagnostics$marginal_p_value)),
            isTRUE(all.equal(recalibrated$p_value,
                              fit$diagnostics$marginal_p_value,
                              check.attributes = FALSE)))
  B <- fit$geometry$smooth[[fit$geometry$target[["global"]]]]$B
  means <- as.numeric(fit$smooth_coefficients$global %*% colMeans(B))
  stopifnot(max(abs(means)) < 1e-9)
  utils::write.csv(fit$diagnostics, file.path(out, paste0(dataset, "-features.csv")),
                   row.names = FALSE)
  utils::write.csv(tested$results, file.path(out, paste0(dataset, "-pairs.csv")),
                   row.names = FALSE)
  saveRDS(list(model = model, fit = fit, test = tested),
          file.path(out, paste0(dataset, "-inla.rds")))
  results[[dataset]] <- data.frame(
    dataset = dataset, observations = ncol(Y), features = nrow(Y),
    mesh_vertices = basis$raw_dimension, setup_seconds = setup_seconds,
    estimate_seconds = fit$timing$elapsed,
    max_abs_observation_mean = max(abs(means)),
    valid_pair_p_values = sum(is.finite(tested$results$p_two_sided))
  )
  print(results[[dataset]])
}
utils::write.csv(do.call(rbind, results), file.path(out, "acceptance-summary.csv"),
                 row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out, "session-info.txt"))
