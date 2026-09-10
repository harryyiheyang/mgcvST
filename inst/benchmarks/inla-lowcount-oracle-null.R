#!/usr/bin/env Rscript

# Mechanism probe for the low-count NB pair-score inflation.
#
# This run is deliberately small and pre-specified: it replays replicates
# 1:100 of nb03_pair_k6, fixes the latent precision and NB size at their DGP
# values, retains the exact observation mean-zero constraint, and evaluates
# the legacy expected-P raw-kernel Davies score.  It is an ablation, not the
# user-target posterior-Vp method and not a definitive type-I calibration
# study.

Sys.setenv(
  OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1", RCPP_PARALLEL_NUM_THREADS = "1"
)

private_library <- normalizePath(".final-library", mustWork = TRUE)
.libPaths(c(private_library, .libPaths()))
suppressPackageStartupMessages({
  library(mgcvST)
  library(mgcv)
})

reps <- as.integer(Sys.getenv("MGCVST_ORACLE_REPS", "100"))
out <- Sys.getenv(
  "MGCVST_ORACLE_OUTPUT",
  "artifacts/lowcount-investigation/oracle"
)
stopifnot(identical(reps, 100L))
dir.create(out, recursive = TRUE, showWarnings = FALSE)

script_file <- "inst/benchmarks/inla-lowcount-oracle-null.R"
operator_file <- "artifacts/constraint-type1/estimated/frozen/inla-constraint-operators.R"
dgp_file <- "artifacts/constraint-type1/estimated/nb03_pair_k6-dgp.rds"
baseline_file <- "artifacts/constraint-type1/estimated/nb03_pair_k6-replicates.csv"
cache_dir <- "artifacts/lowcount-investigation/swaps/cache"
stopifnot(file.exists(script_file), file.exists(operator_file),
          file.exists(dgp_file), file.exists(baseline_file), dir.exists(cache_dir))
file.copy(script_file, file.path(out, "frozen-inla-lowcount-oracle-null.R"),
          overwrite = TRUE)
file.copy(operator_file, file.path(out, "frozen-inla-constraint-operators.R"),
          overwrite = TRUE)

reference <- new.env(parent = asNamespace("mgcvST"))
sys.source(file.path(out, "frozen-inla-constraint-operators.R"), envir = reference)

dgp <- readRDS(dgp_file)
stopifnot(
  identical(dgp$case$name, "nb03_pair_k6"),
  dgp$case$seed == 61000L,
  dgp$case$mean_count == 0.3,
  dgp$case$nb_size == 2,
  dgp$case$kappa == 6,
  is.finite(dgp$tau_truth), dgp$tau_truth > 0
)

d <- dgp$data
n <- nrow(d)
basis <- spde_basis(
  dgp$mesh, as.matrix(d[c("x", "y")]),
  kappa = dgp$case$kappa, project_intercept = TRUE
)
model <- inlaST.set(
  response ~ offset(offset0), d, basis, family = mgcv::nb()
)
spec <- model$inla_spec
raw <- spec$random[[1L]]
X <- spec$fixed$X

# Even fixed positive hyperparameters pass through the package's required
# normal-prior validation.  No PC prior is introduced by this experiment.
oracle_control <- list(
  fixed_precision = dgp$tau_truth,
  nb_size = dgp$case$nb_size,
  precision_prior = list(prior = "normal", param = c(0, 1 / 9), initial = 0),
  nb_size_prior = list(prior = "normal", param = c(0, 1 / 9), initial = 0)
)
validated_control <- mgcvST:::.inlast_control(oracle_control)
stopifnot(
  validated_control$fixed_precision == dgp$tau_truth,
  validated_control$nb_size == dgp$case$nb_size,
  identical(validated_control$precision_prior$prior, "normal"),
  identical(validated_control$nb_size_prior$prior, "normal")
)

baseline_all <- read.csv(baseline_file, stringsAsFactors = FALSE)
baseline <- baseline_all[
  baseline_all$replicate <= reps &
    baseline_all$variant == "raw_kernel_only" &
    baseline_all$calibration == "davies",
  c("replicate", "p_value", "fallback", "tau1", "tau2",
    "nb_size1", "nb_size2", "minimum_working_weight", "error")
]
baseline <- baseline[order(baseline$replicate), ]
stopifnot(
  nrow(baseline) == reps,
  identical(baseline$replicate, seq_len(reps)),
  all(is.finite(baseline$p_value)),
  !any(nzchar(baseline$error[!is.na(baseline$error)]))
)
names(baseline)[names(baseline) != "replicate"] <- paste0(
  "baseline_", names(baseline)[names(baseline) != "replicate"]
)

# The cache contains a later refit of the original responses.  Its INLA modes
# need only agree closely with the earlier CSV because optimizer-level numeric
# differences are expected; the cached counts themselves define the pairing.
cache_fit <- do.call(rbind, lapply(seq_len(reps), function(rep) {
  cached <- readRDS(file.path(cache_dir, sprintf("rep-%04d.rds", rep)))
  data.frame(
    replicate = rep,
    cache_tau1 = as.numeric(cached$inla[[1L]]$tau[1L]),
    cache_tau2 = as.numeric(cached$inla[[2L]]$tau[1L]),
    cache_nb_size1 = as.numeric(cached$inla[[1L]]$family_parameters[1L]),
    cache_nb_size2 = as.numeric(cached$inla[[2L]]$family_parameters[1L])
  )
}))
cache_check <- merge(baseline, cache_fit, by = "replicate", sort = TRUE)
for (metric in c("tau1", "tau2", "nb_size1", "nb_size2")) {
  old <- cache_check[[paste0("baseline_", metric)]]
  cached <- cache_check[[paste0("cache_", metric)]]
  cache_check[[paste0(metric, "_relative_difference")]] <-
    abs(old - cached) / pmax(abs(old), .Machine$double.eps)
}
cache_difference_columns <- grep(
  "_relative_difference$", names(cache_check), value = TRUE
)
max_cache_relative_difference <- max(
  as.matrix(cache_check[cache_difference_columns]), na.rm = TRUE
)
if (!is.finite(max_cache_relative_difference) ||
    max_cache_relative_difference > 0.05) {
  stop("Cached original responses fail the 5% baseline-fit alignment guard.")
}
write.csv(
  cache_check[c("replicate", cache_difference_columns)],
  file.path(out, "cache-baseline-fit-check.csv"), row.names = FALSE
)

simulate_pair <- function(rep) {
  cached <- readRDS(file.path(cache_dir, sprintf("rep-%04d.rds", rep)))
  stopifnot(identical(cached$replicate, rep), length(cached$y) == 2L)
  signal <- cached$signal
  if (max(abs(colMeans(signal))) > 1e-10) {
    stop("DGP latent field violates the observation mean-zero constraint.")
  }
  field_variance <- rowSums(dgp$factor_truth^2)
  beta <- log(dgp$case$mean_count) -
    log(mean(exp(d$offset0 + 0.5 * field_variance)))
  eta <- vapply(seq_len(2L), function(j) {
    beta + 0.25 * (j - 1L) + d$offset0 + signal[, j]
  }, numeric(n))
  y <- do.call(cbind, lapply(cached$y, as.numeric))
  if (!all(dim(y) == c(n, 2L)) || any(y < 0) || any(y != floor(y))) {
    stop("Cached response is not a valid aligned count matrix.")
  }
  list(signal = signal, eta = eta, y = y)
}

failed_row <- function(rep, message, seconds) data.frame(
  replicate = rep, oracle_p_value = NA_real_, oracle_signed_score = NA_real_,
  oracle_fallback = NA, fit_seconds = seconds,
  fitted_mean_error = NA_real_, constraint_residual = NA_real_,
  minimum_working_weight = NA_real_, tau1 = NA_real_, tau2 = NA_real_,
  nb_size1 = NA_real_, nb_size2 = NA_real_,
  eta_rmse1 = NA_real_, eta_rmse2 = NA_real_,
  mean_y1 = NA_real_, mean_y2 = NA_real_, error = message,
  stringsAsFactors = FALSE
)

run_one <- function(rep) {
  started <- proc.time()[["elapsed"]]
  tryCatch({
    sim <- simulate_pair(rep)
    fits <- lapply(seq_len(2L), function(j) {
      mgcvST:::.inlast_fit_feature(
        spec, sim$y[, j], offset = d$offset0, control = oracle_control
      )
    })
    if (!all(vapply(fits, function(fit) isTRUE(fit$converged), logical(1L)))) {
      stop("INLA convergence flag failed.")
    }
    tau <- vapply(fits, function(fit) as.numeric(fit$tau[1L]), numeric(1L))
    size <- vapply(
      fits, function(fit) as.numeric(fit$family_parameters[1L]), numeric(1L)
    )
    if (max(abs(tau - dgp$tau_truth)) > 1e-14 ||
        max(abs(size - dgp$case$nb_size)) > 1e-14) {
      stop("A fixed oracle hyperparameter changed during fitting.")
    }
    fitted_mean_error <- max(abs(unlist(lapply(
      fits, `[[`, "observation_spatial_mean"
    ))))
    constraint_residual <- max(abs(unlist(lapply(
      fits, `[[`, "constraint_residual"
    ))))
    if (fitted_mean_error > 1e-10 || constraint_residual > 1e-10) {
      stop("The fitted latent field violates its mean-zero constraint.")
    }
    refs <- lapply(fits, function(fit) {
      reference$inlast_constraint_reference_states(
        raw$A, raw$Q, raw$constraint, raw$projection,
        working_error = fit$working_error,
        working_variance = fit$working_variance,
        tau = dgp$tau_truth,
        nuisance_X = X,
        nuisance_Vp = fit$nuisance_covariance
      )
    })
    cal <- reference$inlast_constraint_reference_pair(
      refs[[1L]]$states$raw_kernel_only,
      refs[[2L]]$states$raw_kernel_only,
      method = "davies"
    )
    data.frame(
      replicate = rep,
      oracle_p_value = cal$p_two_sided,
      oracle_signed_score = cal$signed_score,
      oracle_fallback = !is.null(cal$liu_parameters),
      fit_seconds = proc.time()[["elapsed"]] - started,
      fitted_mean_error = fitted_mean_error,
      constraint_residual = constraint_residual,
      minimum_working_weight = min(unlist(lapply(
        fits, function(fit) 1 / fit$working_variance
      ))),
      tau1 = tau[1L], tau2 = tau[2L],
      nb_size1 = size[1L], nb_size2 = size[2L],
      eta_rmse1 = sqrt(mean((fits[[1L]]$eta - sim$eta[, 1L])^2)),
      eta_rmse2 = sqrt(mean((fits[[2L]]$eta - sim$eta[, 2L])^2)),
      mean_y1 = mean(sim$y[, 1L]), mean_y2 = mean(sim$y[, 2L]),
      error = NA_character_, stringsAsFactors = FALSE
    )
  }, error = function(e) failed_row(
    rep, conditionMessage(e), proc.time()[["elapsed"]] - started
  ))
}

rows <- vector("list", reps)
for (rep in seq_len(reps)) {
  rows[[rep]] <- run_one(rep)
  if (rep %% 10L == 0L) {
    cat(sprintf("Completed fixed-hyperparameter oracle pairs 1-%d\n", rep))
    flush.console()
  }
}
oracle <- do.call(rbind, rows)
result <- merge(oracle, baseline, by = "replicate", all.x = TRUE, sort = TRUE)
write.csv(result, file.path(out, "replicates.csv"), row.names = FALSE)

summarize_method <- function(p, fallback, method) {
  valid <- is.finite(p) & p >= 0 & p <= 1
  rejected <- sum(p[valid] < 0.05)
  ci <- if (sum(valid)) stats::binom.test(rejected, sum(valid))$conf.int else
    c(NA_real_, NA_real_)
  data.frame(
    method = method, attempted = length(p), valid = sum(valid),
    failed = sum(!valid), fallback = sum(fallback, na.rm = TRUE),
    rejected_005 = rejected,
    rejection_rate_005 = if (sum(valid)) rejected / sum(valid) else NA_real_,
    ci_lower = ci[1L], ci_upper = ci[2L], stringsAsFactors = FALSE
  )
}
summary <- rbind(
  summarize_method(
    result$oracle_p_value, result$oracle_fallback,
    "INLA fixed true tau and NB size; experimental expected-P raw kernel"
  ),
  summarize_method(
    result$baseline_p_value, result$baseline_fallback,
    "INLA estimated tau and NB size; same first 100 pairs"
  )
)
write.csv(summary, file.path(out, "summary.csv"), row.names = FALSE)

paired_valid <- is.finite(result$oracle_p_value) &
  is.finite(result$baseline_p_value)
oracle_reject <- result$oracle_p_value[paired_valid] < 0.05
baseline_reject <- result$baseline_p_value[paired_valid] < 0.05
paired <- data.frame(
  cell = c("both", "oracle_only", "estimated_only", "neither"),
  count = c(
    sum(oracle_reject & baseline_reject),
    sum(oracle_reject & !baseline_reject),
    sum(!oracle_reject & baseline_reject),
    sum(!oracle_reject & !baseline_reject)
  ),
  paired_valid = sum(paired_valid), alpha = 0.05
)
write.csv(paired, file.path(out, "paired-rejections.csv"), row.names = FALSE)

diagnostics <- data.frame(
  metric = c(
    "tau_truth", "nb_size_truth", "max_fitted_mean_error",
    "max_constraint_residual", "minimum_working_weight",
    "median_eta_rmse1", "median_eta_rmse2", "mean_fit_seconds",
    "max_cache_baseline_fit_relative_difference"
  ),
  value = c(
    dgp$tau_truth, dgp$case$nb_size,
    max(result$fitted_mean_error, na.rm = TRUE),
    max(result$constraint_residual, na.rm = TRUE),
    min(result$minimum_working_weight, na.rm = TRUE),
    median(result$eta_rmse1, na.rm = TRUE),
    median(result$eta_rmse2, na.rm = TRUE),
    mean(result$fit_seconds, na.rm = TRUE),
    max_cache_relative_difference
  )
)
write.csv(diagnostics, file.path(out, "diagnostics.csv"), row.names = FALSE)

writeLines(c(
  "Pre-specified mechanism probe: exactly replicates 1:100 of saved nb03_pair_k6.",
  "Responses and centered latent fields are read from the original Snow-worker cache; later cached refit parameters are checked within 5% of the earlier baseline and numeric differences are retained in cache-baseline-fit-check.csv.",
  "This avoids the Mersenne-Twister versus L'Ecuyer-CMRG mismatch and makes the oracle/baseline comparison genuinely paired.",
  sprintf("Latent precision is fixed at DGP tau %.17g; NB size is fixed at 2.", dgp$tau_truth),
  "Both positive-hyperparameter configurations retain log(parameter) ~ N(0,3^2); no PC prior is used.",
  "Every INLA fit retains and validates the exact observation mean-zero latent constraint.",
  "The score uses raw_kernel_only with the constrained fitted null, expected nuisance covariance, and Davies calibration.",
  "This expected-P construction is an experimental legacy ablation, not the user-target method based on INLA's own posterior fixed-effect covariance.",
  "Baseline is raw_kernel_only/Davies from the same first 100 estimated-hyperparameter replicates.",
  "This 100-pair paired ablation explores mechanism; it cannot establish general type-I error control.",
  "The old bam low-count result is not used as an oracle because its fit$Vp nuisance block is under separate expected-versus-observed audit."
), file.path(out, "protocol.txt"))
writeLines(capture.output(sessionInfo()), file.path(out, "session-info.txt"))

print(summary, row.names = FALSE)
print(paired, row.names = FALSE)
print(diagnostics, row.names = FALSE)
