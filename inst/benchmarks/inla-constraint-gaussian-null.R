#!/usr/bin/env Rscript

# Known-hyperparameter Gaussian null experiment for observation-mean
# constraints and covariance-score kernels. No INLA fit is performed.

options(stringsAsFactors = FALSE)
Sys.setenv(
  OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1",
  BLIS_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1",
  RCPP_PARALLEL_NUM_THREADS = "1"
)

parse_args <- function(x) {
  out <- list(
    pairs = 10000L, block = 500L, seed = 18531L,
    output_dir = file.path("artifacts", "constraint-type1", "gaussian-known")
  )
  for (arg in x) {
    z <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1L]]
    z[1L] <- gsub("-", "_", z[1L], fixed = TRUE)
    if (length(z) != 2L || !(z[1L] %in% names(out))) {
      stop("Unknown argument: ", arg)
    }
    out[[z[1L]]] <- z[2L]
  }
  for (name in c("pairs", "block", "seed")) out[[name]] <- as.integer(out[[name]])
  if (out$pairs < 10000L || out$block < 1L || out$seed < 1L) {
    stop("Require pairs >= 10000 and positive block and seed values.")
  }
  out
}

psd_spectrum <- function(H) {
  d <- eigen((H + t(H)) / 2, symmetric = TRUE, only.values = TRUE)$values
  tolerance <- sqrt(.Machine$double.eps) * max(1, max(abs(d)))
  if (min(d) < -100 * tolerance) stop("Score covariance is not positive semidefinite.")
  sort(d[d > tolerance], decreasing = TRUE)
}

davies_two_sided_tail <- function(q, spectrum) {
  fit <- CompQuadForm::davies(
    q, lambda = c(spectrum / 2, -spectrum / 2), acc = 1e-10
  )
  if (!is.finite(fit$Qq) || fit$Qq <= 0 || fit$Qq > 1 || fit$ifault != 0L) {
    stop("Davies failed: ifault=", fit$ifault, ", Qq=", fit$Qq)
  }
  min(1, 2 * fit$Qq)
}

davies_critical <- function(alpha, spectrum) {
  upper <- max(1, 8 * sqrt(sum(spectrum^2)))
  while (davies_two_sided_tail(upper, spectrum) > alpha) upper <- 2 * upper
  uniroot(
    function(q) davies_two_sided_tail(q, spectrum) - alpha,
    interval = c(0, upper), tol = 1e-9
  )$root
}

liu_pvalue <- function(U, spectrum) {
  mgcvST:::.liu_squared_score_moments(
    abs(U), sum(spectrum^2), sum(spectrum^4),
    sum(spectrum^6), sum(spectrum^8)
  )$p_value
}

binomial_interval <- function(rejected, total, level = 0.95) {
  if (!total) return(c(NA_real_, NA_real_))
  stats::binom.test(rejected, total, conf.level = level)$conf.int[1:2]
}

make_locations <- function(n, layout) {
  if (layout == "uniform") {
    return(cbind(runif(n, 0.02, 0.98), runif(n, 0.02, 0.98)))
  }
  group <- rbinom(n, 1L, 0.22)
  center <- cbind(ifelse(group == 0L, 0.22, 0.78),
                  ifelse(group == 0L, 0.28, 0.72))
  loc <- center + matrix(rnorm(2L * n, sd = 0.11), ncol = 2L)
  matrix(pmin(0.98, pmax(0.02, loc)), ncol = 2L)
}

make_geometry <- function(kappa, layout, n = 120L, mesh_side = 8L,
                          tau = 1) {
  vertices <- as.matrix(expand.grid(
    x = seq(0, 1, length.out = mesh_side),
    y = seq(0, 1, length.out = mesh_side)
  ))
  mesh <- list(loc = vertices, graph = list(tv = geometry::delaunayn(vertices)))
  loc <- make_locations(n, layout)
  basis <- spde_basis(mesh, loc, kappa = kappa, project_intercept = TRUE)
  mesh_info <- mgcvST:::.spde_basis_mesh(mesh)
  scaled <- sweep(loc, 2L, mesh_info$transform$center, "-") /
    mesh_info$transform$scale
  A <- as.matrix(mgcvST:::.spde_basis_project(mesh_info, scaled))
  fem <- mgcvST:::.spde_basis_fem(mesh_info)
  Q <- kappa^4 * fem$M0 + 2 * kappa^2 * fem$M1 + fem$M2
  Q <- as.matrix(Matrix::forceSymmetric(Q))

  g <- as.numeric(crossprod(A, rep(1 / n, n)))
  R <- chol(Q)
  L <- backsolve(R, diag(ncol(Q)))
  raw_factor <- A %*% L / sqrt(tau)
  h <- as.numeric(crossprod(L, g))
  h_unit <- h / sqrt(sum(h^2))
  coefficient_projector <- diag(length(h)) - tcrossprod(h_unit)
  centered_factor <- raw_factor %*% coefficient_projector

  Z <- qr.Q(qr(matrix(g, ncol = 1L)), complete = TRUE)[, -1L, drop = FALSE]
  Q_projected <- crossprod(Z, Q %*% Z)
  projected_factor <- A %*% Z %*%
    backsolve(chol(Q_projected), diag(ncol(Q_projected))) / sqrt(tau)
  covariance_error <- max(abs(
    tcrossprod(centered_factor) - tcrossprod(projected_factor)
  ))
  if (covariance_error > 1e-8) {
    stop("Raw constrained and projected field covariances disagree: ", covariance_error)
  }

  D <- 0.09 * exp(1.15 * (loc[, 1L] - 0.5) -
                          0.85 * (loc[, 2L] - 0.5))
  X <- matrix(1, n, 1L, dimnames = list(NULL, "(Intercept)"))
  true_covariance <- tcrossprod(centered_factor) + diag(D)
  raw_covariance <- tcrossprod(raw_factor) + diag(D)
  true_inverse <- chol2inv(chol(true_covariance))
  raw_inverse <- chol2inv(chol(raw_covariance))
  center_vp <- solve(crossprod(X, true_inverse %*% X))
  raw_vp <- solve(crossprod(X, raw_inverse %*% X))
  P_center <- true_inverse - true_inverse %*% X %*% center_vp %*%
    crossprod(X, true_inverse)
  P_raw_keep <- raw_inverse - raw_inverse %*% X %*% center_vp %*%
    crossprod(X, raw_inverse)
  P_raw_recompute <- raw_inverse - raw_inverse %*% X %*% raw_vp %*%
    crossprod(X, raw_inverse)
  P_center <- (P_center + t(P_center)) / 2
  P_raw_keep <- (P_raw_keep + t(P_raw_keep)) / 2
  P_raw_recompute <- (P_raw_recompute + t(P_raw_recompute)) / 2

  extra <- as.numeric(raw_factor %*% h_unit)
  euclidean_residual <- extra - X %*% solve(crossprod(X), crossprod(X, extra))
  diagnostics <- list(
    constraint_covariance_error = covariance_error,
    extra_direction_residual_norm_ratio =
      sqrt(sum(euclidean_residual^2) / sum(extra^2)),
    extra_direction_field_trace_fraction =
      sum(extra^2) / sum(raw_factor^2),
    extra_direction_to_noise_ratio = mean(extra^2) / mean(D),
    heteroskedasticity_ratio = max(D) / min(D),
    center_intercept_annihilation = max(abs(P_center %*% X)),
    raw_keep_intercept_residual = max(abs(P_raw_keep %*% X)),
    raw_recomputed_intercept_annihilation = max(abs(P_raw_recompute %*% X))
  )
  variants <- list(
    centered_null_centered_kernel = list(P = P_center, F = centered_factor),
    centered_null_raw_kernel = list(P = P_center, F = raw_factor),
    raw_null_raw_kernel_keep_centered_nuisance =
      list(P = P_raw_keep, F = raw_factor),
    raw_null_raw_kernel_recompute_nuisance =
      list(P = P_raw_recompute, F = raw_factor)
  )
  variants <- lapply(variants, function(x) {
    H <- crossprod(x$F, x$P %*% x$F)
    H <- (H + t(H)) / 2
    list(K = crossprod(x$F, x$P), H = H, spectrum = psd_spectrum(H))
  })
  reference <- inlast_constraint_reference_states(
    A, Q, g, Z, working_error = numeric(n), working_variance = D,
    tau = tau, nuisance_X = X, nuisance_Vp = center_vp
  )
  reference_names <- c(
    centered_null_centered_kernel = "raw_constrained",
    centered_null_raw_kernel = "raw_kernel_only",
    raw_null_raw_kernel_keep_centered_nuisance = "raw_full_keep_nuisance",
    raw_null_raw_kernel_recompute_nuisance = "raw_full_recompute_nuisance"
  )
  reference_error <- vapply(names(variants), function(name) {
    state <- reference$states[[reference_names[[name]]]]
    max(abs(variants[[name]]$H - state$M)) /
      max(1, max(abs(state$M)))
  }, numeric(1L))
  if (max(reference_error) > 1e-9) {
    stop("Vectorized MC operators disagree with reference states: ",
         max(reference_error))
  }
  diagnostics$reference_state_relative_error <- max(reference_error)
  diagnostics$rank_one_identity_relative_error <-
    reference$diagnostics$rank_one_identity_relative_error
  diagnostics$qinv_g_node_constant_ratio <-
    reference$diagnostics$qinv_g_node_constant_ratio
  list(
    loc = loc, D = D, X = X, centered_factor = centered_factor,
    true_covariance = true_covariance, variants = variants,
    diagnostics = diagnostics, raw_A = A, raw_Q = Q, g = g, Z = Z,
    tau = tau
  )
}

simulate_scores <- function(geometry, pairs, block, beta = c(1, 2)) {
  n <- nrow(geometry$X)
  root <- t(chol(geometry$true_covariance))
  answer <- lapply(geometry$variants, function(x) numeric(pairs))
  starts <- seq.int(1L, pairs, by = block)
  for (first in starts) {
    index <- first:min(pairs, first + block - 1L)
    count <- length(index)
    y1 <- beta[1L] + root %*% matrix(rnorm(n * count), n, count)
    y2 <- beta[2L] + root %*% matrix(rnorm(n * count), n, count)
    for (name in names(geometry$variants)) {
      K <- geometry$variants[[name]]$K
      a1 <- K %*% y1
      a2 <- K %*% y2
      answer[[name]][index] <- colSums(a1 * a2)
    }
  }
  answer
}

summarize_scores <- function(U, spectrum, scenario, variant,
                             alphas = c(0.05, 0.01, 0.001)) {
  finite <- is.finite(U)
  liu <- rep(NA_real_, length(U))
  liu[finite] <- liu_pvalue(U[finite], spectrum)
  out <- list()
  for (alpha in alphas) {
    valid <- finite & is.finite(liu)
    rejected <- sum(liu[valid] <= alpha)
    interval <- binomial_interval(rejected, sum(valid))
    out[[length(out) + 1L]] <- data.frame(
      scenario = scenario, variant = variant,
      calibration = "liu_conditional", alpha = alpha,
      simulations = length(U), valid = sum(valid), failures = sum(!valid),
      failure_rate = mean(!valid), rejected = rejected,
      rejection_rate = rejected / sum(valid), ci_low = interval[1L],
      ci_high = interval[2L], critical_value = NA_real_
    )

    critical <- tryCatch(davies_critical(alpha, spectrum), error = function(e) e)
    if (inherits(critical, "condition")) {
      valid_exact <- rep(FALSE, length(U))
      rejected_exact <- 0L
      interval_exact <- c(NA_real_, NA_real_)
      critical_value <- NA_real_
    } else {
      valid_exact <- finite
      rejected_exact <- sum(abs(U[valid_exact]) >= critical)
      interval_exact <- binomial_interval(rejected_exact, sum(valid_exact))
      critical_value <- critical
    }
    out[[length(out) + 1L]] <- data.frame(
      scenario = scenario, variant = variant,
      calibration = "davies_exact_spectrum_conditional", alpha = alpha,
      simulations = length(U), valid = sum(valid_exact),
      failures = sum(!valid_exact), failure_rate = mean(!valid_exact),
      rejected = rejected_exact,
      rejection_rate = if (any(valid_exact)) rejected_exact / sum(valid_exact) else NA_real_,
      ci_low = interval_exact[1L], ci_high = interval_exact[2L],
      critical_value = critical_value
    )
  }
  do.call(rbind, out)
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
suppressPackageStartupMessages(library(mgcvST))
if (!requireNamespace("geometry", quietly = TRUE) ||
    !requireNamespace("CompQuadForm", quietly = TRUE)) {
  stop("This benchmark requires geometry and CompQuadForm.")
}
operator_script <- file.path("inst", "benchmarks", "inla-constraint-operators.R")
if (!file.exists(operator_script)) stop("Missing reference operator script: ", operator_script)
source(operator_script, local = TRUE)
dir.create(args$output_dir, recursive = TRUE, showWarnings = FALSE)

scenarios <- expand.grid(
  kappa = c(0.7, 6), layout = c("uniform", "clustered"),
  stringsAsFactors = FALSE
)
all_results <- all_diagnostics <- list()
saved <- list()
t0 <- proc.time()[["elapsed"]]
for (i in seq_len(nrow(scenarios))) {
  set.seed(args$seed + i)
  scenario <- paste0("kappa_", scenarios$kappa[i], "_", scenarios$layout[i])
  geometry <- make_geometry(scenarios$kappa[i], scenarios$layout[i])
  U <- simulate_scores(geometry, args$pairs, args$block)
  for (variant in names(U)) {
    all_results[[length(all_results) + 1L]] <- summarize_scores(
      U[[variant]], geometry$variants[[variant]]$spectrum,
      scenario, variant
    )
  }
  all_diagnostics[[i]] <- data.frame(
    scenario = scenario, kappa = scenarios$kappa[i],
    layout = scenarios$layout[i], observations = nrow(geometry$X),
    mesh_vertices = ncol(geometry$raw_A),
    as.data.frame(geometry$diagnostics)
  )
  saved[[scenario]] <- list(
    U = U, diagnostics = geometry$diagnostics,
    spectra = lapply(geometry$variants, `[[`, "spectrum"),
    loc = geometry$loc, D = geometry$D
  )
  message("Completed ", scenario, ".")
}
elapsed <- proc.time()[["elapsed"]] - t0
results <- do.call(rbind, all_results)
diagnostics <- do.call(rbind, all_diagnostics)
rownames(results) <- rownames(diagnostics) <- NULL
results$target_in_ci <- with(
  results, is.finite(ci_low) & alpha >= ci_low & alpha <= ci_high
)
utils::write.csv(results, file.path(args$output_dir, "rejection-rates.csv"), row.names = FALSE)
utils::write.csv(diagnostics, file.path(args$output_dir, "scenario-diagnostics.csv"), row.names = FALSE)
saveRDS(
  list(
    results = results, diagnostics = diagnostics, simulations = saved,
    pairs_per_scenario = args$pairs, seed = args$seed, elapsed = elapsed,
    definition = paste(
      "Independent centered nonzero SPDE fields; cross-feature covariance rho=0;",
      "known hyperparameters; heteroskedastic Gaussian noise; intercept nuisance."
    )
  ),
  file.path(args$output_dir, "simulation-details.rds")
)

display <- results[, c(
  "scenario", "variant", "calibration", "alpha", "rejection_rate",
  "ci_low", "ci_high", "target_in_ci", "failure_rate"
)]
exact05 <- subset(
  results,
  calibration == "davies_exact_spectrum_conditional" & alpha == 0.05
)
rate_range <- function(variant) {
  value <- exact05$rejection_rate[exact05$variant == variant]
  paste(format(range(value), digits = 3), collapse = " to ")
}
lines <- c(
  "# Known-hyperparameter Gaussian covariance-score null experiment",
  "",
  paste0("Generated: ", format(Sys.time(), tz = "UTC", usetz = TRUE)),
  paste0("Pairs per scenario: ", format(args$pairs, big.mark = ",")),
  paste0("Elapsed seconds: ", format(elapsed, digits = 7)),
  "",
  "Each feature has a nonzero SPDE marginal field constrained to have zero mean at the observed locations. Paired features are independent, so the tested cross-feature covariance is zero. This is not a test of zero marginal spatial variance.",
  "",
  "`centered_null_raw_kernel` removes projection only from the tested kernel while retaining the correct centered null covariance. The two `raw_null` variants remove the constraint from both null covariance and kernel; one keeps the centered nuisance covariance and the other recomputes intercept GLS under the raw null.",
  "",
  "Davies critical values are computed once from each conditional signed spectrum. Liu uses the same conditional H matrices and its four-moment approximation.",
  "",
  paste0("At alpha=0.05, exact-spectrum rejection ranges were: current centered test ",
         rate_range("centered_null_centered_kernel"),
         "; raw kernel with the correct centered null ",
         rate_range("centered_null_raw_kernel"),
         "; raw null retaining the centered nuisance covariance ",
         rate_range("raw_null_raw_kernel_keep_centered_nuisance"),
         "; and raw null with nuisance GLS recomputed ",
         rate_range("raw_null_raw_kernel_recompute_nuisance"), "."),
  "",
  "```",
  paste(capture.output(print(display, row.names = FALSE, digits = 5)), collapse = "\n"),
  "```",
  "",
  "Scenario diagnostics, including the removed rank-one direction after intercept projection and its field/noise strength, are in `scenario-diagnostics.csv`."
)
writeLines(lines, file.path(args$output_dir, "report.md"))
print(display, row.names = FALSE)
cat("Elapsed seconds:", elapsed, "\n")
