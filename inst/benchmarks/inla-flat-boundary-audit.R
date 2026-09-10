#!/usr/bin/env Rscript

# Bounded audit of an improper flat log-precision prior at the zero-spatial-
# variance boundary.  This script runs only fixed-hyperparameter profiles; it
# does not enable a flat prior in production or claim posterior propriety.

Sys.setenv(
  OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1", RCPP_PARALLEL_NUM_THREADS = "1"
)

refresh_library <- normalizePath(
  Sys.getenv("MGCVST_FLAT_RESEARCH_LIBRARY",
             "artifacts/flat-research-library"), mustWork = TRUE
)
.libPaths(c(refresh_library, .libPaths()))
suppressPackageStartupMessages({
  library(mgcvST)
  library(mgcv)
})

`%||%` <- function(x, y) if (is.null(x)) y else x

out <- Sys.getenv(
  "MGCVST_FLAT_BOUNDARY_OUTPUT",
  "artifacts/flat-prior-investigation/boundary"
)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
script_file <- "inst/benchmarks/inla-flat-boundary-audit.R"
file.copy(script_file, file.path(out, "frozen-inla-flat-boundary-audit.R"),
          overwrite = TRUE)

theta_grid <- c(-10, -6, -2, 2, 6, 10, 14)
tau_grid <- exp(theta_grid)
nb_size <- 2
cache_files <- file.path(
  "artifacts/lowcount-investigation/validation-500/cache",
  sprintf("rep-%04d-original.rds", 1:2)
)
stopifnot(all(file.exists(cache_files)))
caches <- lapply(cache_files, readRDS)
spec <- caches[[1L]]$spec
offset <- as.numeric(spec$offset)
n <- length(offset)
stopifnot(n == 200L, spec$family == "negative_binomial",
          spec$precision_scale_mode %||% "raw" == "raw")

# Two old validation responses, plus one independently generated pure-noise
# response with the same geometry, offset, mean count and NB size.
responses <- list(
  validation_rep1_feature1 = as.numeric(caches[[1L]]$y[[1L]]),
  validation_rep2_feature1 = as.numeric(caches[[2L]]$y[[1L]])
)
RNGkind("L'Ecuyer-CMRG")
set.seed(771901L)
noise_intercept <- log(0.3) - log(mean(exp(offset)))
responses$pure_noise <- rnbinom(
  n, mu = exp(noise_intercept + offset), size = nb_size
)
stopifnot(all(vapply(responses, length, integer(1L)) == n))

# Laplace boundary for the nested no-spatial NB model with one flat intercept.
# Its absolute normalization need not equal INLA's convention; its finite,
# positive value establishes the nonzero likelihood approached as tau -> Inf.
no_spatial_boundary <- function(y) {
  negative_loglik <- function(beta) {
    -sum(dnbinom(y, mu = exp(beta + offset), size = nb_size, log = TRUE))
  }
  opt <- optimize(negative_loglik, interval = c(-20, 10), tol = 1e-12)
  beta <- opt$minimum
  mu <- exp(beta + offset)
  observed_weight <- nb_size * mu * (y + nb_size) / (nb_size + mu)^2
  curvature <- sum(observed_weight)
  stopifnot(is.finite(curvature), curvature > 0)
  data.frame(
    boundary_beta = beta,
    boundary_loglik = -opt$objective,
    boundary_log_laplace_flat_intercept =
      -opt$objective + 0.5 * log(2 * pi) - 0.5 * log(curvature),
    boundary_curvature = curvature
  )
}

# Independent constrained Laplace profile.  Unlike the extreme INLA mlik
# values, this computation remains directly interpretable at the nested
# no-spatial boundary.
X <- as.matrix(spec$fixed$X)
raw <- spec$random[[1L]]
Z <- as.matrix(raw$projection)
B <- as.matrix(raw$A %*% Z)
Qz <- crossprod(Z, as.matrix(raw$Q %*% Z))
p <- ncol(X)
q <- ncol(B)
logdet_Qz <- 2 * sum(log(diag(chol(Qz))))
# Unit average observation-variance representation used by the paired flat
# pilot.  The physical precision returned by the engine is
# precision_scale * tau_internal, so a common physical initial log-precision
# theta has internal initial theta - log(precision_scale).
Qz_R <- chol(Qz)
Fc <- B %*% backsolve(Qz_R, diag(q))
observation_precision_scale <- mean(rowSums(Fc^2))
stopifnot(is.finite(observation_precision_scale),
          observation_precision_scale > 0)
explicit_laplace <- function(y, tau, start = NULL) {
  # Work with v=sqrt(tau)*R*u, where Qz=R'R. This avoids an ill-scaled
  # high-tau Hessian and makes the N(0,I) latent normalizer explicit.
  scaled_basis <- Fc / sqrt(tau)
  design <- cbind(X, scaled_basis)
  objective <- function(coef) {
    eta <- as.numeric(offset + design %*% coef)
    v <- coef[p + seq_len(q)]
    -sum(dnbinom(y, mu = exp(eta), size = nb_size, log = TRUE)) +
      0.5 * drop(crossprod(v))
  }
  gradient <- function(coef) {
    eta <- as.numeric(offset + design %*% coef)
    mu <- exp(eta)
    residual <- (y + nb_size) * mu / (nb_size + mu) - y
    answer <- as.numeric(crossprod(design, residual))
    v <- coef[p + seq_len(q)]
    answer[p + seq_len(q)] <- answer[p + seq_len(q)] +
      v
    answer
  }
  if (is.null(start)) start <- numeric(p + q)
  opt <- optim(
    start, objective, gradient, method = "BFGS",
    control = list(maxit = 2000, reltol = 1e-12)
  )
  coef <- opt$par
  eta <- as.numeric(offset + design %*% coef)
  mu <- exp(eta)
  w <- (y + nb_size) * nb_size * mu / (nb_size + mu)^2
  H <- crossprod(design, w * design)
  index <- p + seq_len(q)
  H[index, index] <- H[index, index] + diag(q)
  logdet_H <- 2 * sum(log(diag(chol((H + t(H)) / 2))))
  v <- coef[index]
  u <- backsolve(Qz_R, v) / sqrt(tau)
  loglik <- sum(dnbinom(y, mu = mu, size = nb_size, log = TRUE))
  log_prior_v <- -0.5 * q * log(2 * pi) - 0.5 * drop(crossprod(v))
  list(
    log_laplace = loglik + log_prior_v +
      0.5 * (p + q) * log(2 * pi) - 0.5 * logdet_H,
    mode = coef, convergence = opt$convergence,
    maximum_gradient = max(abs(gradient(coef))),
    spatial_mode_norm = sqrt(sum(u^2))
  )
}

failed_row <- function(response_id, theta, message, seconds) data.frame(
  response = response_id, theta = theta, tau = exp(theta),
  valid = FALSE, converged = FALSE, mode_status = NA_integer_,
  mlik_integration = NA_real_, mlik_gaussian = NA_real_,
  normal_logtau_logprior = dnorm(theta, 0, 3, log = TRUE),
  normal_prior_objective = NA_real_, flat_logtau_profile = NA_real_,
  nuisance_Vp = NA_real_, expected_nuisance_Vp = NA_real_,
  Vp_ratio = NA_real_, fixed_mode = NA_real_, spatial_mode_norm = NA_real_,
  maximum_spatial_mode = NA_real_, minimum_working_weight = NA_real_,
  maximum_working_weight = NA_real_, mean_mu = NA_real_, mean_y = NA_real_,
  mean_zero_error = NA_real_, constraint_residual = NA_real_,
  fit_seconds = seconds, error = message, stringsAsFactors = FALSE
)

profile_one <- function(response_id, y, theta) {
  started <- proc.time()[["elapsed"]]
  tryCatch({
    fit <- mgcvST:::.inlast_fit_feature(
      spec, y, offset = offset, diagnostics = TRUE,
      control = list(
        fixed_precision = exp(theta), nb_size = nb_size, keep_fit = TRUE
      )
    )
    mlik <- as.numeric(fit$inla$mlik[, 1L])
    posterior_Vp <- as.numeric(fit$nuisance_covariance[1L, 1L])
    expected_Vp <- as.numeric(fit$expected_nuisance_covariance[1L, 1L])
    spatial <- as.numeric(fit$random_mode[[1L]])
    prior <- dnorm(theta, 0, 3, log = TRUE)
    data.frame(
      response = response_id, theta = theta, tau = exp(theta),
      valid = all(is.finite(c(mlik, posterior_Vp, expected_Vp, spatial))),
      converged = isTRUE(fit$converged), mode_status = fit$mode_status,
      mlik_integration = mlik[1L],
      mlik_gaussian = if (length(mlik) > 1L) mlik[2L] else NA_real_,
      normal_logtau_logprior = prior,
      normal_prior_objective = mlik[1L] + prior,
      flat_logtau_profile = mlik[1L],
      nuisance_Vp = posterior_Vp,
      expected_nuisance_Vp = expected_Vp,
      Vp_ratio = posterior_Vp / expected_Vp,
      fixed_mode = as.numeric(fit$fixed_mode[1L]),
      spatial_mode_norm = sqrt(sum(spatial^2)),
      maximum_spatial_mode = max(abs(spatial)),
      minimum_working_weight = min(1 / fit$working_variance),
      maximum_working_weight = max(1 / fit$working_variance),
      mean_mu = mean(fit$mu), mean_y = mean(y),
      mean_zero_error = max(abs(fit$observation_spatial_mean)),
      constraint_residual = max(abs(fit$constraint_residual)),
      fit_seconds = proc.time()[["elapsed"]] - started,
      error = NA_character_, stringsAsFactors = FALSE
    )
  }, error = function(e) failed_row(
    response_id, theta, conditionMessage(e),
    proc.time()[["elapsed"]] - started
  ))
}

reuse_inla <- identical(Sys.getenv("MGCVST_FLAT_REUSE_INLA"), "1")
if (reuse_inla && file.exists(file.path(out, "profile.csv"))) {
  profile <- read.csv(file.path(out, "profile.csv"), stringsAsFactors = FALSE)
} else {
  rows <- list()
  position <- 0L
  for (response_id in names(responses)) {
    for (theta in theta_grid) {
      position <- position + 1L
      rows[[position]] <- profile_one(response_id, responses[[response_id]], theta)
      cat(sprintf("Completed %s at log(tau)=%g\n", response_id, theta))
      flush.console()
    }
  }
  profile <- do.call(rbind, rows)
  write.csv(profile, file.path(out, "profile.csv"), row.names = FALSE)
}

explicit_rows <- list()
position <- 0L
for (response_id in names(responses)) {
  previous <- NULL
  for (theta in theta_grid) {
    result <- explicit_laplace(
      responses[[response_id]], exp(theta), start = previous
    )
    previous <- result$mode
    position <- position + 1L
    explicit_rows[[position]] <- data.frame(
      response = response_id, theta = theta, tau = exp(theta),
      explicit_log_laplace = result$log_laplace,
      convergence = result$convergence,
      maximum_gradient = result$maximum_gradient,
      spatial_mode_norm = result$spatial_mode_norm
    )
  }
}
explicit_profile <- do.call(rbind, explicit_rows)
write.csv(explicit_profile, file.path(out, "explicit-profile.csv"),
          row.names = FALSE)

boundaries <- do.call(rbind, lapply(names(responses), function(response_id) {
  cbind(response = response_id, no_spatial_boundary(responses[[response_id]]))
}))
write.csv(boundaries, file.path(out, "no-spatial-boundary.csv"),
          row.names = FALSE)

summaries <- do.call(rbind, lapply(split(profile, profile$response), function(z) {
  valid <- z$valid & z$converged & is.finite(z$mlik_integration)
  upper <- z[match(c(6, 10, 14), z$theta), ]
  explicit <- explicit_profile[explicit_profile$response == z$response[1L], ]
  explicit_upper <- explicit[match(c(6, 10, 14), explicit$theta), ]
  boundary <- boundaries[
    boundaries$response == z$response[1L],
    "boundary_log_laplace_flat_intercept"
  ]
  data.frame(
    response = z$response[1L], attempted = nrow(z), valid = sum(valid),
    failed = sum(!valid),
    normal_grid_mode_theta = z$theta[which.max(z$normal_prior_objective)],
    flat_grid_mode_theta = z$theta[which.max(z$flat_logtau_profile)],
    mlik_theta14_minus_theta10 =
      z$mlik_integration[z$theta == 14] - z$mlik_integration[z$theta == 10],
    mlik_theta14_minus_theta6 =
      z$mlik_integration[z$theta == 14] - z$mlik_integration[z$theta == 6],
    upper_tail_slope = unname(coef(lm(mlik_integration ~ theta,
                                      data = upper))[2L]),
    explicit_theta14_minus_theta10 =
      explicit$explicit_log_laplace[explicit$theta == 14] -
      explicit$explicit_log_laplace[explicit$theta == 10],
    explicit_theta14_minus_theta6 =
      explicit$explicit_log_laplace[explicit$theta == 14] -
      explicit$explicit_log_laplace[explicit$theta == 6],
    explicit_upper_tail_slope = unname(coef(lm(
      explicit_log_laplace ~ theta, data = explicit_upper
    ))[2L]),
    explicit_theta14_minus_no_spatial_boundary =
      explicit$explicit_log_laplace[explicit$theta == 14] - boundary,
    max_explicit_gradient = max(explicit$maximum_gradient),
    spatial_norm_theta14 = z$spatial_mode_norm[z$theta == 14],
    Vp_theta14 = z$nuisance_Vp[z$theta == 14],
    max_mean_zero_error = max(z$mean_zero_error, na.rm = TRUE),
    max_constraint_residual = max(z$constraint_residual, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}))
write.csv(summaries, file.path(out, "summary.csv"), row.names = FALSE)

# Free flat-log-precision fits from widely separated initial values.  NB size
# remains fixed, so this isolates the spatial flat tail.  Finite returned modes
# are diagnostics of INLA's optimizer, not evidence that the posterior is
# proper.
flat_initials <- c(-10, 0, 10)
representations <- c(raw = 1, observation_scaled = observation_precision_scale)
if (reuse_inla && file.exists(file.path(out, "flat-multistart.csv"))) {
  flat_multistart <- read.csv(file.path(out, "flat-multistart.csv"),
                              stringsAsFactors = FALSE)
} else {
  flat_rows <- list()
  position <- 0L
  for (response_id in names(responses)) {
    y <- responses[[response_id]]
    for (representation in names(representations)) {
      scale <- unname(representations[[representation]])
      fit_spec <- spec
      fit_spec$random[[1L]]$precision_scale <- scale
      for (physical_initial in flat_initials) {
        internal_initial <- physical_initial - log(scale)
        started <- proc.time()[["elapsed"]]
        position <- position + 1L
        flat_rows[[position]] <- tryCatch({
        fit <- mgcvST:::.inlast_fit_feature(
          fit_spec, y, offset = offset,
          control = list(
            nb_size = nb_size, keep_fit = TRUE,
            precision_prior = list(
              prior = "flat", param = numeric(), initial = internal_initial
            )
          )
        )
        data.frame(
          response = response_id, representation = representation,
          precision_scale = scale,
          physical_initial_logtau = physical_initial,
          internal_initial_logtau = internal_initial,
          returned_tau = as.numeric(fit$tau[1L]),
          returned_tau_internal = as.numeric(fit$tau_internal[1L]),
          returned_logtau = log(as.numeric(fit$tau[1L])),
          returned_internal_logtau = log(as.numeric(fit$tau_internal[1L])),
          converged = isTRUE(fit$converged), mode_status = fit$mode_status,
          mlik = as.numeric(fit$inla$mlik[1L, 1L]),
          spatial_mode_norm = sqrt(sum(fit$random_mode[[1L]]^2)),
          nuisance_Vp = as.numeric(fit$nuisance_covariance[1L, 1L]),
          mean_zero_error = max(abs(fit$observation_spatial_mean)),
          constraint_residual = max(abs(fit$constraint_residual)),
          fit_seconds = proc.time()[["elapsed"]] - started,
          error = NA_character_, stringsAsFactors = FALSE
        )
      }, error = function(e) data.frame(
        response = response_id, representation = representation,
        precision_scale = scale,
        physical_initial_logtau = physical_initial,
        internal_initial_logtau = internal_initial,
        returned_tau = NA_real_, returned_tau_internal = NA_real_,
        returned_logtau = NA_real_, returned_internal_logtau = NA_real_,
        converged = FALSE, mode_status = NA_integer_, mlik = NA_real_,
        spatial_mode_norm = NA_real_, nuisance_Vp = NA_real_,
        mean_zero_error = NA_real_, constraint_residual = NA_real_,
        fit_seconds = proc.time()[["elapsed"]] - started,
        error = conditionMessage(e), stringsAsFactors = FALSE
        ))
      }
    }
  }
  flat_multistart <- do.call(rbind, flat_rows)
  write.csv(flat_multistart, file.path(out, "flat-multistart.csv"),
            row.names = FALSE)
}
flat_summary <- do.call(rbind, lapply(
  split(flat_multistart,
        interaction(flat_multistart$response,
                    flat_multistart$representation, drop = TRUE)), function(z) {
    data.frame(
      response = z$response[1L], representation = z$representation[1L],
      precision_scale = z$precision_scale[1L], attempted = nrow(z),
      valid = sum(is.finite(z$returned_tau)),
      converged = sum(z$converged),
      minimum_returned_logtau = min(z$returned_logtau, na.rm = TRUE),
      maximum_returned_logtau = max(z$returned_logtau, na.rm = TRUE),
      returned_logtau_range = diff(range(z$returned_logtau, na.rm = TRUE)),
      maximum_mean_zero_error = max(z$mean_zero_error, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }
))
write.csv(flat_summary, file.path(out, "flat-multistart-summary.csv"),
          row.names = FALSE)

# Compare physically identical raw and scaled fits at matched physical starts.
raw_compare <- flat_multistart[
  flat_multistart$representation == "raw", , drop = FALSE
]
scaled_compare <- flat_multistart[
  flat_multistart$representation == "observation_scaled", , drop = FALSE
]
representation_comparison <- merge(
  raw_compare, scaled_compare,
  by = c("response", "physical_initial_logtau"),
  suffixes = c("_raw", "_scaled"), sort = TRUE
)
representation_comparison$delta_returned_logtau <-
  representation_comparison$returned_logtau_scaled -
  representation_comparison$returned_logtau_raw
representation_comparison$delta_nuisance_Vp <-
  representation_comparison$nuisance_Vp_scaled -
  representation_comparison$nuisance_Vp_raw
representation_comparison$delta_spatial_mode_norm <-
  representation_comparison$spatial_mode_norm_scaled -
  representation_comparison$spatial_mode_norm_raw
representation_comparison$delta_mlik <-
  representation_comparison$mlik_scaled - representation_comparison$mlik_raw
representation_comparison$generic0_scale_constant <-
  -0.5 * q * log(observation_precision_scale)
representation_comparison$delta_mlik_minus_scale_constant <-
  representation_comparison$delta_mlik -
  representation_comparison$generic0_scale_constant
write.csv(representation_comparison,
          file.path(out, "representation-comparison.csv"), row.names = FALSE)

# Deterministic boundary control: with zero offset and y_i=1, beta=0 and u=0
# exactly solve the latent score equations for every tau because the fitted
# spatial contribution has observation mean zero. A flat log(tau) posterior
# still has its nonintegrable no-spatial tail. These few fits show whether the
# numerical optimizer can report an apparently finite mode and status 0 in
# that exact boundary setting.
balanced_y <- rep.int(1, n)
balanced_offset <- numeric(n)
balanced_rows <- list()
position <- 0L
for (theta in c(2, 6, 10, 14)) {
  position <- position + 1L
  balanced_rows[[position]] <- tryCatch({
    fit <- mgcvST:::.inlast_fit_feature(
      spec, balanced_y, balanced_offset,
      list(fixed_precision = exp(theta), nb_size = nb_size, keep_fit = TRUE)
    )
    data.frame(
      fit_type = "fixed_profile", input_logtau = theta,
      returned_logtau = log(fit$tau[1L]), mode_status = fit$mode_status,
      converged = fit$converged, mlik = fit$inla$mlik[1L, 1L],
      spatial_mode_norm = sqrt(sum(fit$random_mode[[1L]]^2)),
      fixed_mode = fit$fixed_mode[1L], nuisance_Vp = fit$nuisance_covariance[1L, 1L],
      mean_zero_error = max(abs(fit$observation_spatial_mean)),
      error = NA_character_
    )
  }, error = function(e) data.frame(
    fit_type = "fixed_profile", input_logtau = theta,
    returned_logtau = NA_real_, mode_status = NA_integer_, converged = FALSE,
    mlik = NA_real_, spatial_mode_norm = NA_real_, fixed_mode = NA_real_,
    nuisance_Vp = NA_real_, mean_zero_error = NA_real_,
    error = conditionMessage(e)
  ))
}
for (initial in flat_initials) {
  position <- position + 1L
  balanced_rows[[position]] <- tryCatch({
    fit <- mgcvST:::.inlast_fit_feature(
      spec, balanced_y, balanced_offset,
      list(
        nb_size = nb_size, keep_fit = TRUE,
        precision_prior = list(prior = "flat", param = numeric(),
                               initial = initial)
      )
    )
    data.frame(
      fit_type = "free_flat", input_logtau = initial,
      returned_logtau = log(fit$tau[1L]), mode_status = fit$mode_status,
      converged = fit$converged, mlik = fit$inla$mlik[1L, 1L],
      spatial_mode_norm = sqrt(sum(fit$random_mode[[1L]]^2)),
      fixed_mode = fit$fixed_mode[1L], nuisance_Vp = fit$nuisance_covariance[1L, 1L],
      mean_zero_error = max(abs(fit$observation_spatial_mean)),
      error = NA_character_
    )
  }, error = function(e) data.frame(
    fit_type = "free_flat", input_logtau = initial,
    returned_logtau = NA_real_, mode_status = NA_integer_, converged = FALSE,
    mlik = NA_real_, spatial_mode_norm = NA_real_, fixed_mode = NA_real_,
    nuisance_Vp = NA_real_, mean_zero_error = NA_real_,
    error = conditionMessage(e)
  ))
}
balanced_boundary <- do.call(rbind, balanced_rows)
write.csv(balanced_boundary, file.path(out, "balanced-boundary.csv"),
          row.names = FALSE)

writeLines(c(
  "Bounded fixed-hyperparameter audit: 3 responses x 7 log(tau) values; NB size fixed at 2.",
  "The two validation responses are cached feature 1 from independent validation replicates 1 and 2.",
  "pure_noise has no spatial field and was generated once with L'Ecuyer-CMRG seed 771901.",
  "All spatial fits retain the exact observation mean-zero constraint and original Q scale.",
  "Earlier explicit-Laplace and Gaussian audits showed that INLA fixed-hyper mlik matches the likelihood profile up to a nearly theta-invariant constant; it does not contain the varying N(0,9) log-tau prior term.",
  "Therefore flat_logtau_profile equals mlik, while normal_prior_objective adds dnorm(theta,0,3,log=TRUE). No normal-prior subtraction is applied.",
  "A flat improper prior on log(tau) is not made proper by a finite grid. As tau -> Inf the spatial field vanishes and mlik approaches the finite no-spatial NB boundary, yielding a nonintegrable constant upper tail.",
  "Finite-grid convergence flags only describe conditional fits. They do not certify a free flat-prior mode or posterior propriety.",
  "Free-flat spatial-only fits use initial log(tau) -10, 0 and 10 with NB size fixed at 2. Their finite convergence flags do not establish posterior propriety.",
  sprintf("The observation-scaled representation uses precision_scale %.17g.", observation_precision_scale),
  "Raw and observation-scaled free fits use the same physical initial log(tau); the scaled internal initial is physical_initial-log(precision_scale). Differences between representations diagnose numerical optimization or generic0 normalization behavior, not a change in the physical model.",
  "The deterministic boundary control has zero offset and y_i=1. Its exact latent mode is beta=0,u=0 at every fixed tau; free-flat fits therefore audit numerical optimizer status at an exact no-spatial boundary.",
  "No production code or production prior was changed.",
  paste("mgcvST", as.character(packageVersion("mgcvST")),
        "INLA", as.character(packageVersion("INLA")))
), file.path(out, "protocol.txt"))
writeLines(capture.output(sessionInfo()), file.path(out, "session-info.txt"))

print(summaries, row.names = FALSE)
print(boundaries, row.names = FALSE)
print(flat_summary, row.names = FALSE)
print(flat_multistart, row.names = FALSE)
print(representation_comparison[, c(
  "response", "physical_initial_logtau", "delta_returned_logtau",
  "delta_mlik", "generic0_scale_constant",
  "delta_mlik_minus_scale_constant"
)], row.names = FALSE)
print(balanced_boundary, row.names = FALSE)
