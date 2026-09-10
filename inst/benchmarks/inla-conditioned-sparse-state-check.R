#!/usr/bin/env Rscript

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1")
repo <- normalizePath(".", winslash = "/", mustWork = TRUE)
private <- file.path(repo, "artifacts", "flat-test-library")
if (!dir.exists(file.path(private, "mgcvST"))) private <- file.path(repo, ".test-library")
.libPaths(c(private, file.path(repo, ".test-library"), .libPaths()))
suppressPackageStartupMessages(library(mgcvST))
source("inst/benchmarks/inla-conditioned-sparse-state.R")

relative_error <- function(x, y) max(abs(x - y)) / max(1, max(abs(y)))
trace_power <- function(A, k) {
  current <- diag(nrow(A))
  for (j in seq_len(k)) current <- current %*% A
  sum(diag(current))
}
dense_reference <- function(A, Q, Z, e, D, tau, X, Vp) {
  Qp <- crossprod(Z, as.matrix(Q %*% Z))
  R <- chol((Qp + t(Qp)) / 2)
  F <- as.matrix(A %*% Z) %*% backsolve(R, diag(ncol(Z))) / sqrt(tau)
  W <- 1 / D
  WF <- W * F
  woodbury <- chol(diag(ncol(F)) + crossprod(F, WF))
  vsolve <- function(value) {
    WV <- W * value
    WV - WF %*% backsolve(
      woodbury, forwardsolve(t(woodbury), crossprod(F, WV))
    )
  }
  WX <- vsolve(X)
  psolve <- function(value) {
    out <- vsolve(value)
    out - WX %*% Vp %*% crossprod(X, out)
  }
  list(a = as.numeric(crossprod(F, psolve(e))),
       M = crossprod(F, psolve(F)))
}

rows <- list()
position <- 0L
record <- function(case, sparse, reference, sparse_seconds, reference_seconds,
                   production = NULL, production_seconds = NA_real_) {
  sparse_cal <- rkhs_score_calibrate(
    as.numeric(crossprod(sparse$a)), sparse$M, sparse$M, method = "liu"
  )
  reference_cal <- rkhs_score_calibrate(
    as.numeric(crossprod(reference$a)), reference$M, reference$M,
    method = "liu"
  )
  sparse_moments <- vapply(1:4, function(k) trace_power(sparse$M %*% sparse$M, k), numeric(1L))
  reference_moments <- vapply(1:4, function(k) trace_power(reference$M %*% reference$M, k), numeric(1L))
  production_a_error <- production_M_error <- production_p_error <- NA_real_
  if (!is.null(production)) {
    production_a_error <- relative_error(sparse$a, production$a)
    production_M_error <- relative_error(sparse$M, production$M)
    production_cal <- rkhs_score_calibrate(
      as.numeric(crossprod(production$a)), production$M, production$M,
      method = "liu"
    )
    production_p_error <- abs(sparse_cal$p_two_sided - production_cal$p_two_sided)
  }
  position <<- position + 1L
  rows[[position]] <<- data.frame(
    case = case,
    n = sparse$diagnostics$observation_dimension,
    m = sparse$diagnostics$mesh_dimension,
    p = sparse$diagnostics$nuisance_dimension,
    a_relative_error = relative_error(sparse$a, reference$a),
    M_relative_error = relative_error(sparse$M, reference$M),
    score_relative_error = relative_error(
      as.numeric(crossprod(sparse$a)), as.numeric(crossprod(reference$a))
    ),
    four_moment_relative_error = relative_error(sparse_moments, reference_moments),
    p_value_absolute_error = abs(
      sparse_cal$p_two_sided - reference_cal$p_two_sided
    ),
    production_a_relative_error = production_a_error,
    production_M_relative_error = production_M_error,
    production_p_value_absolute_error = production_p_error,
    constraint_solve_error = sparse$diagnostics$constraint_solve_error,
    nuisance_identity_error = sparse$diagnostics$nuisance_identity_error,
    nuisance_identity_relative_error =
      sparse$diagnostics$nuisance_identity_relative_error,
    minimum_M_eigenvalue = sparse$diagnostics$minimum_M_eigenvalue,
    minimum_M_relative_eigenvalue =
      sparse$diagnostics$minimum_M_relative_eigenvalue,
    sparse_seconds = sparse_seconds,
    dense_reference_seconds = reference_seconds,
    current_production_seconds = production_seconds,
    dense_observation_by_mesh_constructed =
      sparse$diagnostics$dense_observation_by_mesh_constructed
  )
}

# Cached low-count working state, m = 36.  Use the retained native posterior Vp.
cached <- readRDS(
  "artifacts/lowcount-investigation/candidates/cache/rep-0001-original.rds"
)
small_fit <- cached$fits[[1L]]
small_spec <- cached$spec
small_random <- small_spec$random[[1L]]
small_Vp <- if (!is.null(cached$posts)) {
  cached$posts[[1L]]$nuisance_covariance
} else small_fit$nuisance_covariance
t0 <- proc.time()[["elapsed"]]
small <- inlast_conditioned_sparse_state(
  small_random$A, small_random$Q, small_random$constraint,
  small_random$projection, small_fit$working_error,
  small_fit$working_variance, small_fit$tau[[1L]],
  small_spec$fixed$X, small_Vp
)
small_seconds <- proc.time()[["elapsed"]] - t0
t0 <- proc.time()[["elapsed"]]
small_ref <- dense_reference(
  small_random$A, small_random$Q, small_random$projection,
  small_fit$working_error, small_fit$working_variance,
  small_fit$tau[[1L]], small_spec$fixed$X, small_Vp
)
small_ref_seconds <- proc.time()[["elapsed"]] - t0
# Reconstruct the compact public fit and exercise the current production state.
dgp <- readRDS("artifacts/constraint-type1/estimated/nb03_pair_k6-dgp.rds")
small_basis <- spde_basis(
  dgp$mesh, as.matrix(dgp$data[c("x", "y")]),
  kappa = 6, project_intercept = TRUE
)
small_model <- inlaST.set(
  response ~ offset(offset0), dgp$data, small_basis, family = mgcv::nb()
)
compact <- list(
  geometry = small_model$geometry,
  dispersion = 1,
  smoothing_parameters = matrix(small_fit$smoothing_parameters, nrow = 1L),
  working_error = matrix(small_fit$working_error, ncol = 1L),
  working_variance = matrix(small_fit$working_variance, ncol = 1L),
  nuisance_covariance = list(small_Vp),
  .mgcvst_fixed_factors = vector("list", length(small_model$geometry$smooth))
)
t0 <- proc.time()[["elapsed"]]
small_production <- mgcvST:::.mgcvst_model_score_state(compact, 1L)
small_production_seconds <- proc.time()[["elapsed"]] - t0
record("cached_lowcount_m36", small, small_ref,
       small_seconds, small_ref_seconds, production = small_production,
       production_seconds = small_production_seconds)

# Real n = 3611, m = 625 geometry with deterministic synthetic working state.
input <- readRDS("artifacts/pathwaylgm-151673/benchmark-input.rds")
real_basis <- spde_basis(
  input$mesh, as.matrix(input$data[c("x", "y")]),
  kappa = 6, project_intercept = TRUE
)
real_model <- inlaST.set(
  response ~ x + offset(offset0), input$data, real_basis, family = mgcv::nb()
)
real_random <- real_model$inla_spec$random[[1L]]
n <- nrow(input$data)
set.seed(150673L)
real_e <- rnorm(n)
real_D <- exp(runif(n, log(0.25), log(3)))
real_X <- real_model$inla_spec$fixed$X
real_tau <- 1.7
# A non-GLS native-like covariance deliberately preserves the P1 mismatch.
W <- 1 / real_D
# Build a controlled non-GLS Vp from the exact marginal information.  The
# preliminary state is excluded from timing and only supplies that p by p block.
preliminary <- inlast_conditioned_sparse_state(
  real_random$A, real_random$Q, real_random$constraint,
  real_random$projection, real_e, real_D, real_tau, real_X,
  solve(crossprod(real_X, W * real_X))
)
real_Vp <- solve(preliminary$audit$nuisance_information) * 0.999
t0 <- proc.time()[["elapsed"]]
real <- inlast_conditioned_sparse_state(
  real_random$A, real_random$Q, real_random$constraint,
  real_random$projection, real_e, real_D, real_tau, real_X, real_Vp
)
real_seconds <- proc.time()[["elapsed"]] - t0
t0 <- proc.time()[["elapsed"]]
real_ref <- dense_reference(
  real_random$A, real_random$Q, real_random$projection,
  real_e, real_D, real_tau, real_X, real_Vp
)
real_ref_seconds <- proc.time()[["elapsed"]] - t0
real_compact <- list(
  geometry = real_model$geometry,
  dispersion = 1,
  smoothing_parameters = matrix(real_tau, nrow = 1L),
  working_error = matrix(real_e, ncol = 1L),
  working_variance = matrix(real_D, ncol = 1L),
  nuisance_covariance = list(real_Vp),
  .mgcvst_fixed_factors = vector("list", length(real_model$geometry$smooth))
)
t0 <- proc.time()[["elapsed"]]
real_production <- mgcvST:::.mgcvst_model_score_state(real_compact, 1L)
real_production_seconds <- proc.time()[["elapsed"]] - t0
record("real_geometry_synthetic_working_m625", real, real_ref,
       real_seconds, real_ref_seconds, production = real_production,
       production_seconds = real_production_seconds)

out <- "artifacts/conditioned-sparse-state"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
checks <- do.call(rbind, rows)
write.csv(checks, file.path(out, "numerical-check.csv"), row.names = FALSE)
saveRDS(list(checks = checks, small = small$diagnostics,
             real = real$diagnostics), file.path(out, "numerical-check.rds"))
capture.output(sessionInfo(), file = file.path(out, "session-info.txt"))
print(checks, row.names = FALSE)
stopifnot(
  max(checks$a_relative_error) < 1e-8,
  max(checks$M_relative_error) < 1e-8,
  max(checks$score_relative_error) < 1e-8,
  max(checks$four_moment_relative_error) < 1e-8,
  max(checks$p_value_absolute_error) < 1e-8,
  !any(checks$dense_observation_by_mesh_constructed)
)
