#!/usr/bin/env Rscript

# Experimental score construction from the full INLA conditional Gaussian
# latent posterior.  This is research code and does not modify production.

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

out <- Sys.getenv(
  "MGCVST_POSTERIOR_SCORE_OUTPUT",
  "artifacts/lowcount-investigation/posterior-score"
)
dir.create(out, recursive = TRUE, showWarnings = FALSE)

posterior_helper_file <- "inst/benchmarks/inla-posterior-vp.R"
operator_file <- "artifacts/constraint-type1/estimated/frozen/inla-constraint-operators.R"
script_file <- "inst/benchmarks/inla-posterior-score.R"
stopifnot(file.exists(posterior_helper_file), file.exists(operator_file),
          file.exists(script_file))
file.copy(script_file, file.path(out, "frozen-inla-posterior-score.R"),
          overwrite = TRUE)
file.copy(posterior_helper_file,
          file.path(out, "frozen-inla-posterior-vp.R"), overwrite = TRUE)
file.copy(operator_file,
          file.path(out, "frozen-inla-constraint-operators.R"), overwrite = TRUE)

posterior_tools <- new.env(parent = asNamespace("mgcvST"))
sys.source(file.path(out, "frozen-inla-posterior-vp.R"),
           envir = posterior_tools)
reference_tools <- new.env(parent = asNamespace("mgcvST"))
sys.source(file.path(out, "frozen-inla-constraint-operators.R"),
           envir = reference_tools)

relative_error <- function(x, y) {
  max(abs(x - y)) / max(1, max(abs(y)))
}

posterior_random_covariance <- function(inla_fit, block, block_size) {
  config <- posterior_tools$.inlast_vp_configuration(inla_fit)
  tag <- paste0(".inlast_r", block)
  index <- vapply(seq_len(block_size), function(within) {
    posterior_tools$.inlast_vp_tag_index(
      config$configs, tag, within, nrow(config$Q)
    )
  }, integer(1L))
  posterior_tools$inlast_posterior_covariance_selected(
    inla_fit, index, labels = paste0(tag, "[", seq_len(block_size), "]")
  )
}

posterior_score_state <- function(fit, spec, block = 1L) {
  stopifnot(!is.null(fit$inla), length(fit$tau) >= block)
  random <- spec$random[[block]]
  A <- as.matrix(random$A)
  Q <- as.matrix(random$Q)
  Z <- as.matrix(random$projection)
  g <- as.numeric(random$constraint)
  uhat <- as.numeric(fit$random_mode[[block]])
  tau <- as.numeric(fit$tau[block])
  m <- ncol(A)
  stopifnot(length(uhat) == m, all(dim(Q) == c(m, m)),
            all(dim(Z) == c(m, m - 1L)), abs(sum(g) - 1) < 1e-10)

  full <- posterior_random_covariance(fit$inla, block, m)
  Vp_u <- full$covariance
  Qc <- crossprod(Z, Q %*% Z)
  Rc <- chol(Qc)
  gamma_hat <- as.numeric(crossprod(Z, uhat))
  a <- as.numeric(sqrt(tau) * Rc %*% gamma_hat)
  Vp_gamma <- crossprod(Z, Vp_u %*% Z)
  Cxi <- tau * Rc %*% Vp_gamma %*% t(Rc)
  Cxi <- (Cxi + t(Cxi)) / 2
  M <- diag(ncol(Z)) - Cxi
  M <- (M + t(M)) / 2

  # Innovation-to-centered-raw transform.  H maps raw mesh coefficients into
  # g' u = 0 along the constant mesh direction because g'1 = 1.
  Rraw <- chol(Q)
  Rraw_inverse <- backsolve(Rraw, diag(m))
  H <- diag(m) - tcrossprod(rep.int(1, m), g)
  T <- Rc %*% crossprod(Z, H %*% Rraw_inverse)
  a_raw <- as.numeric(crossprod(T, a))
  M_raw <- crossprod(T, M %*% T)
  M_raw <- (M_raw + t(M_raw)) / 2

  F_projected <- (A %*% Z) %*% backsolve(Rc, diag(ncol(Z))) /
    sqrt(tau)
  F_raw_no_Z <- A %*% Rraw_inverse
  F_raw_no_Z <- sweep(F_raw_no_Z, 2L, colMeans(F_raw_no_Z), "-") /
    sqrt(tau)
  F_raw_from_T <- F_projected %*% T

  list(
    a = a, M = M, Cxi = Cxi,
    a_raw = a_raw, M_raw = M_raw, T = T,
    F_projected = F_projected, F_raw_no_Z = F_raw_no_Z,
    posterior_u_covariance = Vp_u,
    diagnostics = list(
      fitted_constraint = abs(sum(g * uhat)),
      covariance_constraint = max(abs(crossprod(g, Vp_u))),
      Cxi_min_eigenvalue = min(eigen(Cxi, symmetric = TRUE,
                                     only.values = TRUE)$values),
      Cxi_max_eigenvalue = max(eigen(Cxi, symmetric = TRUE,
                                     only.values = TRUE)$values),
      M_min_eigenvalue = min(eigen(M, symmetric = TRUE,
                                   only.values = TRUE)$values),
      raw_factor_no_Z_relative_error = relative_error(
        F_raw_from_T, F_raw_no_Z
      ),
      posterior_constraint_covariance_error =
        full$diagnostics$constraint_covariance_error
    )
  )
}

observed_working_state <- function(fit, y, spec, posterior_nuisance_Vp) {
  random <- spec$random[[1L]]
  A <- as.matrix(random$A)
  Q <- as.matrix(random$Q)
  Z <- as.matrix(random$projection)
  X <- as.matrix(spec$fixed$X)
  tau <- as.numeric(fit$tau[1L])
  r <- as.numeric(fit$family_parameters[1L])
  mu <- fit$mu
  eta <- fit$eta
  offset <- spec$offset
  w <- r * mu * (y + r) / (r + mu)^2
  z <- eta + (y - mu) * (r + mu) / (mu * (y + r))
  e <- z - offset

  Rc <- chol(crossprod(Z, Q %*% Z))
  F <- (A %*% Z) %*% backsolve(Rc, diag(ncol(Z))) / sqrt(tau)
  WF <- w * F
  middle <- chol(diag(ncol(F)) + crossprod(F, WF))
  Vsolve <- function(value) {
    WV <- w * value
    WV - WF %*% backsolve(
      middle, forwardsolve(t(middle), crossprod(F, WV))
    )
  }
  WX <- Vsolve(X)
  Vp <- as.matrix(posterior_nuisance_Vp)
  P <- function(value) {
    solved <- Vsolve(value)
    solved - WX %*% Vp %*% crossprod(X, solved)
  }
  Pe <- P(matrix(e, ncol = 1L))
  PF <- P(F)
  a <- as.numeric(crossprod(F, Pe))
  M_naive <- crossprod(F, PF)
  M_naive <- (M_naive + t(M_naive)) / 2
  # V = diag(1/w) + F F'; hence PF' V PF is available without n x n.
  M_sandwich <- crossprod(PF, (1 / w) * PF) +
    crossprod(crossprod(F, PF))
  M_sandwich <- (M_sandwich + t(M_sandwich)) / 2
  list(
    a = a, M_naive = M_naive, M_sandwich = M_sandwich,
    P1 = max(abs(P(X))), minimum_weight = min(w),
    formulas = list(
      weight = "r*mu*(y+r)/(r+mu)^2",
      response = "eta+(y-mu)*(r+mu)/(mu*(y+r))"
    )
  )
}

# -------------------------------------------------------------------------
# Four fixed-hyperparameter Gaussian fits: posterior score versus exact P.
# -------------------------------------------------------------------------
dgp <- readRDS("artifacts/constraint-type1/estimated/nb03_pair_k6-dgp.rds")
d <- dgp$data
n <- nrow(d)
basis <- spde_basis(
  dgp$mesh, as.matrix(d[c("x", "y")]),
  kappa = dgp$case$kappa, project_intercept = TRUE
)
gaussian_model <- inlaST.set(
  response ~ offset(offset0), d, basis, family = gaussian()
)
gaussian_spec <- gaussian_model$inla_spec
gaussian_raw <- gaussian_spec$random[[1L]]
gaussian_control <- list(
  fixed_precision = dgp$tau_truth,
  gaussian_precision = 4,
  keep_fit = TRUE
)

RNGkind("L'Ecuyer-CMRG")
gaussian_rows <- vector("list", 4L)
for (rep in seq_len(4L)) {
  set.seed(92000L + rep)
  signal <- as.numeric(dgp$factor_truth %*%
                         rnorm(ncol(dgp$factor_truth)))
  y <- 1 + d$offset0 + signal + rnorm(n, sd = 0.5)
  fit <- mgcvST:::.inlast_fit_feature(
    gaussian_spec, y, offset = d$offset0, control = gaussian_control
  )
  post <- posterior_score_state(fit, gaussian_spec)
  ref <- reference_tools$inlast_constraint_reference_states(
    gaussian_raw$A, gaussian_raw$Q, gaussian_raw$constraint,
    gaussian_raw$projection,
    working_error = fit$working_error,
    working_variance = fit$working_variance,
    tau = fit$tau[1L], nuisance_X = gaussian_spec$fixed$X,
    nuisance_Vp = fit$nuisance_covariance
  )
  gaussian_rows[[rep]] <- data.frame(
    replicate = rep,
    projected_a_relative_error = relative_error(
      post$a, ref$states$projected$a
    ),
    projected_M_relative_error = relative_error(
      post$M, ref$states$projected$M
    ),
    raw_a_relative_error = relative_error(
      post$a_raw, ref$states$raw_kernel_only$a
    ),
    raw_M_relative_error = relative_error(
      post$M_raw, ref$states$raw_kernel_only$M
    ),
    raw_factor_no_Z_relative_error =
      post$diagnostics$raw_factor_no_Z_relative_error,
    fitted_constraint = post$diagnostics$fitted_constraint,
    covariance_constraint = post$diagnostics$covariance_constraint,
    posterior_constraint_covariance_error =
      post$diagnostics$posterior_constraint_covariance_error,
    Cxi_min_eigenvalue = post$diagnostics$Cxi_min_eigenvalue,
    Cxi_max_eigenvalue = post$diagnostics$Cxi_max_eigenvalue,
    M_min_eigenvalue = post$diagnostics$M_min_eigenvalue
  )
}
gaussian_check <- do.call(rbind, gaussian_rows)
write.csv(gaussian_check, file.path(out, "gaussian-fixed-check.csv"),
          row.names = FALSE)

# -------------------------------------------------------------------------
# Four cached NB fits: posterior score versus expected- and observed-working
# score states.  No new NB fit is run.
# -------------------------------------------------------------------------
cache_dir <- "artifacts/lowcount-investigation/candidates/cache"
nb_rows <- list()
position <- 0L
for (rep in seq_len(2L)) {
  cached <- readRDS(file.path(cache_dir, sprintf("rep-%04d-original.rds", rep)))
  stopifnot(length(cached$fits) == 2L, length(cached$posts) == 2L,
            length(cached$y) == 2L)
  for (feature in seq_len(2L)) {
    fit <- cached$fits[[feature]]
    spec <- cached$spec
    post <- posterior_score_state(fit, spec)
    obs <- observed_working_state(
      fit, cached$y[[feature]], spec,
      cached$posts[[feature]]$nuisance_covariance
    )
    expected <- reference_tools$inlast_constraint_reference_states(
      spec$random[[1L]]$A, spec$random[[1L]]$Q,
      spec$random[[1L]]$constraint, spec$random[[1L]]$projection,
      working_error = fit$working_error,
      working_variance = fit$working_variance,
      tau = fit$tau[1L], nuisance_X = spec$fixed$X,
      nuisance_Vp = fit$nuisance_covariance
    )$states$projected
    position <- position + 1L
    nb_rows[[position]] <- data.frame(
      replicate = rep, feature = feature,
      tau = fit$tau[1L], nb_size = fit$family_parameters[1L],
      posterior_vs_expected_a_relative_error = relative_error(
        post$a, expected$a
      ),
      posterior_vs_expected_M_relative_error = relative_error(
        post$M, expected$M
      ),
      posterior_vs_observed_a_relative_error = relative_error(
        post$a, obs$a
      ),
      posterior_vs_observed_M_naive_relative_error = relative_error(
        post$M, obs$M_naive
      ),
      posterior_vs_observed_M_sandwich_relative_error = relative_error(
        post$M, obs$M_sandwich
      ),
      observed_naive_vs_sandwich_relative_error = relative_error(
        obs$M_naive, obs$M_sandwich
      ),
      observed_P1 = obs$P1,
      observed_minimum_weight = obs$minimum_weight,
      raw_factor_no_Z_relative_error =
        post$diagnostics$raw_factor_no_Z_relative_error,
      fitted_constraint = post$diagnostics$fitted_constraint,
      covariance_constraint = post$diagnostics$covariance_constraint,
      Cxi_min_eigenvalue = post$diagnostics$Cxi_min_eigenvalue,
      Cxi_max_eigenvalue = post$diagnostics$Cxi_max_eigenvalue,
      M_min_eigenvalue = post$diagnostics$M_min_eigenvalue
    )
  }
}
nb_check <- do.call(rbind, nb_rows)
write.csv(nb_check, file.path(out, "nb-observed-check.csv"),
          row.names = FALSE)

summary <- data.frame(
  check = c(
    "gaussian max projected a relative error",
    "gaussian max projected M relative error",
    "gaussian max raw a relative error",
    "gaussian max raw M relative error",
    "gaussian max no-Z raw factor relative error",
    "NB median posterior-vs-expected a relative error",
    "NB median posterior-vs-observed a relative error",
    "NB median posterior-vs-observed naive M relative error",
    "NB median posterior-vs-observed sandwich M relative error",
    "NB max observed P1"
  ),
  value = c(
    max(gaussian_check$projected_a_relative_error),
    max(gaussian_check$projected_M_relative_error),
    max(gaussian_check$raw_a_relative_error),
    max(gaussian_check$raw_M_relative_error),
    max(gaussian_check$raw_factor_no_Z_relative_error),
    median(nb_check$posterior_vs_expected_a_relative_error),
    median(nb_check$posterior_vs_observed_a_relative_error),
    median(nb_check$posterior_vs_observed_M_naive_relative_error),
    median(nb_check$posterior_vs_observed_M_sandwich_relative_error),
    max(nb_check$observed_P1)
  )
)
write.csv(summary, file.path(out, "summary.csv"), row.names = FALSE)

writeLines(c(
  "Research validation of the full-latent-posterior score identity.",
  "Projected innovation: xi=sqrt(tau)*Rc*gamma; a=E(xi|y); M=I-Cov(xi|y).",
  "Full random-node covariance is selected from INLA's single empirical-Bayes Gaussian configuration with its exact constraint.",
  "Gaussian check uses four actual fixed-tau, fixed-noise-precision fits and compares with the exact expected-P reference.",
  "NB check uses four already cached actual fits; it runs no new NB fit.",
  "Observed NB formulas are W=r*mu*(y+r)/(r+mu)^2 and z=eta+(y-mu)*(r+mu)/(mu*(y+r)).",
  "NB comparisons are algebraic diagnostics only and do not establish type-I calibration.",
  "Raw centered score uses H=I-1*g' and can be factored without Z as centered columns of A*Rraw^-1/sqrt(tau)."
), file.path(out, "protocol.txt"))
writeLines(capture.output(sessionInfo()), file.path(out, "session-info.txt"))

print(summary, row.names = FALSE)
print(gaussian_check, row.names = FALSE)
print(nb_check, row.names = FALSE)
