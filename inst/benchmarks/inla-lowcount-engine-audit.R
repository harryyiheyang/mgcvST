# Fixed-hyperparameter audit of the low-count INLA working state.
# No production code is modified by this script.
.libPaths(c(".test-library", .libPaths()))
suppressPackageStartupMessages({
  library(mgcvST)
  library(mgcv)
  library(Matrix)
})

source("inst/benchmarks/inla-posterior-vp.R")
out <- "artifacts/lowcount-investigation/engine-audit"
dir.create(out, recursive = TRUE, showWarnings = FALSE)

relative <- function(x, y) max(abs(x - y)) / max(1, max(abs(y)))
matrix_norm <- function(x) sqrt(sum(x^2))

projected_covariance <- function(X, B, Q, weight, tau) {
  design <- cbind(X, B)
  penalty <- bdiag(
    Diagonal(ncol(X), 0), tau * Matrix(Q, sparse = TRUE)
  )
  H <- forceSymmetric(crossprod(sqrt(weight) * design) + penalty)
  as.matrix(solve(H))
}

residual_precision_metrics <- function(X, B, Q, tau, variance, Vp) {
  R <- chol(Q)
  F <- B %*% backsolve(R, diag(ncol(Q))) / sqrt(tau)
  V <- diag(variance) + tcrossprod(F)
  W <- solve(V)
  WX <- W %*% X
  P <- W - WX %*% Vp %*% t(WX)
  c(
    PX_max = max(abs(P %*% X)),
    PVP_relative = matrix_norm(P %*% V %*% P - P) /
      max(matrix_norm(P), .Machine$double.eps)
  )
}

direct_nb_mode <- function(y, offset, X, B, Q, tau, size, start) {
  p <- ncol(X)
  design <- cbind(X, B)
  fn <- function(coef) {
    eta <- offset + as.numeric(design %*% coef)
    mu <- exp(eta)
    -sum(dnbinom(y, mu = mu, size = size, log = TRUE)) +
      0.5 * tau * as.numeric(crossprod(coef[-seq_len(p)],
                                      Q %*% coef[-seq_len(p)]))
  }
  gr <- function(coef) {
    eta <- offset + as.numeric(design %*% coef)
    mu <- exp(eta)
    likelihood_score <- size * (y - mu) / (size + mu)
    answer <- -as.numeric(crossprod(design, likelihood_score))
    answer[-seq_len(p)] <- answer[-seq_len(p)] +
      as.numeric(tau * Q %*% coef[-seq_len(p)])
    answer
  }
  fit <- nlminb(
    start, fn, gradient = gr,
    control = list(eval.max = 1000L, iter.max = 1000L,
                   rel.tol = 1e-12, x.tol = 1e-12)
  )
  list(coefficients = fit$par, objective = fit$objective,
       gradient = gr(fit$par), convergence = fit$convergence,
       message = fit$message)
}

configuration_summary_variance <- function(fit) {
  configs <- fit$misc$configs
  q <- nrow(fit$misc$configs$config[[1L]]$Q)
  answer <- rep(NA_real_, q)
  for (j in seq_along(configs$contents$tag)) {
    tag <- configs$contents$tag[j]
    start <- configs$contents$start[j] - configs$mnpred
    len <- configs$contents$length[j]
    if (start < 1L || start > q) next
    index <- start + seq_len(len) - 1L
    if (tag %in% names(fit$summary.random)) {
      answer[index] <- fit$summary.random[[tag]]$sd^2
    } else if (tag %in% rownames(fit$summary.fixed)) {
      answer[index] <- fit$summary.fixed[tag, "sd"]^2
    }
  }
  if (any(!is.finite(answer))) stop("Could not align INLA summary variances.")
  answer
}

dgp <- readRDS("artifacts/constraint-type1/pilot/nb03_pair_k6-dgp.rds")
case <- dgp$case
data <- dgp$data
basis <- spde_basis(
  dgp$mesh, as.matrix(data[c("x", "y")]),
  kappa = case$kappa, project_intercept = TRUE
)
basis$component <- "global"
basis$score.component <- "global"
size <- 2
tau <- dgp$tau_truth
model <- inlaST.set(
  response ~ offset(offset0), data, basis,
  family = mgcv::nb(theta = size)
)
target <- model$geometry$target[["global"]]
X <- model$geometry$X
B <- model$geometry$smooth[[target]]$B
Q <- model$geometry$smooth[[target]]$penalties[[1L]]

# Five deterministic independent pairs from the existing low-count DGP.
responses <- vector("list", 10L)
pair_id <- feature_in_pair <- integer(10L)
variance_truth <- rowSums(dgp$factor_truth^2)
for (pair in 1:5) {
  set.seed(case$seed + pair)
  signal <- dgp$factor_truth %*%
    matrix(rnorm(ncol(dgp$factor_truth) * 2L), ncol = 2L)
  stopifnot(max(abs(colMeans(signal))) < 1e-10)
  beta <- log(case$mean_count) -
    log(mean(exp(model$offset + 0.5 * variance_truth)))
  for (within in 1:2) {
    j <- 2L * pair - 2L + within
    eta <- beta + 0.25 * (within - 1L) + model$offset + signal[, within]
    responses[[j]] <- rnbinom(nrow(data), mu = exp(eta), size = size)
    pair_id[j] <- pair
    feature_in_pair[j] <- within
  }
}

rows <- vector("list", length(responses))
bam_rows <- vector("list", length(responses))
for (j in seq_along(responses)) {
  y <- responses[[j]]
  engine <- mgcvST:::.inlast_fit_feature(
    model$inla_spec, y, offset = model$offset,
    control = list(
      fixed_precision = tau, nb_size = size,
      fixed_effect_precision = 0, keep_fit = TRUE
    )
  )
  posterior <- inlast_posterior_vp(engine$inla, model$inla_spec)
  config_q <- nrow(engine$inla$misc$configs$config[[1L]]$Q)
  posterior_full <- inlast_posterior_covariance_selected(
    engine$inla, seq_len(config_q)
  )
  summary_variance <- configuration_summary_variance(engine$inla)

  start <- c(engine$fixed_mode, engine$coefficients$global)
  direct <- direct_nb_mode(
    y, model$offset, X, B, Q, tau, size, start
  )
  eta_direct <- model$offset + as.numeric(cbind(X, B) %*% direct$coefficients)
  mu_direct <- exp(eta_direct)
  expected_weight <- 1 / engine$working_variance
  observed_weight <- size * engine$mu * (size + y) /
    (size + engine$mu)^2
  covariance_expected <- projected_covariance(X, B, Q, expected_weight, tau)
  covariance_observed <- projected_covariance(X, B, Q, observed_weight, tau)
  expected_vp <- covariance_expected[seq_len(ncol(X)), seq_len(ncol(X)),
                                     drop = FALSE]
  observed_vp <- covariance_observed[seq_len(ncol(X)), seq_len(ncol(X)),
                                     drop = FALSE]
  p_expected_expected <- residual_precision_metrics(
    X, B, Q, tau, engine$working_variance, expected_vp
  )
  p_posterior_expected <- residual_precision_metrics(
    X, B, Q, tau, engine$working_variance,
    posterior$nuisance_covariance
  )
  p_posterior_observed <- residual_precision_metrics(
    X, B, Q, tau, 1 / observed_weight,
    posterior$nuisance_covariance
  )
  rows[[j]] <- data.frame(
    pair = pair_id[j], feature_in_pair = feature_in_pair[j],
    mean_y = mean(y), zero_fraction = mean(y == 0), tau = tau, nb_size = size,
    inla_converged = engine$converged,
    direct_convergence = direct$convergence,
    mode_max_abs_difference = max(abs(start - direct$coefficients)),
    direct_gradient_max_abs = max(abs(direct$gradient)),
    eta_max_abs_difference = max(abs(engine$eta - eta_direct)),
    mu_max_relative_difference = relative(engine$mu, mu_direct),
    working_D_formula_error = max(abs(
      engine$working_variance - (1 / engine$mu + 1 / size)
    )),
    posterior_vp = posterior$nuisance_covariance[1L, 1L],
    expected_vp = expected_vp[1L, 1L],
    observed_vp = observed_vp[1L, 1L],
    posterior_vs_expected_relative = relative(
      posterior$nuisance_covariance, expected_vp
    ),
    posterior_vs_observed_relative = relative(
      posterior$nuisance_covariance, observed_vp
    ),
    full_posterior_sd2_max_abs_error = max(abs(
      diag(posterior_full$covariance) - summary_variance
    )),
    constraint_covariance_error =
      posterior_full$diagnostics$constraint_covariance_error,
    expected_PX = p_expected_expected[["PX_max"]],
    expected_PVP_relative = p_expected_expected[["PVP_relative"]],
    posterior_expectedD_PX = p_posterior_expected[["PX_max"]],
    posterior_expectedD_PVP_relative =
      p_posterior_expected[["PVP_relative"]],
    posterior_observedD_PX = p_posterior_observed[["PX_max"]],
    posterior_observedD_PVP_relative =
      p_posterior_observed[["PVP_relative"]],
    stringsAsFactors = FALSE
  )

  # Compare bam's conditional Vp with an expected-Fisher reconstruction at
  # bam's own fitted smoothing parameter and working state.
  bam_data <- data
  bam_data$response <- y
  bam_fit <- mgcv::bam(
    response ~ offset(offset0) + s(x, y, bs = "spde", xt = basis),
    data = bam_data, family = mgcv::nb(theta = size), method = "fREML",
    discrete = TRUE, nthreads = 1L
  )
  L_bam <- mgcvST:::.gam_training_lpmatrix(bam_fit)
  geometry_bam <- mgcvST:::.mgcvst_model_geometry(bam_fit, L_bam)
  target_bam <- geometry_bam$target[["global"]]
  X_bam <- geometry_bam$X
  B_bam <- geometry_bam$smooth[[target_bam]]$B
  Q_bam <- geometry_bam$smooth[[target_bam]]$penalties[[1L]]
  tau_bam <- geometry_bam$sp[geometry_bam$smooth[[target_bam]]$sp_index]
  working_bam <- rkhs_extract_working_model(bam_fit)
  nuisance_columns <- setdiff(
    seq_len(ncol(L_bam)), geometry_bam$smooth[[target_bam]]$columns
  )
  bam_vp <- as.matrix(bam_fit$Vp[nuisance_columns, nuisance_columns, drop = FALSE])
  bam_expected_covariance <- projected_covariance(
    X_bam, B_bam, Q_bam, 1 / working_bam$working_variance, tau_bam
  )
  bam_expected_vp <- bam_expected_covariance[
    seq_len(ncol(X_bam)), seq_len(ncol(X_bam)), drop = FALSE
  ]
  bam_p_metrics <- residual_precision_metrics(
    X_bam, B_bam, Q_bam, tau_bam,
    working_bam$working_variance, bam_vp
  )
  bam_expected_metrics <- residual_precision_metrics(
    X_bam, B_bam, Q_bam, tau_bam,
    working_bam$working_variance, bam_expected_vp
  )
  bam_rows[[j]] <- data.frame(
    pair = pair_id[j], feature_in_pair = feature_in_pair[j],
    mean_y = mean(y), tau_bam = tau_bam,
    bam_vp = bam_vp[1L, 1L], bam_expected_vp = bam_expected_vp[1L, 1L],
    bam_vp_expected_relative = relative(bam_vp, bam_expected_vp),
    bam_vp_PX = bam_p_metrics[["PX_max"]],
    bam_vp_PVP_relative = bam_p_metrics[["PVP_relative"]],
    bam_expected_PX = bam_expected_metrics[["PX_max"]],
    bam_expected_PVP_relative = bam_expected_metrics[["PVP_relative"]],
    stringsAsFactors = FALSE
  )
}

engine_audit <- do.call(rbind, rows)
bam_audit <- do.call(rbind, bam_rows)
write.csv(engine_audit, file.path(out, "fixed-hyper-inla.csv"), row.names = FALSE)
write.csv(bam_audit, file.path(out, "bam-vp-vs-expected.csv"), row.names = FALSE)
writeLines(c(
  paste("INLA", as.character(packageVersion("INLA"))),
  paste("mgcv", as.character(packageVersion("mgcv"))),
  "INLA nbinomial internal theta: log(size); from.theta: exp(theta).",
  "NB variance: mu + mu^2/size; expected IRLS D: 1/mu + 1/size.",
  "Every fitted spatial field uses g=A'1/n and g'u=0.",
  "All fixed-hyperparameter fits use tau=dgp$tau_truth and size=2.",
  "No PC prior is used; any estimated positive parameter in the engine retains log(parameter)~N(0,9).",
  "config$Q is the conditional Gaussian latent precision at the single EB configuration.",
  "config$Qinv is a sparse selected inverse and is not used as a dense covariance.",
  "No result in this directory alone proves the cause of type-I inflation."
), file.path(out, "contract.txt"))
saveRDS(
  list(engine = engine_audit, bam = bam_audit,
       session = capture.output(sessionInfo())),
  file.path(out, "engine-audit.rds"), compress = "xz"
)
print(engine_audit)
print(bam_audit)
