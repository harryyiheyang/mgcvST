options(stringsAsFactors = FALSE)
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1", RCPP_PARALLEL_NUM_THREADS = "1")
lib <- normalizePath(Sys.getenv("MGCVST_VALIDATION_LIBRARY",
  "artifacts/inla-bam-validation/library"), mustWork = TRUE)
.libPaths(c(lib, .libPaths()))
library(mgcvST)
library(mgcv)

out <- Sys.getenv("MGCVST_SCALING_OUTPUT", "artifacts/inla-bam-validation/scaling")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
M <- expand.grid(replicate = 1:10, n = c(2000L, 8000L, 32000L))
M$seed <- 950000L + M$n + M$replicate
M$task <- seq_len(nrow(M))
write.csv(M, file.path(out, "manifest.csv"), row.names = FALSE)
args <- commandArgs(trailingOnly = TRUE)
tasks <- if (length(args)) as.integer(sub("^--task=", "", args)) else M$task
stopifnot(all(tasks %in% M$task))
flat <- list(prior = "flat", param = numeric(), initial = 0)
ctl <- list(precision_prior = flat, nb_size_prior = flat)
bp <- BiocParallel::SerialParam()

bam_compact <- function(fits, ids) {
  states <- lapply(fits, function(fit) {
    L <- mgcvST:::.gam_training_lpmatrix(fit)
    H <- mgcvST:::.mgcvst_model_geometry(fit, L)
    N <- mgcvST:::.mgcvst_nuisance_state(fit, H, list(L = L, frozen = TRUE))
    W <- rkhs_extract_working_model(fit)
    list(W = W, H = H, N = N)
  })
  H <- states[[1L]]$H
  H$nuisance_columns <- states[[1L]]$N$columns
  H$nuisance_design <- states[[1L]]$N$design
  H$nuisance_projection <- "conditional_Vp_block"
  sp <- do.call(rbind, lapply(states, function(x) x$H$sp))
  lambda <- setNames(vapply(states, function(x) {
    j <- x$H$target[["global"]]
    x$H$sp[x$H$smooth[[j]]$sp_index]
  }, numeric(1L)), ids)
  structure(list(feature_id = ids,
    working_error = do.call(cbind, lapply(states, function(x) x$W$working_error)),
    working_variance = do.call(cbind, lapply(states, function(x) x$W$working_variance)),
    dispersion = setNames(vapply(states, function(x) x$W$dispersion, numeric(1L)), ids),
    lambda = lambda, component_lambda = matrix(lambda, ncol = 1L,
      dimnames = list(ids, "global")), smoothing_parameters = sp,
    family_parameters = setNames(lapply(states, function(x) x$W$family_parameters), ids),
    nuisance_covariance = setNames(lapply(states, function(x) x$N$covariance), ids),
    geometry = H, row_id = H$row_id, score_components = H$score_components,
    model_setting = "global", test_engine = "single_model",
    diagnostics = data.frame(index = seq_along(ids), feature_id = ids,
      converged = vapply(fits, function(x) isTRUE(x$converged), logical(1L)),
      error_message = NA_character_)),
    class = c("mgcvST_model_fit", "mgcvST_fit", "mgcvST"))
}

for (i in tasks) {
  z <- M[i, ]
  dest <- file.path(out, sprintf("n%d_r%02d", z$n, z$replicate))
  dir.create(dest, recursive = TRUE, showWarnings = FALSE)
  if (file.exists(file.path(dest, "result.rds"))) next
  writeLines("running", file.path(dest, "status.txt"))
  set.seed(z$seed)
  n <- z$n
  t0 <- proc.time()[["elapsed"]]
  d <- data.frame(x = runif(n, .01, .99), y = runif(n, .01, .99))
  d$offset0 <- .35 * (d$x - d$y)
  vertices <- as.matrix(expand.grid(x = seq(0, 1, length.out = 15L),
                                    y = seq(0, 1, length.out = 15L)))
  mesh <- list(loc = vertices, graph = list(tv = geometry::delaunayn(vertices)))
  basis <- spde_basis(mesh, as.matrix(d[c("x", "y")]), kappa = 6,
                       project_intercept = TRUE)
  basis$component <- basis$score.component <- "global"
  F <- basis$B %*% backsolve(chol(basis$Q), diag(ncol(basis$Q)))
  tau <- mean(rowSums(F^2)) / .36
  E <- matrix(rnorm(ncol(F) * 2L), ncol = 2L)
  H <- F %*% E / sqrt(tau)
  v <- rowSums(F^2) / tau
  beta <- log(.3) - log(mean(exp(d$offset0 + .5 * v)))
  Y <- t(vapply(1:2, function(j) rnbinom(n,
    mu = exp(beta + .15 * (j - 1) + d$offset0 + H[, j]), size = 2), numeric(n)))
  ids <- rownames(Y) <- c("feature1", "feature2")
  basis_seconds <- proc.time()[["elapsed"]] - t0
  saveRDS(list(data = d, mesh = mesh, Y = Y, field = H, case = z,
               prior = ctl), file.path(dest, "input.rds"))
  rm(F, E)

  t0 <- proc.time()[["elapsed"]]
  S <- inlaST.set(response ~ offset(offset0), d, basis, family = nb(),
                    control = ctl, score_backend = "sparse")
  inla_setup_seconds <- proc.time()[["elapsed"]] - t0
  t0 <- proc.time()[["elapsed"]]
  I <- inlaST.estimate(Y, S, control = ctl, BPPARAM = bp,
    retain_smooth = TRUE, retain_marginal = TRUE, marginal_args = list(method = "liu"))
  inla_estimator_seconds <- proc.time()[["elapsed"]] - t0
  saveRDS(I, file.path(dest, "inla-fit.rds"))
  if (any(!I$diagnostics$converged) || any(!is.finite(I$working_error))) {
    stop("INLA returned an invalid fit in scaling task ", i, "; checkpoint retained.")
  }
  t0 <- proc.time()[["elapsed"]]
  PI <- inlaST.test(I, pairs = matrix(ids, nrow = 1), calibration = "liu", BPPARAM = bp)
  inla_pair_seconds <- proc.time()[["elapsed"]] - t0

  d$response <- Y[1, ]
  t0 <- proc.time()[["elapsed"]]
  G0 <- bam(response ~ offset(offset0) + s(x, y, bs = "spde", xt = basis),
    data = d, family = nb(), method = "fREML", discrete = TRUE, nthreads = 1L, fit = FALSE)
  bam_setup_seconds <- proc.time()[["elapsed"]] - t0
  ri <- attr(G0$terms, "response")
  fr <- serialize(G0$family, NULL)
  B <- vector("list", 2L)
  t0 <- proc.time()[["elapsed"]]
  for (j in 1:2) {
    G <- G0
    G$y <- Y[j, ]
    G$mf[[ri]] <- Y[j, ]
    G$family <- unserialize(fr)
    B[[j]] <- bam(G = G, method = "fREML", discrete = TRUE, nthreads = 1L)
  }
  bam_fit_seconds <- proc.time()[["elapsed"]] - t0
  t0 <- proc.time()[["elapsed"]]
  MB <- lapply(B, function(b) mgcvST:::taps_score_test(b,
    test.component = 1L, method = "liu", n_threads = 1L))
  bam_marginal_seconds <- proc.time()[["elapsed"]] - t0
  t0 <- proc.time()[["elapsed"]]
  C <- bam_compact(B, ids)
  bam_compact_seconds <- proc.time()[["elapsed"]] - t0
  saveRDS(C, file.path(dest, "bam-fit.rds"))
  if (any(!C$diagnostics$converged) || any(!is.finite(C$working_error))) {
    stop("bam returned an invalid fit in scaling task ", i, "; checkpoint retained.")
  }
  t0 <- proc.time()[["elapsed"]]
  PB <- mgcvST.test(C, pairs = matrix(ids, nrow = 1), calibration = "liu", BPPARAM = bp)
  bam_pair_seconds <- proc.time()[["elapsed"]] - t0
  if (any(!is.finite(c(PI$results$p_two_sided, PB$results$p_two_sided)))) {
    saveRDS(list(inla = PI, bam = PB), file.path(dest, "failed-pair.rds"))
    stop("Invalid pair test in scaling task ", i, "; checkpoint retained.")
  }
  SB <- matrix(NA_real_, n, 2L)
  SI <- I$geometry$smooth[[I$geometry$target[["global"]]]]$B %*%
    t(I$smooth_coefficients$global)
  for (j in 1:2) {
    sm <- B[[j]]$smooth[[1L]]
    L <- mgcvST:::.gam_training_lpmatrix(B[[j]])
    cols <- seq.int(sm$first.para, sm$last.para)
    SB[, j] <- as.numeric(L[, cols, drop = FALSE] %*% coef(B[[j]])[cols])
  }
  R <- data.frame(task = i, n = n, q = ncol(basis$B), replicate = z$replicate,
    seed = z$seed, feature_id = ids, threads = 1L,
    field_correlation = vapply(1:2, function(j) cor(SB[, j], SI[, j]), numeric(1)),
    field_rmse_bam = sqrt(colMeans((SB - H)^2)),
    field_rmse_inla = sqrt(colMeans((SI - H)^2)),
    bam_lambda = C$lambda, inla_lambda = I$lambda,
    bam_marginal_p = vapply(MB, function(x) x$smooth.pvalue, numeric(1)),
    inla_marginal_p = I$diagnostics$marginal_p_value,
    bam_pair_p = PB$results$p_two_sided, inla_pair_p = PI$results$p_two_sided,
    bam_pair_score = PB$results$signed_score, inla_pair_score = PI$results$signed_score,
    basis_seconds = basis_seconds, inla_setup_seconds = inla_setup_seconds,
    inla_estimator_seconds = inla_estimator_seconds,
    inla_marginal_seconds = I$timing$marginal_elapsed,
    inla_pair_seconds = inla_pair_seconds, bam_setup_seconds = bam_setup_seconds,
    bam_fit_seconds = bam_fit_seconds, bam_marginal_seconds = bam_marginal_seconds,
    bam_compact_seconds = bam_compact_seconds, bam_pair_seconds = bam_pair_seconds)
  write.csv(R, file.path(dest, "metrics.csv"), row.names = FALSE)
  saveRDS(list(metrics = R, inla_pair = PI, bam_pair = PB,
    bam_marginal = MB, prior = ctl, package = as.character(packageVersion("mgcvST"))),
    file.path(dest, "result.rds"))
  writeLines("complete", file.path(dest, "status.txt"))
  rm(I, C, B, G0, G, S, basis, SI, SB, L)
  gc()
}
writeLines(capture.output(sessionInfo()), file.path(out, "session-info.txt"))
