#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1", RCPP_PARALLEL_NUM_THREADS = "1")

args <- commandArgs(trailingOnly = TRUE)
task <- if (length(args)) as.integer(sub("^--task=", "", args[1L])) else NA_integer_
out <- Sys.getenv("MGCVST_INFERENCE_OUTPUT",
  "artifacts/inla-bam-validation/inference")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

cases <- data.frame(
  case = c("gaussian_spatial_null", "gaussian_pair_null", "gaussian_positive",
    "gaussian_negative", "nb_spatial_null", "nb_pair_null",
    "nb_positive_fixed_nuisance", "nb_low_count_pair_null",
    "nb_low_count_positive", "nb_positive_random_nuisance"),
  family = c(rep("gaussian", 4L), rep("negative_binomial", 6L)),
  spatial = c(FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE),
  rho = c(0, 0, .7, -.7, 0, 0, .7, 0, .7, .7),
  nuisance = c("none", "none", "none", "none", "none", "none", "fixed", "none", "none", "random"),
  intercept = c(0, 0, 0, 0, log(3), log(3), log(3), log(.3), log(.3), log(3)),
  size = c(NA, NA, NA, NA, 4, 4, 4, 2, 2, 4), stringsAsFactors = FALSE)
M <- merge(cases, data.frame(replicate = 1:10), all = TRUE)
M <- M[order(match(M$case, cases$case), M$replicate), ]
M$seed <- 926000L + 100L * match(M$case, cases$case) + M$replicate
M$task <- seq_len(nrow(M))
write.csv(M, file.path(out, "manifest.csv"), row.names = FALSE)
if (is.na(task)) stop("Supply one fixed manifest row as --task=1,...,100.")
if (length(task) != 1L || task < 1L || task > nrow(M)) stop("Invalid --task value.")

lib <- normalizePath(Sys.getenv("MGCVST_VALIDATION_LIBRARY",
  "artifacts/inla-bam-validation/library"), mustWork = TRUE)
.libPaths(c(lib, .libPaths()))
library(mgcvST)
library(mgcv)
if (!requireNamespace("INLA", quietly = TRUE)) stop("INLA is required.")
if (!requireNamespace("geometry", quietly = TRUE)) stop("geometry is required.")

timed <- function(expr) {
  t0 <- proc.time()[["elapsed"]]
  value <- force(expr)
  list(value = value, seconds = proc.time()[["elapsed"]] - t0)
}

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
    lambda = lambda,
    component_lambda = matrix(lambda, ncol = 1L,
      dimnames = list(ids, "global")),
    smoothing_parameters = sp,
    family_parameters = setNames(lapply(states, function(x) x$W$family_parameters), ids),
    nuisance_covariance = setNames(lapply(states, function(x) x$N$covariance), ids),
    geometry = H, row_id = H$row_id, score_components = H$score_components,
    model_setting = "global", test_engine = "single_model",
    diagnostics = data.frame(index = seq_along(ids), feature_id = ids,
      converged = TRUE, error_message = NA_character_),
    timing = list(backend = "bam", workers = 1L)),
    class = c("mgcvST_model_fit", "mgcvST_fit", "mgcvST"))
}

z <- M[task, ]
set.seed(z$seed)
n <- 200L
vertices <- as.matrix(expand.grid(x = seq(0, 1, length.out = 6L),
  y = seq(0, 1, length.out = 6L)))
mesh <- list(loc = vertices, graph = list(tv = geometry::delaunayn(vertices)))
dat <- data.frame(x = runif(n, .01, .99), y = runif(n, .01, .99),
  z = rnorm(n), offset0 = runif(n, -.2, .2))
basis <- spde_basis(mesh, as.matrix(dat[c("x", "y")]), kappa = 5,
  project_intercept = TRUE)
basis$component <- basis$score.component <- "global"
F <- basis$B %*% backsolve(chol(basis$Q), diag(ncol(basis$Q)))
F <- F / sqrt(mean(rowSums(F^2)))
u <- matrix(rnorm(ncol(F) * 3L), ncol = 3L)
if (z$rho != 0) u[, 2L] <- z$rho * u[, 1L] + sqrt(1 - z$rho^2) * u[, 2L]
eta <- z$intercept + dat$offset0 + if (z$nuisance != "none") .35 * dat$z else 0
eta <- sweep(.45 * F %*% u, 1L, eta, "+")
if (!z$spatial) eta <- matrix(z$intercept + dat$offset0 +
  if (z$nuisance != "none") .35 * dat$z else 0, n, 3L)
if (z$family == "gaussian") {
  Y <- t(vapply(1:3, function(j) eta[, j] + rnorm(n, sd = .7), numeric(n)))
  fam <- gaussian()
} else {
  Y <- t(vapply(1:3, function(j) rnbinom(n, mu = exp(eta[, j]), size = z$size), numeric(n)))
  fam <- nb()
}
ids <- paste0("feature", 1:3)
rownames(Y) <- ids
form <- if (z$nuisance == "fixed") response ~ z + offset(offset0) else if
  (z$nuisance == "random") response ~ s(z, bs = "cr", k = 5) + offset(offset0) else
  response ~ offset(offset0)
flat <- list(prior = "flat", param = numeric(), initial = 0)
ctl <- list(precision_prior = flat, gaussian_precision_prior = flat,
  nb_size_prior = flat)
bp <- BiocParallel::SerialParam()

S <- timed(inlaST.set(form, dat, basis, family = fam))
I <- timed(inlaST.estimate(Y, S$value, retain_smooth = TRUE,
  retain_marginal = TRUE, BPPARAM = bp, control = ctl,
  marginal_args = list(method = "liu")))

dat$response <- Y[1L, ]
bform <- update(form, . ~ . + s(x, y, bs = "spde", xt = basis))
G0 <- timed(bam(bform, data = dat, family = fam, method = "fREML",
  discrete = TRUE, nthreads = 1L, fit = FALSE))
ri <- attr(G0$value$terms, "response")
fr <- serialize(G0$value$family, NULL)
B <- timed(lapply(1:3, function(j) {
  G <- G0$value
  G$y <- Y[j, ]
  G$mf[[ri]] <- Y[j, ]
  G$family <- unserialize(fr)
  bam(G = G, method = "fREML", discrete = TRUE, nthreads = 1L)
}))
C <- timed(bam_compact(B$value, ids))

pairs <- t(combn(ids, 2L))
rows <- list()
fits <- list(inla = I$value, bam = C$value)
for (method in names(fits)) for (cal in c("liu", "davies")) {
  P <- timed(mgcvST.test(fits[[method]], pairs = pairs, calibration = cal,
    BPPARAM = bp))
  R <- P$value$results
  R$backend <- method
  R$calibration <- cal
  R$seconds <- P$seconds
  rows[[length(rows) + 1L]] <- R
}
pair <- do.call(rbind, rows)
marg <- list()
for (cal in c("liu", "davies")) {
  X <- timed(mgcvST.marginal(I$value, calibration = cal, BPPARAM = bp))
  X$value$backend <- "inla"
  X$value$calibration <- cal
  X$value$seconds <- X$seconds
  marg[[length(marg) + 1L]] <- X$value
}
for (cal in c("liu", "davies")) {
  t0 <- proc.time()[["elapsed"]]
  X <- do.call(rbind, lapply(B$value, function(b)
    as.data.frame(mgcvST:::taps_score_test(b, test.component = 1L,
      method = cal, n_threads = 1L))))
  X$feature_id <- ids
  X$backend <- "bam"
  X$calibration <- cal
  X$seconds <- proc.time()[["elapsed"]] - t0
  marg[[length(marg) + 1L]] <- X
}

dir <- file.path(out, sprintf("task-%03d-%s-r%02d", task, z$case, z$replicate))
dir.create(dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(list(manifest = z, data = dat, mesh = mesh, basis = basis, Y = Y),
  file.path(dir, "input.rds"), compress = FALSE)
saveRDS(I$value, file.path(dir, "inla-fit.rds"), compress = FALSE)
saveRDS(C$value, file.path(dir, "bam-compact-fit.rds"), compress = FALSE)
write.csv(pair, file.path(dir, "pair-tests.csv"), row.names = FALSE)
write.csv(data.table::rbindlist(marg, fill = TRUE),
  file.path(dir, "marginal-tests.csv"), row.names = FALSE)
write.csv(data.frame(backend = c("inla", "bam"), setup_seconds = c(S$seconds, G0$seconds),
  fit_seconds = c(I$seconds, B$seconds), compact_seconds = c(0, C$seconds)),
  file.path(dir, "timings.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(dir, "session-info.txt"))
