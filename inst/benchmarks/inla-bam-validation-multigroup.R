#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
Sys.setenv(LC_ALL = "C", OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1", RCPP_PARALLEL_NUM_THREADS = "1")
args <- commandArgs(trailingOnly = TRUE)
replicate <- if (length(args)) as.integer(sub("^--replicate=", "", args[1L])) else NA_integer_
if (length(replicate) != 1L || is.na(replicate) || replicate < 1L || replicate > 10L) {
  stop("Supply --replicate=1,...,10.")
}
out <- Sys.getenv("MGCVST_MULTIGROUP_OUTPUT",
  "artifacts/inla-bam-validation/multigroup")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
lib <- normalizePath(Sys.getenv("MGCVST_VALIDATION_LIBRARY",
  "artifacts/inla-bam-validation/validated-library"), mustWork = TRUE)
.libPaths(c(lib, .libPaths()))
library(mgcvST)
library(mgcv)
flat <- list(prior = "flat", param = numeric(), initial = 0)
ctl <- list(precision_prior = flat, gaussian_precision_prior = flat,
  nb_size_prior = flat)
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
  rownames(sp) <- ids
  components <- names(H$target)
  component_lambda <- do.call(rbind, lapply(states, function(x) {
    vapply(components, function(component) {
      j <- x$H$target[[component]]
      x$H$sp[x$H$smooth[[j]]$sp_index]
    }, numeric(1L))
  }))
  dimnames(component_lambda) <- list(ids, components)
  structure(list(feature_id = ids,
    working_error = do.call(cbind, lapply(states, function(x) x$W$working_error)),
    working_variance = do.call(cbind, lapply(states, function(x) x$W$working_variance)),
    dispersion = setNames(vapply(states, function(x) x$W$dispersion, numeric(1L)), ids),
    lambda = component_lambda[, "global"], component_lambda = component_lambda,
    smoothing_parameters = sp,
    family_parameters = setNames(lapply(states, function(x) x$W$family_parameters), ids),
    nuisance_covariance = setNames(lapply(states, function(x) x$N$covariance), ids),
    geometry = H, row_id = H$row_id, score_components = components,
    model_setting = "global_local", test_engine = NULL,
    diagnostics = data.frame(index = seq_along(ids), feature_id = ids,
      converged = vapply(fits, function(x) isTRUE(x$converged), logical(1L)),
      error_message = NA_character_), timing = list(backend = "bam", workers = 1L)),
    class = c("mgcvST_model_fit", "mgcvST_fit", "mgcvST"))
}

seed <- 928000L + replicate
set.seed(seed)
n <- 200L
vertices <- as.matrix(expand.grid(x = seq(0, 1, length.out = 6L),
  y = seq(0, 1, length.out = 6L)))
mesh <- list(loc = vertices, graph = list(tv = geometry::delaunayn(vertices)))
dat <- data.frame(x = runif(n, .01, .99), y = runif(n, .01, .99),
  u = runif(n, .01, .99), v = runif(n, .01, .99), z = rnorm(n),
  offset0 = runif(n, -.15, .15))
global <- spde_basis(mesh, as.matrix(dat[c("x", "y")]), kappa = 4,
  project_intercept = TRUE)
local <- spde_basis(mesh, as.matrix(dat[c("u", "v")]), kappa = 9,
  project_intercept = TRUE)
global$component <- global$score.component <- "global"
local$component <- local$score.component <- "local"
Fg <- global$B %*% backsolve(chol(global$Q), diag(ncol(global$Q)))
Fl <- local$B %*% backsolve(chol(local$Q), diag(ncol(local$Q)))
Fg <- Fg / sqrt(mean(rowSums(Fg^2)))
Fl <- Fl / sqrt(mean(rowSums(Fl^2)))
ug <- matrix(rnorm(ncol(Fg) * 3L), ncol = 3L)
ul <- matrix(rnorm(ncol(Fl) * 3L), ncol = 3L)
ug[, 2L] <- .65 * ug[, 1L] + sqrt(1 - .65^2) * ug[, 2L]
ul[, 3L] <- -.6 * ul[, 1L] + sqrt(1 - .6^2) * ul[, 3L]
eta <- .3 * Fg %*% ug + .25 * Fl %*% ul
eta <- sweep(eta, 1L, .3 * dat$z + dat$offset0, "+")
Y <- t(vapply(1:3, function(j) rnbinom(n, mu = exp(log(2) + eta[, j]), size = 4),
  numeric(n)))
ids <- paste0("feature", 1:3)
rownames(Y) <- ids
dat$response <- Y[1L, ]
form <- response ~ z + offset(offset0) +
  s(x, y, bs = "spde", xt = global) +
  s(u, v, bs = "spde", xt = local)
bp <- BiocParallel::SerialParam()

S <- timed(inlaST.set(form, dat, nb(), control = ctl))
I <- timed(inlaST.estimate(Y, S$value, feature_id = ids, retain_smooth = TRUE,
  retain_marginal = TRUE, BPPARAM = bp, control = ctl))
G0 <- timed(bam(form, data = dat, family = nb(), method = "fREML",
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

wrows <- list()
for (backend in c("inla", "bam")) for (group in list("global", "local", c("global", "local"))) {
  fit <- if (backend == "inla") I$value else C$value
  W <- timed(mgcvST.wgcna(fit, ids, group = group,
    wgcna.para = list(minClusterSize = 2L, deepSplit = 0L)))
  X <- W$value$modules
  X$backend <- backend
  X$group <- paste(group, collapse = "+")
  X$seconds <- W$seconds
  wrows[[length(wrows) + 1L]] <- X
}

dout <- file.path(out, sprintf("replicate-%02d", replicate))
dir.create(dout, recursive = TRUE, showWarnings = FALSE)
saveRDS(list(seed = seed, data = dat, mesh = mesh, global = global,
  local = local, Y = Y), file.path(dout, "input.rds"), compress = FALSE)
saveRDS(I$value, file.path(dout, "inla-fit.rds"), compress = FALSE)
saveRDS(C$value, file.path(dout, "bam-compact-fit.rds"), compress = FALSE)
write.csv(do.call(rbind, wrows), file.path(dout, "wgcna.csv"), row.names = FALSE)
write.csv(data.frame(backend = c("inla", "bam"), setup_seconds = c(S$seconds, G0$seconds),
  fit_seconds = c(I$seconds, B$seconds), compact_seconds = c(0, C$seconds),
  test_engine_available = c(!is.null(I$value$test_engine), !is.null(C$value$test_engine)),
  marginal_components_available = FALSE), file.path(dout, "status.csv"), row.names = FALSE)
writeLines(c("Current model-score dispatch implements only test_engine='single_model'.",
  "The global_local fits correctly record test_engine=NULL; pair all/omit/part is unavailable.",
  "The marginal API retains one feature-level result and does not expose separate global/local tests."),
  file.path(dout, "api-boundary.txt"))
writeLines(capture.output(sessionInfo()), file.path(dout, "session-info.txt"))
