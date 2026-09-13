# Run from the canonical checkout; argument is the extracted transfer directory.
library(Matrix)
library(INLA)
library(fmesher)
options(warn = 2)
args <- commandArgs(TRUE)
if (length(args) != 1L) stop("Supply the extracted spde3d transfer directory.")
src <- normalizePath(args[1L], winslash = "/", mustWork = TRUE)
out <- "artifacts/inla3d-transfer/generic-flat"
dir.create(file.path(out, "fits"), recursive = TRUE, showWarnings = FALSE)
if (file.exists(file.path(out, "started.txt"))) stop("This experiment has already started; inspect its saved outputs.")
p <- file.path(src, "data", "magic")
spec <- jsonlite::read_json(file.path(p, "mesh_contract.json"), simplifyVector = TRUE)
xyz <- as.matrix(read.csv(file.path(p, "nodes.csv")))
tv <- as.matrix(read.csv(file.path(p, "tetrahedra.csv")))
mesh <- fm_mesh_3d(loc = xyz, tv = tv)
k <- spec$kappa_fixed
spde <- inla.spde2.matern(mesh, alpha = 2, constr = TRUE,
  B.tau = matrix(c(-0.5 * log(8 * pi * k), -1), 1, 2),
  B.kappa = matrix(c(log(k), 0), 1, 2),
  theta.prior.mean = 0, theta.prior.prec = 1)
flat <- list(prior = "flat", param = numeric(), initial = 0, fixed = FALSE)
stopifnot(spde$n.theta == 1L)
d <- read.delim(gzfile(file.path(p, "aligned_points.tsv.gz")), check.names = FALSE)
g <- read.delim(gzfile(file.path(p, "Snap25.tsv.gz")), check.names = FALSE)
ix <- match(d$point_id, g$point_id)
stopifnot(!anyNA(ix), !anyDuplicated(d$point_id), !anyDuplicated(g$point_id))
d$count <- g$count[ix]
stopifnot(all(d$total_umi == g$total_umi[ix]), all(d$total_umi > 0),
  all(d$count >= 0), all(d$count == round(d$count)))
loc <- sweep(sweep(as.matrix(d[, c("x_aligned", "y_aligned", "z")]),
  2, spec$source_origin, "-"), 2, spec$source_units_per_mesh_unit, "/")
A <- inla.spde.make.A(mesh, loc = loc)
err <- max(abs(as.matrix(A %*% xyz) - loc))
stopifnot(max(abs(rowSums(A) - 1)) < 1e-7, err < 1e-7,
  max(rowSums(A != 0)) <= 4L, min(A@x) > -1e-10)
Q <- inla.spde.precision(spde, theta = 0)
C <- spde$f$extraconstr
F <- fm_fem(mesh, order = 2L)
Q1 <- forceSymmetric((k^4 * F$c0 + 2 * k^2 * F$g1 + F$g2) / (8 * pi * k))
qerr <- max(abs(Q - Q1)) / max(abs(Q))
stopifnot(qerr < 1e-10)
E <- d$total_umi / 10000
sel <- unique(round(seq(1, nrow(d), length.out = 5000L)))
write.csv(data.frame(n = nrow(d), nodes = mesh$n, tetrahedra = nrow(tv),
  nnz_A = nnzero(A), A_MiB = as.numeric(object.size(A)) / 2^20,
  coordinate_error = err, native_generic_Q_relative_error = qerr,
  kappa = k, practical_range = 2 / k), file.path(out, "geometry.csv"), row.names = FALSE)
writeLines(c(capture.output(sessionInfo()), paste("source", system("git rev-parse HEAD", intern = TRUE))),
  file.path(out, "session.txt"))
saveRDS(list(spatial = flat, nb_size = flat, constraint = spde$f$extraconstr,
  kappa = k, source = src, integration = "eb", latent = "gaussian"), file.path(out, "config.rds"))
writeLines("started", file.path(out, "started.txt"))
R1 <- list()
# Two real-data sizes, then ten NB response replicates on the full real design.
for (r in 0:11) {
  j <- if (r == 0L) sel else seq_len(nrow(d))
  n <- length(j)
  if (r < 2L) {
    y <- d$count[j]
  } else {
    set.seed(20260913L + r)
    y <- rnbinom(n, mu = mu0, size = size0)
  }
  st <- inla.stack(data = list(y = y), A = list(1, A[j, ]),
    effects = list(data.frame(b0 = rep(1, n)), list(spatial = seq_len(mesh$n))),
    compress = TRUE, remove.unused = FALSE)
  # Q(sigma) = Q(1) / sigma^2; flat log precision is flat log sigma.
  f <- y ~ -1 + b0 + f(spatial, model = "generic0", Cmatrix = Q,
    constr = FALSE, extraconstr = C, rankdef = 1L, diagonal = 0,
    hyper = list(prec = flat))
  writeLines(paste("fit", r), file.path(out, "current-task.txt"))
  t0 <- proc.time()[["elapsed"]]
  fit <- inla(f, family = "nbinomial", data = inla.stack.data(st), E = E[j],
    control.family = list(variant = 0, hyper = list(theta = flat)),
    control.fixed = list(mean = 0, prec = 0, mean.intercept = 0, prec.intercept = 0),
    control.predictor = list(A = inla.stack.A(st), compute = FALSE),
    control.compute = list(config = FALSE, dic = FALSE, waic = FALSE, cpo = FALSE,
      mlik = FALSE, return.marginals = FALSE, return.marginals.predictor = FALSE),
    control.inla = list(strategy = "gaussian", int.strategy = "eb",
      control.vb = list(enable = FALSE)), inla.mode = "experimental",
    num.threads = "2:1", verbose = FALSE, safe = FALSE, keep = FALSE)
  elapsed <- proc.time()[["elapsed"]] - t0
  u <- fit$summary.random$spatial$mean
  b <- fit$summary.fixed["b0", "mean"]
  eta <- as.numeric(b + A[j, ] %*% u)
  mu <- E[j] * exp(eta)
  cr <- as.numeric(spde$f$extraconstr$A %*% u)
  th <- fit$mode$theta
  stopifnot(length(th) == 2L, all(is.finite(th)), all(is.finite(mu)),
    all(mu > 0), abs(cr) < 1e-6)
  size <- exp(th[1L])
  sigma <- exp(-th[2L] / 2)
  stopifnot(is.finite(size), size > 0, is.finite(sigma), sigma > 0)
  R1[[r + 1L]] <- data.frame(kind = if (r < 2L) "real" else "simulation",
    replicate = if (r < 2L) 0L else r - 1L, n = n, nodes = mesh$n,
    seconds = elapsed, mode_status = fit$mode$mode.status,
    nb_size = size, sigma = sigma, beta = b, constraint_error = cr,
    warnings = length(fit$misc$warnings), count_rmse = sqrt(mean((y - mu)^2)),
    count_correlation = cor(y, mu),
    eta_rmse = if (r < 2L) NA_real_ else sqrt(mean((eta - eta0)^2)),
    eta_correlation = if (r < 2L) NA_real_ else cor(eta, eta0))
  saveRDS(list(metrics = R1[[r + 1L]], mode = fit$mode, hyper = fit$summary.hyperpar,
    fixed = fit$summary.fixed, u = u, eta = eta, y = y, cpu = fit$cpu.used,
    warnings = fit$misc$warnings, seed = if (r < 2L) NA_integer_ else 20260913L + r),
    file.path(out, "fits", paste0("fit-", r, ".rds")))
  write.csv(do.call(rbind, R1), file.path(out, "fits.csv"), row.names = FALSE)
  if (r == 1L) {
    eta0 <- eta
    mu0 <- mu
    size0 <- size
    saveRDS(list(eta = eta0, mu = mu0, size = size0, E = E,
      point_id = d$point_id, description = "Conditional parametric recovery on MAGIC real geometry"),
      file.path(out, "simulation-truth.rds"))
  }
  rm(fit, st)
  gc()
}
writeLines("complete", file.path(out, "current-task.txt"))
