# Run from the package checkout after inla3d-mesh.py.
library(Matrix)
library(INLA)
library(fmesher)
options(warn = 2)

base <- "artifacts/inla3d"
out <- file.path(base, "eb")
dir.create(file.path(out, "fits"), recursive = TRUE, showWarnings = FALSE)
args <- commandArgs(TRUE)
nr <- if (length(args)) as.integer(args[1L]) else 10L
if (length(nr) != 1L || is.na(nr) || nr < 1L) stop("Replicates must be positive.")
manifest <- read.csv(file.path(base, "mesh", "manifest.csv"))
set.seed(20260913)
n <- 30000L
loc <- cbind(runif(n, 0.05, 2.95), runif(n, 0.05, 1.95),
             rep(seq(0.05, 0.95, length.out = 10L), each = 3000L))
offset <- rnorm(n, 0, 0.15)
r2 <- rowSums(sweep(loc, 2L, c(1, 1, 0.5))^2)
roi <- r2 < 0.35^2
truth <- cbind(broad = 0.45 * sin(pi * loc[, 1L] / 1.5) *
                 cos(pi * loc[, 2L] / 2) + 0.25 * sin(2 * pi * loc[, 3L]),
               focal = 0.3 * sin(pi * loc[, 1L] / 1.5) +
                 1.2 * exp(-r2 / (2 * 0.12^2)))
truth <- sweep(truth, 2L, colMeans(truth))
Y <- array(NA_integer_, c(n, 2L, nr))
for (s in seq_len(2L)) {
  for (r in seq_len(nr)) {
    set.seed(20260913L + 1000L * s + r)
    Y[, s, r] <- rnbinom(n, mu = exp(log(6) + offset + truth[, s]), size = 8)
  }
}
saveRDS(list(loc = loc, offset = offset, truth = truth, Y = Y, roi = roi,
             seed = 20260913L, size = 8), file.path(out, "input.rds"))
writeLines(c(capture.output(sessionInfo()),
             paste("source", system("git rev-parse HEAD", intern = TRUE))),
           file.path(out, "session.txt"))
flat <- list(prior = "flat", param = numeric(), initial = 0, fixed = FALSE)
kappa <- 5
geometry <- vector("list", nrow(manifest))
G <- list()
for (m in seq_len(nrow(manifest))) {
  label <- manifest$mesh[m]
  writeLines(paste0("geometry-", label), file.path(out, "current-task.txt"))
  xyz <- as.matrix(read.csv(file.path(base, "mesh", paste0(label, "-vertices.csv"))))
  tv <- as.matrix(read.csv(file.path(base, "mesh", paste0(label, "-tetrahedra.csv"))))
  mesh <- fm_mesh_3d(loc = xyz, tv = tv)
  if (ncol(tv) != 4L || mesh$n >= 3000L) stop("Expected a tetrahedral mesh below 3000 nodes.")
  t0 <- proc.time()[["elapsed"]]
  E <- fm_fem(mesh, order = 2L)
  fem_seconds <- proc.time()[["elapsed"]] - t0
  if (any(E$ta <= 0) || any(E$va <= 0)) stop("Non-positive tetrahedron or vertex volume.")
  if (abs(sum(E$ta) - 6) > 1e-8) stop("Mesh does not cover the synthetic volume.")
  if (max(abs(as.numeric(E$g1 %*% rep(1, mesh$n)))) > 1e-8) stop("Stiffness constant test failed.")
  t0 <- proc.time()[["elapsed"]]
  A <- fm_basis(mesh, loc = loc)
  A_seconds <- proc.time()[["elapsed"]] - t0
  row_error <- max(abs(rowSums(A) - 1))
  affine_error <- max(abs(as.matrix(A %*% xyz) - loc))
  if (any(!is.finite(A@x)) || min(A@x) < -1e-10 || row_error > 1e-9 ||
      affine_error > 1e-9 || any(rowSums(A != 0) > 4L)) stop("3D interpolation check failed.")
  t0 <- proc.time()[["elapsed"]]
  Q <- forceSymmetric(kappa^4 * E$c0 + 2 * kappa^2 * E$g1 + E$g2)
  Q_seconds <- proc.time()[["elapsed"]] - t0
  # alpha=2 in R3: nu=1/2. Q is scaled only by the fitted precision.
  t0 <- proc.time()[["elapsed"]]
  L <- Cholesky(Q, LDL = FALSE, perm = TRUE)
  chol_seconds <- proc.time()[["elapsed"]] - t0
  g <- as.numeric(colMeans(A))
  geometry[[m]] <- list(A = A, Q = Q, g = g)
  G[[m]] <- data.frame(mesh = label, n = n, nodes = mesh$n,
    tetrahedra = nrow(tv), fem_seconds = fem_seconds, A_seconds = A_seconds,
    Q_seconds = Q_seconds, chol_seconds = chol_seconds, nnz_A = nnzero(A),
    nnz_Q = nnzero(Q), nnz_L = nnzero(expand(L)$L),
    A_MiB = as.numeric(object.size(A)) / 2^20,
    Q_MiB = as.numeric(object.size(Q)) / 2^20,
    row_error = row_error, affine_error = affine_error, volume = sum(E$ta))
}
write.csv(do.call(rbind, G), file.path(out, "geometry.csv"), row.names = FALSE)
saveRDS(geometry, file.path(out, "geometry.rds"))
R1 <- list()
i <- 0L
for (s in seq_len(2L)) {
  for (r in seq_len(nr)) {
    for (m in seq_len(nrow(manifest))) {
      label <- manifest$mesh[m]
      task <- paste(colnames(truth)[s], r, label, sep = "-")
      writeLines(task, file.path(out, "current-task.txt"))
      t0 <- proc.time()[["elapsed"]]
      A <- geometry[[m]]$A
      Q <- geometry[[m]]$Q
      C <- list(A = matrix(geometry[[m]]$g, nrow = 1L), e = 0)
      st <- inla.stack(data = list(y = Y[, s, r]), A = list(1, A),
        effects = list(data.frame(b0 = rep(1, n), off = offset),
                       field = seq_len(ncol(A))))
      stack_seconds <- proc.time()[["elapsed"]] - t0
      f <- y ~ -1 + b0 + offset(off) + f(field, model = "generic0", Cmatrix = Q,
        constr = FALSE, rankdef = 1L, extraconstr = C, hyper = list(prec = flat))
      t0 <- proc.time()[["elapsed"]]
      fit <- inla(f, family = "nbinomial", data = inla.stack.data(st),
        control.predictor = list(A = inla.stack.A(st), compute = FALSE),
        control.family = list(hyper = list(theta = flat)),
        control.fixed = list(mean = 0, prec = 0, mean.intercept = 0, prec.intercept = 0),
        control.inla = list(strategy = "gaussian", int.strategy = "eb",
                            control.vb = list(enable = FALSE)),
        control.compute = list(config = FALSE, return.marginals = FALSE),
        num.threads = "1:1", verbose = FALSE)
      fit_seconds <- proc.time()[["elapsed"]] - t0
      u <- fit$summary.random$field$mean
      est <- as.numeric(A %*% u)
      e <- est - truth[, s]
      if (any(!is.finite(est))) stop("Non-finite fitted field.")
      cr <- sum(geometry[[m]]$g * u)
      if (abs(cr) > 1e-6) stop("Observation-mean-zero constraint failed.")
      i <- i + 1L
      R1[[i]] <- data.frame(scenario = colnames(truth)[s], replicate = r, mesh = label,
        nodes = ncol(A), stack_seconds = stack_seconds, fit_seconds = fit_seconds,
        mode_status = fit$mode$mode.status, correlation = cor(est, truth[, s]),
        rmse = sqrt(mean(e^2)), roi_rmse = sqrt(mean(e[roi]^2)),
        outside_rmse = sqrt(mean(e[!roi]^2)), constraint_error = cr,
        beta = fit$summary.fixed["b0", "mean"],
        nb_size_mode = unname(exp(fit$mode$theta[1L])),
        precision_mode = unname(exp(fit$mode$theta[2L])))
      saveRDS(list(metrics = R1[[i]], field = est,
                   fixed = fit$summary.fixed, hyper = fit$summary.hyperpar,
                   mode = fit$mode, cpu = fit$cpu.used, flat = flat,
                   integration = "eb", latent = "gaussian", kappa = kappa),
              file.path(out, "fits", paste0(task, ".rds")))
      write.csv(do.call(rbind, R1), file.path(out, "fits.csv"), row.names = FALSE)
      print(R1[[i]])
      rm(fit, st)
      gc()
    }
  }
}
writeLines("complete", file.path(out, "current-task.txt"))
