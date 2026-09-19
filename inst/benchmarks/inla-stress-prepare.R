library(Matrix)
library(INLA)
library(mgcvST)
options(warn = 2)
out <- "artifacts/inla-stress-calibration"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
p <- "artifacts/inla3d-transfer/spde3d_transfer_2026-09-13/data/magic"
spec <- jsonlite::read_json(file.path(p, "mesh_contract.json"), simplifyVector = TRUE)
xyz <- as.matrix(read.csv(file.path(p, "nodes.csv")))
tv <- as.matrix(read.csv(file.path(p, "tetrahedra.csv")))
mesh <- fmesher::fm_mesh_3d(loc = xyz, tv = tv)
k <- spec$kappa_fixed
spde <- inla.spde2.matern(mesh, alpha = 2, constr = TRUE,
  B.tau = matrix(c(-0.5 * log(8 * pi * k), -1), 1, 2),
  B.kappa = matrix(c(log(k), 0), 1, 2))
Q <- inla.spde.precision(spde, theta = 0)
g <- as.numeric(spde$f$extraconstr$A)
d <- read.delim(gzfile(file.path(p, "aligned_points.tsv.gz")))
loc <- sweep(sweep(as.matrix(d[, c("x_aligned", "y_aligned", "z")]),
  2, spec$source_origin, "-"), 2, spec$source_units_per_mesh_unit, "/")
A <- inla.spde.make.A(mesh, loc = loc)
E <- d$total_umi / 10000
L <- expand(Cholesky(Q, LDL = FALSE, perm = TRUE))
Qg <- as.numeric(solve(Q, g))
set.seed(2026091401)
u <- as.matrix(t(L$P) %*% solve(t(L$L), matrix(rnorm(mesh$n * 40L), mesh$n, 40L)))
u <- 0.4 * (u - tcrossprod(Qg, as.numeric(crossprod(g, u)) / sum(g * Qg)))
stopifnot(max(abs(crossprod(g, u))) < 1e-8)
Y <- matrix(NA_integer_, nrow(A), 40L)
for (j in 1:40) {
  set.seed(2026091401L + j)
  Y[, j] <- rnbinom(nrow(A), mu = E * exp(0.75 + as.numeric(A %*% u[, j])), size = 15)
}
saveRDS(list(A = A, Q = Q, g = g, Y = Y, E = E, u = u,
  kappa = k, sigma = 0.4, size = 15, seed = 2026091401L),
  file.path(out, "stress-input.rds"), compress = FALSE)
set.seed(2026091402)
n <- 200L
loc <- matrix(runif(n * 2L, 0.01, 0.99), ncol = 2L)
v <- as.matrix(expand.grid(x = seq(0, 1, length.out = 6L), y = seq(0, 1, length.out = 6L)))
m <- list(loc = v, graph = list(tv = geometry::delaunayn(v)))
d <- data.frame(x = loc[, 1], y = loc[, 2], off = runif(n, -0.2, 0.2))
B <- spde_basis(m, loc, kappa = 5, project_intercept = TRUE)
B$component <- B$score.component <- "global"
F <- B$B %*% backsolve(chol(B$Q), diag(ncol(B$Q)))
F <- 0.45 * F / sqrt(mean(rowSums(F^2)))
S <- inlaST.set(response ~ offset(off), d, B, family = mgcv::nb())
saveRDS(list(S = S, F = F, d = d, B = B), file.path(out, "null2d-input.rds"))
capture.output(sessionInfo(), file = file.path(out, "session.txt"))
