# One isolated stress job. Arguments: label, independent/joint, n, genes, threads, replicate.
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
library(Matrix)
library(INLA)
options(warn = 2)
a <- commandArgs(TRUE)
stopifnot(length(a) == 6L)
label <- a[1L]
kind <- a[2L]
n <- as.integer(a[3L])
p <- as.integer(a[4L])
nt <- as.integer(a[5L])
r <- as.integer(a[6L])
out <- file.path("artifacts/inla-stress-calibration/stress", label)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
d <- readRDS("artifacts/inla-stress-calibration/stress-input.rds")
j <- unique(round(seq(1, nrow(d$A), length.out = n)))
A <- d$A[j, ]
Q <- d$Q
g <- d$g
E <- d$E[j]
q <- ncol(A)
flat <- list(prior = "flat", param = numeric(), initial = 0)
ids <- if (kind == "independent") r else seq_len(p)
stopifnot(all(ids <= ncol(d$Y)))
if (kind == "joint") {
  AA <- kronecker(Diagonal(p), A)
  QQ <- kronecker(Diagonal(p), Q)
  C <- list(A = as.matrix(kronecker(Diagonal(p), matrix(g, nrow = 1L))), e = rep(0, p))
  X <- kronecker(Diagonal(p), Matrix(1, nrow(A), 1, sparse = TRUE))
  fx <- as.data.frame(diag(p))
  names(fx) <- paste0("b", seq_len(p))
  y <- as.numeric(d$Y[j, ids])
  st <- inla.stack(data = list(y = y), A = list(X, AA),
    effects = list(fx, list(field = seq_len(q * p))), compress = TRUE, remove.unused = FALSE)
  form <- as.formula(paste("y ~ -1 +", paste(names(fx), collapse = " + "),
    "+ f(field, model='generic0', Cmatrix=QQ, constr=FALSE, rankdef=p, extraconstr=C, hyper=list(prec=flat))"))
  exposure <- rep(E, p)
} else {
  C <- list(A = matrix(g, 1), e = 0)
  y <- d$Y[j, ids]
  st <- inla.stack(data = list(y = y), A = list(1, A),
    effects = list(data.frame(b0 = rep(1, nrow(A))), list(field = seq_len(q))))
  form <- y ~ -1 + b0 + f(field, model = "generic0", Cmatrix = Q,
    constr = FALSE, rankdef = 1L, extraconstr = C, hyper = list(prec = flat))
  exposure <- E
}
rm(d)
gc()
t0 <- proc.time()[["elapsed"]]
fit <- inla(form, family = "nbinomial", data = inla.stack.data(st), E = exposure,
  control.family = list(hyper = list(theta = flat)),
  control.fixed = list(mean = 0, prec = 0, mean.intercept = 0, prec.intercept = 0),
  control.predictor = list(A = inla.stack.A(st), compute = FALSE),
  control.inla = list(strategy = "gaussian", int.strategy = "eb", control.vb = list(enable = FALSE)),
  control.compute = list(config = FALSE, return.marginals = FALSE),
  num.threads = paste0(nt, ":1"), safe = FALSE, verbose = FALSE)
elapsed <- proc.time()[["elapsed"]] - t0
u <- fit$summary.random$field$mean
th <- fit$mode$theta
cr <- max(abs(C$A %*% u))
stopifnot(length(th) == 2L, all(is.finite(th)), all(is.finite(u)), cr < 1e-6)
z <- data.frame(label = label, kind = kind, n = n, genes = length(ids),
  nodes_per_gene = q, threads = nt, replicate = r, seconds = elapsed,
  mode_status = fit$mode$mode.status, nb_size = exp(th[1]),
  sigma = exp(-th[2] / 2), constraint_error = cr,
  warnings = length(fit$misc$warnings))
saveRDS(list(metrics = z, u = u, fixed = fit$summary.fixed, theta = th),
  file.path(out, sprintf("job-%03d.rds", r)))
write.csv(z, file.path(out, sprintf("job-%03d.csv", r)), row.names = FALSE)
