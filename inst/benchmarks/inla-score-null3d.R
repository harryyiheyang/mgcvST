# Standalone 3D pair-score experiment; not a public 3D package API.
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
library(Matrix)
library(INLA)
library(mgcvST)
options(warn = 2)
a <- commandArgs(TRUE)
stopifnot(length(a) == 4L)
kind <- a[1L]
mu0 <- as.numeric(a[2L])
rr <- seq.int(as.integer(a[3L]), as.integer(a[4L]))
out <- file.path("artifacts/inla-stress-calibration/null3d-paired", paste(kind, mu0, sep = "-"))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
set.seed(2026091403)
xyz <- as.matrix(expand.grid(x = seq(0, 1, length.out = 4L),
  y = seq(0, 1, length.out = 4L), z = seq(0, 1, length.out = 4L)))
tv <- geometry::delaunayn(xyz)
mesh <- fmesher::fm_mesh_3d(loc = xyz, tv = tv)
n <- 200L
loc <- matrix(runif(n * 3L, 0.01, 0.99), ncol = 3L)
A <- inla.spde.make.A(mesh, loc = loc)
spde <- inla.spde2.matern(mesh, alpha = 2, constr = TRUE,
  B.tau = matrix(c(-0.5 * log(8 * pi * 5), -1), 1, 2),
  B.kappa = matrix(c(log(5), 0), 1, 2))
Q <- inla.spde.precision(spde, theta = 0)
g <- as.numeric(spde$f$extraconstr$A)
Z <- qr.Q(qr(matrix(g, ncol = 1L)), complete = TRUE)[, -1L]
F <- as.matrix(A %*% Z) %*% backsolve(chol(crossprod(Z, as.matrix(Q %*% Z))), diag(ncol(Z)))
scale <- mean(rowSums(F^2)) / 0.45^2
Q <- Q * scale
F <- F / sqrt(scale)
off <- runif(n, -0.2, 0.2)
flat <- list(prior = "flat", param = numeric(), initial = 0)
p <- if (kind == "joint") 40L else 2L
q <- ncol(A)
if (kind == "joint") {
  AA <- kronecker(Diagonal(p), A)
  QQ <- kronecker(Diagonal(p), Q)
  C <- list(A = as.matrix(kronecker(Diagonal(p), matrix(g, nrow = 1))), e = rep(0, p))
  X <- kronecker(Diagonal(p), Matrix(1, n, 1, sparse = TRUE))
  fx <- as.data.frame(diag(p))
  names(fx) <- paste0("b", seq_len(p))
  form <- as.formula(paste("y ~ -1 +", paste(names(fx), collapse = " + "),
    "+ f(field, model='generic0', Cmatrix=QQ, constr=FALSE, rankdef=p, extraconstr=C, hyper=list(prec=flat))"))
} else {
  C <- list(A = matrix(g, nrow = 1L), e = 0)
  form <- y ~ -1 + b0 + f(field, model = "generic0", Cmatrix = Q,
    constr = FALSE, rankdef = 1L, extraconstr = C, hyper = list(prec = flat))
}
for (r in rr) {
  file <- file.path(out, sprintf("rep-%04d.rds", r))
  if (file.exists(file)) stop("Replicate already exists: ", file)
  seed <- 202650000L + r + if (mu0 < 1) 0L else 10000L
  beta <- log(mu0) - log(mean(exp(off + 0.5 * rowSums(F^2))))
  Y <- matrix(NA_integer_, n, p)
  for (j in seq_len(p)) {
    set.seed(seed + 100000L * j)
    eta <- as.numeric(F %*% rnorm(ncol(F)))
    Y[, j] <- rnbinom(n, mu = exp(beta + off + eta), size = 2)
  }
  states <- list()
  meta <- list()
  t0 <- proc.time()[["elapsed"]]
  for (h in seq_len(if (kind == "joint") 1L else 2L)) {
    if (kind == "joint") {
      st <- inla.stack(data = list(y = as.numeric(Y)), A = list(X, AA),
        effects = list(fx, list(field = seq_len(q * p))), compress = TRUE, remove.unused = FALSE)
      exposure <- rep(exp(off), p)
    } else {
      st <- inla.stack(data = list(y = Y[, h]), A = list(1, A),
        effects = list(data.frame(b0 = rep(1, n)), list(field = seq_len(q))))
      exposure <- exp(off)
    }
    fit <- inla(form, family = "nbinomial", data = inla.stack.data(st), E = exposure,
      control.family = list(hyper = list(theta = flat)),
      control.fixed = list(mean = 0, prec = 0, mean.intercept = 0, prec.intercept = 0),
      control.predictor = list(A = inla.stack.A(st), compute = FALSE),
      control.inla = list(strategy = "gaussian", int.strategy = "eb", control.vb = list(enable = FALSE)),
      control.compute = list(config = FALSE, return.marginals = FALSE),
      num.threads = "1:1", safe = FALSE, verbose = FALSE)
    th <- fit$mode$theta
    size <- exp(th[1L])
    tau <- exp(th[2L])
    u <- fit$summary.random$field$mean
    cr <- max(abs(C$A %*% u))
    meta[[h]] <- data.frame(fit = h, mode_status = fit$mode$mode.status,
      log_size = th[1L], log_precision = th[2L], constraint_error = cr,
      warnings = length(fit$misc$warnings))
    stopifnot(all(is.finite(th)), all(is.finite(u)), cr < 1e-6)
    # Conditional Gaussian fixed-effect variance; never replace with expected Fisher Vp.
    for (j in if (kind == "joint") 1:2 else h) {
      uj <- if (kind == "joint") u[(j - 1L) * q + seq_len(q)] else u
      tag <- if (kind == "joint") paste0("b", j) else "b0"
      b <- fit$summary.fixed[tag, "mean"]
      Vp <- fit$summary.fixed[tag, "sd"]^2
      eta1 <- b + as.numeric(A %*% uj) + off
      mu <- exp(eta1)
      D <- 1 / mu + exp(-th[1L])
      e <- eta1 + (Y[, j] - mu) / mu - off
      FF <- F * exp(-th[2L] / 2)
      V <- diag(D) + tcrossprod(FF)
      Vi <- chol2inv(chol(V))
      vx <- rowSums(Vi)
      P <- Vi - tcrossprod(vx) * Vp
      states[[j]] <- list(a = as.numeric(crossprod(FF, P %*% e)),
        M = crossprod(FF, P %*% FF), Vp = Vp, Vp_expected = 1 / sum(vx),
        G = crossprod(FF, Vi %*% FF), w = as.numeric(crossprod(FF, vx)),
        F = FF, D = D, mu = mu, e = e)
    }
    rm(fit, st)
    gc()
  }
  rows <- list()
  U <- sum(states[[1L]]$a * states[[2L]]$a)
  mineig <- min(vapply(states, function(x) min(eigen(x$M, symmetric = TRUE, only.values = TRUE)$values), numeric(1L)))
  valid <- is.finite(U) && mineig >= -1e-10
  for (cal in c("davies", "liu")) {
    z <- if (valid) rkhs_score_calibrate(U, states[[1L]]$M, states[[2L]]$M, method = cal) else
      list(p_two_sided = NA_real_, information = NA_real_)
    rows[[cal]] <- data.frame(case = paste(kind, mu0, sep = "-"), replicate = r,
      seed = seed, calibration = cal, p_value = z$p_two_sided, statistic = U,
      information = z$information, minimum_M_eigenvalue = mineig,
      error = if (valid) NA_character_ else "Non-PSD plug-in score matrix; retained as invalid")
  }
  z <- do.call(rbind, rows)
  saveRDS(list(rows = z, hyper = do.call(rbind, meta), seconds = proc.time()[["elapsed"]] - t0,
    Y = Y, states = states, geometry = list(nodes = q, observations = n, genes = p,
      kappa = 5, constraint = "volume integral", sigma2 = 0.45^2)), file)
  write.csv(z, sub("rds$", "csv", file), row.names = FALSE)
}
