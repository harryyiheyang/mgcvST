Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
library(Matrix)
library(INLA)
options(warn = 2)

out <- "artifacts/inla-stress-calibration/replicate-pilot"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
d <- readRDS("artifacts/inla-stress-calibration/stress-input.rds")
n <- 500L
p <- 3L
j <- unique(round(seq(1, nrow(d$A), length.out = n)))
A <- d$A[j, ]
Q <- d$Q
g <- d$g
E <- d$E[j]
Y <- d$Y[j, seq_len(p), drop = FALSE]
q <- ncol(A)
flat <- list(prior = "flat", param = numeric(), initial = 0)

AA <- kronecker(Diagonal(p), A)
QQ <- kronecker(Diagonal(p), Q)
C <- list(A = as.matrix(kronecker(Diagonal(p), matrix(g, nrow = 1L))), e = rep(0, p))
X <- kronecker(Diagonal(p), Matrix(1, nrow(A), 1, sparse = TRUE))
fx <- as.data.frame(diag(p))
names(fx) <- paste0("b", seq_len(p))
y <- as.numeric(Y)
exposure <- rep(E, p)

st_block <- inla.stack(data = list(y = y), A = list(X, AA),
  effects = list(fx, list(field = seq_len(q * p))), compress = TRUE,
  remove.unused = FALSE)
form_block <- as.formula(paste("y ~ -1 +", paste(names(fx), collapse = " + "),
  "+ f(field, model='generic0', Cmatrix=QQ, constr=FALSE, rankdef=p, extraconstr=C, hyper=list(prec=flat))"))
t0 <- proc.time()[["elapsed"]]
fit_block <- inla(form_block, family = "nbinomial", data = inla.stack.data(st_block),
  E = exposure, control.family = list(hyper = list(theta = flat)),
  control.fixed = list(mean = 0, prec = 0, mean.intercept = 0, prec.intercept = 0),
  control.predictor = list(A = inla.stack.A(st_block), compute = FALSE),
  control.inla = list(strategy = "gaussian", int.strategy = "eb",
    control.vb = list(enable = FALSE)),
  control.compute = list(config = FALSE, return.marginals = FALSE),
  num.threads = "1:1", safe = FALSE, verbose = FALSE)
seconds_block <- proc.time()[["elapsed"]] - t0
u_block <- matrix(fit_block$summary.random$field$mean, q, p)
b_block <- fit_block$summary.fixed$mean
theta_block <- fit_block$mode$theta
eta_block <- sweep(as.matrix(A %*% u_block), 2, b_block, "+")
ll_block <- sum(dnbinom(Y, mu = E * exp(eta_block), size = exp(theta_block[1L]), log = TRUE))
constraint_block <- as.numeric(crossprod(g, u_block))
block <- list(u = u_block, beta = b_block, theta = theta_block,
  constraint = constraint_block, conditional_loglik = ll_block,
  marginal_loglik = fit_block$mlik, seconds = seconds_block,
  mode_status = fit_block$mode$mode.status, warnings = fit_block$misc$warnings)
saveRDS(block, file.path(out, "block.rds"))
rm(fit_block, st_block)
gc()

field <- rep(seq_len(q), p)
rep_id <- rep(seq_len(p), each = q)
C1 <- list(A = matrix(g, nrow = 1L), e = 0)
st_replicate <- inla.stack(data = list(y = y), A = list(X, AA),
  effects = list(fx, data.frame(field = field, rep_id = rep_id)),
  compress = TRUE, remove.unused = FALSE)
form_replicate <- as.formula(paste("y ~ -1 +", paste(names(fx), collapse = " + "),
  "+ f(field, model='generic0', Cmatrix=Q, replicate=rep_id, nrep=p, constr=FALSE, rankdef=1L, extraconstr=C1, hyper=list(prec=flat))"))
t0 <- proc.time()[["elapsed"]]
fit_replicate <- inla(form_replicate, family = "nbinomial",
  data = inla.stack.data(st_replicate), E = exposure,
  control.family = list(hyper = list(theta = flat)),
  control.fixed = list(mean = 0, prec = 0, mean.intercept = 0, prec.intercept = 0),
  control.predictor = list(A = inla.stack.A(st_replicate), compute = FALSE),
  control.inla = list(strategy = "gaussian", int.strategy = "eb",
    control.vb = list(enable = FALSE)),
  control.compute = list(config = FALSE, return.marginals = FALSE),
  num.threads = "1:1", safe = FALSE, verbose = FALSE)
seconds_replicate <- proc.time()[["elapsed"]] - t0
u_replicate <- matrix(fit_replicate$summary.random$field$mean, q, p)
b_replicate <- fit_replicate$summary.fixed$mean
theta_replicate <- fit_replicate$mode$theta
eta_replicate <- sweep(as.matrix(A %*% u_replicate), 2, b_replicate, "+")
ll_replicate <- sum(dnbinom(Y, mu = E * exp(eta_replicate),
  size = exp(theta_replicate[1L]), log = TRUE))
constraint_replicate <- as.numeric(crossprod(g, u_replicate))
replicate_fit <- list(u = u_replicate, beta = b_replicate,
  theta = theta_replicate, constraint = constraint_replicate,
  conditional_loglik = ll_replicate, marginal_loglik = fit_replicate$mlik,
  seconds = seconds_replicate, mode_status = fit_replicate$mode$mode.status,
  warnings = fit_replicate$misc$warnings)
saveRDS(replicate_fit, file.path(out, "replicate.rds"))

z <- data.frame(n = n, genes = p, nodes_per_gene = q,
  block_seconds = seconds_block, replicate_seconds = seconds_replicate,
  max_abs_u = max(abs(u_block - u_replicate)),
  max_abs_beta = max(abs(b_block - b_replicate)),
  max_abs_theta = max(abs(theta_block - theta_replicate)),
  max_abs_constraint_block = max(abs(constraint_block)),
  max_abs_constraint_replicate = max(abs(constraint_replicate)),
  conditional_loglik_block = ll_block,
  conditional_loglik_replicate = ll_replicate,
  conditional_loglik_difference = ll_replicate - ll_block,
  marginal_loglik_block = as.numeric(block$marginal_loglik[1L]),
  marginal_loglik_replicate = as.numeric(fit_replicate$mlik[1L]),
  marginal_loglik_difference = as.numeric(fit_replicate$mlik[1L] - block$marginal_loglik[1L]),
  block_mode_status = block$mode_status,
  replicate_mode_status = fit_replicate$mode$mode.status,
  block_warnings = length(block$warnings),
  replicate_warnings = length(fit_replicate$misc$warnings))
write.csv(z, file.path(out, "comparison.csv"), row.names = FALSE)
saveRDS(list(input_rows = j, input_genes = seq_len(p), block = block,
  replicate = replicate_fit, comparison = z), file.path(out, "comparison.rds"))
