# Public 2D score calibration. Arguments: pair/marginal, mean, first replicate, last replicate.
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
library(mgcvST)
library(mgcv)
options(warn = 2)
a <- commandArgs(TRUE)
stopifnot(length(a) == 4L, as.character(packageVersion("mgcvST")) == "0.0.1.9006")
kind <- a[1L]
mu0 <- as.numeric(a[2L])
rr <- seq.int(as.integer(a[3L]), as.integer(a[4L]))
case <- paste(kind, mu0, sep = "-")
out <- file.path("artifacts/inla-stress-calibration/null2d", case)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
d <- readRDS("artifacts/inla-stress-calibration/null2d-input.rds")
flat <- list(prior = "flat", param = numeric(), initial = 0)
ctl <- list(precision_prior = flat, nb_size_prior = flat, num_threads = 1L)
n <- nrow(d$F)
p <- if (kind == "pair") 2L else 1L
ids <- paste0("g", seq_len(p))
for (r in rr) {
  file <- file.path(out, sprintf("rep-%04d.rds", r))
  if (file.exists(file)) stop("Replicate already exists: ", file)
  seed <- 202610000L + if (kind == "pair") 0L else 10000L
  seed <- seed + if (mu0 < 1) 0L else 20000L
  seed <- seed + r
  set.seed(seed)
  eta <- matrix(0, n, p)
  if (kind == "pair") eta <- d$F %*% matrix(rnorm(ncol(d$F) * p), ncol = p)
  beta <- log(mu0) - log(mean(exp(d$d$off + if (kind == "pair") 0.5 * rowSums(d$F^2) else 0)))
  Y <- matrix(NA_integer_, p, n, dimnames = list(ids, NULL))
  for (j in seq_len(p)) Y[j, ] <- rnbinom(n, mu = exp(beta + d$d$off + eta[, j]), size = 2)
  t0 <- proc.time()[["elapsed"]]
  I <- inlaST.estimate(Y, d$S, control = ctl, diagnostics = TRUE,
    retain_marginal = TRUE)
  elapsed <- proc.time()[["elapsed"]] - t0
  rows <- list()
  for (cal in c("davies", "liu")) {
    if (kind == "pair") {
      z <- mgcvST.test(I, pairs = matrix(ids, nrow = 1L), calibration = cal)$results
      rows[[cal]] <- data.frame(case = case, replicate = r, seed = seed,
        calibration = cal, p_value = z$p_two_sided, statistic = z$signed_score,
        information = z$information, error = z$error_message)
    } else {
      z <- mgcvST.marginal(I, calibration = cal)
      rows[[cal]] <- data.frame(case = case, replicate = r, seed = seed,
        calibration = cal, p_value = z$p_value, statistic = z$statistic,
        information = NA_real_, error = z$error_message)
    }
  }
  z <- do.call(rbind, rows)
  H <- data.frame(feature = ids, converged = I$diagnostics$converged,
    error = I$diagnostics$error_message,
    tau = I$smoothing_parameters[, 1L] / I$dispersion,
    size = vapply(I$family_parameters, function(x) if (length(x)) as.numeric(x[1L]) else NA_real_, numeric(1L)))
  saveRDS(list(rows = z, hyper = H, seconds = elapsed, Y = Y), file)
  write.csv(z, sub("rds$", "csv", file), row.names = FALSE)
  rm(I)
  gc()
}
