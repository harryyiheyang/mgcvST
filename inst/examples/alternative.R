library(BiocParallel)
library(mgcv)
library(mgcvST)

workers <- as.integer(Sys.getenv("MGCVST_EXAMPLE_WORKERS", "20"))
B.power <- as.integer(Sys.getenv("MGCVST_EXAMPLE_POWER_REPS", "100"))
out.dir <- Sys.getenv("MGCVST_EXAMPLE_OUTPUT", ".")
if (!dir.exists(out.dir)) dir.create(out.dir, recursive = TRUE)
rho <- seq(0, 0.5, length.out = 100L)
kappa <- 0.1
cutoff <- 0.999
ctrl <- gam.control(maxit = 50L, nthreads = 1L, ncv.threads = 1L)

if (!is.finite(workers) || workers < 1L || !is.finite(B.power) ||
    B.power < 1L || !dir.exists(out.dir)) {
  stop("workers, power simulation count, and output directory are invalid.")
}

fit_pair <- function(D, gene1, gene2, mesh, pc = FALSE) {
  D1 <- D$covariates
  D2 <- D$covariates
  D1$response <- D$expression[[gene1]]
  D2$response <- D$expression[[gene2]]
  if (pc) {
    f1 <- gam(
      response ~ offset(offset0) + s(
        x, y, bs = "spdePC", xt = list(mesh = mesh, pc_cutoff = cutoff),
        sp = c(-1, kappa)
      ), data = D1, family = nb(link = "log"), method = "REML", control = ctrl
    )
    f2 <- gam(
      response ~ offset(offset0) + s(
        x, y, bs = "spdePC", xt = list(mesh = mesh, pc_cutoff = cutoff),
        sp = c(-1, kappa)
      ), data = D2, family = nb(link = "log"), method = "REML", control = ctrl
    )
  } else {
    f1 <- gam(
      response ~ offset(offset0) + s(
        x, y, bs = "spde", xt = list(mesh = mesh), sp = c(-1, kappa)
      ), data = D1, family = nb(link = "log"), method = "REML", control = ctrl
    )
    f2 <- gam(
      response ~ offset(offset0) + s(
        x, y, bs = "spde", xt = list(mesh = mesh), sp = c(-1, kappa)
      ), data = D2, family = nb(link = "log"), method = "REML", control = ctrl
    )
  }
  list(fit1 = f1, fit2 = f2)
}

psd_root <- function(H) {
  H <- (H + t(H)) / 2
  E <- CppMatrix::matrixEigen(H)
  d <- as.numeric(E$values)
  V <- as.matrix(E$vectors)
  tol <- sqrt(.Machine$double.eps) * max(1, max(abs(d)))
  if (min(d) < -tol) stop("Score covariance matrix is not positive semidefinite.")
  R <- CppMatrix::matrixMultiply(sweep(V, 2L, sqrt(pmax(d, 0)), "*"), t(V))
  (R + t(R)) / 2
}

liu_p <- function(U, moments) {
  mgcvST:::.liu_squared_score_moments(
    U, moments[1L], moments[2L], moments[3L], moments[4L]
  )$p_value
}

power_chunk <- function(ids, R1, R2, moments, rho, B.power) {
  set.seed(2026082600L + ids[1L])
  q <- nrow(R1)
  out <- vector("list", length(ids))
  for (j in seq_along(ids)) {
    i <- ids[j]
    Z1 <- matrix(rnorm(q * B.power), q, B.power)
    Z2 <- matrix(rnorm(q * B.power), q, B.power)
    U <- colSums(
      (R1 %*% Z1) *
        (R2 %*% (rho[i] * Z1 + sqrt(1 - rho[i]^2) * Z2))
    )
    out[[j]] <- data.frame(
      rho_index = i, rho = rho[i], replicate = seq_len(B.power),
      p_value = liu_p(U, moments)
    )
  }
  do.call(rbind, out)
}

init_worker <- function(i) {
  Sys.setenv(
    OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1", BLIS_NUM_THREADS = "1",
    VECLIB_MAXIMUM_THREADS = "1", RCPP_PARALLEL_NUM_THREADS = "1"
  )
  library(mgcvST)
  data.frame(worker = i, pid = Sys.getpid())
}

data(MISO_E13, package = "mgcvST")
data(Visium_B, package = "mgcvST")
S <- list(
  MISO_spde = c(list(slice = "MISO_E13", representation = "spde",
                     pair = "Mapt--Map1b"),
                fit_pair(MISO_E13, "Mapt", "Map1b", MISO_E13$meshes$spde)),
  MISO_spdePC_G999 = c(list(slice = "MISO_E13", representation = "spdePC_G999",
                            pair = "Mapt--Map1b"),
                       fit_pair(MISO_E13, "Mapt", "Map1b",
                                MISO_E13$meshes$spdePC_g999, pc = TRUE)),
  Visium_spde = c(list(slice = "Visium_B", representation = "spde",
                       pair = "mt-co3--BRAFhuman"),
                  fit_pair(Visium_B, "mt_co3", "BRAFhuman", Visium_B$meshes$spde)),
  Visium_spdePC_G999 = c(list(slice = "Visium_B", representation = "spdePC_G999",
                              pair = "mt-co3--BRAFhuman"),
                         fit_pair(Visium_B, "mt_co3", "BRAFhuman",
                                  Visium_B$meshes$spdePC_g999, pc = TRUE))
)

R.geometry <- list()
R.power <- list()
for (nm in names(S)) {
  z <- rkhs_covariance_score_irls(S[[nm]]$fit1, S[[nm]]$fit2, method = "liu")
  S[[nm]]$R1 <- psd_root(z$summary1$H)
  S[[nm]]$R2 <- psd_root(z$summary2$H)
  S[[nm]]$moments <- z$moments
  sm <- S[[nm]]$fit1$smooth[[1L]]
  Q <- if (inherits(sm, "spdePC.smooth")) sm$pc_score_Q else sm$S[[1L]]
  R.geometry[[nm]] <- data.frame(
    scenario = nm, slice = S[[nm]]$slice, representation = S[[nm]]$representation,
    pair = S[[nm]]$pair, n_spots = length(S[[nm]]$fit1$y),
    fit_dimension = sm$last.para - sm$first.para + 1L,
    score_dimension = ncol(z$summary1$H), score_Q_rank = as.integer(Matrix::rankMatrix(Q)),
    information = z$information, effective_rank = z$effective_rank,
    kappa = kappa, pc_cutoff = if (inherits(sm, "spdePC.smooth")) sm$pc_cutoff else NA_real_
  )
}

BPPARAM <- SnowParam(
  workers = workers, type = "SOCK", tasks = 0L, stop.on.error = TRUE,
  progressbar = FALSE, RNGseed = 20260817L
)
BPPARAM <- bpstart(BPPARAM)
W <- do.call(rbind, bplapply(seq_len(workers), init_worker, BPPARAM = BPPARAM))
if (length(unique(W$pid)) != workers) stop("Worker initialization did not reach every worker.")

chunks <- split(seq_along(rho), rep(seq_len(workers), length.out = length(rho)))
for (nm in names(S)) {
  R.power[[nm]] <- do.call(rbind, bplapply(
    chunks, power_chunk, R1 = S[[nm]]$R1, R2 = S[[nm]]$R2,
    moments = S[[nm]]$moments, rho = rho, B.power = B.power, BPPARAM = BPPARAM
  ))
  R.power[[nm]]$scenario <- nm
  R.power[[nm]]$slice <- S[[nm]]$slice
  R.power[[nm]]$representation <- S[[nm]]$representation
  R.power[[nm]]$pair <- S[[nm]]$pair
}
BPPARAM <- bpstop(BPPARAM)

R.geometry <- do.call(rbind, R.geometry)
R.power <- do.call(rbind, R.power)
R.alternative <- list()
k <- 0L
for (nm in names(S)) {
  d <- R.power[R.power$scenario == nm, ]
  for (a in c(0.05, 0.01, 0.001)) {
    for (r in rho) {
      p <- d$p_value[d$rho == r]
      R.alternative[[k <- k + 1L]] <- data.frame(
        scenario = nm, slice = S[[nm]]$slice, representation = S[[nm]]$representation,
        pair = S[[nm]]$pair, rho = r, alpha = a, replicates = length(p),
        power = mean(p <= a)
      )
    }
  }
}
R.alternative <- do.call(rbind, R.alternative)

write.csv(R.geometry, file.path(out.dir, "alternative_geometry.csv"), row.names = FALSE)
write.csv(R.alternative, file.path(out.dir, "alternative_power.csv"), row.names = FALSE)
saveRDS(R.power, file.path(out.dir, "alternative_pvalues.rds"))
saveRDS(
  list(
    MISO_spdePC_G999 = S$MISO_spdePC_G999$fit1$smooth[[1L]]$pc_score_Q,
    Visium_spdePC_G999 = S$Visium_spdePC_G999$fit1$smooth[[1L]]$pc_score_Q
  ), file.path(out.dir, "alternative_score_Q.rds")
)
