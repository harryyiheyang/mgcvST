#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1", BLIS_NUM_THREADS = "1")
if (!requireNamespace("pkgload", quietly = TRUE)) stop("pkgload is required to load the current checkout.")
pkgload::load_all(".", export_all = FALSE, quiet = TRUE)
library(BiocParallel)

out <- Sys.getenv("MGCVST_MAGIC_PAIR_OUTPUT", "artifacts/magic-inla-pairwise-liu")
run.fit <- identical(Sys.getenv("MGCVST_MAGIC_RUN_FITS", "no"), "yes")
n.fit <- as.integer(Sys.getenv("MGCVST_MAGIC_N_FIT", "3"))
if (is.na(n.fit) || n.fit < 1L || n.fit > 16L) stop("MGCVST_MAGIC_N_FIT must be an integer from 1 through 16.")
if (n.fit > 3L && !identical(Sys.getenv("MGCVST_MAGIC_APPROVED_GT3"), "yes")) stop("More than three new INLA fits require MGCVST_MAGIC_APPROVED_GT3=yes.")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out, "matrix-f32"), recursive = TRUE, showWarnings = FALSE)

magic.file <- "artifacts/datasets/MAGIC/MAGIC.rds"
csv.file <- "C:/Users/yxy1234/Downloads/magicST/magic_inlast_bench/expression_3d_100genes.csv"
fit.dir <- "artifacts/magic-slide-exploration/fits"
stress.file <- "artifacts/inla-public-parallel-stress/serial/fit.rds"
if (!file.exists(magic.file) || !file.exists(csv.file) || !file.exists(stress.file)) stop("One or more registered inputs are unavailable.")
desc <- read.dcf("DESCRIPTION")[1L, "Version"]
if (as.character(utils::packageVersion("mgcvST")) != desc || !exists("mgcvst_pair_trace_powers_cpp", envir = asNamespace("mgcvST"), inherits = FALSE)) stop("The loaded mgcvST namespace does not match the current checkout.")

MAGIC <- readRDS(magic.file)
d <- MAGIC$covariates
mesh <- fmesher::fm_mesh_3d(loc = MAGIC$meshes$native3d$loc, tv = MAGIC$meshes$native3d$tv)
kappa <- MAGIC$meshes$native3d$contract$kappa_fixed
ids <- colnames(MAGIC$expression)
Y <- t(MAGIC$expression[, ids, drop = FALSE])
rownames(Y) <- ids
if (n.fit > length(ids)) {
  X <- data.table::fread(csv.file, check.names = FALSE, data.table = FALSE)
  gene <- as.character(X[[1L]])
  point.id <- names(X)[-1L]
  ix <- match(d$point_id, point.id)
  if (anyNA(ix)) stop("The 100-gene CSV does not cover every MAGIC point_id.")
  X <- as.matrix(X[, -1L, drop = FALSE])
  storage.mode(X) <- "double"
  rownames(X) <- gene
  X <- X[, ix, drop = FALSE]
  keep <- setdiff(rownames(X), ids)
  keep <- keep[order(rowMeans(X[keep, , drop = FALSE]))]
  take <- keep[unique(round(seq(1L, length(keep), length.out = n.fit - length(ids))))]
  Y <- rbind(Y, X[take, , drop = FALSE])
}
Y <- Y[seq_len(n.fit), , drop = FALSE]
write.csv(data.frame(source = c("MAGIC.rds", "expression_3d_100genes.csv", "historical standalone fits", "Visium interface fit"),
  path = c(normalizePath(magic.file), normalizePath(csv.file), normalizePath(fit.dir), normalizePath(stress.file)),
  role = c("current data and native 3D mesh", "additional candidate counts", "integral constraint: context only, never production M", "2D q=298 interface-only check, never MAGIC M")), file.path(out, "source-hierarchy.csv"), row.names = FALSE)

flat <- list(prior = "flat", param = numeric(), initial = 0)
S <- mgcvST::inlaST.set(response ~ offset(log(exposure)), d, mesh = mesh, kappa = kappa,
  coordinates = c("x_mm", "y_mm", "z_mm"), family = mgcv::nb())
if (!identical(S$mean_constraint, "observation")) stop("The current native model did not retain the observation-mean constraint.")
write.csv(data.frame(observations = nrow(d), mesh_nodes = ncol(S$inla_spec$random[[1L]]$Q), kappa = kappa,
  constraint = "g = A^T 1/n", nuisance = "expected-curvature -U Vp U^T", disk_budget_bytes = 1000000000,
  io_scope = "sub-1GB files; OS-cache first/repeat read only; no cache clearing; no 154GB cold-read extrapolation"), file.path(out, "production-contract.csv"), row.names = FALSE)
if (!run.fit) quit(save = "no", status = 0L)

ctl <- list(precision_prior = flat, nb_size_prior = flat, num_threads = 1L, keep_fit = FALSE)
F <- mgcvST::inlaST.estimate(Y, S, feature_id = rownames(Y), control = ctl, diagnostics = TRUE,
  retain_smooth = TRUE, BPPARAM = SerialParam(), threads = 1L)
if (!all(F$diagnostics$converged)) stop("At least one current MAGIC INLA fit did not converge.")
z <- mgcvST:::.inlast_sparse_batch(F, seq_len(nrow(Y)), threads = 1L)
bad <- vapply(z, function(x) !is.null(x$error) && length(x$error) == 1L && nzchar(x$error), logical(1L))
if (any(bad)) stop("Sparse score construction failed: ", paste(F$feature_id[bad], vapply(z[bad], `[[`, character(1L), "error"), collapse = " | "))
M <- lapply(z, `[[`, "M")
names(M) <- F$feature_id
a <- lapply(z, `[[`, "a")
q <- ncol(M[[1L]])
W <- 1 / F$working_variance
checkpoint.file <- file.path(out, "production-double-checkpoint.rds")
saveRDS(list(feature_id = F$feature_id, M = M, a = a, W = W,
  diagnostics = F$diagnostics, lambda = F$lambda, dispersion = F$dispersion,
  family_parameters = F$family_parameters), checkpoint.file, compress = FALSE)
checkpoint.bytes <- file.info(checkpoint.file)$size
if (checkpoint.bytes > 1000000000) stop("The double score checkpoint exceeded the 1 GB disk budget.")
R1 <- list()
for (j in seq_len(ncol(W))) {
  alpha <- sum(W[, j] * W[, 1L]) / sum(W[, 1L]^2)
  R1[[j]] <- data.frame(feature_id = colnames(W)[j], reference = colnames(W)[1L], alpha = alpha,
    relative_frobenius_residual = norm(matrix(W[, j] - alpha * W[, 1L], ncol = 1L), "F") / norm(matrix(W[, j], ncol = 1L), "F"))
}
write.csv(do.call(rbind, R1), file.path(out, "working-curvature-proportionality.csv"), row.names = FALSE)

t0 <- proc.time()[["elapsed"]]
for (j in seq_along(M)) {
  con <- file(file.path(out, "matrix-f32", paste0(names(M)[j], ".bin")), "wb")
  writeBin(as.vector(M[[j]]), con, size = 4L)
  close(con)
}
write.seconds <- proc.time()[["elapsed"]] - t0
ff <- file.path(out, "matrix-f32", paste0(names(M), ".bin"))
bytes <- sum(file.info(ff)$size)
slot.dir <- file.path(out, "io-slots-repeated-content")
dir.create(slot.dir, recursive = TRUE, showWarnings = FALSE)
slot <- file.path(slot.dir, sprintf("slot-%02d.bin", 1:16))
t0 <- proc.time()[["elapsed"]]
for (j in seq_along(slot)) file.copy(ff[(j - 1L) %% length(ff) + 1L], slot[j], overwrite = TRUE)
slot.write.seconds <- proc.time()[["elapsed"]] - t0
slot.bytes <- sum(file.info(slot)$size)
if (bytes + slot.bytes + checkpoint.bytes > 1000000000) stop("The checkpoint, float32 matrix and repeated I/O slots exceeded the 1 GB disk budget.")
read.m <- function(id) {
  con <- file(file.path(out, "matrix-f32", paste0(id, ".bin")), "rb")
  x <- readBin(con, numeric(), n = q * q, size = 4L)
  close(con)
  matrix(x, q, q)
}
read.pass <- function(order) {
  t0 <- proc.time()[["elapsed"]]
  for (j in order) {
    con <- file(slot[j], "rb")
    readBin(con, numeric(), n = q * q, size = 4L)
    close(con)
  }
  proc.time()[["elapsed"]] - t0
}
set.seed(20260922L)
ord <- sample(seq_along(slot))
seq.first <- read.pass(seq_along(slot)); seq.repeat <- read.pass(seq_along(slot))
random.first <- read.pass(ord); random.repeat <- read.pass(ord)
M32 <- lapply(names(M), read.m)
names(M32) <- names(M)
write.csv(data.frame(operation = c("write_actual_matrices", "write_16_repeated_content_slots", "sequential_first_application_read", "sequential_repeat_application_read", "random_first_application_read_after_sequential", "random_repeat_application_read"), bytes = c(bytes, slot.bytes, rep(slot.bytes, 4L)),
  seconds = c(write.seconds, slot.write.seconds, seq.first, seq.repeat, random.first, random.repeat),
  gib_per_second = c(bytes, slot.bytes, rep(slot.bytes, 4L)) / c(write.seconds, slot.write.seconds, seq.first, seq.repeat, random.first, random.repeat) / 1024^3,
  files = c(length(ff), rep(length(slot), 5L)), content = c("independent fitted matrices", rep("16 repeated-content I/O slots", 5L)),
  cold_device_io_isolated = FALSE), file.path(out, "matrix-io.csv"), row.names = FALSE)

R2 <- list()
for (j in seq_along(M)) {
  ev <- eigen(M[[j]], symmetric = TRUE, only.values = TRUE)$values
  alpha <- sum(M[[j]] * M[[1L]]) / sum(M[[1L]]^2)
  R2[[j]] <- data.frame(feature_id = names(M)[j], reference = names(M)[1L], alpha = alpha,
    relative_frobenius_residual = norm(M[[j]] - alpha * M[[1L]], "F") / norm(M[[j]], "F"),
    trace = sum(ev), effective_rank = sum(ev)^2 / sum(ev^2), minimum_eigenvalue = min(ev))
}
write.csv(do.call(rbind, R2), file.path(out, "matrix-spectrum.csv"), row.names = FALSE)

R3 <- list()
for (test in seq_along(M)) {
  train <- setdiff(seq_along(M), test)
  E <- eigen(Reduce(`+`, M[train]) / length(train), symmetric = TRUE)$vectors
  for (rank in c(64L, 128L, 256L)) {
    V <- E[, seq_len(min(rank, q)), drop = FALSE]
    B <- crossprod(V, M[[test]] %*% V)
    R3[[length(R3) + 1L]] <- data.frame(train = paste(names(M)[train], collapse = ";"), test = names(M)[test], rank = ncol(V),
      trace_tail = sum(diag(M[[test]])) - sum(diag(B)), trace_tail_ratio = 1 - sum(diag(B)) / sum(diag(M[[test]])),
      frobenius_tail = sqrt(sum(M[[test]]^2) - sum(B^2)), frobenius_tail_ratio = sqrt(sum(M[[test]]^2) - sum(B^2)) / norm(M[[test]], "F"),
      projected_offdiagonal_energy = (sum(B^2) - sum(diag(B)^2)) / sum(B^2))
  }
}
write.csv(do.call(rbind, R3), file.path(out, "common-basis-leave-one-out.csv"), row.names = FALSE)

P <- t(combn(seq_along(M), 2L))
R4 <- list()
for (j in seq_len(nrow(P))) {
  i1 <- P[j, 1L]; i2 <- P[j, 2L]
  C <- M[[i1]] %*% M[[i2]] - M[[i2]] %*% M[[i1]]
  R4[[j]] <- data.frame(feature1 = names(M)[i1], feature2 = names(M)[i2], commutator_relative_frobenius = norm(C, "F") / norm(M[[i1]] %*% M[[i2]], "F"))
}
write.csv(do.call(rbind, R4), file.path(out, "commutators.csv"), row.names = FALSE)

T0 <- mgcvST:::mgcvst_pair_trace_powers_cpp(M, P, maxPower = 4L, threads = 1L)
Tq <- mgcvST:::mgcvst_pair_trace_powers_cpp(M32, P, maxPower = 4L, threads = 1L)
R5 <- list()
for (j in seq_len(nrow(P))) {
  u <- sum(a[[P[j, 1L]]] * a[[P[j, 2L]]])
  U <- c(u, 4 * sqrt(T0[j, 1L]), 6 * sqrt(T0[j, 1L]), 8 * sqrt(T0[j, 1L]))
  L0 <- mgcvST:::.liu_squared_score_moments(U, T0[j, 1L], T0[j, 2L], T0[j, 3L], T0[j, 4L])
  Lq <- mgcvST:::.liu_squared_score_moments(U, Tq[j, 1L], Tq[j, 2L], Tq[j, 3L], Tq[j, 4L])
  R5[[j]] <- data.frame(feature1 = names(M)[P[j, 1L]], feature2 = names(M)[P[j, 2L]], U_kind = c("observed", "z4", "z6", "z8"), U = U,
    c1_relative_error = (Tq[j, 1L] - T0[j, 1L]) / T0[j, 1L], c2_relative_error = (Tq[j, 2L] - T0[j, 2L]) / T0[j, 2L],
    c3_relative_error = (Tq[j, 3L] - T0[j, 3L]) / T0[j, 3L], c4_relative_error = (Tq[j, 4L] - T0[j, 4L]) / T0[j, 4L],
    liu_relative_error = (Lq$p_value - L0$p_value) / L0$p_value, delta_log10_p = log10(Lq$p_value) - log10(L0$p_value))
}
write.csv(do.call(rbind, R5), file.path(out, "float32-liu-error.csv"), row.names = FALSE)

R6 <- list()
for (r in 1:5) for (power in if (r %% 2L) c(1L, 2L, 4L) else c(4L, 2L, 1L)) {
  t0 <- proc.time()[["elapsed"]]
  mgcvST:::mgcvst_pair_trace_powers_cpp(M, P, maxPower = power, threads = 1L)
  R6[[length(R6) + 1L]] <- data.frame(repetition = r, max_power = power, seconds = proc.time()[["elapsed"]] - t0)
}
R6 <- do.call(rbind, R6)
write.csv(R6, file.path(out, "trace-power-repeats.csv"), row.names = FALSE)
T <- aggregate(seconds ~ max_power, R6, median)
T$c2_reduction_median_seconds <- ifelse(T$max_power == 2L, T$seconds - T$seconds[T$max_power == 1L], NA_real_)
T$c3_c4_median_seconds <- ifelse(T$max_power == 4L, T$seconds - T$seconds[T$max_power == 2L], NA_real_)
write.csv(T, file.path(out, "trace-power-timing.csv"), row.names = FALSE)
g <- gc()
p <- as.numeric(system2("powershell", c("-NoProfile", "-Command", paste0("(Get-Process -Id ", Sys.getpid(), ").PeakWorkingSet64")), stdout = TRUE))
write.csv(data.frame(double_checkpoint_bytes = checkpoint.bytes, matrix_float32_bytes = bytes, repeated_io_slot_bytes = slot.bytes,
  output_disk_bytes = sum(file.info(list.files(out, recursive = TRUE, full.names = TRUE))$size),
  r_heap_gc_max_mib_approx = sum(g[, "max used"] * c(56, 8)) / 1024^2,
  windows_peak_working_set_bytes = p), file.path(out, "resource-use.csv"), row.names = FALSE)
capture.output(sessionInfo(), file = file.path(out, "session-info.txt"))
