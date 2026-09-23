#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
if (!requireNamespace("pkgload", quietly = TRUE)) stop("pkgload is required to load the current checkout.")
pkgload::load_all(".", export_all = FALSE, quiet = TRUE)

out <- "artifacts/magic-inla-pairwise-liu"
checkpoint.file <- file.path(out, "production-double-checkpoint.rds")
csv.file <- "C:/Users/yxy1234/Downloads/magicST/magic_inlast_bench/expression_3d_100genes.csv"
magic.file <- "artifacts/datasets/MAGIC/MAGIC.rds"
if (!file.exists(checkpoint.file) || !file.exists(csv.file) || !file.exists(magic.file)) stop("The checkpoint or registered candidate input is unavailable.")
Z <- readRDS(checkpoint.file)
M <- Z$M
a <- Z$a
if (length(M) != 3L || length(a) != 3L) stop("This full-basis check requires exactly three retained MAGIC score states.")
q <- ncol(M[[1L]])
if (any(vapply(M, function(x) !all(dim(x) == c(q, q)), logical(1L)))) stop("The retained matrices have incompatible dimensions.")

P <- t(combn(seq_along(M), 2L))
V1 <- lapply(M, function(x) eigen(x, symmetric = TRUE)$vectors)
Vmean <- eigen(Reduce(`+`, M) / length(M), symmetric = TRUE)$vectors
R1 <- list()
DD <- vector("list", nrow(P))
for (basis in c(paste0("single_gene_", Z$feature_id), "third_gene_external", "three_gene_mean_optimistic")) for (j in seq_len(nrow(P))) {
  i1 <- P[j, 1L]
  i2 <- P[j, 2L]
  if (basis == "third_gene_external") {
    train <- setdiff(seq_along(M), c(i1, i2))
    V <- V1[[train]]
  } else if (basis == "three_gene_mean_optimistic") {
    train <- seq_along(M)
    V <- Vmean
  } else {
    train <- match(sub("^single_gene_", "", basis), Z$feature_id)
    V <- V1[[train]]
  }
  B1 <- crossprod(V, M[[i1]] %*% V)
  B2 <- crossprod(V, M[[i2]] %*% V)
  d1 <- diag(B1)
  d2 <- diag(B2)
  T0 <- as.numeric(mgcvST:::mgcvst_pair_trace_powers_cpp(M[c(i1, i2)], matrix(c(1L, 2L), nrow = 1L), maxPower = 4L, threads = 1L))
  Td <- vapply(1:4, function(k) sum((d1 * d2)^k), numeric(1L))
  U <- c(sum(a[[i1]] * a[[i2]]), -3 * sqrt(T0[1L]), -6 * sqrt(T0[1L]),
    -8 * sqrt(T0[1L]), 3 * sqrt(T0[1L]), 6 * sqrt(T0[1L]), 8 * sqrt(T0[1L]))
  L0 <- mgcvST:::.liu_squared_score_moments(U, T0[1L], T0[2L], T0[3L], T0[4L])
  Ld <- mgcvST:::.liu_squared_score_moments(U, Td[1L], Td[2L], Td[3L], Td[4L])
  p0.pos <- ifelse(U >= 0, L0$p_value / 2, 1 - L0$p_value / 2)
  p0.neg <- ifelse(U <= 0, L0$p_value / 2, 1 - L0$p_value / 2)
  pd.pos <- ifelse(U >= 0, Ld$p_value / 2, 1 - Ld$p_value / 2)
  pd.neg <- ifelse(U <= 0, Ld$p_value / 2, 1 - Ld$p_value / 2)
  R1[[length(R1) + 1L]] <- data.frame(
    basis = basis, training_features = paste(Z$feature_id[train], collapse = ";"),
    feature1 = Z$feature_id[i1], feature2 = Z$feature_id[i2],
    basis_rank = q, basis_spectral_tail = 0,
    diagonal_energy_1 = sum(d1^2), offdiagonal_energy_1 = pmax(0, sum(B1^2) - sum(d1^2)),
    diagonal_energy_2 = sum(d2^2), offdiagonal_energy_2 = pmax(0, sum(B2^2) - sum(d2^2)),
    offdiagonal_frobenius_ratio_1 = sqrt(pmax(0, sum(B1^2) - sum(d1^2))) / norm(B1, "F"),
    offdiagonal_frobenius_ratio_2 = sqrt(pmax(0, sum(B2^2) - sum(d2^2))) / norm(B2, "F"),
    c1_relative_error = (Td[1L] - T0[1L]) / T0[1L],
    c2_relative_error = (Td[2L] - T0[2L]) / T0[2L],
    c3_relative_error = (Td[3L] - T0[3L]) / T0[3L],
    c4_relative_error = (Td[4L] - T0[4L]) / T0[4L],
    U_kind = c("observed", "zminus3", "zminus6", "zminus8", "z3", "z6", "z8"), U = U,
    delta_log10_p_two_sided = log10(Ld$p_value) - log10(L0$p_value),
    delta_log10_p_positive = log10(pd.pos) - log10(p0.pos),
    delta_log10_p_negative = log10(pd.neg) - log10(p0.neg))
  if (basis == "third_gene_external") DD[[j]] <- list(d1 = d1, d2 = d2, T0 = T0, U = U)
}
write.csv(do.call(rbind, R1), file.path(out, "common-diagonalization-full.csv"), row.names = FALSE)

reps <- 1000L
t0 <- proc.time()[["elapsed"]]
for (r in seq_len(reps)) for (j in seq_len(nrow(P))) {
  Td <- vapply(1:4, function(k) sum((DD[[j]]$d1 * DD[[j]]$d2)^k), numeric(1L))
  mgcvST:::.liu_squared_score_moments(DD[[j]]$U, Td[1L], Td[2L], Td[3L], Td[4L])
}
diag.seconds <- (proc.time()[["elapsed"]] - t0) / reps
exact <- read.csv(file.path(out, "trace-power-repeats.csv"))
exact.seconds <- median(exact$seconds[exact$max_power == 4L])
write.csv(data.frame(method = c("exact_maxPower4", "external_basis_diagonal_reduction"),
  seconds_3_pairs = c(exact.seconds, diag.seconds),
  seconds_per_pair = c(exact.seconds / 3, diag.seconds / 3),
  note = c("existing three-pair single-thread median", "O(q) moments plus Liu after the full basis transforms; 1000 repetitions normalized to one three-pair pass")),
  file.path(out, "common-diagonalization-timing.csv"), row.names = FALSE)

MAGIC <- readRDS(magic.file)
d <- MAGIC$covariates
X <- data.table::fread(csv.file, check.names = FALSE, data.table = FALSE)
gene <- as.character(X[[1L]])
point.id <- names(X)[-1L]
ix <- match(d$point_id, point.id)
if (anyNA(ix)) stop("The 100-gene CSV does not cover every MAGIC point_id.")
X <- as.matrix(X[, -1L, drop = FALSE])
storage.mode(X) <- "double"
X <- X[, ix, drop = FALSE]
mu <- rowMeans(X)
vv <- rowMeans((X - mu)^2)
poisson.like <- vv <= mu
size.proxy <- rep(Inf, length(mu))
size.proxy[!poisson.like] <- mu[!poisson.like]^2 / (vv[!poisson.like] - mu[!poisson.like])
rate <- mu / mean(d$exposure)
W.proxy <- 1 / sweep(1 / pmax(rate %o% d$exposure, .Machine$double.xmin), 1L, 1 / size.proxy, "+")
W.cv <- apply(W.proxy, 1L, stats::sd) / rowMeans(W.proxy)
sec <- split(seq_len(nrow(d)), d$slice_order)
sv <- vapply(seq_len(nrow(X)), function(j) {
  m <- vapply(sec, function(ii) mean(X[j, ii]), numeric(1L))
  stats::var(m) / vv[j]
}, numeric(1L))
cut3 <- function(x) cut(rank(x, ties.method = "first"), breaks = quantile(rank(x, ties.method = "first"), c(0, 1/3, 2/3, 1)), include.lowest = TRUE, labels = c("low", "middle", "high"))
R2 <- data.frame(gene = gene, mean_count = mu, nb_size_proxy = size.proxy, poisson_like_proxy = poisson.like,
  working_curvature_proxy_cv = W.cv, section_mean_variance_fraction = sv,
  mean_stratum = cut3(mu), size_proxy_stratum = cut3(size.proxy),
  working_curvature_stratum = cut3(W.cv), section_stratum = cut3(sv),
  selection_role = "future 12-16 fit design: include low/middle/high and edge combinations of expression, NB-size proxy, W heterogeneity and section load; verify fitted NB size only after fitting")
write.csv(R2, file.path(out, "candidate-strata-design.csv"), row.names = FALSE)
capture.output(sessionInfo(), file = file.path(out, "common-diagonalization-session-info.txt"))
