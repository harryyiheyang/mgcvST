# Step 2 of the INLA validation of the PCAlearning saddlepoint route (run
# spa-validation-inla-prepare.R first). The 1,225 held-out pairs of 50 MAGIC
# genes have stored exact spectra (inla_check_exact_spectra.rds, q = 1404).
# For each pair and each reference level -log10 p in 5, 8, 12, 20, 30 the
# score x at which the full-spectrum saddlepoint of the exact spectrum has that
# level is fed to the production PCAlearning kernels (rank 20, k = 50), and
# the deviation of the route's -log10 p from the level is reported. The route
# must stay within 0.02 at 1e-12 and 0.05 at 1e-20.
#
# Rscript spa-validation-inla.R <library> <prepare.rds> [threads] [study_dir] [fit.rds]
args <- commandArgs(TRUE)
if (nzchar(args[1L])) .libPaths(c(args[1L], .libPaths()))
prepared_file <- args[2L]
threads <- if (length(args) >= 3L) as.integer(args[3L]) else 2L
study <- if (length(args) >= 4L) args[4L] else
  "C:/Users/yxy1234/Downloads/magicST/paper_workspace/05_analysis"
f <- if (length(args) >= 5L) args[5L] else
  "C:/Users/yxy1234/Downloads/magicST/downstream_data/output-no-celltype/inlaST-estimate.rds"
suppressMessages(library(mgcvST))
cat("mgcvST", as.character(packageVersion("mgcvST")), "\n")

SP <- readRDS(file.path(study, "pair_spa_shared_basis", "inla_check_exact_spectra.rds"))
old <- readRDS(prepared_file)
fit <- readRDS(f)
fit$mu_bar <- old$mu_bar
pairs <- cbind(SP$test[SP$pairs[, 1L]], SP$test[SP$pairs[, 2L]])
stopifnot(all(pairs[, 1L] < pairs[, 2L]))

t0 <- proc.time()[["elapsed"]]
fit <- mgcvST:::.inlast_sparse_prepare(fit)
used <- sort(unique(as.vector(pairs)))
prep <- mgcvST:::.mgcvst_pca_prepare(
  fit, used, old$basis, q = old$basis$rank, rank = 20L, n_per_cell = 3L, seed = 1L,
  k = 50L, threads = threads, verbose = TRUE)
cat("prepared in", round(proc.time()[["elapsed"]] - t0), "s; max |B'B - I| =",
    format(prep$basis_check, digits = 3), "\n")

# Reference: the full-spectrum saddlepoint of the exact spectrum, from the
# shared kernel (S = the whole spectrum, no remainder).
reference_nlp <- function(s, x) {
  tm <- vapply(1:4, function(r) sum(s^(2 * r)), numeric(1L))
  z <- mgcvST:::mgcvst_spa_cpp(x, matrix(s, ncol = 1L), matrix(tm, 4L, 1L), 4L, 1L)
  -z[, "log_p_two_sided"] / log(10)
}
level_x <- function(s, level) {
  hi <- max(s) * (level * log(10) + 60) + 20 * sqrt(sum(s^2))
  stats::uniroot(function(x) reference_nlp(s, x) - level,
                 c(1e-8 * sqrt(sum(s^2)), hi), tol = 1e-12)$root
}
levels <- c(5, 8, 12, 20, 30)
loc <- match(pairs, prep$used)
dim(loc) <- dim(pairs)
rows <- vector("list", nrow(pairs) * length(levels))
n <- 0L
for (p in seq_len(nrow(pairs))) {
  s <- SP$spectra[[p]]
  s <- s[s > 0]
  i <- loc[p, 1L]
  j <- loc[p, 2L]
  for (level in levels) {
    x <- level_x(s, level)
    # a_i = a_j = sqrt(x) e_1 gives the raw score U = x.
    A <- matrix(0, nrow(prep$A), 2L)
    A[1L, ] <- sqrt(x)
    out <- mgcvST:::mgcvst_pca_spa_pairs_cpp(
      A, prep$C[c(i, j), , drop = FALSE], prep$K2[, c(i, j), drop = FALSE],
      prep$T2, prep$R[, c(i, j), drop = FALSE], prep$scale[c(i, j)],
      1L, 2L, 1L)
    n <- n + 1L
    rows[[n]] <- c(pair = p, level = level, x = x,
                   delta = -out[1L, "logp_two_sided"] / log(10) - level,
                   kind = out[1L, "remainder_kind"], status = out[1L, "status"])
  }
}
D <- as.data.frame(do.call(rbind, rows))
cat("\nmax |delta -log10 p| of the PCAlearning route (r = 20, k = 50) against the",
    "exact-spectrum saddlepoint, by level:\n")
tab <- do.call(rbind, lapply(levels, function(l) {
  z <- D[D$level == l, ]
  data.frame(level = l, pairs = nrow(z), max_abs = max(abs(z$delta)),
             mean_abs = mean(abs(z$delta)),
             failed = sum(z$status != 0), one_node = sum(z$kind == 1),
             gaussian = sum(z$kind == 3))
}))
print(tab, row.names = FALSE, digits = 3)
cat("\nlimits: 0.02 at level 12, 0.05 at level 20 ->",
    if (tab$max_abs[tab$level == 12] <= 0.02 && tab$max_abs[tab$level == 20] <= 0.05)
      "met" else "NOT met", "\n")
