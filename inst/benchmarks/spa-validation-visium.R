# Validation of the exact-moment saddlepoint route on the Visium-B pair
# universe of the shared-basis study (1,025 genes, q = 298, 25,853 pairs).
#
# Reference: nlp_full of pair_scheme_results.rds, the full-spectrum saddlepoint
# -log10 p-value of each pair at its stored score. The exact route (shared basis
# of k = 20 leading eigenvectors, four-moment two-node remainder) must stay
# within 0.002 of the reference up to p = 1e-30; k = q must reproduce it.
#
# Rscript spa-validation-visium.R <library> [threads] [study_dir]
args <- commandArgs(TRUE)
lib <- args[1L]
threads <- if (length(args) >= 2L) as.integer(args[2L]) else 2L
study <- if (length(args) >= 3L) args[3L] else
  "C:/Users/yxy1234/Downloads/magicST/paper_workspace/05_analysis"
if (!is.na(lib) && nzchar(lib)) .libPaths(c(lib, .libPaths()))
suppressMessages(library(mgcvST))
d <- file.path(study, "pair_spa_shared_basis")

PS <- readRDS(file.path(d, "pair_scheme_results.rds"))
GU <- readRDS(file.path(d, "gene_universe.rds"))
SF <- readRDS(file.path(study, "pair_liu_davies_tail", "spectrum_factors.rds"))
EF <- readRDS(file.path(d, "extra_gene_factors.rds"))
Fac <- c(SF$Fac, EF$Fac)
rm(SF, EF)
stopifnot(identical(length(Fac), length(GU$genes)))
q <- nrow(Fac[[1L]])
H <- lapply(Fac, function(F) tcrossprod(F))
rm(Fac)
n <- length(H)
scale <- vapply(H, function(M) max(abs(M)), numeric(1L))

left <- match(PS$feature1, GU$genes)
right <- match(PS$feature2, GU$genes)
stopifnot(!anyNA(left), !anyNA(right))
ord <- order(left, right)
x <- PS$q / sqrt(scale[left] * scale[right])
a <- matrix(0, q, n)

shared_basis <- function(k) {
  S <- mgcvST:::mgcvst_pair_basis_sum_cpp(H)
  eigen(S, symmetric = TRUE)$vectors[, seq_len(k), drop = FALSE]
}
run <- function(V, rows) {
  G <- mgcvST:::mgcvst_pair_basis_cpp(H, V, threads)
  z <- mgcvST:::mgcvst_pair_spa_cpp(H, G, a, left[rows], right[rows], threads, 4L, x[rows])
  -z$log_p_two_sided / log(10)
}

report <- function(delta, reference, label) {
  bins <- cut(reference, c(0, 4, 6, 8, 10, 12, 15, 20, 30, Inf), include.lowest = TRUE)
  tab <- do.call(rbind, lapply(split(seq_along(delta), bins), function(i) {
    if (!length(i)) return(NULL)
    data.frame(bin = as.character(bins[i[1L]]), n = length(i),
               max_abs = max(abs(delta[i])))
  }))
  cat("\n", label, "\n", sep = "")
  print(tab, row.names = FALSE, digits = 3)
  below30 <- reference <= 30
  cat("max |delta -log10 p| for p >= 1e-30:", signif(max(abs(delta[below30])), 3), "\n")
  invisible(tab)
}

t0 <- proc.time()[["elapsed"]]
rows <- ord
nlp <- numeric(nrow(PS))
nlp[rows] <- run(shared_basis(20L), rows)
cat("exact route, k = 20:", nrow(PS), "pairs in",
    round(proc.time()[["elapsed"]] - t0), "s,", threads, "threads\n")
report(nlp - PS$nlp_full, PS$nlp_full, "k = 20 against nlp_full (max |delta| by reference bin)")

# k = q: V = I, G_g = H_g^{1/2}; the pair spectrum is complete and the remainder
# vanishes, so the route reproduces the full-spectrum saddlepoint.
set.seed(1)
sub <- sort(sample(nrow(PS), 600L))
rows_q <- sub[order(left[sub], right[sub])]
full <- numeric(nrow(PS))
full[rows_q] <- run(diag(q), rows_q)
d_full <- full[rows_q] - PS$nlp_full[rows_q]
cat("\nk = q = ", q, " on ", length(rows_q), " pairs: max |delta| = ",
    signif(max(abs(d_full)), 3), "\n", sep = "")
