# Step 1 of the INLA validation of the PCAlearning saddlepoint route. Run it
# with the mgcvST build that produced the stored MAGIC fit and the stored exact
# spectra (the 0.995-coverage observation basis, q = 1404): it saves that
# basis and the mean of the fitted mean of every gene (mu_bar), which a fit of
# the current version stores at estimation.
#
# Rscript spa-validation-inla-prepare.R <old library> <out.rds> [threads] [fit.rds]
args <- commandArgs(TRUE)
if (nzchar(args[1L])) .libPaths(c(args[1L], .libPaths()))
out <- args[2L]
threads <- if (length(args) >= 3L) as.integer(args[3L]) else 2L
f <- if (length(args) >= 4L) args[4L] else
  "C:/Users/yxy1234/Downloads/magicST/downstream_data/output-no-celltype/inlaST-estimate.rds"
suppressMessages(library(mgcvST))
cat("mgcvST", as.character(packageVersion("mgcvST")), "\n")
fit <- mgcvST:::.inlast_sparse_prepare(readRDS(f))
basis <- mgcvST:::.inlast_sparse_observation_basis(fit)
cat("observation basis:", basis$kind, " rank", basis$rank, "\n")
scales <- mgcvST:::.mgcvst_pca_scales(fit, threads = threads)
saveRDS(list(basis = basis[c("coordinate", "basis", "rank", "kind")],
             mu_bar = scales$mu_bar, library = as.character(packageVersion("mgcvST"))),
        out)
cat("saved", out, "\n")
