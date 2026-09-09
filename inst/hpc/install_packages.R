options(repos = c(CRAN = "https://cloud.r-project.org"), timeout = 1200)
if ("mgcvST" %in% loadedNamespaces()) stop("Run the installer in a fresh R session.")
path <- find.package("mgcvST", quiet = TRUE)
lib <- if (length(path)) dirname(path) else .libPaths()[[1L]]
if (file.access(lib, 2L) != 0L) stop("Package library is not writable: ", lib)
Sys.setenv(R_LIBS_USER = paste(.libPaths(), collapse = .Platform$path.sep), MAKEFLAGS = "-j1", OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
cran <- c("BiocManager", "Rcpp", "RcppArmadillo", "Matrix", "mgcv", "data.table",
  "dynamicTreeCut", "fastcluster", "mclust", "jsonlite", "ps")
missing <- cran[!vapply(cran, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing, lib = lib, Ncpus = 4L,
  dependencies = c("Depends", "Imports", "LinkingTo"))
bioc <- c("BiocParallel", "impute", "preprocessCore", "GO.db", "AnnotationDbi")
missing <- bioc[!vapply(bioc, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) BiocManager::install(missing, lib = lib, ask = FALSE,
  update = FALSE, Ncpus = 4L)
if (!requireNamespace("WGCNA", quietly = TRUE)) install.packages("WGCNA", lib = lib,
  Ncpus = 4L, dependencies = c("Depends", "Imports", "LinkingTo"))
if (!requireNamespace("CppMatrix", quietly = TRUE)) {
  install.packages("packages/CppMatrix_0.1.0.tar.gz", repos = NULL, type = "source", lib = lib)
}
install.packages("packages/mgcvST_0.0.1.tar.gz", repos = NULL, type = "source", lib = lib)
needed <- c(cran, bioc, "WGCNA", "CppMatrix", "mgcvST")
ok <- vapply(needed, requireNamespace, logical(1), quietly = TRUE)
print(ok)
if (!all(ok)) stop("Some packages failed to install; inspect the installation error above.")
if (!exists("model.set", asNamespace("mgcvST"))) {
  stop("The installed mgcvST lacks model.set; inspect the installation log and restart R.")
}
print("Installation check PASS")
