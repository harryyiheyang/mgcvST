# Run in a fresh Rscript session, from the HPC payload directory.
if ("mgcvST" %in% loadedNamespaces()) {
  stop("Restart R before replacing mgcvST; do not reuse loaded namespaces or workers.")
}
args <- commandArgs(trailingOnly = TRUE)
src <- if (length(args)) args[[1L]] else "packages/mgcvST_0.0.1.tar.gz"
if (!file.exists(src)) stop("Source archive not found: ", src)
path <- find.package("mgcvST", quiet = TRUE)
lib <- if (length(path)) dirname(path) else .libPaths()[[1L]]
if (file.access(lib, 2L) != 0L) stop("Package library is not writable: ", lib)
print(list(source = normalizePath(src), overwrite_library = lib))
Sys.setenv(MAKEFLAGS = "-j4")
install.packages(src, repos = NULL, type = "source", lib = lib)

# install.packages can warn rather than stop after a failed source installation.
if (!dir.exists(file.path(lib, "mgcvST"))) stop("mgcvST installation is missing.")
print("Installation command finished. Check its DONE message, then verify in a fresh R session.")
