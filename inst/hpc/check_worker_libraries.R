# This check creates two fresh workers but fits no models.
Sys.setenv(R_LIBS_USER = paste(.libPaths(), collapse = .Platform$path.sep))
pkg <- c("mgcvST", "mgcv", "CppMatrix", "BiocParallel")
paths <- vapply(pkg, find.package, character(1L))
versions <- vapply(pkg, function(p) as.character(utils::packageVersion(p)), character(1L))
print(data.frame(package = pkg, path = paths, version = versions))
BP <- BiocParallel::bpstart(BiocParallel::SnowParam(workers = 2L, type = "SOCK"))
workers <- BiocParallel::bplapply(1:2, function(i, pkg) {
  list(paths = vapply(pkg, find.package, character(1L)),
    versions = vapply(pkg, function(p) as.character(utils::packageVersion(p)), character(1L)),
    model_api = exists("model.set", asNamespace("mgcvST"), inherits = FALSE))
}, pkg = pkg, BPPARAM = BP)
BiocParallel::bpstop(BP)
print(workers)
stopifnot(all(vapply(workers, function(x) identical(x$paths, paths) &&
  identical(x$versions, versions) && x$model_api, logical(1L))))
print("Worker package consistency PASS")
