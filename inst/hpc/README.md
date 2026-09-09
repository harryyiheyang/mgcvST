# Replace mgcvST on HPC

Replace the existing package in its original library. Do not create another
project `Rlib` merely to preserve the old mgcvST. Run installation in a fresh R
session, then restart the analysis and its worker pool. Never overwrite a package
while active jobs are using that installation.

From the WGCNA payload directory, with this installer copied there:

```bash
env -u R_LIBS_USER Rscript --vanilla install_mgcvST.R packages/mgcvST_0.0.1.tar.gz
env -u R_LIBS_USER Rscript --vanilla check_worker_libraries.R
env -u R_LIBS_USER Rscript --vanilla check.R
```

The `env -u` commands undo the temporary project-Rlib override used in the earlier
deployment. If the site's normal configuration requires a custom `R_LIBS_USER`,
restore that original value instead. The installer prints the exact original
library selected and stops if it is not writable; it does not silently create a
second installation. Inspect `DONE (mgcvST)` before running the checks.

The worker check starts two fresh SOCK workers without fitting any genes. It
requires identical package paths and versions in the master and workers, and
checks availability of `model.set`. The six-gene numerical preflight is a
separate test; do not launch the full simulation until it passes.

Before starting a SOCK pool, pass the selected library paths to fresh sessions:

```r
Sys.setenv(R_LIBS_USER = paste(.libPaths(), collapse = .Platform$path.sep))
```

Do not prepend the old project `Rlib` in the analysis loader. A master's
`.libPaths()` change alone is not inherited by SOCK workers. A local reproduction
showed that the legacy mesh-based SPDE prediction method, applied to a newer
basis-backed model, fails with `array(STATS, dims[perm])` because its expected
`object$transform$center` is NULL. This matches the reported HPC error; worker
paths and the numerical preflight must still confirm the repair on that host.

These changes concern installation and worker startup only. No covariance,
WGCNA, simulation, or fitting parameters are changed.
