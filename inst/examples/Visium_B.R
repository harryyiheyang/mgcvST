library(BiocParallel)
library(mgcv)
library(mgcvST)

workers <- 3L
kappa <- 0.1
out.dir <- Sys.getenv("MGCVST_EXAMPLE_OUTPUT", ".")
if (!dir.exists(out.dir)) dir.create(out.dir, recursive = TRUE)

data(Visium_B, package = "mgcvST")
genes <- c("mt_co3", "mt_co2", "BRAFhuman")
D <- Visium_B$covariates
D$response <- Visium_B$expression[[genes[1L]]]
G <- gam(
  response ~ offset(offset0) +
    s(x, y, bs = "spde", xt = list(mesh = Visium_B$meshes$spde),
      sp = c(-1, kappa)),
  data = D, family = nb(link = "log"), method = "REML", fit = FALSE,
  control = gam.control(nthreads = 1L, ncv.threads = 1L)
)
Y <- t(as.matrix(Visium_B$expression[, genes, drop = FALSE]))

BPPARAM <- SnowParam(
  workers = workers, type = "SOCK", tasks = 0L, stop.on.error = TRUE,
  progressbar = FALSE
)
BPPARAM <- bpstart(BPPARAM)
fitmgcvST <- mgcvST.estimate(Y, G, feature_id = genes, BPPARAM = BPPARAM)
if (nrow(fitmgcvST$failures)) stop("At least one Visium_B marginal fit failed.")

pairs <- t(combn(genes, 2L))
testmgcvST <- mgcvST.test(
  fitmgcvST, pairs = pairs, BPPARAM = BPPARAM,
  calibration = "liu", threads = workers
)
BPPARAM <- bpstop(BPPARAM)

write.csv(fitmgcvST$diagnostics,
          file.path(out.dir, "Visium_B_mgcvST_estimate.csv"), row.names = FALSE)
write.csv(testmgcvST$results,
          file.path(out.dir, "Visium_B_mgcvST_test.csv"), row.names = FALSE)
saveRDS(fitmgcvST, file.path(out.dir, "Visium_B_mgcvST_estimate.rds"))
saveRDS(testmgcvST, file.path(out.dir, "Visium_B_mgcvST_test.rds"))
