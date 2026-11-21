library(BiocParallel)
library(mgcv)
library(spdeST)

workers <- 3L
kappa <- 0.1
out.dir <- Sys.getenv("SPDEST_EXAMPLE_OUTPUT", ".")
if (!dir.exists(out.dir)) dir.create(out.dir, recursive = TRUE)

data(Visium_B, package = "spdeST")
genes <- c("mt_co3", "mt_co2", "BRAFhuman")
D <- Visium_B$covariates
basis <- spde_basis(
  Visium_B$meshes$spde, as.matrix(D[, c("x", "y")]), kappa = kappa
)
model <- spdeST.set(
  response ~ offset(offset0) + s(x, y, bs = "spde", xt = basis),
  D, family = nb(link = "log")
)
Y <- t(as.matrix(Visium_B$expression[, genes, drop = FALSE]))

BPPARAM <- SnowParam(
  workers = workers, type = "SOCK", tasks = 0L, stop.on.error = TRUE,
  progressbar = FALSE
)
BPPARAM <- bpstart(BPPARAM)
scan <- spdeST.scan(Y, model, feature_id = genes, BPPARAM = BPPARAM)
fit <- spdeST.estimate(scan, Y, features = "all", BPPARAM = BPPARAM)
diagnostics <- fit$diagnostics
failed <- with(diagnostics,
  is.na(converged) | !converged |
    !is.na(error_class) | !is.na(error_message)
)
if (any(failed)) {
  stop("At least one Visium_B spatial fit failed: ",
       paste(diagnostics$feature_id[failed], collapse = ", "))
}

pairs <- t(combn(genes, 2L))
test <- spdeST.test(fit, pairs = pairs, threads = workers, moments = "exact")
BPPARAM <- bpstop(BPPARAM)

write.csv(fit$diagnostics,
          file.path(out.dir, "Visium_B_spdeST_estimate.csv"), row.names = FALSE)
results <- cbind(
  feature1 = test$feature_id[test$results$i],
  feature2 = test$feature_id[test$results$j],
  test$results
)
write.csv(results, file.path(out.dir, "Visium_B_spdeST_test.csv"), row.names = FALSE)
saveRDS(fit, file.path(out.dir, "Visium_B_spdeST_estimate.rds"))
saveRDS(test, file.path(out.dir, "Visium_B_spdeST_test.rds"))
