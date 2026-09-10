# Public-API smoke check on the two lowest-mean genes of the saved gold slice.
# This is an execution check, not a type-I experiment or a speed benchmark.
Sys.setenv(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
suppressPackageStartupMessages(library(mgcvST))
input <- readRDS("artifacts/pathwaylgm-151673/benchmark-input.rds")
chosen <- order(rowMeans(input$Y))[1:2]
Y <- input$Y[chosen,,drop=FALSE]
basis <- spde_basis(input$mesh, as.matrix(input$data[c("x","y")]),
                    kappa=.7, project_intercept=TRUE)
out <- "artifacts/lowcount-investigation/production-real-smoke"
dir.create(out, recursive=TRUE, showWarnings=FALSE)
ans <- lapply(c("raw","observation"), function(scale) {
  model <- inlaST.set(response ~ offset(offset0), input$data, basis,
                      family=mgcv::nb(), precision_scale=scale)
  fit <- inlaST.estimate(Y, model, diagnostics=TRUE,
                         BPPARAM=BiocParallel::SerialParam())
  stopifnot(all(fit$diagnostics$converged),
            max(abs(fit$observation_spatial_mean)) < 1e-10,
            all(vapply(fit$nuisance_covariance,is.matrix,logical(1))),
            all(vapply(fit$nuisance_covariance,function(x) all(is.finite(x)),logical(1))))
  score <- mgcvST.test(fit, pairs=matrix(rownames(Y),nrow=1),
                       calibration="davies", BPPARAM=BiocParallel::SerialParam())
  stopifnot(is.finite(score$results$p_two_sided))
  saveRDS(list(diagnostics=fit$diagnostics,inla=fit$inla_diagnostics,
               nuisance_covariance=fit$nuisance_covariance,
               score=score$results),file.path(out,paste0(scale,".rds")))
  data.frame(scale=scale,observations=ncol(Y),mesh_vertices=basis$raw_dimension,
             gene=rownames(Y),mean_count=rowMeans(Y),
             max_mean_error=max(abs(fit$observation_spatial_mean)),
             posterior_Vp=vapply(fit$nuisance_covariance,function(x)x[1,1],numeric(1)),
             pair_p=score$results$p_two_sided)
})
ans <- do.call(rbind,ans)
write.csv(ans,file.path(out,"summary.csv"),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(out,"session-info.txt"))
print(ans,row.names=FALSE)
