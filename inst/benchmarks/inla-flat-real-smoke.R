# Real-data execution check; these observed pairs are not null simulations.
Sys.setenv(OMP_NUM_THREADS="1",OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="1")
suppressPackageStartupMessages(library(mgcvST))
input <- readRDS("artifacts/pathwaylgm-151673/benchmark-input.rds")
chosen <- order(rowMeans(input$Y))[1:2]
Y <- input$Y[chosen,,drop=FALSE]
basis <- spde_basis(input$mesh,as.matrix(input$data[c("x","y")]),kappa=.7,project_intercept=TRUE)
model <- inlaST.set(response~offset(offset0),input$data,basis,family=mgcv::nb())
flat <- list(prior="flat",param=numeric(),initial=0)
historical_normal <- list(prior="normal",param=c(0,1/9),initial=0)
out <- "artifacts/flat-prior-investigation/real-slice"
dir.create(out,recursive=TRUE,showWarnings=FALSE)
rows <- lapply(c("flat_spatial","flat_both"),function(variant) {
  ctl <- list(precision_prior=flat,nb_size_prior=historical_normal)
  if(variant=="flat_both") ctl$nb_size_prior <- flat
  fit <- inlaST.estimate(Y,model,diagnostics=TRUE,BPPARAM=BiocParallel::SerialParam(),control=ctl)
  score <- mgcvST.test(fit,pairs=matrix(rownames(Y),nrow=1),calibration="davies",BPPARAM=BiocParallel::SerialParam())
  saveRDS(list(fit=fit,score=score,raw=model$inla_spec$random[[1]],Y=Y),file.path(out,paste0(variant,".rds")))
  data.frame(variant=variant,gene=rownames(Y),mean_count=rowMeans(Y),
    converged=fit$diagnostics$converged,mean_error=apply(abs(fit$observation_spatial_mean),1,max),
    tau=as.numeric(fit$lambda/fit$dispersion),
    nb_size=vapply(fit$family_parameters,function(x)x[1],numeric(1)),
    pair_p=score$results$p_two_sided,
    fit_error=fit$diagnostics$error_message)
})
rows <- do.call(rbind,rows)
write.csv(rows,file.path(out,"summary.csv"),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(out,"session-info.txt"))
print(rows,row.names=FALSE)
