Sys.setenv(OMP_NUM_THREADS="1",OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="1")
suppressPackageStartupMessages(library(mgcvST))
out <- Sys.getenv("MGCVST_CANDIDATE_OUT","artifacts/lowcount-investigation/candidates")
reps <- as.integer(Sys.getenv("MGCVST_CANDIDATE_REPS","100"))
workers <- as.integer(Sys.getenv("MGCVST_CANDIDATE_WORKERS","3"))
seed_base <- as.integer(Sys.getenv("MGCVST_CANDIDATE_SEED","61000"))
mean_count <- as.numeric(Sys.getenv("MGCVST_CANDIDATE_MEAN","0.3"))
rho <- as.numeric(Sys.getenv("MGCVST_CANDIDATE_RHO","0"))
fixed_truth <- identical(Sys.getenv("MGCVST_CANDIDATE_FIXED_TRUTH","false"),"true")
stopifnot(is.finite(rho),abs(rho)<1)
scale_names<-strsplit(Sys.getenv("MGCVST_CANDIDATE_SCALES","original,unit_mean_variance"),",",fixed=TRUE)[[1]]
working_names<-strsplit(Sys.getenv("MGCVST_CANDIDATE_WORKING","legacy_expected,posterior_expected,posterior_observed"),",",fixed=TRUE)[[1]]
kernel_names<-strsplit(Sys.getenv("MGCVST_CANDIDATE_KERNELS","conditioned,raw,raw_centered"),",",fixed=TRUE)[[1]]
dir.create(out,recursive=TRUE,showWarnings=FALSE)
cache<-file.path(out,"cache");dir.create(cache,showWarnings=FALSE)
frozen<-file.path(out,"frozen");dir.create(frozen,showWarnings=FALSE)
for(f in c("inla-lowcount-candidate-worker.R","inla-posterior-vp.R"))
  stopifnot(file.copy(file.path("inst/benchmarks",f),file.path(frozen,f),overwrite=TRUE))
dgp<-readRDS("artifacts/constraint-type1/estimated/nb03_pair_k6-dgp.rds")
d<-dgp$data
basis<-spde_basis(dgp$mesh,as.matrix(d[c("x","y")]),kappa=6,project_intercept=TRUE)
model<-inlaST.set(response~offset(offset0),d,basis,family=mgcv::nb())
raw<-model$inla_spec$random[[1]];A<-as.matrix(raw$A);Q<-as.matrix(raw$Q);Z<-raw$projection
Fc<-(A%*%Z)%*%backsolve(chol(crossprod(Z,Q%*%Z)),diag(ncol(Z)))
Fr<-A%*%backsolve(chol(Q),diag(ncol(Q)))
Q_scale<-mean(rowSums(Fc^2))
payload<-list(spec=model$inla_spec,X=model$inla_spec$fixed$X,offset=model$offset,
  truth=dgp$factor_truth,Fconditioned=Fc,Fraw=Fr,Fraw_centered=sweep(Fr,2,colMeans(Fr)),
  Q_scale=Q_scale,seed_base=seed_base,mean_count=mean_count,
  rho=rho,fixed_truth=fixed_truth,tau_truth=dgp$tau_truth,
  scale_names=scale_names,working_names=working_names,kernel_names=kernel_names,
  cache=normalizePath(cache,winslash="/"),frozen=normalizePath(frozen,winslash="/"))
writeLines(c(sprintf("Prespecified %d replicates, RNGkind L'Ecuyer-CMRG, seed base %d; NB mean %g, size 2.",reps,seed_base,mean_count),
  "Every INLA fit has mean-zero constraint and log positive hyperparameter N(0,9); no PC prior.",
  sprintf("Latent cross-correlation rho=%g; fixed true hyperparameters=%s.",rho,fixed_truth),
  sprintf("Unit mean observation variance Q scale c=%.12g; Q_new=c Q, tau_original=c tau_new.",Q_scale),
  "Scaling changes the physical prior by declaring the same N(0,9) prior on standardized precision; it is not the same prior on original tau.",
  "Posterior Vp is extracted directly from the actual INLA conditional posterior precision with exact constraints.",
  "Kernel conditioned, raw, and raw_centered are reported separately: posterior Vp plus expected D need not annihilate the intercept.",
  paste("scales",paste(scale_names,collapse=",")),paste("working",paste(working_names,collapse=",")),
  paste("kernels",paste(kernel_names,collapse=","))),file.path(out,"protocol.txt"))
entry<-function(indices,payload) {
  e<-new.env(parent=asNamespace("mgcvST"))
  for(f in c("inla-posterior-vp.R","inla-lowcount-candidate-worker.R")) sys.source(file.path(payload$frozen,f),envir=e)
  e$lowcount_candidate_block(indices,payload)
}
pool<-if(workers>1) BiocParallel::SnowParam(workers=workers,type="SOCK") else BiocParallel::SerialParam()
ans<-BiocParallel::bplapply(split(seq_len(reps),ceiling(seq_len(reps)/10)),entry,payload=payload,BPPARAM=pool)
ans<-unlist(ans,recursive=FALSE)
rows<-do.call(rbind,lapply(ans,`[[`,"rows"));meta<-do.call(rbind,lapply(ans,`[[`,"metadata"))
fail<-do.call(rbind,lapply(ans,`[[`,"failure"))
write.csv(rows,file.path(out,"replicates.csv"),row.names=FALSE)
write.csv(meta,file.path(out,"parameters.csv"),row.names=FALSE)
if(!is.null(fail)) {write.csv(fail,file.path(out,"fit-failures.csv"),row.names=FALSE);print(fail)}
if(!is.null(rows)) {
  summary<-do.call(rbind,lapply(split(rows,interaction(rows$scale,rows$working,rows$kernel,drop=TRUE)),function(z) {
    valid<-is.finite(z$p_value);nv<-sum(valid);k<-sum(z$p_value[valid]<.05)
    ci<-if(nv) binom.test(k,nv)$conf.int else c(NA,NA)
    data.frame(scale=z$scale[1],working=z$working[1],kernel=z$kernel[1],attempted=reps,
      valid=nv,failed=reps-nv,rejected=k,rate=if(nv)k/nv else NA_real_,ci_lower=ci[1],ci_upper=ci[2],
      fallback=sum(z$fallback,na.rm=TRUE),maximum_P1=max(c(NA,z$P1),na.rm=TRUE))
  }))
  write.csv(summary,file.path(out,"summary.csv"),row.names=FALSE);print(summary,row.names=FALSE)
}
