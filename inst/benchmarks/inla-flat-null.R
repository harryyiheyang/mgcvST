# Paired flat-log-hyperprior probe; retains every attempted pair and fit failure.
Sys.setenv(OMP_NUM_THREADS="1",OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="1")
suppressPackageStartupMessages(library(mgcvST))
out <- Sys.getenv("MGCVST_FLAT_OUT","artifacts/flat-prior-investigation/pilot")
reps <- as.integer(Sys.getenv("MGCVST_FLAT_REPS","10"))
workers <- as.integer(Sys.getenv("MGCVST_FLAT_WORKERS","2"))
seed_base <- as.integer(Sys.getenv("MGCVST_FLAT_SEED","161000"))
mean_count <- as.numeric(Sys.getenv("MGCVST_FLAT_MEAN",".3"))
rho <- as.numeric(Sys.getenv("MGCVST_FLAT_RHO","0"))
variants <- strsplit(Sys.getenv("MGCVST_FLAT_VARIANTS","flat_spatial,flat_both,flat_both_scaled"),",",fixed=TRUE)[[1]]
stopifnot(reps>0,workers>0,all(variants%in%c("flat_spatial","flat_both","flat_both_scaled")))
dir.create(out,recursive=TRUE,showWarnings=FALSE)
cache <- file.path(out,"cache"); dir.create(cache,showWarnings=FALSE)
results <- file.path(out,"completed"); dir.create(results,showWarnings=FALSE)
frozen <- file.path(out,"frozen"); dir.create(frozen,showWarnings=FALSE)
for(f in c("inla-flat-null-worker.R","inla-lowcount-candidate-worker.R"))
  stopifnot(file.copy(file.path("inst/benchmarks",f),file.path(frozen,f),overwrite=TRUE))
ns <- asNamespace("mgcvST")
engine_names <- ls(ns,all.names=TRUE)
dump(engine_names[grepl("^\\.inlast_",engine_names)],envir=ns,
     file=file.path(frozen,"installed-inla-functions.R"))
dgp <- readRDS("artifacts/constraint-type1/estimated/nb03_pair_k6-dgp.rds")
d <- dgp$data
basis <- spde_basis(dgp$mesh,as.matrix(d[c("x","y")]),kappa=6,project_intercept=TRUE)
model <- inlaST.set(response~offset(offset0),d,basis,family=mgcv::nb())
raw <- model$inla_spec$random[[1]]; A<-as.matrix(raw$A);Q<-as.matrix(raw$Q);Z<-raw$projection
Fc <- (A%*%Z)%*%backsolve(chol(crossprod(Z,Q%*%Z)),diag(ncol(Z)))
Fr <- A%*%backsolve(chol(Q),diag(ncol(Q)))
payload <- list(spec=model$inla_spec,X=model$inla_spec$fixed$X,offset=model$offset,
  truth=dgp$factor_truth,Fconditioned=Fc,Fraw=Fr,Fraw_centered=sweep(Fr,2,colMeans(Fr)),
  Q_scale=mean(rowSums(Fc^2)),seed_base=seed_base,mean_count=mean_count,rho=rho,
  variants=variants,kernels=c("conditioned","raw_centered"),
  results=normalizePath(results,winslash="/"),
  cache=normalizePath(cache,winslash="/"),frozen=normalizePath(frozen,winslash="/"),
  old_cache=normalizePath("artifacts/lowcount-investigation/validation-500/cache",winslash="/"))
saveRDS(payload,file.path(out,"payload.rds"))
writeLines(c(sprintf("Prespecified %d pairs, seed base %d, L'Ecuyer-CMRG, mean %g, rho %g.",reps,seed_base,mean_count,rho),
  "All fits retain observation mean=0, native INLA conditional Gaussian Vp, expected working states; no PC prior.",
  "flat_spatial: flat log spatial precision, normal(0,9) log NB size.",
  "flat_both: flat log spatial precision and flat log NB size.",
  "flat_both_scaled: same flat objective, Q scaled to unit average observation variance.",
  "Finite returned hyperparameter values or mode status 0 do not certify interior maxima or proper hyperposteriors.",
  "All invalid scores retained; no boundary clipping, p-value substitution, or renormalization.",
  paste("Variants:",paste(variants,collapse=","))),file.path(out,"protocol.txt"))
entry <- function(indices,payload) {
  env <- new.env(parent=asNamespace("mgcvST"))
  for(f in c("inla-lowcount-candidate-worker.R","inla-flat-null-worker.R"))
    sys.source(file.path(payload$frozen,f),envir=env)
  env$flat_null_block(indices,payload)
}
bp <- if(workers>1) BiocParallel::SnowParam(workers,type="SOCK") else BiocParallel::SerialParam()
ans <- tryCatch(BiocParallel::bplapply(split(seq_len(reps),ceiling(seq_len(reps)/5)),entry,
                   payload=payload,BPPARAM=bp),finally=BiocParallel::bpstop(bp))
ans <- unlist(ans,recursive=FALSE)
rows <- do.call(rbind,lapply(ans,`[[`,"rows")); meta<-do.call(rbind,lapply(ans,`[[`,"metadata"))
write.csv(rows,file.path(out,"replicates.csv"),row.names=FALSE)
write.csv(meta,file.path(out,"parameters.csv"),row.names=FALSE)
summary <- do.call(rbind,lapply(split(rows,interaction(rows$variant,rows$kernel,drop=TRUE)),function(z) {
  ok<-is.finite(z$p_value); nv<-sum(ok); k<-sum(z$p_value[ok]<.05)
  ci<-if(nv)binom.test(k,nv)$conf.int else c(NA,NA)
  data.frame(variant=z$variant[1],kernel=z$kernel[1],attempted=reps,valid=nv,invalid=reps-nv,
    rejected=k,rate=if(nv)k/nv else NA_real_,ci_lower=ci[1],ci_upper=ci[2],
    all_attempt_lower=k/reps,all_attempt_upper=(k+reps-nv)/reps,
    fallback=sum(z$fallback,na.rm=TRUE),mean_pair_fit_seconds=mean(z$elapsed,na.rm=TRUE))
}))
write.csv(summary,file.path(out,"summary.csv"),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(out,"session-info.txt"))
print(summary,row.names=FALSE)
