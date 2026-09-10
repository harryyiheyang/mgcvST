Sys.setenv(OMP_NUM_THREADS="1",OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="1")
suppressPackageStartupMessages({library(mgcvST);library(mgcv)})
out <- Sys.getenv("MGCVST_SWAP_OUT","artifacts/lowcount-investigation/swaps")
reps <- as.integer(Sys.getenv("MGCVST_SWAP_REPS","100"))
workers <- as.integer(Sys.getenv("MGCVST_SWAP_WORKERS","3"))
dir.create(out,recursive=TRUE,showWarnings=FALSE)
cache <- file.path(out,"cache"); dir.create(cache,showWarnings=FALSE)
frozen <- file.path(out,"frozen");dir.create(frozen,showWarnings=FALSE)
source_files <- c("inla-lowcount-swap-worker.R","inla-raw-kernel-sparse.R","inla-constraint-operators.R")
for(f in source_files) file.copy(file.path("inst/benchmarks",f),file.path(frozen,f),overwrite=TRUE)
frozen <- normalizePath(frozen,winslash="/")
dgp <- readRDS("artifacts/constraint-type1/estimated/nb03_pair_k6-dgp.rds")
d <- dgp$data; d$response <- 0
basis <- spde_basis(dgp$mesh,as.matrix(d[c("x","y")]),kappa=6,project_intercept=TRUE)
basis$component <- basis$score.component <- "global"
G <- mgcv::bam(response~offset(offset0)+s(x,y,bs="spde",xt=basis),data=d,
              family=mgcv::nb(),method="fREML",discrete=TRUE,fit=FALSE,nthreads=1L)
model <- inlaST.set(response~offset(offset0),d,basis,family=mgcv::nb())
raw <- model$inla_spec$random[[1]]
payload <- list(d=d,truth=dgp$factor_truth,offset=d$offset0,G=G,
 response_index=attr(G$terms,"response"),family_raw=serialize(G$family,NULL),
 spec=model$inla_spec,A=raw$A,Q=raw$Q,g=raw$constraint,Z=raw$projection,
 X=model$inla_spec$fixed$X,cache=normalizePath(cache,winslash="/"),frozen=frozen)
writeLines(c(sprintf("Exploratory matched component swaps; prespecified first %d replicates, seeds 61000+i.",reps),
 "All fitted fields mean zero; INLA positive hyperparameters log N(0,9); no PC priors.",
 "Eta/tau/size swaps are diagnostic synthetic states, not validated alternative estimators.",
 "Every swapped state reconstructs its coherent expected-Fisher nuisance covariance.",
 "bam_original_Vp preserves original package nuisance block as a separate comparison."),file.path(out,"protocol.txt"))
entry <- function(indices,payload) {
 e <- new.env(parent=asNamespace("mgcvST"))
 for(f in c("inla-constraint-operators.R","inla-raw-kernel-sparse.R","inla-lowcount-swap-worker.R"))
   sys.source(file.path(payload$frozen,f),envir=e)
 e$lowcount_swap_block(indices,payload)
}
pool <- if(workers>1L) BiocParallel::SnowParam(workers=workers,type="SOCK") else BiocParallel::SerialParam()
ans <- BiocParallel::bplapply(split(seq_len(reps),ceiling(seq_len(reps)/10)),entry,payload=payload,BPPARAM=pool)
ans <- unlist(ans,recursive=FALSE)
rows <- do.call(rbind,lapply(ans,`[[`,"rows"))
metadata <- do.call(rbind,lapply(ans,`[[`,"metadata"))
failures <- do.call(rbind,lapply(ans,`[[`,"failure"))
write.csv(rows,file.path(out,"replicates.csv"),row.names=FALSE)
write.csv(metadata,file.path(out,"parameters.csv"),row.names=FALSE)
if(!is.null(failures)) write.csv(failures,file.path(out,"failures.csv"),row.names=FALSE)
if(!is.null(rows)) {
 summary <- do.call(rbind,lapply(split(rows,rows$variant),function(z) {
   valid <- is.finite(z$p_value); k <- sum(z$p_value[valid]<.05); nv <- sum(valid)
   ci <- binom.test(k,nv)$conf.int
   data.frame(variant=z$variant[1],attempted=reps,valid=nv,rejected=k,rate=k/nv,
     ci_lower=ci[1],ci_upper=ci[2],fallback=sum(z$fallback,na.rm=TRUE),
     max_nuisance_leakage=max(z$nuisance_leakage))
 }))
 write.csv(summary,file.path(out,"summary.csv"),row.names=FALSE);print(summary,row.names=FALSE)
}
if(!is.null(failures)) print(failures)
