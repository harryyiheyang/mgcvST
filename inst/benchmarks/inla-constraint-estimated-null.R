# Estimated-hyperparameter type-I study; every actual INLA fit remains centered.
Sys.setenv(OMP_NUM_THREADS="1",OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="1")
suppressPackageStartupMessages(library(mgcvST))
reps <- as.integer(Sys.getenv("MGCVST_MC_REPS","500"))
workers <- as.integer(Sys.getenv("MGCVST_MC_WORKERS","6"))
out <- Sys.getenv("MGCVST_MC_OUTPUT","artifacts/constraint-type1/estimated")
stopifnot(reps>=2L,workers>=1L)
dir.create(out,recursive=TRUE,showWarnings=FALSE)
variants <- c("projected","raw_constrained","raw_kernel_only",
              "raw_full_keep_nuisance","raw_full_recompute_nuisance")
cases <- data.frame(
  name=c("gaussian_pair_k6","nb3_pair_k6","nb03_pair_k6","nb3_pair_k07",
         "gaussian_marginal","nb3_marginal","nb03_marginal"),
  test=c(rep("pair",4),rep("marginal",3)),
  family=c("gaussian",rep("negative_binomial",3),"gaussian",rep("negative_binomial",2)),
  mean_count=c(NA,3,.3,3,NA,3,.3), kappa=c(6,6,6,.7,.7,.7,.7),
  phi=.25,nb_size=2,seed=31000L+seq_len(7L)*10000L,
  stringsAsFactors=FALSE)
selected <- Sys.getenv("MGCVST_MC_CASES","")
if (nzchar(selected)) cases <- cases[cases$name %in% strsplit(selected,",",fixed=TRUE)[[1]],,drop=FALSE]
stopifnot(nrow(cases)>0)
write.csv(cases,file.path(out,"design.csv"),row.names=FALSE)
writeLines(c(
  sprintf("Pre-specified replicates per case: %d; INLA workers: %d, one thread each",reps,workers),
  "n=200 irregular observations, 6x6 mesh; 70% Beta(4,1) x Beta(2,5) cluster, 30% uniform.",
  "Pair null: two independent, nonzero observation-mean-zero Gaussian SPDE fields (rho=0).",
  "Each field has mean pointwise variance 0.36; Gaussian noise variance 0.25; NB true size=2.",
  "NB first-feature unconditional mean across points is exactly mean_count; second feature is mean_count*exp(.25).",
  "Marginal null: one feature, spatial field exactly zero, same known offset.",
  "All latent precisions and observation parameters are estimated anew by INLA for each feature.",
  "Every fitted spatial field has mean zero; log-positive-parameter priors are N(0,3^2), no PC priors.",
  "Five pair test constructions and two marginal kernels are applied to identical fits per replicate.",
  "Pair working model and calibration use existing package conventions. Marginal uses null D with intercept projection.",
  "Liu and Davies (with explicitly recorded Liu fallback) are reported separately.",
  "Nominal levels .05, .01, .001; exact binomial CIs; failures retained, not silently removed.",
  "The all-attempts rate treating failures as nonrejects is only a lower bound, not a validated type-I estimate.",
  "Extended-family marginal masks weights <=1e-12 exactly as the package; minimum weights are recorded.",
  "Production/reference pair operator equivalence was validated in artifacts/constraint-type1/operator-check.csv.",
  "This is a finite-sample simulation at specified meshes/count regimes, not a guarantee for arbitrary data."
),file.path(out,"protocol.txt"))

set.seed(203909)
n <- 200L
cluster <- seq_len(n)<=140L
d <- data.frame(x=runif(n,.01,.99),y=runif(n,.01,.99))
d$x[cluster] <- .01+.98*rbeta(sum(cluster),4,1)
d$y[cluster] <- .01+.98*rbeta(sum(cluster),2,5)
d$offset0 <- .5*(d$x-d$y)
vertices <- as.matrix(expand.grid(x=seq(0,1,length.out=6),y=seq(0,1,length.out=6)))
mesh <- list(loc=vertices,graph=list(tv=geometry::delaunayn(vertices)))
source_dir <- normalizePath("inst/benchmarks",winslash="/",mustWork=TRUE)
frozen <- file.path(out,"frozen")
dir.create(frozen,showWarnings=FALSE)
for (f in c("inla-constraint-estimated-worker.R","inla-constraint-operators.R")) {
  file.copy(file.path(source_dir,f),file.path(frozen,f),overwrite=TRUE)
}
worker_file <- normalizePath(file.path(frozen,"inla-constraint-estimated-worker.R"),winslash="/")
operator_file <- normalizePath(file.path(frozen,"inla-constraint-operators.R"),winslash="/")
pool <- if (workers==1L) BiocParallel::SerialParam() else
  BiocParallel::SnowParam(workers=workers,type="SOCK",stop.on.error=TRUE)
BiocParallel::bpstart(pool)

worker_entry <- function(indices,payload) {
  Sys.setenv(OMP_NUM_THREADS="1",OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="1")
  e <- new.env(parent=asNamespace("mgcvST"))
  sys.source(payload$operator_file,envir=e)
  sys.source(payload$worker_file,envir=e)
  e$inlast_constraint_mc_block(indices,payload)
}
summarize_rows <- function(rows) {
  groups <- split(rows,interaction(rows$case,rows$variant,rows$calibration,drop=TRUE))
  do.call(rbind,lapply(groups,function(g) {
    valid <- is.finite(g$p_value) & g$p_value>=0 & g$p_value<=1
    do.call(rbind,lapply(c(.05,.01,.001),function(alpha) {
      nv <- sum(valid); nr <- sum(g$p_value[valid]<alpha)
      ci <- if(nv) binom.test(nr,nv)$conf.int else c(NA_real_,NA_real_)
      data.frame(case=g$case[1],variant=g$variant[1],calibration=g$calibration[1],
                 alpha=alpha,attempted=nrow(g),valid=nv,failed=nrow(g)-nv,
                 fallback=sum(g$fallback,na.rm=TRUE),rejected=nr,
                 rejection_rate=if(nv) nr/nv else NA_real_,
                 ci_lower=ci[1],ci_upper=ci[2],
                 rejection_rate_lower_bound_all_attempts=nr/nrow(g))
    }))
  }))
}
all_rows <- list()
for (k in seq_len(nrow(cases))) {
  case <- as.list(cases[k,,drop=FALSE])
  cat(sprintf("Starting %s: %d replicates\n",case$name,reps)); flush.console()
  basis <- spde_basis(mesh,as.matrix(d[c("x","y")]),kappa=case$kappa,project_intercept=TRUE)
  family <- if(case$family=="gaussian") gaussian() else mgcv::nb()
  model <- inlaST.set(response~offset(offset0),d,basis,family=family)
  raw <- model$inla_spec$random[[1]]
  Q <- as.matrix(raw$Q); A <- as.matrix(raw$A); Z <- raw$projection
  factor_raw <- A %*% backsolve(chol(Q),diag(ncol(Q)))
  factor_centered <- (A %*% Z) %*% backsolve(chol(crossprod(Z,Q%*%Z)),diag(ncol(Z)))
  tau_truth <- mean(rowSums(factor_centered^2))/.36
  payload <- list(case=case,spec=model$inla_spec,A=raw$A,Q=raw$Q,g=raw$constraint,
                  Z=Z,X=model$inla_spec$fixed$X,offset=model$offset,
                  factor_truth=factor_centered/sqrt(tau_truth),
                  factor_centered=factor_centered,factor_raw=factor_raw,
                  variants=variants,worker_file=worker_file,operator_file=operator_file)
  saveRDS(list(case=case,data=d,mesh=mesh,tau_truth=tau_truth,
               factor_truth=payload$factor_truth),file.path(out,paste0(case$name,"-dgp.rds")))
  chunks <- split(seq_len(reps),ceiling(seq_len(reps)/25L))
  evaluated <- BiocParallel::bplapply(chunks,worker_entry,payload=payload,BPPARAM=pool)
  rows <- do.call(rbind,evaluated)
  write.csv(rows,file.path(out,paste0(case$name,"-replicates.csv")),row.names=FALSE)
  all_rows[[case$name]] <- rows
  combined <- do.call(rbind,all_rows)
  write.csv(summarize_rows(combined),file.path(out,"type1-summary.csv"),row.names=FALSE)
  cat(sprintf("Finished %s: %d valid / %d test results\n",case$name,sum(is.finite(rows$p_value)),nrow(rows)))
  flush.console()
}
BiocParallel::bpstop(pool)
writeLines(capture.output(sessionInfo()),file.path(out,"session-info.txt"))
print(subset(summarize_rows(do.call(rbind,all_rows)),alpha==.05),row.names=FALSE)
