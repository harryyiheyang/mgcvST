#!/usr/bin/env Rscript
options(stringsAsFactors=FALSE)
Sys.setenv(OMP_NUM_THREADS="1",OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="1",
           RCPP_PARALLEL_NUM_THREADS="1")
suppressPackageStartupMessages({library(mgcvST);library(mgcv)})

out<-Sys.getenv("MGCVST_SCORE_PAIR_OUTPUT",
  "artifacts/flat-prior-investigation/large-benchmark/sparse-score-paired-n32000-m225")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
n<-32000L;mesh_side<-15L;repeats<-3L
timed<-function(expr){gc(FALSE);t<-proc.time()[["elapsed"]];v<-force(expr);
  list(value=v,seconds=proc.time()[["elapsed"]]-t)}
set.seed(121015L)
d<-data.frame(x=runif(n,.01,.99),y=runif(n,.01,.99))
d$offset0<-.35*(d$x-d$y)
vertices<-as.matrix(expand.grid(x=seq(0,1,length.out=mesh_side),
                                y=seq(0,1,length.out=mesh_side)))
mesh<-list(loc=vertices,graph=list(tv=geometry::delaunayn(vertices)))
basis<-mgcvST:::spde_basis(mesh,as.matrix(d[c("x","y")]),kappa=6,project_intercept=TRUE)
basis$component<-basis$score.component<-"global"
F<-basis$B%*%backsolve(chol(basis$Q),diag(ncol(basis$Q)))
tau<-mean(rowSums(F^2))/.36
signal<-F%*%matrix(rnorm(ncol(F)*2L),ncol=2L)/sqrt(tau)
variance<-rowSums(F^2)/tau
beta<-log(.3)-log(mean(exp(d$offset0+.5*variance)))
Y<-t(vapply(1:2,function(j)rnbinom(n,mu=exp(beta+.15*(j-1)+d$offset0+signal[,j]),size=2),numeric(n)))
rownames(Y)<-c("feature1","feature2")
setting<-inlaST.set(response~offset(offset0),d,basis,family=mgcv::nb())
fit_time<-timed(inlaST.estimate(Y,setting,retain_smooth=FALSE,
  BPPARAM=BiocParallel::SerialParam(),score_backend="auto",
  control=list(precision_prior=list(prior="flat",param=numeric(),initial=0),
               nb_size_prior=list(prior="flat",param=numeric(),initial=0)),
  marginal_args=list(method="liu")))
fit<-fit_time$value
rows<-list();k<-0L
for(iteration in seq_len(repeats)) for(backend in c("sparse","dense")) {
  candidate<-fit
  candidate$score_backend<-backend
  candidate$.mgcvst_fixed_factors<-NULL
  candidate$.mgcvst_state_cache<-NULL
  if(backend=="dense") candidate$score_sparse<-NULL
  ans<-timed(mgcvST.test(candidate,pairs=matrix(c("feature1","feature2"),1),
    calibration="liu",BPPARAM=BiocParallel::SerialParam()))
  k<-k+1L
  rows[[k]]<-data.frame(iteration=iteration,warm=iteration>1L,backend=backend,
    pair_score_seconds=ans$seconds,p_two_sided=ans$value$results$p_two_sided[1],
    statistic=ans$value$results$statistic[1],valid=is.finite(ans$value$results$p_two_sided[1]))
  rm(candidate,ans);gc(FALSE)
}
result<-do.call(rbind,rows)
write.csv(result,file.path(out,"paired-score.csv"),row.names=FALSE)
equiv<-merge(result[result$backend=="sparse",],result[result$backend=="dense",],by="iteration")
equiv$p_abs_difference<-abs(equiv$p_two_sided.x-equiv$p_two_sided.y)
equiv$statistic_abs_difference<-abs(equiv$statistic.x-equiv$statistic.y)
write.csv(equiv,file.path(out,"equivalence.csv"),row.names=FALSE)
writeLines(c(sprintf("mgcvST library: %s",find.package("mgcvST")),
 sprintf("mgcvST version: %s",packageVersion("mgcvST")),
 sprintf("single shared fit seconds: %.6f",fit_time$seconds),
 "n=32000, mesh=225, two NB features; single-threaded fit and score.",
 "Flat precision and NB-size priors; observation spatial mean-zero constraint active.",
 "Each backend starts with fixed-factor and state caches cleared; dense also removes score_sparse.",
 "Runs occurred during the root four-worker Monte Carlo, so use within-fit paired ratios."
),file.path(out,"protocol.txt"))
writeLines(capture.output(sessionInfo()),file.path(out,"session-info.txt"))
print(result);print(equiv)
