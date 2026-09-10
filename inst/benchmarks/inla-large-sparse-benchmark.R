#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE)
Sys.setenv(OMP_NUM_THREADS="1",OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="1",
           RCPP_PARALLEL_NUM_THREADS="1")
suppressPackageStartupMessages({library(mgcvST);library(mgcv)})

sizes <- as.integer(strsplit(Sys.getenv("MGCVST_LARGE_N","2000,8000,32000"),",",fixed=TRUE)[[1]])
mesh_sides <- as.integer(strsplit(Sys.getenv("MGCVST_LARGE_MESH_SIDE","15"),",",fixed=TRUE)[[1]])
repeats <- as.integer(Sys.getenv("MGCVST_LARGE_REPEATS","3"))
backends <- strsplit(Sys.getenv("MGCVST_LARGE_BACKENDS","inla,bam"),",",fixed=TRUE)[[1]]
prior_label <- Sys.getenv("MGCVST_LARGE_PRIOR_LABEL","normal_log_N0_9")
out <- Sys.getenv("MGCVST_LARGE_OUTPUT","artifacts/flat-prior-investigation/large-benchmark")
stopifnot(all(sizes>=100L),all(mesh_sides>=4L),repeats>=2L,
          all(backends%in%c("inla","bam")))
dir.create(out,recursive=TRUE,showWarnings=FALSE)
flat_prior <- list(prior="flat",param=numeric(),initial=0)
historical_normal <- list(prior="normal",param=c(0,1/9),initial=0)
inla_control <- if(prior_label=="flat_both") list(
  precision_prior=flat_prior,
  nb_size_prior=flat_prior
) else if(prior_label=="flat_spatial_only") list(
  precision_prior=flat_prior,
  nb_size_prior=historical_normal
) else list(
  precision_prior=historical_normal,
  nb_size_prior=historical_normal
)

timed <- function(expr) {gc(FALSE);t<-proc.time()[["elapsed"]];v<-force(expr);
  list(value=v,seconds=proc.time()[["elapsed"]]-t)}

bam_compact <- function(fits,ids) {
  states<-lapply(fits,function(fit) {
    L<-mgcvST:::.gam_training_lpmatrix(fit)
    geometry<-mgcvST:::.mgcvst_model_geometry(fit,L)
    nuisance<-mgcvST:::.mgcvst_nuisance_state(fit,geometry,list(L=L,frozen=TRUE))
    W<-rkhs_extract_working_model(fit)
    list(W=W,geometry=geometry,nuisance=nuisance)
  })
  geometry<-states[[1]]$geometry
  geometry$nuisance_columns<-states[[1]]$nuisance$columns
  geometry$nuisance_design<-states[[1]]$nuisance$design
  geometry$nuisance_projection<-"conditional_Vp_block"
  lambda<-setNames(vapply(states,function(x) {
    j<-x$geometry$target[["global"]]
    x$geometry$sp[x$geometry$smooth[[j]]$sp_index]
  },numeric(1)),ids)
  structure(list(feature_id=ids,
    working_error=do.call(cbind,lapply(states,function(x)x$W$working_error)),
    working_variance=do.call(cbind,lapply(states,function(x)x$W$working_variance)),
    dispersion=setNames(vapply(states,function(x)x$W$dispersion,numeric(1)),ids),
    lambda=lambda,component_lambda=matrix(lambda,ncol=1,dimnames=list(ids,"global")),
    smoothing_parameters=do.call(rbind,lapply(states,function(x)x$geometry$sp)),
    nuisance_covariance=setNames(lapply(states,function(x)x$nuisance$covariance),ids),
    geometry=geometry,row_id=geometry$row_id,score_components=geometry$score_components,
    model_setting="global",test_engine="single_model",
    diagnostics=data.frame(index=1:2,feature_id=ids,converged=TRUE,error_message=NA_character_)),
    class=c("mgcvST_model_fit","mgcvST_fit","mgcvST"))
}

make_case <- function(n,mesh_side,seed) {
  set.seed(seed)
  d<-data.frame(x=runif(n,.01,.99),y=runif(n,.01,.99))
  d$offset0<-.35*(d$x-d$y)
  vertices<-as.matrix(expand.grid(x=seq(0,1,length.out=mesh_side),
                                  y=seq(0,1,length.out=mesh_side)))
  mesh<-list(loc=vertices,graph=list(tv=geometry::delaunayn(vertices)))
  basis<-spde_basis(mesh,as.matrix(d[c("x","y")]),kappa=6,project_intercept=TRUE)
  basis$component<-basis$score.component<-"global"
  F<-basis$B%*%backsolve(chol(basis$Q),diag(ncol(basis$Q)))
  tau<-mean(rowSums(F^2))/.36
  innovation<-matrix(rnorm(ncol(F)*2L),ncol=2L)
  signal<-F%*%innovation/sqrt(tau)
  variance<-rowSums(F^2)/tau
  beta<-log(.3)-log(mean(exp(d$offset0+.5*variance)))
  Y<-t(vapply(1:2,function(j)rnbinom(n,mu=exp(beta+.15*(j-1)+d$offset0+signal[,j]),size=2),numeric(n)))
  rownames(Y)<-c("feature1","feature2")
  list(d=d,mesh=mesh,basis=basis,Y=Y)
}

rows<-list();sparsity<-list();case_index<-0L
for(mesh_side in mesh_sides) for(n in sizes) {
  # Optional 900-node stress is bounded to n=8000 unless explicitly requested alone.
  if(mesh_side>=30L && n!=8000L && length(sizes)>1L) next
  case_index<-case_index+1L
  made<-timed(make_case(n,mesh_side,74000L+1000L*mesh_side+n))
  obj<-made$value;d<-obj$d;basis<-obj$basis;Y<-obj$Y
  inla_setup<-timed(inlaST.set(response~offset(offset0),d,basis,family=mgcv::nb()))
  spec<-inla_setup$value$inla_spec;A<-spec$random[[1]]$A;Q<-spec$random[[1]]$Q
  cholQ<-Matrix::Cholesky(Matrix::forceSymmetric(Q),LDL=FALSE)
  cholQ_sparse<-as(cholQ,"sparseMatrix")
  sparsity[[case_index]]<-data.frame(n=n,mesh_side=mesh_side,mesh_vertices=ncol(A),
    A_stored_entries=length(A@x),A_matrix_nnzero=Matrix::nnzero(A),
    A_density=Matrix::nnzero(A)/prod(dim(A)),
    Q_stored_triangle_entries=length(Q@x),Q_matrix_nnzero=Matrix::nnzero(Q),
    Q_density=Matrix::nnzero(Q)/prod(dim(Q)),
    chol_factor_stored_entries=length(cholQ_sparse@x),
    chol_factor_matrix_nnzero=Matrix::nnzero(cholQ_sparse),
    projected_B_bytes=as.numeric(object.size(basis$B)),
    one_dense_score_factor_bytes=8*n*(ncol(A)-1L),basis_dgp_seconds=made$seconds)

  if("inla"%in%backends) for(iteration in seq_len(repeats)) {
    fit<-tryCatch(timed(inlaST.estimate(Y,inla_setup$value,retain_smooth=FALSE,
      BPPARAM=BiocParallel::SerialParam(),control=inla_control,marginal_args=list(method="liu"))),
      error=function(e)e)
    if(inherits(fit,"condition")) {
      rows[[length(rows)+1L]]<-data.frame(n=n,mesh_side=mesh_side,iteration=iteration,
        warm=iteration>1L,backend="inla",prior=prior_label,
        setup_seconds=inla_setup$seconds,pure_fit_seconds=NA_real_,
        estimator_total_seconds=NA_real_,compact_only_seconds=NA_real_,compact_marginal_seconds=NA_real_,
        marginal_seconds=NA_real_,pair_score_seconds=NA_real_,valid_pair=FALSE,
        max_mean_error=NA_real_,peak_memory_note=paste("FIT FAILURE:",conditionMessage(fit)))
      next
    }
    engine_seconds<-sum(fit$value$diagnostics$fit_seconds,na.rm=TRUE)
    marginal_seconds<-fit$value$timing$marginal_elapsed
    pair<-timed(mgcvST.test(fit$value,pairs=matrix(c("feature1","feature2"),1),
      calibration="liu",BPPARAM=BiocParallel::SerialParam()))
    rows[[length(rows)+1L]]<-data.frame(n=n,mesh_side=mesh_side,iteration=iteration,
      warm=iteration>1L,backend="inla",prior=prior_label,
      setup_seconds=inla_setup$seconds,pure_fit_seconds=engine_seconds,
      estimator_total_seconds=fit$seconds,
      compact_only_seconds=max(0,fit$seconds-engine_seconds-marginal_seconds),
      compact_marginal_seconds=max(0,fit$seconds-engine_seconds),
      marginal_seconds=marginal_seconds,pair_score_seconds=pair$seconds,
      valid_pair=is.finite(pair$value$results$p_two_sided[1]),
      max_mean_error=max(abs(fit$value$observation_spatial_mean),na.rm=TRUE),
      peak_memory_note="R object proxies only")
    rm(fit,pair);gc(FALSE)
  }

  if("bam"%in%backends) {
    data_setup<-d;data_setup$response<-Y[1,]
    formula<-response~offset(offset0)+s(x,y,bs="spde",xt=basis)
    setup_bam<-timed(mgcv::bam(formula,data=data_setup,family=mgcv::nb(),
      method="fREML",discrete=TRUE,nthreads=1L,fit=FALSE))
    response_index<-attr(setup_bam$value$terms,"response")
    family_raw<-serialize(setup_bam$value$family,NULL)
    for(iteration in seq_len(repeats)) {
      fitted<-timed(lapply(1:2,function(j) {
        G<-setup_bam$value;G$y<-Y[j,];G$mf[[response_index]]<-Y[j,]
        G$family<-unserialize(family_raw)
        mgcv::bam(G=G,method="fREML",discrete=TRUE,nthreads=1L)
      }))
      marginal<-timed(lapply(fitted$value,function(x)
        mgcvST:::taps_score_test(x,test.component=1L,method="liu",n_threads=1L)))
      compact<-timed(bam_compact(fitted$value,rownames(Y)))
      pair<-timed(mgcvST.test(compact$value,pairs=matrix(c("feature1","feature2"),1),
        calibration="liu",BPPARAM=BiocParallel::SerialParam()))
      mean_error<-max(abs(vapply(fitted$value,function(x) {
        sm<-which(vapply(x$smooth,inherits,logical(1),"spde.smooth"))[1]
        cols<-x$smooth[[sm]]$first.para:x$smooth[[sm]]$last.para
        L<-mgcvST:::.gam_training_lpmatrix(x)
        mean(as.numeric(L[,cols,drop=FALSE]%*%coef(x)[cols]))
      },numeric(1))))
      rows[[length(rows)+1L]]<-data.frame(n=n,mesh_side=mesh_side,iteration=iteration,
        warm=iteration>1L,backend="bam_fREML_discrete",prior="fREML",
        setup_seconds=setup_bam$seconds,pure_fit_seconds=fitted$seconds,
        estimator_total_seconds=fitted$seconds+marginal$seconds+compact$seconds,
        compact_only_seconds=compact$seconds,
        compact_marginal_seconds=marginal$seconds+compact$seconds,
        marginal_seconds=marginal$seconds,pair_score_seconds=pair$seconds,
        valid_pair=is.finite(pair$value$results$p_two_sided[1]),
        max_mean_error=mean_error,peak_memory_note="R object proxies only")
      rm(fitted,marginal,compact,pair);gc(FALSE)
    }
  }
  rm(obj,basis,Y,spec,A,Q,cholQ,cholQ_sparse);gc(FALSE)
}

timings<-do.call(rbind,rows);sparse<-do.call(rbind,sparsity)
write.csv(timings,file.path(out,"timings.csv"),row.names=FALSE)
write.csv(sparse,file.path(out,"sparsity.csv"),row.names=FALSE)
warm<-timings[timings$warm,]
if(nrow(warm)) {
  metric_names<-c("pure_fit_seconds","estimator_total_seconds","compact_only_seconds",
                  "marginal_seconds","pair_score_seconds")
  warm_summary<-aggregate(warm[metric_names],
    warm[c("n","mesh_side","backend","prior")],median,na.rm=TRUE)
  write.csv(warm_summary,file.path(out,"warm-summary.csv"),row.names=FALSE)
}
writeLines(c(sprintf("mgcvST library: %s",find.package("mgcvST")),
  sprintf("mgcvST version: %s",packageVersion("mgcvST")),
  sprintf("INLA version: %s",packageVersion("INLA")),
  sprintf("INLA prior label supplied by runner: %s",prior_label),
  "Single thread within every fit; two features; NB mean 0.3 and true size 2.",
  "All spatial fields use the active observation mean-zero constraint.",
  "Warm timing is iteration>1. Setup, engine fit, compact/marginal, and pair score are separate.",
  "No n-by-n matrix is constructed. Dense projected B and score-factor byte proxies are reported."
),file.path(out,"protocol.txt"))
writeLines(capture.output(sessionInfo()),file.path(out,"session-info.txt"))
print(timings);print(sparse)
