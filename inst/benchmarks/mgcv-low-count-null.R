#!/usr/bin/env Rscript

Sys.setenv(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
           MKL_NUM_THREADS="1", RCPP_PARALLEL_NUM_THREADS="1")
suppressPackageStartupMessages({library(mgcvST); library(mgcv)})
source("inst/benchmarks/inla-constraint-operators.R", local=TRUE)
RNGkind("L'Ecuyer-CMRG")

reps <- as.integer(Sys.getenv("MGCVST_BAM_REPS", "500"))
workers <- as.integer(Sys.getenv("MGCVST_BAM_WORKERS", "1"))
oracle_reps <- as.integer(Sys.getenv("MGCVST_BAM_ORACLE_REPS", "5"))
out <- Sys.getenv("MGCVST_BAM_OUTPUT",
                  "artifacts/constraint-type1/mgcv-low-count")
mean_count <- as.numeric(Sys.getenv("MGCVST_BAM_MEAN","0.3"))
rho <- as.numeric(Sys.getenv("MGCVST_BAM_RHO","0"))
seed_base <- as.integer(Sys.getenv("MGCVST_BAM_SEED","61000"))
cache_dir <- Sys.getenv("MGCVST_BAM_CACHE",
                        "artifacts/lowcount-investigation/swaps/cache")
stopifnot(reps >= 2L, workers >= 1L, oracle_reps >= 0L)
stopifnot(is.finite(mean_count),mean_count>0,is.finite(rho),abs(rho)<1)
dir.create(out, recursive=TRUE, showWarnings=FALSE)

dgp <- readRDS("artifacts/constraint-type1/estimated/nb03_pair_k6-dgp.rds")
stopifnot(identical(dgp$case$name, "nb03_pair_k6"), dgp$case$seed == 61000L,
          dgp$case$mean_count == .3, dgp$case$nb_size == 2,
          dgp$case$kappa == 6)
d <- dgp$data; n <- nrow(d)
basis <- spde_basis(dgp$mesh, as.matrix(d[c("x","y")]), kappa=6,
                    project_intercept=TRUE)
basis$component <- basis$score.component <- "global"
data_setup <- d; data_setup$response <- 0
formula <- response ~ offset(offset0) + s(x, y, bs="spde", xt=basis)
G0 <- mgcv::bam(formula, data=data_setup, family=mgcv::nb(), method="fREML",
                discrete=TRUE, nthreads=1L, fit=FALSE)
response_index <- attr(G0$terms, "response")
family_raw <- serialize(G0$family, NULL)
inla_model <- inlaST.set(response ~ offset(offset0), d, basis,
                         family=mgcv::nb())
raw <- inla_model$inla_spec$random[[1L]]

compact_pair <- function(fits) {
  ids <- c("feature1","feature2")
  states <- lapply(fits, function(fit) {
    L <- mgcvST:::.gam_training_lpmatrix(fit)
    geometry <- mgcvST:::.mgcvst_model_geometry(fit, L)
    nuisance <- mgcvST:::.mgcvst_nuisance_state(fit, geometry,
                                                 list(L=L, frozen=TRUE))
    if (is.null(nuisance)) stop("bam did not return conditional nuisance Vp")
    list(W=rkhs_extract_working_model(fit), geometry=geometry,
         nuisance=nuisance)
  })
  geometry <- states[[1L]]$geometry
  geometry$nuisance_columns <- states[[1L]]$nuisance$columns
  geometry$nuisance_design <- states[[1L]]$nuisance$design
  geometry$nuisance_projection <- "conditional_Vp_block"
  lambda <- setNames(vapply(states,function(x) {
    j <- x$geometry$target[["global"]]
    x$geometry$sp[x$geometry$smooth[[j]]$sp_index]
  },numeric(1)),ids)
  structure(list(
    feature_id=ids,
    working_error=do.call(cbind,lapply(states,function(x)x$W$working_error)),
    working_variance=do.call(cbind,lapply(states,function(x)x$W$working_variance)),
    dispersion=setNames(vapply(states,function(x)x$W$dispersion,numeric(1)),ids),
    lambda=lambda,
    component_lambda=matrix(lambda,ncol=1L,
                            dimnames=list(ids,"global")),
    smoothing_parameters=do.call(rbind,lapply(states,function(x)x$geometry$sp)),
    nuisance_covariance=setNames(lapply(states,function(x)x$nuisance$covariance),ids),
    geometry=geometry, row_id=geometry$row_id,
    score_components=geometry$score_components,
    model_setting="global", test_engine="single_model",
    diagnostics=data.frame(index=1:2,feature_id=ids,converged=TRUE,
                           error_message=NA_character_)
  ), class=c("mgcvST_model_fit","mgcvST_fit","mgcvST"))
}

simulate_response <- function(i) {
  set.seed(seed_base + i)
  innovation <- matrix(rnorm(ncol(dgp$factor_truth)*2L), ncol=2L)
  if(rho!=0) innovation[,2] <- rho*innovation[,1]+sqrt(1-rho^2)*innovation[,2]
  signal <- dgp$factor_truth %*% innovation
  variance <- rowSums(dgp$factor_truth^2)
  beta <- log(mean_count) - log(mean(exp(d$offset0 + .5*variance)))
  y <- vapply(1:2, function(j) {
    eta <- beta + .25*(j-1L) + d$offset0 + signal[,j]
    rnbinom(n, mu=exp(eta), size=2)
  }, numeric(n))
  list(y=y, signal=signal, beta=beta)
}

failed_rows <- function(i, msg, seconds) data.frame(
  replicate=i, variant=c("conditioned_spde","raw_G_after_P_intercept"),
  calibration="davies", p_value=NA_real_, signed_score=NA_real_,
  fallback=NA, fit_seconds=seconds, mean_y1=NA_real_,mean_y2=NA_real_,
  y_checksum1=NA_real_,y_checksum2=NA_real_,response_source=NA_character_,
  theta1=NA_real_,theta2=NA_real_,sp1=NA_real_,sp2=NA_real_,
  spatial_mean_error=NA_real_,oracle_p_error=NA_real_,
  error=msg, stringsAsFactors=FALSE)

run_one <- function(i) {
  t0 <- proc.time()[["elapsed"]]
  tryCatch({
    sim <- simulate_response(i)
    cache_candidates <- file.path(cache_dir,c(sprintf("rep-%04d.rds",i),
                              sprintf("rep-%04d-unit_mean_variance.rds",i)))
    cache_file <- cache_candidates[file.exists(cache_candidates)][1L]
    response_source <- "regenerated_L'Ecuyer-CMRG"
    if (length(cache_file) && !is.na(cache_file) && file.exists(cache_file)) {
      cached <- readRDS(cache_file)$y
      cached_y <- do.call(cbind,cached)
      if (!identical(unname(sim$y),unname(cached_y)))
        stop("regenerated response does not match swap cache")
      sim$y <- cached_y
      response_source <- "verified_swap_cache"
    }
    fits <- lapply(1:2, function(j) {
      G <- G0; G$y <- sim$y[,j]; G$mf[[response_index]] <- sim$y[,j]
      G$family <- unserialize(family_raw)
      mgcv::bam(G=G, method="fREML", discrete=TRUE, nthreads=1L)
    })
    if (!all(vapply(fits,function(x)isTRUE(x$converged),logical(1))))
      stop("bam convergence failure")
    compact <- compact_pair(fits)
    prod <- mgcvST.test(compact, pairs=matrix(c("feature1","feature2"),1L),
                        calibration="davies", BPPARAM=BiocParallel::SerialParam())
    refs <- lapply(1:2,function(j) {
      phi <- compact$dispersion[j]
      lambda <- compact$smoothing_parameters[j,
        compact$geometry$smooth[[compact$geometry$target[["global"]]]]$sp_index]
      inlast_constraint_reference_states(
        raw$A,raw$Q,raw$constraint,raw$projection,
        compact$working_error[,j],compact$working_variance[,j],
        tau=lambda/phi,nuisance_X=compact$geometry$nuisance_design,
        nuisance_Vp=compact$nuisance_covariance[[j]])
    })
    condition <- inlast_constraint_reference_pair(
      refs[[1]]$states$projected,refs[[2]]$states$projected,"davies")
    raw_score <- inlast_constraint_reference_pair(
      refs[[1]]$states$raw_kernel_only,refs[[2]]$states$raw_kernel_only,"davies")
    oracle_error <- abs(condition$p_two_sided-prod$results$p_two_sided[1L])
    if (i <= oracle_reps && (!is.finite(oracle_error) || oracle_error > 1e-7))
      stop("production conditioned oracle mismatch: ",oracle_error)
    spatial_mean <- max(abs(vapply(fits,function(fit) {
      sm <- which(vapply(fit$smooth,inherits,logical(1),"spde.smooth"))[1L]
      cols <- fit$smooth[[sm]]$first.para:fit$smooth[[sm]]$last.para
      L <- mgcvST:::.gam_training_lpmatrix(fit)
      mean(as.numeric(L[,cols,drop=FALSE] %*% coef(fit)[cols]))
    },numeric(1))))
    rows <- data.frame(
      replicate=i,variant=c("conditioned_spde","raw_G_after_P_intercept"),
      calibration="davies",
      p_value=c(condition$p_two_sided,raw_score$p_two_sided),
      signed_score=c(condition$signed_score,raw_score$signed_score),
      fallback=c(!is.null(condition$liu_parameters),!is.null(raw_score$liu_parameters)),
      fit_seconds=proc.time()[["elapsed"]]-t0,
      mean_y1=mean(sim$y[,1]),mean_y2=mean(sim$y[,2]),
      y_checksum1=sum(seq_len(n)*sim$y[,1]),
      y_checksum2=sum(seq_len(n)*sim$y[,2]),
      response_source=response_source,
      theta1=fits[[1]]$family$getTheta(TRUE),theta2=fits[[2]]$family$getTheta(TRUE),
      sp1=fits[[1]]$sp[1],sp2=fits[[2]]$sp[1],
      spatial_mean_error=spatial_mean,oracle_p_error=oracle_error,
      error=NA_character_,stringsAsFactors=FALSE)
    rows
  },error=function(e) failed_rows(i,conditionMessage(e),
                                  proc.time()[["elapsed"]]-t0))
}

index <- seq_len(reps)
if (workers == 1L) {
  result <- do.call(rbind,lapply(index,run_one))
} else {
  pool <- BiocParallel::SnowParam(workers=workers,type="SOCK",stop.on.error=TRUE)
  result <- do.call(rbind,BiocParallel::bplapply(index,run_one,BPPARAM=pool))
}
write.csv(result,file.path(out,"replicates.csv"),row.names=FALSE)
summary <- do.call(rbind,lapply(split(result,result$variant),function(g) {
  valid <- is.finite(g$p_value); nr <- sum(g$p_value[valid] < .05)
  ci <- if(sum(valid)) binom.test(nr,sum(valid))$conf.int else c(NA,NA)
  data.frame(variant=g$variant[1],attempted=nrow(g),valid=sum(valid),
             failures=sum(!valid),fallback=sum(g$fallback,na.rm=TRUE),
             rejected_005=nr,rejection_rate_005=nr/sum(valid),
             ci_lower=ci[1],ci_upper=ci[2],
             mean_fit_seconds=mean(g$fit_seconds),
             max_oracle_p_error=max(g$oracle_p_error,na.rm=TRUE))
}))
write.csv(summary,file.path(out,"summary.csv"),row.names=FALSE)
writeLines(c(
  sprintf("Saved nb03_pair_k6 geometry; mean=%g, rho=%g, seed base=%d, RNGkind L'Ecuyer-CMRG.",mean_count,rho,seed_base),
  "Each replicate refits both features with bam(method='fREML', discrete=TRUE).",
  "The fitted SPDE basis retains the observation mean-zero constraint.",
  "conditioned_spde is checked against the package production pair score for the requested oracle replicates.",
  "raw_G_after_P_intercept retains the same fitted null P and QR-removes its intercept component."
),file.path(out,"protocol.txt"))
print(summary)
