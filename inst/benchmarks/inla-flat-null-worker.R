flat_null_one <- function(i,payload) {
  RNGkind("L'Ecuyer-CMRG"); set.seed(payload$seed_base+i)
  innovation <- matrix(rnorm(ncol(payload$truth)*2),ncol=2)
  if(payload$rho!=0) innovation[,2] <- payload$rho*innovation[,1]+sqrt(1-payload$rho^2)*innovation[,2]
  signal <- payload$truth%*%innovation
  beta <- log(payload$mean_count)-log(mean(exp(payload$offset+.5*rowSums(payload$truth^2))))
  ys <- lapply(1:2,function(j) rnbinom(nrow(payload$X),
        mu=exp(beta+.25*(j-1)+payload$offset+signal[,j]),size=2))
  if(payload$seed_base==161000 && payload$rho==0 && payload$mean_count==.3 && i<=2) {
    old <- readRDS(file.path(payload$old_cache,sprintf("rep-%04d-original.rds",i)))
    stopifnot(identical(ys,old$y))
  }
  rows <- list(); metadata <- list()
  for(variant in payload$variants) {
    result <- tryCatch({
      spec <- payload$spec
      if(grepl("scaled$",variant)) spec$random[[1]]$precision_scale <- payload$Q_scale
      flat <- list(prior="flat",param=numeric(),initial=0)
      historical_normal <- list(prior="normal",param=c(0,1/9),initial=0)
      ctl <- list(keep_fit=TRUE,precision_prior=flat,
                  nb_size_prior=historical_normal)
      if(grepl("both",variant)) ctl$nb_size_prior <- flat
      t0 <- proc.time()[["elapsed"]]
      fits <- lapply(ys,function(y) mgcvST:::.inlast_fit_feature(spec,y,payload$offset,ctl))
      elapsed <- proc.time()[["elapsed"]]-t0
      for(j in 1:2) {
        fit <- fits[[j]]
        metadata[[length(metadata)+1]] <- data.frame(replicate=i,variant=variant,feature=j,
          tau=fit$tau[1],tau_internal=fit$tau_internal[1],nb_size=fit$family_parameters[1],
          average_field_variance=payload$Q_scale/fit$tau[1],
          mean_error=max(abs(fit$observation_spatial_mean)),converged=fit$converged,
          mode_status=fit$mode_status,Vp=fit$nuisance_covariance[1,1],fit_seconds=fit$fit_seconds)
      }
      if(i<=2) saveRDS(list(y=ys,fits=fits,spec=spec),
          file.path(payload$cache,sprintf("rep-%04d-%s.rds",i,variant)))
      if(!all(vapply(fits,function(x)isTRUE(x$converged),logical(1)))) stop("INLA mode status not converged")
      if(max(vapply(fits,function(x)max(abs(x$observation_spatial_mean)),numeric(1)))>1e-8) stop("mean-zero constraint failed")
      do.call(rbind,lapply(payload$kernels,function(kernel) {
        tryCatch({
          states <- lapply(1:2,function(j) lowcount_candidate_score(fits[[j]],ys[[j]],payload,
             as.numeric(fits[[j]]$tau[1]),fits[[j]]$nuisance_covariance,"posterior_expected",kernel))
          U <- sum(states[[1]]$a*states[[2]]$a)
          cal <- mgcvST::rkhs_score_calibrate(U,states[[1]]$M,states[[2]]$M,method="davies")
          data.frame(replicate=i,variant=variant,kernel=kernel,p_value=cal$p_two_sided,
            information=cal$information,fallback=!is.null(cal$liu_parameters),elapsed=elapsed,
            maximum_P1=max(vapply(states,function(x)x$diagnostics$P1,numeric(1))),
            error=if(is.finite(cal$p_two_sided)) NA_character_ else
              "score information non-finite or below package threshold")
        },error=function(e) data.frame(replicate=i,variant=variant,kernel=kernel,p_value=NA_real_,
              information=NA_real_,fallback=NA,elapsed=elapsed,maximum_P1=NA_real_,error=conditionMessage(e)))
      }))
    },error=function(e) data.frame(replicate=i,variant=variant,kernel=payload$kernels,
       p_value=NA_real_,information=NA_real_,fallback=NA,elapsed=NA_real_,maximum_P1=NA_real_,error=conditionMessage(e)))
    rows[[length(rows)+1]] <- result
  }
  answer <- list(rows=do.call(rbind,rows),metadata=do.call(rbind,metadata))
  saveRDS(answer,file.path(payload$results,sprintf("rep-%04d.rds",i)))
  answer
}

flat_null_block <- function(indices,payload) lapply(indices,function(i) flat_null_one(i,payload))
