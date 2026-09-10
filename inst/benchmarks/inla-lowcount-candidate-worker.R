lowcount_candidate_score <- function(fit, y, payload, tau_original, Vp, working, kernel) {
  mu <- fit$mu; r <- fit$family_parameters[1]; eta <- fit$eta
  if(working=="posterior_observed") {
    w <- r*mu*(y+r)/(r+mu)^2
    e <- eta+(y-mu)*(r+mu)/(mu*(y+r))-payload$offset
  } else {
    w <- 1/(1/mu+1/r)
    e <- eta+(y-mu)/mu-payload$offset
  }
  F0 <- payload$Fconditioned/sqrt(tau_original)
  X <- payload$X
  wf <- w*F0
  middle <- chol(diag(ncol(F0))+crossprod(F0,wf))
  Vsolve <- function(B) {
    WB <- w*B
    WB-wf%*%backsolve(middle,forwardsolve(t(middle),crossprod(F0,WB)))
  }
  WX <- Vsolve(X); glsVp <- solve(crossprod(X,WX))
  usedVp <- if(working=="legacy_expected") glsVp else Vp
  P <- function(B) {
    WB <- Vsolve(B)
    WB-WX%*%usedVp%*%crossprod(X,WB)
  }
  F <- switch(kernel,conditioned=payload$Fconditioned,
              raw=payload$Fraw,raw_centered=payload$Fraw_centered)/sqrt(tau_original)
  a <- as.numeric(crossprod(F,P(matrix(e,ncol=1))))
  M <- crossprod(F,P(F));M<-(M+t(M))/2
  # The sandwich is retained for diagnostics; it is not silently substituted.
  PF <- P(F)
  M_sandwich <- crossprod(PF,(1/w)*PF)+crossprod(crossprod(F0,PF))
  list(a=a,M=M,diagnostics=list(P1=max(abs(P(X))),
    Vp_ratio=Vp[1,1]/glsVp[1,1],minimum_eigenvalue=min(eigen(M,symmetric=TRUE,only.values=TRUE)$values),
    sandwich_relative_error=max(abs(M-M_sandwich))/max(1,max(abs(M_sandwich)))))
}

lowcount_candidate_one <- function(i,payload) {
  RNGkind("L'Ecuyer-CMRG");set.seed(payload$seed_base+i)
  n <- nrow(payload$X)
  innovation <- matrix(rnorm(ncol(payload$truth)*2L),ncol=2L)
  if(payload$rho!=0) innovation[,2]<-payload$rho*innovation[,1]+sqrt(1-payload$rho^2)*innovation[,2]
  signal <- payload$truth%*%innovation
  beta <- log(payload$mean_count)-log(mean(exp(payload$offset+.5*rowSums(payload$truth^2))))
  ys <- lapply(1:2,function(j) rnbinom(n,mu=exp(beta+.25*(j-1)+payload$offset+signal[,j]),size=2))
  rows <- list(); metadata <- list(); position<-0L
  for(scale_name in payload$scale_names) {
    scale <- if(scale_name=="original") 1 else payload$Q_scale
    spec <- payload$spec;spec$random[[1]]$Q <- scale*spec$random[[1]]$Q
    historical_normal <- list(prior="normal",param=c(0,1/9),initial=0)
    fit_control <- list(keep_fit=TRUE,
      precision_prior=historical_normal,nb_size_prior=historical_normal)
    if(payload$fixed_truth) {
      fit_control$fixed_precision <- payload$tau_truth/scale
      fit_control$nb_size <- 2
    }
    fits <- lapply(ys,function(y) mgcvST:::.inlast_fit_feature(spec,y,offset=payload$offset,
                                       control=fit_control))
    stopifnot(all(vapply(fits,function(x)isTRUE(x$converged),logical(1))))
    posts <- lapply(fits,function(x) inlast_posterior_vp(x$inla,spec))
    for(j in 1:2) {
      if(max(abs(fits[[j]]$observation_spatial_mean))>1e-8) stop("mean constraint failure")
      metadata[[length(metadata)+1L]]<-data.frame(replicate=i,scale=scale_name,feature=j,
        tau_internal=fits[[j]]$tau[1],tau_original=scale*fits[[j]]$tau[1],Q_scale=scale,
        nb_size=fits[[j]]$family_parameters[1],mean_error=max(abs(fits[[j]]$observation_spatial_mean)),
        Vp_posterior=posts[[j]]$nuisance_covariance[1,1],
        Vp_legacy=fits[[j]]$nuisance_covariance[1,1])
    }
    if(i<=2) saveRDS(list(fits=fits,posts=posts,y=ys,spec=spec),
       file.path(payload$cache,sprintf("rep-%04d-%s.rds",i,scale_name)))
    for(working in payload$working_names) for(kernel in payload$kernel_names) {
      position<-position+1L
      rows[[position]] <- tryCatch({
        states<-lapply(1:2,function(j) lowcount_candidate_score(fits[[j]],ys[[j]],payload,
                  scale*as.numeric(fits[[j]]$tau[1]),posts[[j]]$nuisance_covariance,working,kernel))
        U <- sum(states[[1]]$a*states[[2]]$a)
        cal <- mgcvST::rkhs_score_calibrate(U,states[[1]]$M,states[[2]]$M,method="davies")
        data.frame(replicate=i,scale=scale_name,working=working,kernel=kernel,
          p_value=cal$p_two_sided,U=U,information=cal$information,
          fallback=!is.null(cal$liu_parameters),P1=max(vapply(states,function(x)x$diagnostics$P1,numeric(1))),
          Vp_ratio_min=min(vapply(states,function(x)x$diagnostics$Vp_ratio,numeric(1))),
          Vp_ratio_max=max(vapply(states,function(x)x$diagnostics$Vp_ratio,numeric(1))),
          minimum_eigenvalue=min(vapply(states,function(x)x$diagnostics$minimum_eigenvalue,numeric(1))),
          sandwich_relative_error=max(vapply(states,function(x)x$diagnostics$sandwich_relative_error,numeric(1))),
          error=NA_character_)
      },error=function(e) data.frame(replicate=i,scale=scale_name,working=working,kernel=kernel,
           p_value=NA_real_,U=NA_real_,information=NA_real_,fallback=NA,P1=NA_real_,
           Vp_ratio_min=NA_real_,Vp_ratio_max=NA_real_,minimum_eigenvalue=NA_real_,
           sandwich_relative_error=NA_real_,error=conditionMessage(e)))
    }
  }
  list(rows=do.call(rbind,rows),metadata=do.call(rbind,metadata))
}

lowcount_candidate_block <- function(indices,payload) {
  lapply(indices,function(i) tryCatch(lowcount_candidate_one(i,payload),
      error=function(e) list(failure=data.frame(replicate=i,error=conditionMessage(e)))))
}
