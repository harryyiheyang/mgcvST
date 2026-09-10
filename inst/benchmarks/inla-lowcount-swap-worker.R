lowcount_swap_one <- function(i, payload) {
  RNGkind("L'Ecuyer-CMRG")
  set.seed(61000L+i)
  n <- nrow(payload$d)
  signal <- payload$truth %*% matrix(rnorm(ncol(payload$truth)*2L),ncol=2L)
  beta <- log(.3)-log(mean(exp(payload$offset+.5*rowSums(payload$truth^2))))
  ys <- lapply(1:2,function(j) rnbinom(n,mu=exp(beta+.25*(j-1)+payload$offset+signal[,j]),size=2))
  historical_normal <- list(prior="normal",param=c(0,1/9),initial=0)
  inla <- lapply(ys,function(y) mgcvST:::.inlast_fit_feature(
    payload$spec,y,offset=payload$offset,
    control=list(precision_prior=historical_normal,
                 nb_size_prior=historical_normal)))
  bam <- lapply(ys,function(y) {
    G <- payload$G; G$y <- y; G$mf[[payload$response_index]] <- y
    G$family <- unserialize(payload$family_raw)
    mgcv::bam(G=G,method="fREML",discrete=TRUE,nthreads=1L)
  })
  stopifnot(all(vapply(inla,function(x)isTRUE(x$converged),logical(1))),
            all(vapply(bam,function(x)isTRUE(x$converged),logical(1))))
  extracts <- lapply(1:2,function(j) list(
    inla=list(eta=inla[[j]]$eta,tau=as.numeric(inla[[j]]$tau[1]),
              size=inla[[j]]$family_parameters[1],Vp=inla[[j]]$nuisance_covariance),
    bam=list(eta=as.numeric(bam[[j]]$linear.predictors),tau=as.numeric(bam[[j]]$sp[1])/bam[[j]]$sig2,
             size=bam[[j]]$family$getTheta(TRUE),Vp=bam[[j]]$Vp[1,1,drop=FALSE])))
  if (max(abs(unlist(lapply(inla,`[[`,"observation_spatial_mean"))))>1e-8)
    stop("INLA mean constraint failure")
  for(j in 1:2) {
    L <- mgcvST:::.gam_training_lpmatrix(bam[[j]])
    field <- L[,-1,drop=FALSE] %*% coef(bam[[j]])[-1]
    if(abs(mean(field))>1e-8) stop("bam mean constraint failure")
  }
  saveRDS(list(replicate=i,y=ys,extracts=extracts,signal=signal,
              inla=inla),file.path(payload$cache,sprintf("rep-%04d.rds",i)))
  design <- expand.grid(eta=c("inla","bam"),tau=c("inla","bam"),
                        size=c("inla","bam"),stringsAsFactors=FALSE)
  states <- vector("list",nrow(design)); rows <- vector("list",nrow(design)+1L)
  row_from_pair <- function(pair,variant,eta,tau,size,leak=0) data.frame(
    replicate=i,variant=variant,eta_source=eta,tau_source=tau,size_source=size,
    p_value=pair$p_two_sided,U=pair$signed_score,
    information=pair$information,fallback=!is.null(pair$liu_parameters),
    nuisance_leakage=leak,error=NA_character_)
  for(k in seq_len(nrow(design))) {
    v <- design[k,]
    states[[k]] <- lapply(1:2,function(j) {
      e <- extracts[[j]][[v$eta]]$eta; mu <- exp(e)
      r <- extracts[[j]][[v$size]]$size; tau <- extracts[[j]][[v$tau]]$tau
      inlast_raw_kernel_sparse_state(payload$A,payload$Q,payload$g,
                e+(ys[[j]]-mu)/mu-payload$offset,1/mu+1/r,tau,payload$X)
    })
    pair <- inlast_raw_kernel_sparse_pair(states[[k]][[1]],states[[k]][[2]],"davies")
    rows[[k]] <- row_from_pair(pair,paste(v,collapse="_"),v$eta,v$tau,v$size)
  }
  brefs <- lapply(1:2,function(j) {
    z <- extracts[[j]]$bam; mu <- exp(z$eta)
    inlast_constraint_reference_states(payload$A,payload$Q,payload$g,payload$Z,
      z$eta+(ys[[j]]-mu)/mu-payload$offset,1/mu+1/z$size,z$tau,
      nuisance_X=payload$X,nuisance_Vp=z$Vp)$states$raw_kernel_only
  })
  bp <- inlast_constraint_reference_pair(brefs[[1]],brefs[[2]],"davies")
  rows[[nrow(design)+1L]] <- row_from_pair(bp,"bam_original_Vp","bam","bam","bam",
     max(vapply(brefs,`[[`,numeric(1),"nuisance_annihilation_error")))
  answer <- do.call(rbind,rows)
  meta <- do.call(rbind,lapply(1:2,function(j) {
    bi <- which(design$eta=="bam" & design$tau=="bam" & design$size=="bam")
    data.frame(replicate=i,feature=j,tau_inla=extracts[[j]]$inla$tau,
      tau_bam=extracts[[j]]$bam$tau,size_inla=extracts[[j]]$inla$size,
      size_bam=extracts[[j]]$bam$size,
      eta_rmse=sqrt(mean((extracts[[j]]$inla$eta-extracts[[j]]$bam$eta)^2)),
      Vp_bam=extracts[[j]]$bam$Vp[1,1],
      Vp_bam_expected=states[[bi]][[j]]$nuisance_covariance[1,1],
      Vp_inla=extracts[[j]]$inla$Vp[1,1])
  }))
  list(rows=answer,metadata=meta)
}

lowcount_swap_block <- function(indices,payload) {
  lapply(indices,function(i) tryCatch(lowcount_swap_one(i,payload),
     error=function(e) list(failure=data.frame(replicate=i,error=conditionMessage(e)))))
}
