#!/usr/bin/env Rscript

Sys.setenv(OMP_NUM_THREADS="1",OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="1")
suppressPackageStartupMessages(library(mgcvST))
reps <- as.integer(Sys.getenv("MGCVST_HYPER_REPS","3"))
out <- Sys.getenv("MGCVST_HYPER_OUTPUT",
                  "artifacts/lowcount-investigation/hyperprofile")
stopifnot(reps>=1L,reps<=3L)
dir.create(out,recursive=TRUE,showWarnings=FALSE)

dgp <- readRDS("artifacts/constraint-type1/estimated/nb03_pair_k6-dgp.rds")
d <- dgp$data; n <- nrow(d); size <- 2
basis <- spde_basis(dgp$mesh,as.matrix(d[c("x","y")]),kappa=6,
                    project_intercept=TRUE)
model <- inlaST.set(response~offset(offset0),d,basis,family=mgcv::nb())
raw <- model$inla_spec$random[[1L]]
Z <- raw$projection; B <- as.matrix(raw$A %*% Z)
Qz <- crossprod(Z,as.matrix(raw$Q %*% Z)); r <- ncol(B)
X <- as.matrix(model$inla_spec$fixed$X); p <- ncol(X)
Rz <- chol(Qz); logdetQz <- 2*sum(log(diag(Rz)))
unit_variance_scale <- mean(rowSums((B %*% backsolve(Rz,diag(r)))^2))
scaled_spec <- model$inla_spec
scaled_spec$random[[1L]]$Q <- Matrix::forceSymmetric(
  unit_variance_scale*scaled_spec$random[[1L]]$Q)

simulate_first <- function(i) {
  set.seed(dgp$case$seed+i)
  signal <- dgp$factor_truth %*%
    matrix(rnorm(ncol(dgp$factor_truth)*2L),ncol=2L)
  variance <- rowSums(dgp$factor_truth^2)
  beta <- log(.3)-log(mean(exp(d$offset0+.5*variance)))
  eta <- beta+d$offset0+signal[,1L]
  rnbinom(n,mu=exp(eta),size=size)
}

laplace_at <- function(y,tau,start=NULL) {
  design <- cbind(X,B)
  objective <- function(z) {
    eta <- as.numeric(d$offset0+design%*%z); u <- z[p+seq_len(r)]
    -sum(dnbinom(y,mu=exp(eta),size=size,log=TRUE))+
      .5*tau*drop(crossprod(u,Qz%*%u))
  }
  gradient <- function(z) {
    eta <- as.numeric(d$offset0+design%*%z); mu <- exp(eta)
    residual <- (y+size)*mu/(size+mu)-y
    ans <- as.numeric(crossprod(design,residual))
    u <- z[p+seq_len(r)]
    ans[p+seq_len(r)] <- ans[p+seq_len(r)]+tau*as.numeric(Qz%*%u)
    ans
  }
  if(is.null(start)) start <- numeric(p+r)
  opt <- optim(start,objective,gradient,method="BFGS",
               control=list(maxit=1000,reltol=1e-11))
  z <- opt$par; eta <- as.numeric(d$offset0+design%*%z); mu <- exp(eta)
  w <- (y+size)*size*mu/(size+mu)^2
  H <- crossprod(design,w*design)
  H[p+seq_len(r),p+seq_len(r)] <-
    H[p+seq_len(r),p+seq_len(r)]+tau*Qz
  RH <- chol((H+t(H))/2); logdetH <- 2*sum(log(diag(RH)))
  ll <- sum(dnbinom(y,mu=mu,size=size,log=TRUE))
  log_prior_u <- .5*r*log(tau)+.5*logdetQz-.5*r*log(2*pi)-
    .5*tau*drop(crossprod(z[p+seq_len(r)],Qz%*%z[p+seq_len(r)]))
  log_laplace <- ll+log_prior_u+.5*(p+r)*log(2*pi)-.5*logdetH
  list(value=log_laplace,mode=z,convergence=opt$convergence,
       gradient=max(abs(gradient(z))))
}

theta_grid <- seq(-8,2,by=.5)
rows <- list(); modes <- vector("list",reps)
for(i in seq_len(reps)) {
  y <- simulate_first(i)
  eb <- mgcvST:::.inlast_fit_feature(
    model$inla_spec,y,offset=model$offset,
    control=list(nb_size=size,keep_fit=TRUE))
  eb_scaled <- mgcvST:::.inlast_fit_feature(
    scaled_spec,y,offset=model$offset,
    control=list(nb_size=size,keep_fit=TRUE))
  previous <- NULL
  for(theta in theta_grid) {
    tau <- exp(theta)
    explicit <- laplace_at(y,tau,previous); previous <- explicit$mode
    fixed <- mgcvST:::.inlast_fit_feature(
      model$inla_spec,y,offset=model$offset,
      control=list(nb_size=size,fixed_precision=tau,keep_fit=TRUE))
    mlik <- as.numeric(fixed$inla$mlik[,1L]); names(mlik)<-rownames(fixed$inla$mlik)
    prior <- dnorm(theta,0,3,log=TRUE)
    scaled_prior <- dnorm(theta-log(unit_variance_scale),0,3,log=TRUE)
    rows[[length(rows)+1L]] <- data.frame(
      replicate=i,theta=theta,tau=tau,eb_tau=as.numeric(eb$tau[1L]),
      unit_variance_scale=unit_variance_scale,
      eb_scaled_tau=as.numeric(eb_scaled$tau[1L]),
      eb_scaled_implied_original_tau=
        unit_variance_scale*as.numeric(eb_scaled$tau[1L]),
      inla_mlik_integration=mlik[1L],
      inla_mlik_gaussian=if(length(mlik)>1L)mlik[2L] else NA_real_,
      explicit_laplace=explicit$value,logtau_prior=prior,
      scaled_logtau_prior_on_original_grid=scaled_prior,
      inla_objective=mlik[1L]+prior,
      explicit_objective=explicit$value+prior,
      scaled_explicit_objective=explicit$value+scaled_prior,
      explicit_convergence=explicit$convergence,
      explicit_max_gradient=explicit$gradient,
      fitted_constraint_error=max(abs(fixed$observation_spatial_mean)),
      stringsAsFactors=FALSE)
  }
}
ans <- do.call(rbind,rows)
ans$inla_minus_explicit <- ans$inla_mlik_integration-ans$explicit_laplace
write.csv(ans,file.path(out,"tau-profile.csv"),row.names=FALSE)

summary <- do.call(rbind,lapply(split(ans,ans$replicate),function(g) {
  delta <- g$inla_minus_explicit
  slope <- coef(lm(delta~theta,data=g))[2L]
  data.frame(replicate=g$replicate[1L],eb_tau=g$eb_tau[1L],
    grid_inla_tau=g$tau[which.max(g$inla_objective)],
    grid_explicit_tau=g$tau[which.max(g$explicit_objective)],
    eb_scaled_tau=g$eb_scaled_tau[1L],
    eb_scaled_implied_original_tau=g$eb_scaled_implied_original_tau[1L],
    grid_scaled_implied_original_tau=
      g$tau[which.max(g$scaled_explicit_objective)],
    delta_range=max(delta)-min(delta),delta_theta_slope=unname(slope),
    max_explicit_gradient=max(g$explicit_max_gradient),
    max_constraint_error=max(g$fitted_constraint_error))
}))
write.csv(summary,file.path(out,"profile-summary.csv"),row.names=FALSE)

# Gaussian exact restricted likelihood audit with the same constrained field.
set.seed(91811)
u_gaussian <- backsolve(chol(dgp$tau_truth*Qz),rnorm(r))
y_gaussian <- as.numeric(1+d$offset0+B%*%u_gaussian+rnorm(n,sd=.5))
gaussian_model <- inlaST.set(response~offset(offset0),d,basis,family=gaussian())
gaussian_eb <- mgcvST:::.inlast_fit_feature(
  gaussian_model$inla_spec,y_gaussian,offset=gaussian_model$offset,
  control=list(gaussian_precision=4))
K <- B%*%solve(Qz,t(B))
gaussian_rows <- lapply(theta_grid,function(theta) {
  tau <- exp(theta); V <- diag(.25,n)+K/tau; RV <- chol(V)
  VinvX <- backsolve(RV,forwardsolve(t(RV),X))
  XtVinvX <- crossprod(X,VinvX)
  beta <- solve(XtVinvX,crossprod(VinvX,y_gaussian-d$offset0))
  residual <- y_gaussian-d$offset0-as.numeric(X%*%beta)
  Vinvr <- backsolve(RV,forwardsolve(t(RV),residual))
  reml <- -.5*((n-p)*log(2*pi)+2*sum(log(diag(RV)))+
                 as.numeric(determinant(XtVinvX,logarithm=TRUE)$modulus)+
                 sum(residual*Vinvr))
  fixed <- mgcvST:::.inlast_fit_feature(
    gaussian_model$inla_spec,y_gaussian,offset=gaussian_model$offset,
    control=list(gaussian_precision=4,fixed_precision=tau,keep_fit=TRUE))
  data.frame(theta=theta,tau=tau,eb_tau=as.numeric(gaussian_eb$tau[1L]),
             inla_mlik=as.numeric(fixed$inla$mlik[1L,1L]),
             exact_restricted_loglik=reml,logtau_prior=dnorm(theta,0,3,log=TRUE))
})
gaussian_rows <- do.call(rbind,gaussian_rows)
gaussian_rows$inla_minus_exact <-
  gaussian_rows$inla_mlik-gaussian_rows$exact_restricted_loglik
write.csv(gaussian_rows,file.path(out,"gaussian-profile.csv"),row.names=FALSE)
gaussian_summary <- data.frame(
  eb_tau=gaussian_eb$tau[1L],
  grid_inla_tau=gaussian_rows$tau[which.max(gaussian_rows$inla_mlik+
                                             gaussian_rows$logtau_prior)],
  grid_exact_tau=gaussian_rows$tau[which.max(gaussian_rows$exact_restricted_loglik+
                                              gaussian_rows$logtau_prior)],
  delta_range=diff(range(gaussian_rows$inla_minus_exact)),
  delta_theta_slope=unname(coef(lm(inla_minus_exact~theta,
                                    data=gaussian_rows))[2L]))
write.csv(gaussian_summary,file.path(out,"gaussian-summary.csv"),row.names=FALSE)
writeLines(c(
  "Fixed NB size=2 mechanism profile; latent tau remains governed by log(tau)~N(0,3^2) in the free EB fit.",
  "Explicit Laplace integrates flat fixed effects and the m-1 observation-mean-zero SPDE coordinates.",
  "A missing or duplicated rank adjustment would appear as an approximately +/-0.5 slope of INLA-minus-explicit against log(tau).",
  paste("INLA version",as.character(packageVersion("INLA")))
),file.path(out,"protocol.txt"))
print(summary)
print(gaussian_summary)
