# Worker functions for the estimated-hyperparameter constraint/null study.
# Source together with inla-constraint-operators.R into one environment.

.inlast_mc_marginal <- function(engine, payload, raw_kernel, method) {
  F <- if (raw_kernel) payload$factor_raw else payload$factor_centered
  Dinv <- 1 / engine$working_variance
  X <- payload$X
  e <- engine$working_error
  if (payload$case$family == "negative_binomial") {
    keep <- is.finite(Dinv) & Dinv > 1e-12
    F <- F[keep,,drop=FALSE]; X <- X[keep,,drop=FALSE]
    e <- e[keep]; Dinv <- Dinv[keep]
  }
  WX <- Dinv * X
  Vp <- solve(crossprod(X, WX))
  Papply <- function(value) {
    weighted <- Dinv * value
    weighted - WX %*% Vp %*% crossprod(X, weighted)
  }
  a <- as.numeric(crossprod(F, Papply(matrix(e, ncol=1))))
  H <- crossprod(F, Papply(F))
  spectrum <- eigen((H + t(H))/2, symmetric=TRUE, only.values=TRUE)$values
  if (min(spectrum) < -1e-8 * max(1, abs(spectrum))) stop("Negative marginal eigenvalue")
  spectrum <- spectrum[spectrum > 1e-12 * max(spectrum)]
  statistic <- sum(a^2)
  moments <- vapply(1:4, function(k) sum(spectrum^k), numeric(1L))
  fallback <- FALSE
  p <- if (method == "liu") {
    mgcvST:::.mgcvst_marginal_liu(statistic, moments)
  } else {
    answer <- CompQuadForm::davies(statistic, lambda=spectrum, lim=10000L, acc=1e-6)
    if (answer$ifault != 0L || !is.finite(answer$Qq) || answer$Qq <= 0 || answer$Qq > 1) {
      fallback <- TRUE
      mgcvST:::.mgcvst_marginal_liu(statistic, moments)
    } else answer$Qq
  }
  list(statistic=statistic, p=p, fallback=fallback)
}

.inlast_mc_one <- function(rep, payload) {
  case <- payload$case
  set.seed(case$seed + rep)
  nf <- if (case$test == "pair") 2L else 1L
  n <- nrow(payload$A)
  signal <- if (case$test == "pair") {
    payload$factor_truth %*% matrix(rnorm(ncol(payload$factor_truth)*nf), ncol=nf)
  } else matrix(0, n, nf)
  if (max(abs(colMeans(signal))) > 1e-10) stop("DGP spatial mean is not zero")
  response <- lapply(seq_len(nf), function(j) {
    beta <- if (case$family == "gaussian") 1 else {
      variance <- if (case$test == "pair") rowSums(payload$factor_truth^2) else rep(0,n)
      log(case$mean_count) - log(mean(exp(payload$offset + 0.5*variance)))
    }
    eta <- beta + 0.25*(j-1) + payload$offset + signal[,j]
    if (case$family == "gaussian") eta + rnorm(n, sd=sqrt(case$phi)) else
      rnbinom(n, mu=exp(eta), size=case$nb_size)
  })
  historical_normal <- list(prior="normal",param=c(0,1/9),initial=0)
  fits <- lapply(response, function(y) {
    mgcvST:::.inlast_fit_feature(
      payload$spec, y, offset=payload$offset,
      control=list(precision_prior=historical_normal,
                   nb_size_prior=historical_normal)
    )
  })
  if (!all(vapply(fits, function(z) isTRUE(z$converged), logical(1L)))) stop("INLA convergence flag failed")
  mean_error <- max(abs(unlist(lapply(fits, `[[`, "observation_spatial_mean"))))
  if (mean_error > 1e-8) stop("Fitted mean-zero constraint failed")
  parameter <- function(name, i) {
    if (length(fits) < i) return(NA_real_)
    v <- fits[[i]][[name]]
    if (!length(v)) NA_real_ else as.numeric(v[1])
  }
  metadata <- list(case=case$name, replicate=rep, test=case$test, family=case$family,
                   n=n, mesh_vertices=ncol(payload$A), kappa=case$kappa,
                   mean_count=case$mean_count, fitted_mean_error=mean_error,
                   minimum_working_weight=min(unlist(lapply(fits,function(z) 1/z$working_variance))),
                   tau1=parameter("tau",1),tau2=parameter("tau",2),
                   dispersion1=parameter("dispersion",1),dispersion2=parameter("dispersion",2),
                   nb_size1=parameter("family_parameters",1),nb_size2=parameter("family_parameters",2))
  if (case$test == "pair") {
    refs <- lapply(fits, function(z) {
      inlast_constraint_reference_states(
        payload$A, payload$Q, payload$g, payload$Z,
        working_error=z$working_error, working_variance=z$working_variance,
        tau=as.numeric(z$tau[1]), nuisance_X=payload$X,
        nuisance_Vp=z$nuisance_covariance)
    })
    equivalence <- inlast_constraint_reference_equivalence(refs[[1]], refs[[2]], tolerance=1e-7)
    if (!equivalence$passed) stop("Projected/raw-constrained invariant mismatch")
    rows <- lapply(payload$variants, function(variant) {
      do.call(rbind, lapply(c("liu","davies"), function(method) {
        cal <- tryCatch(inlast_constraint_reference_pair(refs[[1]]$states[[variant]],
                                                refs[[2]]$states[[variant]],method=method),
                        error=function(e) e)
        if (inherits(cal,"condition")) return(data.frame(
          variant=variant,calibration=method,p_value=NA_real_,statistic=NA_real_,
          fallback=NA,equivalence_error=max(equivalence$relative_errors),
          nuisance_leakage=NA_real_,error=conditionMessage(cal)))
        data.frame(variant=variant, calibration=method,
                   p_value=cal$p_two_sided, statistic=cal$statistic,
                   fallback=method=="davies" && !is.null(cal$liu_parameters),
                   equivalence_error=max(equivalence$relative_errors),
                   nuisance_leakage=max(refs[[1]]$states[[variant]]$nuisance_annihilation_error,
                                        refs[[2]]$states[[variant]]$nuisance_annihilation_error),
                   error=NA_character_)
      }))
    })
    answer <- do.call(rbind,rows)
  } else {
    answer <- do.call(rbind,lapply(c(FALSE,TRUE), function(raw_kernel) {
      do.call(rbind,lapply(c("liu","davies"), function(method) {
        cal <- tryCatch(.inlast_mc_marginal(fits[[1]],payload,raw_kernel,method),
                        error=function(e) e)
        if (inherits(cal,"condition")) return(data.frame(
          variant=if(raw_kernel) "raw_kernel_only" else "projected",calibration=method,
          p_value=NA_real_,statistic=NA_real_,fallback=NA,equivalence_error=NA_real_,
          nuisance_leakage=NA_real_,error=conditionMessage(cal)))
        data.frame(variant=if(raw_kernel) "raw_kernel_only" else "projected",
                   calibration=method,p_value=cal$p,statistic=cal$statistic,
                   fallback=cal$fallback,equivalence_error=NA_real_,
                   nuisance_leakage=NA_real_,error=NA_character_)
      }))
    }))
  }
  cbind(as.data.frame(metadata),answer)
}

.inlast_mc_failed <- function(rep, payload, error) {
  variants <- if(payload$case$test=="pair") payload$variants else c("projected","raw_kernel_only")
  cells <- expand.grid(variant=variants,calibration=c("liu","davies"),stringsAsFactors=FALSE)
  data.frame(case=payload$case$name,replicate=rep,test=payload$case$test,
             family=payload$case$family,n=nrow(payload$A),mesh_vertices=ncol(payload$A),
             kappa=payload$case$kappa,mean_count=payload$case$mean_count,
             fitted_mean_error=NA_real_,minimum_working_weight=NA_real_,tau1=NA_real_,tau2=NA_real_,
             dispersion1=NA_real_,dispersion2=NA_real_,nb_size1=NA_real_,nb_size2=NA_real_,
             cells,p_value=NA_real_,statistic=NA_real_,fallback=NA,
             equivalence_error=NA_real_,nuisance_leakage=NA_real_,
             error=conditionMessage(error))
}

inlast_constraint_mc_block <- function(indices,payload) {
  output <- lapply(indices,function(i) {
    tryCatch(.inlast_mc_one(i,payload),error=function(e) .inlast_mc_failed(i,payload,e))
  })
  cat(sprintf("Completed %s replicates %d-%d\n",payload$case$name,min(indices),max(indices)))
  do.call(rbind,output)
}
