# Same-fit diagnostic: do not refit, change priors, or replace production Vp.
library(mgcvST)
out <- Sys.getenv("MGCVST_PROJECTION_OUT", "artifacts/score-projection-check")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
fit <- readRDS("artifacts/nb-flat-ten/default-smoke-fit.rds")
dense <- fit
dense$score_backend <- "dense"
dense$.mgcvst_fixed_factors <- NULL
dense$.mgcvst_state_cache <- NULL
n <- nrow(fit$working_error)
center <- function(F) sweep(F, 2L, colMeans(F), "-")
Fraw0 <- mgcvST:::.mgcvst_spde_factor_base(fit$score_sparse$A, fit$score_sparse$Q)
rows <- list(); states <- list(); comparisons <- list()
for (vp in c("native", "expected_diagnostic_only")) {
  used <- dense
  if (vp != "native") used$nuisance_covariance <- fit$expected_nuisance_covariance
  for (j in seq_len(ncol(fit$working_error))) {
    op <- mgcvST:::.mgcvst_model_operator_vp(used, j)
    P <- mgcvST:::.mgcvst_model_apply_P(op$operator, diag(n))
    eig <- eigen((P + t(P))/2, symmetric = TRUE)
    eigenP <- eig$values
    # Verify the user's square-root identity directly only when P is PSD.
    # No eigenvalue clipping or covariance replacement is used here.
    sqrtP <- if (min(eigenP) >= 0) {
      sweep(eig$vectors, 2L, sqrt(eigenP), "*") %*% t(eig$vectors)
    } else NULL
    raw <- Fraw0 * sqrt(fit$dispersion[j]/fit$smoothing_parameters[j,1])
    conditioned <- op$target[[1]]
    factors <- list(raw = raw, raw_centered = center(raw),
                    conditioned = conditioned, conditioned_centered = center(conditioned))
    for (kernel in names(factors)) {
      F <- factors[[kernel]]
      PF <- P %*% F
      a <- drop(crossprod(F, P %*% fit$working_error[,j]))
      M <- crossprod(F, PF); M <- (M + t(M))/2
      states[[paste(vp,j,kernel,sep="/")]] <- list(a=a,M=M)
      rows[[length(rows)+1L]] <- data.frame(Vp=vp,feature=j,kernel=kernel,
        max_P1=max(abs(P %*% rep(1,n))), min_eigen_P=min(eigenP),
        max_asymmetry_P=max(abs(P-t(P))), max_factor_column_mean=max(abs(colMeans(F))),
        factor_rank=qr(F)$rank, norm_score_vector=sqrt(sum(a*a)),
        trace_score_M=sum(diag(M)), min_eigen_M=min(eigen(M,symmetric=TRUE,only.values=TRUE)$values),
        max_FtP1=max(abs(crossprod(F,P %*% rep(1,n)))))
    }
    getstate <- function(k) states[[paste(vp,j,k,sep="/")]]
    for (kernel in c("raw", "conditioned")) {
      s <- getstate(kernel); sc <- getstate(paste0(kernel,"_centered"))
      comparisons[[length(comparisons)+1L]] <- data.frame(Vp=vp,feature=j,kernel=kernel,
        max_a_change=max(abs(s$a-sc$a)),max_M_change=max(abs(s$M-sc$M)),
        relative_M_change=max(abs(s$M-sc$M))/max(abs(s$M)),
        max_sqrtP_G_sqrtP_change=if (is.null(sqrtP)) NA_real_ else {
          F <- factors[[kernel]]; Fc <- factors[[paste0(kernel,"_centered")]]
          max(abs(tcrossprod(sqrtP %*% F) - tcrossprod(sqrtP %*% Fc)))
        })
    }
  }
}
pairs <- list()
for (vp in c("native", "expected_diagnostic_only")) for (kernel in names(factors)) {
  a <- states[[paste(vp,1,kernel,sep="/")]]
  b <- states[[paste(vp,2,kernel,sep="/")]]
  score <- sum(a$a*b$a)
  cal <- rkhs_score_calibrate(score,a$M,b$M,method="davies")
  pairs[[length(pairs)+1L]] <- data.frame(Vp=vp,kernel=kernel,score=score,
    information=cal$information,p_value=cal$p_two_sided)
}
sparse <- lapply(1:2,function(j) mgcvST:::.mgcvst_model_score_state(fit,j))
sparse_check <- do.call(rbind,lapply(1:2,function(j) {
  d <- states[[paste("native",j,"conditioned",sep="/")]]
  data.frame(feature=j,max_a_difference=max(abs(d$a-sparse[[j]]$a)),
             max_M_difference=max(abs(d$M-sparse[[j]]$M)))
}))
diagnostics <- do.call(rbind,rows)
comparison <- do.call(rbind,comparisons)
pair_results <- do.call(rbind,pairs)
stopifnot(max(abs(fit$observation_spatial_mean)) < 1e-10,
          max(sparse_check$max_a_difference) < 1e-9,
          max(sparse_check$max_M_difference) < 1e-9,
          max(comparison$max_M_change[comparison$kernel=="conditioned"]) < 1e-9,
          max(comparison$max_a_change[comparison$kernel=="conditioned"]) < 1e-9)
write.csv(diagnostics,file.path(out,"operator-diagnostics.csv"),row.names=FALSE)
write.csv(comparison,file.path(out,"centering-comparison.csv"),row.names=FALSE)
write.csv(pair_results,file.path(out,"pair-results.csv"),row.names=FALSE)
write.csv(sparse_check,file.path(out,"sparse-dense-check.csv"),row.names=FALSE)
saveRDS(list(diagnostics=diagnostics,comparison=comparison,pairs=pair_results,
             sparse_check=sparse_check,source_fit="artifacts/nb-flat-ten/default-smoke-fit.rds",
             session=sessionInfo()),file.path(out,"results.rds"))
print(diagnostics, digits=8, row.names=FALSE)
print(comparison, digits=8, row.names=FALSE)
print(pair_results, digits=10, row.names=FALSE)
print(sparse_check, digits=4, row.names=FALSE)
