dir.create("artifacts/inla-openmp", recursive = TRUE, showWarnings = FALSE)
library(mgcvST)
set.seed(31)
n <- 512L
q <- 256L
A <- Matrix::sparseMatrix(i = rep(seq_len(n), each = 2L),
  j = as.vector(rbind((seq_len(n) - 1L) %% q + 1L, seq_len(n) %% q + 1L)), x = 0.5)
Q <- Matrix::bandSparse(q, k = c(-1L, 0L, 1L), diagonals = list(rep(-1, q-1L), rep(3,q), rep(-1,q-1L)))
fit <- list(score_sparse=list(A=A,Q=Q,constraint=rep(1,q),target="global",sp_index=1L),
  dispersion=c(1,1), smoothing_parameters=matrix(c(1,2),2,1),
  working_error=matrix(rnorm(n*2),n,2),working_variance=matrix(runif(n*2,0.5,2),n,2),
  geometry=list(nuisance_design=matrix(1,n,1)))
fit <- mgcvST:::.inlast_sparse_prepare(fit)
build <- system.time(units <- mgcvST:::.inlast_sparse_units(fit,1:2,threads=2L))[["elapsed"]]
saveRDS(units,"artifacts/inla-openmp/units-small.rds",compress=FALSE)
restore <- system.time(units2 <- readRDS("artifacts/inla-openmp/units-small.rds"))[["elapsed"]]
elapsed <- numeric(3)
for (i in 1:3) elapsed[i] <- system.time(M <- mgcvST:::.inlast_sparse_materialize(fit,units2,threads=2L))[["elapsed"]]
result <- data.frame(nodes=q,features=2L,unit_bytes=file.info("artifacts/inla-openmp/units-small.rds")$size,
  dense_M_bytes=sum(vapply(M,function(z)as.numeric(object.size(z$M)),numeric(1))),
  build_seconds=build,read_seconds=restore,materialize_seconds=median(elapsed))
print(result)
stopifnot(all(vapply(units2,function(z)is.null(z$M),logical(1))))
write.csv(result,"artifacts/inla-openmp/reconstruction.csv",row.names=FALSE)
