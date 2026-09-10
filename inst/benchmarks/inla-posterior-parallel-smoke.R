# Verify current posterior-Vp and observation-scale objects on fresh SOCK workers.
Sys.setenv(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
suppressPackageStartupMessages(library(mgcvST))
set.seed(1402)
vertices <- as.matrix(expand.grid(x=seq(0,1,length.out=5),y=seq(0,1,length.out=5)))
mesh <- list(loc=vertices,graph=list(tv=geometry::delaunayn(vertices)))
d <- data.frame(x=runif(80,.01,.99),y=runif(80,.01,.99),offset0=0)
basis <- spde_basis(mesh,as.matrix(d[c("x","y")]),kappa=.7,project_intercept=TRUE)
model <- inlaST.set(response~1,d,basis,family=mgcv::nb(),precision_scale="observation")
Y <- rbind(a=rnbinom(80,mu=exp(sin(4*d$x)),size=2),
           b=rnbinom(80,mu=exp(cos(4*d$y)),size=2))
a <- inlaST.estimate(Y,model,BPPARAM=BiocParallel::SerialParam())
bp <- BiocParallel::SnowParam(2,type="SOCK")
b <- tryCatch(inlaST.estimate(Y,model,BPPARAM=bp,chunk_size=1),
              finally=BiocParallel::bpstop(bp))
# Estimated INLA modes vary slightly across processes; compare at 1e-5.
print(all.equal(a$nuisance_covariance,b$nuisance_covariance,tolerance=1e-5))
print(all.equal(a$working_error,b$working_error,tolerance=1e-5))
saveRDS(list(serial=a,parallel=b),"artifacts/posterior-parallel-smoke.rds")
stopifnot(all(a$diagnostics$converged),all(b$diagnostics$converged),
          isTRUE(all.equal(a$nuisance_covariance,b$nuisance_covariance,tolerance=1e-5)),
          isTRUE(all.equal(a$working_error,b$working_error,tolerance=1e-5)),
          max(abs(b$observation_spatial_mean))<1e-10)
print(b$diagnostics)
cat("Serial/SOCK native posterior Vp and working-error agreement: passed\n")
