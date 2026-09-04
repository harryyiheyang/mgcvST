# Run with the candidate mgcvST installed. Output location is supplied by caller.
# No mesh-search backend changes or new C++ interpolation kernels are used.
suppressPackageStartupMessages(library(mgcvST))
out <- Sys.getenv("MGCVST_BENCH_OUT", tempdir())
dir.create(out, showWarnings=FALSE, recursive=TRUE)
set.seed(20260904)
data("MISO_E13",package="mgcvST")
basis <- spde_basis(MISO_E13$meshes$spde,
                    as.matrix(MISO_E13$covariates[,c("x","y")]),kappa=.1,
                    pc_cutoff=as.numeric(Sys.getenv("MGCVST_PC_CUTOFF","0.999")))
ns <- asNamespace("mgcvST")
xy <- basis$mesh_vertices
tv <- basis$mesh_triangles
timed <- function(fun,reps=5L) {
  fun() # warmup outside measurement
  if (requireNamespace("microbenchmark",quietly=TRUE)) {
    # Nanosecond timer for short stages on Windows; return seconds.
    return(median(microbenchmark::microbenchmark(fun(),times=20L)$time)/1e9)
  }
  elapsed <- numeric(reps)
  for (r in seq_len(reps)) {
    # Average repeated calls to make sub-millisecond stages measurable with
    # Windows' elapsed-time clock. Identical workload in every repetition.
    gc(FALSE)
    t0 <- proc.time()[["elapsed"]]
    for (j in 1:5) value <- fun()
    elapsed[r] <- (proc.time()[["elapsed"]]-t0)/5
  }
  median(elapsed)
}
results <- list()
for (n in c(100L,1000L,10000L,100000L)) {
  tri <- sample.int(nrow(tv),n,replace=TRUE)
  w <- matrix(rexp(3*n),n,3)
  w <- w/rowSums(w)
  scaled0 <- w[,1]*xy[tv[tri,1],]+w[,2]*xy[tv[tri,2],]+w[,3]*xy[tv[tri,3],]
  raw <- sweep(scaled0*basis$transform$scale,2,basis$transform$center,"+")
  transform <- function() sweep(raw,2,basis$transform$center,"-")/basis$transform$scale
  scaled <- transform()
  searchn <- function() geometry::tsearchn(xy,tv,scaled)
  search2 <- function() geometry::tsearch(xy[,1],xy[,2],tv,scaled[,1],scaled[,2],bary=TRUE,method="quadtree")
  hit <- searchn()
  stopifnot(identical(hit,search2()),!anyNA(hit$idx))
  build_A <- function() Matrix::sparseMatrix(i=rep(seq_len(n),each=3L),
    j=as.vector(t(tv[hit$idx,,drop=FALSE])),x=as.vector(t(hit$p)),dims=c(n,nrow(xy)))
  A <- build_A()
  full <- as.matrix(A %*% basis$projection)
  direct_pc <- as.matrix(A %*% basis$pc_mesh_projection)
  manual_pc <- full %*% basis$pc_vectors[,seq_len(basis$pc_cached_dimension),drop=FALSE]
  stopifnot(isTRUE(all.equal(direct_pc,manual_pc,tolerance=1e-12)))
  result <- data.frame(n_new=n,vertices=nrow(xy),triangles=nrow(tv),
    pc_cutoff=basis$pc_cutoff,
    full_dimension=ncol(basis$B),pc_dimension=basis$pc_cached_dimension,
    coordinate_transform=timed(transform), triangle_lookup_tsearchn=timed(searchn),
    triangle_lookup_quadtree=timed(search2), A_construction=timed(build_A),
    A_times_Z=timed(function() as.matrix(A %*% basis$projection)),
    A_times_ZV=timed(function() as.matrix(A %*% basis$pc_mesh_projection)),
    total_full_prediction=timed(function() ns$.spde_basis_at(basis,raw)),
    total_pc_prediction=timed(function() ns$.spde_basis_at(basis,raw,pc=TRUE)),
    tsearch_identical=identical(hit,search2()),
    pc_reference_max_abs=max(abs(direct_pc-manual_pc)),
    avoided_full_intermediate_bytes=as.numeric(object.size(full)))
  results[[as.character(n)]] <- result
  print(result)
}
write.csv(do.call(rbind,results),file.path(out,"projector-benchmark.csv"),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(out,"projector-session.txt"))
