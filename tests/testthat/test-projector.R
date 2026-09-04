test_that("recomputed training predictors agree with the old implementation", {
  old <- old_st()
  for (pc in c(FALSE, TRUE)) {
    f <- st_fixture(pc = pc)
    sm <- f$G$smooth[[1L]]
    method <- if (pc) "Predict.matrix.spdePC.smooth" else "Predict.matrix.spde.smooth"
    old_pred <- old[[method]](sm,f$data)
    new_pred <- get(method,asNamespace("mgcvST"))(sm,f$data)
    expect_equal(new_pred,old_pred,tolerance=1e-12)
    expect_equal(mgcvST:::.spde_basis_at(f$basis, f$basis$coordinates, pc), old_pred,tolerance=1e-12)
    if (pc) {
      expect_identical(sm$score_basis, f$basis$B)
      expect_identical(ncol(sm$score_basis), ncol(f$basis$projection))
      expect_lt(ncol(new_pred), ncol(sm$score_basis))
    }
  }
})

test_that("new, reordered, subset and empty coordinates use the correct space", {
  f <- st_fixture()
  basis <- f$basis
  new <- rbind(c(.123,.456), c(.891,.654))
  mesh <- list(xy=basis$mesh_vertices,tv=basis$mesh_triangles)
  A <- mgcvST:::.spde_basis_project(mesh,new)
  Z <- basis$projection
  q <- basis$pc_cached_dimension
  V <- basis$pc_vectors[,seq_len(q),drop=FALSE]
  expect_equal(mgcvST:::.spde_basis_at(basis,new),as.matrix(A %*% Z),tolerance=1e-13)
  expect_equal(mgcvST:::.spde_basis_at(basis,new,TRUE),as.matrix((A %*% Z) %*% V),tolerance=1e-12)
  expect_equal(basis$pc_mesh_projection, Z %*% V, tolerance=1e-13)
  for (pc in c(FALSE,TRUE)) {
    reference <- if (pc) basis$pc_training_basis else basis$B
    for (rows in list(2L,c(8L,1L),rev(seq_len(nrow(basis$B))))) {
      expect_equal(mgcvST:::.spde_basis_at(basis,basis$coordinates[rows,,drop=FALSE],pc),
                       reference[rows,,drop=FALSE],tolerance=1e-12)
    }
    expect_identical(dim(mgcvST:::.spde_basis_at(basis,matrix(numeric(),0,2),pc)),c(0L,ncol(reference)))
    expect_error(mgcvST:::.spde_basis_at(basis,matrix(c(-10,10),1,2),pc),"outside")
    expect_error(mgcvST:::.spde_basis_at(basis,matrix(c(NA,1),1,2),pc),"finite")
  }
})

test_that("raw new coordinates use the saved nontrivial transform", {
  f <- st_fixture()
  transform <- list(center=c(125,-370),scale=31)
  mesh <- structure(list(mesh=list(loc=f$basis$mesh_vertices,
                                   graph=list(tv=f$basis$mesh_triangles)),
                         transform=transform),class="spde_mesh")
  raw <- sweep(f$basis$coordinates*transform$scale,2,transform$center,"+")
  basis <- spde_basis(mesh,raw,kappa=.7,pc_cutoff=.95)
  scaled <- rbind(c(.21,.41),c(.52,.18),c(.62,.82))
  raw_new <- sweep(scaled*transform$scale,2,transform$center,"+")
  A <- mgcvST:::.spde_basis_project(list(xy=basis$mesh_vertices,tv=basis$mesh_triangles),scaled)
  expect_equal(mgcvST:::.spde_basis_at(basis,raw_new),as.matrix(A %*% basis$projection),tolerance=1e-12)
  expect_equal(mgcvST:::.spde_basis_at(basis,raw_new,TRUE),as.matrix(A %*% basis$pc_mesh_projection),tolerance=1e-12)
})

test_that("prediction never rebuilds FEM, precision or eigensystems", {
  f <- st_fixture()
  testthat::local_mocked_bindings(.spde_basis_fem=function(...) stop("FEM recomputed"),
                                 .package="mgcvST")
  testthat::local_mocked_bindings(matrixEigen=function(...) stop("eigendecomposition recomputed"),
                                 matrixSolve=function(...) stop("precision recomputed"),
                                 .package="CppMatrix")
  new <- matrix(c(.33,.71),1,2)
  expect_true(is.matrix(mgcvST:::.spde_basis_at(f$basis,new)))
  expect_true(is.matrix(mgcvST:::.spde_basis_at(f$basis,new,TRUE)))
  project <- mgcvST:::.spde_basis_project
  calls <- 0L
  testthat::local_mocked_bindings(.spde_basis_project=function(...) {
    calls <<- calls + 1L
    project(...)
  }, .package="mgcvST")
  expect_identical(mgcvST:::.spde_basis_at(f$basis,f$basis$coordinates),f$basis$B)
  expect_equal(mgcvST:::.spde_basis_at(f$basis,f$basis$coordinates,TRUE),f$basis$pc_training_basis,tolerance=1e-12)
  expect_equal(mgcvST:::.spde_basis_at(f$basis,f$basis$coordinates[1:2,],TRUE),f$basis$pc_training_basis[1:2,],tolerance=1e-12)
  expect_identical(calls, 3L)
})

test_that("predict.gam works on new locations and prediction blocks", {
  for (pc in c(FALSE,TRUE)) {
    f <- st_fixture(pc=pc)
    fit <- mgcv::gam(G=f$G,method="REML")
    new <- data.frame(x=c(.132,.254,.815),y=c(.271,.724,.158),offset0=c(.1,0,-.1))
    L <- predict(fit,newdata=new,type="lpmatrix",block.size=1L)
    B <- mgcvST:::.spde_basis_at(f$basis,as.matrix(new[,c("x","y")]),pc)
    cols <- fit$smooth[[1]]$first.para:fit$smooth[[1]]$last.para
    expect_equal(unname(L[,cols,drop=FALSE]),unname(B),tolerance=1e-12)
    expect_equal(as.numeric(predict(fit,newdata=new,type="link",block.size=1L)),
                 as.numeric(L %*% coef(fit))+new$offset0,tolerance=1e-12)
    expect_equal(as.numeric(predict(fit,newdata=new,type="response")),
                 fit$family$linkinv(as.numeric(L %*% coef(fit))+new$offset0),tolerance=1e-12)
    training <- predict(fit,type="lpmatrix",block.size=11L)
    expect_equal(as.numeric(training[,cols,drop=FALSE]),
                     as.numeric(if(pc) f$basis$pc_training_basis else f$basis$B),tolerance=1e-12)
    new$x[1] <- -5
    expect_error(predict(fit,newdata=new,type="lpmatrix"),"outside")
  }
})

test_that("old prepared bases can be upgraded without recomputing PCs", {
  f <- st_fixture(pc=TRUE)
  legacy <- f$basis
  legacy$pc_cached_dimension <- legacy$pc_training_basis <- legacy$pc_mesh_projection <- NULL
  legacy$coordinate_keys <- NULL
  testthat::local_mocked_bindings(matrixEigen=function(...) stop("eigen"),.package="CppMatrix")
  expect_equal(mgcvST:::.spde_basis_at(legacy,legacy$coordinates,TRUE),f$basis$pc_training_basis,tolerance=1e-12)
  changed <- f$basis
  changed$pc_cutoff <- 1
  updated <- mgcvST:::.spde_basis_pc_cache(changed)
  expect_identical(ncol(updated$pc_training_basis),ncol(changed$B))
})

test_that("2D tsearchn and explicit quadtree share results, including edges", {
  f <- st_fixture()
  xy <- f$basis$mesh_vertices
  tv <- f$basis$mesh_triangles
  points <- rbind(xy,(xy[tv[,1],]+xy[tv[,2],])/2,c(.31,.68),c(-1,3))
  a <- geometry::tsearchn(xy,tv,points)
  b <- geometry::tsearch(xy[,1],xy[,2],tv,points[,1],points[,2],bary=TRUE,method="quadtree")
  expect_identical(a,b)
})
