test_that("Stage 1 p-values are finite", {
  f <- st_fixture()
  fit <- spdeST:::mgcvST.estimate_spatial(
    mgcvST.estimate(f$Y, f$model),
    f$Y, features = "all")
  expect_true(all(is.finite(fit$diagnostics$marginal_p_value)))
})

test_that("Snow workers give the pair results of serial chunks", {
  skip_on_cran()
  f <- st_fixture()
  bp <- BiocParallel::SnowParam(workers=2,type="SOCK",progressbar=FALSE)
  on.exit(BiocParallel::bpstop(bp),add=TRUE)
  Y <- f$Y[rep(1:3, 2),,drop=FALSE]
  rownames(Y) <- paste0("snow", seq_len(nrow(Y)))
  fit <- spdeST:::mgcvST.estimate_spatial(
    mgcvST.estimate(Y,f$model,BPPARAM=bp,chunk_size=1),
    Y, features = "all", BPPARAM=bp, chunk_size=1)
  expect_true(all(vapply(fit$nuisance_covariance, is.matrix, logical(1L))))
  pairs <- t(combn(1:3,2))
  serial <- mgcvST.test(fit,pairs=pairs, moments = "exact")
  blocks <- mgcvST.test(fit,pairs=pairs,chunk_size=1, moments = "exact")
  expect_identical(serial$results,blocks$results)
})
