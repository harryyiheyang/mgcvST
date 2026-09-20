test_that("nuisance rankdef travels with the block spec", {
  Q <- Matrix::Diagonal(x = c(1, 1, 1, 1e-18, 1e-18))
  spec <- list(
    n = 8L, family = "gaussian",
    fixed = list(X = matrix(1, 8L, 1L), names = "intercept"),
    random = list(list(name = "nuisance_1",
      A = Matrix::Matrix(matrix(rnorm(40L), 8L, 5L), sparse = TRUE),
      Q = Q, kind = "nuisance", target = FALSE, rankdef = 2L,
      sp_index = 1L)), offset = rep(0, 8L)
  )
  z <- mgcvST:::.inlast_validate_spec(spec, rnorm(8L), rep(0, 8L))
  expect_identical(z$random[[1L]]$rankdef, 2L)
})

test_that("an unsupplied rankdef defaults to zero", {
  Q <- 1e-10 * Matrix::Diagonal(5L)
  spec <- list(
    n = 8L, family = "gaussian",
    fixed = list(X = matrix(1, 8L, 1L), names = "intercept"),
    random = list(list(name = "nuisance_1",
      A = Matrix::Matrix(matrix(rnorm(40L), 8L, 5L), sparse = TRUE),
      Q = Q, kind = "nuisance", target = FALSE,
      sp_index = 1L)), offset = rep(0, 8L)
  )
  z <- mgcvST:::.inlast_validate_spec(spec, rnorm(8L), rep(0, 8L))
  expect_identical(z$random[[1L]]$rankdef, 0L)
})
