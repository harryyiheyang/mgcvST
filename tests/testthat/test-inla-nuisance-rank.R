test_that("factorable singular nuisance penalties keep their algebraic rank", {
  Q <- Matrix::Diagonal(x = c(1, 1, 1, 1e-18, 1e-18))
  spec <- list(
    n = 8L, family = "gaussian",
    fixed = list(X = matrix(1, 8L, 1L), names = "intercept"),
    random = list(list(name = "nuisance_1",
      A = Matrix::Matrix(matrix(rnorm(40L), 8L, 5L), sparse = TRUE),
      Q = Q, kind = "nuisance", target = FALSE, rankdef = 2L,
      sp_index = 1L)), offset = rep(0, 8L)
  )
  expect_no_error(mgcvST:::.inlast_validate_spec(spec, rnorm(8L), rep(0, 8L)))
})

test_that("small full-rank nuisance penalties retain rankdef zero", {
  Q <- 1e-10 * Matrix::Diagonal(5L)
  spec <- list(
    n = 8L, family = "gaussian",
    fixed = list(X = matrix(1, 8L, 1L), names = "intercept"),
    random = list(list(name = "nuisance_1",
      A = Matrix::Matrix(matrix(rnorm(40L), 8L, 5L), sparse = TRUE),
      Q = Q, kind = "nuisance", target = FALSE, rankdef = 0L,
      sp_index = 1L)), offset = rep(0, 8L)
  )
  expect_no_error(mgcvST:::.inlast_validate_spec(spec, rnorm(8L), rep(0, 8L)))
})

test_that("zero nuisance penalties remain invalid", {
  Q <- Matrix::Diagonal(x = rep(0, 4L))
  spec <- list(
    n = 8L, family = "gaussian",
    fixed = list(X = matrix(1, 8L, 1L), names = "intercept"),
    random = list(list(name = "nuisance_1",
      A = Matrix::Matrix(matrix(rnorm(32L), 8L, 4L), sparse = TRUE),
      Q = Q, kind = "nuisance", target = FALSE, rankdef = 4L,
      sp_index = 1L)), offset = rep(0, 8L)
  )
  expect_error(mgcvST:::.inlast_validate_spec(spec, rnorm(8L), rep(0, 8L)),
    "nonzero positive-semidefinite")
})

test_that("indefinite nuisance penalties remain invalid", {
  Q <- Matrix::Diagonal(x = c(1, 1, -1, 0))
  spec <- list(
    n = 8L, family = "gaussian",
    fixed = list(X = matrix(1, 8L, 1L), names = "intercept"),
    random = list(list(name = "nuisance_1",
      A = Matrix::Matrix(matrix(rnorm(32L), 8L, 4L), sparse = TRUE),
      Q = Q, kind = "nuisance", target = FALSE, rankdef = 1L,
      sp_index = 1L)), offset = rep(0, 8L)
  )
  expect_error(mgcvST:::.inlast_validate_spec(spec, rnorm(8L), rep(0, 8L)),
    "positive-semidefinite")
})
