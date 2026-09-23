test_that("batched landmark traces agree with exact noncommuting pairs", {
  genes <- list(
    diag(c(1, 2, 3)),
    crossprod(matrix(c(1, 2, 0, 0, 1, 1, 2, 0, 1), 3, 3)) + diag(3) / 4,
    tcrossprod(matrix(c(1, 0, 2, 1, 2, 1, 0, 1, 1), 3, 3)) + diag(3) / 3
  )
  references <- list(
    crossprod(matrix(c(1, 1, 0, 0, 2, 1, 1, 0, 2), 3, 3)) + diag(3) / 2,
    diag(c(2, 1, 4))
  )
  expect_gt(max(abs(genes[[1L]] %*% references[[1L]] -
                    references[[1L]] %*% genes[[1L]])), 0)

  pairs <- cbind(
    rep(seq_along(genes), each = length(references)),
    length(genes) + rep(seq_along(references), times = length(genes))
  )
  expected <- mgcvST:::mgcvst_pair_trace_powers_cpp(
    c(genes, references), pairs, 4L, 1L
  )
  double <- mgcvST:::mgcvst_landmark_trace_cpp(genes, references, 1L, FALSE)
  float <- mgcvST:::mgcvst_landmark_trace_cpp(genes, references, 1L, TRUE)

  expect_identical(dim(double), c(length(genes) * length(references), 4L))
  expect_equal(double, expected, tolerance = 1e-12)
  expect_equal(float, double, tolerance = 1e-5)
  expect_identical(
    mgcvST:::mgcvst_landmark_trace_cpp(genes, references, 2L, FALSE), double
  )
  expect_identical(
    mgcvST:::mgcvst_landmark_trace_cpp(genes, references, 2L, TRUE), float
  )
})

test_that("batched landmark traces validate dimensions and thread count", {
  M <- diag(2)
  expect_error(mgcvST:::mgcvst_landmark_trace_cpp(list(), list(M)),
               "non-empty")
  expect_error(mgcvST:::mgcvst_landmark_trace_cpp(list(M), list()),
               "non-empty")
  expect_error(mgcvST:::mgcvst_landmark_trace_cpp(list(M), list(diag(3))),
               "same square dimension")
  expect_error(mgcvST:::mgcvst_landmark_trace_cpp(list(M), list(M), 0L),
               "positive integer")
})
