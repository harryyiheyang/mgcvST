# Multiply dense matrices through the package's single-thread numerical backend.
.magic_mm <- function(A, B, transA = FALSE, transB = FALSE) {
  CppMatrix::matrixMultiply(
    as.matrix(A), as.matrix(B), transA = transA, transB = transB
  )
}

