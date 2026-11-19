# Multiply dense matrices through the package's single-thread numerical backend.
.magic_mm <- function(A, B, transA = FALSE, transB = FALSE) {
  CppMatrix::matrixMultiply(
    as.matrix(A), as.matrix(B), transA = transA, transB = transB
  )
}

# Normalize and validate a finite numeric matrix input.
.as_numeric_matrix <- function(x, name) {
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  if (length(dim(x)) != 2L || any(!is.finite(x))) {
    stop(name, " must be a finite numeric matrix.")
  }
  x
}
