# Pack one symmetric score matrix without retaining a duplicated triangle.
.mgcvst_pack_symmetric <- function(M) {
  M <- as.matrix(M)
  if (nrow(M) != ncol(M)) stop("A score-state matrix must be square.")
  list(n = nrow(M), upper = M[upper.tri(M, diag = TRUE)])
}

.mgcvst_unpack_symmetric <- function(x) {
  n <- as.integer(x$n)
  M <- matrix(0, n, n)
  keep <- upper.tri(M, diag = TRUE)
  if (length(x$upper) != sum(keep)) stop("Packed score state has an invalid size.")
  M[keep] <- x$upper
  M[lower.tri(M)] <- t(M)[lower.tri(M)]
  M
}

.mgcvst_pack_score_state <- function(state) {
  list(a = as.numeric(state$a), M = .mgcvst_pack_symmetric(state$M),
       width = state$width)
}

.mgcvst_unpack_score_state <- function(unit) {
  list(a = unit$a, M = .mgcvst_unpack_symmetric(unit$M), width = unit$width)
}

# Create and remove only a verified child of the R session temporary directory.
.mgcvst_dense_temp_dir <- function() {
  path <- tempfile("mgcvst-dense-cache-", tmpdir = tempdir())
  if (!dir.create(path)) stop("Could not create the dense score cache directory.")
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

.mgcvst_dense_cleanup <- function(path) {
  root <- paste0(normalizePath(tempdir(), winslash = "/", mustWork = TRUE), "/")
  target <- paste0(normalizePath(path, winslash = "/", mustWork = TRUE), "/")
  if (!startsWith(target, root) || identical(target, root)) {
    stop("Refusing to remove an unverified dense score cache directory.")
  }
  unlink(sub("/$", "", target), recursive = TRUE, force = TRUE)
}
