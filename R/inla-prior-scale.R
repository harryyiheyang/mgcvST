# Normalize the spatial precision prior on the observation scale. This uses
# the already prepared constrained basis, so no subtraction of two enormous
# near-intercept covariance terms is needed for small kappa. Triangular solves
# are blocked over observations; no covariance inverse or n-by-n matrix forms.
.inlast_observation_precision_scale <- function(prepared, block_size = 256L) {
  B <- prepared$basis$B
  Q <- prepared$basis$Q
  R <- chol(Q)
  total <- 0
  n <- nrow(B)
  for (start in seq.int(1L, n, by = block_size)) {
    index <- seq.int(start, min(n, start + block_size - 1L))
    whitened <- forwardsolve(t(R), t(B[index, , drop = FALSE]))
    total <- total + sum(whitened * whitened)
  }
  scale <- total / n
  if (!is.finite(scale) || scale <= 0) {
    stop("The constrained SPDE has no positive observation-variance scale.")
  }
  as.numeric(scale)
}
