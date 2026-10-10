# Benchmark: fused C++ Liu pair kernel (mgcvst_pair_liu_cpp) vs the old
# mgcvst_pair_trace_powers_cpp + R Liu path. Not run under testthat; run
# interactively with
#   Rscript tests/benchmark/bench-liu-exact.R

suppressPackageStartupMessages({
  library(mgcvST)
})

cat("== Part 1: fused Liu pair kernel vs old trace-powers + R Liu ==\n")

bench_kernel <- function(q, K, threads = 4L, seed = 1L) {
  set.seed(seed)
  H <- lapply(seq_len(K), function(k) {
    z <- matrix(rnorm(q * q), q, q)
    crossprod(z) + diag(q) * 0.1
  })
  a <- matrix(rnorm(q * K), q, K)
  idx <- which(upper.tri(matrix(0, K, K)), arr.ind = TRUE)
  left <- idx[, 2L]
  right <- idx[, 1L]
  ord <- order(left)
  left <- left[ord]
  right <- right[ord]
  pairs <- cbind(left, right)

  t_old <- system.time({
    old_score <- colSums(a[, left, drop = FALSE] * a[, right, drop = FALSE])
    old_moments <- mgcvST:::mgcvst_pair_trace_powers_cpp(
      H, pairs, maxPower = 4L, threads = threads
    )
    old_liu <- mgcvST:::.liu_squared_score_moments(
      abs(old_score), old_moments[, 1L], old_moments[, 2L],
      old_moments[, 3L], old_moments[, 4L]
    )
  })["elapsed"]

  t_new <- system.time({
    new <- mgcvST:::mgcvst_pair_liu_cpp(H, a, left, right, threads = threads)
  })["elapsed"]

  data.frame(q = q, K = K, pairs = nrow(pairs),
             old_seconds = as.numeric(t_old), new_seconds = as.numeric(t_new),
             speedup = as.numeric(t_old) / as.numeric(t_new))
}

kernel_table <- do.call(rbind, lapply(c(30L, 60L, 100L), function(q) {
  bench_kernel(q, K = 400L, threads = 4L)
}))
print(kernel_table, row.names = FALSE)
