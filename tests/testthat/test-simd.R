test_that("the SIMD kernels are chosen at run time", {
  expect_true(spdeST:::mgcvst_simd_cpp() %in% c("generic", "avx2", "avx512"))
})

test_that("the full-spectrum pair kernel equals the single-pair reference", {
  set.seed(11)
  q <- 12L
  H <- lapply(1:4, function(g) {
    Z <- matrix(rnorm(q * (q - (g == 2) * 4L)), q)
    tcrossprod(Z) * 10^(g - 2)
  })
  a <- sapply(H, function(M) drop(M %*% rnorm(q)) / sqrt(max(abs(M))))
  G <- spdeST:::mgcvst_pair_basis_cpp(H, diag(q), 1L)
  pr <- t(combn(4L, 2L))
  z <- spdeST:::mgcvst_pair_spa_cpp(H, G, a, pr[, 1L], pr[, 2L], 1L, 4L)
  expect_true(all(z$status == 0L))
  expect_true(all(z$remainder_kind == 0L))
  for (t in seq_len(nrow(pr))) {
    ref <- rkhs_score_calibrate(z$score[t], H[[pr[t, 1L]]], H[[pr[t, 2L]]])
    expect_equal(z$log_p_two_sided[t], ref$log_p_two_sided, tolerance = 1e-8)
  }
})
