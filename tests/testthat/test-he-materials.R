he_states <- function(seed, q = 12L, n = 4L) {
  set.seed(seed)
  H <- lapply(seq_len(n), function(i) {
    X <- matrix(stats::rnorm(q * q), q)
    crossprod(X) / q
  })
  list(H = H, a = matrix(stats::rnorm(q * n), q))
}

test_that("the exact pair kernel returns t1 = tr(M_i M_j) in public units", {
  s <- he_states(9049L)
  q <- nrow(s$a)
  V <- qr.Q(qr(matrix(stats::rnorm(q * 5L), q)))
  G <- spdeST:::mgcvst_pair_basis_cpp(s$H, V, 1L)
  res <- spdeST:::mgcvst_pair_spa_cpp(s$H, G, s$a, c(1L, 1L, 2L), c(2L, 3L, 4L),
                                      1L, 4L)
  ref <- c(sum(s$H[[1]] * s$H[[2]]), sum(s$H[[1]] * s$H[[3]]),
           sum(s$H[[2]] * s$H[[4]]))
  expect_equal(res$t1, ref, tolerance = 1e-12)
  expect_equal(res$score, colSums(s$a[, c(1, 1, 2)] * s$a[, c(2, 3, 4)]))
  # The full-spectrum path (k = q) reports the same t1.
  Gq <- spdeST:::mgcvst_pair_basis_cpp(s$H, diag(q), 1L)
  full <- spdeST:::mgcvst_pair_spa_cpp(s$H, Gq, s$a, c(1L, 1L, 2L),
                                       c(2L, 3L, 4L), 1L, 4L)
  expect_equal(full$t1, ref, tolerance = 1e-12)
})

test_that("the PCAlearning materialization reports tr(M) and tr(M^2)", {
  s <- he_states(9050L)
  z <- spdeST:::mgcvst_pca_dense_cpp(s$H, s$a, NULL, rep(FALSE, 4L), NULL, 1L)
  expect_equal(z$trace, vapply(s$H, function(M) sum(diag(M)), 0), tolerance = 1e-12)
  expect_equal(z$fro2, vapply(s$H, function(M) sum(M * M), 0), tolerance = 1e-10)
})

test_that("pair results carry t1 and the test keeps per-feature traces", {
  skip_if_not_installed("arrow")
  .pair_pipeline_mocks()
  fit <- .pair_pipeline_fit(6L)
  exact <- spdeST:::.mgcvst_pair_pipeline(fit, NULL, 1L, 100L, FALSE, k = 2L)
  pca <- spdeST:::.mgcvst_pair_pcalearning(fit, NULL, 1L, 100L, FALSE,
                                           rank = 2L, k = 2L)
  a <- do.call(rbind, lapply(exact$shards, spdeST:::.mgcvst_read_shard))
  a <- a[order(a$i, a$j), ]
  b <- .pca_pairs(pca)
  expect_true("t1" %in% names(a))
  expect_true(all(is.finite(a$t1) & a$t1 > 0))
  # The PCAlearning t1 is c_i' c_j of the projected matrices.
  expect_equal(b$t1, a$t1, tolerance = 0.05)
  expect_identical(exact$features$feature_index, pca$features$feature_index)
  expect_equal(exact$features$trace_M, pca$features$trace_M, tolerance = 1e-10)
  expect_equal(exact$features$trace_M2, pca$features$trace_M2, tolerance = 1e-8)
  expect_identical(spdeST:::.mgcvst_contract("exact")$schema, "compact_v2")
})
