.score_scale_fixture <- function() {
  list(H = list(diag(c(0.2, 0.5, 1)), diag(c(0.7, 0.4, 0.8))),
       a = cbind(c(0.6, 0.3, 1), c(0.7, -0.1, 0.8)))
}

.score_pair_spa <- function(H, a) {
  q <- nrow(H[[1L]])
  G <- spdeST:::mgcvst_pair_basis_cpp(H, diag(q), 1L)
  spdeST:::mgcvst_pair_spa_cpp(H, G, a, 1L, 2L)
}

test_that("the fused exact path normalizes covariance and score together", {
  f <- .score_scale_fixture()
  ref <- .score_pair_spa(f$H, f$a)
  y <- .pair_ref_log_p(sum(f$a[, 1L] * f$a[, 2L]), f$H[[1L]], f$H[[2L]])
  scales <- rbind(c(1e-240, 1e-240), c(1e240, 1e240), c(1e-240, 1e240),
                  c(1e-20, 1e-20), c(1e20, 1e20))
  for (j in seq_len(nrow(scales))) {
    s <- scales[j, ]
    units <- sqrt(s[1L]) * sqrt(s[2L])
    H <- list(f$H[[1L]] * s[1L], f$H[[2L]] * s[2L])
    a <- sweep(f$a, 2L, sqrt(s), "*")
    x <- .score_pair_spa(H, a)
    expect_identical(x$status, 0L)
    expect_equal(x$score / units, ref$score, tolerance = 1e-12)
    expect_equal(x$log_p_two_sided, ref$log_p_two_sided, tolerance = 1e-10)
    expect_equal(x$log_p_two_sided, y, tolerance = 1e-9)
    G <- spdeST:::mgcvst_pair_basis_cpp(rev(H), diag(3L), 1L)
    swapped <- spdeST:::mgcvst_pair_spa_cpp(rev(H), G, a[, 2:1], 1L, 2L, 2L)
    expect_equal(swapped$log_p_two_sided, x$log_p_two_sided, tolerance = 1e-10)
    a[, 2L] <- -a[, 2L]
    negative <- .score_pair_spa(H, a)
    expect_equal(negative$score, -x$score, tolerance = 1e-12)
    expect_equal(negative$log_p_positive, x$log_p_negative, tolerance = 1e-10)
  }
  zero <- .score_pair_spa(list(matrix(0, 3L, 3L), f$H[[2L]]), f$a)
  expect_identical(zero$status, 1L)
  expect_true(is.na(zero$log_p_two_sided))
})

test_that("the seventeen historical positive-information failures are recovered", {
  fixture <- readRDS(test_path("fixtures", "score-calibration-history.rds"))
  expect_length(fixture, 17L)
  for (f in fixture) {
    expect_true(is.na(f$original_liu) && is.na(f$original_davies))
    L <- .pair_ref_log_p(f$U, f$H[[1L]], f$H[[2L]])
    C <- .score_pair_spa(f$H, f$a)
    # The pair kernel calibrates the pair that Liu moment matching could not,
    # and agrees with the full-spectrum reference.
    expect_true(is.finite(L))
    expect_identical(C$status, 0L)
    expect_equal(C$log_p_two_sided, L, tolerance = 1e-6)
    expect_equal(C$score / f$U, 1, tolerance = 1e-12)
  }
})
