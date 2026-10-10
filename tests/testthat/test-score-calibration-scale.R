.score_scale_fixture <- function() {
  list(H = list(diag(c(0.2, 0.5, 1)), diag(c(0.7, 0.4, 0.8))),
       a = cbind(c(0.6, 0.3, 1), c(0.7, -0.1, 0.8)))
}

.score_pair_spa <- function(H, a) {
  q <- nrow(H[[1L]])
  G <- mgcvST:::mgcvst_pair_basis_cpp(H, diag(q), 1L)
  mgcvST:::mgcvst_pair_spa_cpp(H, G, a, 1L, 2L)
}

test_that("R saddlepoint calibration is invariant to covariance units", {
  f <- .score_scale_fixture()
  U <- sum(f$a[, 1L] * f$a[, 2L])
  scales <- rbind(c(1e-240, 1e-240), c(1e240, 1e240),
                  c(1e-240, 1e240), c(1e-20, 1e-20), c(1e20, 1e20))
  ref <- rkhs_score_calibrate(U, f$H[[1L]], f$H[[2L]])
  for (j in seq_len(nrow(scales))) {
    s <- scales[j, ]
    units <- sqrt(s[1L]) * sqrt(s[2L])
    x <- rkhs_score_calibrate(U * units, f$H[[1L]] * s[1L], f$H[[2L]] * s[2L])
    expect_true(is.finite(x$p_two_sided))
    expect_equal(x$p_two_sided, ref$p_two_sided, tolerance = 1e-10)
    expect_equal(x$log_p_two_sided, ref$log_p_two_sided, tolerance = 1e-10)
    expect_equal(x$effective_rank, ref$effective_rank, tolerance = 1e-12)
    expect_equal(x$information, (ref$information * units) * units, tolerance = 1e-12)
  }
  negative <- rkhs_score_calibrate(-U, f$H[[1L]], f$H[[2L]])
  swapped <- rkhs_score_calibrate(U, f$H[[2L]], f$H[[1L]])
  expect_equal(negative$p_two_sided, ref$p_two_sided)
  expect_equal(negative$p_positive, ref$p_negative)
  expect_equal(swapped$p_two_sided, ref$p_two_sided, tolerance = 1e-10)
  expect_equal(ref$p_two_sided, min(1, 2 * min(ref$p_positive, ref$p_negative)))
  expect_false("singular_values" %in% names(ref))
  expect_false("liu_parameters" %in% names(ref))
  expect_false("method" %in% names(formals(rkhs_score_calibrate)))
  # The calibration is the saddlepoint of the full singular spectrum.
  expect_identical(ref$spa$k, 3L)
  expect_identical(ref$spa$remainder_kind, 0)
  s <- mgcvST:::.rkhs_score_spectrum(f$H[[1L]], f$H[[2L]])
  expect_equal(ref$log_p_two_sided, .spa_ref(abs(U), s), tolerance = 1e-8)
})

test_that("public moments, cumulants and singular values retain input units", {
  f <- .score_scale_fixture()
  U <- sum(f$a[, 1L] * f$a[, 2L])
  base <- mgcvST:::.rkhs_score_spectrum(f$H[[1L]], f$H[[2L]])
  ref <- rkhs_score_calibrate(U, f$H[[1L]], f$H[[2L]])
  for (s in c(1e-20, 1, 1e20)) {
    H <- lapply(f$H, function(z) z * s)
    moments <- mgcvST:::.rkhs_score_moments(H[[1L]], H[[2L]])
    x <- rkhs_score_calibrate(U * s, H[[1L]], H[[2L]])
    expect_equal(x$moments / moments, rep(1, 4L), tolerance = 1e-12)
    expect_equal(x$information / moments[1L], 1, tolerance = 1e-12)
    expect_equal(x$p_two_sided, ref$p_two_sided, tolerance = 1e-10)
    expect_equal(mgcvST:::.rkhs_score_spectrum(H[[1L]], H[[2L]]) / base,
                 rep(s, length(base)), tolerance = 1e-12)
  }
  # The trace moments are those of the singular spectrum.
  expect_equal(ref$moments, vapply(1:4, function(r) sum(base^(2 * r)), numeric(1L)),
               tolerance = 1e-12)
})

test_that("the PSD tolerance and zero-state contract are unchanged", {
  accepted <- diag(c(-1e-11, 1e-5))
  other <- diag(c(1e-4, 1e-4))
  expect_silent(mgcvST:::.psd_factor(accepted))
  expect_true(length(mgcvST:::.rkhs_score_spectrum(accepted, other)) > 0L)
  rejected <- diag(c(-1e-9, 1e-5))
  expect_error(mgcvST:::.rkhs_score_spectrum(rejected, other), "not positive semidefinite")
  zero <- rkhs_score_calibrate(0, matrix(0, 2L, 2L), other)
  expect_true(is.na(zero$p_two_sided))
  expect_equal(zero$information, 0)
  expect_equal(zero$moments, rep(0, 4L))
  expect_null(zero$spa)
})

test_that("the fused exact path normalizes covariance and score together", {
  f <- .score_scale_fixture()
  ref <- .score_pair_spa(f$H, f$a)
  scales <- rbind(c(1e-240, 1e-240), c(1e240, 1e240), c(1e-240, 1e240),
                  c(1e-20, 1e-20), c(1e20, 1e20))
  for (j in seq_len(nrow(scales))) {
    s <- scales[j, ]
    units <- sqrt(s[1L]) * sqrt(s[2L])
    H <- list(f$H[[1L]] * s[1L], f$H[[2L]] * s[2L])
    a <- sweep(f$a, 2L, sqrt(s), "*")
    x <- .score_pair_spa(H, a)
    y <- rkhs_score_calibrate(sum(a[, 1L] * a[, 2L]), H[[1L]], H[[2L]])
    expect_identical(x$status, 0L)
    expect_equal(x$score / units, ref$score, tolerance = 1e-12)
    expect_equal(x$log_p_two_sided, ref$log_p_two_sided, tolerance = 1e-10)
    expect_equal(x$log_p_two_sided, y$log_p_two_sided, tolerance = 1e-9)
    G <- mgcvST:::mgcvst_pair_basis_cpp(rev(H), diag(3L), 1L)
    swapped <- mgcvST:::mgcvst_pair_spa_cpp(rev(H), G, a[, 2:1], 1L, 2L, 2L)
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
    L <- rkhs_score_calibrate(f$U, f$H[[1L]], f$H[[2L]])
    C <- .score_pair_spa(f$H, f$a)
    # Both paths calibrate the pair that Liu moment matching could not, and
    # they agree.
    expect_true(is.finite(L$log_p_two_sided))
    expect_identical(C$status, 0L)
    expect_equal(C$log_p_two_sided, L$log_p_two_sided, tolerance = 1e-6)
    expect_equal(L$information / f$information, 1, tolerance = 1e-12)
    expect_equal(C$score / f$U, 1, tolerance = 1e-12)
  }
})
