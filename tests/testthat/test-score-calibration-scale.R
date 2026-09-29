.score_scale_fixture <- function() {
  list(H = list(diag(c(0.2, 0.5, 1)), diag(c(0.7, 0.4, 0.8))),
       a = cbind(c(0.6, 0.3, 1), c(0.7, -0.1, 0.8)))
}

test_that("R Liu and Davies calibration are invariant to covariance units", {
  skip_if_not_installed("CompQuadForm")
  f <- .score_scale_fixture()
  U <- sum(f$a[, 1L] * f$a[, 2L])
  scales <- rbind(c(1e-240, 1e-240), c(1e240, 1e240),
                  c(1e-240, 1e240), c(1e-20, 1e-20), c(1e20, 1e20))
  for (method in c("liu", "davies")) {
    ref <- rkhs_score_calibrate(U, f$H[[1L]], f$H[[2L]], method)
    if (method == "davies") {
      s <- rkhs_score_singular_values(f$H[[1L]], f$H[[2L]])
      raw <- CompQuadForm::davies(abs(U), lambda = c(s / 2, -s / 2))
      expect_equal(raw$ifault, 0)
      expect_lt(abs(ref$p_two_sided - min(1, 2 * raw$Qq)), 2e-4)
    }
    for (j in seq_len(nrow(scales))) {
      s <- scales[j, ]
      units <- sqrt(s[1L]) * sqrt(s[2L])
      x <- rkhs_score_calibrate(U * units, f$H[[1L]] * s[1L], f$H[[2L]] * s[2L], method)
      expect_true(is.finite(x$p_two_sided))
      expect_equal(x$p_two_sided, ref$p_two_sided, tolerance = 1e-10)
      expect_equal(x$effective_rank, ref$effective_rank, tolerance = 1e-12)
      expect_equal(x$information, (ref$information * units) * units, tolerance = 1e-12)
    }
    negative <- rkhs_score_calibrate(-U, f$H[[1L]], f$H[[2L]], method)
    swapped <- rkhs_score_calibrate(U, f$H[[2L]], f$H[[1L]], method)
    expect_equal(negative$p_two_sided, ref$p_two_sided)
    expect_equal(negative$p_positive, ref$p_negative)
    expect_equal(swapped$p_two_sided, ref$p_two_sided, tolerance = 1e-10)
    expect_equal(ref$p_two_sided, min(1, 2 * min(ref$p_positive, ref$p_negative)))
  }
})

test_that("public moments, cumulants and singular values retain input units", {
  f <- .score_scale_fixture()
  U <- sum(f$a[, 1L] * f$a[, 2L])
  for (s in c(1e-20, 1, 1e20)) {
    H <- lapply(f$H, function(z) z * s)
    moments <- mgcvST:::.rkhs_score_moments(H[[1L]], H[[2L]])
    old <- mgcvST:::.liu_squared_score_moments(U * s, moments[1L], moments[2L], moments[3L], moments[4L])
    x <- rkhs_score_calibrate(U * s, H[[1L]], H[[2L]], "liu")
    expect_equal(x$moments / moments, rep(1, 4L), tolerance = 1e-12)
    expect_equal(x$information / moments[1L], 1, tolerance = 1e-12)
    expect_equal(unlist(x$liu_parameters[paste0("c", 1:4)]) /
                   unlist(old[paste0("c", 1:4)]), rep(1, 4L), ignore_attr = TRUE,
                 tolerance = 1e-12)
    expect_equal(x$p_two_sided, old$p_value, tolerance = 1e-12)
    y <- rkhs_score_calibrate(U * s, H[[1L]], H[[2L]], "davies")
    expected <- rkhs_score_singular_values(H[[1L]], H[[2L]])
    expect_equal(y$singular_values / expected, rep(1, length(expected)), tolerance = 1e-12)
  }
})

test_that("Davies retains the original PSD tolerance and zero-state contract", {
  accepted <- diag(c(-1e-11, 1e-5))
  other <- diag(c(1e-4, 1e-4))
  expect_silent(mgcvST:::.psd_factor(accepted))
  expect_true(is.finite(rkhs_score_calibrate(1e-6, accepted, other, "davies")$p_two_sided))
  rejected <- diag(c(-1e-9, 1e-5))
  expect_error(rkhs_score_calibrate(1e-6, rejected, other, "davies"), "not positive semidefinite")
  for (method in c("liu", "davies")) {
    zero <- rkhs_score_calibrate(0, matrix(0, 2L, 2L), other, method)
    expect_true(is.na(zero$p_two_sided))
    expect_equal(zero$information, 0)
    expect_equal(zero$moments, rep(0, 4L))
  }
})

test_that("the fused Liu path normalizes covariance and score together", {
  f <- .score_scale_fixture()
  ref <- mgcvST:::mgcvst_pair_liu_cpp(f$H, f$a, 1L, 2L)
  scales <- rbind(c(1e-240, 1e-240), c(1e240, 1e240), c(1e-240, 1e240),
                  c(1e-20, 1e-20), c(1e20, 1e20))
  for (j in seq_len(nrow(scales))) {
    s <- scales[j, ]
    units <- sqrt(s[1L]) * sqrt(s[2L])
    H <- list(f$H[[1L]] * s[1L], f$H[[2L]] * s[2L])
    a <- sweep(f$a, 2L, sqrt(s), "*")
    x <- mgcvST:::mgcvst_pair_liu_cpp(H, a, 1L, 2L)
    y <- rkhs_score_calibrate(sum(a[, 1L] * a[, 2L]), H[[1L]], H[[2L]], "liu")
    expect_identical(x$status, 0L)
    expect_equal(x$score / units, ref$score, tolerance = 1e-12)
    expect_equal(x$information, (ref$information * units) * units, tolerance = 1e-12)
    expect_equal(x$effective_rank, ref$effective_rank, tolerance = 1e-12)
    expect_equal(x$log_p_two_sided, ref$log_p_two_sided, tolerance = 1e-12)
    expect_equal(exp(x$log_p_two_sided), y$p_two_sided, tolerance = 1e-12)
    swapped <- mgcvST:::mgcvst_pair_liu_cpp(rev(H), a[, 2:1], 1L, 2L, 2L)
    expect_equal(swapped, x, tolerance = 1e-12)
    a[, 2L] <- -a[, 2L]
    negative <- mgcvST:::mgcvst_pair_liu_cpp(H, a, 1L, 2L)
    expect_equal(negative$score, -x$score, tolerance = 1e-12)
    expect_equal(negative$log_p_positive, x$log_p_negative, tolerance = 1e-12)
  }
  zero <- mgcvST:::mgcvst_pair_liu_cpp(list(matrix(0, 3L, 3L), f$H[[2L]]), f$a, 1L, 2L)
  expect_identical(zero$status, 1L)
  expect_true(is.na(zero$log_p_two_sided))
})

test_that("the seventeen historical positive-information failures are recovered", {
  skip_if_not_installed("CompQuadForm")
  fixture <- readRDS(test_path("fixtures", "score-calibration-history.rds"))
  expect_length(fixture, 17L)
  for (f in fixture) {
    expect_true(is.na(f$original_liu) && is.na(f$original_davies))
    L <- rkhs_score_calibrate(f$U, f$H[[1L]], f$H[[2L]], "liu")
    D <- rkhs_score_calibrate(f$U, f$H[[1L]], f$H[[2L]], "davies")
    C <- mgcvST:::mgcvst_pair_liu_cpp(f$H, f$a, 1L, 2L)
    expect_equal(L$p_two_sided, f$reference_liu, tolerance = 1e-12)
    expect_equal(D$p_two_sided, f$reference_davies, tolerance = 1e-8)
    expect_identical(C$status, 0L)
    expect_equal(exp(C$log_p_two_sided), f$reference_liu, tolerance = 1e-12)
    expect_equal(L$information / f$information, 1, tolerance = 1e-12)
    expect_equal(C$information / f$information, 1, tolerance = 1e-12)
    expect_equal(C$score / f$U, 1, tolerance = 1e-12)
  }
})

test_that("normalized calibration does not resume old pair results", {
  root <- tempfile("score-scale-checkpoint-")
  dir.create(root)
  store <- list(path = root, temporary = FALSE)
  index <- matrix(1:2, 1L)
  result <- data.frame(pair_index = 1L, score = 1e-6, information = 1e-11,
    effective_rank = 0, p_value = NA_real_, log_p_two_sided = NA_real_,
    log_p_positive = NA_real_, log_p_negative = NA_real_, error_message = NA_character_)
  for (calibration in list(NULL, "davies")) {
    inputs <- list(version = 1L, index = index, pair_index = 1L)
    if (!is.null(calibration)) inputs$calibration <- calibration
    old <- file.path(root, paste0("pairs-", digest::digest(inputs, algo = "sha256")))
    dir.create(old)
    mgcvST:::.mgcvst_pair_checkpoint_write(old, 1L, 1L, result)
    path <- mgcvST:::.mgcvst_pair_checkpoint(store, index, 1L, calibration)
    expect_false(identical(path, old))
    expect_null(mgcvST:::.mgcvst_pair_checkpoint_read(path, 1L, 1L))
    repaired <- result
    repaired$p_value <- 0.5
    repaired$log_p_two_sided <- log(0.5)
    repaired$log_p_positive <- log(0.25)
    repaired$log_p_negative <- log(0.75)
    mgcvST:::.mgcvst_pair_checkpoint_write(path, 1L, 1L, repaired)
    expect_identical(mgcvST:::.mgcvst_pair_checkpoint(store, index, 1L, calibration), path)
    expect_identical(mgcvST:::.mgcvst_pair_checkpoint_read(path, 1L, 1L)$result, repaired)
    expect_identical(mgcvST:::.mgcvst_pair_checkpoint_read(old, 1L, 1L)$result, result)
  }
})
