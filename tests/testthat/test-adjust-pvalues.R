.ref_step_up <- function(lp, by) {
  m <- length(lp)
  o <- order(lp)
  k <- seq_len(m)
  raw <- lp[o] + log(m) + (if (by) log(sum(1 / k)) else 0) - log(k)
  q <- numeric(m)
  q[o] <- pmin(0, rev(cummin(rev(raw))))
  q
}

test_that("log-space BH and BY match stats::p.adjust within the double range", {
  withr::local_seed(9L)
  for (n in c(1L, 2L, 7L, 500L, 20000L)) {
    p <- runif(n)^3
    lp <- log(p)
    for (method in c("BH", "BY")) {
      lq <- mgcvST:::.mgcvst_log_adjust(lp, method)
      expect_equal(lq$n, n)
      expect_equal(exp(lq$log_q), p.adjust(p, method), tolerance = 1e-12)
      expect_identical(lq$log_q <= log(0.05), p.adjust(p, method) <= 0.05)
    }
  }
  # Ties receive the same adjusted value, as in p.adjust.
  p <- c(0.01, 0.01, 0.01, 0.2, 0.2, 0.9)
  for (method in c("BH", "BY")) {
    expect_equal(exp(mgcvST:::.mgcvst_log_adjust(log(p), method)$log_q),
                 p.adjust(p, method), tolerance = 1e-14)
  }
})

test_that("BY uses c(m) = sum(1 / i) and BH does not", {
  lp <- log(c(1e-4, 1e-3, 0.02, 0.3))
  by <- mgcvST:::.mgcvst_log_adjust(lp, "BY")$log_q
  bh <- mgcvST:::.mgcvst_log_adjust(lp, "BH")$log_q
  expect_equal(exp(by[1L]) / exp(bh[1L]), sum(1 / 1:4), tolerance = 1e-12)
  expect_equal(exp(by), p.adjust(exp(lp), "BY"), tolerance = 1e-12)
  expect_equal(exp(bh), p.adjust(exp(lp), "BH"), tolerance = 1e-12)
})

test_that("single-step Sidak matches 1 - (1 - p)^m and keeps tiny p-values", {
  withr::local_seed(10L)
  p <- runif(300)^4
  lq <- mgcvST:::.mgcvst_log_adjust(log(p), "Sidak")$log_q
  expect_equal(exp(lq), 1 - (1 - p)^length(p), tolerance = 1e-12)
  deep <- mgcvST:::.mgcvst_log_adjust(c(-5000, -1e-3, -2), "Sidak")$log_q
  expect_equal(deep[1L], -5000 + log(3), tolerance = 1e-12)
  expect_true(all(deep <= 0))
  expect_equal(mgcvST:::.mgcvst_log_adjust(c(0, -Inf), "Sidak")$log_q, c(0, -Inf))
})

test_that("none returns the log p-values unchanged", {
  lp <- c(-5000, -3, -1e-9, NA, 0)
  expect_identical(mgcvST:::.mgcvst_log_adjust(lp, "none")$log_q, lp)
})

test_that("log p-values down to exp(-5000) keep their order and adjusted values", {
  lp <- c(-5000, -3000, -1, -800, -2, -4999)
  for (method in c("BH", "BY")) {
    lq <- mgcvST:::.mgcvst_log_adjust(lp, method)$log_q
    expect_true(all(is.finite(lq)))
    expect_equal(lq, .ref_step_up(lp, method == "BY"), tolerance = 1e-12)
    # Adjustment is monotone in the p-value.
    expect_false(is.unsorted(lq[order(lp)]))
  }
  expect_equal(mgcvST:::.mgcvst_log_adjust(c(-5000, -3000, -1), "BH")$log_q,
               c(-5000 + log(3), -3000 + log(1.5), -1), tolerance = 1e-12)
  # Decisions are defined where the p-values underflow to zero.
  expect_true(all(exp(lp[c(1L, 2L, 4L)]) == 0))
  expect_true(mgcvST:::.mgcvst_log_adjust(lp, "BY")$log_q[1L] <= log(1e-300))
  big <- sort(c(-5000, -4000, -3000, -2000, -1000, seq(-50, -1, length.out = 45)))
  by <- mgcvST:::.mgcvst_log_adjust(big, "BY")$log_q
  expect_equal(by, .ref_step_up(big, TRUE), tolerance = 1e-12)
  expect_true(all(by[1:5] < log(1e-100)))
})

test_that("zero p-values, rounding above zero and missing values are handled", {
  lq <- mgcvST:::.mgcvst_log_adjust(c(-Inf, -2, 1e-12, NA, NaN, Inf, -1), "BH")
  expect_equal(lq$n, 4)
  expect_identical(is.na(lq$log_q), c(FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE))
  expect_identical(lq$log_q[1L], -Inf)
  expect_true(all(lq$log_q[!is.na(lq$log_q)] <= 0))
  # m counts only the usable p-values: the same as dropping the missing ones.
  expect_equal(lq$log_q[c(1L, 2L, 3L, 7L)],
               mgcvST:::.mgcvst_log_adjust(c(-Inf, -2, 0, -1), "BH")$log_q)
  empty <- mgcvST:::.mgcvst_log_adjust(numeric(), "BY")
  expect_length(empty$log_q, 0L)
  allna <- mgcvST:::.mgcvst_log_adjust(c(NA, NaN), "BY")
  expect_true(all(is.na(allna$log_q)) && allna$n == 0)
  expect_error(mgcvST:::mgcvst_log_adjust_cpp(0, "holm"), "BY")
  expect_error(mgcvST:::.mgcvst_log_adjust(0, "holm"), "should be one of")
})

test_that("compact result shards carry one adjustment and the sign split", {
  withr::local_seed(11L)
  dir <- tempfile("mgcvst-finalize-")
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  contract <- mgcvST:::.mgcvst_contract("exact")
  pair_dir <- mgcvST:::.mgcvst_pairs_open(dir, list(test = 1L), contract)
  make <- function(i, j, score, p2) {
    mgcvST:::.mgcvst_pairs_frame(i, j, score, log(p2),
      log_p_positive = ifelse(score > 0, log(p2 / 2), log1p(-p2 / 2)),
      log_p_negative = ifelse(score < 0, log(p2 / 2), log1p(-p2 / 2)))
  }
  a <- make(c(1L, 1L, 1L), c(2L, 3L, 4L), c(2.5, -3, 0.1), c(1e-8, 1e-6, 0.4))
  b <- make(c(2L, 2L, 3L), c(3L, 4L, 4L), c(-1, 4, 0), c(0.7, 1e-300, 1))
  extra <- mgcvST:::.mgcvst_pairs_frame(c(4L, 5L), c(5L, 6L), status = 3L)
  files <- mgcvST:::.mgcvst_shard_file(pair_dir, 1:2)
  mgcvST:::.mgcvst_write_parquet(a, files[1L])
  mgcvST:::.mgcvst_write_parquet(b, files[2L])
  lp <- c(a$log_p_two_sided, b$log_p_two_sided)
  by <- mgcvST:::.mgcvst_log_adjust(lp, "BY")$log_q

  out <- mgcvST:::.mgcvst_pairs_finalize(pair_dir, files, c(3L, 3L), extra, "BY",
                                         0.05, temporary = FALSE)
  expect_true(out$adjustment$computed)
  expect_equal(out$n_pairs, 8)
  expect_length(out$shards, 3L)
  expect_true(all(file.exists(out$shards)) && all(file.exists(files)))
  r <- out$results
  expect_named(r, c("i", "j", "score", "log_p_two_sided", "log_p_positive",
                    "log_p_negative", "log_q", "remainder_kind", "status"))
  expect_false(is.unsorted(r$i * 10L + r$j))
  key <- function(x) paste(x$i, x$j)
  expect_equal(r$log_q[match(c(key(a), key(b)), key(r))], by, tolerance = 1e-14)
  expect_true(all(is.na(r$log_q[r$status == 3L])))
  # Discoveries are the adjusted two-sided hits, split by the score sign.
  hit <- !is.na(by) & by <= log(0.05)
  score <- c(a$score, b$score)
  expect_equal(out$discoveries$pairs_discovered, sum(hit))
  expect_equal(out$discoveries$pairs_discovered_positive, sum(hit & score > 0))
  expect_equal(out$discoveries$pairs_discovered_negative, sum(hit & score < 0))
  expect_equal(sum(hit), 3)
  expect_equal(out$discoveries$pairs_discovered,
               out$discoveries$pairs_discovered_positive +
                 out$discoveries$pairs_discovered_negative)
  expect_equal(out$discoveries$pairs_requested, 8)
  expect_equal(out$discoveries$pairs_tested, 6)
  expect_equal(out$discoveries$pairs_with_p_value, 6)
  expect_equal(out$threshold$log_p_threshold, max(lp[hit]))

  # A repeated finalize rewrites the shards; a temporary run removes raw shards.
  again <- mgcvST:::.mgcvst_pairs_finalize(pair_dir, files, c(3L, 3L), extra, "BH",
                                           0.05, temporary = TRUE)
  expect_false(any(file.exists(files)))
  expect_true(all(file.exists(again$shards)))
  expect_equal(again$results$log_q[match(key(a), key(again$results))],
               mgcvST:::.mgcvst_log_adjust(lp, "BH")$log_q[1:3], tolerance = 1e-14)
})

test_that("the memory guard skips the adjustment and the in-memory table", {
  dir <- tempfile("mgcvst-guard-")
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  contract <- mgcvST:::.mgcvst_contract("exact")
  pair_dir <- mgcvST:::.mgcvst_pairs_open(dir, list(test = 2L), contract)
  frame <- mgcvST:::.mgcvst_pairs_frame(1:3, 2:4, c(1, -1, 2), log(c(0.1, 0.2, 0.3)))
  file <- mgcvST:::.mgcvst_shard_file(pair_dir, 1L)
  mgcvST:::.mgcvst_write_parquet(frame, file)
  testthat::local_mocked_bindings(
    .mgcvst_memory_probe = function(...) list(available = 100), .package = "mgcvST")
  expect_warning(out <- mgcvST:::.mgcvst_pairs_finalize(
    pair_dir, file, 3L, NULL, "BY", 0.05, temporary = FALSE), "adjustment skipped")
  expect_false(out$adjustment$computed)
  expect_match(out$adjustment$reason, "safe memory line")
  expect_null(out$results)
  expect_length(out$shards, 1L)
  z <- mgcvST:::.mgcvst_read_shard(out$shards)
  expect_true(all(is.na(z$log_q)))
  expect_identical(z$log_p_two_sided, frame$log_p_two_sided)
  expect_true(is.na(out$discoveries$pairs_discovered))
  # Unadjusted p-values need no adjustment memory.
  none <- mgcvST:::.mgcvst_pairs_finalize(pair_dir, file, 3L, NULL, "none", 0.15,
                                          temporary = FALSE)
  expect_true(none$adjustment$computed)
  expect_identical(mgcvST:::.mgcvst_read_shard(none$shards)$log_q, frame$log_p_two_sided)
  expect_equal(none$discoveries$pairs_discovered, 1)
})

test_that("adjustment arguments are validated and shared by both tests", {
  f <- st_fixture()
  fit <- mgcvST.estimate(f$Y, f$model, BPPARAM = BiocParallel::SerialParam())
  expect_error(mgcvST.test(fit, adjust = "holm"), "should be one of")
  expect_error(mgcvST.test(fit, q.value = 0), "q.value")
  expect_error(mgcvST.test(fit, threads = 0L), "threads")
  expect_error(mgcvST.test(fit, chunk_size = 0), "chunk_size")
  expect_error(mgcvST.test(fit, verbose = NA), "verbose")
  expect_error(mgcvST.test(fit, checkpoint_dir = c("a", "b")), "checkpoint_dir")
  by <- mgcvST.test(fit)
  for (adjust in c("BH", "Sidak", "none")) {
    x <- mgcvST.test(fit, adjust = adjust)
    expect_equal(x$results$log_p_two_sided, by$results$log_p_two_sided)
    expect_identical(x$threshold$adjust, adjust)
    expect_equal(x$results$log_q, mgcvST:::.mgcvst_log_adjust(x$results$log_p_two_sided,
                                                              adjust)$log_q,
                 tolerance = 1e-14)
  }
  expect_equal(by$results$log_q, mgcvST:::.mgcvst_log_adjust(by$results$log_p_two_sided,
                                                             "BY")$log_q)
  expect_s3_class(by, "mgcvST_test")
  expect_output(print(by), "adjustment: BY")
})
