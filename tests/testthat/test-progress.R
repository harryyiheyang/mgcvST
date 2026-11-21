test_that("progress lines appear at each tenth and at the end", {
  withr::local_options(spdeST.progress = TRUE)
  msgs <- testthat::capture_messages({
    p <- spdeST:::.mgcvst_progress("null", 100, "features")
    for (i in 1:20) p$add(5)
    p$done()
  })
  expect_match(msgs[1L], "^\\[null\\] 100 features")
  lines <- msgs[-1L]
  expect_length(lines, 10L)
  expect_match(lines[1L], "10 / 100 features \\(10%\\)")
  expect_match(lines[10L], "100 / 100 features \\(100%\\) .*ETA 0s")
})

test_that("restored work counts as done and the option silences progress", {
  withr::local_options(spdeST.progress = TRUE)
  msgs <- testthat::capture_messages({
    p <- spdeST:::.mgcvst_progress("pairs", 1e6, "pairs", resumed = 5e5)
    p$add(1e5)
    p$add(4e5, fresh = FALSE)
    p$done()
  })
  expect_match(msgs[1L], "1,000,000 pairs, 500,000 restored")
  expect_match(msgs[2L], "600,000 / 1,000,000 pairs \\(60%\\)")
  expect_length(msgs, 3L)
  withr::local_options(spdeST.progress = FALSE)
  expect_silent({
    p <- spdeST:::.mgcvst_progress("null", 10)
    p$add(10)
    p$done()
  })
})

test_that("a heartbeat line appears between tenths", {
  withr::local_options(spdeST.progress = TRUE)
  msgs <- testthat::capture_messages({
    p <- spdeST:::.mgcvst_progress("spatial", 1000, heartbeat = 0)
    p$add(1)
    p$add(1)
  })
  expect_length(msgs, 3L)
  expect_match(msgs[3L], "2 / 1,000 features \\(0%\\)")
})

test_that("the chunk runner reports features as chunks complete", {
  withr::local_options(spdeST.progress = TRUE)
  groups <- split(1:30, rep(1:3, each = 10))
  worker <- function(payload) as.list(payload$index)
  msgs <- testthat::capture_messages(
    out <- spdeST:::.mgcvst_run_chunks(groups, function(index) list(index = index),
      "null", NULL, "digest", rep(FALSE, 30), BiocParallel::SerialParam(),
      worker)
  )
  expect_identical(unlist(out$results), 1:30)
  expect_match(msgs[1L], "\\[null\\] 30 features")
  expect_match(msgs[length(msgs)], "30 / 30 features \\(100%\\)")
})
