test_that("chunk workers hold BLAS to one thread", {
  seen <- 0L
  testthat::local_mocked_bindings(
    .mgcvst_thread_limit = function() { seen <<- seen + 1L; invisible(NULL) },
    .inlast_null_chunk = function(payload) payload * 2,
    .package = "spdeST"
  )
  withr::local_libpaths(.libPaths())
  # The library that holds the loaded namespace, so the task keeps it loaded.
  lib <- dirname(getNamespaceInfo(asNamespace("spdeST"), "path"))
  task <- spdeST:::.inlast_chunk_task()
  # A worker (another process id) limits its threads before the work.
  expect_identical(task(3, ".inlast_null_chunk", list(), lib, manager = -1L), 6)
  expect_identical(seen, 1L)
  # The manager itself (a serial backend) is limited by the entry point.
  expect_identical(task(3, ".inlast_null_chunk", list(), lib,
                        manager = Sys.getpid()), 6)
  expect_identical(seen, 1L)
})
