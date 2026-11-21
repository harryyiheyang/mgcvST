test_that("score-state shards round-trip exactly", {
  store <- spdeST:::.mgcvst_store_open(c("gene/one", "gene:two"))
  on.exit(spdeST:::.mgcvst_store_cleanup(store), add = TRUE)
  M <- matrix(c(0.75, 0.125, 0.125, 1.25), 2L)
  state <- list(a = c(1.125, -0.375), M = M, width = c(global = 2L))
  expect_false(spdeST:::.mgcvst_store_has(store, "gene/one"))
  spdeST:::.mgcvst_store_write(store, "gene/one", state)
  expect_true(spdeST:::.mgcvst_store_has(store, 1L))
  actual <- spdeST:::.mgcvst_store_read(store, "gene/one")
  expect_identical(actual$a, state$a)
  expect_equal(actual$M, M, tolerance = 0)
  expect_identical(actual$width, state$width)
  expect_error(spdeST:::.mgcvst_store_write(store, 1L, state),
               "already exists")
  expect_identical(names(formals(spdeST:::.mgcvst_store_open)), "feature_ids")
})

test_that("damaged shards fail clearly", {
  store <- spdeST:::.mgcvst_store_open(c("a", "b"))
  on.exit(spdeST:::.mgcvst_store_cleanup(store), add = TRUE)
  spdeST:::.mgcvst_store_write(store, "a", list(a = 1, M = matrix(2)))
  expect_equal(spdeST:::.mgcvst_store_read(store, "a")$M, matrix(2))
  writeBin(as.raw(c(1L, 2L, 3L)),
           file.path(store$path, "feature-0000000001.rds"))
  expect_error(spdeST:::.mgcvst_store_read(store, "a"), "unreadable")
  expect_error(spdeST:::.mgcvst_store_read(store, "b"), "has not been written")
})

test_that("failed features are isolated and the store is a temporary directory that cleanup removes", {
  store <- spdeST:::.mgcvst_store_open(c("bad", "good"))
  spdeST:::.mgcvst_store_write(store, "bad", list(error = "fit failed"))
  spdeST:::.mgcvst_store_write(store, "good",
    list(a = c(1, 2), M = diag(2)))
  expect_identical(spdeST:::.mgcvst_store_read(store, "bad"),
                   list(error = "fit failed"))
  expect_equal(spdeST:::.mgcvst_store_read(store, "good")$M, diag(2))
  expect_true(startsWith(normalizePath(store$path, winslash = "/"),
                         normalizePath(tempdir(), winslash = "/")))
  expect_true(dir.exists(store$path))
  expect_true(spdeST:::.mgcvst_store_cleanup(store))
  expect_false(dir.exists(store$path))
})
