test_that("score-state shards round-trip at both storage precisions", {
  for (storage in c("double", "float32")) {
    path <- tempfile("mgcvst-score-store-")
    on.exit(unlink(path, recursive = TRUE), add = TRUE)
    store <- mgcvST:::.mgcvst_store_open(
      path, signature = list(fit = "example", basis_rank = 2L),
      feature_ids = c("gene/one", "gene:two"), storage = storage
    )
    M <- matrix(c(0.75, 0.125, 0.125, 1.25), 2L)
    state <- list(a = c(1.125, -0.375), M = M, width = c(global = 2L))
    expect_false(mgcvST:::.mgcvst_store_has(store, "gene/one"))
    mgcvST:::.mgcvst_store_write(store, "gene/one", state)
    expect_true(mgcvST:::.mgcvst_store_has(store, 1L))
    actual <- mgcvST:::.mgcvst_store_read(store, "gene/one")
    expect_identical(actual$a, state$a)
    expect_equal(actual$M, M, tolerance = if (storage == "double") 0 else 1e-7)
    expect_identical(actual$width, state$width)
    expect_error(mgcvST:::.mgcvst_store_write(store, 1L, state),
                 "already exists")
  }
})

test_that("score-state stores resume only with identical metadata", {
  path <- tempfile("mgcvst-score-resume-")
  on.exit(unlink(path, recursive = TRUE), add = TRUE)
  signature <- list(fit = "a", projection = list(coverage = 0.995, rank = 2L))
  store <- mgcvST:::.mgcvst_store_open(
    path, signature = signature, feature_ids = c("a", "b")
  )
  state <- list(a = c(1, 2), M = diag(c(0.5, 1.5)))
  mgcvST:::.mgcvst_store_write(store, "a", state)
  reopened <- mgcvST:::.mgcvst_store_open(
    path, signature = signature, feature_ids = c("a", "b"), resume = TRUE
  )
  expect_true(mgcvST:::.mgcvst_store_has(reopened, "a"))
  expect_false(mgcvST:::.mgcvst_store_has(reopened, "b"))
  expect_equal(mgcvST:::.mgcvst_store_read(reopened, "a")$M, state$M)
  expect_error(mgcvST:::.mgcvst_store_open(
    path, signature = list(fit = "other"), feature_ids = c("a", "b")
  ), "do not match")
  expect_error(mgcvST:::.mgcvst_store_open(
    path, signature = signature, feature_ids = c("b", "a")
  ), "do not match")
  expect_error(mgcvST:::.mgcvst_store_open(
    path, signature = signature, feature_ids = c("a", "b"), storage = "float32"
  ), "do not match")
  expect_error(mgcvST:::.mgcvst_store_open(
    path, signature = signature, feature_ids = c("a", "b"), resume = FALSE
  ), "already exists")
})

test_that("incomplete writes are ignored and damaged shards fail clearly", {
  path <- tempfile("mgcvst-score-partial-")
  on.exit(unlink(path, recursive = TRUE), add = TRUE)
  store <- mgcvST:::.mgcvst_store_open(
    path, signature = list(fit = "partial"), feature_ids = c("a", "b")
  )
  writeBin(as.raw(c(1L, 2L, 3L)),
           file.path(path, "feature-0000000001.rds-interrupted.tmp"))
  reopened <- mgcvST:::.mgcvst_store_open(
    path, signature = list(fit = "partial"), feature_ids = c("a", "b")
  )
  expect_false(mgcvST:::.mgcvst_store_has(reopened, "a"))
  mgcvST:::.mgcvst_store_write(reopened, "a", list(a = 1, M = matrix(2)))
  expect_equal(mgcvST:::.mgcvst_store_read(reopened, "a")$M, matrix(2))
  writeBin(as.raw(c(1L, 2L, 3L)),
           file.path(path, "feature-0000000001.rds"))
  expect_error(mgcvST:::.mgcvst_store_read(reopened, "a"), "unreadable")
  expect_error(mgcvST:::.mgcvst_store_write(reopened, "a",
    list(a = 1, M = matrix(2))), "already exists")
})

test_that("failed features are isolated and only temporary stores are cleaned", {
  path <- tempfile("mgcvst-score-errors-")
  on.exit(unlink(path, recursive = TRUE), add = TRUE)
  store <- mgcvST:::.mgcvst_store_open(
    path, signature = list(fit = "errors"), feature_ids = c("bad", "good")
  )
  mgcvST:::.mgcvst_store_write(store, "bad", list(error = "fit failed"))
  mgcvST:::.mgcvst_store_write(store, "good",
    list(a = c(1, 2), M = diag(2)))
  expect_identical(mgcvST:::.mgcvst_store_read(store, "bad"),
                   list(error = "fit failed"))
  expect_equal(mgcvST:::.mgcvst_store_read(store, "good")$M, diag(2))
  expect_false(mgcvST:::.mgcvst_store_cleanup(store))
  expect_true(dir.exists(path))

  temporary <- mgcvST:::.mgcvst_store_open(
    signature = list(fit = "temp"), feature_ids = "a"
  )
  temp_path <- temporary$path
  expect_true(dir.exists(temp_path))
  expect_true(mgcvST:::.mgcvst_store_cleanup(temporary))
  expect_false(dir.exists(temp_path))
})
