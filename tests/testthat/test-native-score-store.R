test_that("RDS score stores preserve the format-one manifest and roundtrip", {
  path <- tempfile("mgcvst-rds-state-")
  on.exit(unlink(path, recursive = TRUE), add = TRUE)
  signature <- list(version = 9L, fit = "fixture")
  ids <- c("g1", "g2")
  store <- mgcvST:::.mgcvst_store_open(path, signature, ids)
  expect_identical(store$encoding, "rds")
  expect_identical(readRDS(file.path(path, "manifest.rds")),
    list(format = 1L, signature = signature, feature_ids = ids,
         storage = "double"))
  state <- list(a = c(1, 2), M = diag(2), width = 2L)
  mgcvST:::.mgcvst_store_write(store, "g1", state)
  mgcvST:::.mgcvst_store_write(store, "g2", list(error = "fit failed"))
  expect_true(file.exists(file.path(path, "feature-0000000001.rds")))
  expect_identical(mgcvST:::.mgcvst_store_read(store, "g1"), state)
  expect_identical(mgcvST:::.mgcvst_store_read(store, "g2"),
                   list(error = "fit failed"))
})

test_that("native score stores keep double states and validate identity", {
  skip_if(!exists("mgcvst_state_write_cpp", asNamespace("mgcvST"),
                  inherits = FALSE),
          "native score-state I/O has not been exported yet")
  path <- tempfile("mgcvst-native-state-")
  on.exit(unlink(path, recursive = TRUE), add = TRUE)
  signature <- list(version = 9L, fit = "fixture")
  ids <- c("g1", "g2")
  store <- mgcvST:::.mgcvst_store_open(
    path, signature, ids, encoding = "native"
  )
  expect_identical(store$encoding, "native")
  expect_identical(readRDS(file.path(path, "manifest.rds")),
    list(format = 2L, encoding = "native", signature = signature,
         feature_ids = ids, storage = "double"))
  state <- list(a = c(1.25, -2.5), M = matrix(c(2, .25, .25, 3), 2L),
                width = 2L)
  mgcvST:::.mgcvst_store_write(store, "g1", state)
  mgcvST:::.mgcvst_store_write(store, "g2", list(error = "fit failed"))
  expect_true(file.exists(file.path(path, "feature-0000000001.bin")))
  restored <- mgcvST:::.mgcvst_store_read(store, "g1")
  expect_equal(restored$a, state$a, tolerance = 0)
  expect_equal(restored$M, state$M, tolerance = 0)
  expect_equal(restored$width, state$width)
  expect_identical(mgcvST:::.mgcvst_store_read(store, "g2"),
                   list(error = "fit failed"))
  expect_error(mgcvST:::mgcvst_state_read_cpp(
    file.path(path, "feature-0000000001.bin"), "wrong-signature", "g1"
  ))
  expect_error(mgcvST:::.mgcvst_store_open(
    tempfile("mgcvst-invalid-native-"), signature, ids,
    storage = "float32", encoding = "native"
  ), "storage = 'double'")
})

test_that("approximate summary RDS uses format three and remains resumable", {
  path <- tempfile("mgcvst-summary-")
  on.exit(unlink(path, recursive = TRUE), add = TRUE)
  signature <- list(fit = "fixture", references = c(1L, 2L))
  ids <- c("g1", "g2")
  store <- mgcvST:::.mgcvst_approx_summary_open(path, signature, ids, TRUE)
  expect_identical(readRDS(file.path(path, "manifest.rds")),
    list(format = 3L, signature = signature, feature_ids = ids))
  z <- list(a = c(1, 2), self = c(1, 2, 3, 4),
            cross = matrix(1:8, 2L, 4L))
  mgcvST:::.mgcvst_approx_summary_write(store, 1L, z)
  expect_identical(mgcvST:::.mgcvst_approx_summary_read(store, 1L, 2L),
    c(z, list(feature_id = "g1", signature = signature)))
  expect_null(mgcvST:::.mgcvst_approx_summary_read(store, 2L, 2L))
  expect_identical(mgcvST:::.mgcvst_approx_summary_open(
    path, signature, ids, TRUE
  )$feature_ids, ids)
})

test_that("native trace summaries attach R metadata and prefer an RDS shard", {
  skip_if(!exists("mgcvst_trace_read_cpp", asNamespace("mgcvST"),
                  inherits = FALSE),
          "native trace-summary I/O has not been exported yet")
  path <- tempfile("mgcvst-native-summary-")
  on.exit(unlink(path, recursive = TRUE), add = TRUE)
  signature <- list(fit = "fixture", references = 1L)
  ids <- c("g1", "g2")
  store <- mgcvST:::.mgcvst_approx_summary_open(path, signature, ids, TRUE)
  native_path <- mgcvST:::.mgcvst_approx_summary_native_file(store, 2L)
  file.create(native_path)
  calls <- new.env(parent = emptyenv())
  calls$count <- 0L
  testthat::local_mocked_bindings(
    mgcvst_trace_read_cpp = function(path, signature, feature_id, n_ref) {
      calls$count <- calls$count + 1L
      expect_identical(feature_id, "g2")
      expect_identical(n_ref, 2L)
      list(a = c(3, 4), self = NULL, cross = matrix(2, 2L, 4L))
    },
    .package = "mgcvST"
  )
  native <- mgcvST:::.mgcvst_approx_summary_read(store, 2L, 2L)
  expect_identical(native$feature_id, "g2")
  expect_identical(native$signature, signature)
  expect_null(native$self)
  expect_identical(calls$count, 1L)
  expect_error(mgcvST:::.mgcvst_approx_summary_write(store, 2L, native),
               "already exists")

  rds_value <- list(a = c(5, 6), self = NULL,
                    cross = matrix(7, 2L, 4L))
  mgcvST:::.mgcvst_approx_summary_write(store, 1L, rds_value)
  file.create(mgcvST:::.mgcvst_approx_summary_native_file(store, 1L))
  calls$count <- 1L
  from_rds <- mgcvST:::.mgcvst_approx_summary_read(store, 1L, 2L)
  expect_equal(from_rds$a, rds_value$a)
  expect_identical(calls$count, 1L)
})
