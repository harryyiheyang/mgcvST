.unit_store_fixture <- function(q = 2L) {
  list(
    a = as.numeric(seq_len(q)),
    expected_vp = matrix(1, 1L, 1L),
    K = Matrix::Diagonal(q, as.numeric(seq_len(q))),
    U = matrix(numeric(), q, 0L),
    H_L = Matrix::Diagonal(q),
    H_D = rep(1, q),
    H_perm = seq_len(q),
    hinv_g = rep(1, q),
    hden = 1,
    tau = 2,
    width = q,
    normalization = q - 1L
  )
}

test_that("sparse unit shards roundtrip, resume, and validate metadata", {
  path <- tempfile("mgcvst-unit-store-")
  on.exit(unlink(path, recursive = TRUE), add = TRUE)
  ids <- c("g1", "g2")
  signature <- list(version = 3L, fit = "fixture")
  store <- mgcvST:::.mgcvst_unit_store_open(path, signature, ids)
  unit <- .unit_store_fixture()
  mgcvST:::.mgcvst_unit_store_write(store, "g1", unit)
  mgcvST:::.mgcvst_unit_store_write(store, "g2", list(error = "fit failed"))
  expect_true(mgcvST:::.mgcvst_unit_store_has(store, "g1"))
  expect_identical(mgcvST:::.mgcvst_unit_store_read(store, "g1"), unit)
  expect_identical(mgcvST:::.mgcvst_unit_store_read(store, 2L),
                   list(error = "fit failed"))
  expect_identical(mgcvST:::.mgcvst_unit_store_open(
    path, signature, ids, resume = TRUE
  )$feature_ids, ids)
  expect_identical(store$kind, "sparse")
  expect_error(mgcvST:::.mgcvst_unit_store_open(
    path, signature, ids, kind = "dense"
  ), "do not match")
  expect_error(mgcvST:::.mgcvst_unit_store_open(
    path, list(version = 4L), ids
  ), "do not match")
  expect_error(mgcvST:::.mgcvst_unit_store_open(
    path, signature, rev(ids)
  ), "do not match")
  expect_error(mgcvST:::.mgcvst_unit_store_write(store, "g1", unit),
               "already exists")
  expect_error(mgcvST:::.mgcvst_unit_store_read(store, "missing"),
               "not in this unit store")

  dense_path <- tempfile("mgcvst-dense-unit-store-")
  on.exit(unlink(dense_path, recursive = TRUE), add = TRUE)
  dense <- mgcvST:::.mgcvst_unit_store_open(
    dense_path, signature, ids, kind = "dense"
  )
  dense_unit <- list(operator = diag(2L), target = "global", a = c(2, 3))
  mgcvST:::.mgcvst_unit_store_write(dense, "g1", dense_unit)
  expect_identical(mgcvST:::.mgcvst_unit_store_read(dense, "g1"), dense_unit)
})

test_that("sparse pair batches materialize supplied units without rebuilding", {
  fit <- list(feature_id = c("g1", "g2", "g3"), score_backend = "sparse")
  ids <- c(2L, 1L, 3L)
  units <- list(g2 = list(error = "feature two failed"),
                g1 = .unit_store_fixture(), g3 = .unit_store_fixture())
  units$g1$a <- c(11, 12)
  units$g3$a <- c(31, 32)
  seen <- new.env(parent = emptyenv())
  seen$ids <- character()
  testthat::local_mocked_bindings(
    .inlast_sparse_units = function(...) stop("units were rebuilt"),
    .inlast_sparse_materialize_reduced = function(fit, units, basis,
                                                   threads = 1L) {
      seen$ids <- names(units)
      lapply(units, function(unit) list(
        a = unit$a + 1,
        M = diag(unit$a[1L], length(unit$a)),
        width = length(unit$a)
      ))
    },
    .package = "mgcvST"
  )
  built <- mgcvST:::.mgcvst_pair_build_batch(
    fit, ids, threads = 2L, basis = list(), mode = "sparse",
    sparse_units = units
  )
  expect_identical(seen$ids, c("g1", "g3"))
  expect_identical(built[[1L]], list(error = "feature two failed"))
  expect_equal(built[[2L]]$a, c(12, 13))
  expect_equal(built[[3L]]$a, c(32, 33))
  expect_equal(built[[2L]]$M, diag(11, 2L))
  expect_equal(built[[3L]]$M, diag(31, 2L))
  expect_error(mgcvST:::.mgcvst_pair_build_batch(
    fit, ids, threads = 1L, basis = list(), mode = "sparse",
    sparse_units = units[c(2L, 1L, 3L)]
  ), "same order")
})
