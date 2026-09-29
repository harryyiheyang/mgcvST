.fp16_scale_fixture <- function(scales = c(1, 1), signs = rep(1, length(scales))) {
  A <- as(Matrix::Diagonal(4L), "generalMatrix")
  X <- Matrix::Matrix(1, 4L, 1L, sparse = TRUE)
  Z <- diag(4L)[, 2:4, drop = FALSE]
  a <- cbind(c(0, 1, 1, 0.2), c(0, 0.8, 0.4, -0.2))
  cache <- mgcvST:::mgcvst_fp16_cache_cpp(length(scales), 3L)
  for (j in seq_along(scales)) {
    z <- Z * sqrt(scales[j])
    mgcvST:::mgcvst_fp16_build_cpp(
      cache, as.integer(j), A, A, c(1, 0, 0, 0), X,
      matrix(0, 4L, 1L), matrix(0, 1L, 1L), matrix(0, 4L, 1L),
      0L, NA_real_, 1, 1, signs[j] * a[, 1L + (j - 1L) %% 2L, drop = FALSE],
      z, z, NA_real_, 1L
    )
  }
  list(cache = cache, a = a[2:4, , drop = FALSE],
       M = diag(0.5, 3L) - matrix(0.1, 3L, 3L))
}

test_that("fp16 Liu is invariant to per-gene curvature units", {
  scales <- c(1, 1, 1e-240, 1e-240, 1e240, 1e240,
              1e-240, 1e240, 1e-6, 1e6)
  f <- .fp16_scale_fixture(scales)
  on.exit(mgcvST:::mgcvst_fp16_cache_release_cpp(f$cache))
  info <- mgcvST:::mgcvst_fp16_cache_info_cpp(f$cache)
  expect_identical(info$state, rep(1L, length(scales)))
  expect_equal(as.numeric(info$scale) / scales, rep(0.4, length(scales)))
  expect_equal(info$gene_bytes, 2 * 6 + 8 * 4)
  left <- seq.int(1L, 9L, 2L)
  out <- mgcvST:::mgcvst_fp16_pairs_cpp(
    f$cache, seq_along(scales), 1L, length(scales), left, left + 1L, 1L
  )
  U <- sum(f$a[, 1L] * f$a[, 2L])
  expect_true(all(is.finite(out$mlog10p)))
  expect_equal(out$mlog10p, rep(out$mlog10p[1L], length(left)), tolerance = 1e-10)
  units <- sqrt(scales[left]) * sqrt(scales[left + 1L])
  expect_equal(out$score / units, rep(U, length(left)), tolerance = 1e-12)
  ref <- rkhs_score_calibrate(U, f$M, f$M, method = "liu")
  expect_equal(10^-out$mlog10p, rep(ref$p_two_sided, length(left)), tolerance = 1e-5)
  threaded <- mgcvST:::mgcvst_fp16_pairs_cpp(
    f$cache, seq_along(scales), 1L, length(scales), left, left + 1L, 2L
  )
  expect_equal(threaded, out, tolerance = 1e-12)
})

test_that("scaled fp16 states preserve score sign and report zero curvature", {
  positive <- .fp16_scale_fixture()
  negative <- .fp16_scale_fixture(signs = c(1, -1))
  zero <- .fp16_scale_fixture(c(0, 1))
  on.exit({
    mgcvST:::mgcvst_fp16_cache_release_cpp(positive$cache)
    mgcvST:::mgcvst_fp16_cache_release_cpp(negative$cache)
    mgcvST:::mgcvst_fp16_cache_release_cpp(zero$cache)
  })
  p <- mgcvST:::mgcvst_fp16_pairs_cpp(positive$cache, 1:2, 1L, 1L)
  n <- mgcvST:::mgcvst_fp16_pairs_cpp(negative$cache, 1:2, 1L, 1L)
  expect_equal(n$score, -p$score)
  expect_equal(n$mlog10p, p$mlog10p)
  info <- mgcvST:::mgcvst_fp16_cache_info_cpp(zero$cache)
  expect_identical(info$state, c(2L, 1L))
  expect_match(info$error[1L], "identically zero")
  z <- mgcvST:::mgcvst_fp16_pairs_cpp(zero$cache, 1:2, 1L, 1L)
  expect_true(is.na(z$score))
  expect_true(is.na(z$mlog10p))
})

test_that("scaled fp16 states round trip and reject earlier binary formats", {
  f <- .fp16_scale_fixture(c(1e-240, 1e240))
  on.exit(mgcvST:::mgcvst_fp16_cache_release_cpp(f$cache))
  path <- tempfile(fileext = ".bin")
  mgcvST:::mgcvst_fp16_write_cpp(f$cache, 1:2, 1:2, c("a", "b"), "scale-test", path)
  target <- mgcvST:::mgcvst_fp16_cache_cpp(2L, 3L)
  on.exit(mgcvST:::mgcvst_fp16_cache_release_cpp(target), add = TRUE)
  mgcvST:::mgcvst_fp16_read_cpp(target, 1:2, 1:2, c("a", "b"), "scale-test", path, NA_real_)
  expect_equal(mgcvST:::mgcvst_fp16_cache_info_cpp(target),
               mgcvST:::mgcvst_fp16_cache_info_cpp(f$cache))
  expect_identical(mgcvST:::mgcvst_fp16_pairs_cpp(target, 1:2, 1L, 1L),
                   mgcvST:::mgcvst_fp16_pairs_cpp(f$cache, 1:2, 1L, 1L))
  con <- file(path, "r+b")
  seek(con, 8L, rw = "write")
  writeBin(1L, con, size = 4L)
  close(con)
  old <- mgcvST:::mgcvst_fp16_cache_cpp(2L, 3L)
  on.exit(mgcvST:::mgcvst_fp16_cache_release_cpp(old), add = TRUE)
  expect_error(mgcvST:::mgcvst_fp16_read_cpp(
    old, 1:2, 1:2, c("a", "b"), "scale-test", path, NA_real_
  ), "incompatible header")
  expect_identical(mgcvST:::mgcvst_fp16_cache_info_cpp(old)$state, c(0L, 0L))
})

test_that("scaled fp16 pair checkpoints bind the gene-state signature", {
  skip_if_not_installed("arrow")
  f <- .fp16_scale_fixture()
  on.exit(mgcvST:::mgcvst_fp16_cache_release_cpp(f$cache))
  state <- list(cache = f$cache, used = 1:2, n = 2L,
                root = tempfile("fp16-explicit-"), signature = "current-state")
  dir.create(state$root)
  pairs <- cbind(i = 1L, j = 2L)
  out <- mgcvST:::.mgcvst_inla_fp16_explicit_shard(state, pairs, 1L, TRUE, FALSE)
  again <- mgcvST:::.mgcvst_inla_fp16_explicit_shard(state, pairs, 1L, TRUE, FALSE)
  expect_identical(lapply(again$result, unname), lapply(out$result, unname))
  manifest <- readRDS(file.path(state$root, "pairs", "manifest.rds"))
  expect_identical(manifest$version, 2L)
  expect_identical(manifest$state_signature, state$signature)
  state$signature <- "earlier-state"
  expect_error(mgcvST:::.mgcvst_inla_fp16_explicit_shard(state, pairs, 1L, TRUE, FALSE),
               "different state format, gene states")
  state$root <- tempfile("fp16-stream-")
  dir.create(state$root)
  out <- mgcvST:::.mgcvst_inla_fp16_stream(state, 1L, 1L, FALSE)
  again <- mgcvST:::.mgcvst_inla_fp16_stream(state, 1L, 1L, FALSE)
  expect_identical(again$resumed_pairs, 1)
  state$signature <- "current-state"
  expect_error(mgcvST:::.mgcvst_inla_fp16_stream(state, 1L, 1L, FALSE),
               "different state format, gene states")
})
