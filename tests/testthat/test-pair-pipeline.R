test_that("bounded pair blocks stream every pair once whatever the state budget", {
  .pair_pipeline_mocks()
  fit <- .pair_pipeline_fit()
  pairs <- t(utils::combn(4L, 2L))
  original <- mgcvST:::.mgcvst_spa_pairs
  seen <- list()
  testthat::local_mocked_bindings(
    .mgcvst_spa_pairs = function(index, active, states, G, threads, order = 4L) {
      seen[[length(seen) + 1L]] <<- c(rows = nrow(index), features = length(active))
      original(index, active, states, G, threads, order)
    },
    .package = "mgcvST"
  )
  # The resident-state budget is 0.7 of the available memory less a reserve;
  # the mocked probe sets it to `cache_bytes`. The reserve repeats the formula
  # of the pipeline for these two-column states, six pairs and the resident
  # pair bases of four features (8 * 2 * 2 * 4 bytes).
  run <- function(index, cache_bytes, chunk_size = 9L) {
    reserve <- 4 * 8 * 2^2 + 2 * (8 * (2^2 + 2) + 2048) +
      256 * min(6, chunk_size) + 64 * 1024^2 + 8 * 2 * 2 * 4
    testthat::local_mocked_bindings(
      .mgcvst_memory_probe = function(...) list(available = (cache_bytes + reserve) / 0.7),
      .package = "mgcvST")
    z <- mgcvST:::.mgcvst_pair_pipeline(fit, index, threads = 1L,
      chunk_size = chunk_size, verbose = FALSE)
    out <- do.call(rbind, lapply(z$shards, mgcvST:::.mgcvst_read_shard))
    expect_identical(sum(z$rows), nrow(out))
    out[order(out$i, out$j), ]
  }
  for (index in list(NULL, cbind(pairs[, 1L], pairs[, 2L]))) {
    seen <- list()
    full <- run(index, 1e5)
    seen <- list()
    bounded <- run(index, 4200)
    expect_identical(nrow(bounded), nrow(pairs))
    expect_equal(bounded, full, tolerance = 1e-12)
    expect_identical(unname(as.matrix(bounded[, c("i", "j")])), unname(pairs))
    expect_gt(length(seen), 1L)
    expect_true(all(vapply(seen, `[[`, numeric(1L), "features") <= 2L))
    expect_true(all(full$status == 0L))
    expect_true(all(is.finite(full$log_p_two_sided)))
  }
  small <- run(NULL, 1e5, chunk_size = 2L)
  expect_equal(small, run(NULL, 1e5), tolerance = 1e-12)
  expect_true(all(vapply(seen, `[[`, numeric(1L), "rows") <= 9L))
})

test_that("a pair of a failed state is returned with status 3, never a p-value", {
  .pair_pipeline_mocks()
  testthat::local_mocked_bindings(
    .mgcvst_pair_build_batch = function(fit, ids, threads, native) {
      lapply(ids, function(i) {
        if (i == 2L) list(error = "state failed") else
          list(a = c(i, i + 0.25), M = diag(c(i + 0.5, i + 1)), width = 2L)
      })
    }, .package = "mgcvST"
  )
  z <- mgcvST:::.mgcvst_pair_pipeline(.pair_pipeline_fit(3L), NULL, 1L, 10L, FALSE)
  out <- mgcvST:::.mgcvst_read_shard(z$shards[1L])
  bad <- out$i == 2L | out$j == 2L
  expect_true(all(out$status[bad] == 3L))
  expect_true(all(is.na(out$score[bad]) & is.na(out$log_p_two_sided[bad])))
  expect_identical(out$status[!bad], 0L)
  expect_identical(z$failed, data.frame(feature_id = "g2", error = "state failed"))
})

test_that("shared preparation keeps dense native state contracts", {
  fit <- list(
    working_error = matrix(0, 3L, 2L),
    working_variance = matrix(1, 3L, 2L),
    dispersion = c(1, 2), smoothing_parameters = matrix(c(2, 4), 2L),
    nuisance_covariance = list(matrix(1), matrix(1)),
    geometry = list(X = matrix(1, 3L, 1L))
  )
  seen <- new.env(parent = emptyenv())
  seen$scale <- NULL
  native_batch <- function(T0, variance, error, scale, X, nuisance, threads) {
    seen$scale <- scale
    lapply(seq_along(scale), function(i) list(a = c(i, i + 1),
                                              H = diag(c(i, i + 1))))
  }
  testthat::local_mocked_bindings(
    mgcvst_dense_score_batch_cpp = native_batch,
    .package = "mgcvST"
  )
  native <- list(T0 = matrix(1, 3L, 2L), X = matrix(1, 3L, 1L),
                 sp_index = 1L, width = c(global = 2L))
  model <- mgcvST:::.mgcvst_pair_build_batch(fit, 1:2, 2L, native)
  expect_equal(seen$scale, c(0.5, 0.5))
  expect_identical(model[[1L]]$a, c(1, 2))
  expect_equal(model[[2L]]$M, diag(c(2, 3)))
  expect_identical(model[[1L]]$width, c(global = 2L))
})

test_that("the fused exact pair kernel matches the trace powers and the full-spectrum saddlepoint", {
  set.seed(20260924)
  q <- 6L
  K <- 8L
  H <- lapply(seq_len(K), function(k) {
    z <- matrix(rnorm(q * q), q, q)
    crossprod(z) + diag(q) * 0.1
  })
  a <- matrix(rnorm(q * K), q, K)

  idx <- which(upper.tri(matrix(0, K, K)), arr.ind = TRUE)
  left <- idx[, 2L]
  right <- idx[, 1L]
  ord <- order(left)
  left <- left[ord]
  right <- right[ord]

  G <- mgcvST:::mgcvst_pair_basis_cpp(H, diag(q), 1L)
  new <- mgcvST:::mgcvst_pair_spa_cpp(H, G, a, left, right, threads = 1L)
  old_score <- colSums(a[, left, drop = FALSE] * a[, right, drop = FALSE])
  expect_equal(new$score, old_score, tolerance = 1e-12)
  expect_true(all(new$status == 0L))

  # The saddlepoint of each pair from the trace powers of the product, which
  # are the power sums of the singular spectrum.
  pairs <- cbind(left, right)
  moments <- mgcvST:::mgcvst_pair_trace_powers_cpp(H, pairs, maxPower = 4L, threads = 1L)
  spectra <- lapply(seq_along(left), function(p) {
    svd(.spa_sqrtm(H[[left[p]]]) %*% .spa_sqrtm(H[[right[p]]]))$d
  })
  for (p in seq_along(left)) {
    expect_equal(.spa_powers(spectra[[p]]), moments[p, ], tolerance = 1e-9)
  }
  reference <- vapply(seq_along(left), function(p) {
    .spa_ref(abs(old_score[p]), spectra[[p]])
  }, numeric(1L))
  expect_equal(new$log_p_two_sided, reference, tolerance = 1e-7)

  # A constructed strong pair whose p underflows to 0 has a finite
  # log_p_two_sided well below log(1e-300).
  strong_H <- lapply(1:2, function(k) diag(rep(100, q)))
  strong_a <- cbind(rep(600, q), rep(600, q))
  strong_G <- mgcvST:::mgcvst_pair_basis_cpp(strong_H, diag(q), 1L)
  strong <- mgcvST:::mgcvst_pair_spa_cpp(strong_H, strong_G, strong_a, 1L, 2L,
                                         threads = 1L)
  expect_identical(strong$status, 0L)
  expect_true(is.finite(strong$log_p_two_sided))
  expect_lt(strong$log_p_two_sided, log(1e-300))
  expect_equal(exp(strong$log_p_two_sided), 0)
})

test_that("an invalid exact-route p-value is missing and never enters the adjustment", {
  testthat::local_mocked_bindings(
    mgcvst_pair_spa_cpp = function(H, G, avec, left, right, threads, order) {
      list(score = c(1, 2), log_p_two_sided = c(-Inf, -3),
           log_p_positive = c(-Inf, -3.7), log_p_negative = c(-Inf, -0.02),
           remainder_kind = c(0L, 2L), status = c(2L, 0L))
    }, .package = "mgcvST")
  states <- lapply(1:3, function(i) list(a = c(i, 1), M = diag(2)))
  out <- mgcvST:::.mgcvst_spa_pairs(rbind(c(1L, 2L), c(1L, 3L)), 1:3, states,
                                    vector("list", 3L), 1L)
  expect_identical(out$status, c(2L, 0L))
  expect_equal(out$score, c(1, 2))
  expect_identical(out$remainder_kind, c(0L, 2L))
  expect_true(all(is.na(out[1L, c("log_p_two_sided", "log_p_positive",
                                  "log_p_negative")])))
  expect_equal(out$log_p_two_sided[2L], -3)
  adjusted <- mgcvST:::.mgcvst_log_adjust(out$log_p_two_sided, "BY")
  expect_equal(adjusted$n, 1)
  expect_true(is.na(adjusted$log_q[1L]))
})

test_that("stale pair results are refused before work and the pair directory follows the states", {
  .pair_pipeline_mocks()
  built <- 0L
  testthat::local_mocked_bindings(
    .mgcvst_pair_build_batch = function(fit, ids, threads, native) {
      built <<- built + length(ids)
      lapply(ids, function(i) list(a = c(i, i + 0.25), M = diag(c(i + 0.5, i + 1)),
                                   width = 2L))
    }, .package = "mgcvST")
  fit <- .pair_pipeline_fit()
  dir <- tempfile("mgcvst-pair-open-")
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  mgcvST:::.mgcvst_pair_pipeline(fit, NULL, 1L, 10L, FALSE, checkpoint_dir = dir)
  expect_identical(built, 4L)
  # Leave one missing state and a pair directory of another contract.
  unlink(list.files(dir, "^pairs-", full.names = TRUE), recursive = TRUE)
  unlink(file.path(dir, "feature-0000000004.rds"))
  stale <- file.path(dir, "pairs-0123")
  dir.create(stale)
  saveRDS(list(), file.path(stale, "block-0000000001.rds"))
  built <- 0L
  expect_error(mgcvST:::.mgcvst_pair_pipeline(fit, NULL, 1L, 10L, FALSE,
                                              checkpoint_dir = dir),
               "different algorithm contract")
  expect_identical(built, 0L)
  expect_false(file.exists(file.path(dir, "feature-0000000004.rds")))
  unlink(stale, recursive = TRUE)

  # A run whose state preparation fails leaves no pair directory behind.
  testthat::local_mocked_bindings(
    .mgcvst_pair_build_batch = function(fit, ids, threads, native) stop("state failure"),
    .package = "mgcvST")
  expect_error(mgcvST:::.mgcvst_pair_pipeline(fit, NULL, 1L, 10L, FALSE,
                                              checkpoint_dir = dir),
               "state failure")
  expect_length(list.files(dir, "^pairs-"), 0L)
})
