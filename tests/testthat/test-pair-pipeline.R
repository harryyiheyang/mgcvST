test_that("bounded pair blocks stream every pair once whatever the state budget", {
  .pair_pipeline_mocks()
  fit <- .pair_pipeline_fit()
  pairs <- t(utils::combn(4L, 2L))
  original <- spdeST:::.mgcvst_spa_pairs
  seen <- list()
  testthat::local_mocked_bindings(
    .mgcvst_spa_pairs = function(index, active, states, G, threads, order = 4L) {
      seen[[length(seen) + 1L]] <<- c(rows = nrow(index), features = length(active))
      original(index, active, states, G, threads, order)
    },
    .package = "spdeST"
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
      .package = "spdeST")
    z <- spdeST:::.mgcvst_pair_pipeline(fit, index, threads = 1L,
      chunk_size = chunk_size, verbose = FALSE, basis = .pair_pipeline_basis)
    out <- do.call(rbind, lapply(z$shards, spdeST:::.mgcvst_read_shard))
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
  .pair_pipeline_mocks(materialize = function(ids) lapply(ids, function(i) {
    if (i == 2L) list(error = "state failed") else .pair_pipeline_state(i)
  }))
  z <- spdeST:::.mgcvst_pair_pipeline(.pair_pipeline_fit(3L), NULL, 1L, 10L, FALSE,
                                      basis = .pair_pipeline_basis)
  out <- spdeST:::.mgcvst_read_shard(z$shards[1L])
  bad <- out$i == 2L | out$j == 2L
  expect_true(all(out$status[bad] == 3L))
  expect_true(all(is.na(out$score[bad]) & is.na(out$log_p_two_sided[bad])))
  expect_identical(out$status[!bad], 0L)
  expect_identical(z$failed, data.frame(feature_id = "g2", error = "state failed"))
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

  G <- spdeST:::mgcvst_pair_basis_cpp(H, diag(q), 1L)
  new <- spdeST:::mgcvst_pair_spa_cpp(H, G, a, left, right, threads = 1L)
  old_score <- colSums(a[, left, drop = FALSE] * a[, right, drop = FALSE])
  expect_equal(new$score, old_score, tolerance = 1e-12)
  expect_true(all(new$status == 0L))

  # The saddlepoint of each pair from the trace powers of the product, which
  # are the power sums of the singular spectrum.
  pairs <- cbind(left, right)
  moments <- spdeST:::mgcvst_pair_trace_powers_cpp(H, pairs, maxPower = 4L, threads = 1L)
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
  strong_G <- spdeST:::mgcvst_pair_basis_cpp(strong_H, diag(q), 1L)
  strong <- spdeST:::mgcvst_pair_spa_cpp(strong_H, strong_G, strong_a, 1L, 2L,
                                         threads = 1L)
  expect_identical(strong$status, 0L)
  expect_true(is.finite(strong$log_p_two_sided))
  expect_lt(strong$log_p_two_sided, log(1e-300))
  expect_equal(exp(strong$log_p_two_sided), 0)
})

test_that("an invalid exact-route p-value is p = 1, keeps its status and stays in the adjustment", {
  testthat::local_mocked_bindings(
    mgcvst_pair_spa_cpp = function(H, G, avec, left, right, threads, order) {
      list(score = c(1, 2), log_p_two_sided = c(-Inf, -3),
           log_p_positive = c(-Inf, -3.7), log_p_negative = c(-Inf, -0.02),
           remainder_kind = c(0L, 2L), status = c(2L, 0L), t1 = c(1, 1))
    }, .package = "spdeST")
  states <- lapply(1:3, function(i) list(a = c(i, 1), M = diag(2)))
  out <- spdeST:::.mgcvst_spa_pairs(rbind(c(1L, 2L), c(1L, 3L)), 1:3, states,
                                    vector("list", 3L), 1L)
  expect_identical(out$status, c(2L, 0L))
  expect_equal(out$score, c(1, 2))
  expect_identical(out$remainder_kind, c(0L, 2L))
  expect_identical(unlist(out[1L, c("log_p_two_sided", "log_p_positive",
                                    "log_p_negative")], use.names = FALSE), c(0, 0, 0))
  expect_equal(out$log_p_two_sided[2L], -3)
  adjusted <- spdeST:::.mgcvst_log_adjust(out$log_p_two_sided, "BY")
  expect_equal(adjusted$n, 2)
  expect_equal(exp(adjusted$log_q), stats::p.adjust(exp(out$log_p_two_sided), "BY"),
               tolerance = 1e-12)
})

test_that("a run whose state preparation fails leaves no file behind", {
  fit <- .pair_pipeline_fit()
  before <- list.files(tempdir(), all.files = TRUE, no.. = TRUE)
  .pair_pipeline_mocks(materialize = function(ids) stop("state failure"))
  expect_error(spdeST:::.mgcvst_pair_pipeline(fit, NULL, 1L, 10L, FALSE,
                                              basis = .pair_pipeline_basis),
               "state failure")
  expect_identical(setdiff(list.files(tempdir(), all.files = TRUE, no.. = TRUE), before),
                   character())
})
