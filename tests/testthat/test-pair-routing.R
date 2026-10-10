test_that("auto takes the exact route when its pair phase, bases and store fit the budgets", {
  resolve <- function(q, n_used, n_pairs, threads = 4L, rank = 20L, memory = 64 * 1024^3,
                      moments = "auto") {
    mgcvST:::.mgcvst_route_resolve(moments, q, n_used, n_pairs, threads, rank,
                                   k_exact = 20L, k_pca = 50L, available_memory = memory)
  }
  # Exact cost: 5 ms * (q / 298)^3 per pair and thread.
  expect_equal(mgcvST:::.mgcvst_exact_pair_seconds(298), 5e-3)
  expect_equal(mgcvST:::.mgcvst_exact_pair_seconds(596), 8 * 5e-3)
  small <- resolve(q = 298L, n_used = 1025L, n_pairs = 525000)
  expect_identical(small$moments, "exact")
  expect_identical(small$k, 20L)
  expect_equal(small$seconds, 525000 * 5e-3 / 4)
  expect_match(small$reason, "within 2.0 h")
  # 9.2 million pairs at q = 298 take 3.2 hours on 4 threads: PCAlearning.
  large <- resolve(q = 298L, n_used = 4300L, n_pairs = 9.2e6)
  expect_identical(large$moments, "pcalearning")
  expect_identical(large$k, 50L)
  expect_match(large$reason, "exact pair phase would take an estimated")
  expect_lt(large$seconds, large$seconds_exact)
  # More threads bring the same case under the budget.
  expect_identical(resolve(q = 298L, n_used = 4300L, n_pairs = 9.2e6, threads = 8L)$moments,
                   "exact")
  # A large score dimension is PCAlearning whatever the number of pairs.
  expect_identical(resolve(q = 1404L, n_used = 200L, n_pairs = 19900)$moments, "exact")
  expect_identical(resolve(q = 1404L, n_used = 2000L, n_pairs = 2e6)$moments, "pcalearning")
  # Pair bases above 30% of the memory, or a store above 64 GB, are PCAlearning.
  memory <- resolve(q = 298L, n_used = 4300L, n_pairs = 1e5, memory = 4300 * 298 * 20 * 8 / 0.5)
  expect_identical(memory$moments, "pcalearning")
  expect_match(memory$reason, "pair bases need")
  disk <- resolve(q = 1000L, n_used = 9000L, n_pairs = 1e4)
  expect_identical(disk$moments, "pcalearning")
  expect_match(disk$reason, "on disk")
  # With no more genes than the PCAlearning rank there is nothing to learn.
  few <- resolve(q = 1404L, n_used = 15L, n_pairs = 1e9)
  expect_identical(few$moments, "exact")
  expect_match(few$reason, "do not exceed the PCAlearning rank")
  # k never exceeds q, and an explicit choice is kept.
  expect_identical(resolve(q = 12L, n_used = 100L, n_pairs = 10)$k, 12L)
  expect_identical(resolve(q = 298L, n_used = 4300L, n_pairs = 9.2e6, moments = "exact")$moments,
                   "exact")
  expect_identical(resolve(q = 12L, n_used = 3L, n_pairs = 3, moments = "pcalearning")$moments,
                   "pcalearning")
  expect_null(resolve(q = 12L, n_used = 3L, n_pairs = 3, moments = "exact")$reason)
  expect_identical(mgcvST:::.mgcvst_format_duration(30), "30 s")
  expect_identical(mgcvST:::.mgcvst_format_duration(7200), "2.0 h")
  expect_identical(mgcvST:::.mgcvst_format_duration(3 * 86400), "3.0 days")
})

test_that("mgcvST.test reports its route and the estimated time, and the contract records it", {
  f <- st_fixture()
  fit <- suppressWarnings(mgcvST.estimate(f$Y, f$model, BPPARAM = BiocParallel::SerialParam(),
                                          spatial = "all"))
  q <- mgcvST:::.mgcvst_state_width(fit)
  messages <- testthat::capture_messages(auto <- mgcvST.test(fit, verbose = TRUE))
  expect_match(paste(messages, collapse = ""),
               "Pair test: exact moments \\(k = 20\\) on q = 24, 3 pairs and 1 thread; estimated pair time 0 s \\(auto: 3 genes do not exceed the PCAlearning rank 20\\)")
  expect_identical(auto$moments, "exact")
  expect_identical(auto$calibration, "saddlepoint")
  expect_identical(auto$timing$route$moments, "exact")
  expect_identical(auto$timing$route$q, q)
  contract <- auto$contract
  expect_identical(contract$calibration_contract, "spa_v1")
  expect_identical(contract$route, "exact")
  expect_identical(contract$k, 20L)
  expect_identical(contract$remainder_order, 4L)
  expect_match(contract$basis_sha, "^[0-9a-f]{64}$")
  expect_identical(contract$kernel_version, 1L)
  expect_true(all(auto$results$remainder_kind %in% 0:3))

  # k = q is the full-spectrum saddlepoint; the default k = 20 differs from it
  # only by the remainder.
  full <- mgcvST.test(fit, moments = "exact", k = 1000L)
  expect_identical(full$contract$k, q)
  expect_equal(full$results$log_p_two_sided, auto$results$log_p_two_sided, tolerance = 1e-5)
  for (seen in c("exact", "pcalearning")) {
    messages <- testthat::capture_messages(
      explicit <- mgcvST.test(fit, moments = seen, rank = 2L, k = 3L, verbose = TRUE))
    expect_identical(explicit$moments, seen)
    expect_match(paste(messages, collapse = ""),
                 if (seen == "exact") "exact moments \\(k = 3\\)" else
                   "PCAlearning \\(rank 2, k = 3\\)")
    expect_false(grepl("auto:", paste(messages, collapse = ""), fixed = TRUE))
  }
  expect_identical(explicit$contract$k, 3L)
  expect_identical(explicit$contract$route, "pcalearning")
  expect_identical(explicit$contract$remainder_order, 2L)
  expect_false(is.null(explicit$pca_learning))
  expect_equal(explicit$results$log_p_two_sided, full$results$log_p_two_sided, tolerance = 0.3)
  expect_error(mgcvST.test(fit, moments = "pcalearning", rank = 4L), "achievable PCAlearning rank")
  expect_error(mgcvST.test(fit, moments = "approximate"), "should be one of")
})

test_that("a checkpoint directory records its route and a resumed auto run follows it", {
  f <- st_fixture()
  fit <- suppressWarnings(mgcvST.estimate(f$Y, f$model, BPPARAM = BiocParallel::SerialParam(),
                                          spatial = "all"))
  dir <- tempfile("mgcvst-route-")
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  first <- mgcvST.test(fit, checkpoint_dir = dir, chunk_size = 2L)
  expect_identical(first$moments, "exact")
  expect_identical(readRDS(file.path(dir, "route.rds"))$moments, "exact")
  # With a budget that would send auto to PCAlearning, the resumed run still
  # follows the recorded route.
  testthat::local_mocked_bindings(
    .mgcvst_route_limits = function() {
      list(exact_seconds = -1, store_bytes = 0, resident_fraction = 0)
    }, .package = "mgcvST")
  fresh <- mgcvST.test(fit, moments = "auto", rank = 2L, k = 3L)
  expect_identical(fresh$moments, "pcalearning")
  messages <- testthat::capture_messages(
    again <- mgcvST.test(fit, checkpoint_dir = dir, chunk_size = 2L, verbose = TRUE))
  expect_identical(again$moments, "exact")
  expect_identical(again$timing$pair_pipeline$resumed_pairs, 3)
  expect_match(paste(messages, collapse = ""), "the route of the checkpoint directory")
  expect_identical(again$results, first$results)
  # An explicit other route in the same directory is refused by that route.
  expect_error(mgcvST.test(fit, checkpoint_dir = dir, moments = "pcalearning", rank = 2L),
               "holds the score states of the exact route")
})

test_that("an INLA fit takes the exact route at small q and agrees with the full spectrum", {
  skip_on_cran()
  fit <- .pca_nb_fit()
  pairs <- t(combn(fit$feature_id, 2L))
  reference <- .pca_exact_reference(fit, pairs)
  auto <- inlaST.test(fit, pairs = pairs, threads = 2L)
  expect_identical(auto$moments, "exact")
  expect_identical(auto$contract$route, "exact")
  exact <- inlaST.test(fit, pairs = pairs, threads = 2L, moments = "exact", k = 1000L)
  expect_identical(exact$contract$k, reference$q)
  expect_equal(exact$results$score, reference$score, tolerance = 1e-10)
  expect_true(all(exact$results$status == 0L))
  expect_equal(-exact$results$log_p_two_sided / log(10), -reference$log_p / log(10),
               tolerance = 1e-6)
  # The default k = 20 differs from the full spectrum by the remainder only.
  expect_equal(-auto$results$log_p_two_sided / log(10), -reference$log_p / log(10),
               tolerance = 1e-3)
  expect_identical(auto$timing$pair_pipeline$preparation_backend, "sparse_reduced")
  # The same states serve the PCAlearning route.
  pca <- inlaST.test(fit, pairs = pairs, threads = 2L, moments = "pcalearning",
                     rank = length(fit$feature_id), k = 1000L)
  expect_equal(-pca$results$log_p_two_sided / log(10), -reference$log_p / log(10),
               tolerance = 1e-3)
  # One and four threads give identical exact results.
  one <- inlaST.test(fit, pairs = pairs, threads = 1L, moments = "exact")
  four <- inlaST.test(fit, pairs = pairs, threads = 4L, moments = "exact")
  expect_identical(one$results, four$results)
  expect_identical(one$contract$basis_sha, four$contract$basis_sha)
})
