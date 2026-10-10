.routing_fit <- function() {
  f <- st_fixture()
  suppressWarnings(mgcvST.estimate(f$Y, f$model, BPPARAM = BiocParallel::SerialParam(),
                                   spatial = "all"))
}

test_that("the route is chosen by the user: moments has no default and no automatic value", {
  fit <- .routing_fit()
  test <- function(...) mgcvST.test(fit, ...)
  expect_error(test(), "moments must be given")
  expect_error(test(verbose = TRUE), "moments must be given")
  expect_error(test(pairs = rbind(c(1L, 2L))), "moments must be given")
  expect_error(test(moments = "auto"), "moments must be one of")
  expect_error(test(moments = c("exact", "pcalearning")), "moments must be one of")
  expect_error(test(moments = NA_character_), "moments must be one of")
  expect_false(exists(".mgcvst_route_resolve", asNamespace("mgcvST"), inherits = FALSE))
  expect_false(exists(".mgcvst_route_limits", asNamespace("mgcvST"), inherits = FALSE))
})

test_that("mgcvST.test reports its route and k, and the contract records them", {
  fit <- .routing_fit()
  q <- mgcvST:::.mgcvst_state_width(fit)
  messages <- testthat::capture_messages(exact <- mgcvST.test(fit, moments = "exact",
                                                              verbose = TRUE))
  expect_match(paste(messages, collapse = ""),
               "Pair test: exact moments \\(k = 20\\) on q = 24, 3 pairs and 1 thread\\.")
  expect_false(grepl("estimated", paste(messages, collapse = "")))
  expect_identical(exact$moments, "exact")
  expect_identical(exact$calibration, "saddlepoint")
  expect_identical(exact$timing$route, list(moments = "exact", k = 20L, q = q))
  contract <- exact$contract
  expect_identical(contract$calibration_contract, "spa_v1")
  expect_identical(contract$route, "exact")
  expect_identical(contract$k, 20L)
  expect_identical(contract$remainder_order, 4L)
  expect_match(contract$basis_sha, "^[0-9a-f]{64}$")
  expect_identical(contract$kernel_version, 2L)
  expect_true(all(exact$results$remainder_kind %in% 0:3))
  # The diagnostic of the compression is one aggregate count in the metadata.
  expect_true(is.numeric(exact$timing$pair_pipeline$nodes_above_leading))
  expect_length(exact$timing$pair_pipeline$nodes_above_leading, 1L)
  expect_false("nodes_above_leading" %in% names(exact$results))

  # k = q is the full-spectrum saddlepoint; the default k = 20 differs from it
  # only by the remainder.
  full <- mgcvST.test(fit, moments = "exact", k = 1000L)
  expect_identical(full$contract$k, q)
  expect_equal(full$results$log_p_two_sided, exact$results$log_p_two_sided, tolerance = 1e-5)
  messages <- testthat::capture_messages(
    pca <- mgcvST.test(fit, moments = "pcalearning", rank = 2L, k = 3L, verbose = TRUE))
  expect_match(paste(messages, collapse = ""), "PCAlearning \\(rank 2, k = 3\\) on q = 24")
  expect_identical(pca$moments, "pcalearning")
  expect_identical(pca$contract$k, 3L)
  expect_identical(pca$contract$route, "pcalearning")
  expect_identical(pca$contract$remainder_order, 2L)
  expect_false(is.null(pca$pca_learning))
  expect_true(is.numeric(pca$timing$pcalearning$nodes_above_leading))
  expect_equal(pca$results$log_p_two_sided, full$results$log_p_two_sided, tolerance = 0.3)
  # Too few genes for the rank: the message points to the exact route.
  expect_error(mgcvST.test(fit, moments = "pcalearning", rank = 4L),
               "achievable PCAlearning rank .*moments = \"exact\"")
})

test_that("a checkpoint directory records its route and a resume with another route stops", {
  fit <- .routing_fit()
  dir <- tempfile("mgcvst-route-")
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  first <- mgcvST.test(fit, checkpoint_dir = dir, chunk_size = 2L, moments = "exact")
  expect_identical(readRDS(file.path(dir, "route.rds"))$moments, "exact")
  expect_error(mgcvST.test(fit, checkpoint_dir = dir, moments = "pcalearning", rank = 2L),
               "written with moments = \"exact\".*asks for moments = \"pcalearning\"")
  # The refusal comes before any work: nothing was added to the directory.
  expect_false(file.exists(file.path(dir, "pca-basis.rds")))
  again <- mgcvST.test(fit, checkpoint_dir = dir, chunk_size = 2L, moments = "exact")
  expect_identical(again$timing$pair_pipeline$resumed_pairs, 3)
  expect_identical(again$results, first$results)

  other <- tempfile("mgcvst-route-pca-")
  on.exit(unlink(other, recursive = TRUE), add = TRUE)
  pca <- mgcvST.test(fit, checkpoint_dir = other, moments = "pcalearning", rank = 2L, k = 3L)
  expect_identical(readRDS(file.path(other, "route.rds"))$moments, "pcalearning")
  expect_error(mgcvST.test(fit, checkpoint_dir = other, moments = "exact"),
               "written with moments = \"pcalearning\".*asks for moments = \"exact\"")
  expect_identical(mgcvST.test(fit, checkpoint_dir = other, moments = "pcalearning",
                               rank = 2L, k = 3L)$results, pca$results)
})

test_that("an INLA fit takes either route and the exact route agrees with the full spectrum", {
  skip_on_cran()
  fit <- .pca_nb_fit()
  pairs <- t(combn(fit$feature_id, 2L))
  reference <- .pca_exact_reference(fit, pairs)
  exact <- inlaST.test(fit, pairs = pairs, threads = 2L, moments = "exact", k = 1000L)
  expect_identical(exact$contract$k, reference$q)
  expect_identical(exact$contract$route, "exact")
  expect_equal(exact$results$score, reference$score, tolerance = 1e-10)
  expect_true(all(exact$results$status == 0L))
  expect_equal(-exact$results$log_p_two_sided / log(10), -reference$log_p / log(10),
               tolerance = 1e-6)
  # The default k = 20 differs from the full spectrum by the remainder only.
  default <- inlaST.test(fit, pairs = pairs, threads = 2L, moments = "exact")
  expect_equal(-default$results$log_p_two_sided / log(10), -reference$log_p / log(10),
               tolerance = 1e-3)
  expect_identical(default$timing$pair_pipeline$preparation_backend, "sparse_reduced")
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
