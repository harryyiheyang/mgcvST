.wgcna_test <- local({
  cached <- NULL
  function() {
    if (is.null(cached)) cached <<- spdeST.test(.pca_nb_fit(), moments = "exact")
    cached
  }
})

test_that("the similarity of two genes is U / t1 of their pair", {
  skip_on_cran()
  test <- .wgcna_test()
  r <- test$results
  used <- c(3L, 1L, 6L, 2L)
  pairs <- spdeST:::.mgcvst_wgcna_pairs(test, used)
  expect_identical(nrow(pairs), 6L)
  position <- integer(length(test$feature_id))
  position[used] <- seq_along(used)
  z <- spdeST:::mgcvst_he_correlation_cpp(pairs$i, pairs$j, pairs$score, pairs$t1,
                                          position, length(used))
  expect_identical(diag(z$rho), rep(1, 4L))
  expect_identical(z$missing, 0)
  for (a in 1:3) for (b in (a + 1L):4L) {
    row <- pair_row(test, min(used[a], used[b]), max(used[a], used[b]))
    expect_identical(z$rho[a, b], row$score / row$t1)
    expect_identical(z$rho[b, a], z$rho[a, b])
  }
  # The shards hold the same pairs as the in-memory table.
  from_shards <- spdeST.test(.pca_nb_fit(), moments = "exact",
                             checkpoint_dir = withr::local_tempdir())
  expect_true(length(from_shards$shards) > 0L)
  from_shards$results <- NULL
  shard_pairs <- spdeST:::.mgcvst_wgcna_pairs(from_shards, used)
  expect_equal(shard_pairs[order(shard_pairs$i, shard_pairs$j), ],
               pairs[order(pairs$i, pairs$j), ], ignore_attr = TRUE)

  # An estimate outside [-1, 1] enters the similarity at the limit and is
  # counted; a pair without information is missing.
  z <- spdeST:::mgcvst_he_correlation_cpp(
    c(1L, 1L, 2L), c(2L, 3L, 3L), c(3, -0.5, -4), c(2, 1, 2), 1:3, 3L)
  expect_identical(z$rho[1, 2], 1.5)
  expect_identical(z$similarity[1, 2], 1)
  expect_identical(z$similarity[3, 2], -1)
  expect_identical(z$similarity[1, 3], -0.5)
  expect_identical(z$limited, 2)
  z <- spdeST:::mgcvst_he_correlation_cpp(
    c(1L, 1L), c(2L, 3L), c(1, NA), c(2, 1), 1:3, 3L)
  expect_identical(z$missing, 2)
  expect_error(spdeST:::mgcvst_he_correlation_cpp(1L, 9L, 1, 1, 1:3, 3L),
               "outside the fit")
})

test_that("WGCNA parameters and gene blocks have explicit contracts", {
  skip_if_not_installed("WGCNA")
  skip_if_not_installed("dynamicTreeCut")
  skip_if_not_installed("fastcluster")
  skip_on_cran()

  parameters <- spdeST:::.mgcvst_wgcna_parameters
  defaults <- parameters(NULL)
  expect_identical(defaults, parameters(list()))
  expect_identical(parameters(list(power = 4))$power, 4)
  expect_error(parameters(list(powers = 4)), "Unknown wgcna.para")
  expect_error(parameters(list(power = NULL)), "one finite number")
  expect_error(parameters(list(4)), "unique non-empty name")
  expect_error(parameters(list(minClusterSize = 2.5)), "integer of at least 2")

  test <- .wgcna_test()
  expect_error(spdeST.wgcna(.pca_nb_fit(), 1:3), "result of spdeST.test")
  expect_error(spdeST.wgcna(test), "indices must explicitly")
  expect_error(spdeST.wgcna(test, c(1.1, 2)), "valid integer feature positions")
  expect_error(spdeST.wgcna(test, c("g1", "missing")), "Unknown feature ID")
  expect_error(spdeST.wgcna(test, c("g1", "g1")),
               "at least two distinct features")
  expect_error(spdeST.wgcna(test, list(c(1, 2))), "uniquely named")
  # Every pair of a block must have been tested.
  part <- spdeST.test(.pca_nb_fit(), pairs = rbind(c("g1", "g2"), c("g1", "g3")),
                      moments = "exact")
  expect_error(spdeST.wgcna(part, c("g1", "g2", "g3")), "1 pairs of the selected features")
})

test_that("the network of a block is built from U / t1 and matches a hand network", {
  skip_if_not_installed("WGCNA")
  skip_if_not_installed("dynamicTreeCut")
  skip_if_not_installed("fastcluster")
  skip_on_cran()

  test <- .wgcna_test()
  before <- serialize(test, NULL)
  ids <- test$feature_id[c(3L, 1L, 2L, 5L)]
  para <- list(power = 4, minClusterSize = 2L, deepSplit = 0L)
  W <- spdeST.wgcna(test, ids, wgcna.para = para)

  used <- match(ids, test$feature_id)
  R <- diag(4L)
  for (a in 1:3) for (b in (a + 1L):4L) {
    row <- pair_row(test, min(used[a], used[b]), max(used[a], used[b]))
    R[a, b] <- R[b, a] <- row$score / row$t1
  }
  dimnames(R) <- list(ids, ids)
  similarity <- pmax(pmin(R, 1), -1)
  adj <- WGCNA::adjacency.fromSimilarity(similarity, type = "signed", power = 4)
  TOM <- WGCNA::TOMsimilarity(adj, TOMType = "signed", verbose = 0)
  dimnames(adj) <- dimnames(TOM) <- list(ids, ids)
  H <- fastcluster::hclust(stats::as.dist(1 - TOM), method = "average")
  labels <- as.integer(dynamicTreeCut::cutreeDynamic(
    H, distM = 1 - TOM, minClusterSize = 2L, deepSplit = 0L, verbose = 0
  ))

  expect_identical(serialize(test, NULL), before)
  expect_identical(W$networks$selected$correlation, R)
  expect_identical(W$networks$selected$limited, sum(abs(R[upper.tri(R)]) > 1) + 0)
  expect_equal(W$networks$selected$adjacency, adj, tolerance = 1e-12)
  expect_equal(W$networks$selected$TOM, TOM, tolerance = 1e-12)
  expect_identical(unname(W$networks$selected$labels), labels)
  expect_identical(W$networks$selected$feature_id, ids)
  expect_null(W$networks$selected$covariance)
  expect_null(W$score)
  expect_output(print(W), "spdeST WGCNA modules")
})

test_that("WGCNA preserves overlapping block order", {
  skip_if_not_installed("WGCNA")
  skip_if_not_installed("dynamicTreeCut")
  skip_if_not_installed("fastcluster")
  skip_on_cran()

  test <- .wgcna_test()
  ids <- test$feature_id
  blocks <- list(second = ids[c(3L, 1L, 2L)], first = ids[c(2L, 1L, 3L)])
  W <- spdeST.wgcna(test, blocks)
  expect_identical(names(W$networks), names(blocks))
  expect_identical(W$networks$second$feature_id, blocks$second)
  expect_identical(W$networks$first$feature_id, blocks$first)
  expect_identical(W$networks$second$correlation[blocks$first, blocks$first],
                   W$networks$first$correlation)
})
