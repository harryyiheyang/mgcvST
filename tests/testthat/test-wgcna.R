test_that("WGCNA parameters and gene blocks have explicit contracts", {
  skip_if_not_installed("WGCNA")
  skip_if_not_installed("dynamicTreeCut")
  skip_if_not_installed("fastcluster")

  parameters <- spdeST:::.mgcvst_wgcna_parameters
  defaults <- parameters(NULL)
  expect_identical(defaults, parameters(list()))
  expect_identical(parameters(list(power = 4))$power, 4)
  expect_error(parameters(list(powers = 4)), "Unknown wgcna.para")
  expect_error(parameters(list(power = NULL)), "one finite number")
  expect_error(parameters(list(4)), "unique non-empty name")
  expect_error(parameters(list(minClusterSize = 2.5)), "integer of at least 2")

  fit <- .pca_nb_fit()
  expect_error(spdeST.wgcna(fit), "indices must explicitly")
  expect_error(spdeST.wgcna(fit, c(1.1, 2)), "valid integer feature positions")
  expect_error(spdeST.wgcna(fit, c("g1", "missing")), "Unknown feature ID")
  expect_error(spdeST.wgcna(fit, c("g1", "g1")),
               "at least two distinct features")
  expect_error(spdeST.wgcna(fit, list(c(1, 2))), "uniquely named")
})

test_that("WGCNA uses current fitted score states and matches a hand network", {
  skip_if_not_installed("WGCNA")
  skip_if_not_installed("dynamicTreeCut")
  skip_if_not_installed("fastcluster")

  skip_on_cran()
  fit <- .pca_nb_fit()
  before <- serialize(fit, NULL)
  ids <- fit$feature_id[c(3L, 1L, 2L)]
  para <- list(power = 4, minClusterSize = 2L, deepSplit = 0L)
  W <- spdeST.wgcna(fit, ids, wgcna.para = para)

  used <- match(ids, fit$feature_id)
  basis <- spdeST:::.spdest_sparse_observation_basis(spdeST:::.spdest_sparse_prepare(fit))
  A <- crossprod(basis$coordinate, fit$score_a[, used, drop = FALSE])
  dimnames(A) <- list(NULL, ids)
  S <- crossprod(A) / nrow(A)
  dimnames(S) <- list(ids, ids)
  R <- stats::cov2cor(S)
  adj <- WGCNA::adjacency.fromSimilarity(R, type = "signed", power = 4)
  TOM <- WGCNA::TOMsimilarity(adj, TOMType = "signed", verbose = 0)
  dimnames(TOM) <- list(ids, ids)
  H <- fastcluster::hclust(stats::as.dist(1 - TOM), method = "average")
  labels <- as.integer(dynamicTreeCut::cutreeDynamic(
    H, distM = 1 - TOM, minClusterSize = 2L, deepSplit = 0L, verbose = 0
  ))

  expect_identical(serialize(fit, NULL), before)
  expect_equal(W$score$A, A, tolerance = 1e-10)
  expect_equal(W$networks$selected$covariance, S, tolerance = 1e-12)
  expect_equal(W$networks$selected$correlation, R, tolerance = 1e-12)
  expect_equal(W$networks$selected$adjacency, adj, tolerance = 1e-12)
  expect_equal(W$networks$selected$TOM, TOM, tolerance = 1e-12)
  expect_identical(unname(W$networks$selected$labels), labels)
  expect_identical(W$networks$selected$feature_id, ids)
})

test_that("WGCNA preserves overlapping block order", {
  skip_if_not_installed("WGCNA")
  skip_if_not_installed("dynamicTreeCut")
  skip_if_not_installed("fastcluster")
  skip_on_cran()
  skip_if_not_installed("geometry")

  fit <- .pca_nb_fit()
  ids <- fit$feature_id
  blocks <- list(second = ids[c(3L, 1L, 2L)], first = ids[c(2L, 1L, 3L)])

  W <- spdeST.wgcna(fit, blocks)

  prepared <- spdeST:::.spdest_sparse_prepare(fit)
  basis <- spdeST:::.spdest_sparse_observation_basis(prepared)
  used <- match(W$score$feature_id, ids)
  A <- crossprod(basis$coordinate, fit$score_a[, used, drop = FALSE])
  dimnames(A) <- list(NULL, W$score$feature_id)

  expect_identical(names(W$networks), names(blocks))
  expect_identical(W$networks$second$feature_id, blocks$second)
  expect_identical(W$networks$first$feature_id, blocks$first)
  expect_identical(W$score$feature_id, unique(unlist(blocks, use.names = FALSE)))
  expect_equal(W$score$A, A, tolerance = 1e-12)
})

test_that("WGCNA scores are the test's projected observation-kernel scores", {
  skip_on_cran()
  skip_if_not_installed("WGCNA")
  skip_if_not_installed("dynamicTreeCut")
  skip_if_not_installed("fastcluster")
  skip_if_not_installed("geometry")
  set.seed(1731L)
  n <- 60L
  vertices <- as.matrix(expand.grid(
    x = seq(0, 1, length.out = 5L), y = seq(0, 1, length.out = 5L)
  ))
  mesh <- list(loc = vertices, graph = list(tv = geometry::delaunayn(vertices)))
  data <- data.frame(x = runif(n, 0.02, 0.98), y = runif(n, 0.02, 0.98),
                     z = seq(-1, 1, length.out = n))
  basis <- spde_basis(mesh, as.matrix(data[c("x", "y")]), kappa = 1.2,
                      project_intercept = TRUE)
  eta <- 0.3 + 0.2 * data$z + 0.15 * sin(2 * pi * data$x)
  Y <- rbind(
    a = rpois(n, exp(1 + eta)), b = rpois(n, exp(1 + eta)),
    c = rpois(n, exp(1 + eta)), d = rpois(n, exp(1 + eta))
  )
  model <- spdeST.set(response ~ z, data, basis, family = poisson())
  fit <- spdeST:::.spdest_estimate(
    .spdest_scan(
    Y, model, BPPARAM = BiocParallel::SerialParam()),
    Y, features = "all", BPPARAM = BiocParallel::SerialParam())
  ids <- rownames(Y)
  used <- match(ids, fit$feature_id)

  W <- spdeST.wgcna(fit, ids)
  prepared <- spdeST:::.spdest_sparse_prepare(fit)
  basis_r <- spdeST:::.spdest_sparse_observation_basis(prepared)
  A <- crossprod(basis_r$coordinate, fit$score_a[, used, drop = FALSE])
  dimnames(A) <- list(NULL, ids)
  expect_equal(W$score$A, A, tolerance = 1e-12)

  pairs <- t(utils::combn(ids, 2L))
  T <- spdeST.test(fit, pairs = pairs, rank = 1L, moments = "exact")
  i <- match(fit$feature_id[T$results$i], ids)
  j <- match(fit$feature_id[T$results$j], ids)
  expect_equal(T$results$score,
               unname(colSums(A[, i, drop = FALSE] * A[, j, drop = FALSE])),
               tolerance = 1e-10)
})
