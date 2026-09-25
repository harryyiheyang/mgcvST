test_that("WGCNA parameters and gene blocks have explicit contracts", {
  skip_if_not_installed("WGCNA")
  skip_if_not_installed("dynamicTreeCut")
  skip_if_not_installed("fastcluster")

  parameters <- mgcvST:::.mgcvst_wgcna_parameters
  defaults <- parameters(NULL)
  expect_identical(defaults, parameters(list()))
  expect_identical(parameters(list(power = 4))$power, 4)
  expect_error(parameters(list(powers = 4)), "Unknown wgcna.para")
  expect_error(parameters(list(power = NULL)), "one finite number")
  expect_error(parameters(list(4)), "unique non-empty name")
  expect_error(parameters(list(minClusterSize = 2.5)), "integer of at least 2")

  f <- st_fixture(n = 45L, family = gaussian())
  fit <- mgcvST.estimate(
    f$Y, f$model, diagnostics = FALSE,
    BPPARAM = BiocParallel::SerialParam()
  )
  expect_error(mgcvST.wgcna(fit), "indices must explicitly")
  expect_error(mgcvST.wgcna(fit, c(1.1, 2)), "valid integer feature positions")
  expect_error(mgcvST.wgcna(fit, c("response", "missing")), "Unknown feature ID")
  expect_error(mgcvST.wgcna(fit, c("response", "response")),
               "at least two distinct features")
  expect_error(mgcvST.wgcna(fit, list(c(1, 2))), "uniquely named")
})

test_that("WGCNA uses current fitted score states and matches a hand network", {
  skip_if_not_installed("WGCNA")
  skip_if_not_installed("dynamicTreeCut")
  skip_if_not_installed("fastcluster")

  f <- st_fixture(n = 54L, family = gaussian(), nuisance = TRUE)
  fit <- mgcvST.estimate(
    f$Y, f$model, diagnostics = FALSE,
    BPPARAM = BiocParallel::SerialParam()
  )
  before <- serialize(fit, NULL)
  ids <- fit$feature_id[c(3L, 1L, 2L)]
  para <- list(power = 4, minClusterSize = 2L, deepSplit = 0L)
  W <- mgcvST.wgcna(fit, ids, wgcna.para = para)

  A <- vapply(match(ids, fit$feature_id), function(i) {
    mgcvST:::.mgcvst_model_score_state(fit, i)$a
  }, numeric(nrow(W$score$A)))
  colnames(A) <- ids
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

test_that("mgcv WGCNA scores agree with mgcvST.test()'s pairwise scores", {
  skip_if_not_installed("WGCNA")
  skip_if_not_installed("dynamicTreeCut")
  skip_if_not_installed("fastcluster")

  f <- st_fixture(n = 54L, family = gaussian(), nuisance = TRUE)
  fit <- mgcvST.estimate(
    f$Y, f$model, diagnostics = FALSE,
    BPPARAM = BiocParallel::SerialParam()
  )
  ids <- fit$feature_id[c(3L, 1L, 2L)]
  W <- mgcvST.wgcna(fit, ids)
  pairs <- t(utils::combn(ids, 2L))
  T <- mgcvST.test(fit, pairs = pairs, calibration = "liu")
  i <- match(T$results$feature1, colnames(W$score$A))
  j <- match(T$results$feature2, colnames(W$score$A))
  G <- crossprod(W$score$A)
  expect_equal(T$results$signed_score, G[cbind(i, j)], tolerance = 1e-10)
})

test_that("WGCNA preserves overlapping block order", {
  skip_if_not_installed("WGCNA")
  skip_if_not_installed("dynamicTreeCut")
  skip_if_not_installed("fastcluster")
  skip_on_cran()
  skip_if_not_installed("geometry")

  f <- st_fixture(family = gaussian())
  fit <- mgcvST.estimate(f$Y, f$G, diagnostics = FALSE,
                         BPPARAM = BiocParallel::SerialParam())
  ids <- fit$feature_id
  blocks <- list(second = ids[c(3L, 1L, 2L)], first = ids[c(2L, 1L, 3L)])

  W <- mgcvST.wgcna(fit, blocks)

  used <- match(W$score$feature_id, ids)
  T0 <- mgcvST:::.mgcvst_legacy_shared_score_factor(fit$geometry)
  field_scale <- mgcvST:::.mgcvst_field_scale(fit)
  z <- mgcvST:::mgcvst_dense_score_batch_cpp(
    T0, fit$working_variance[, used, drop = FALSE],
    fit$working_error[, used, drop = FALSE], field_scale[used],
    fit$geometry$X, list(), 1L, score_only = TRUE
  )
  A <- vapply(z, `[[`, numeric(ncol(T0)), "a")
  colnames(A) <- W$score$feature_id

  expect_identical(names(W$networks), names(blocks))
  expect_identical(W$networks$second$feature_id, blocks$second)
  expect_identical(W$networks$first$feature_id, blocks$first)
  expect_identical(W$score$feature_id, unique(unlist(blocks, use.names = FALSE)))
  expect_identical(W$score$group, "global")
  expect_identical(unname(W$score$width), ncol(T0))
  expect_equal(W$score$A, A, tolerance = 1e-10)
})

test_that("WGCNA retains legacy compact B-Q-X score semantics", {
  skip_if_not_installed("WGCNA")
  skip_if_not_installed("dynamicTreeCut")
  skip_if_not_installed("fastcluster")

  set.seed(811L)
  n <- 18L
  p <- 4L
  B <- matrix(rnorm(n * 3L), n, 3L)
  Q <- diag(c(1, 2, 3))
  X <- matrix(1, n, 1L)
  ids <- paste0("old", seq_len(p))
  fit <- structure(list(
    feature_id = ids,
    working_error = matrix(rnorm(n * p), n, p),
    working_variance = matrix(runif(n * p, 0.8, 1.3), n, p),
    dispersion = seq(0.9, 1.2, length.out = p),
    lambda = seq(1, 1.6, length.out = p),
    geometry = list(B = B, Q = Q, X = X)
  ), class = c("mgcvST_fit", "mgcvST"))
  colnames(fit$working_error) <- colnames(fit$working_variance) <- ids
  selected <- c("old4", "old2", "old1")
  W <- mgcvST.wgcna(fit, selected)
  A <- vapply(match(selected, ids), function(i) {
    op <- rkhs_score_operator(
      B, Q, fit$working_variance[, i], X,
      field_scale = fit$dispersion[i] / fit$lambda[i]
    )
    rkhs_score_summary(fit$working_error[, i], op)$a
  }, numeric(ncol(B)))
  colnames(A) <- selected

  expect_equal(W$score$A, A, tolerance = 1e-10)
  expect_identical(W$score$group, "global")
  expect_identical(W$score$width, c(global = ncol(B)))
})

test_that("inlaST.wgcna rejects an mgcv fit and names the right entry point", {
  skip_on_cran()
  skip_if_not_installed("geometry")
  f <- st_fixture()
  fit <- mgcvST.estimate(f$Y, f$G)
  expect_error(
    inlaST.wgcna(fit, rownames(f$Y)),
    "requires a fit returned by inlaST.estimate"
  )
})

test_that("INLA WGCNA scores are the test's projected observation-kernel scores", {
  skip_on_cran()
  skip_if_not_installed("INLA")
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
    a = eta + rnorm(n, sd = 0.3), b = eta + rnorm(n, sd = 0.3),
    c = eta + rnorm(n, sd = 0.3), d = eta + rnorm(n, sd = 0.3)
  )
  model <- inlaST.set(response ~ z, data, basis, family = gaussian())
  fit <- inlaST.estimate(
    Y, model, BPPARAM = BiocParallel::SerialParam(),
    control = list(fixed_precision = 1.7, gaussian_precision = 1 / 0.09)
  )
  ids <- rownames(Y)
  used <- match(ids, fit$feature_id)

  W <- inlaST.wgcna(fit, ids)
  prepared <- mgcvST:::.inlast_sparse_prepare(fit)
  basis_r <- mgcvST:::.inlast_sparse_observation_basis(prepared)
  A <- crossprod(basis_r$coordinate, fit$score_a[, used, drop = FALSE])
  dimnames(A) <- list(NULL, ids)
  expect_equal(W$score$A, A, tolerance = 1e-12)

  pairs <- t(utils::combn(ids, 2L))
  T <- inlaST.test(fit, pairs = pairs, calibration = "liu")
  i <- match(fit$feature_id[T$result$i], ids)
  j <- match(fit$feature_id[T$result$j], ids)
  expect_equal(T$result$score,
               unname(colSums(A[, i, drop = FALSE] * A[, j, drop = FALSE])),
               tolerance = 1e-10)

  expect_error(
    mgcvST.test(fit),
    "mgcvST.test\\(\\) does not accept inlaST.estimate\\(\\) fits; use inlaST.test\\(\\)."
  )
  expect_error(
    mgcvST.wgcna(fit, ids),
    "mgcvST.wgcna\\(\\) does not accept inlaST.estimate\\(\\) fits; use inlaST.wgcna\\(\\)."
  )
})
