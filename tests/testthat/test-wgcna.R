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

test_that("WGCNA preserves overlapping block order", {
  skip_if_not_installed("WGCNA")
  skip_if_not_installed("dynamicTreeCut")
  skip_if_not_installed("fastcluster")

  set.seed(719L)
  n <- 16L
  p <- 6L
  ids <- paste0("g", seq_len(p))
  global <- list(
    B = matrix(rnorm(n * 3L), n, 3L), fixed = FALSE,
    score_component = "global", penalties = list(diag(c(1, 2, 3))),
    sp_index = 1L
  )
  nuisance <- list(
    B = cbind(1, seq_len(n) / n), fixed = FALSE,
    score_component = NULL, penalties = list(diag(c(0, 2))), sp_index = 2L
  )
  fit <- structure(list(
    feature_id = ids,
    working_error = matrix(rnorm(n * p), n, p),
    working_variance = matrix(runif(n * p, 0.7, 1.4), n, p),
    dispersion = seq(0.8, 1.2, length.out = p),
    smoothing_parameters = cbind(rep(1.1, p), rep(0.8, p)),
    geometry = list(
      X = matrix(1, n, 1L), smooth = list(global, nuisance),
      target = list(global = 1L)
    ),
    score_components = "global"
  ), class = c("mgcvST_model_fit", "mgcvST_fit", "mgcvST"))
  colnames(fit$working_error) <- colnames(fit$working_variance) <- ids
  rownames(fit$smoothing_parameters) <- ids
  colnames(fit$smoothing_parameters) <- c("global", "s(z)")
  fit$.mgcvst_fixed_factors <- mgcvST:::.mgcvst_model_fixed_factors(fit)
  blocks <- list(second = c("g6", "g2", "g4"),
                 first = c("g4", "g1", "g2"))

  W <- mgcvST.wgcna(fit, blocks)
  states <- lapply(match(W$score$feature_id, ids), function(i) {
    mgcvST:::.mgcvst_model_score_state(fit, i)
  })
  A <- vapply(states, function(z) z$a[1:3], numeric(3L))
  colnames(A) <- W$score$feature_id

  expect_identical(names(W$networks), names(blocks))
  expect_identical(W$networks$second$feature_id, blocks$second)
  expect_identical(W$networks$first$feature_id, blocks$first)
  expect_identical(W$score$feature_id, unique(unlist(blocks, use.names = FALSE)))
  expect_identical(W$score$group, "global")
  expect_identical(W$score$width, c(global = 3L))
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
