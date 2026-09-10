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

test_that("WGCNA consumes conditioned sparse INLA score states", {
  skip_if_not_installed("WGCNA")
  skip_if_not_installed("dynamicTreeCut")
  skip_if_not_installed("fastcluster")

  set.seed(9021)
  n <- 18L
  p <- 3L
  g <- rep(1 / n, n)
  Z <- qr.Q(qr(matrix(g, ncol = 1L)), complete = TRUE)[, -1L, drop = FALSE]
  D <- seq(0.7, 1.4, length.out = n)
  X <- matrix(1, n, 1L, dimnames = list(NULL, "(Intercept)"))
  base <- mgcvST:::.rkhs_score_operator_factor(
    Z, D, matrix(numeric(), n, 0L), field_scale = 1, B = NULL, Q = NULL
  )
  WX <- mgcvST:::.mgcvst_model_vsolve(base, X)
  Vp <- 0.973 * solve(crossprod(X, WX))
  ids <- paste0("gene", seq_len(p))
  fit <- structure(list(
    score_backend = "sparse",
    score_sparse = list(
      A = Matrix::Diagonal(n), Q = Matrix::Diagonal(n), constraint = g,
      projection = Z, coefficient_factor = Z, target = "global", sp_index = 1L
    ),
    feature_id = ids,
    dispersion = stats::setNames(rep(1, p), ids),
    smoothing_parameters = matrix(1, p, 1L,
      dimnames = list(ids, "global")),
    working_error = sapply(seq_len(p), function(j) {
      sin(seq_len(n) / (j + 1)) + seq_len(n) / (40 + j)
    }),
    working_variance = matrix(D, n, p),
    nuisance_covariance = stats::setNames(rep(list(Vp), p), ids),
    geometry = list(
      nuisance_design = X, X = X, target = c(global = 1L),
      smooth = list(list(
        fixed = FALSE, score_component = "global", B = Z,
        penalties = list(diag(ncol(Z))), sp_index = 1L
      ))
    )
  ), class = c("inlaST_fit", "mgcvST_model_fit", "mgcvST_fit", "mgcvST"))
  dimnames(fit$working_error) <- dimnames(fit$working_variance) <- list(NULL, ids)

  expected <- vapply(seq_len(p), function(i) {
    mgcvST:::.mgcvst_model_score_state(fit, i)$a
  }, numeric(ncol(Z)))
  colnames(expected) <- ids
  for (i in seq_len(p)) {
    full <- mgcvST:::.mgcvst_model_sparse_score_state(fit, i)
    score_only <- mgcvST:::.mgcvst_model_sparse_score_state(
      fit, i, score_only = TRUE
    )
    expect_equal(score_only$a, full$a, tolerance = 2e-12)
    expect_identical(score_only$width, full$width)
    expect_null(score_only$M)
  }
  testthat::local_mocked_bindings(
    .mgcvst_model_operator = function(...) {
      stop("dense score operator was called")
    },
    .package = "mgcvST"
  )
  W <- mgcvST.wgcna(
    fit, ids, wgcna.para = list(minClusterSize = 2L, deepSplit = 0L)
  )

  expect_equal(W$score$A, expected, tolerance = 2e-12)
  expect_identical(W$score$feature_id, ids)
  expect_identical(W$score$width, c(global = ncol(Z)))
})

test_that("WGCNA preserves overlapping block and score-group order", {
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
  local <- list(
    B = matrix(rnorm(n * 2L), n, 2L), fixed = FALSE,
    score_component = "local", penalties = list(diag(c(1, 2))),
    sp_index = 2L
  )
  nuisance <- list(
    B = cbind(1, seq_len(n) / n), fixed = FALSE,
    score_component = NULL, penalties = list(diag(c(0, 2))), sp_index = 3L
  )
  fit <- structure(list(
    feature_id = ids,
    working_error = matrix(rnorm(n * p), n, p),
    working_variance = matrix(runif(n * p, 0.7, 1.4), n, p),
    dispersion = seq(0.8, 1.2, length.out = p),
    smoothing_parameters = cbind(rep(1.1, p), rep(1.3, p), rep(0.8, p)),
    geometry = list(
      X = matrix(1, n, 1L), smooth = list(global, local, nuisance),
      target = list(global = 1L, local = 2L)
    ),
    score_components = c("global", "local")
  ), class = c("mgcvST_model_fit", "mgcvST_fit", "mgcvST"))
  colnames(fit$working_error) <- colnames(fit$working_variance) <- ids
  rownames(fit$smoothing_parameters) <- ids
  colnames(fit$smoothing_parameters) <- c("global", "local", "s(z)")
  fit$.mgcvst_fixed_factors <- mgcvST:::.mgcvst_model_fixed_factors(fit)
  blocks <- list(second = c("g6", "g2", "g4"),
                 first = c("g4", "g1", "g2"))

  expect_error(mgcvST.wgcna(fit, blocks), "group must explicitly select")
  W <- mgcvST.wgcna(fit, blocks, group = c("local", "global"))
  states <- lapply(match(W$score$feature_id, ids), function(i) {
    mgcvST:::.mgcvst_model_score_state(fit, i)
  })
  A <- vapply(states, function(z) c(z$a[4:5], z$a[1:3]), numeric(5L))
  colnames(A) <- W$score$feature_id

  expect_identical(names(W$networks), names(blocks))
  expect_identical(W$networks$second$feature_id, blocks$second)
  expect_identical(W$networks$first$feature_id, blocks$first)
  expect_identical(W$score$feature_id, unique(unlist(blocks, use.names = FALSE)))
  expect_identical(W$score$group, c("local", "global"))
  expect_identical(W$score$width, c(local = 2L, global = 3L))
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

test_that("WGCNA runs after a real sparse INLA estimate", {
  skip_if_not_installed("INLA")
  skip_if_not_installed("geometry")
  skip_if_not_installed("WGCNA")
  skip_if_not_installed("dynamicTreeCut")
  skip_if_not_installed("fastcluster")

  set.seed(1741L)
  n <- 36L
  vertices <- as.matrix(expand.grid(
    x = seq(0, 1, length.out = 4L), y = seq(0, 1, length.out = 4L)
  ))
  mesh <- list(loc = vertices, graph = list(tv = geometry::delaunayn(vertices)))
  dat <- data.frame(
    x = runif(n, 0.02, 0.98), y = runif(n, 0.02, 0.98),
    z = seq(-1, 1, length.out = n), offset0 = runif(n, -0.1, 0.1)
  )
  basis <- spde_basis(
    mesh, as.matrix(dat[c("x", "y")]), kappa = 1.1,
    project_intercept = TRUE
  )
  eta <- 0.4 + 0.2 * dat$z + dat$offset0
  Y <- rbind(
    a = eta + rnorm(n, sd = 0.3),
    b = eta + 0.1 * sin(2 * pi * dat$x) + rnorm(n, sd = 0.3),
    c = eta - 0.1 * cos(2 * pi * dat$y) + rnorm(n, sd = 0.3)
  )
  model <- inlaST.set(
    response ~ z + offset(offset0), dat, basis, family = gaussian(),
    score_backend = "sparse"
  )
  fit <- inlaST.estimate(
    Y, model, score_backend = "sparse", diagnostics = FALSE,
    BPPARAM = BiocParallel::SerialParam(),
    control = list(fixed_precision = 2, gaussian_precision = 1 / 0.09)
  )
  W <- mgcvST.wgcna(
    fit, rownames(Y),
    wgcna.para = list(minClusterSize = 2L, deepSplit = 0L)
  )
  A <- vapply(seq_len(nrow(Y)), function(i) {
    mgcvST:::.mgcvst_model_score_state(fit, i)$a
  }, numeric(nrow(W$score$A)))
  colnames(A) <- rownames(Y)

  expect_s3_class(fit, "inlaST_fit")
  expect_identical(fit$score_backend, "sparse")
  expect_equal(W$score$A, A, tolerance = 1e-9)
  expect_true(all(is.finite(W$networks$selected$TOM)))
})
