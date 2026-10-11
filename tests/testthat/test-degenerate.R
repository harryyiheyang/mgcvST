# Degenerate spatial fits: a gene whose spatial smooth has (almost) no effective
# degrees of freedom has p = 1 in every pair (status 4); a pair the kernel cannot
# evaluate also has p = 1 and stays in the adjustment family.

.deg_field <- function(basis, v, seed) {
  set.seed(seed)
  f <- as.numeric(basis$B %*% rnorm(ncol(basis$B)))
  f <- f - mean(f)
  f * sqrt(v / stats::var(f))
}

# Eleven genes on a 6 x 6 mesh with 1500 spots: strong, medium and weak (v = 0.03)
# fields, and eight pure-noise genes. The spatial smooth of a pure-noise gene is
# penalized to the boundary or to a small edf in most fits; the tests use the
# genes that this fit sends to the boundary (`none4`, `none7`, `none8`).
.deg_mgcv_ids <- c("strong", "medium", "weak", "none4", "none7")
.deg_mgcv_pairs <- function() t(combn(.deg_mgcv_ids, 2L))
.deg_mgcv_fit <- local({
  cached <- NULL
  function() {
    skip_if_not_installed("geometry")
    if (!is.null(cached)) return(cached)
    set.seed(5)
    n <- 1500L
    xy <- as.matrix(expand.grid(x = seq(0, 1, length.out = 6), y = seq(0, 1, length.out = 6)))
    mesh <- list(loc = xy, graph = list(tv = geometry::delaunayn(xy)))
    data <- data.frame(x = runif(n, .01, .99), y = runif(n, .01, .99),
                       offset0 = runif(n, -.2, .2), z = runif(n))
    basis <- spde_basis(mesh, as.matrix(data[, c("x", "y")]), kappa = 1.5,
                        project_intercept = TRUE)
    mu0 <- exp(1.2 + data$offset0)
    rnb <- function(mu) rnbinom(length(mu), mu = mu, size = 6)
    set.seed(11)
    Y <- rbind(strong = rnb(mu0 * exp(.deg_field(basis, 0.30, 1))),
               medium = rnb(mu0 * exp(.deg_field(basis, 0.10, 2))),
               weak = rnb(mu0 * exp(.deg_field(basis, 0.03, 3))),
               none1 = rnb(mu0), none2 = rnb(mu0), none3 = rnb(mu0), none4 = rnb(mu0),
               none5 = rnb(mu0), none6 = rnb(mu0), none7 = rnb(mu0), none8 = rnb(mu0))
    model <- model.set(response ~ offset(offset0) + z, cbind(response = 1, data), basis,
                       family = mgcv::nb())
    cached <<- suppressWarnings(mgcvST.estimate(
      Y, model, BPPARAM = BiocParallel::SerialParam(), spatial = "all"))
    cached
  }
})

test_that("the effective degrees of freedom flag a pure-noise gene and not a weak field", {
  skip_on_cran()
  fit <- .deg_mgcv_fit()
  d <- fit$diagnostics
  rownames(d) <- d$feature_id
  expect_identical(fit$format, 3L)
  expect_true(is.numeric(mgcvST:::.mgcvst_edf_min) && length(mgcvST:::.mgcvst_edf_min) == 1L)
  expect_identical(d$spatial_degenerate, d$edf_spatial < mgcvST:::.mgcvst_edf_min)
  # The pure-noise genes sit at the boundary; the weak field (v = 0.03) keeps
  # several effective degrees of freedom, above the minimum.
  expect_lt(max(d[c("none4", "none7", "none8"), "edf_spatial"]), 0.01)
  expect_true(all(d[c("none4", "none7"), "spatial_degenerate"]))
  expect_false(any(d[c("strong", "medium", "weak"), "spatial_degenerate"]))
  expect_gt(d["weak", "edf_spatial"], 2 * mgcvST:::.mgcvst_edf_min)
  expect_true(all(diff(d[c("strong", "medium", "weak"), "edf_spatial"]) < 0))
  expect_identical(mgcvST:::.mgcvst_degenerate_features(fit), d$spatial_degenerate)
})

test_that("a pair with a degenerate gene has p = 1 and status 4 in both routes", {
  skip_on_cran()
  fit <- .deg_mgcv_fit()
  P <- .deg_mgcv_pairs()
  routes <- list(
    exact = mgcvST.test(fit, pairs = P, moments = "exact", adjust = "BY"),
    pcalearning = mgcvST.test(fit, pairs = P, moments = "pcalearning", rank = 2L, k = 3L,
                              adjust = "BY"))
  noise <- match(c("none4", "none7"), fit$feature_id)
  for (route in names(routes)) {
    z <- routes[[route]]
    r <- z$results
    deg <- r$i %in% noise | r$j %in% noise
    expect_identical(sum(deg), 7L)
    expect_true(all(r$status[deg] == 4L), info = route)
    expect_true(all(r$log_p_two_sided[deg] == 0 & r$log_p_positive[deg] == 0 &
                      r$log_p_negative[deg] == 0), info = route)
    expect_true(all(is.na(r$score[deg])), info = route)
    expect_true(all(r$remainder_kind[deg] == 0L), info = route)
    expect_true(all(r$status[!deg] == 0L), info = route)
    expect_true(all(is.finite(r$log_p_two_sided[!deg]) & r$log_p_two_sided[!deg] < 0), info = route)
    expect_true(all(c("none4", "none7") %in% z$degenerate), info = route)
    expect_true(all(z$degenerate %in% fit$feature_id[fit$diagnostics$spatial_degenerate]))
    # The degenerate pairs stay in the family: all ten pairs are adjusted and
    # BY equals stats::p.adjust on the p-values including the ones.
    expect_equal(z$adjustment$n_adjusted, 10, info = route)
    expect_equal(exp(r$log_q), stats::p.adjust(exp(r$log_p_two_sided), "BY"),
                 tolerance = 1e-12, info = route)
    expect_equal(z$discoveries$pairs_tested, 10)
    expect_equal(z$discoveries$pairs_with_p_value, 10)
    expect_identical(z$contract$kernel_version, 3L)
  }
  # The pairs of the other genes do not depend on the degenerate genes: the same
  # three p-values from a test restricted to them.
  ids <- fit$feature_id
  only <- mgcvST.test(fit, pairs = rbind(c("strong", "medium"), c("strong", "weak"),
                                         c("medium", "weak")), moments = "exact")
  full <- routes$exact$results
  expect_equal(only$results$log_p_two_sided,
               full$log_p_two_sided[full$status == 0L], tolerance = 1e-9)
})

test_that("a test whose genes are all degenerate gives p = 1 for every pair", {
  skip_on_cran()
  fit <- .deg_mgcv_fit()
  z <- mgcvST.test(fit, pairs = rbind(c("none4", "none7")), moments = "exact")
  expect_identical(z$results$status, 4L)
  expect_identical(z$results$log_p_two_sided, 0)
  expect_match(z$contract$basis_sha, "^[0-9a-f]{64}$")
  expect_true(all(c("none4", "none7") %in% z$degenerate))
  # One degenerate and one usable gene.
  mixed <- mgcvST.test(fit, pairs = rbind(c("strong", "none4"), c("strong", "medium")),
                       moments = "exact")
  expect_identical(mixed$results$status[order(mixed$results$j)], c(0L, 4L))
  expect_identical(mixed$results$log_p_negative[mixed$results$status == 4L], 0)
})

test_that("a pair the kernel cannot evaluate has p = 1 and stays in the family", {
  # exact route: the kernel returns the invalid statuses 2 and 1
  testthat::local_mocked_bindings(
    mgcvst_pair_spa_cpp = function(H, G, avec, left, right, threads, order) {
      list(score = c(1, 2, 3), log_p_two_sided = c(NA, -3, -Inf),
           log_p_positive = c(NA, -3.7, -Inf), log_p_negative = c(NA, -0.02, -Inf),
           remainder_kind = c(0L, 2L, 0L), status = c(2L, 0L, 1L),
           nodes_above_leading = 0)
    }, .package = "mgcvST")
  states <- lapply(1:3, function(i) list(a = c(i, 1), M = diag(2)))
  out <- mgcvST:::.mgcvst_spa_pairs(rbind(c(1L, 2L), c(1L, 3L), c(2L, 3L)), 1:3, states,
                                    vector("list", 3L), 1L)
  expect_identical(out$status, c(2L, 0L, 1L))
  expect_identical(out$log_p_two_sided, c(0, -3, 0))
  expect_identical(out$log_p_positive, c(0, -3.7, 0))
  expect_identical(out$log_p_negative, c(0, -0.02, 0))
  expect_equal(out$score, c(1, 2, 3))
  adjusted <- mgcvST:::.mgcvst_log_adjust(out$log_p_two_sided, "BY")
  expect_equal(adjusted$n, 3)
  expect_equal(exp(adjusted$log_q), stats::p.adjust(exp(out$log_p_two_sided), "BY"),
               tolerance = 1e-12)
  # a failed state is still status 3 without a p-value, whatever the other gene
  states[[2L]] <- list(error = "state failed")
  bad <- mgcvST:::.mgcvst_spa_pairs(rbind(c(1L, 2L), c(1L, 3L)), 1:3, states,
                                    vector("list", 3L), 1L,
                                    degenerate = c(FALSE, FALSE, TRUE))
  expect_identical(bad$status, c(3L, 4L))
  expect_true(is.na(bad$log_p_two_sided[1L]))
  expect_identical(bad$log_p_two_sided[2L], 0)
})

test_that("a PCAlearning pair the kernel cannot evaluate has p = 1", {
  skip_on_cran()
  fit <- .deg_mgcv_fit()
  pairs_kernel <- mgcvST:::mgcvst_pca_spa_pairs_cpp
  testthat::local_mocked_bindings(
    mgcvst_pca_spa_pairs_cpp = function(...) {
      out <- pairs_kernel(...)
      out[1L, c("logp_two_sided", "logp_positive", "logp_negative")] <- NA_real_
      out[1L, "status"] <- 2
      out
    }, .package = "mgcvST")
  z <- mgcvST.test(fit, pairs = .deg_mgcv_pairs(), moments = "pcalearning", rank = 2L, k = 3L,
                   adjust = "BY")
  r <- z$results
  expect_identical(sum(r$status == 2L), 1L)
  bad <- r$status == 2L
  expect_identical(c(r$log_p_two_sided[bad], r$log_p_positive[bad], r$log_p_negative[bad]),
                   c(0, 0, 0))
  expect_true(is.finite(r$score[bad]))
  expect_equal(z$adjustment$n_adjusted, 10)
  expect_equal(exp(r$log_q), stats::p.adjust(exp(r$log_p_two_sided), "BY"), tolerance = 1e-12)
})

test_that("fits and checkpoints written before the edf diagnostics are refused", {
  skip_on_cran()
  fit <- .deg_mgcv_fit()
  old <- fit
  old$format <- 2L
  old$diagnostics$edf_spatial <- NULL
  old$diagnostics$spatial_degenerate <- NULL
  expect_error(mgcvST.test(old, moments = "exact"),
               "estimated before mgcvST 0.0.1.9034.*effective degrees of freedom")
  expect_error(mgcvST.estimate_spatial(old, matrix(1, 5L, 3L), features = "all"),
               "before mgcvST 0.0.1.9034")
  # an estimation checkpoint of the earlier format
  dir <- tempfile("mgcvst-edf-")
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  dir.create(dir)
  saveRDS(list(format = 2L, kind = "mgcv", signature = list()), file.path(dir, "estimation-manifest.rds"))
  expect_error(mgcvST:::.mgcvst_chunk_store(dir, "mgcv", list(), TRUE),
               "before 0.0.1.9034")
  # the pair kernel version keys the pair directories
  expect_identical(mgcvST:::.mgcvst_contract("exact")$kernel_version, 3L)
})

.deg_inla_fit <- local({
  cached <- list()
  function(mesh_side = 5L, n = 150L, seed = 1701L) {
    skip_if_not_installed("INLA")
    skip_if_not_installed("geometry")
    key <- paste(mesh_side, n, seed)
    if (!is.null(cached[[key]])) return(cached[[key]])
    set.seed(seed)
    vertices <- as.matrix(expand.grid(x = seq(0, 1, length.out = mesh_side),
                                      y = seq(0, 1, length.out = mesh_side)))
    mesh <- list(loc = vertices, graph = list(tv = geometry::delaunayn(vertices)))
    data <- data.frame(x = runif(n, 0.02, 0.98), y = runif(n, 0.02, 0.98),
                       z = seq(-1, 1, length.out = n), exposure = runif(n, 0.8, 1.3))
    data$offset0 <- log(data$exposure)
    basis <- spde_basis(mesh, as.matrix(data[c("x", "y")]), kappa = if (mesh_side > 10) 3 else 1.2,
                        project_intercept = TRUE)
    eta0 <- 1 + 0.25 * data$z + data$offset0
    set.seed(seed + 1L)
    rnb <- function(f) rnbinom(n, mu = exp(eta0 + f), size = 6)
    Y <- rbind(strong = rnb(.deg_field(basis, 0.30, 1)), weak = rnb(.deg_field(basis, 0.03, 3)),
               medium = rnb(.deg_field(basis, 0.10, 2)),
               none1 = rnb(0), none2 = rnb(0))
    model <- inlaST.set(response ~ z + offset(offset0), data, basis, family = mgcv::nb())
    fit <- suppressMessages(inlaST.estimate(Y, model, BPPARAM = BiocParallel::SerialParam(),
                                            spatial = "all"))
    cached[[key]] <<- fit
    fit
  }
})

test_that("INLA: the field edf is the trace of the reduced curvature, and flags noise genes", {
  skip_on_cran()
  fit <- .deg_inla_fit()
  d <- fit$diagnostics
  expect_identical(fit$format, 3L)
  # With the full-rank observation basis the trace of the reduced curvature is
  # the field block of the hat matrix: the same number as edf_spatial (exact
  # for a small field).
  prepared <- mgcvST:::.inlast_sparse_prepare(fit)
  basis <- mgcvST:::.inlast_sparse_observation_basis(prepared)
  units <- mgcvST:::.inlast_sparse_units(prepared, seq_along(fit$feature_id), threads = 1L)
  reduced <- mgcvST:::.inlast_sparse_materialize_reduced(prepared, units, basis, threads = 1L)
  trace <- vapply(reduced, function(z) sum(diag(z$M)), numeric(1L))
  expect_equal(d$edf_spatial, trace, tolerance = 1e-8)
  expect_identical(d$spatial_degenerate, d$edf_spatial < mgcvST:::.mgcvst_edf_min)
  expect_true(all(d$spatial_degenerate[4:5]))
  expect_false(any(d$spatial_degenerate[1:3]))
  expect_gt(d$edf_spatial[2L], mgcvST:::.mgcvst_edf_min)   # the weak field keeps its degrees of freedom
  for (moments in c("exact", "pcalearning")) {
    z <- inlaST.test(fit, moments = moments, rank = 2L, k = 3L, threads = 2L)
    r <- z$results
    deg <- r$i %in% 4:5 | r$j %in% 4:5
    expect_true(all(r$status[deg] == 4L), info = moments)
    expect_true(all(r$log_p_two_sided[deg] == 0 & r$log_p_positive[deg] == 0 &
                      r$log_p_negative[deg] == 0), info = moments)
    expect_true(all(r$status[!deg] == 0L), info = moments)
    expect_identical(z$degenerate, c("none1", "none2"))
  }
})

test_that("INLA: the trace estimate of a large field agrees with the exact trace", {
  skip_on_cran()
  fit <- .deg_inla_fit(mesh_side = 15L, n = 500L, seed = 99L)
  expect_gt(ncol(fit$score_sparse$Q), 200L)
  prepared <- mgcvST:::.inlast_sparse_prepare(fit)
  basis <- mgcvST:::.inlast_sparse_observation_basis(prepared)
  units <- mgcvST:::.inlast_sparse_units(prepared, seq_along(fit$feature_id), threads = 1L)
  reduced <- mgcvST:::.inlast_sparse_materialize_reduced(prepared, units, basis, threads = 1L)
  trace <- vapply(reduced, function(z) sum(diag(z$M)), numeric(1L))
  edf <- fit$diagnostics$edf_spatial
  # Hutchinson with 128 fixed probes: standard error at most sqrt(2 edf / 128).
  expect_lt(max(abs(edf - trace) / sqrt(2 * pmax(trace, 1e-3) / 128)), 5)
  expect_equal(edf[trace > 10], trace[trace > 10], tolerance = 0.1)
  expect_identical(fit$diagnostics$spatial_degenerate[trace < 0.5], rep(TRUE, sum(trace < 0.5)))
})
