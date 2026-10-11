# Features that the Stage 1 test of the fit did not select. spatial = "all" fits
# a spatial model for every feature; every pair that contains a feature with a
# Stage 1 q-value above q.value (or missing) then has p = 1 and status 4 and
# stays in the adjustment family. A fit with spatial = "discoveries" or with
# feature IDs has no such feature. A pair the kernel cannot evaluate also has
# p = 1 and stays in the family.

.deg_field <- function(basis, v, seed) {
  set.seed(seed)
  f <- as.numeric(basis$B %*% rnorm(ncol(basis$B)))
  f <- f - mean(f)
  f * sqrt(v / stats::var(f))
}

# Eleven genes on a 6 x 6 mesh with 1500 spots: strong, medium and weak
# (v = 0.03) fields, and eight pure-noise genes. At q = 0.05 (BY) Stage 1
# selects the three genes with a field and none of the noise genes.
.deg_ids <- c("strong", "medium", "weak", "none4", "none7")
.deg_pairs <- function() t(combn(.deg_ids, 2L))
.deg_data <- local({
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
    cached <<- list(Y = Y, model = model)
    cached
  }
})
.deg_mgcv_fit <- local({
  cache <- list()
  function(spatial = "all") {
    key <- paste(spatial, collapse = ",")
    if (!is.null(cache[[key]])) return(cache[[key]])
    d <- .deg_data()
    cache[[key]] <<- suppressWarnings(mgcvST.estimate(
      d$Y, d$model, BPPARAM = BiocParallel::SerialParam(), spatial = spatial))
    cache[[key]]
  }
})

test_that("the Stage 1 test of the fit decides which features are set aside under spatial = 'all'", {
  skip_on_cran()
  all_fit <- .deg_mgcv_fit("all")
  d <- all_fit$diagnostics
  expect_identical(all_fit$format, 2L)
  expect_identical(d$spatial_route, rep("all", 11L))
  expect_false("edf_spatial" %in% names(d))
  selected <- !is.na(d$marginal_q_value) & d$marginal_q_value <= all_fit$stage1$q.value
  expect_identical(all_fit$feature_id[selected], c("strong", "medium", "weak"))
  unselected <- mgcvST:::.mgcvst_stage1_unselected(all_fit)
  expect_identical(unselected, !selected)
  # A missing Stage 1 q-value is not a selection.
  missing_q <- all_fit
  missing_q$diagnostics$marginal_q_value[1L] <- NA_real_
  expect_identical(which(mgcvST:::.mgcvst_stage1_unselected(missing_q)), c(1L, 4:11))
  # The threshold is that of the fit.
  strict <- all_fit
  strict$stage1$q.value <- 0.035
  expect_identical(which(mgcvST:::.mgcvst_stage1_unselected(strict)), c(2L, 4:11))
  # A feature without a spatial model is never set aside by this rule.
  unfitted <- all_fit
  unfitted$diagnostics$spatial_fitted[4L] <- FALSE
  expect_false(mgcvST:::.mgcvst_stage1_unselected(unfitted)[4L])

  # Feature IDs are the user's selection, whatever their Stage 1 q-values.
  ids_fit <- .deg_mgcv_fit(.deg_ids)
  expect_identical(ids_fit$diagnostics$spatial_route[ids_fit$diagnostics$spatial_fitted],
                   rep("user", 5L))
  expect_true(all(is.na(ids_fit$diagnostics$spatial_route[!ids_fit$diagnostics$spatial_fitted])))
  expect_false(any(mgcvST:::.mgcvst_stage1_unselected(ids_fit)))
  logical_fit <- .deg_mgcv_fit(c(rep(TRUE, 4L), rep(FALSE, 7L)))
  expect_identical(sum(logical_fit$diagnostics$spatial_route %in% "user"), 4L)
  expect_false(any(mgcvST:::.mgcvst_stage1_unselected(logical_fit)))

  # The default selects the discoveries and sets nothing aside.
  disc <- .deg_mgcv_fit("discoveries")
  expect_identical(disc$feature_id[disc$diagnostics$spatial_fitted],
                   c("strong", "medium", "weak"))
  expect_identical(disc$diagnostics$spatial_route[1:3], rep("discoveries", 3L))
  expect_false(any(mgcvST:::.mgcvst_stage1_unselected(disc)))
})

test_that("estimate_spatial records the route of the features it adds", {
  skip_on_cran()
  d <- .deg_data()
  disc <- .deg_mgcv_fit("discoveries")
  sp <- BiocParallel::SerialParam()
  added <- suppressWarnings(mgcvST.estimate_spatial(disc, d$Y, features = "all", BPPARAM = sp))
  expect_identical(added$diagnostics$spatial_route, c(rep("discoveries", 3L), rep("all", 8L)))
  expect_identical(which(mgcvST:::.mgcvst_stage1_unselected(added)), 4:11)
  expect_identical(added$format, 2L)
  # Feature IDs added later are the user's selection.
  by_id <- suppressWarnings(mgcvST.estimate_spatial(disc, d$Y, features = c("none4", "none7"),
                                                    BPPARAM = sp))
  expect_identical(by_id$diagnostics$spatial_route[c(1:3, 7L, 10L)],
                   c(rep("discoveries", 3L), "user", "user"))
  expect_false(any(mgcvST:::.mgcvst_stage1_unselected(by_id)))
  # A looser Stage 1 rule given to the add-later call is the user's selection.
  loose <- suppressWarnings(mgcvST.estimate_spatial(
    disc, d$Y, features = "discoveries", adjust = "none", q.value = 0.9, BPPARAM = sp))
  newly <- setdiff(which(loose$diagnostics$spatial_fitted), 1:3)
  expect_gt(length(newly), 0L)
  expect_true(all(loose$diagnostics$marginal_q_value[newly] > disc$stage1$q.value))
  expect_identical(loose$diagnostics$spatial_route[newly], rep("discoveries", length(newly)))
  expect_false(any(mgcvST:::.mgcvst_stage1_unselected(loose)))
  # The original fit is not changed.
  expect_identical(disc$diagnostics$spatial_route[4:11], rep(NA_character_, 8L))
})

test_that("a pair with a feature that Stage 1 did not select has p = 1 and status 4 in both routes", {
  skip_on_cran()
  fit <- .deg_mgcv_fit("all")
  P <- .deg_pairs()
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
    expect_true(all(is.finite(r$log_p_two_sided[!deg]) & r$log_p_two_sided[!deg] < 0),
                info = route)
    expect_setequal(z$degenerate, paste0("none", 1:8))
    # The pairs stay in the family: all ten are adjusted, and BY equals
    # stats::p.adjust on the written p-values including the ones.
    expect_equal(z$adjustment$n_adjusted, 10, info = route)
    expect_equal(exp(r$log_q), stats::p.adjust(exp(r$log_p_two_sided), "BY"),
                 tolerance = 1e-12, info = route)
    expect_equal(z$discoveries$pairs_tested, 10)
    expect_equal(z$discoveries$pairs_with_p_value, 10)
    expect_identical(z$contract$kernel_version, 3L)
  }
  # The three pairs of the other genes do not depend on the set-aside genes:
  # the same p-values from a test restricted to them.
  only <- mgcvST.test(fit, pairs = rbind(c("strong", "medium"), c("strong", "weak"),
                                         c("medium", "weak")), moments = "exact")
  full <- routes$exact$results
  expect_equal(only$results$log_p_two_sided,
               full$log_p_two_sided[full$status == 0L], tolerance = 1e-9)

  # The same genes with the user's own selection are tested normally.
  ids_fit <- .deg_mgcv_fit(.deg_ids)
  for (moments in c("exact", "pcalearning")) {
    z <- mgcvST.test(ids_fit, pairs = P, moments = moments, rank = 2L, k = 3L)
    expect_true(all(z$results$status == 0L), info = moments)
    expect_true(all(is.finite(z$results$score)), info = moments)
    expect_true(all(z$results$log_p_two_sided < 0), info = moments)
    expect_identical(z$degenerate, character(0), info = moments)
  }

  # The default selection sets nothing aside.
  disc <- .deg_mgcv_fit("discoveries")
  for (moments in c("exact", "pcalearning")) {
    z <- mgcvST.test(disc, moments = moments, rank = 2L, k = 3L)
    expect_identical(nrow(z$results), 3L)
    expect_true(all(z$results$status == 0L), info = moments)
    expect_identical(z$degenerate, character(0), info = moments)
  }
})

test_that("a test whose genes are all set aside gives p = 1 for every pair", {
  skip_on_cran()
  fit <- .deg_mgcv_fit("all")
  z <- mgcvST.test(fit, pairs = rbind(c("none4", "none7")), moments = "exact")
  expect_identical(z$results$status, 4L)
  expect_identical(z$results$log_p_two_sided, 0)
  expect_match(z$contract$basis_sha, "^[0-9a-f]{64}$")
  expect_true(all(c("none4", "none7") %in% z$degenerate))
  # One gene set aside and one usable gene.
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
  fit <- .deg_mgcv_fit("all")
  pairs_kernel <- mgcvST:::mgcvst_pca_spa_pairs_cpp
  testthat::local_mocked_bindings(
    mgcvst_pca_spa_pairs_cpp = function(...) {
      out <- pairs_kernel(...)
      out[1L, c("logp_two_sided", "logp_positive", "logp_negative")] <- NA_real_
      out[1L, "status"] <- 2
      out
    }, .package = "mgcvST")
  z <- mgcvST.test(fit, pairs = .deg_pairs(), moments = "pcalearning", rank = 2L, k = 3L,
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

test_that("fits and estimation checkpoints of format 2 are accepted", {
  skip_on_cran()
  d <- .deg_data()
  expect_identical(mgcvST:::.mgcvst_fit_format, 2L)
  # A fit of 0.0.1.9032 or 0.0.1.9033 has no route: the selection is the user's.
  old <- .deg_mgcv_fit("discoveries")
  old$diagnostics$spatial_route <- NULL
  expect_identical(old$format, 2L)
  expect_false(any(mgcvST:::.mgcvst_stage1_unselected(old)))
  for (moments in c("exact", "pcalearning")) {
    z <- mgcvST.test(old, moments = moments, rank = 2L, k = 3L)
    expect_true(all(z$results$status == 0L), info = moments)
    expect_identical(z$degenerate, character(0), info = moments)
  }
  # ... even when every feature has a spatial model.
  old_all <- .deg_mgcv_fit("all")
  old_all$diagnostics$spatial_route <- NULL
  z <- mgcvST.test(old_all, pairs = .deg_pairs(), moments = "exact")
  expect_true(all(z$results$status == 0L))
  expect_identical(z$degenerate, character(0))
  # It is extended by the add-later function; only the new features get a route.
  ext <- suppressWarnings(mgcvST.estimate_spatial(
    old, d$Y, features = "all", BPPARAM = BiocParallel::SerialParam()))
  expect_identical(ext$diagnostics$spatial_route, c(rep(NA_character_, 3L), rep("all", 8L)))
  expect_identical(which(mgcvST:::.mgcvst_stage1_unselected(ext)), 4:11)

  # An estimation checkpoint of format 2 is resumed.
  dir <- tempfile("mgcvst-format2-")
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  first <- suppressWarnings(mgcvST.estimate(
    d$Y, d$model, BPPARAM = BiocParallel::SerialParam(), spatial = "discoveries",
    chunk_size = 4L, checkpoint_dir = dir))
  manifest <- file.path(dir, "estimation-manifest.rds")
  expect_identical(readRDS(manifest)$format, 2L)
  again <- suppressWarnings(mgcvST.estimate(
    d$Y, d$model, BPPARAM = BiocParallel::SerialParam(), spatial = "discoveries",
    chunk_size = 4L, checkpoint_dir = dir))
  expect_gt(again$timing$resumed_null_chunks, 0L)
  expect_identical(again$diagnostics$marginal_p_value, first$diagnostics$marginal_p_value)
  # The pair kernel version keys the pair directories.
  expect_identical(mgcvST:::.mgcvst_contract("exact")$kernel_version, 3L)
})

.deg_inla_fit <- local({
  cache <- list()
  function(spatial = "all") {
    skip_if_not_installed("INLA")
    skip_if_not_installed("geometry")
    key <- paste(spatial, collapse = ",")
    if (!is.null(cache[[key]])) return(cache[[key]])
    n <- 150L
    set.seed(1701L)
    vertices <- as.matrix(expand.grid(x = seq(0, 1, length.out = 5L),
                                      y = seq(0, 1, length.out = 5L)))
    mesh <- list(loc = vertices, graph = list(tv = geometry::delaunayn(vertices)))
    data <- data.frame(x = runif(n, 0.02, 0.98), y = runif(n, 0.02, 0.98),
                       z = seq(-1, 1, length.out = n), exposure = runif(n, 0.8, 1.3))
    data$offset0 <- log(data$exposure)
    basis <- spde_basis(mesh, as.matrix(data[c("x", "y")]), kappa = 1.2,
                        project_intercept = TRUE)
    eta0 <- 1 + 0.25 * data$z + data$offset0
    set.seed(1702L)
    rnb <- function(f) rnbinom(n, mu = exp(eta0 + f), size = 6)
    Y <- rbind(strong = rnb(.deg_field(basis, 0.30, 1)), weak = rnb(.deg_field(basis, 0.03, 3)),
               medium = rnb(.deg_field(basis, 0.10, 2)),
               none1 = rnb(0), none2 = rnb(0))
    model <- inlaST.set(response ~ z + offset(offset0), data, basis, family = mgcv::nb())
    cache[[key]] <<- list(
      fit = suppressMessages(inlaST.estimate(Y, model, BPPARAM = BiocParallel::SerialParam(),
                                             spatial = spatial)),
      Y = Y)
    cache[[key]]
  }
})

test_that("INLA: the Stage 1 test of the fit decides which features are set aside", {
  skip_on_cran()
  fit <- .deg_inla_fit("all")$fit
  d <- fit$diagnostics
  expect_identical(fit$format, 2L)
  expect_identical(d$spatial_route, rep("all", 5L))
  # Stage 1 selects strong and weak; medium, none1 and none2 are set aside.
  expect_identical(which(mgcvST:::.mgcvst_stage1_unselected(fit)), 3:5)
  for (moments in c("exact", "pcalearning")) {
    z <- inlaST.test(fit, moments = moments, rank = 2L, k = 3L, threads = 2L, adjust = "BY")
    r <- z$results
    deg <- r$i %in% 3:5 | r$j %in% 3:5
    expect_identical(sum(deg), 9L)
    expect_true(all(r$status[deg] == 4L), info = moments)
    expect_true(all(r$log_p_two_sided[deg] == 0 & r$log_p_positive[deg] == 0 &
                      r$log_p_negative[deg] == 0), info = moments)
    expect_true(all(is.na(r$score[deg])), info = moments)
    expect_identical(r$status[!deg], 0L, info = moments)
    expect_true(is.finite(r$log_p_two_sided[!deg]) && r$log_p_two_sided[!deg] < 0)
    expect_setequal(z$degenerate, c("medium", "none1", "none2"))
    expect_equal(z$adjustment$n_adjusted, 10, info = moments)
    expect_equal(exp(r$log_q), stats::p.adjust(exp(r$log_p_two_sided), "BY"),
                 tolerance = 1e-12, info = moments)
  }
  # The user's own selection of the same features is tested normally.
  ids <- .deg_inla_fit(c("strong", "weak", "medium", "none1", "none2"))$fit
  expect_identical(ids$diagnostics$spatial_route, rep("user", 5L))
  z <- inlaST.test(ids, moments = "exact", threads = 2L)
  expect_true(all(z$results$status == 0L))
  expect_identical(z$degenerate, character(0))
  # The default selection sets nothing aside, and the add-later function
  # records the route of the features it adds.
  disc <- .deg_inla_fit("discoveries")
  expect_identical(disc$fit$diagnostics$spatial_route,
                   c("discoveries", "discoveries", rep(NA_character_, 3L)))
  expect_false(any(mgcvST:::.mgcvst_stage1_unselected(disc$fit)))
  added <- suppressMessages(inlaST.estimate_spatial(
    disc$fit, disc$Y, features = "all", BPPARAM = BiocParallel::SerialParam()))
  expect_identical(added$diagnostics$spatial_route,
                   c("discoveries", "discoveries", rep("all", 3L)))
  expect_identical(which(mgcvST:::.mgcvst_stage1_unselected(added)), 3:5)
  # A fit of format 2 without a route is accepted and sets nothing aside.
  old <- fit
  old$diagnostics$spatial_route <- NULL
  z <- inlaST.test(old, moments = "exact", threads = 2L)
  expect_true(all(z$results$status == 0L))
  expect_identical(z$degenerate, character(0))
})
