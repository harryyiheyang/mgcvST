# Two-step estimation: Stage 1 null fits for every feature, spatial fits for the
# selected features only, in both branches.

# ---- helpers ----------------------------------------------------------------

# Any numeric vector or array with at least `n` elements anywhere in `x`.
.ts_has_long_numeric <- function(x, n) {
  found <- FALSE
  walk <- function(z) {
    if (is.numeric(z) && length(z) >= n) {
      found <<- TRUE
    } else if (is.list(z)) {
      for (element in z) walk(element)
    }
  }
  walk(x)
  found
}

.ts_chunks <- function(dir, step) {
  files <- list.files(dir, paste0("^", step, "-.*[.]rds$"), full.names = TRUE)
  lapply(files, function(file) readRDS(file)$result)
}

# INLA fits are not reproducible from run to run, so the INLA tests memoize the
# per-feature fits: every call of the engine with the same response, model and
# controls returns the first result.
.ts_memo <- new.env(parent = emptyenv())

.ts_memoize_fits <- function(env = parent.frame()) {
  real <- mgcvST:::.inlast_fit_feature
  testthat::local_mocked_bindings(
    .inlast_fit_feature = function(spec, y, offset = NULL, control = list(),
                                   diagnostics = FALSE) {
      key <- digest::digest(list(y, spec$family,
        vapply(spec$random, function(z) z$name, ""), offset, control, diagnostics))
      if (is.null(.ts_memo[[key]])) {
        .ts_memo[[key]] <- real(spec, y, offset = offset, control = control,
                                diagnostics = diagnostics)
      }
      .ts_memo[[key]]
    },
    .package = "mgcvST", .env = env
  )
}

.ts_inla <- local({
  cached <- NULL
  function() {
    skip_if_not_installed("INLA")
    skip_if_not_installed("geometry")
    if (!is.null(cached)) return(cached)
    withr::local_seed(2026L)
    n <- 60L
    vertices <- as.matrix(expand.grid(x = seq(0, 1, length.out = 4L),
                                      y = seq(0, 1, length.out = 4L)))
    mesh <- list(loc = vertices, graph = list(tv = geometry::delaunayn(vertices)))
    data <- data.frame(x = runif(n, 0.02, 0.98), y = runif(n, 0.02, 0.98),
                       z = seq(-1, 1, length.out = n), exposure = runif(n, 0.8, 1.3))
    data$offset0 <- log(data$exposure)
    basis <- spde_basis(mesh, as.matrix(data[c("x", "y")]), kappa = 1.2,
                        project_intercept = TRUE)
    G <- 6L
    Y <- t(vapply(seq_len(G), function(g) {
      eta <- 1 + 0.25 * data$z + data$offset0 + 0.5 * sin(2 * pi * (data$x + g / G))
      rnbinom(n, mu = exp(eta), size = 3 + g)
    }, numeric(n)))
    dimnames(Y) <- list(paste0("g", seq_len(G)), NULL)
    model <- inlaST.set(response ~ z + offset(offset0), data, basis, family = mgcv::nb())
    cached <<- list(Y = Y, model = model, n = n,
                    m = ncol(model$inla_spec$random[[1L]]$A))
    cached
  }
})

# Identical estimates: every field the tests read, without timings.
.ts_same_estimates <- function(a, b, fields) {
  for (name in fields) expect_identical(a[[name]], b[[name]], info = name)
  keep <- !grepl("_seconds$", names(a$diagnostics))
  expect_identical(a$diagnostics[, keep], b$diagnostics[, keep])
}

# ---- Stage 1 selection ------------------------------------------------------

test_that("spatial selection resolves discoveries, all, none, IDs, indices and logicals", {
  ids <- paste0("f", 1:5)
  q <- c(0.001, 0.2, NA, 0.04, 0.5)
  sel <- function(x, qv = 0.05) mgcvST:::.mgcvst_select_spatial(x, ids, q, qv)
  expect_identical(sel("discoveries"), c(1L, 4L))
  expect_identical(sel("discoveries", 0.3), c(1L, 2L, 4L))
  expect_identical(sel("all"), 1:5)
  expect_identical(sel("none"), integer())
  expect_identical(sel(c("f4", "f2")), c(2L, 4L))
  expect_identical(sel(c(5, 1)), c(1L, 5L))
  expect_identical(sel(c(TRUE, FALSE, FALSE, TRUE, FALSE)), c(1L, 4L))
  expect_error(sel("f9"), "unknown feature IDs")
  expect_error(sel(7), "valid one-based")
  expect_error(sel(1.5), "valid one-based")
  expect_error(sel(c(TRUE, FALSE)), "one non-missing value per feature")
  expect_error(sel(NULL), "must not be NULL")
  expect_error(sel(list(1)), "spatial must be")
})

test_that("Stage 1 q-values are the log-space adjustment of the null p-values", {
  p <- c(1e-6, 0.3, NA, 0.02, 1e-3, 0.9)
  for (adjust in c("BY", "BH", "none")) {
    expect_equal(mgcvST:::.mgcvst_stage1_q(p, adjust), p.adjust(p, adjust),
                 tolerance = 1e-12)
  }
  expect_true(is.na(mgcvST:::.mgcvst_stage1_q(p, "BY")[3L]))
  expect_equal(mgcvST:::.mgcvst_stage1_q(c(0, 0.5), "BH"), c(0, 0.5))
  expect_error(mgcvST:::.mgcvst_check_q_value(0), "q.value")
  expect_error(mgcvST:::.mgcvst_check_q_value(1.5), "q.value")
  expect_error(mgcvST:::.mgcvst_check_q_value(c(0.1, 0.2)), "q.value")
  expect_identical(mgcvST:::.mgcvst_check_q_value(1), 1)
})

# ---- chunk checkpoints ------------------------------------------------------

test_that("chunk checkpoints are resumed, and refused when they do not match", {
  dir <- tempfile("mgcvst-chunks-")
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  store <- mgcvST:::.mgcvst_chunk_store(dir, "inla", "signature-a", TRUE)
  digest <- c("a", "b", "c", "d")
  calls <- 0L
  work <- function(payload, scale) {
    calls <<- calls + length(payload$index)
    out <- as.list(payload$index * scale)
    mgcvST:::.mgcvst_chunk_save(payload$chunk_file, payload$chunk_key, out)
    out
  }
  payloads <- list(list(index = 1:2), list(index = 3:4))
  sp <- BiocParallel::SerialParam()
  first <- mgcvST:::.mgcvst_run_chunks(payloads, "null", store, digest, sp, work, scale = 10)
  expect_identical(first$resumed, 0L)
  expect_identical(unlist(first$results), c(10, 20, 30, 40))
  expect_identical(calls, 4L)
  again <- mgcvST:::.mgcvst_run_chunks(payloads, "null", store, digest, sp, work, scale = 10)
  expect_identical(again$resumed, 2L)
  expect_identical(calls, 4L)
  expect_identical(again$results, first$results)
  # A chunk is keyed by its features and their response digests.
  other <- mgcvST:::.mgcvst_run_chunks(payloads, "null", store, c("a", "b", "x", "d"),
                                       sp, work, scale = 10)
  expect_identical(other$resumed, 1L)
  # A step has its own chunks.
  spatial <- mgcvST:::.mgcvst_run_chunks(payloads, "spatial", store, digest, sp, work,
                                         scale = 10)
  expect_identical(spatial$resumed, 0L)

  # A damaged chunk is reported, not recomputed silently.
  files <- list.files(dir, "^null-", full.names = TRUE)
  writeBin(as.raw(1:20), files[1L])
  expect_error(mgcvST:::.mgcvst_run_chunks(payloads, "null", store, digest, sp, work,
                                           scale = 10), "damaged")

  # The manifest ties the directory to one estimator, format and signature.
  expect_error(mgcvST:::.mgcvst_chunk_store(dir, "mgcv", "signature-a", TRUE),
               "another estimator or by a version before 0.0.1.9032")
  expect_error(mgcvST:::.mgcvst_chunk_store(dir, "inla", "signature-b", TRUE),
               "different model, offset or controls")
  expect_error(mgcvST:::.mgcvst_chunk_store(dir, "inla", "signature-a", FALSE),
               "already exists")
  manifest <- file.path(dir, "estimation-manifest.rds")
  record <- readRDS(manifest)
  earlier <- record
  earlier$format <- 1L
  saveRDS(earlier, manifest)
  expect_error(mgcvST:::.mgcvst_chunk_store(dir, "inla", "signature-a", TRUE),
               "before 0.0.1.9032")
  saveRDS(record, manifest)
  expect_identical(mgcvST:::.mgcvst_chunk_store(dir, "inla", "signature-a", TRUE)$kind,
                   "inla")
  loose <- tempfile("mgcvst-loose-")
  dir.create(loose)
  on.exit(unlink(loose, recursive = TRUE), add = TRUE)
  writeLines("x", file.path(loose, "other.txt"))
  expect_error(mgcvST:::.mgcvst_chunk_store(loose, "inla", "signature-a", TRUE),
               "no manifest")
  expect_null(mgcvST:::.mgcvst_chunk_store(NULL, "inla", "signature-a", TRUE))
  expect_error(mgcvST:::.mgcvst_chunk_store(c("a", "b"), "inla", "s", TRUE),
               "checkpoint_dir")
})

test_that("a fit estimated before the two-step estimators is refused only where it lacks data", {
  inla_old <- list(estimator = "INLA")
  expect_error(mgcvST:::.mgcvst_check_fit_format(inla_old), "before mgcvST 0.0.1.9032")
  inla_old$format <- 1L
  expect_error(mgcvST:::.mgcvst_check_fit_format(inla_old), "re-run inlaST.estimate")
  inla_new <- list(estimator = "INLA", format = 2L, mu_bar = c(1, 2))
  expect_identical(mgcvST:::.mgcvst_check_fit_format(inla_new), inla_new)
  expect_error(mgcvST:::.mgcvst_check_fit_format(list(estimator = "INLA", format = 2L)),
               "re-run inlaST.estimate")
  # An mgcv fit of the earlier format holds everything the mgcv tests read.
  expect_silent(mgcvST:::.mgcvst_check_fit_format(list(feature_id = "a")))
  expect_error(mgcvST:::.mgcvst_check_fit_format(list(format = 99L)), "newer version")
})

# ---- mgcv branch ------------------------------------------------------------

test_that("mgcvST.estimate fits the spatial model of the Stage 1 discoveries only", {
  f <- st_fixture()
  sp <- BiocParallel::SerialParam()
  all <- suppressWarnings(mgcvST.estimate(f$Y, f$model, BPPARAM = sp, spatial = "all"))
  q <- all$diagnostics$marginal_q_value
  expect_equal(q, p.adjust(all$diagnostics$marginal_p_value, "BY"), tolerance = 1e-12)
  expect_true(all(all$diagnostics$spatial_selected & all$diagnostics$spatial_fitted))
  expect_true(all(mgcvST:::.mgcvst_feature_available(all)))
  # A threshold between the second and the third q-value selects two features.
  ordered <- sort(q)
  cut <- sqrt(ordered[2L] * ordered[3L])
  fit <- suppressWarnings(mgcvST.estimate(f$Y, f$model, BPPARAM = sp, q.value = cut))
  chosen <- which(q <= cut)
  expect_length(chosen, 2L)
  rest <- setdiff(seq_along(q), chosen)
  expect_identical(fit$diagnostics$spatial_selected, q <= cut)
  expect_identical(fit$diagnostics$spatial_fitted, q <= cut)
  expect_identical(fit$diagnostics$marginal_p_value, all$diagnostics$marginal_p_value)
  expect_identical(fit$stage1, list(adjust = "BY", q.value = cut))
  # Unselected features carry no working model and are unavailable.
  expect_true(all(is.na(fit$working_error[, rest])))
  expect_true(all(is.na(fit$working_variance[, rest])))
  expect_true(all(is.na(fit$dispersion[rest])) && all(is.na(fit$lambda[rest])))
  expect_identical(unname(mgcvST:::.mgcvst_feature_available(fit)), q <= cut)
  expect_false(fit$diagnostics$converged[rest])
  # Selected features are the fits of the one-step run.
  expect_identical(fit$working_error[, chosen], all$working_error[, chosen])
  expect_identical(fit$working_variance[, chosen], all$working_variance[, chosen])
  expect_identical(fit$dispersion[chosen], all$dispersion[chosen])
  expect_identical(fit$smoothing_parameters[chosen, ], all$smoothing_parameters[chosen, ])
  expect_identical(fit$nuisance_covariance[chosen], all$nuisance_covariance[chosen])
  expect_output(print(fit), "spatial models: 2 of 3 features")
  # pairs = NULL tests the pairs among the fitted features; the others are not failures.
  tested <- mgcvST.test(fit)
  expect_identical(nrow(tested$results), 1L)
  expect_identical(c(tested$results$i, tested$results$j), chosen)
  expect_identical(nrow(tested$failed), 0L)
  expect_identical(tested$results$log_p_two_sided,
                   mgcvST.test(all, pairs = rbind(chosen))$results$log_p_two_sided)
  # An explicit pair with an unselected feature is reported, not tested.
  bad <- mgcvST.test(fit, pairs = rbind(c(chosen[1L], rest)))
  expect_identical(bad$results$status, 3L)
  expect_match(bad$failed$error, "not selected in step 2")
  expect_identical(bad$failed$feature_id, fit$feature_id[rest])
  expect_error(mgcvST.wgcna(fit, indices = fit$feature_id), "no spatial fit")

  none <- suppressWarnings(mgcvST.estimate(f$Y, f$model, BPPARAM = sp, spatial = "none"))
  expect_false(any(none$diagnostics$spatial_fitted))
  expect_true(is.null(none$geometry))
  expect_error(mgcvST.test(none), "no spatial model")
  one_model <- suppressWarnings(mgcvST.estimate(f$Y, f$model, BPPARAM = sp,
                                                spatial = "response2"))
  expect_error(mgcvST.test(one_model), "At least two available features")
  by_id <- suppressWarnings(mgcvST.estimate(f$Y, f$model, BPPARAM = sp,
                                            spatial = c("response3", "response")))
  expect_identical(by_id$diagnostics$spatial_fitted, c(TRUE, FALSE, TRUE))
  by_flag <- suppressWarnings(mgcvST.estimate(f$Y, f$model, BPPARAM = sp,
                                              spatial = c(FALSE, TRUE, FALSE)))
  expect_identical(by_flag$diagnostics$spatial_fitted, c(FALSE, TRUE, FALSE))
  expect_error(mgcvST.estimate(f$Y, f$model, spatial = "response9"), "unknown feature IDs")
  expect_error(mgcvST.estimate(f$Y, f$model, q.value = 0), "q.value")
  expect_error(mgcvST.estimate(f$Y, f$model, adjust = "holm"), "should be one of")
})

test_that("mgcvST.estimate_spatial adds spatial models to a step 1 fit", {
  f <- st_fixture()
  sp <- BiocParallel::SerialParam()
  all <- suppressWarnings(mgcvST.estimate(f$Y, f$model, BPPARAM = sp, spatial = "all"))
  none <- suppressWarnings(mgcvST.estimate(f$Y, f$model, BPPARAM = sp, spatial = "none"))
  one <- suppressWarnings(mgcvST.estimate_spatial(none, f$Y, "response2", BPPARAM = sp))
  # The supplied fit is not changed.
  expect_false(any(none$diagnostics$spatial_fitted))
  expect_true(is.null(none$geometry))
  expect_identical(one$diagnostics$spatial_fitted, c(FALSE, TRUE, FALSE))
  expect_identical(one$diagnostics$spatial_selected, c(FALSE, TRUE, FALSE))
  expect_false(is.null(one$geometry))
  # A feature that has a spatial model is skipped; the rest are added.
  full <- suppressWarnings(mgcvST.estimate_spatial(one, f$Y, "all", BPPARAM = sp))
  fields <- c("working_error", "working_variance", "dispersion", "lambda",
              "component_lambda", "smoothing_parameters", "nuisance_covariance",
              "family_parameters", "offset", "row_id", "linear_design",
              "score_components", "y_digest")
  .ts_same_estimates(full, all, fields)
  # The shared geometry carries the smoothing parameter of the feature that
  # established it; everything else in it is the same.
  geometry <- function(x) {
    x$geometry$sp <- NULL
    x$geometry
  }
  expect_identical(geometry(full), geometry(all))
  expect_identical(mgcvST.test(full)$results, mgcvST.test(all)$results)
  # Discoveries are selected from the stored Stage 1 p-values.
  q <- all$diagnostics$marginal_q_value
  cut <- sqrt(sort(q)[2L] * sort(q)[3L])
  some <- suppressWarnings(mgcvST.estimate_spatial(none, f$Y, q.value = cut, BPPARAM = sp))
  expect_identical(some$diagnostics$spatial_fitted, q <= cut)
  same_again <- suppressWarnings(mgcvST.estimate_spatial(full, f$Y, "all", BPPARAM = sp))
  expect_identical(same_again$working_error, full$working_error)

  expect_error(mgcvST.estimate_spatial(none, f$Y[, rev(seq_len(ncol(f$Y)))], "all"),
               "Y differs from the responses of step 1 for response")
  expect_error(mgcvST.estimate_spatial(none, f$Y[, -1L], "all"), "Y must be the feature")
  expect_error(mgcvST.estimate_spatial(none, f$Y, "nothing"), "unknown feature IDs")
  expect_error(mgcvST.estimate_spatial(list(), f$Y), "must be returned by mgcvST.estimate")
  stale <- none
  stale$estimation_context <- NULL
  expect_error(mgcvST.estimate_spatial(stale, f$Y, "all"), "cannot be extended")
})

test_that("mgcv null chunks hold no observation-length vector and resume exactly", {
  f <- st_fixture()
  n <- ncol(f$Y)
  sp <- BiocParallel::SerialParam()
  dir <- tempfile("mgcvst-two-step-")
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  first <- suppressWarnings(mgcvST.estimate(
    f$Y, f$model, BPPARAM = sp, spatial = "all", chunk_size = 1L,
    checkpoint_dir = dir))
  expect_identical(first$timing$resumed_null_chunks, 0L)
  nulls <- .ts_chunks(dir, "null")
  spatial <- .ts_chunks(dir, "spatial")
  expect_length(nulls, 3L)
  expect_length(spatial, 3L)
  # What a worker returns of a null fit is a handful of scalars.
  for (chunk in nulls) {
    expect_false(.ts_has_long_numeric(chunk, n))
    expect_lt(as.numeric(object.size(chunk)), 5000)
  }
  expect_true(all(vapply(spatial, .ts_has_long_numeric, logical(1L), n = n)))

  again <- suppressWarnings(mgcvST.estimate(
    f$Y, f$model, BPPARAM = sp, spatial = "all", chunk_size = 1L,
    checkpoint_dir = dir))
  expect_identical(again$timing$resumed_null_chunks, 3L)
  expect_identical(again$timing$resumed_spatial_chunks, 3L)
  fields <- c("working_error", "working_variance", "dispersion", "lambda",
              "smoothing_parameters", "nuisance_covariance", "family_parameters")
  .ts_same_estimates(again, first, fields)

  # A lost chunk is recomputed, and only that one.
  unlink(list.files(dir, "^spatial-", full.names = TRUE)[2L])
  unlink(list.files(dir, "^null-", full.names = TRUE)[1L])
  partial <- suppressWarnings(mgcvST.estimate(
    f$Y, f$model, BPPARAM = sp, spatial = "all", chunk_size = 1L,
    checkpoint_dir = dir))
  expect_identical(partial$timing$resumed_null_chunks, 2L)
  expect_identical(partial$timing$resumed_spatial_chunks, 2L)
  .ts_same_estimates(partial, first, fields)

  # Step 2 of a finished step 1 reuses the null chunks of the same directory.
  none_dir <- tempfile("mgcvst-two-step-none-")
  on.exit(unlink(none_dir, recursive = TRUE), add = TRUE)
  none <- suppressWarnings(mgcvST.estimate(
    f$Y, f$model, BPPARAM = sp, spatial = "none", chunk_size = 1L,
    checkpoint_dir = none_dir))
  expect_length(list.files(none_dir, "^spatial-"), 0L)
  added <- suppressWarnings(mgcvST.estimate_spatial(
    none, f$Y, "all", BPPARAM = sp, chunk_size = 1L, checkpoint_dir = none_dir))
  expect_length(list.files(none_dir, "^spatial-"), 3L)
  .ts_same_estimates(added, first, fields)

  # Another model, offset or control is refused rather than mixed in.
  expect_error(suppressWarnings(mgcvST.estimate(
    f$Y, f$model, BPPARAM = sp, spatial = "all", chunk_size = 1L,
    offset = rep(0.1, n), checkpoint_dir = dir)), "different model, offset or controls")
  expect_error(suppressWarnings(mgcvST.estimate(
    f$Y, f$model, BPPARAM = sp, spatial = "all", chunk_size = 1L,
    checkpoint_dir = dir, resume = FALSE)), "already exists")
  inla_dir <- tempfile("mgcvst-two-step-inla-")
  on.exit(unlink(inla_dir, recursive = TRUE), add = TRUE)
  mgcvST:::.mgcvst_chunk_store(inla_dir, "inla", "x", TRUE)
  expect_error(suppressWarnings(mgcvST.estimate(
    f$Y, f$model, BPPARAM = sp, checkpoint_dir = inla_dir)), "another estimator")
  # Other responses recompute their chunks; none of the old ones is reused.
  shifted <- f$Y
  shifted[2L, ] <- rev(shifted[2L, ])
  changed <- suppressWarnings(mgcvST.estimate(
    shifted, f$model, BPPARAM = sp, spatial = "all", chunk_size = 1L,
    checkpoint_dir = dir))
  expect_identical(changed$timing$resumed_null_chunks, 2L)
  expect_identical(changed$timing$resumed_spatial_chunks, 2L)
})

# ---- INLA branch ------------------------------------------------------------

test_that("inlaST.estimate fits the null model of every feature and the spatial model of the discoveries", {
  skip_on_cran()
  d <- .ts_inla()
  .ts_memoize_fits()
  sp <- BiocParallel::SerialParam()
  all <- inlaST.estimate(d$Y, d$model, BPPARAM = sp, spatial = "all")
  expect_identical(all$format, 2L)
  expect_s3_class(all, "inlaST_fit")
  p <- all$diagnostics$marginal_p_value
  expect_true(all(is.finite(p)))
  expect_equal(all$diagnostics$marginal_q_value, p.adjust(p, "BY"), tolerance = 1e-12)
  expect_true(all(all$diagnostics$spatial_selected & all$diagnostics$spatial_fitted))
  expect_true(all(is.finite(all$mu_bar)))

  q <- all$diagnostics$marginal_q_value
  cut <- sqrt(sort(q)[3L] * sort(q)[4L])
  fit <- inlaST.estimate(d$Y, d$model, BPPARAM = sp, q.value = cut)
  chosen <- which(q <= cut)
  rest <- setdiff(seq_along(q), chosen)
  expect_length(chosen, 3L)
  expect_identical(fit$diagnostics$spatial_selected, q <= cut)
  expect_identical(fit$diagnostics$spatial_fitted, q <= cut)
  expect_identical(fit$diagnostics$marginal_p_value, p)
  # Every feature keeps its null state; only the selected ones have a spatial model.
  expect_true(all(vapply(fit$null_state, is.list, logical(1L))))
  expect_true(all(is.na(fit$score_a[, rest])))
  expect_true(all(is.na(fit$dispersion[rest])) && all(is.na(fit$mu_bar[rest])))
  expect_true(all(is.na(fit$target_coefficients[, rest])))
  expect_identical(unname(mgcvST:::.mgcvst_feature_available(fit)), q <= cut)
  for (name in c("score_a", "target_coefficients", "nuisance_coefficients")) {
    expect_identical(fit[[name]][, chosen], all[[name]][, chosen], info = name)
  }
  expect_identical(fit$mu_bar[chosen], all$mu_bar[chosen])
  expect_null(fit$working_error)
  expect_null(fit$working_variance)
  expect_false(.ts_has_long_numeric(fit$null_state, d$n))
  expect_output(print(fit), "spatial models: 3 of 6 features")

  # The test covers the fitted features; the others are not failures.
  tested <- inlaST.test(fit, rank = 3L, seed = 4L)
  expect_identical(nrow(tested$results), 3L)
  expect_true(all(c(tested$results$i, tested$results$j) %in% chosen))
  expect_identical(nrow(tested$failed), 0L)
  bad <- inlaST.test(fit, pairs = rbind(c(chosen[1L], rest[1L])), rank = 3L, seed = 4L)
  expect_identical(bad$results$status, 3L)
  expect_match(bad$failed$error, "not selected in step 2")
  expect_error(inlaST.wgcna(fit, indices = fit$feature_id), "no spatial fit")

  none <- inlaST.estimate(d$Y, d$model, BPPARAM = sp, spatial = "none")
  expect_false(any(none$diagnostics$spatial_fitted))
  expect_true(all(is.na(none$score_a)))
  expect_identical(none$diagnostics$marginal_p_value, p)
  by_id <- inlaST.estimate(d$Y, d$model, BPPARAM = sp, spatial = c("g6", "g1"))
  expect_identical(which(by_id$diagnostics$spatial_fitted), c(1L, 6L))
  expect_error(inlaST.estimate(d$Y, d$model, spatial = "g99"), "unknown feature IDs")
  expect_error(inlaST.estimate(d$Y, d$model, q.value = 2), "q.value")
  expect_error(inlaST.estimate(d$Y, d$model, adjust = "holm"), "should be one of")
})

test_that("inlaST.estimate_spatial adds spatial models to a step 1 fit", {
  skip_on_cran()
  d <- .ts_inla()
  .ts_memoize_fits()
  sp <- BiocParallel::SerialParam()
  all <- inlaST.estimate(d$Y, d$model, BPPARAM = sp, spatial = "all")
  none <- inlaST.estimate(d$Y, d$model, BPPARAM = sp, spatial = "none")
  two <- inlaST.estimate_spatial(none, d$Y, c("g2", "g5"), BPPARAM = sp)
  expect_false(any(none$diagnostics$spatial_fitted))
  expect_identical(which(two$diagnostics$spatial_fitted), c(2L, 5L))
  expect_false(is.null(two$estimation))
  full <- inlaST.estimate_spatial(two, d$Y, "all", BPPARAM = sp)
  fields <- c("score_a", "target_coefficients", "nuisance_coefficients", "mu_bar",
              "dispersion", "lambda", "component_lambda", "smoothing_parameters",
              "family_parameters", "null_state", "constraint_residual",
              "observation_spatial_mean", "feature_family", "y_digest")
  .ts_same_estimates(full, all, fields)
  expect_identical(inlaST.test(full, rank = 3L, seed = 4L)$results,
                   inlaST.test(all, rank = 3L, seed = 4L)$results)
  q <- all$diagnostics$marginal_q_value
  cut <- sqrt(sort(q)[3L] * sort(q)[4L])
  some <- inlaST.estimate_spatial(none, d$Y, q.value = cut, BPPARAM = sp)
  expect_identical(some$diagnostics$spatial_fitted, q <= cut)

  expect_error(inlaST.estimate_spatial(none, d$Y[, rev(seq_len(ncol(d$Y)))], "g1"),
               "Y differs from the responses of step 1 for g1")
  expect_error(inlaST.estimate_spatial(none, d$Y[, -1L], "g1"), "Y must be the feature")
  expect_error(inlaST.estimate_spatial(none, d$Y, "g99"), "unknown feature IDs")
  expect_error(inlaST.estimate_spatial(list(), d$Y), "must be returned by inlaST.estimate")
})

test_that("INLA workers return compact results and a checkpoint resumes exactly", {
  skip_on_cran()
  d <- .ts_inla()
  .ts_memoize_fits()
  sp <- BiocParallel::SerialParam()
  dir <- tempfile("mgcvst-two-step-inla-")
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  first <- inlaST.estimate(d$Y, d$model, BPPARAM = sp, spatial = "all",
                           chunk_size = 2L, checkpoint_dir = dir)
  expect_identical(first$timing$resumed_null_chunks, 0L)
  nulls <- .ts_chunks(dir, "null")
  spatial <- .ts_chunks(dir, "spatial")
  expect_length(nulls, 3L)
  expect_length(spatial, 3L)
  # Nothing observation-length reaches the manager from either step.
  for (chunk in c(nulls, spatial)) {
    expect_false(.ts_has_long_numeric(chunk, d$n))
    expect_true(all(vapply(chunk, function(z) is.null(z$error), logical(1L))))
  }
  for (z in unlist(nulls, recursive = FALSE)) {
    expect_null(z$null$working_error)
    expect_null(z$null$working_variance)
    expect_null(z$null$eta)
    expect_null(z$null$mu)
  }
  expect_true(all(vapply(unlist(spatial, recursive = FALSE),
                         function(z) all(is.finite(z$score_a)) && is.finite(z$mu_bar),
                         logical(1L))))

  fields <- c("score_a", "target_coefficients", "nuisance_coefficients", "mu_bar",
              "dispersion", "lambda", "smoothing_parameters", "family_parameters",
              "null_state")
  again <- inlaST.estimate(d$Y, d$model, BPPARAM = sp, spatial = "all",
                           chunk_size = 2L, checkpoint_dir = dir)
  expect_identical(again$timing$resumed_null_chunks, 3L)
  expect_identical(again$timing$resumed_spatial_chunks, 3L)
  .ts_same_estimates(again, first, fields)

  unlink(list.files(dir, "^spatial-", full.names = TRUE)[1L])
  partial <- inlaST.estimate(d$Y, d$model, BPPARAM = sp, spatial = "all",
                             chunk_size = 2L, checkpoint_dir = dir)
  expect_identical(partial$timing$resumed_null_chunks, 3L)
  expect_identical(partial$timing$resumed_spatial_chunks, 2L)
  .ts_same_estimates(partial, first, fields)

  # Step 2 reuses the directory of step 1.
  later <- tempfile("mgcvst-two-step-later-")
  on.exit(unlink(later, recursive = TRUE), add = TRUE)
  none <- inlaST.estimate(d$Y, d$model, BPPARAM = sp, spatial = "none",
                          chunk_size = 2L, checkpoint_dir = later)
  added <- inlaST.estimate_spatial(none, d$Y, "all", BPPARAM = sp, chunk_size = 2L,
                                   checkpoint_dir = later)
  expect_length(list.files(later, "^spatial-"), 3L)
  .ts_same_estimates(added, first, fields)
  resumed <- inlaST.estimate_spatial(none, d$Y, "all", BPPARAM = sp, chunk_size = 2L,
                                     checkpoint_dir = later)
  expect_identical(resumed$timing$resumed_spatial_chunks, 3L)

  # A different offset or control is refused; a changed response recomputes.
  expect_error(inlaST.estimate(d$Y, d$model, BPPARAM = sp, spatial = "all",
    chunk_size = 2L, offset = rep(0.1, d$n), checkpoint_dir = dir),
    "different model, offset or controls")
  expect_error(inlaST.estimate(d$Y, d$model, BPPARAM = sp, spatial = "all",
    chunk_size = 2L, control = list(control.inla = list(tolerance = 1e-3)),
    checkpoint_dir = dir), "different model, offset or controls")
  expect_error(inlaST.estimate(d$Y, d$model, BPPARAM = sp, spatial = "all",
    chunk_size = 2L, checkpoint_dir = dir, resume = FALSE), "already exists")
  mgcv_dir <- tempfile("mgcvst-two-step-mgcv-")
  on.exit(unlink(mgcv_dir, recursive = TRUE), add = TRUE)
  mgcvST:::.mgcvst_chunk_store(mgcv_dir, "mgcv", "x", TRUE)
  expect_error(inlaST.estimate(d$Y, d$model, BPPARAM = sp, checkpoint_dir = mgcv_dir),
               "another estimator")
})

test_that("mu_bar is stored at estimation and read by the PCAlearning scales", {
  skip_on_cran()
  d <- .ts_inla()
  .ts_memoize_fits()
  fit <- inlaST.estimate(d$Y, d$model, BPPARAM = BiocParallel::SerialParam(),
                         spatial = "all")
  prepared <- mgcvST:::.inlast_sparse_prepare(fit)
  state <- mgcvST:::.inlast_working_state(prepared, seq_along(fit$feature_id))
  expect_equal(unname(fit$mu_bar), colMeans(state$mu), tolerance = 1e-10)
  scales <- mgcvST:::.mgcvst_pca_scales(fit)
  expect_identical(scales$mu_bar, unname(fit$mu_bar))
  nb <- fit$diagnostics$family_used == "negative_binomial"
  theta <- vapply(fit$family_parameters, function(x) x[1L], numeric(1L))
  expect_equal(scales$sigma_e2, ifelse(nb, 1 + fit$mu_bar / theta, 1), tolerance = 1e-12)

  # The scales read the stored mean: an altered mu_bar changes them, and a fit
  # that lacks it is refused instead of being recomputed silently.
  altered <- fit
  altered$mu_bar[] <- 2 * altered$mu_bar
  expect_equal(mgcvST:::.mgcvst_pca_scales(altered)$mu_bar, 2 * unname(fit$mu_bar))
  missing <- fit
  missing$mu_bar <- NULL
  expect_error(mgcvST:::.mgcvst_pca_scales(missing), "does not store mu_bar")
  for (old in list(missing, local({ x <- fit; x$format <- NULL; x }),
                   local({ x <- fit; x$format <- 1L; x }))) {
    expect_error(inlaST.test(old, rank = 2L), "before mgcvST 0.0.1.9032|re-run inlaST.estimate")
    expect_error(inlaST.wgcna(old, indices = fit$feature_id[1:3]), "re-run inlaST.estimate")
    expect_error(inlaST.estimate_spatial(old, d$Y), "re-run inlaST.estimate")
  }
})

test_that("the observation basis is full rank and enters the pair signature", {
  skip_on_cran()
  d <- .ts_inla()
  .ts_memoize_fits()
  fit <- inlaST.estimate(d$Y, d$model, BPPARAM = BiocParallel::SerialParam(),
                         spatial = "all")
  prepared <- mgcvST:::.inlast_sparse_prepare(fit)
  basis <- mgcvST:::.inlast_sparse_observation_basis(prepared)
  expect_identical(basis$kind, "full_rank")
  expect_identical(basis$rank, d$m - 1L)
  expect_identical(ncol(basis$coordinate), d$m - 1L)
  # The basis is built once and is the one the test reports.
  tested <- inlaST.test(fit, rank = 3L, seed = 4L)
  expect_identical(tested$timing$inla_projection$q, d$m)
  expect_identical(tested$timing$inla_projection$r, d$m - 1L)
  expect_identical(tested$timing$inla_projection$basis_kind, "full_rank")
  # WGCNA uses the same basis, and its normalizer is q - 1.
  scores <- mgcvST:::.mgcvst_inla_wgcna_scores(fit, 1:3, 1L, FALSE)
  expect_identical(scores$normalization, d$m - 1L)
  expect_equal(unname(scores$A),
               unname(crossprod(basis$coordinate, fit$score_a[, 1:3])),
               tolerance = 1e-12)
  # The estimation records the basis; a test or WGCNA run with another one fails.
  expect_identical(fit$basis_spec, list(kind = "full_rank", rank = d$m - 1L))
  changed <- fit
  changed$basis_spec$kind <- "truncated"
  expect_error(inlaST.test(changed, rank = 3L, seed = 4L), "differs from the one recorded")
  expect_error(mgcvST:::.mgcvst_inla_wgcna_scores(changed, 1:3, 1L, FALSE),
               "differs from the one recorded")
  expect_error(inlaST.test(local({ x <- fit; x$basis_spec <- NULL; x }), rank = 3L),
               "recorded by the estimation (none)", fixed = TRUE)
  # The pair signature, and with it every pair checkpoint, depends on the basis.
  signature <- mgcvST:::.mgcvst_pair_signature(prepared, basis)
  other <- basis
  other$kind <- "truncated"
  expect_false(identical(signature, mgcvST:::.mgcvst_pair_signature(prepared, other)))
  fewer <- basis
  fewer$rank <- basis$rank - 1L
  expect_false(identical(signature, mgcvST:::.mgcvst_pair_signature(prepared, fewer)))
})
