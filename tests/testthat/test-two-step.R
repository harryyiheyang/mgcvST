# Two-step estimation: Stage 1 null fits for every feature, spatial fits for
# the selected features only.

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

.ts_fit <- local({
  cached <- NULL
  function() {
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
    model <- spdeST.set(response ~ z + offset(offset0), data, basis, family = mgcv::nb())
    cached <<- list(Y = Y, model = model, n = n,
                    m = ncol(model$spde_spec$random[[1L]]$A))
    cached
  }
})

# Identical estimates: every field the tests read, without timings.
# Wall-clock fields (`*_seconds`) differ from run to run, also inside the
# per-feature null states (the spde.scan null fits are not memoized, they are
# deterministic apart from their timing).
.ts_strip_seconds <- function(x) {
  if (is.list(x) && !is.data.frame(x)) {
    if (!is.null(names(x))) x <- x[!grepl("_seconds$", names(x))]
    x[] <- lapply(x, .ts_strip_seconds)
  }
  x
}

.ts_same_estimates <- function(a, b, fields) {
  for (name in fields) {
    expect_identical(.ts_strip_seconds(a[[name]]), .ts_strip_seconds(b[[name]]),
                     info = name)
  }
  keep <- !grepl("_seconds$", names(a$diagnostics))
  expect_identical(a$diagnostics[, keep], b$diagnostics[, keep])
}

# ---- Stage 1 selection ------------------------------------------------------

test_that("spatial selection resolves discoveries, all, none, IDs, indices and logicals", {
  ids <- paste0("f", 1:5)
  q <- c(0.001, 0.2, NA, 0.04, 0.5)
  sel <- function(x, qv = 0.05) spdeST:::.mgcvst_select_spatial(x, ids, q, qv)
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
    expect_equal(spdeST:::.mgcvst_stage1_q(p, adjust), p.adjust(p, adjust),
                 tolerance = 1e-12)
  }
  expect_true(is.na(spdeST:::.mgcvst_stage1_q(p, "BY")[3L]))
  expect_equal(spdeST:::.mgcvst_stage1_q(c(0, 0.5), "BH"), c(0, 0.5))
  expect_error(spdeST:::.mgcvst_check_q_value(0), "q.value")
  expect_error(spdeST:::.mgcvst_check_q_value(1.5), "q.value")
  expect_error(spdeST:::.mgcvst_check_q_value(c(0.1, 0.2)), "q.value")
  expect_identical(spdeST:::.mgcvst_check_q_value(1), 1)
})

# ---- chunk checkpoints ------------------------------------------------------

test_that("chunk checkpoints are resumed, and refused when they do not match", {
  dir <- tempfile("mgcvst-chunks-")
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  store <- spdeST:::.mgcvst_chunk_store(dir, "spde", "signature-a", TRUE)
  digest <- c("a", "b", "c", "d")
  calls <- 0L
  work <- function(payload, scale) {
    calls <<- calls + length(payload$index)
    out <- as.list(payload$index * scale)
    spdeST:::.mgcvst_chunk_save(payload$chunk_file, payload$chunk_key, out)
    out
  }
  groups <- list(1:2, 3:4)
  make_payload <- function(index) list(index = index)
  sp <- BiocParallel::SerialParam()
  run <- function(step, dig, route, grp) {
    spdeST:::.mgcvst_run_chunks(grp, make_payload, step, store, dig, route,
                                sp, work, scale = 10)
  }
  digest_a <- digest
  groups_a <- groups
  first <- run("null", digest_a, NULL, groups_a)
  expect_identical(first$resumed, 0L)
  expect_identical(unlist(first$results), c(10, 20, 30, 40))
  expect_identical(calls, 4L)
  again <- run("null", digest_a, NULL, groups_a)
  expect_identical(again$resumed, 2L)
  expect_identical(calls, 4L)
  expect_identical(again$results, first$results)
  # A chunk is keyed by its features and their response digests.
  other <- run("null", c("a", "b", "x", "d"), NULL, groups_a)
  expect_identical(other$resumed, 1L)
  # ... and by the routing of its features.
  plain <- run("null", digest_a, rep(FALSE, 4L), groups_a)
  expect_identical(plain$resumed, 0L)
  routed <- run("null", digest_a, c(FALSE, FALSE, TRUE, FALSE), groups_a)
  expect_identical(routed$resumed, 1L)
  # A chunk is found by its key whatever its position in the list.
  moved <- run("null", digest_a, NULL, list(3:4, 1:2))
  expect_identical(moved$resumed, 2L)
  expect_identical(unlist(moved$results), c(30, 40, 10, 20))
  # A step has its own chunks.
  spatial <- run("spatial", digest_a, NULL, groups_a)
  expect_identical(spatial$resumed, 0L)

  # A damaged chunk is reported, not recomputed silently.
  victim <- spdeST:::.mgcvst_chunk_file(store, "null",
    spdeST:::.mgcvst_chunk_key("null", groups_a[[1L]], digest_a, NULL))
  expect_true(file.exists(victim))
  writeBin(as.raw(1:20), victim)
  expect_error(run("null", digest_a, NULL, groups_a), "damaged")
  unlink(victim)

  # The manifest ties the directory to one estimator, format and signature.
  expect_error(spdeST:::.mgcvst_chunk_store(dir, "mgcv", "signature-a", TRUE),
               "another estimator or by a version before 0.0.1.9032")
  expect_error(spdeST:::.mgcvst_chunk_store(dir, "spde", "signature-b", TRUE),
               "different model, offset or controls")
  expect_error(spdeST:::.mgcvst_chunk_store(dir, "spde", "signature-a", FALSE),
               "already exists")
  manifest <- file.path(dir, "estimation-manifest.rds")
  record <- readRDS(manifest)
  earlier <- record
  earlier$format <- 1L
  saveRDS(earlier, manifest)
  expect_error(spdeST:::.mgcvst_chunk_store(dir, "spde", "signature-a", TRUE),
               "before 0.0.1.9032")
  saveRDS(record, manifest)
  expect_identical(spdeST:::.mgcvst_chunk_store(dir, "spde", "signature-a", TRUE)$kind,
                   "spde")
  # Files left by an interrupted write are removed when a run resumes.
  stray <- file.path(dir, c("chunk-1a2b.tmp", "manifest-3c4d.tmp"))
  writeLines("x", stray[1L])
  writeLines("y", stray[2L])
  spdeST:::.mgcvst_chunk_store(dir, "spde", "signature-a", TRUE)
  expect_false(any(file.exists(stray)))
  expect_true(file.exists(manifest))
  loose <- tempfile("mgcvst-loose-")
  dir.create(loose)
  on.exit(unlink(loose, recursive = TRUE), add = TRUE)
  writeLines("x", file.path(loose, "other.txt"))
  expect_error(spdeST:::.mgcvst_chunk_store(loose, "spde", "signature-a", TRUE),
               "no manifest")
  expect_null(spdeST:::.mgcvst_chunk_store(NULL, "spde", "signature-a", TRUE))
  expect_error(spdeST:::.mgcvst_chunk_store(c("a", "b"), "spde", "s", TRUE),
               "checkpoint_dir")
})

test_that("chunk payloads are built lazily, only for the chunks that are computed", {
  dir <- tempfile("mgcvst-lazy-")
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  store <- spdeST:::.mgcvst_chunk_store(dir, "spde", "signature", TRUE)
  digest <- letters[1:6]
  groups <- list(1:2, 3:4, 5:6)
  built <- integer()
  make_payload <- function(index) {
    built <<- c(built, index[1L])
    list(index = index)
  }
  work <- function(payload, scale) {
    out <- as.list(payload$index * scale)
    spdeST:::.mgcvst_chunk_save(payload$chunk_file, payload$chunk_key, out)
    out
  }
  run <- function() {
    spdeST:::.mgcvst_run_chunks(groups, make_payload, "null", store, digest, NULL,
                                BiocParallel::SerialParam(), work, scale = 1)
  }
  first <- run()
  expect_identical(built, c(1L, 3L, 5L))
  expect_identical(first$built, 3L)
  # Nothing is built for chunks that are on disk.
  built <- integer()
  again <- run()
  expect_length(built, 0L)
  expect_identical(c(again$built, again$resumed), c(0L, 3L))
  # A lost chunk builds one payload.
  unlink(list.files(dir, "^null-", full.names = TRUE)[2L])
  built <- integer()
  third <- run()
  expect_identical(third$built, 1L)
  expect_length(built, 1L)
  expect_identical(unlist(third$results), as.numeric(1:6))

  # The payloads are made on demand when workers are free: with two workers the
  # constructor is called in the manager, one chunk at a time.
  skip_on_cran()
  sink_dir <- tempfile("mgcvst-lazy-snow-")
  on.exit(unlink(sink_dir, recursive = TRUE), add = TRUE)
  store2 <- spdeST:::.mgcvst_chunk_store(sink_dir, "spde", "signature", TRUE)
  calls <- new.env()
  calls$n <- 0L
  make2 <- function(index) {
    calls$n <- calls$n + 1L
    list(index = index)
  }
  wide <- split(seq_len(8L), rep(1:8, each = 1L))
  fun <- function(payload) {
    Sys.sleep(0.1)
    list(payload$index)
  }
  environment(fun) <- baseenv()
  res <- spdeST:::.mgcvst_run_chunks(wide, make2, "null", store2, as.character(1:8), NULL,
    BiocParallel::SnowParam(2L, type = "SOCK", progressbar = FALSE), fun)
  expect_identical(calls$n, 8L)
  expect_identical(unlist(res$results), 1:8)
})

test_that("chunk results are placed by key, whatever the order bpiterate returns them in", {
  groups <- split(1:8, rep(1:4, each = 2L))
  digest <- letters[1:8]
  work <- function(payload, scale) as.list(payload$index * scale)
  # An old BiocParallel returned results in completion order. The mock runs
  # the tasks and returns them reversed, which is what a fully reversed
  # completion looks like.
  reversed <- function(ITER, FUN, ..., REDUCE, init, BPPARAM) {
    out <- list()
    while (!is.null(payload <- ITER())) out[[length(out) + 1L]] <- FUN(payload, ...)
    out <- rev(out)
    if (missing(REDUCE)) out else Reduce(REDUCE, out, init)
  }
  res <- testthat::with_mocked_bindings(
    spdeST:::.mgcvst_run_chunks(groups, function(i) list(index = i), "spatial",
                                NULL, digest, NULL, BiocParallel::SerialParam(),
                                work, scale = 10),
    bpiterate = reversed, .package = "BiocParallel")
  expect_identical(unlist(res$results), as.numeric(1:8) * 10)
  expect_identical(lapply(res$results, unlist),
                   unname(lapply(groups, function(i) as.numeric(i) * 10)))

  # Checkpointed chunks: the file of each chunk holds the result of its own
  # features and a resumed run reads the same values back.
  dir <- tempfile("mgcvst-keyed-")
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  store <- spdeST:::.mgcvst_chunk_store(dir, "spde", "signature", TRUE)
  saving <- function(payload, scale) {
    out <- as.list(payload$index * scale)
    spdeST:::.mgcvst_chunk_save(payload$chunk_file, payload$chunk_key, out)
    out
  }
  testthat::with_mocked_bindings(
    spdeST:::.mgcvst_run_chunks(groups, function(i) list(index = i), "null",
                                store, digest, NULL, BiocParallel::SerialParam(),
                                saving, scale = 10),
    bpiterate = reversed, .package = "BiocParallel")
  again <- spdeST:::.mgcvst_run_chunks(groups, function(i) list(index = i), "null",
                                       store, digest, NULL,
                                       BiocParallel::SerialParam(), saving, scale = 10)
  expect_identical(again$resumed, 4L)
  expect_identical(unlist(again$results), as.numeric(1:8) * 10)

  # A result that cannot be matched to a dispatched chunk is an error.
  place <- spdeST:::.mgcvst_place_by_key
  ok <- list(list(key = "b", result = 2), list(key = "a", result = 1))
  expect_identical(place(ok, c("a", "b")), list(1, 2))
  expect_error(place(list(list(key = "a", result = 1)), c("a", "b")), "No result returned")
  expect_error(place(c(ok, list(list(key = "a", result = 3))), c("a", "b")), "same key")
  expect_error(place(list(list(key = "z", result = 1), ok[[2L]]), c("a", "b")), "unknown key")
  expect_error(place(list(list(result = 1), ok[[2L]]), c("a", "b")), "without the key")
  expect_error(place(list(simpleError("x"), ok[[2L]]), c("a", "b")), "without the key")
  expect_error(
    spdeST:::.mgcvst_run_chunks(list(1:2, 1:2), function(i) list(index = i),
                                "null", NULL, digest, NULL,
                                BiocParallel::SerialParam(), work, scale = 1),
    "distinct")

  # Two workers whose chunks finish in reverse order: every feature gets its
  # own result.
  skip_on_cran()
  many <- split(1:6, 1:6)
  slow <- function(payload) {
    Sys.sleep(0.12 * (7L - payload$index))
    list(payload$index * 100L)
  }
  environment(slow) <- baseenv()
  snow <- BiocParallel::SnowParam(2L, type = "SOCK", progressbar = FALSE)
  on.exit(BiocParallel::bpstop(snow), add = TRUE)
  res <- spdeST:::.mgcvst_run_chunks(many, function(i) list(index = i), "null",
                                     NULL, as.character(1:6), NULL, snow, slow)
  expect_identical(unlist(res$results), (1:6) * 100L)
  expect_identical(res$built, 6L)
})

test_that("a chunk that holds a failed feature is computed again on resume", {
  dir <- tempfile("mgcvst-failed-")
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  store <- spdeST:::.mgcvst_chunk_store(dir, "spde", "signature", TRUE)
  digest <- letters[1:4]
  groups <- list(1:2, 3:4)
  attempts <- 0L
  work <- function(payload) {
    attempts <<- attempts + 1L
    out <- lapply(payload$index, function(j) {
      if (j == 3L && attempts <= 2L) list(error = list(class = "e", message = "m", call = ""))
      else list(value = j)
    })
    spdeST:::.mgcvst_chunk_save(payload$chunk_file, payload$chunk_key, out)
    out
  }
  run <- function() spdeST:::.mgcvst_run_chunks(groups, function(i) list(index = i),
    "spatial", store, digest, NULL, BiocParallel::SerialParam(), work)
  first <- run()
  expect_identical(attempts, 2L)
  expect_false(is.null(first$results[[2L]][[1L]]$error))
  second <- run()
  expect_identical(c(second$resumed, second$built), c(1L, 1L))
  expect_identical(attempts, 3L)
  expect_null(second$results[[2L]][[1L]]$error)
  third <- run()
  expect_identical(c(third$resumed, third$built), c(2L, 0L))
  # The failed null step of an mgcv feature counts as a failure as well.
  expect_true(spdeST:::.mgcvst_chunk_has_error(list(list(marginal_error = list(message = "m")))))
  expect_false(spdeST:::.mgcvst_chunk_has_error(list(list(marginal_error = NULL, value = 1))))
})

test_that("the response matrix is checked in blocks and kept as given", {
  check <- spdeST:::.mgcvst_check_response_matrix
  Y <- matrix(as.integer(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)), 4L, 3L)
  expect_identical(check(Y, counts = TRUE), Y)
  expect_identical(typeof(check(Y)), "integer")
  D <- Y + 0
  expect_identical(check(D, counts = TRUE), D)
  bad <- D
  bad[4L, 2L] <- -1
  expect_error(check(bad, counts = TRUE), "Count responses must be non-negative integers")
  expect_identical(check(bad, counts = FALSE), bad)
  bad[4L, 2L] <- 1.5
  expect_error(check(bad, counts = TRUE), "Count responses must be non-negative integers")
  for (value in c(NA, NaN, Inf)) {
    bad[1L, 3L] <- value
    expect_error(check(bad), "finite numeric feature-by-observation matrix")
  }
  expect_error(check(matrix(numeric(), 0L, 3L)), "non-empty")
  expect_error(check(matrix("a", 2L, 2L)), "numeric")
  # The last block of a multi-block pass is checked too.
  wide <- matrix(1, 7L, 5L)
  wide[7L, 5L] <- -2
  expect_error(check(wide, counts = TRUE, block_elements = 10), "non-negative")
  expect_identical(check(matrix(1, 7L, 5L), counts = TRUE, block_elements = 10), matrix(1, 7L, 5L))
  # A block of rows is converted to double only when a payload is made.
  block <- spdeST:::.mgcvst_double_rows(Y, 2:3)
  expect_identical(typeof(block), "double")
  expect_equal(block, Y[2:3, ] + 0)
  expect_identical(spdeST:::.mgcvst_double_rows(D, 2:3), D[2:3, , drop = FALSE])
  # The digest of a row does not depend on names or on the storage type.
  named <- D
  dimnames(named) <- list(letters[1:4], LETTERS[1:3])
  expect_identical(spdeST:::.mgcvst_row_digests(named), spdeST:::.mgcvst_row_digests(D))
  expect_identical(spdeST:::.mgcvst_row_digests(Y), spdeST:::.mgcvst_row_digests(D))
})

test_that("chunk_size is one positive integer in every entry point", {
  check <- spdeST:::.mgcvst_check_chunk_size
  expect_null(check(NULL))
  expect_identical(check(3), 3L)
  expect_identical(check(3L), 3L)
  for (bad in list(0, -1, 1.5, 2.9, NA, Inf, c(1, 2), "2", 3e10)) {
    expect_error(check(bad), "chunk_size must be one positive integer")
  }
  # No truncation anywhere: the estimators and the add-later functions agree.
  expect_identical(
    names(formals(.spdest_estimate)),
    c("fit", "Y", "features", "adjust", "q.value", "BPPARAM", "chunk_size",
      "checkpoint_dir", "resume", "threads"))
})

test_that("step-1 chunks are the same whatever the number of workers when chunk_size is given", {
  features <- seq_len(11L)
  serial <- BiocParallel::SerialParam()
  many <- BiocParallel::SnowParam(5L, type = "SOCK", progressbar = FALSE)
  for (checkpoint in c(FALSE, TRUE)) {
    expect_identical(spdeST:::.mgcvst_feature_chunks(features, 3L, serial, checkpoint),
                     spdeST:::.mgcvst_feature_chunks(features, 3L, many, checkpoint))
  }
  # Without chunk_size they follow the workers, which is why a resumable run
  # passes it explicitly.
  expect_false(identical(spdeST:::.mgcvst_feature_chunks(features, NULL, serial),
                         spdeST:::.mgcvst_feature_chunks(features, NULL, many)))
  digest <- as.character(features)
  key <- function(index, route = NULL) spdeST:::.mgcvst_chunk_key("null", index, digest, route)
  expect_identical(key(1:3), key(1:3))
  expect_false(identical(key(1:3), key(1:4)))
  expect_false(identical(key(1:3, rep(FALSE, 11L)), key(1:3, c(FALSE, TRUE, rep(FALSE, 9L)))))
})

test_that("a fit estimated before the two-step estimators is refused only where it lacks data", {
  check <- spdeST:::.mgcvst_check_fit_format
  for (old in list(list(feature_id = "a"), list(format = 1L, mu_bar = c(1, 2)),
                   list(format = 2L))) {
    expect_error(check(old), "re-run spdeST.estimate")
  }
  new <- list(format = 2L, mu_bar = c(1, 2))
  expect_identical(check(new), new)
  expect_error(check(list(format = 99L)), "newer version")
})

# ---- estimation ----------------------------------------------------------

test_that("the scan step fits the null model of every feature and the spatial model of the discoveries", {
  skip_on_cran()
  d <- .ts_fit()
  sp <- BiocParallel::SerialParam()
  all <- spdeST:::.spdest_estimate(
    .spdest_scan(d$Y, d$model, BPPARAM = sp),
    d$Y, features = "all", BPPARAM = sp)
  expect_identical(all$format, 2L)
  expect_s3_class(all, "spdeST_scan")
  p <- all$diagnostics$marginal_p_value
  expect_true(all(is.finite(p)))
  expect_equal(all$diagnostics$marginal_q_value, p.adjust(p, "BY"), tolerance = 1e-12)
  expect_true(all(all$diagnostics$spatial_selected & all$diagnostics$spatial_fitted))
  expect_true(all(is.finite(all$mu_bar)))

  q <- all$diagnostics$marginal_q_value
  cut <- sqrt(sort(q)[3L] * sort(q)[4L])
  fit <- spdeST:::.spdest_estimate(
    .spdest_scan(d$Y, d$model, BPPARAM = sp, q.value = cut),
    d$Y, features = "discoveries", BPPARAM = sp)
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
  expect_identical(unname(spdeST:::.mgcvst_feature_available(fit)), q <= cut)
  for (name in c("score_a", "target_coefficients", "nuisance_coefficients")) {
    expect_identical(fit[[name]][, chosen], all[[name]][, chosen], info = name)
  }
  expect_identical(fit$mu_bar[chosen], all$mu_bar[chosen])
  expect_null(fit$working_error)
  expect_null(fit$working_variance)
  expect_false(.ts_has_long_numeric(fit$null_state, d$n))
  expect_output(print(fit), "spatial models: 3 of 6 features")

  # The test covers the fitted features; the others are not failures.
  tested <- spdeST.test(fit, rank = 3L, seed = 4L, moments = "exact")
  expect_identical(nrow(tested$results), 3L)
  expect_true(all(c(tested$results$i, tested$results$j) %in% chosen))
  expect_identical(nrow(tested$failed), 0L)
  bad <- spdeST.test(fit, pairs = rbind(c(chosen[1L], rest[1L])), rank = 3L, seed = 4L, moments = "exact")
  expect_identical(bad$results$status, 3L)
  expect_match(bad$failed$error, "not selected in step 2")

  none <- .spdest_scan(d$Y, d$model, BPPARAM = sp)
  expect_false(any(none$diagnostics$spatial_fitted))
  expect_true(all(is.na(none$score_a)))
  expect_identical(none$diagnostics$marginal_p_value, p)
  by_id <- spdeST:::.spdest_estimate(
    .spdest_scan(d$Y, d$model, BPPARAM = sp),
    d$Y, features = c("g6", "g1"), BPPARAM = sp)
  expect_identical(which(by_id$diagnostics$spatial_fitted), c(1L, 6L))
  expect_error(spdeST:::.spdest_estimate(
    .spdest_scan(d$Y, d$model),
    d$Y, features = "g99"), "unknown feature IDs")
  expect_error(spdeST:::.spdest_estimate(
    .spdest_scan(d$Y, d$model, q.value = 2),
    d$Y, features = "discoveries"), "q.value")
  expect_error(spdeST:::.spdest_estimate(
    .spdest_scan(d$Y, d$model, adjust = "holm"),
    d$Y, features = "discoveries"), "should be one of")
})

test_that("the estimation step adds spatial models to a step 1 fit", {
  skip_on_cran()
  d <- .ts_fit()
  sp <- BiocParallel::SerialParam()
  all <- spdeST:::.spdest_estimate(
    .spdest_scan(d$Y, d$model, BPPARAM = sp),
    d$Y, features = "all", BPPARAM = sp)
  none <- .spdest_scan(d$Y, d$model, BPPARAM = sp)
  two <- .spdest_estimate(none, d$Y, c("g2", "g5"), BPPARAM = sp)
  expect_false(any(none$diagnostics$spatial_fitted))
  expect_identical(which(two$diagnostics$spatial_fitted), c(2L, 5L))
  expect_false(is.null(two$estimation))
  full <- .spdest_estimate(two, d$Y, "all", BPPARAM = sp)
  fields <- c("score_a", "target_coefficients", "nuisance_coefficients", "mu_bar",
              "dispersion", "lambda", "component_lambda", "smoothing_parameters",
              "family_parameters", "null_state", "constraint_residual",
              "observation_spatial_mean", "feature_family", "y_digest")
  .ts_same_estimates(full, all, fields)
  expect_identical(spdeST.test(full, rank = 3L, seed = 4L, moments = "exact")$results,
                   spdeST.test(all, rank = 3L, seed = 4L, moments = "exact")$results)
  q <- all$diagnostics$marginal_q_value
  cut <- sqrt(sort(q)[3L] * sort(q)[4L])
  some <- .spdest_estimate(none, d$Y, q.value = cut, BPPARAM = sp)
  expect_identical(some$diagnostics$spatial_fitted, q <= cut)

  expect_error(.spdest_estimate(none, d$Y[, rev(seq_len(ncol(d$Y)))], "g1"),
               "Y differs from the responses of step 1 for g1")
  expect_error(.spdest_estimate(none, d$Y[, -1L], "g1"), "Y must be the feature")
  expect_error(.spdest_estimate(none, d$Y, "g99"), "unknown feature IDs")
  expect_error(.spdest_estimate(list(), d$Y), "must be returned by spdeST.scan")
})

test_that("workers return compact results and a checkpoint resumes exactly", {
  skip_on_cran()
  d <- .ts_fit()
  sp <- BiocParallel::SerialParam()
  dir <- tempfile("spdest-two-step-")
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  first <- spdeST:::.spdest_estimate(
    .spdest_scan(d$Y, d$model, BPPARAM = sp,
                           chunk_size = 2L, checkpoint_dir = dir),
    d$Y, features = "all", BPPARAM = sp, chunk_size = 2L, checkpoint_dir = dir)
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
  again <- spdeST:::.spdest_estimate(
    .spdest_scan(d$Y, d$model, BPPARAM = sp,
                           chunk_size = 2L, checkpoint_dir = dir),
    d$Y, features = "all", BPPARAM = sp, chunk_size = 2L, checkpoint_dir = dir)
  expect_identical(again$timing$resumed_null_chunks, 3L)
  expect_identical(again$timing$resumed_spatial_chunks, 3L)
  .ts_same_estimates(again, first, fields)

  unlink(list.files(dir, "^spatial-", full.names = TRUE)[1L])
  partial <- spdeST:::.spdest_estimate(
    .spdest_scan(d$Y, d$model, BPPARAM = sp,
                             chunk_size = 2L, checkpoint_dir = dir),
    d$Y, features = "all", BPPARAM = sp, chunk_size = 2L, checkpoint_dir = dir)
  expect_identical(partial$timing$resumed_null_chunks, 3L)
  expect_identical(partial$timing$resumed_spatial_chunks, 2L)
  .ts_same_estimates(partial, first, fields)

  # Step 2 reuses the directory of step 1.
  later <- tempfile("mgcvst-two-step-later-")
  on.exit(unlink(later, recursive = TRUE), add = TRUE)
  none <- .spdest_scan(d$Y, d$model, BPPARAM = sp,
                          chunk_size = 2L, checkpoint_dir = later)
  added <- .spdest_estimate(none, d$Y, "all", BPPARAM = sp, chunk_size = 2L,
                                   checkpoint_dir = later)
  expect_length(list.files(later, "^spatial-"), 3L)
  .ts_same_estimates(added, first, fields)
  resumed <- .spdest_estimate(none, d$Y, "all", BPPARAM = sp, chunk_size = 2L,
                                     checkpoint_dir = later)
  expect_identical(resumed$timing$resumed_spatial_chunks, 3L)

  # A different offset or control is refused; a changed response recomputes.
  expect_error(spdeST:::.spdest_estimate(
    .spdest_scan(d$Y, d$model, BPPARAM = sp,
    chunk_size = 2L, offset = rep(0.1, d$n), checkpoint_dir = dir),
    d$Y, features = "all", BPPARAM = sp, chunk_size = 2L, checkpoint_dir = dir),
    "different model, offset or controls")
  expect_error(spdeST:::.spdest_estimate(
    .spdest_scan(d$Y, d$model, BPPARAM = sp,
    chunk_size = 2L, control = list(spde = list(tol = 1e-6)),
    checkpoint_dir = dir),
    d$Y, features = "all", BPPARAM = sp, chunk_size = 2L, checkpoint_dir = dir), "different model, offset or controls")
  expect_error(spdeST:::.spdest_estimate(
    .spdest_scan(d$Y, d$model, BPPARAM = sp,
    chunk_size = 2L, checkpoint_dir = dir, resume = FALSE),
    d$Y, features = "all", BPPARAM = sp, chunk_size = 2L, checkpoint_dir = dir, resume = FALSE), "already exists")
  mgcv_dir <- tempfile("mgcvst-two-step-mgcv-")
  on.exit(unlink(mgcv_dir, recursive = TRUE), add = TRUE)
  spdeST:::.mgcvst_chunk_store(mgcv_dir, "mgcv", "x", TRUE)
  expect_error(spdeST:::.spdest_estimate(
    .spdest_scan(d$Y, d$model, BPPARAM = sp, checkpoint_dir = mgcv_dir),
    d$Y, features = "discoveries", BPPARAM = sp, checkpoint_dir = mgcv_dir),
               "another estimator")
})

test_that("mu_bar is stored at estimation and read by the PCAlearning scales", {
  skip_on_cran()
  d <- .ts_fit()
  fit <- spdeST:::.spdest_estimate(
    .spdest_scan(d$Y, d$model, BPPARAM = BiocParallel::SerialParam()),
    d$Y, features = "all", BPPARAM = BiocParallel::SerialParam())
  scales <- spdeST:::.mgcvst_pca_scales(fit)
  expect_identical(scales$mu_bar, unname(fit$mu_bar))
  nb <- fit$diagnostics$family_used == "negative_binomial"
  theta <- vapply(fit$family_parameters, function(x) x[1L], numeric(1L))
  expect_equal(scales$sigma_e2, ifelse(nb, 1 + fit$mu_bar / theta, 1), tolerance = 1e-12)

  # The scales read the stored mean: an altered mu_bar changes them, and a fit
  # that lacks it is refused instead of being recomputed silently.
  altered <- fit
  altered$mu_bar[] <- 2 * altered$mu_bar
  expect_equal(spdeST:::.mgcvst_pca_scales(altered)$mu_bar, 2 * unname(fit$mu_bar))
  missing <- fit
  missing$mu_bar <- NULL
  expect_error(spdeST:::.mgcvst_pca_scales(missing), "does not store mu_bar")
  for (old in list(missing, local({ x <- fit; x$format <- NULL; x }),
                   local({ x <- fit; x$format <- 1L; x }))) {
    expect_error(spdeST.test(old, rank = 2L, moments = "exact"), "re-run spdeST.estimate")
    expect_error(.spdest_estimate(old, d$Y), "re-run spdeST.estimate")
  }
})

test_that("the observation basis is full rank and enters the pair signature", {
  skip_on_cran()
  d <- .ts_fit()
  fit <- spdeST:::.spdest_estimate(
    .spdest_scan(d$Y, d$model, BPPARAM = BiocParallel::SerialParam()),
    d$Y, features = "all", BPPARAM = BiocParallel::SerialParam())
  prepared <- spdeST:::.spdest_sparse_prepare(fit)
  basis <- spdeST:::.spdest_sparse_observation_basis(prepared)
  expect_identical(basis$kind, "full_rank")
  expect_identical(basis$rank, d$m - 1L)
  expect_identical(ncol(basis$coordinate), d$m - 1L)
  # The basis is built once and is the one the test reports.
  tested <- spdeST.test(fit, rank = 3L, seed = 4L, moments = "exact")
  expect_identical(tested$timing$score_projection$q, d$m)
  expect_identical(tested$timing$route$q, d$m - 1L)
  expect_identical(tested$timing$score_projection$r, d$m - 1L)
  expect_identical(tested$timing$score_projection$basis_kind, "full_rank")
  # The estimation records the basis; a test with another one fails.
  expect_identical(fit$basis_spec, list(kind = "full_rank", rank = d$m - 1L))
  changed <- fit
  changed$basis_spec$kind <- "truncated"
  expect_error(spdeST.test(changed, rank = 3L, seed = 4L, moments = "exact"), "differs from the one recorded")
  expect_error(spdeST.test(local({ x <- fit; x$basis_spec <- NULL; x }), rank = 3L, moments = "exact"),
               "recorded by the estimation (none)", fixed = TRUE)
  # The pair signature, and with it every pair checkpoint, depends on the basis.
  signature <- spdeST:::.mgcvst_pair_signature(prepared, basis)
  other <- basis
  other$kind <- "truncated"
  expect_false(identical(signature, spdeST:::.mgcvst_pair_signature(prepared, other)))
  fewer <- basis
  fewer$rank <- basis$rank - 1L
  expect_false(identical(signature, spdeST:::.mgcvst_pair_signature(prepared, fewer)))
})

# ---- resume, lazy payloads and multi-process runs ----------------------------

.ts_fields_mgcv <- c("working_error", "working_variance", "dispersion", "lambda",
                     "smoothing_parameters", "nuisance_covariance", "family_parameters")
.ts_fields <- c("score_a", "target_coefficients", "nuisance_coefficients",
                     "mu_bar", "dispersion", "lambda", "smoothing_parameters",
                     "family_parameters", "null_state")

# Replace the first feature of a saved chunk by a failure record.
.ts_fail_chunk <- function(file) {
  z <- readRDS(file)
  first <- z$result[[1L]]
  z$result[[1L]] <- list(error = list(class = "simpleError", message = "transient",
                                      call = ""),
                         index = first$index, feature_id = first$feature_id)
  saveRDS(z, file)
}

test_that("a resumed run builds no payload, repeats failed chunks and cleans strays", {
  skip_on_cran()
  d <- .ts_fit()
  sp <- BiocParallel::SerialParam()
  dir <- tempfile("spdest-resume-")
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  estimate <- function() spdeST:::.spdest_estimate(
    .spdest_scan(d$Y, d$model, BPPARAM = sp, chunk_size = 2L,
                    checkpoint_dir = dir),
    d$Y, features = "all", BPPARAM = sp, chunk_size = 2L, checkpoint_dir = dir)
  first <- estimate()
  expect_identical(c(first$timing$null_payloads, first$timing$spatial_payloads), c(3L, 3L))
  stray <- file.path(dir, "chunk-dead.tmp")
  writeLines("x", stray)
  again <- estimate()
  expect_false(file.exists(stray))
  expect_identical(c(again$timing$null_payloads, again$timing$spatial_payloads), c(0L, 0L))
  .ts_same_estimates(again, first, .ts_fields)

  .ts_fail_chunk(list.files(dir, "^spatial-", full.names = TRUE)[2L])
  retried <- estimate()
  expect_identical(retried$timing$resumed_spatial_chunks, 2L)
  expect_identical(retried$timing$spatial_payloads, 1L)
  expect_true(all(retried$diagnostics$spatial_fitted))
  .ts_same_estimates(retried, first, .ts_fields)

  .ts_fail_chunk(list.files(dir, "^null-", full.names = TRUE)[1L])
  null_again <- estimate()
  expect_identical(null_again$timing$resumed_null_chunks, 2L)
  expect_true(all(is.finite(null_again$diagnostics$marginal_p_value)))
  .ts_same_estimates(null_again, first, .ts_fields)
})

test_that("another selection in the same directory recomputes step 2 and never mixes", {
  skip_on_cran()
  d <- .ts_fit()
  sp <- BiocParallel::SerialParam()
  all <- spdeST:::.spdest_estimate(
    .spdest_scan(d$Y, d$model, BPPARAM = sp),
    d$Y, features = "all", BPPARAM = sp)
  dir <- tempfile("spdest-selection-")
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  estimate <- function(spatial) spdeST:::.spdest_estimate(
    .spdest_scan(d$Y, d$model, BPPARAM = sp, chunk_size = 1L,
                    checkpoint_dir = dir),
    d$Y, features = spatial, BPPARAM = sp, chunk_size = 1L, checkpoint_dir = dir)
  one <- estimate(c("g1", "g2"))
  expect_identical(which(one$diagnostics$spatial_fitted), 1:2)
  two <- estimate(c("g2", "g3", "g1"))
  expect_identical(two$timing$resumed_null_chunks, 6L)
  expect_identical(two$timing$resumed_spatial_chunks, 2L)
  expect_identical(two$timing$spatial_payloads, 1L)
  expect_identical(which(two$diagnostics$spatial_fitted), 1:3)
  expect_identical(two$score_a[, 1:3], all$score_a[, 1:3])
  expect_true(all(is.na(two$score_a[, 4:6])))
  q <- all$diagnostics$marginal_q_value
  cut <- sqrt(sort(q)[2L] * sort(q)[3L])
  dir2 <- tempfile("spdest-qvalue-")
  on.exit(unlink(dir2, recursive = TRUE), add = TRUE)
  small <- spdeST:::.spdest_estimate(
    .spdest_scan(d$Y, d$model, BPPARAM = sp, q.value = cut, chunk_size = 1L,
                           checkpoint_dir = dir2),
    d$Y, features = "discoveries", BPPARAM = sp, chunk_size = 1L, checkpoint_dir = dir2)
  large <- spdeST:::.spdest_estimate(
    .spdest_scan(d$Y, d$model, BPPARAM = sp, q.value = 1, chunk_size = 1L,
                           checkpoint_dir = dir2),
    d$Y, features = "discoveries", BPPARAM = sp, chunk_size = 1L, checkpoint_dir = dir2)
  expect_identical(sum(small$diagnostics$spatial_fitted), 2L)
  expect_true(all(large$diagnostics$spatial_fitted))
  expect_identical(large$timing$resumed_spatial_chunks, 2L)
  expect_identical(large$score_a, all$score_a)
})

test_that("Y may be integer or unnamed, and chunk_size is validated everywhere", {
  skip_on_cran()
  d <- .ts_fit()
  sp <- BiocParallel::SerialParam()
  reference <- spdeST:::.spdest_estimate(
    .spdest_scan(d$Y, d$model, BPPARAM = sp),
    d$Y, features = "all", BPPARAM = sp)
  integer_Y <- d$Y
  storage.mode(integer_Y) <- "integer"
  from_integer <- spdeST:::.spdest_estimate(
    .spdest_scan(integer_Y, d$model, BPPARAM = sp),
    integer_Y, features = "all", BPPARAM = sp)
  .ts_same_estimates(from_integer, reference, .ts_fields)
  named <- d$Y
  colnames(named) <- paste0("obs", seq_len(ncol(named)))
  none <- .spdest_scan(named, d$model, BPPARAM = sp)
  added <- .spdest_estimate(none, unname(d$Y), "all", BPPARAM = sp)
  .ts_same_estimates(added, reference, .ts_fields)
  added_integer <- .spdest_estimate(none, integer_Y, "all", BPPARAM = sp)
  .ts_same_estimates(added_integer, reference, .ts_fields)
  expect_error(spdeST:::.spdest_estimate(
    .spdest_scan(d$Y, d$model, chunk_size = 1.5),
    d$Y, features = "discoveries", chunk_size = 1.5),
               "chunk_size must be one positive integer")
  expect_error(.spdest_estimate(none, d$Y, "all", chunk_size = 2.5),
               "chunk_size must be one positive integer")
  bad <- d$Y
  bad[2L, 3L] <- -1
  expect_error(spdeST:::.spdest_estimate(
    .spdest_scan(bad, d$model),
    bad, features = "discoveries"), "Count responses must be non-negative integers")
  bad[2L, 3L] <- 0.5
  expect_error(spdeST:::.spdest_estimate(
    .spdest_scan(bad, d$model),
    bad, features = "discoveries"), "Count responses must be non-negative integers")
  bad[2L, 3L] <- NA
  expect_error(spdeST:::.spdest_estimate(
    .spdest_scan(bad, d$model),
    bad, features = "discoveries"), "finite numeric feature-by-observation matrix")
})

test_that("a multi-process run with checkpoint_dir resumes under another backend", {
  skip_on_cran()
  d <- .ts_fit()
  Y <- d$Y[1:3, , drop = FALSE]
  dir <- tempfile("spdest-snow-")
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  snow <- BiocParallel::SnowParam(2L, type = "SOCK", progressbar = FALSE)
  first <- spdeST:::.spdest_estimate(
    .spdest_scan(Y, d$model, BPPARAM = snow, chunk_size = 1L,
                           checkpoint_dir = dir),
    Y, features = "all", BPPARAM = snow, chunk_size = 1L, checkpoint_dir = dir)
  expect_length(list.files(dir, "^null-"), 3L)
  expect_length(list.files(dir, "^spatial-"), 3L)
  expect_length(list.files(dir, "[.]tmp$"), 0L)
  expect_true(all(first$diagnostics$spatial_fitted))
  again <- spdeST:::.spdest_estimate(
    .spdest_scan(Y, d$model, BPPARAM = BiocParallel::SerialParam(), chunk_size = 1L, checkpoint_dir = dir),
    Y, features = "all", BPPARAM = BiocParallel::SerialParam(), chunk_size = 1L, checkpoint_dir = dir)
  expect_identical(c(again$timing$resumed_null_chunks, again$timing$resumed_spatial_chunks),
                   c(3L, 3L))
  expect_identical(c(again$timing$null_payloads, again$timing$spatial_payloads), c(0L, 0L))
  .ts_same_estimates(again, first, .ts_fields)
})
