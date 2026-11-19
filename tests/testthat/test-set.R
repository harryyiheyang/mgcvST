test_that("set expands factor interactions once for shared BAM null and full designs", {
  for (family in list(mgcv::nb())) {
    f <- st_fixture(family = family)
    d <- f$data
    d$a <- factor(rep(1:3, length.out = nrow(d)))
    d$b <- factor(rep(1:2, length.out = nrow(d)))
    d$x <- d$x + 1e-6
    d$y <- d$y - 5e-7
    d$w <- runif(nrow(d))
    basis <- f$basis
    saved <- serialize(basis, NULL)
    form <- response ~ a * z + b + w + offset(offset0) +
      s(x, y, bs = "spde", xt = basis)
    model <- mgcvST.set(form, d, family)
    expect_true(model$shared_design)
    expect_identical(serialize(basis, NULL), saved)
    external <- mgcvST.set(G = model$G)
    expect_equal(external$L, model$L, tolerance = 1e-8)
    expect_true(all(unlist(model$timing) >= 0))
    expect_false(grepl("spde", paste(deparse(model$null_formula), collapse = " ")))
    sm <- model$G$smooth[[1L]]
    expect_true(sm$timing$basis_calls >= 1L)
    expected <- spdeST:::.spde_basis_at(basis, as.matrix(d[, c("x", "y")]))
    expect_equal(unname(model$L[, sm$first.para:sm$last.para]), expected, tolerance = 1e-8)
    offset <- matrix(seq(-.15, .2, length.out = length(f$Y)), nrow(f$Y))
    fit <- testthat::with_mocked_bindings(
      spdeST:::mgcvST.estimate_spatial(
    mgcvST.estimate(f$Y, model, offset = offset,
        BPPARAM = BiocParallel::SerialParam(), chunk_size = 2L,
        # This test verifies the frozen design against direct mgcv fits with
        # the requested family; disable the Poisson prescreen so near-Poisson
        # synthetic genes are not routed away from nb().
        control = local({
          k <- mgcv::gam.control(nthreads = 1L)
          k$ncv.threads <- 1L
          k
        })),
    f$Y, features = "all", BPPARAM = BiocParallel::SerialParam(), chunk_size = 2L),
      .gam_training_lpmatrix = function(...) stop("L rebuilt during estimate"),
      .mgcvst_prescreen_route = function(Y, X, offset, active, threads = 1L) {
        p <- nrow(as.matrix(Y))
        list(phi = rep(NA_real_, p), score = rep(NA_real_, p),
             poisson = rep(FALSE, p))
      },
      .package = "spdeST"
    )
    expect_true(all(is.finite(fit$diagnostics$marginal_p_value)))
    expect_equal(fit$offset, sweep(offset, 2L, model$offset, "+"))
    family_raw <- serialize(model$G$family, NULL)
    for (i in seq_len(nrow(f$Y))) {
      full_data <- model$full_data
      full_data[[model$null_response]] <- f$Y[i, ]
      control <- mgcv::gam.control(nthreads = 1L)
      control$ncv.threads <- 1L
      direct <- mgcv::bam(formula = model$full_formula, data = full_data,
        family = unserialize(family_raw), offset = offset[i, ], method = "fREML",
        discrete = TRUE, nthreads = 1L, control = control)
      W <- rkhs_extract_working_model(direct)
      expect_equal(fit$working_error[, i], unname(W$working_error), tolerance = 1e-7)
      expect_equal(fit$working_variance[, i], unname(W$working_variance), tolerance = 1e-7)
    }
    pair <- mgcvST.test(fit, pairs = matrix(c(1L, 2L), 1L), moments = "exact")
    expect_true(all(is.finite(pair$results$log_p_two_sided)))
  }
})

test_that("set rejects gene-specific designs and invalid offsets", {
  f <- st_fixture()
  model <- mgcvST.set(G = f$G)
  expect_error(mgcvST.set(G = f$G, data = f$data), "Supply G alone")
  expect_error(mgcvST.set(response ~ response + x, f$data), "response cannot")
  expect_error(spdeST:::mgcvST.estimate_spatial(
    mgcvST.estimate(f$Y, model, data = f$data),
    f$Y, features = "all"), "Do not supply")
  expect_true(all(is.finite(spdeST:::mgcvST.estimate_spatial(
    mgcvST.estimate(
    f$Y, model, BPPARAM = BiocParallel::SerialParam()),
    f$Y, features = "all", BPPARAM = BiocParallel::SerialParam())$working_error)))
  expect_error(spdeST:::mgcvST.estimate_spatial(
    mgcvST.estimate(f$Y, model, offset = matrix(0, 2, 3)),
    f$Y, features = "all"), "offset must")
  expect_error(spdeST:::mgcvST.estimate_spatial(
    mgcvST.estimate(f$Y, model, offset = rep(NA_real_, ncol(f$Y))),
    f$Y, features = "all"), "offset must")
})

test_that("null formula keeps leading tensor nuisance smooths", {
  for (term in list(quote(te(z, x)), quote(ti(z, x)), quote(t2(z, x)))) {
    full <- spdeST:::.mgcvst_rebuild_formula(
      response ~ 1, list(term, quote(s(x, y, bs = "spde", xt = basis)))
    )
    null <- spdeST:::.mgcvst_null_formula(full, 2L)
    expect_match(paste(deparse(null), collapse = " "), as.character(term[[1L]]))
    expect_false(grepl("spde", paste(deparse(null), collapse = " ")))
  }
})

test_that("both constructors and predictions evaluate every supplied coordinate set", {
  f <- st_fixture()
  basis <- f$basis
  s <- mgcv::s
  for (delta in c(0, 1e-12, 1e-6)) {
    d <- f$data
    d$x <- d$x + delta
    form <- response ~ s(x, y, bs = "spde", xt = basis)
    G <- mgcv::gam(form, data = d, fit = FALSE)
    fit <- mgcv::gam(G = G)
    L <- spdeST:::.gam_training_lpmatrix(fit)
    expect_equal(as.numeric(G$X), as.numeric(L), tolerance = 1e-14)
    sm <- fit$smooth[[1L]]
    calls <- sm$timing$prediction_calls
    invisible(mgcv::PredictMat(sm, d))
    expect_gt(sm$timing$prediction_calls, calls)
    expect_gte(sm$timing$prediction_seconds, 0)
    expect_gte(sm$timing$basis_seconds, 0)
  }
})

test_that("set shares L with SOCK workers and gene offsets", {
  f <- st_fixture()
  model <- mgcvST.set(G = f$G)
  offset <- matrix(seq(-.1, .1, length.out = length(f$Y)), nrow(f$Y))
  a <- spdeST:::mgcvST.estimate_spatial(
    mgcvST.estimate(f$Y, model, offset = offset,
    BPPARAM = BiocParallel::SerialParam()),
    f$Y, features = "all", BPPARAM = BiocParallel::SerialParam())
  b <- spdeST:::mgcvST.estimate_spatial(
    mgcvST.estimate(f$Y, model, offset = offset,
    BPPARAM = BiocParallel::SnowParam(2L, type = "SOCK")),
    f$Y, features = "all", BPPARAM = BiocParallel::SnowParam(2L, type = "SOCK"))
  expect_equal(a$working_error, b$working_error, tolerance = 1e-12)
  expect_equal(a$nuisance_covariance, b$nuisance_covariance, tolerance = 1e-12)
  expect_equal(a$diagnostics$marginal_p_value, b$diagnostics$marginal_p_value, tolerance = 1e-12)
  restored <- unserialize(serialize(model, NULL))
  c <- spdeST:::mgcvST.estimate_spatial(
    mgcvST.estimate(f$Y, restored, offset = offset,
    BPPARAM = BiocParallel::SerialParam()),
    f$Y, features = "all", BPPARAM = BiocParallel::SerialParam())
  expect_equal(a$working_error, c$working_error, tolerance = 1e-12)
})
