test_that("set expands factor interactions once for shared BAM null and full designs", {
  for (family in list(gaussian(), mgcv::nb())) for (pc in c(FALSE, TRUE)) {
    f <- st_fixture(family = family)
    d <- f$data
    d$a <- factor(rep(1:3, length.out = nrow(d)))
    d$b <- factor(rep(1:2, length.out = nrow(d)))
    d$x <- d$x + 1e-6
    d$y <- d$y - 5e-7
    basis <- f$basis
    saved <- serialize(basis, NULL)
    form <- if (pc) response ~ a * z + b + offset(offset0) + s(z, k = 5) +
      s(x, y, bs = "spdePC", xt = basis) else
      response ~ a * z + b + offset(offset0) + s(z, k = 5) + s(x, y, bs = "spde", xt = basis)
    model <- mgcvST.set(form, d, family)
    expect_true(model$shared_design)
    expect_identical(serialize(basis, NULL), saved)
    external <- mgcvST.set(G = model$G)
    expect_equal(external$L, model$L, tolerance = 1e-8)
    expect_true(all(unlist(model$timing) >= 0))
    expect_false(grepl("spde", paste(deparse(model$null_formula), collapse = " ")))
    expect_match(paste(deparse(model$null_formula), collapse = " "), "s\\(z")
    sm <- model$G$smooth[[2L]]
    expect_true(sm$timing$basis_calls >= 1L)
    expected <- mgcvST:::.spde_basis_at(basis, as.matrix(d[, c("x", "y")]), pc)
    expect_equal(unname(model$L[, sm$first.para:sm$last.para]), expected, tolerance = 1e-8)
    if (pc) expect_equal(sm$score_basis,
      mgcvST:::.spde_basis_at(basis, as.matrix(d[, c("x", "y")])) , tolerance = 1e-8)
    offset <- matrix(seq(-.15, .2, length.out = length(f$Y)), nrow(f$Y))
    fit <- testthat::with_mocked_bindings(
      mgcvST.estimate(f$Y, model, offset = offset,
        BPPARAM = BiocParallel::SerialParam(), chunk_size = 2L,
        # This test verifies the frozen design against direct mgcv fits with
        # the requested family; disable the Poisson prescreen so near-Poisson
        # synthetic genes are not routed away from nb().
        control = local({
          k <- mgcv::gam.control(nthreads = 1L)
          k$ncv.threads <- 1L
          k$poisson_screen_phi <- 0
          k
        }),
        marginal_args = list(method = "liu"), retain_marginal = TRUE),
      .gam_training_lpmatrix = function(...) stop("L rebuilt during estimate"),
      .package = "mgcvST"
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
      null_data <- model$null_data
      null_data[[model$null_response]] <- f$Y[i, ]
      null_fit <- mgcv::bam(
        formula = model$null_formula, data = null_data,
        family = unserialize(family_raw), offset = offset[i, ],
        method = "fREML", discrete = TRUE, nthreads = 1L, control = control
      )
      setup <- mgcvST:::.mgcvst_null_score_setup(
        model$G, which(vapply(model$G$smooth, function(s) identical(s$score.component, "global"), logical(1L))),
        list(formula = model$null_formula, data = model$null_data,
             response = model$null_response, X0 = model$null_X)
      )
      p <- mgcvST:::.mgcvst_null_score_test(null_fit, setup, method = "liu")
      expect_equal(unname(fit$diagnostics$marginal_p_value[i]), p$smooth.pvalue, tolerance = 1e-7)
    }
    retained <- mgcvST.marginal(fit, calibration = "liu", BPPARAM = BiocParallel::SerialParam())
    expect_equal(retained$p_value, fit$diagnostics$marginal_p_value, tolerance = 1e-7)
    pair <- mgcvST.test(fit, pairs = matrix(c(1L, 2L), 1L), calibration = "liu")
    expect_true(all(is.finite(pair$results$p_two_sided)))
  }
})

test_that("set rejects gene-specific designs and invalid offsets", {
  f <- st_fixture()
  model <- mgcvST.set(G = f$G)
  expect_error(mgcvST.set(G = f$G, data = f$data), "Supply G alone")
  expect_error(mgcvST.set(response ~ response + x, f$data), "response cannot")
  expect_error(mgcvST.estimate(f$Y, model, data = f$data), "Do not supply")
  expect_true(all(is.finite(mgcvST.estimate(
    f$Y, model, worker_init = function() NULL,
    BPPARAM = BiocParallel::SerialParam()
  )$working_error)))
  expect_error(mgcvST.estimate(f$Y, model, offset = matrix(0, 2, 3)), "offset must")
  expect_error(mgcvST.estimate(f$Y, model, offset = rep(NA_real_, ncol(f$Y))), "offset must")
})

test_that("null formula keeps leading tensor nuisance smooths", {
  for (term in list(quote(te(z, x)), quote(ti(z, x)), quote(t2(z, x)))) {
    full <- mgcvST:::.mgcvst_rebuild_formula(
      response ~ 1, list(term, quote(s(x, y, bs = "spde", xt = basis)))
    )
    null <- mgcvST:::.mgcvst_null_formula(full, 2L)
    expect_match(paste(deparse(null), collapse = " "), as.character(term[[1L]]))
    expect_false(grepl("spde", paste(deparse(null), collapse = " ")))
  }
})

test_that("both constructors and predictions evaluate every supplied coordinate set", {
  for (pc in c(FALSE, TRUE)) {
    f <- st_fixture(pc = pc)
    basis <- f$basis
    s <- mgcv::s
    for (delta in c(0, 1e-12, 1e-6)) {
      d <- f$data
      d$x <- d$x + delta
      form <- if (pc) response ~ s(x, y, bs = "spdePC", xt = basis) else
        response ~ s(x, y, bs = "spde", xt = basis)
      G <- mgcv::gam(form, data = d, fit = FALSE)
      fit <- mgcv::gam(G = G)
      L <- mgcvST:::.gam_training_lpmatrix(fit)
      expect_equal(as.numeric(G$X), as.numeric(L), tolerance = 1e-14)
      sm <- fit$smooth[[1L]]
      calls <- sm$timing$prediction_calls
      invisible(mgcv::PredictMat(sm, d))
      expect_gt(sm$timing$prediction_calls, calls)
      expect_gte(sm$timing$prediction_seconds, 0)
      expect_gte(sm$timing$basis_seconds, 0)
    }
  }
})

test_that("set shares L with SOCK workers and gene offsets", {
  f <- st_fixture(pc = TRUE)
  model <- mgcvST.set(G = f$G)
  offset <- matrix(seq(-.1, .1, length.out = length(f$Y)), nrow(f$Y))
  a <- mgcvST.estimate(f$Y, model, offset = offset,
    marginal_args = list(method = "liu"), BPPARAM = BiocParallel::SerialParam())
  b <- mgcvST.estimate(f$Y, model, offset = offset,
    marginal_args = list(method = "liu"), BPPARAM = BiocParallel::SnowParam(2L, type = "SOCK"))
  expect_equal(a$working_error, b$working_error, tolerance = 1e-12)
  expect_equal(a$nuisance_covariance, b$nuisance_covariance, tolerance = 1e-12)
  expect_equal(a$diagnostics$marginal_p_value, b$diagnostics$marginal_p_value, tolerance = 1e-12)
  restored <- unserialize(serialize(model, NULL))
  c <- mgcvST.estimate(f$Y, restored, offset = offset,
    marginal_args = list(method = "liu"), BPPARAM = BiocParallel::SerialParam())
  expect_equal(a$working_error, c$working_error, tolerance = 1e-12)
})
