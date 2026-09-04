test_that("set freezes a formal design with two factors, nuisance and changed coordinates", {
  for (family in list(gaussian(), mgcv::nb())) for (pc in c(FALSE, TRUE)) {
    f <- st_fixture(family = family)
    d <- f$data
    d$a <- factor(rep(1:3, length.out = nrow(d)))
    d$b <- factor(rep(1:2, length.out = nrow(d)))
    d$x <- d$x + 1e-6
    d$y <- d$y - 5e-7
    basis <- f$basis
    saved <- serialize(basis, NULL)
    form <- if (pc) response ~ a + b + offset(offset0) + s(z, k = 5) +
      s(x, y, bs = "spdePC", xt = basis) else
      response ~ a + b + offset(offset0) + s(z, k = 5) + s(x, y, bs = "spde", xt = basis)
    model <- mgcvST.set(form, d, family)
    expect_true(model$shared_design)
    expect_identical(serialize(basis, NULL), saved)
    external <- mgcvST.set(G = model$G)
    expect_identical(external$L, model$L)
    expect_true(all(unlist(model$timing) >= 0))
    sm <- model$G$smooth[[2L]]
    expect_true(sm$timing$basis_calls >= 1L)
    expect_true(sm$timing$prediction_calls >= 1L)
    expected <- mgcvST:::.spde_basis_at(basis, as.matrix(d[, c("x", "y")]), pc)
    expect_equal(unname(model$L[, sm$first.para:sm$last.para]), expected, tolerance = 1e-12)
    if (pc) expect_equal(sm$score_basis,
      mgcvST:::.spde_basis_at(basis, as.matrix(d[, c("x", "y")])) , tolerance = 1e-12)
    offset <- matrix(seq(-.15, .2, length.out = length(f$Y)), nrow(f$Y))
    fit <- testthat::with_mocked_bindings(
      mgcvST.estimate(f$Y, model, offset = offset,
        BPPARAM = BiocParallel::SerialParam(), chunk_size = 2L,
        marginal_args = list(method = "liu"), retain_marginal = TRUE),
      .gam_training_lpmatrix = function(...) stop("L rebuilt during estimate"),
      .package = "mgcvST"
    )
    expect_true(all(is.finite(fit$diagnostics$marginal_p_value)))
    expect_equal(fit$offset, sweep(offset, 2L, model$offset, "+"))
    expect_identical(fit$geometry$nuisance_projection, "conditional_Vp_block")
    expect_true(all(vapply(fit$nuisance_covariance, is.matrix, logical(1))))
    family_raw <- serialize(model$G$family, NULL)
    for (i in seq_len(nrow(f$Y))) {
      G <- model$G
      G$family <- unserialize(family_raw)
      G$y <- f$Y[i, ]
      G$mf[[1L]] <- as.numeric(G$y)
      G$offset <- model$offset + offset[i, ]
      control <- mgcv::gam.control(nthreads = 1L)
      control$ncv.threads <- 1L
      direct <- mgcv::gam(G = G, method = "REML", control = control)
      expect_identical(model$L, mgcvST:::.gam_training_lpmatrix(direct))
      W <- rkhs_extract_working_model(direct)
      expect_equal(fit$working_error[, i], unname(W$working_error), tolerance = 1e-10)
      expect_equal(fit$working_variance[, i], unname(W$working_variance), tolerance = 1e-10)
      cols <- fit$geometry$nuisance_columns
      expect_equal(fit$nuisance_covariance[[i]], direct$Vp[cols, cols, drop = FALSE], tolerance = 1e-10)
      p <- mgcvST:::taps_score_test(direct, test.component = 2L, method = "liu", lpmatrix = model$L)
      expect_equal(unname(fit$diagnostics$marginal_p_value[i]), p$smooth.pvalue, tolerance = 1e-10)
    }
    retained <- mgcvST.marginal(fit, calibration = "liu", BPPARAM = BiocParallel::SerialParam())
    expect_equal(retained$p_value, fit$diagnostics$marginal_p_value, tolerance = 1e-10)
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
  expect_error(mgcvST.estimate(f$Y, model, worker_init = function() NULL), "fixes the shared design")
  expect_error(mgcvST.estimate(f$Y, model, offset = matrix(0, 2, 3)), "offset must")
  expect_error(mgcvST.estimate(f$Y, model, offset = rep(NA_real_, ncol(f$Y))), "offset must")
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
