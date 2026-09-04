test_that("inline TAPS reuses the formal matrix for both estimation routes", {
  for (pc in c(FALSE, TRUE)) {
    f <- st_fixture(pc = pc, nuisance = TRUE)
    for (G in list(f$G, f$model)) {
      old_options <- options(mgcvST.test_predict_count = 0L)
      trace("predict.gam", where = asNamespace("mgcv"), print = FALSE,
            tracer = quote(if (identical(type, "lpmatrix")) {
              options(mgcvST.test_predict_count =
                        getOption("mgcvST.test_predict_count") + 1L)
            }))
      fit <- mgcvST.estimate(G = G, Y = f$Y,
                             marginal_args = list(method = "liu"),
                             BPPARAM = BiocParallel::SerialParam())
      untrace("predict.gam", where = asNamespace("mgcv"))
      count <- getOption("mgcvST.test_predict_count")
      options(old_options)
      expect_identical(count, 1L)
      expect_true(all(is.finite(fit$diagnostics$marginal_p_value)))
      expect_null(fit$gam)
      expect_null(fit$marginal_data)
    }
  }
})

test_that("cached TAPS matches direct TAPS with two factors and a nuisance smooth", {
  f <- st_fixture()
  d <- f$data
  d$a <- factor(rep(1:3, length.out = nrow(d)))
  d$b <- factor(rep(1:2, each = 3, length.out = nrow(d)))
  s <- mgcv::s
  for (family in list(gaussian(), mgcv::nb())) {
    G <- mgcv::gam(response ~ a + b + offset(offset0) + s(z, k = 5) +
                    s(x, y, bs = "spde", xt = f$basis),
                  data = d, family = family, fit = FALSE)
    g <- mgcv::gam(G = G, method = "REML")
    L <- predict(g, type = "lpmatrix")
    model <- model.set(response ~ a + b + offset(offset0) + s(z, k = 5),
                       d, f$basis, family = family)
    for (method in c("liu", "davies")) {
      a <- mgcvST:::taps_score_test(g, test.component = 2L, method = method)
      b <- mgcvST:::taps_score_test(g, test.component = 2L, method = method,
                                     lpmatrix = L)
      expect_identical(a, b)
      out <- mgcvST.estimate(f$Y[1L, , drop = FALSE], model,
                             marginal_args = list(method = method),
                             BPPARAM = BiocParallel::SerialParam())
      expect_equal(out$diagnostics$marginal_p_value, a$smooth.pvalue,
                   tolerance = 1e-12)
    }
  }
})

test_that("inline marginal tests agree on Serial and SOCK workers", {
  skip_on_cran()
  f <- st_fixture(nuisance = TRUE)
  bp <- BiocParallel::SnowParam(2L, type = "SOCK", progressbar = FALSE)
  on.exit(BiocParallel::bpstop(bp), add = TRUE)
  BiocParallel::bpstart(bp)
  for (G in list(f$G, f$model)) {
    for (method in c("liu", "davies")) {
      a <- mgcvST.estimate(f$Y, G, marginal_args = list(method = method),
                           BPPARAM = BiocParallel::SerialParam())
      b <- mgcvST.estimate(f$Y, G, marginal_args = list(method = method),
                           BPPARAM = bp, chunk_size = 1L)
      expect_equal(a$diagnostics$marginal_p_value,
                   b$diagnostics$marginal_p_value, tolerance = 1e-12)
      expect_true(all(is.finite(b$diagnostics$marginal_p_value)))
      for (field in c("marginal_requested_method", "marginal_method", "marginal_fallback")) {
        expect_identical(a$diagnostics[[field]], b$diagnostics[[field]])
      }
    }
  }
})

test_that("both estimation routes retain mixed Davies and Liu outcomes", {
  f <- st_fixture(nuisance = TRUE)
  for (G in list(f$G, f$model)) {
    liu <- mgcvST.estimate(f$Y, G, marginal_args = list(method = "liu"),
                           BPPARAM = BiocParallel::SerialParam())
    local({
      calls <- 0L
      local_mocked_bindings(davies = function(...) {
        calls <<- calls + 1L
        if (calls == 1L) list(ifault = 0L, Qq = .222) else
          list(ifault = 1L, Qq = .333)
      }, .package = "CompQuadForm")
      fit <- mgcvST.estimate(f$Y, G, BPPARAM = BiocParallel::SerialParam())
      expect_identical(fit$diagnostics$marginal_requested_method, rep("davies", 3L))
      expect_identical(fit$diagnostics$marginal_method, c("davies", "liu", "liu"))
      expect_identical(fit$diagnostics$marginal_fallback, c(FALSE, TRUE, TRUE))
      expect_identical(fit$diagnostics$marginal_p_value[1L], .222)
      expect_equal(fit$diagnostics$marginal_p_value[-1L],
                   liu$diagnostics$marginal_p_value[-1L], tolerance = 1e-12)
    })
    expect_identical(liu$diagnostics$marginal_method, rep("liu", 3L))
    expect_false(any(liu$diagnostics$marginal_fallback))
    custom <- mgcvST.estimate(f$Y, G, marginal_test = function(...) list(smooth.pvalue = .5),
                              BPPARAM = BiocParallel::SerialParam())
    expect_true(all(is.na(custom$diagnostics$marginal_method)))
    expect_true(all(is.na(custom$diagnostics$marginal_fallback)))
  }
})
