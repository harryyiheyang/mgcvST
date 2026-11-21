.spdest_sparse_fixture <- function(n = 72L, seed = 1701L) {
  skip_if_not_installed("geometry")
  set.seed(seed)
  vertices <- as.matrix(expand.grid(
    x = seq(0, 1, length.out = 5L),
    y = seq(0, 1, length.out = 5L)
  ))
  mesh <- list(loc = vertices, graph = list(tv = geometry::delaunayn(vertices)))
  data <- data.frame(
    x = runif(n, 0.02, 0.98), y = runif(n, 0.02, 0.98),
    z = seq(-1, 1, length.out = n), exposure = runif(n, 0.8, 1.3)
  )
  data$offset0 <- log(data$exposure)
  basis <- spde_basis(
    mesh, as.matrix(data[c("x", "y")]), kappa = 1.2,
    project_intercept = TRUE
  )
  eta <- 0.4 + 0.25 * data$z + data$offset0 + 0.2 * sin(2 * pi * data$x)
  Y <- rbind(
    a = rpois(n, exp(1 + eta)),
    b = rpois(n, exp(1 + eta + 0.1 * cos(2 * pi * data$y))),
    c = rpois(n, exp(1 + eta - 0.1 * sin(2 * pi * data$y)))
  )
  list(data = data, basis = basis, Y = Y, mesh = mesh)
}

test_that("the pair test takes no BiocParallel backend and agrees across OpenMP counts", {
  skip_on_cran()
  f <- .spdest_sparse_fixture(n = 64L, seed = 1711L)
  model <- spdeST.set(
    response ~ z + offset(offset0), f$data, f$basis, family = poisson()
  )
  fit <- spdeST:::.spdest_estimate(
    .spdest_scan(
    f$Y, model,
    BPPARAM = BiocParallel::SerialParam()),
    f$Y, features = "all", BPPARAM = BiocParallel::SerialParam())
  expect_true(all(is.finite(fit$diagnostics$marginal_p_value)))
  pairs <- t(combn(rownames(f$Y), 2L))
  serial <- spdeST.test(fit, pairs = pairs, rank = 1L, chunk_size = 1L, moments = "exact")
  threaded <- spdeST.test(fit, pairs = pairs, rank = 1L, threads = 2L,
                          chunk_size = 1L, moments = "exact")
  expect_identical(nrow(serial$results), nrow(pairs))
  expect_equal(threaded$results, serial$results, tolerance = 1e-10)
  bp <- BiocParallel::SnowParam(2L, type = "SOCK", progressbar = FALSE)
  expect_error(spdeST.test(fit, pairs = pairs, BPPARAM = bp, moments = "exact"), "unused argument")
})

test_that("a nuisance smooth is rejected at set() by both setups", {
  skip_on_cran()
  f <- .spdest_sparse_fixture(n = 56L, seed = 1721L)
  s <- mgcv::s
  # The Stage 1 null is a GLM on both paths.
  expected <- "The null of spdeST.set\\(\\) is a GLM"
  expect_error(
    spdeST.set(
      response ~ s(z, k = 5) + offset(offset0),
      f$data, f$basis, family = poisson()
    ),
    expected
  )
  expect_error(
    spdeST.set(
      response ~ s(z, k = 5) + offset(offset0), data = f$data,
      family = poisson(), mesh = f$mesh, kappa = 1.2,
      coordinates = c("x", "y")
    ),
    expected
  )
})

test_that("parametric covariates remain fully supported by both setups", {
  skip_on_cran()
  f <- .spdest_sparse_fixture(n = 56L, seed = 1721L)
  f$data$lu <- as.numeric(scale(f$data$exposure))
  legacy <- spdeST.set(
    response ~ lu + offset(offset0), f$data, f$basis, family = poisson()
  )
  expect_length(legacy$spde_spec$random, 1L)
  expect_true(spdeST:::.spdest_sparse_score_capability(legacy)$eligible)
  expect_identical(legacy$spde_spec$fixed$names, c("(Intercept)", "lu"))
  expect_equal(unname(as.matrix(legacy$spde_spec$fixed$X)),
               unname(as.matrix(legacy$spde_spec$nuisance_design)))
  fit <- spdeST:::.spdest_estimate(
    .spdest_scan(
    f$Y, legacy, diagnostics = TRUE, BPPARAM = BiocParallel::SerialParam()),
    f$Y, features = "all", BPPARAM = BiocParallel::SerialParam())
  expect_true(all(fit$diagnostics$converged))
  expect_identical(ncol(fit$geometry$nuisance_design), 2L)

  # A factor covariate is a parametric term too.
  f$data$grp <- factor(rep(c("a", "b", "c"), length.out = nrow(f$data)))
  factored <- spdeST.set(
    response ~ grp + offset(offset0), f$data, f$basis, family = poisson()
  )
  expect_identical(ncol(factored$spde_spec$fixed$X), 3L)
  expect_length(factored$spde_spec$random, 1L)
})

test_that("a second spatial SPDE term is still rejected at set()", {
  skip_on_cran()
  f <- .spdest_sparse_fixture(n = 56L, seed = 1723L)
  s <- mgcv::s
  basis <- f$basis
  # Complete-formula route: two spatial SPDE terms would need two random
  # blocks, which the sparse kernel cannot carry.
  expect_error(
    spdeST.set(
      response ~ s(x, y, bs = "spde", xt = basis) +
        s(x, y, bs = "spde", xt = basis) + offset(offset0),
      data = f$data, family = poisson()
    ),
    "exactly one spatial SPDE term"
  )
})

test_that("a wide parametric nuisance design is rejected with a dedicated message", {
  expect_error(spdeST:::.spdest_check_nuisance_width(1001L),
               "dedicated implementation")
  expect_silent(spdeST:::.spdest_check_nuisance_width(1000L))
  skip_on_cran()
  f <- .spdest_sparse_fixture(n = 1100L, seed = 1725L)
  # A factor with many levels expands to a wide dense nuisance design.
  f$data$grp <- factor(rep(seq_len(1010L), length.out = nrow(f$data)))
  expect_error(
    spdeST.set(
      response ~ grp + offset(offset0), f$data, f$basis, family = poisson()
    ),
    "dedicated implementation"
  )
})

