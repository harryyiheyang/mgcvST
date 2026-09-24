.cct_fit <- local({
  cached <- NULL
  function() {
    skip_if_not_installed("INLA")
    skip_if_not_installed("geometry")
    if (!is.null(cached)) return(cached)
    vertices <- as.matrix(expand.grid(x = seq(0, 1, length.out = 5L),
                                      y = seq(0, 1, length.out = 5L)))
    mesh <- list(loc = vertices, graph = list(tv = geometry::delaunayn(vertices)))
    data <- expand.grid(x = seq(0.04, 0.96, length.out = 8L),
                        y = seq(0.04, 0.96, length.out = 8L))
    data$z <- seq(-1, 1, length.out = nrow(data))
    altitude <- datasets::volcano[cbind(1L + round(86 * data$x),
                                        1L + round(60 * data$y))]
    Y <- rbind(
      gene_a = altitude / 100 + sin(17 * data$x + 11 * data$y) / 10,
      gene_b = altitude / 105 + cos(13 * data$x - 9 * data$y) / 10,
      gene_c = altitude / 95 + sin(7 * data$x - 15 * data$y) / 10,
      gene_d = altitude / 90 + cos(5 * data$x + 3 * data$y) / 10
    )
    basis <- spde_basis(mesh, as.matrix(data[c("x", "y")]), kappa = 1.2,
                        project_intercept = TRUE)
    model <- inlaST.set(response ~ z, data, basis, family = gaussian())
    cached <<- inlaST.estimate(
      Y, model, BPPARAM = BiocParallel::SerialParam(),
      control = list(fixed_precision = 1.7, gaussian_precision = 1 / 0.09)
    )
    cached
  }
})

test_that("conditional Cauchy runs end to end and returns raw results", {
  fit <- .cct_fit()
  out <- inlaST.test(fit, pairwise_method = "conditional_cauchy", FDR = FALSE)
  r <- out$results
  expect_s3_class(out, "mgcvST_test")
  expect_identical(names(r), c(
    "feature1", "feature2", "signed_score", "statistic", "p_two_sided",
    "log_p_two_sided", "p_1_given_2", "log_p_1_given_2", "p_2_given_1",
    "log_p_2_given_1", "p_adjusted", "log_p_adjusted", "discovered"))
  expect_identical(nrow(r), 6L)
  expect_identical(out$discoveries$pairs_tested, 6L)
  expect_identical(out$discoveries$pairs_discovered, sum(r$discovered))
  expect_identical(out$threshold$adjustment_method, "none")
  expect_identical(out$calibration, "conditional_cauchy")
  expect_equal(r$statistic, r$signed_score^2)
  expect_equal(r$p_two_sided, exp(r$log_p_two_sided))
  expect_equal(r$p_1_given_2, exp(r$log_p_1_given_2))
  expect_equal(r$p_2_given_1, exp(r$log_p_2_given_1))
  expect_true(all(r$p_two_sided >= pmin(r$p_1_given_2, r$p_2_given_1)))
  cauchy <- (tan((0.5 - r$p_1_given_2) * pi) + tan((0.5 - r$p_2_given_1) * pi)) / 2
  expect_equal(r$p_two_sided, 0.5 - atan(cauchy) / pi, tolerance = 1e-10)
  expect_identical(r$p_adjusted, r$p_two_sided)
  expect_identical(r$discovered, r$p_two_sided <= 0.05)
})

test_that("conditional Cauchy multiple testing is the caller's choice", {
  fit <- .cct_fit()
  pairs <- rbind(c("gene_a", "gene_b"), c("gene_a", "gene_c"), c("gene_b", "gene_d"))
  raw <- inlaST.test(fit, pairwise_method = "conditional_cauchy", pairs = pairs,
                     FDR = FALSE)
  by <- inlaST.test(fit, pairwise_method = "conditional_cauchy", pairs = pairs,
                    method = "BY", q.value = 0.2)
  bh <- inlaST.test(fit, pairwise_method = "conditional_cauchy", pairs = pairs)
  core <- c("feature1", "feature2", "signed_score", "statistic", "p_two_sided",
            "log_p_two_sided", "p_1_given_2", "p_2_given_1")
  expect_identical(by$results[core], raw$results[core])
  expect_identical(bh$results[core], raw$results[core])
  expect_identical(raw$results$feature1, pairs[, 1L])
  lp <- raw$results$log_p_two_sided
  expect_equal(by$results$log_p_adjusted, mgcvST:::.mgcvst_log_by(lp))
  expect_equal(by$results$p_adjusted, p.adjust(exp(lp), "BY"))
  expect_identical(by$results$discovered, by$results$log_p_adjusted <= log(0.2))
  expect_equal(bh$results$p_adjusted, p.adjust(exp(lp), "BH"))
  expect_identical(by$threshold$adjustment_method, "BY")
  expect_identical(bh$threshold$adjustment_method, "BH")
})

test_that("conditional Cauchy checkpoints resume to identical results", {
  fit <- .cct_fit()
  path <- withr::local_tempdir()
  first <- inlaST.test(fit, pairwise_method = "conditional_cauchy",
                       checkpoint_dir = path)
  again <- inlaST.test(fit, pairwise_method = "conditional_cauchy",
                       checkpoint_dir = path)
  expect_identical(again$results, first$results)
  expect_identical(again$timing$variance_rows_reused, 4L)
  expect_error(inlaST.test(fit, pairwise_method = "conditional_cauchy",
                           checkpoint_dir = path, resume = FALSE),
               "already exists")
})

test_that("inlaST.test validates the split argument set", {
  fit <- .cct_fit()
  pair <- matrix(c("gene_a", "gene_b"), 1L)
  expect_error(inlaST.test(fit, pairwise_method = "conditional_cauchy",
                           highlight = pair), "highlight")
  expect_error(inlaST.test(fit, pairwise_method = "conditional_cauchy",
                           liu_approximation = "pca_learning"),
               "liu_approximation")
  expect_error(inlaST.test(fit, pairs = pair, conditional_precision = "float32"),
               "conditional_precision")
  # mgcvST.test() does not accept inlaST.estimate() fits (Task E1); every
  # call below now fails at that guard, before its own argument checks.
  expect_error(mgcvST.test(fit, pairs = pair, pairwise_method = "score_liu"),
               "mgcvST.test\\(\\) does not accept inlaST.estimate\\(\\) fits; use inlaST.test\\(\\).")
  expect_error(inlaST.test(fit, pairs = pair, pairwise_method = "liu"),
               "score_liu")
  exact <- inlaST.test(fit, pairs = pair, pairwise_method = "score_liu",
                       liu_approximation = "exact")
  expect_error(mgcvST.test(fit, pairs = pair),
               "mgcvST.test\\(\\) does not accept inlaST.estimate\\(\\) fits; use inlaST.test\\(\\).")
})
