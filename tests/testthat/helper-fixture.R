st_fixture <- function(n = 90L, family = mgcv::nb(), pc = FALSE, nuisance = FALSE) {
  skip_if_not_installed("geometry")
  set.seed(81)
  xy <- as.matrix(expand.grid(x = seq(0, 1, length.out = 5), y = seq(0, 1, length.out = 5)))
  mesh <- list(loc = xy, graph = list(tv = geometry::delaunayn(xy)))
  data <- data.frame(x = runif(n, .01, .99), y = runif(n, .01, .99),
                     offset0 = runif(n, -.2, .2), z = runif(n))
  basis <- spde_basis(mesh, as.matrix(data[, c("x", "y")]), kappa = .7, pc_cutoff = .95)
  data$response <- rpois(n, exp(1 + data$x + sin(5 * data$y)))
  data$response2 <- rpois(n, exp(.8 + data$x))
  data$response3 <- rpois(n, exp(.7 + data$y))
  s <- mgcv::s
  f <- if (nuisance) response ~ offset(offset0) + s(z, k = 5) else response ~ offset(offset0)
  model <- model.set(f, data, basis, family = family)
  f2 <- if (pc) response ~ offset(offset0) + s(x, y, bs = "spdePC", xt = basis) else
    response ~ offset(offset0) + s(x, y, bs = "spde", xt = basis)
  G <- mgcv::gam(f2, data = data, family = family, fit = FALSE)
  Y <- t(as.matrix(data[, c("response", "response2", "response3")]))
  list(data = data, basis = basis, model = model, G = G, Y = Y)
}

old_st <- function() {
  repo <- Sys.getenv("MGCVST_BASELINE")
  skip_if(!dir.exists(file.path(repo, "R")), "Set MGCVST_BASELINE for pinned-source equivalence checks")
  e <- new.env(parent = asNamespace("mgcvST"))
  for (f in list.files(file.path(repo, "R"), full.names = TRUE, pattern = "\\.R$")) sys.source(f, e)
  e
}

strip_elapsed <- function(x) {
  x$call <- x$timing <- NULL
  x$nuisance_covariance <- NULL
  if (!is.null(x$geometry)) {
    x$geometry$nuisance_columns <- NULL
    x$geometry$nuisance_design <- NULL
    x$geometry$nuisance_projection <- NULL
  }
  if (!is.null(x$diagnostics)) {
    keep <- !grepl("_seconds$", names(x$diagnostics))
    x$diagnostics <- x$diagnostics[, keep, drop = FALSE]
  }
  x
}

expect_numerically_equivalent_test <- function(x, y, tolerance = 1e-10) {
  numeric <- c("signed_score", "statistic", "information", "effective_rank",
               "p_two_sided", "p_positive", "p_negative", "p_adjusted",
               "p_positive_adjusted", "p_negative_adjusted")
  expect_identical(names(x), names(y))
  expect_identical(names(x$results), names(y$results))
  for (name in numeric) {
    expect_equal(x$results[[name]], y$results[[name]], tolerance = tolerance)
    y$results[[name]] <- x$results[[name]]
  }
  expect_equal(x$threshold$raw_p_threshold, y$threshold$raw_p_threshold,
               tolerance = tolerance)
  y$threshold$raw_p_threshold <- x$threshold$raw_p_threshold
  expect_identical(strip_elapsed(x), strip_elapsed(y))
}
