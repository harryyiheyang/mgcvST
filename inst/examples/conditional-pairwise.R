# Small deterministic INLA example with three spatial features.
vertices <- as.matrix(expand.grid(
  x = seq(0, 1, length.out = 5L),
  y = seq(0, 1, length.out = 5L)
))
mesh <- list(loc = vertices, graph = list(tv = geometry::delaunayn(vertices)))
data <- expand.grid(
  x = seq(0.04, 0.96, length.out = 8L),
  y = seq(0.04, 0.96, length.out = 8L)
)
data$z <- seq(-1, 1, length.out = nrow(data))
altitude <- datasets::volcano[cbind(
  1L + round(86 * data$x), 1L + round(60 * data$y)
)]
Y <- rbind(
  gene_a = altitude / 100 + sin(17 * data$x + 11 * data$y) / 10,
  gene_b = altitude / 105 + cos(13 * data$x - 9 * data$y) / 10,
  gene_c = altitude / 95 + sin(7 * data$x - 15 * data$y) / 10
)
basis <- spde_basis(
  mesh, as.matrix(data[c("x", "y")]), kappa = 1.2,
  project_intercept = TRUE
)
model <- inlaST.set(response ~ z, data, basis, family = gaussian())
fit <- inlaST.estimate(
  Y, model, BPPARAM = BiocParallel::SerialParam(),
  control = list(fixed_precision = 1.7, gaussian_precision = 1 / 0.09)
)
result <- mgcvST::inlaST.test(
  fit, pairwise_method = "conditional_cauchy", method = "BY", threads = 2L,
  checkpoint_dir = file.path(tempdir(), "mgcvst-conditional-example")
)
print(result$results)
