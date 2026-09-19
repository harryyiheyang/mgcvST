# Small reproducible demonstration of the current public native 3D API.
library(mgcvST)
library(BiocParallel)
library(ggplot2)
library(patchwork)

out <- "artifacts/inla3d-visualization"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
set.seed(20260919)

xyz <- as.matrix(expand.grid(x = seq(0, 1, length.out = 4L),
                             y = seq(0, 1, length.out = 4L),
                             z = seq(0, 1, length.out = 4L)))
tv <- geometry::delaunayn(xyz)
mesh <- fmesher::fm_mesh_3d(loc = xyz, tv = tv)
n <- 600L
d <- data.frame(x = runif(n, 0.02, 0.98),
                y = runif(n, 0.02, 0.98),
                z = runif(n, 0.02, 0.98))
H <- 0.45 * sin(2 * pi * d$x) * cos(2 * pi * d$y) +
  0.25 * (d$z - 0.5)
H <- H - mean(H)
d$response <- 0
Y <- matrix(0.7 + H + rnorm(n, sd = 0.2), nrow = 1L,
            dimnames = list("simulated_feature", NULL))

m <- inlaST.set(response ~ 1, data = d, family = gaussian(), mesh = mesh,
                kappa = 3, coordinates = c("x", "y", "z"),
                control = list(fixed_precision = 4, gaussian_precision = 25))
fit <- inlaST.estimate(Y, m, BPPARAM = SerialParam(), retain_smooth = TRUE,
                      diagnostics = TRUE, control = list(poisson_screen_phi = 0))
if (!isTRUE(fit$diagnostics$converged[1L])) {
  stop("The public native 3D demonstration did not converge.")
}

A <- m$inla_spec$random[[1L]]$A
u <- as.numeric(fit$smooth_coefficients$global[1L, ])
Hhat <- as.numeric(A %*% u)
R <- data.frame(d[, c("x", "y", "z")], truth = H, fitted_field = Hhat,
                observed = as.numeric(Y))
R$axon_x <- R$x + 0.42 * R$z
R$axon_y <- R$y + 0.24 * R$z
lim <- max(abs(c(R$truth, R$fitted_field)))

theme_set(theme_classic(base_size = 7, base_family = "Arial") +
            theme(plot.title = element_text(size = 7.5, face = "bold"),
                  plot.subtitle = element_text(size = 6.2),
                  plot.tag = element_text(size = 9, face = "bold")))
p1 <- ggplot(R, aes(axon_x, axon_y, colour = truth)) +
  geom_point(size = 0.8, alpha = 0.8) + coord_equal() +
  scale_colour_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
                         limits = c(-lim, lim), midpoint = 0, name = "True field") +
  labs(title = "Simulated spatial field", subtitle = "All 600 observations",
       x = "Axonometric x", y = "Axonometric y")
p2 <- ggplot(R, aes(axon_x, axon_y, colour = fitted_field)) +
  geom_point(size = 0.8, alpha = 0.8) + coord_equal() +
  scale_colour_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
                         limits = c(-lim, lim), midpoint = 0, name = "Fitted field") +
  labs(title = "Current public native fit",
       subtitle = "Observation-mean constrained 3D SPDE",
       x = "Axonometric x", y = "Axonometric y")
p3 <- ggplot(R, aes(truth, fitted_field)) +
  geom_point(size = 0.8, alpha = 0.55, colour = "#3182BD") +
  geom_abline(slope = 1, intercept = 0, colour = "#B2182B", linewidth = 0.4) +
  coord_equal(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
  labs(title = "Field recovery",
       subtitle = sprintf("Pearson r = %.3f; RMSE = %.3f",
                          cor(R$truth, R$fitted_field),
                          sqrt(mean((R$truth - R$fitted_field)^2))),
       x = "True field", y = "Fitted field")
fig <- p1 | p2 | p3 + plot_annotation(tag_levels = "a")

ggsave(file.path(out, "public-native-3d-demo.png"), fig,
       width = 183, height = 62, units = "mm", dpi = 600, bg = "white")
grDevices::cairo_pdf(file.path(out, "public-native-3d-demo.pdf"),
                     width = 183 / 25.4, height = 62 / 25.4, family = "Arial")
print(fig)
dev.off()
saveRDS(list(data = R, mesh = mesh, model = m, fit = fit,
             truth = list(intercept = 0.7, field = H, error_sd = 0.2)),
        file.path(out, "public-native-3d-demo.rds"), compress = "xz")
write.csv(data.frame(
  observations = n, mesh_nodes = nrow(xyz), tetrahedra = nrow(tv),
  field_correlation = cor(R$truth, R$fitted_field),
  field_rmse = sqrt(mean((R$truth - R$fitted_field)^2)),
  fitted_observation_spatial_mean = mean(Hhat),
  constraint_residual = fit$constraint_residual[1L, 1L]
), file.path(out, "public-native-3d-demo-summary.csv"), row.names = FALSE)
