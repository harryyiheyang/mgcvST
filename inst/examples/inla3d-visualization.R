# MAGIC Snap25 3D visualization from the saved full-data INLA fit.
library(ggplot2)
library(patchwork)
library(plotly)
library(htmlwidgets)

src <- "artifacts/inla3d-transfer/spde3d_transfer_2026-09-13/data/magic"
fit_file <- "artifacts/inla3d-transfer/generic-flat/fits/fit-1.rds"
out <- "artifacts/inla3d-visualization"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
dir.create("man/figures", recursive = TRUE, showWarnings = FALSE)

d <- read.delim(gzfile(file.path(src, "aligned_points.tsv.gz")), check.names = FALSE)
g <- read.delim(gzfile(file.path(src, "Snap25.tsv.gz")), check.names = FALSE)
ix <- match(d$point_id, g$point_id)
if (anyNA(ix) || anyDuplicated(d$point_id) || anyDuplicated(g$point_id)) {
  stop("MAGIC coordinates and Snap25 counts do not have a one-to-one point match.")
}
d$count <- g$count[ix]
if (!all(d$total_umi == g$total_umi[ix])) {
  stop("MAGIC library sizes disagree after point matching.")
}

fit <- readRDS(fit_file)
if (length(fit$eta) != nrow(d) || length(fit$y) != nrow(d) ||
    !all(fit$y == d$count)) {
  stop("The saved full-data fit does not match the MAGIC Snap25 observations.")
}

spec <- jsonlite::read_json(file.path(src, "mesh_contract.json"), simplifyVector = TRUE)
xyz <- as.matrix(d[, c("x_aligned", "y_aligned", "z")])
xyz <- sweep(sweep(xyz, 2L, spec$source_origin, "-"),
             2L, spec$source_units_per_mesh_unit, "/")
colnames(xyz) <- c("x", "y", "z")
E <- d$total_umi / 10000
beta <- as.numeric(fit$fixed["b0", "mean"])

D <- data.frame(
  point_id = d$point_id,
  slice_id = d$slice_id,
  x = xyz[, 1L], y = xyz[, 2L], z = xyz[, 3L],
  count = d$count,
  exposure = E,
  observed_rate = d$count / E,
  fitted_eta = as.numeric(fit$eta),
  spatial_field = as.numeric(fit$eta) - beta,
  fitted_rate = exp(as.numeric(fit$eta)),
  spatial_fold = exp(as.numeric(fit$eta) - beta),
  fitted_mean_count = E * exp(as.numeric(fit$eta)),
  stringsAsFactors = FALSE
)

rx <- diff(range(D$x))
ry <- diff(range(D$y))
rz <- diff(range(D$z))
cx <- stats::median(D$x)
cy <- stats::median(D$y)
uz <- sort(unique(D$z))
cz <- uz[which.min(abs(uz - stats::median(D$z)))]
sx <- abs(D$x - cx) <= 0.04 * rx
sy <- abs(D$y - cy) <= 0.04 * ry
sz <- D$z == cz

lim_count <- c(0, stats::quantile(D$count, 0.98, names = FALSE))
lim_rate <- c(0, max(stats::quantile(D$observed_rate, 0.98, names = FALSE),
                     stats::quantile(D$fitted_rate, 0.98, names = FALSE)))
lim_fold <- c(0, stats::quantile(D$spatial_fold, 0.98, names = FALSE))
blue_scale <- c("#FFFFFF", "#6BAED6", "#08306B")
plotly_blue_scale <- list(c(0, "#FFFFFF"), c(0.5, "#6BAED6"), c(1, "#08306B"))
scale_blue <- function(limits, name) {
  scale_colour_gradientn(
    colours = blue_scale, limits = limits, oob = scales::squish, name = name,
    guide = guide_colourbar(
      direction = "horizontal", position = "bottom",
      barwidth = grid::unit(35, "mm"), barheight = grid::unit(2.5, "mm"),
      title.theme = element_text(size = 5.2),
      label.theme = element_text(size = 5)
    )
  )
}
plotly_colorbar <- function(title) {
  list(
    title = list(text = title, side = "top", font = list(size = 10)),
    orientation = "h", x = 0.5, xanchor = "center", y = -0.22,
    yanchor = "top", thickness = 10, len = 0.5, tickfont = list(size = 9)
  )
}
D$axon_x <- D$x + 0.42 * D$z
D$axon_y <- D$y + 0.24 * D$z
z_show <- uz[unique(round(seq(1, length(uz), length.out = 6L)))]
Dz <- D[D$z %in% z_show, ]
Dz$z_panel <- factor(sprintf("z = %.2f", Dz$z),
                     levels = sprintf("z = %.2f", z_show))

theme_set(
  theme_classic(base_size = 7, base_family = "Arial") +
    theme(axis.line = element_line(linewidth = 0.3),
          axis.ticks = element_line(linewidth = 0.3),
          legend.title = element_text(size = 5.2),
          legend.text = element_text(size = 5),
          legend.position = "bottom",
          legend.box = "horizontal",
          plot.title = element_text(size = 6.2, face = "bold"),
          plot.subtitle = element_text(size = 5.2),
          plot.tag = element_text(size = 9, face = "bold"))
)

p1 <- ggplot(D, aes(axon_x, axon_y, colour = count)) +
  geom_point(size = 0.12, alpha = 0.68) +
  scale_blue(lim_count, "Raw count") +
  coord_equal() +
  labs(title = "Observed Snap25 count",
       subtitle = sprintf("All 97,830 observations; upper colors capped at count %.0f (98th percentile)",
                          lim_count[2L]),
       x = "Axonometric x", y = "Axonometric y")

p2 <- ggplot(D, aes(axon_x, axon_y, colour = observed_rate)) +
  geom_point(size = 0.12, alpha = 0.68) +
  scale_blue(lim_rate, "Count per\n10,000 UMI") +
  coord_equal() +
  labs(title = "Observed normalized expression",
       subtitle = "count/(total UMI/10,000); shared rate scale",
       x = "Axonometric x", y = "Axonometric y")

p3 <- ggplot(D, aes(axon_x, axon_y, colour = fitted_rate)) +
  geom_point(size = 0.12, alpha = 0.68) +
  scale_blue(lim_rate, "Fitted count per\n10,000 UMI") +
  coord_equal() +
  labs(title = "Fitted expression rate",
       subtitle = "exp(eta); shared rate scale",
       x = "Axonometric x", y = "Axonometric y")

p4 <- ggplot(D[sx, ], aes(y, z, colour = spatial_fold)) +
  geom_point(size = 0.3, alpha = 0.75) +
  scale_blue(lim_fold, "Spatial\nfold-change") +
  coord_equal() +
  labs(title = "Central X slab", subtitle = sprintf("|x - %.2f| <= %.2f; n = %s", cx, 0.04 * rx,
                                                    format(sum(sx), big.mark = ",")),
       x = "y", y = "z")

p5 <- ggplot(D[sy, ], aes(x, z, colour = spatial_fold)) +
  geom_point(size = 0.3, alpha = 0.75) +
  scale_blue(lim_fold, "Spatial\nfold-change") +
  coord_equal() +
  labs(title = "Central Y slab", subtitle = sprintf("|y - %.2f| <= %.2f; n = %s", cy, 0.04 * ry,
                                                   format(sum(sy), big.mark = ",")),
       x = "x", y = "z")

p6 <- ggplot(D[sz, ], aes(x, y, colour = spatial_fold)) +
  geom_point(size = 0.3, alpha = 0.75) +
  scale_blue(lim_fold, "Spatial\nfold-change") +
  coord_equal() +
  labs(title = "Central Z plane", subtitle = sprintf("z = %.2f; n = %s", cz,
                                                    format(sum(sz), big.mark = ",")),
       caption = "Spatial fold-change 1 is the intercept-only reference.",
       x = "x", y = "y")

p7 <- ggplot(Dz, aes(x, y, colour = observed_rate)) +
  geom_point(size = 0.16, alpha = 0.65) +
  facet_wrap(~z_panel, nrow = 1) + coord_equal() +
  scale_blue(lim_rate, "Count per\n10,000 UMI") +
  labs(title = "Observed normalized expression across selected Z planes",
       subtitle = "Six sections spanning z; every point in each plane",
       x = "x", y = "y")

p8 <- ggplot(Dz, aes(x, y, colour = fitted_rate)) +
  geom_point(size = 0.16, alpha = 0.65) +
  facet_wrap(~z_panel, nrow = 1) + coord_equal() +
  scale_blue(lim_rate, "Fitted count per\n10,000 UMI") +
  labs(title = "Fitted expression rate across the same Z planes",
       subtitle = "Same units and color scale as the observed row",
       x = "x", y = "y")

fig <- ((p1 | p2 | p3) / (p4 | p5 | p6) / p7 / p8 +
  plot_layout(heights = c(1.05, 1, 0.72, 0.72), guides = "collect") +
  plot_annotation(tag_levels = "a")) &
  theme(legend.position = "bottom", legend.box = "horizontal")

ggsave(file.path(out, "magic-snap25-inla3d.png"), fig,
       width = 183, height = 255, units = "mm", dpi = 600, bg = "white")
svglite::svglite(file.path(out, "magic-snap25-inla3d.svg"),
                 width = 183 / 25.4, height = 255 / 25.4)
print(fig)
dev.off()
grDevices::cairo_pdf(file.path(out, "magic-snap25-inla3d.pdf"),
                     width = 183 / 25.4, height = 255 / 25.4, family = "Arial")
print(fig)
dev.off()
ragg::agg_tiff(file.path(out, "magic-snap25-inla3d.tiff"),
               width = 183 / 25.4, height = 255 / 25.4,
               units = "in", res = 600, background = "white")
print(fig)
dev.off()

hover <- paste0(
  D$point_id, "<br>x=", signif(D$x, 4), ", y=", signif(D$y, 4),
  ", z=", signif(D$z, 4), "<br>count=", D$count,
  "<br>fitted mean=", signif(D$fitted_mean_count, 4),
  "<br>observed rate=", signif(D$observed_rate, 4),
  "<br>fitted rate=", signif(D$fitted_rate, 4),
  "<br>spatial fold-change=", signif(D$spatial_fold, 4)
)
p3d <- plot_ly()
p3d <- add_trace(p3d, data = D, x = ~x, y = ~y, z = ~z, text = hover,
                 hoverinfo = "text", type = "scatter3d", mode = "markers",
                 name = "Raw count", visible = TRUE,
                 marker = list(size = 1.25, opacity = 0.72, color = D$count,
                               colorscale = plotly_blue_scale, cmin = lim_count[1L], cmax = lim_count[2L],
                               colorbar = plotly_colorbar("Raw count")))
p3d <- add_trace(p3d, data = D, x = ~x, y = ~y, z = ~z, text = hover,
                 hoverinfo = "text", type = "scatter3d", mode = "markers",
                 name = "Observed rate", visible = FALSE,
                 marker = list(size = 1.25, opacity = 0.72, color = D$observed_rate,
                               colorscale = plotly_blue_scale, cmin = lim_rate[1L], cmax = lim_rate[2L],
                               colorbar = plotly_colorbar("Count per<br>10,000 UMI")))
p3d <- add_trace(p3d, data = D, x = ~x, y = ~y, z = ~z, text = hover,
                 hoverinfo = "text", type = "scatter3d", mode = "markers",
                 name = "Fitted rate", visible = FALSE,
                 marker = list(size = 1.25, opacity = 0.72, color = D$fitted_rate,
                               colorscale = plotly_blue_scale, cmin = lim_rate[1L], cmax = lim_rate[2L],
                               colorbar = plotly_colorbar("Fitted count per<br>10,000 UMI")))
p3d <- add_trace(p3d, data = D, x = ~x, y = ~y, z = ~z, text = hover,
                 hoverinfo = "text", type = "scatter3d", mode = "markers",
                 name = "Spatial fold-change", visible = FALSE,
                 marker = list(size = 1.25, opacity = 0.72, color = D$spatial_fold,
                               colorscale = plotly_blue_scale,
                               cmin = lim_fold[1L], cmax = lim_fold[2L],
                               colorbar = plotly_colorbar("Spatial fold-change<br>(1 = intercept reference)")))
p3d <- layout(
  p3d,
  title = "MAGIC Snap25: all 97,830 observations",
  showlegend = FALSE,
  margin = list(t = 90, b = 85),
  scene = list(xaxis = list(title = "x"), yaxis = list(title = "y"),
               zaxis = list(title = "z"), aspectmode = "data"),
  updatemenus = list(list(type = "dropdown", x = 0.02, y = 1.08,
    buttons = list(
      list(method = "update", args = list(list(visible = c(TRUE, FALSE, FALSE, FALSE)),
           list(title = "MAGIC Snap25: raw count (all observations)")), label = "Raw count"),
      list(method = "update", args = list(list(visible = c(FALSE, TRUE, FALSE, FALSE)),
           list(title = "MAGIC Snap25: observed count per 10,000 UMI")), label = "Observed rate"),
      list(method = "update", args = list(list(visible = c(FALSE, FALSE, TRUE, FALSE)),
           list(title = "MAGIC Snap25: fitted count per 10,000 UMI")), label = "Fitted rate"),
      list(method = "update", args = list(list(visible = c(FALSE, FALSE, FALSE, TRUE)),
           list(title = "MAGIC Snap25: spatial relative fold-change")), label = "Spatial fold-change")
    )
  ))
)
saveWidget(p3d, file.path(out, "magic-snap25-inla3d-interactive.html"),
           selfcontained = FALSE, title = "MAGIC Snap25 3D INLA fit")

S <- rbind(
  data.frame(view = "central_x_slab", center = cx, half_width = 0.04 * rx,
             observations = sum(sx)),
  data.frame(view = "central_y_slab", center = cy, half_width = 0.04 * ry,
             observations = sum(sy)),
  data.frame(view = "central_z_plane", center = cz, half_width = 0,
             observations = sum(sz))
)
write.csv(S, file.path(out, "slice-definition.csv"), row.names = FALSE)
write.csv(data.frame(
  observations = nrow(D), mesh_nodes = length(fit$u),
  nb_size = as.numeric(fit$metrics$nb_size), spatial_sigma = as.numeric(fit$metrics$sigma),
  intercept = beta, count_correlation = cor(D$count, D$fitted_mean_count),
  rate_correlation = cor(D$observed_rate, D$fitted_rate),
  count_rmse = sqrt(mean((D$count - D$fitted_mean_count)^2)),
  spatial_fold_scale_min = lim_fold[1L], spatial_fold_p98 = lim_fold[2L],
  constraint_error = as.numeric(fit$metrics$constraint_error)
), file.path(out, "fit-summary.csv"), row.names = FALSE)

Z <- aggregate(spatial_fold ~ z, D, function(x) {
  c(n = length(x), median = stats::median(x), q25 = stats::quantile(x, 0.25),
    q75 = stats::quantile(x, 0.75), p10 = stats::quantile(x, 0.10),
    p90 = stats::quantile(x, 0.90))
})
Z <- data.frame(z = Z$z, Z$spatial_fold, row.names = NULL)
write.csv(Z, file.path(out, "spatial-fold-by-z.csv"), row.names = FALSE)

readme_fig <- ((p2 + theme(legend.position = "bottom") |
  p3 + theme(legend.position = "bottom") |
  p6 + theme(legend.position = "bottom")) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a")) &
  theme(legend.position = "bottom", legend.box = "horizontal")
ggsave("man/figures/inla3d-snap25.png", readme_fig,
       width = 180, height = 62, units = "mm", dpi = 300, bg = "white")
saveRDS(D, file.path(out, "plotting-data.rds"), compress = "xz")
