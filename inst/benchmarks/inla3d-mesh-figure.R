# Actual z=0.5 sections of all four volume meshes; no synthetic display edges.
library(ggplot2)
library(patchwork)
options(warn = 2)

out <- "artifacts/inla3d/mesh"
fig <- "inst/validation/inla3d"
dir.create(fig, recursive = TRUE, showWarnings = FALSE)
M <- read.csv(file.path(out, "manifest.csv"))
M <- M[c(1, 3, 2, 4), ]
E <- t(combn(1:4, 2))
S <- list()
P <- list()
theta <- seq(0, 2 * pi, length.out = 201L)
circle <- data.frame(x = 1 + 0.35 * cos(theta), y = 1 + 0.35 * sin(theta))
for (m in seq_len(nrow(M))) {
  xyz <- as.matrix(read.csv(file.path(out, paste0(M$mesh[m], "-vertices.csv"))))
  tv <- as.matrix(read.csv(file.path(out, paste0(M$mesh[m], "-tetrahedra.csv"))))
  R1 <- list()
  j <- 0L
  for (tt in seq_len(nrow(tv))) {
    v <- xyz[tv[tt, ], , drop = FALSE]
    z <- v[, 3L] - 0.5
    if (min(z) > 0 || max(z) < 0) next
    p <- v[z == 0, 1:2, drop = FALSE]
    for (edge in seq_len(nrow(E))) {
      a <- E[edge, 1L]
      b <- E[edge, 2L]
      if (z[a] * z[b] < 0) {
        p <- rbind(p, v[a, 1:2] + (-z[a] / (z[b] - z[a])) * (v[b, 1:2] - v[a, 1:2]))
      }
    }
    p <- unique(round(p, 12L))
    if (nrow(p) < 3L) next
    mid <- colMeans(p)
    idx <- order(atan2(p[, 2L] - mid[2L], p[, 1L] - mid[1L]))
    p <- p[c(idx, idx[1L]), , drop = FALSE]
    j <- j + 1L
    R1[[j]] <- data.frame(x = p[, 1L], y = p[, 2L], polygon = tt, mesh = M$mesh[m])
  }
  S[[m]] <- do.call(rbind, R1)
  title <- paste0(if (M$kind[m] == "uniform") "Uniform" else "Adaptive", ", ", M$nodes[m], " nodes")
  P[[m]] <- ggplot(S[[m]], aes(x, y, group = polygon)) +
    geom_path(linewidth = 0.13, colour = "#5A6874") +
    geom_path(data = circle, aes(x, y), inherit.aes = FALSE,
              linewidth = 0.45, colour = "#B6573E", linetype = "dashed") +
    coord_equal(xlim = c(0, 3), ylim = c(0, 2), expand = FALSE) +
    labs(title = title, x = "x (mm)", y = "y (mm)") +
    theme_classic(base_size = 8, base_family = "Arial") +
    theme(plot.title = element_text(size = 9), axis.text = element_text(size = 7),
          axis.ticks = element_line(linewidth = 0.3), axis.line = element_line(linewidth = 0.3))
}
p <- wrap_plots(P, ncol = 2) + plot_annotation(
  title = "Geometry-driven refinement below 3,000 nodes",
  subtitle = "True tetrahedral sections at z = 0.5 mm; dashed circle: prespecified region",
  tag_levels = "a", theme = theme(text = element_text(family = "Arial", size = 9),
    plot.title = element_text(size = 11), plot.subtitle = element_text(size = 8)))
write.csv(do.call(rbind, S), file.path(fig, "mesh-section-source.csv"), row.names = FALSE)
grDevices::cairo_pdf(file.path(fig, "mesh-sections.pdf"),
                    width = 183 / 25.4, height = 145 / 25.4, family = "Arial")
print(p)
dev.off()
svglite::svglite(file.path(fig, "mesh-sections.svg"),
                 width = 183 / 25.4, height = 145 / 25.4)
print(p)
dev.off()
ggsave(file.path(fig, "mesh-sections.png"), p, device = ragg::agg_png,
       width = 183, height = 145, units = "mm", dpi = 300)
