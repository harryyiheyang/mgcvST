# MAGIC 3D overview: orthographic isometric projections on a two-dimensional page.
library(ggplot2)
library(patchwork)
library(svglite)
library(ragg)

src <- "artifacts/datasets/MAGIC/MAGIC.rds"
fit_dir <- "artifacts/magic-slide-exploration/fits"
audit_dir <- "artifacts/magic-slide-exploration/3d-overview"
fig_dir <- "man/figures"
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

magic <- readRDS(src)
X <- magic$covariates
mesh <- magic$meshes$native3d
fit <- readRDS(file.path(fit_dir, "Tfap2b-spatial.rds"))
cfg <- readRDS(file.path(fit_dir, "configuration.rds"))

if (nrow(X) != 97830L || nrow(magic$slices) != 93L) {
  stop("MAGIC.rds must contain the 97,830-observation, 93-section retained set.")
}
if (nrow(mesh$loc) != 1962L || nrow(mesh$tv) != 7676L ||
    !identical(mesh$contract$tetrahedra_index_base, 1L)) {
  stop("MAGIC.rds does not contain the expected 1,962-node, 7,676-tetrahedron mesh.")
}
priors <- c(cfg$spatial_prior$prior, cfg$nb_prior$prior,
            cfg$slide_prior$prec$prior, cfg$ou_prior$prec$prior,
            cfg$ou_prior$phi$prior)
if (!all(priors == "flat") ||
    !identical(fit$gene, "Tfap2b") || !identical(fit$model, "spatial")) {
  stop("The overview requires the saved flat-prior Tfap2b spatial baseline fit.")
}
ix <- match(X$point_id, fit$point_id)
if (anyNA(ix) || anyDuplicated(X$point_id) || length(fit$mu) != nrow(X)) {
  stop("The Tfap2b fitted means do not have a one-to-one match to MAGIC observations.")
}

# This orthographic isometric map preserves the physical length of each x, y, and z axis.
iso <- function(x, y, z) {
  data.frame(u = sqrt(3) / 2 * (x - y), v = z + 0.5 * (x + y))
}
P <- iso(X$x_mm, X$y_mm, X$z_mm)
N <- iso(mesh$loc[, 1L], mesh$loc[, 2L], mesh$loc[, 3L])
D <- data.frame(u = P$u, v = P$v, z = X$z_mm,
                depth = -X$x_mm - X$y_mm + X$z_mm,
                rate = fit$mu[ix] / (X$total_umi / 10000))
lim <- c(0, stats::quantile(D$rate, 0.98, names = FALSE))
if (!is.finite(lim[2L]) || lim[2L] <= 0) stop("Tfap2b fitted-rate P98 is not positive.")

# Extract actual tetrahedron boundary faces. A face occurs once only if it is not shared.
tv <- mesh$tv
F <- rbind(tv[, c(1L, 2L, 3L)], tv[, c(1L, 2L, 4L)],
           tv[, c(1L, 3L, 4L)], tv[, c(2L, 3L, 4L)])
F <- t(apply(F, 1L, sort))
key <- apply(F, 1L, paste, collapse = "-")
keep <- key %in% names(table(key))[table(key) == 1L]
B <- F[keep, , drop = FALSE]
centroid_x <- rowMeans(matrix(mesh$loc[B, 1L], ncol = 3L))
cut_x <- stats::median(mesh$loc[, 1L])
B <- B[centroid_x <= cut_x, , drop = FALSE]
B <- B[order(rowMeans(matrix(-mesh$loc[B, 1L] - mesh$loc[B, 2L] + mesh$loc[B, 3L],
                             ncol = 3L))), , drop = FALSE]
M <- data.frame(face = rep(seq_len(nrow(B)), each = 3L),
                u = as.vector(t(matrix(N$u[B], ncol = 3L))),
                v = as.vector(t(matrix(N$v[B], ncol = 3L))))
N$depth <- -mesh$loc[, 1L] - mesh$loc[, 2L] + mesh$loc[, 3L]

all_u <- range(c(D$u, N$u))
all_v <- range(c(D$v, N$v))
pad_u <- diff(all_u) * 0.03
pad_v <- diff(all_v) * 0.03
xlim <- all_u + c(-pad_u, pad_u)
ylim <- all_v + c(-pad_v, pad_v)
base_theme <- theme_classic(base_size = 7, base_family = "Arial") +
  theme(
    axis.line = element_line(linewidth = 0.3),
    axis.ticks = element_line(linewidth = 0.3),
    axis.text = element_text(size = 5.5),
    axis.title = element_text(size = 6),
    plot.title = element_text(size = 7.4, face = "bold"),
    plot.subtitle = element_text(size = 5.3, margin = margin(b = 2)),
    plot.caption = element_text(size = 5, hjust = 0),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.title = element_text(size = 5.4),
    legend.text = element_text(size = 5.1),
    plot.margin = margin(2, 2, 2, 2)
  )
theme_set(base_theme)

p1 <- ggplot(D[order(D$depth), ], aes(u, v)) +
  geom_point(colour = "#2B6CA3", alpha = 0.20, size = 0.10) +
  geom_segment(data = data.frame(axis = c("x", "y", "z"),
                                 u = xlim[1L] + 0.75, v = ylim[1L] + 0.55,
                                 uend = xlim[1L] + c(1.14, 0.36, 0.75),
                                 vend = ylim[1L] + c(0.775, 0.775, 1.00)),
               aes(x = u, y = v, xend = uend, yend = vend), inherit.aes = FALSE,
               colour = "#263B4E", linewidth = 0.35,
               arrow = grid::arrow(length = grid::unit(1.2, "mm"))) +
  geom_text(data = data.frame(axis = c("x", "y", "z"),
                              u = xlim[1L] + c(1.23, 0.27, 0.75),
                              v = ylim[1L] + c(0.82, 0.82, 1.11)),
            aes(u, v, label = axis), inherit.aes = FALSE,
            colour = "#263B4E", size = 2.3, fontface = "bold") +
  coord_equal(xlim = xlim, ylim = ylim, expand = FALSE) +
  labs(title = "a  Observations",
       subtitle = "97,830 points across 93 z layers",
       x = "u (mm)", y = "v (mm)")

p2 <- ggplot() +
  geom_polygon(data = M, aes(u, v, group = face), fill = NA,
               colour = "#627B96", linewidth = 0.10, alpha = 0.72) +
  geom_point(data = N[order(N$depth), ], aes(u, v), colour = "#165B8E", alpha = 0.42, size = 0.18) +
  coord_equal(xlim = xlim, ylim = ylim, expand = FALSE) +
  labs(title = "b  Tetrahedral mesh",
       subtitle = sprintf("1,962 nodes; 7,676 tetrahedra; %d faces", nrow(B)),
       x = "u (mm)", y = "v (mm)")

p3 <- ggplot(D[order(D$depth), ], aes(u, v, colour = rate)) +
  geom_point(alpha = 0.70, size = 0.10) +
  coord_equal(xlim = xlim, ylim = ylim, expand = FALSE) +
  scale_colour_gradientn(
    colours = c("#FFFFFF", "#6BAED6", "#08306B"), limits = lim,
    oob = scales::squish, name = "Tfap2b fitted count / 10,000 UMI",
    guide = guide_colourbar(direction = "horizontal", position = "bottom",
                            barwidth = grid::unit(24, "mm"),
                            barheight = grid::unit(2.4, "mm"),
                            title.theme = element_text(size = 5.4),
                            label.theme = element_text(size = 5.1))
  ) +
  labs(title = "c  Tfap2b prediction",
       subtitle = sprintf("Flat spatial baseline; P98 cap = %.2f", lim[2L]),
       x = "u (mm)", y = "v (mm)")

p <- (p1 | p2 | p3) +
  plot_annotation(
    title = "MAGIC: stacked sections and 3D prediction",
    theme = theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5))
  )

png <- file.path(fig_dir, "magic-3d-overview.png")
pdf <- file.path(audit_dir, "magic-3d-overview.pdf")
svg <- file.path(audit_dir, "magic-3d-overview.svg")
ggsave(png, p, width = 183, height = 88, units = "mm", dpi = 600, bg = "white")
grDevices::cairo_pdf(pdf, width = 183 / 25.4, height = 88 / 25.4, family = "Arial")
print(p)
dev.off()
svglite(svg, width = 183 / 25.4, height = 88 / 25.4)
print(p)
dev.off()

S <- data.frame(
  retained_observations = nrow(D), z_layers = nrow(magic$slices),
  mesh_nodes = nrow(mesh$loc), tetrahedra = nrow(mesh$tv),
  boundary_faces_total = nrow(F[keep, , drop = FALSE]),
  boundary_faces_displayed = nrow(B), cutaway_rule = sprintf("face centroid x <= %.6f mm", cut_x),
  tfap2b_rate_p98 = lim[2L], projection = "u=sqrt(3)/2*(x-y); v=z+0.5*(x+y)",
  stringsAsFactors = FALSE
)
write.csv(S, file.path(audit_dir, "magic-3d-overview-source-audit.csv"), row.names = FALSE)
writeLines(c(
  "# MAGIC 3D overview source audit",
  "",
  "- Input: `artifacts/datasets/MAGIC/MAGIC.rds`, including its native 3D mesh.",
  sprintf("- Retained input: %s observations over %s registered z sections.", nrow(D), nrow(magic$slices)),
  sprintf("- Mesh: %s nodes and %s tetrahedra; %s boundary faces are shown after the declared cutaway.",
          nrow(mesh$loc), nrow(mesh$tv), nrow(B)),
  sprintf("- Cutaway rule: retain a true boundary face only when its centroid x <= %.6f mm.", cut_x),
  "  Boundary faces are tetrahedron faces seen once; the script does not compute a new Delaunay mesh.",
  "- Panel c input: saved `Tfap2b-spatial.rds`, matched one-to-one by `point_id`.",
  "  The shared configuration was checked to use flat spatial, negative-binomial, slide, and OU priors.",
  sprintf("- Panel c rate: fitted count divided by total UMI / 10,000; linear 0--P98 scale, P98 = %.6f.", lim[2L]),
  "  Zero is pure white (`#FFFFFF`); values above P98 are capped for display.",
  "- Projection: `u = sqrt(3)/2*(x-y)` and `v = z + 0.5*(x+y)`, with `coord_equal()`.",
  "  The x/y/z direction triad in panel a follows those projected vectors.",
  "  Points and displayed faces are drawn far-to-near by depth `-x - y + z`, toward the right-hand normal `(-1, -1, 1)` of the projection plane.",
  "  It is a 3D-to-2D orthographic view, not a planar slice or a volume interpolation.",
  "- Reproducible rendering sources: `magic-3d-overview.pdf` and `magic-3d-overview.svg`."
), file.path(audit_dir, "magic-3d-overview-source-audit.md"))
