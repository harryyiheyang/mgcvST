# MAGIC Foxp1 and Tfap2b observation-only visualization on shared z planes.
library(ggplot2)
library(svglite)
library(ragg)

src <- "artifacts/inla3d-transfer/spde3d_transfer_2026-09-13/data/magic"
gene_dir <- "artifacts/datasets/MAGIC/additional-genes"
out <- "artifacts/magic-slide-exploration/markers"
dir.create(out, recursive = TRUE, showWarnings = FALSE)

d <- read.delim(gzfile(file.path(src, "aligned_points.tsv.gz")), check.names = FALSE)
spec <- jsonlite::read_json(file.path(src, "mesh_contract.json"), simplifyVector = TRUE)
xyz <- as.matrix(d[, c("x_aligned", "y_aligned", "z")])
xyz <- sweep(sweep(xyz, 2L, spec$source_origin, "-"),
             2L, spec$source_units_per_mesh_unit, "/")
colnames(xyz) <- c("x", "y", "z")
uz <- sort(unique(xyz[, "z"]))
z_show <- uz[unique(round(seq(1L, length(uz), length.out = 6L)))]

theme_set(
  theme_classic(base_size = 7, base_family = "Arial") +
    theme(
      axis.line = element_line(linewidth = 0.3),
      axis.ticks = element_line(linewidth = 0.3),
      legend.title = element_text(size = 5.2),
      legend.text = element_text(size = 5),
      plot.title = element_text(size = 7.5, face = "bold"),
      plot.subtitle = element_text(size = 6.2)
    )
)

blue_scale <- c("#FFFFFF", "#6BAED6", "#08306B")
if (!identical(blue_scale[[1L]], "#FFFFFF")) stop("The zero-expression color must be pure white.")
S <- vector("list", 2L)
for (i in seq_along(c("Foxp1", "Tfap2b"))) {
  gene <- c("Foxp1", "Tfap2b")[i]
  g <- read.delim(gzfile(file.path(gene_dir, paste0(gene, ".tsv.gz"))),
                  check.names = FALSE)
  ix <- match(d$point_id, g$point_id)
  if (anyNA(ix) || anyDuplicated(d$point_id) || anyDuplicated(g$point_id)) {
    stop(gene, " counts and aligned coordinates do not have a one-to-one point match.")
  }
  if (!all(g$total_umi[ix] == d$total_umi) || !all(g$gene_symbol == gene)) {
    stop(gene, " input does not match the retained MAGIC observations.")
  }
  D <- data.frame(x = xyz[, "x"], y = xyz[, "y"], z = xyz[, "z"],
                  count = g$count[ix], total_umi = g$total_umi[ix])
  D$rate <- D$count / (D$total_umi / 10000)
  lim <- c(0, stats::quantile(D$rate, 0.98, names = FALSE))
  Dz <- D[D$z %in% z_show, ]
  Dz$z_panel <- factor(sprintf("z = %.2f mm", Dz$z),
                       levels = sprintf("z = %.2f mm", z_show))
  p <- ggplot(Dz, aes(x, y, colour = rate)) +
    geom_point(size = 0.16, alpha = 0.7) +
    facet_wrap(~z_panel, nrow = 1) +
    coord_equal() +
    scale_colour_gradientn(
      colours = blue_scale, limits = lim, oob = scales::squish,
      name = "Count per\n10,000 UMI",
      guide = guide_colourbar(
        barwidth = grid::unit(2.5, "mm"), barheight = grid::unit(12, "mm"),
        title.theme = element_text(size = 5.2),
        label.theme = element_text(size = 5)
      )
    ) +
    labs(
      title = paste0(gene, ": observed expression across six z planes"),
      subtitle = sprintf(
        "All points in the same six planes; linear count/10,000 UMI; colors capped at %.2f (gene-wide P98)",
        lim[2L]
      ),
      x = "x (mm)", y = "y (mm)"
    )
  base <- file.path(out, paste0(tolower(gene), "-six-z-observed"))
  ggsave(paste0(base, ".png"), p, width = 180, height = 55,
         units = "mm", dpi = 600, bg = "white")
  svglite(paste0(base, ".svg"), width = 180 / 25.4, height = 55 / 25.4)
  print(p)
  dev.off()
  grDevices::cairo_pdf(paste0(base, ".pdf"), width = 180 / 25.4,
                        height = 55 / 25.4, family = "Arial")
  print(p)
  dev.off()
  ragg::agg_tiff(paste0(base, ".tiff"), width = 180 / 25.4,
                 height = 55 / 25.4, units = "in", res = 600,
                 background = "white")
  print(p)
  dev.off()
  S[[i]] <- data.frame(
    gene = gene, retained_observations = nrow(D), displayed_observations = nrow(Dz),
    detected_observations = sum(D$count > 0L), p98_count_per_10000 = lim[2L],
    cap_count_per_10000 = lim[2L], stringsAsFactors = FALSE
  )
}
write.csv(do.call(rbind, S), file.path(out, "marker-plot-summary.csv"), row.names = FALSE)
