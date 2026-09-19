# MAGIC Foxp1 and Tfap2b observed-versus-spatial-baseline interactive views.
library(plotly)
library(htmlwidgets)
library(jsonlite)

src <- "artifacts/inla3d-transfer/spde3d_transfer_2026-09-13/data/magic"
gene_dir <- "artifacts/datasets/MAGIC/additional-genes"
fit_dir <- "artifacts/magic-slide-exploration/fits"
out <- "artifacts/inla3d-visualization"
dir.create(out, recursive = TRUE, showWarnings = FALSE)

d <- read.delim(gzfile(file.path(src, "aligned_points.tsv.gz")), check.names = FALSE)
spec <- read_json(file.path(src, "mesh_contract.json"), simplifyVector = TRUE)
xyz <- as.matrix(d[, c("x_aligned", "y_aligned", "z")])
xyz <- sweep(sweep(xyz, 2L, spec$source_origin, "-"),
             2L, spec$source_units_per_mesh_unit, "/")
colnames(xyz) <- c("x", "y", "z")

plotly_blue_scale <- list(c(0, "#FFFFFF"), c(0.5, "#6BAED6"), c(1, "#08306B"))
plotly_colorbar <- function(title) {
  list(
    title = list(text = title, font = list(size = 10)),
    thickness = 15, len = 0.5, tickfont = list(size = 9)
  )
}

configuration <- readRDS(file.path(fit_dir, "configuration.rds"))
priors <- c(
  configuration$spatial_prior$prior,
  configuration$nb_prior$prior,
  configuration$slide_prior$prec$prior,
  configuration$ou_prior$prec$prior,
  configuration$ou_prior$phi$prior
)
if (length(priors) != 5L || anyNA(priors) || !all(priors == "flat")) {
  stop("The marker comparison requires flat priors in every saved configuration block.")
}

p3d <- plot_ly()
for (i in seq_along(c("Foxp1", "Tfap2b"))) {
  gene <- c("Foxp1", "Tfap2b")[i]
  g <- read.delim(gzfile(file.path(gene_dir, paste0(gene, ".tsv.gz"))),
                  check.names = FALSE)
  fit <- readRDS(file.path(fit_dir, paste0(gene, "-spatial.rds")))
  ix <- match(d$point_id, g$point_id)
  jx <- match(d$point_id, fit$point_id)
  if (anyNA(ix) || anyNA(jx) || anyDuplicated(d$point_id) ||
      anyDuplicated(g$point_id) || anyDuplicated(fit$point_id)) {
    stop(gene, " input, fit and aligned coordinates must have one-to-one point IDs.")
  }
  if (!all(g$total_umi[ix] == d$total_umi) || !all(g$gene_symbol == gene)) {
    stop(gene, " counts do not match the retained MAGIC observations.")
  }
  if (length(fit$mu) != length(fit$point_id) || any(!is.finite(fit$mu)) || any(fit$mu < 0)) {
    stop(gene, " spatial-baseline fit must provide non-negative finite mu values.")
  }
  observed <- g$count[ix] / (g$total_umi[ix] / 10000)
  fitted <- fit$mu[jx] / (g$total_umi[ix] / 10000)
  lim <- c(0, max(stats::quantile(observed, 0.98, names = FALSE),
                    stats::quantile(fitted, 0.98, names = FALSE)))
  hover_observed <- paste0(
    d$point_id, "<br>gene=", gene, "<br>observed count=", g$count[ix],
    "<br>observed count/10,000 UMI=", signif(observed, 4)
  )
  hover_fitted <- paste0(
    d$point_id, "<br>gene=", gene,
    "<br>spatial-baseline fitted count/10,000 UMI=", signif(fitted, 4),
    "<br>spatial field=", signif(fit$spatial[jx], 4)
  )
  p3d <- add_trace(
    p3d, x = xyz[, "x"], y = xyz[, "y"], z = xyz[, "z"], text = hover_observed,
    hoverinfo = "text", type = "scatter3d", mode = "markers",
    name = paste(gene, "observed"), visible = i == 1L,
    marker = list(size = 1.25, opacity = 0.72, color = observed,
                  colorscale = plotly_blue_scale, cmin = lim[1L], cmax = lim[2L],
                  colorbar = plotly_colorbar(paste0(gene, " observed<br>count per 10,000 UMI")))
  )
  p3d <- add_trace(
    p3d, x = xyz[, "x"], y = xyz[, "y"], z = xyz[, "z"], text = hover_fitted,
    hoverinfo = "text", type = "scatter3d", mode = "markers",
    name = paste(gene, "fitted"), visible = FALSE,
    marker = list(size = 1.25, opacity = 0.72, color = fitted,
                  colorscale = plotly_blue_scale, cmin = lim[1L], cmax = lim[2L],
                  colorbar = plotly_colorbar(paste0(gene, " spatial-baseline fit<br>count per 10,000 UMI")))
  )
}

p3d <- layout(
  p3d,
  title = list(text = "MAGIC markers: 97,830 observations", x = 0.5, xanchor = "center"),
  showlegend = FALSE,
  margin = list(t = 90),
  scene = list(xaxis = list(title = "x (mm)"), yaxis = list(title = "y (mm)"),
               zaxis = list(title = "z (mm)"), aspectmode = "data"),
  updatemenus = list(list(type = "dropdown", x = 0.02, y = 0.90, buttons = list(
    list(method = "update", args = list(list(visible = c(TRUE, FALSE, FALSE, FALSE)),
         list(title = list(text = "MAGIC markers: 97,830 observations", x = 0.5, xanchor = "center"))), label = "Foxp1 observed"),
    list(method = "update", args = list(list(visible = c(FALSE, TRUE, FALSE, FALSE)),
         list(title = list(text = "MAGIC markers: 97,830 observations", x = 0.5, xanchor = "center"))), label = "Foxp1 fitted"),
    list(method = "update", args = list(list(visible = c(FALSE, FALSE, TRUE, FALSE)),
         list(title = list(text = "MAGIC markers: 97,830 observations", x = 0.5, xanchor = "center"))), label = "Tfap2b observed"),
    list(method = "update", args = list(list(visible = c(FALSE, FALSE, FALSE, TRUE)),
         list(title = list(text = "MAGIC markers: 97,830 observations", x = 0.5, xanchor = "center"))), label = "Tfap2b fitted")
  )))
)
saveWidget(p3d, file.path(out, "magic-markers-inla3d-interactive.html"),
           selfcontained = FALSE, title = "MAGIC marker 3D spatial-baseline fits")
