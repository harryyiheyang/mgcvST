# Inspect a completed benchmark; keep this outside the timed experiment.
suppressPackageStartupMessages(library(mgcvST))
out <- Sys.getenv("MGCVST_INLA_OUTPUT", "artifacts/pathwaylgm-151673")
z <- readRDS(file.path(out, "benchmark-fits.rds"))
a <- z$inla_fit
j <- a$geometry$target[["global"]]
spatial_inla <- a$geometry$smooth[[j]]$B %*% t(a$smooth_coefficients$global)
spatial_bam <- vapply(z$bam_fits, function(fit) {
  s <- which(vapply(fit$smooth, inherits, logical(1L), "spde.smooth"))
  columns <- seq.int(fit$smooth[[s]]$first.para, fit$smooth[[s]]$last.para)
  L <- mgcvST:::.gam_training_lpmatrix(fit)
  as.numeric(L[, columns, drop = FALSE] %*% coef(fit)[columns])
}, numeric(ncol(z$Y)))
stopifnot(all(a$diagnostics$converged),
          all(vapply(z$bam_fits, function(f) isTRUE(f$converged), logical(1L))),
          all(is.finite(a$diagnostics$marginal_p_value)),
          all(is.finite(z$inla_score$results$p_two_sided)),
          all(is.finite(z$bam_score$results$p_two_sided)),
          max(abs(colMeans(spatial_inla))) < 1e-9,
          max(abs(colMeans(spatial_bam))) < 1e-9)
agreement <- data.frame(
  feature_id = rownames(z$Y),
  spatial_correlation = vapply(seq_len(nrow(z$Y)), function(k)
    cor(spatial_inla[, k], spatial_bam[, k]), numeric(1L)),
  spatial_rmse = sqrt(colMeans((spatial_inla - spatial_bam)^2)),
  inla_spatial_sd = apply(spatial_inla, 2L, sd),
  bam_spatial_sd = apply(spatial_bam, 2L, sd),
  inla_mean = colMeans(spatial_inla), bam_mean = colMeans(spatial_bam)
)
utils::write.csv(agreement, file.path(out, "spatial-agreement.csv"), row.names = FALSE)
print(agreement)
print(a$diagnostics)

png(file.path(out, "spatial-comparison.png"), width = 1900, height = 650, res = 150)
par(mfrow = c(1, 4), mar = c(1.8, 1.8, 3, .5))
xy <- z$data[c("x", "y")]
labels <- factor(z$labels)
plot(xy, asp = 1, pch = 16, cex = .32,
     col = hcl.colors(nlevels(labels), "Dark 3")[as.integer(labels)],
     axes = FALSE, xlab = "", ylab = "", main = "151673: manual cortical layers")
legend("bottomleft", levels(labels), col = hcl.colors(nlevels(labels), "Dark 3"),
       pch = 16, bty = "n", cex = .58)
k <- match("ENSG00000197971", rownames(z$Y))
if (is.na(k)) k <- 1L
v1 <- spatial_inla[, k]
v2 <- spatial_bam[, k]
lim <- range(c(v1, v2))
show_field <- function(v, limits, title) {
  color <- hcl.colors(101L, "Blue-Red 3")
  bins <- pmin(101L, pmax(1L, 1L + floor(100 * (v - limits[1L]) / diff(limits))))
  plot(xy, asp = 1, pch = 16, cex = .32, col = color[bins],
       axes = FALSE, xlab = "", ylab = "", main = title)
  mtext(sprintf("log-mean spatial component: %.2f to %.2f", limits[1L], limits[2L]),
        side = 1, cex = .63)
}
show_field(v1, lim, "MBP: INLA")
show_field(v2, lim, "MBP: bam fREML / discrete")
delta <- max(abs(v1 - v2))
show_field(v1 - v2, c(-delta, delta), "MBP: INLA minus bam")
dev.off()
writeLines(capture.output(sessionInfo()), file.path(out, "session-info.txt"))
