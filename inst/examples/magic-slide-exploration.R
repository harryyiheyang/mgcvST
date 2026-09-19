# Exploratory slide adjustment on the retained MAGIC geometry.
library(Matrix)
library(INLA)
library(fmesher)
MAGIC <- readRDS("artifacts/datasets/MAGIC/MAGIC.rds")
d <- MAGIC$covariates
s <- MAGIC$slices
out <- "artifacts/magic-slide-exploration/fits"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
mesh <- fm_mesh_3d(loc = MAGIC$meshes$native3d$loc, tv = MAGIC$meshes$native3d$tv)
k <- MAGIC$meshes$native3d$contract$kappa_fixed
spde <- inla.spde2.matern(mesh, alpha = 2, constr = TRUE,
  B.tau = matrix(c(-0.5 * log(8 * pi * k), -1), 1, 2),
  B.kappa = matrix(c(log(k), 0), 1, 2),
  theta.prior.mean = 0, theta.prior.prec = 1)
Q <- inla.spde.precision(spde, theta = 0)
C <- spde$f$extraconstr
A <- inla.spde.make.A(mesh, loc = as.matrix(d[, c("x_mm", "y_mm", "z_mm")]))
Z <- sparseMatrix(i = seq_len(nrow(d)), j = d$slice_order, x = 1,
  dims = c(nrow(d), nrow(s)))
stopifnot(max(abs(rowSums(A) - 1)) < 1e-7,
  max(abs(as.matrix(A %*% mesh$loc) - as.matrix(d[, c("x_mm", "y_mm", "z_mm")]))) < 1e-7)
flat <- list(prior = "flat", param = numeric(), initial = 0, fixed = FALSE)
hp <- list(prec = list(prior = "flat", param = numeric(), initial = log(100), fixed = FALSE))
hu <- c(hp, list(phi = list(prior = "flat", param = numeric(), initial = log(20), fixed = FALSE)))
saveRDS(list(spatial_prior = flat, nb_prior = flat, slide_prior = hp,
  ou_prior = hu, constraint = C, kappa = k, coordinates = "millimetres",
  ou_covariance = "exp(-phi * abs(z_s-z_t)) / precision",
  comparison = "exploratory in-sample comparison; no covariance tests or model selection",
  provenance = MAGIC$provenance), file.path(out, "configuration.rds"))
R1 <- list()
R2 <- list()
j <- 0L
for (gene in colnames(MAGIC$expression)) {
  y <- MAGIC$expression[, gene]
  for (model in c("spatial", "iid", "ou")) {
    j <- j + 1L
    f <- y ~ -1 + b0 + f(spatial, model = "generic0", Cmatrix = Q,
      constr = FALSE, extraconstr = C, rankdef = 1L, diagonal = 0,
      hyper = list(prec = flat))
    aa <- list(1, A)
    ee <- list(data.frame(b0 = rep(1, nrow(d))), list(spatial = seq_len(mesh$n)))
    if (model == "iid") {
      f <- update(f, . ~ . + f(slide, model = "iid", constr = FALSE, hyper = hp))
      aa <- c(aa, list(Z))
      ee <- c(ee, list(list(slide = seq_len(nrow(s)))))
    }
    if (model == "ou") {
      f <- update(f, . ~ . + f(slide, model = "ou", values = s$z / 1000,
        constr = FALSE, hyper = hu))
      aa <- c(aa, list(Z))
      ee <- c(ee, list(list(slide = s$z / 1000)))
    }
    st <- inla.stack(data = list(y = y), A = aa, effects = ee,
      compress = TRUE, remove.unused = FALSE)
    cat(gene, model, "started\n")
    t0 <- proc.time()[["elapsed"]]
    fit <- inla(f, family = "nbinomial", data = inla.stack.data(st), E = d$exposure,
      control.family = list(variant = 0, hyper = list(theta = flat)),
      control.fixed = list(mean = 0, prec = 0, mean.intercept = 0, prec.intercept = 0),
      control.predictor = list(A = inla.stack.A(st), compute = FALSE),
      control.compute = list(config = FALSE, dic = FALSE, waic = FALSE, cpo = FALSE,
        mlik = FALSE, return.marginals = FALSE, return.marginals.predictor = FALSE),
      control.inla = list(strategy = "gaussian", int.strategy = "eb",
        control.vb = list(enable = FALSE)), inla.mode = "experimental",
      num.threads = "2:1", verbose = FALSE, safe = FALSE, keep = FALSE)
    seconds <- proc.time()[["elapsed"]] - t0
    th <- setNames(fit$mode$theta, rownames(fit$summary.hyperpar))
    size <- exp(th[grep("^size", names(th), ignore.case = TRUE)])
    tau <- exp(th[grep("^Precision for spatial$", names(th))])
    u <- fit$summary.random$spatial$mean
    b <- fit$summary.fixed["b0", "mean"]
    v <- rep(0, nrow(s))
    slide_sd <- NA_real_
    phi <- NA_real_
    if (model != "spatial") {
      ids <- if (model == "ou") s$z / 1000 else seq_len(nrow(s))
      ix <- match(ids, fit$summary.random$slide$ID)
      stopifnot(!anyNA(ix))
      v <- fit$summary.random$slide$mean[ix]
      slide_sd <- exp(-th[grep("^Precision for slide$", names(th))] / 2)
    }
    if (model == "ou") phi <- exp(th[grep("^Phi for slide$", names(th))])
    spatial <- as.numeric(A %*% u)
    eta <- b + spatial + v[d$slice_order]
    mu <- d$exposure * exp(eta)
    res <- (y - mu) / sqrt(mu + mu^2 / size)
    stopifnot(length(size) == 1L, length(tau) == 1L, length(phi) == 1L,
      length(slide_sd) == 1L, all(is.finite(eta)), all(is.finite(res)),
      abs(as.numeric(C$A %*% u)) < 1e-6)
    ztab <- aggregate(cbind(y, mu, res, depth = d$total_umi, spatial),
      by = list(slice_order = d$slice_order), FUN = mean)
    ztab$slice_id <- s$slice_id[ztab$slice_order]
    ztab$z_um <- s$z[ztab$slice_order]
    ztab$observed_fitted_ratio <- ztab$y / ztab$mu
    ztab$slide_effect <- v[ztab$slice_order]
    ztab$gene <- gene
    ztab$model <- model
    R2[[j]] <- ztab
    R1[[j]] <- data.frame(gene = gene, model = model, n = length(y), seconds = seconds,
      mode_status = fit$mode$mode.status, warnings = length(fit$misc$warnings),
      nb_size = size, spatial_sd = 1 / sqrt(tau), slide_sd = slide_sd, phi_per_mm = phi,
      correlation_at_60um = if (model == "ou") exp(-phi * 0.06) else NA_real_,
      count_correlation = cor(y, mu), rate_correlation = cor(y / d$exposure, exp(eta)),
      pearson_mean = mean(res), pearson_rms = sqrt(mean(res^2)),
      slice_mean_residual_rms = sqrt(mean(ztab$res^2)),
      adjacent_slice_residual_correlation = cor(head(ztab$res, -1), tail(ztab$res, -1)),
      spatial_slide_mean_correlation = if (model == "spatial") NA_real_ else cor(ztab$spatial, v),
      slice_ratio_min = min(ztab$observed_fitted_ratio),
      slice_ratio_max = max(ztab$observed_fitted_ratio), row.names = NULL)
    saveRDS(list(gene = gene, model = model, point_id = d$point_id, metrics = R1[[j]],
      mode = fit$mode, hyper = fit$summary.hyperpar, fixed = fit$summary.fixed,
      spatial = spatial, slide_effect = v, u = u, eta = eta, mu = mu,
      warnings = fit$misc$warnings), file.path(out, paste0(gene, "-", model, ".rds")))
    write.csv(do.call(rbind, R1), file.path(out, "comparison.csv"), row.names = FALSE)
    write.csv(do.call(rbind, R2), file.path(out, "slice-diagnostics.csv"), row.names = FALSE)
    print(R1[[j]])
  }
}
writeLines(capture.output(sessionInfo()), file.path(out, "session.txt"))
