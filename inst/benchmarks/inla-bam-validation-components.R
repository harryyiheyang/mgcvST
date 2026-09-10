#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
Sys.setenv(LC_ALL = "C", OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1", BLIS_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1", RCPP_PARALLEL_NUM_THREADS = "1")

root <- normalizePath(Sys.getenv("MGCVST_REAL_DATA_ROOT",
  "C:/Users/yxy1234/Downloads/magicST"), mustWork = TRUE)
out <- Sys.getenv("MGCVST_COMPONENT_OUTPUT",
  "artifacts/inla-bam-validation/components")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
lib <- normalizePath(Sys.getenv("MGCVST_VALIDATION_LIBRARY",
  "artifacts/inla-bam-validation/library"), mustWork = TRUE)
.libPaths(c(lib, .libPaths()))
Sys.setenv(R_LIBS_USER = paste(.libPaths(), collapse = .Platform$path.sep))

library(BiocParallel)
library(data.table)
library(hdf5r)
library(Matrix)
library(mgcv)
library(mgcvST)
setDTthreads(1L)
if (!requireNamespace("INLA", quietly = TRUE)) stop("INLA is required.")
if (!requireNamespace("WGCNA", quietly = TRUE)) stop("WGCNA is required.")
write.csv(data.frame(package = "mgcvST", version = as.character(packageVersion("mgcvST")),
  library = find.package("mgcvST"), package_source_commit = "f2ba738"),
  file.path(out, "package-provenance.csv"), row.names = FALSE)

workers <- 4L
flat <- list(prior = "flat", param = numeric(), initial = 0)
ctl <- list(precision_prior = flat, nb_size_prior = flat,
  num_threads = 1L, keep_fit = FALSE)

timed <- function(expr) {
  t0 <- proc.time()[["elapsed"]]
  value <- force(expr)
  list(value = value, seconds = proc.time()[["elapsed"]] - t0)
}

memory_row <- function(dataset, stage, object = NULL) {
  g <- gc()
  data.frame(dataset = dataset, stage = stage,
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S%z"),
    gc_used_mb = sum(g[, 2L]),
    object_mb = if (is.null(object)) NA_real_ else
      as.numeric(object.size(object)) / 1024^2,
    seconds = NA_real_)
}

bam_compact <- function(fits, ids) {
  states <- lapply(fits, function(fit) {
    L <- mgcvST:::.gam_training_lpmatrix(fit)
    H <- mgcvST:::.mgcvst_model_geometry(fit, L)
    N <- mgcvST:::.mgcvst_nuisance_state(fit, H, list(L = L, frozen = TRUE))
    W <- rkhs_extract_working_model(fit)
    list(W = W, H = H, N = N)
  })
  H <- states[[1L]]$H
  H$nuisance_columns <- states[[1L]]$N$columns
  H$nuisance_design <- states[[1L]]$N$design
  H$nuisance_projection <- "conditional_Vp_block"
  sp <- do.call(rbind, lapply(states, function(x) x$H$sp))
  rownames(sp) <- ids
  target.sp <- vapply(H$target, function(j) H$smooth[[j]]$sp_index,
    integer(1L))
  component.lambda <- sp[, target.sp, drop = FALSE]
  colnames(component.lambda) <- names(target.sp)
  structure(list(feature_id = ids,
    working_error = do.call(cbind, lapply(states, function(x) x$W$working_error)),
    working_variance = do.call(cbind, lapply(states, function(x) x$W$working_variance)),
    dispersion = setNames(vapply(states, function(x) x$W$dispersion,
      numeric(1L)), ids),
    lambda = component.lambda[, "global"],
    component_lambda = component.lambda, smoothing_parameters = sp,
    nuisance_covariance = setNames(lapply(states, function(x) x$N$covariance), ids),
    geometry = H, row_id = H$row_id, score_components = H$score_components,
    model_setting = "global", test_engine = "single_model",
    diagnostics = data.frame(index = seq_along(ids), feature_id = ids,
      converged = vapply(fits, function(x) isTRUE(x$converged), logical(1L)),
      error_message = NA_character_),
    timing = list(backend = "bam", workers = workers)),
    class = c("mgcvST_model_fit", "mgcvST_fit", "mgcvST"))
}

bam_payload_compact <- function(payload, ids) {
  H <- payload[[1L]]$H
  if (is.null(H)) stop("The first BAM payload must retain shared geometry.")
  H$nuisance_columns <- payload[[1L]]$N$columns
  H$nuisance_design <- payload[[1L]]$N$design
  H$nuisance_projection <- "conditional_Vp_block"
  sp <- do.call(rbind, lapply(payload, `[[`, "sp"))
  rownames(sp) <- ids
  target.sp <- vapply(H$target, function(j) H$smooth[[j]]$sp_index,
    integer(1L))
  component.lambda <- sp[, target.sp, drop = FALSE]
  colnames(component.lambda) <- names(target.sp)
  structure(list(feature_id = ids,
    working_error = do.call(cbind, lapply(payload, function(x) x$W$working_error)),
    working_variance = do.call(cbind, lapply(payload, function(x) x$W$working_variance)),
    dispersion = setNames(vapply(payload, function(x) x$W$dispersion,
      numeric(1L)), ids), lambda = component.lambda[, "global"],
    component_lambda = component.lambda, smoothing_parameters = sp,
    nuisance_covariance = setNames(lapply(payload, function(x) x$N$covariance), ids),
    geometry = H, row_id = H$row_id, score_components = H$score_components,
    model_setting = "global", test_engine = "single_model",
    diagnostics = data.frame(index = seq_along(ids), feature_id = ids,
      converged = vapply(payload, `[[`, logical(1L), "converged"),
      error_message = NA_character_),
    timing = list(backend = "bam", workers = workers)),
    class = c("mgcvST_model_fit", "mgcvST_fit", "mgcvST"))
}

obj <- readRDS(file.path(root,
  "research/results/gse159709_B_manual_two_hole_mesh_relax3.rds"))
spots <- obj$retained_spots
mesh <- obj$mesh
if (nrow(spots) != 2125L || mesh$mesh$n - 1L != 298L) {
  stop("The validated Visium-B geometry must have 2,125 spots and 298 coefficients.")
}

h5.file <- file.path(root, "research/data/miso_zebrafish_melanoma",
  "GSM4838132_Visium_Sample_B_filtered_feature_bc_matrix.h5")
h5 <- H5File$new(h5.file, mode = "r")
g <- h5[["matrix"]]
dim0 <- as.integer(g[["shape"]][])
counts <- sparseMatrix(i = as.integer(g[["indices"]][]) + 1L,
  p = as.integer(g[["indptr"]][]), x = as.numeric(g[["data"]][]), dims = dim0)
barcodes <- as.character(g[["barcodes"]][])
feature.id <- as.character(g[["features"]][["id"]][])
feature.name <- as.character(g[["features"]][["name"]][])
h5$close_all()
if (anyDuplicated(feature.id)) stop("The H5 feature identifiers must be unique.")
ii <- match(spots$barcode, barcodes)
if (anyNA(ii)) stop("Retained spot barcodes are absent from the H5 matrix.")
counts <- counts[, ii, drop = FALSE]
nz <- Matrix::rowSums(counts > 0)
sm <- Matrix::rowSums(counts)
ss <- Matrix::rowSums(counts^2)
vr <- (ss - sm^2 / ncol(counts)) / (ncol(counts) - 1L)
keep <- which(nz >= 10L & vr > 0)
if (length(keep) != 12965L) stop("RNA QC must retain exactly 12,965 genes.")
lib.size <- Matrix::colSums(counts)
if (any(lib.size <= 0)) stop("Every retained spot must have positive library size.")

basis <- spde_basis(mesh, as.matrix(data.frame(x = spots$px, y = -spots$py)),
  kappa = 0.1, project_intercept = TRUE)
if (ncol(basis$B) != 298L) stop("The common basis must have 298 coefficients.")
basis$component <- basis$score.component <- "global"

cell.file <- file.path(root,
  "research/results/gse159709_B_rctd_spot_cell_type.csv")
cell <- fread(cell.file)
ci <- match(spots$barcode, cell$barcode)
if (anyNA(ci) || anyDuplicated(cell$barcode)) {
  stop("RCTD cell types must match retained spots one-to-one.")
}
cell.type <- factor(cell$cell_type_factor[ci],
  levels = c("muscle", "tumor", "other"))
if (anyNA(cell.type) || any(table(cell.type) == 0L)) {
  stop("The cell-type model requires muscle, tumor and other levels.")
}

datasets <- strsplit(Sys.getenv("MGCVST_COMPONENT_DATASETS",
  "unadjusted,celltype"), ",", fixed = TRUE)[[1L]]
if (!length(datasets) || any(!datasets %in% c("unadjusted", "celltype"))) {
  stop("MGCVST_COMPONENT_DATASETS must select unadjusted and/or celltype.")
}
node.files <- c(
  unadjusted = "research/results/mgcvst_gse159709_B_rna_mesh298_component_nodes.csv",
  celltype = "research/results/mgcvst_gse159709_B_cell_type_component_nodes.csv")
edge.files <- c(
  unadjusted = "research/results/mgcvst_gse159709_B_rna_mesh298_component_edges.csv",
  celltype = "research/results/mgcvst_gse159709_B_cell_type_component_edges.csv")
expected <- c(unadjusted = 635L, celltype = 491L)

for (dataset in datasets) {
  d.out <- file.path(out, dataset)
  dir.create(d.out, recursive = TRUE, showWarnings = FALSE)
  stage.file <- file.path(d.out, "stages.csv")
  stage_seconds <- function(stage) {
    if (!file.exists(stage.file)) return(NA_real_)
    z <- fread(stage.file)
    target <- stage
    value <- z[stage == target & is.finite(seconds), seconds]
    if (length(value)) value[1L] else NA_real_
  }
  z <- memory_row(dataset, "input_start", counts)
  write.table(z, stage.file, sep = ",", row.names = FALSE,
    col.names = !file.exists(stage.file), append = file.exists(stage.file))

  nodes <- fread(file.path(root, node.files[[dataset]]))
  ids <- nodes[direction == "positive" & component_id == "component_001",
    feature_id]
  if (length(ids) != expected[[dataset]] || anyDuplicated(ids)) {
    stop("Unexpected largest component membership for ", dataset, ".")
  }
  jj <- match(ids, feature.id)
  if (anyNA(jj) || any(!jj %in% keep)) {
    stop("Component genes do not match the original QC opportunity set.")
  }
  Y <- as.matrix(counts[jj, , drop = FALSE])
  storage.mode(Y) <- "double"
  rownames(Y) <- ids
  dat <- data.frame(response = rep(1, nrow(spots)), x = spots$px,
    y = -spots$py, offset0 = log(lib.size / 1e4), row.names = spots$barcode)
  form <- response ~ offset(offset0)
  if (dataset == "celltype") {
    dat$cell_type <- cell.type
    form <- response ~ offset(offset0) + cell_type
  }
  input.file <- file.path(d.out, "input.rds")
  if (!file.exists(input.file)) saveRDS(list(dataset = dataset, feature_id = ids,
    feature_name = feature.name[jj], Y = Y, data = dat, basis = basis,
    component_nodes = nodes[direction == "positive" &
      component_id == "component_001"]), input.file, compress = FALSE)

  bp <- SnowParam(workers, type = "SOCK", tasks = 0L, stop.on.error = TRUE,
    progressbar = FALSE)
  bp <- bpstart(bp)

  inla.model.file <- file.path(d.out, "inla-model.rds")
  if (file.exists(inla.model.file)) {
    inla.model <- readRDS(inla.model.file)
    inla.setup.seconds <- NA_real_
  } else {
    S <- timed(inlaST.set(form, dat, basis, family = nb(link = "log"),
      control = ctl, score_backend = "sparse"))
    inla.model <- S$value
    inla.setup.seconds <- S$seconds
    saveRDS(inla.model, inla.model.file, compress = FALSE)
  }
  inla.file <- file.path(d.out, "inla-fit.rds")
  if (file.exists(inla.file)) {
    inla.fit <- readRDS(inla.file)
    inla.seconds <- stage_seconds("inla_fit_saved")
  } else {
    inla.bp <- if (dataset == "celltype") SerialParam() else bp
    I <- timed(inlaST.estimate(Y, inla.model, feature_id = ids,
      retain_marginal = TRUE, marginal_args = list(method = "liu"),
      diagnostics = TRUE, BPPARAM = inla.bp,
      chunk_size = if (dataset == "celltype") 1L else ceiling(length(ids) / workers),
      control = ctl, score_backend = "sparse"))
    inla.fit <- I$value
    inla.seconds <- I$seconds
    saveRDS(inla.fit, inla.file, compress = FALSE)
  }
  z <- memory_row(dataset, "inla_fit_saved", inla.fit)
  z$seconds <- inla.seconds
  write.table(z, stage.file, sep = ",", row.names = FALSE,
    col.names = FALSE, append = TRUE)

  dat$response <- Y[1L, ]
  bform <- update(form, . ~ . + s(x, y, bs = "spde", xt = basis))
  bam.setup.file <- file.path(d.out, "bam-setup.rds")
  if (file.exists(bam.setup.file)) {
    G <- readRDS(bam.setup.file)
    bam.setup.seconds <- NA_real_
  } else {
    B0 <- timed(bam(bform, data = dat, family = nb(link = "log"),
      method = "fREML", discrete = TRUE, nthreads = 1L, fit = FALSE))
    G <- B0$value
    bam.setup.seconds <- B0$seconds
    saveRDS(G, bam.setup.file, compress = FALSE)
  }
  bam.raw.file <- file.path(d.out, "bam-fits.rds")
  bam.payload.file <- file.path(d.out, "bam-payloads.rds")
  if (dataset == "celltype" && file.exists(bam.payload.file)) {
    payload <- readRDS(bam.payload.file)
    bam.seconds <- stage_seconds("bam_fit_saved")
  } else if (dataset == "celltype") {
    ri <- attr(G$terms, "response")
    family.raw <- serialize(G$family, NULL)
    B <- timed(bplapply(seq_along(ids), function(j, Y, G, ids, ri, family.raw) {
      Gj <- G
      Gj$y <- Y[j, ]
      Gj$mf[[ri]] <- Y[j, ]
      Gj$family <- unserialize(family.raw)
      fitj <- mgcv::bam(G = Gj, method = "fREML", discrete = TRUE,
        nthreads = 1L)
      L <- mgcvST:::.gam_training_lpmatrix(fitj)
      H <- mgcvST:::.mgcvst_model_geometry(fitj, L)
      N <- mgcvST:::.mgcvst_nuisance_state(fitj, H,
        list(L = L, frozen = TRUE))
      W <- rkhs_extract_working_model(fitj)
      marginal <- as.data.frame(mgcvST:::taps_score_test(
        fitj, test.component = 1L, method = "liu", n_threads = 1L
      ))
      marginal$feature_id <- ids[j]
      list(W = W, N = N, sp = H$sp, H = if (j == 1L) H else NULL,
        marginal = marginal, converged = isTRUE(fitj$converged),
        response_identical = identical(as.numeric(fitj$y), as.numeric(Y[j, ])))
    }, Y = Y, G = G, ids = ids, ri = ri, family.raw = family.raw,
    BPPARAM = bp))
    payload <- B$value
    bam.seconds <- B$seconds
    saveRDS(payload, bam.payload.file, compress = FALSE)
  } else if (file.exists(bam.raw.file) &&
      (!file.exists(file.path(d.out, "bam-compact-fit.rds")) ||
       !file.exists(file.path(d.out, "bam-marginal-tests.rds")))) {
    bam.fits <- readRDS(bam.raw.file)
    bam.seconds <- stage_seconds("bam_fit_saved")
  } else if (file.exists(file.path(d.out, "bam-compact-fit.rds")) &&
      file.exists(file.path(d.out, "bam-marginal-tests.rds"))) {
    bam.seconds <- stage_seconds("bam_fit_saved")
  } else {
    ri <- attr(G$terms, "response")
    family.raw <- serialize(G$family, NULL)
    B <- timed(bplapply(seq_along(ids), function(j) {
      Gj <- G
      Gj$y <- Y[j, ]
      Gj$mf[[ri]] <- Y[j, ]
      Gj$family <- unserialize(family.raw)
      mgcv::bam(G = Gj, method = "fREML", discrete = TRUE, nthreads = 1L)
    }, BPPARAM = bp))
    bam.fits <- B$value
    bam.seconds <- B$seconds
    saveRDS(bam.fits, bam.raw.file, compress = FALSE)
  }
  bam.file <- file.path(d.out, "bam-compact-fit.rds")
  if (file.exists(bam.file)) {
    bam.fit <- readRDS(bam.file)
    compact.seconds <- NA_real_
  } else {
    if (!exists("bam.fits", inherits = FALSE) &&
        !exists("payload", inherits = FALSE)) {
      stop("The BAM compact checkpoint is missing its required fit payload.")
    }
    C <- timed(if (dataset == "celltype")
      bam_payload_compact(payload, ids) else bam_compact(bam.fits, ids))
    bam.fit <- C$value
    compact.seconds <- C$seconds
    saveRDS(bam.fit, bam.file, compress = FALSE)
  }
  z <- memory_row(dataset, "bam_fit_saved", bam.fit)
  z$seconds <- bam.seconds
  write.table(z, stage.file, sep = ",", row.names = FALSE,
    col.names = FALSE, append = TRUE)

  bam.marginal.file <- file.path(d.out, "bam-marginal-tests.rds")
  if (dataset == "celltype" && !file.exists(bam.marginal.file)) {
    same.y <- vapply(payload, `[[`, logical(1L), "response_identical")
    converged <- vapply(payload, `[[`, logical(1L), "converged")
    write.csv(data.frame(feature_id = ids, response_identical = same.y,
      converged = converged), file.path(d.out, "bam-fit-integrity.csv"),
      row.names = FALSE)
    if (!all(same.y)) stop("A BAM fit response is not aligned with its feature ID.")
    if (!all(converged)) stop("Nonconverged BAM features: ",
      paste(ids[!converged], collapse = ", "), ".")
    saveRDS(do.call(rbind, lapply(payload, `[[`, "marginal")),
      bam.marginal.file, compress = FALSE)
    write.csv(data.frame(dataset = dataset, engine = "bam",
      stage = "fit_and_marginal_in_worker", seconds = bam.seconds,
      exclusive_seconds = FALSE),
      file.path(d.out, "bam-marginal-timing.csv"), row.names = FALSE)
  } else if (dataset != "celltype" && !file.exists(bam.marginal.file)) {
    if (length(bam.fits) != length(ids)) stop("The BAM fit count changed.")
    same.y <- vapply(seq_along(ids), function(j) {
      identical(as.numeric(bam.fits[[j]]$y), as.numeric(Y[j, ]))
    }, logical(1L))
    converged <- vapply(bam.fits, function(x) isTRUE(x$converged), logical(1L))
    write.csv(data.frame(feature_id = ids, response_identical = same.y,
      converged = converged), file.path(d.out, "bam-fit-integrity.csv"),
      row.names = FALSE)
    if (!all(same.y)) stop("A BAM fit response is not aligned with its feature ID.")
    if (!all(converged)) stop("Nonconverged BAM features: ",
      paste(ids[!converged], collapse = ", "), ".")
    M <- timed(lapply(seq_along(ids), function(j) {
      z <- as.data.frame(mgcvST:::taps_score_test(
        bam.fits[[j]], test.component = 1L, method = "liu", n_threads = 1L
      ))
      z$feature_id <- ids[j]
      z
    }))
    saveRDS(do.call(rbind, M$value), bam.marginal.file, compress = FALSE)
    write.csv(data.frame(dataset = dataset, engine = "bam",
      stage = "marginal", seconds = M$seconds),
      file.path(d.out, "bam-marginal-timing.csv"), row.names = FALSE)
  }
  if (exists("bam.fits", inherits = FALSE)) rm(bam.fits)
  if (exists("payload", inherits = FALSE)) rm(payload)
  if (exists("B", inherits = FALSE)) rm(B)
  if (exists("C", inherits = FALSE)) rm(C)
  rm(G)
  gc()
  pairs <- t(combn(ids, 2L))
  fits <- list(bam = bam.fit, inla = inla.fit)
  for (engine in names(fits)) {
    pair.file <- file.path(d.out, paste0(engine, "-pair-tests.rds"))
    if (!file.exists(pair.file)) {
      P <- timed(mgcvST.test(fits[[engine]], pairs = pairs,
        calibration = "liu", BPPARAM = bp,
        chunk_size = ceiling(nrow(pairs) / workers)))
      saveRDS(P$value, pair.file, compress = FALSE)
      write.csv(data.frame(dataset = dataset, engine = engine,
        stage = "pair", seconds = P$seconds),
        file.path(d.out, paste0(engine, "-pair-timing.csv")), row.names = FALSE)
    }
    marginal.file <- file.path(d.out, paste0(engine, "-marginal-tests.rds"))
    if (engine == "inla" && !file.exists(marginal.file)) {
      M <- timed(mgcvST.marginal(fits[[engine]], features = ids,
        calibration = "liu", BPPARAM = bp))
      saveRDS(M$value, marginal.file, compress = FALSE)
      write.csv(data.frame(dataset = dataset, engine = engine,
        stage = "marginal", seconds = M$seconds),
        file.path(d.out, paste0(engine, "-marginal-timing.csv")), row.names = FALSE)
    }
    wgcna.file <- file.path(d.out, paste0(engine, "-wgcna.rds"))
    if (!file.exists(wgcna.file)) {
      W <- timed(mgcvST.wgcna(fits[[engine]], indices = ids, group = "global"))
      saveRDS(W$value, wgcna.file, compress = FALSE)
      write.csv(data.frame(dataset = dataset, engine = engine,
        stage = "wgcna", seconds = W$seconds),
        file.path(d.out, paste0(engine, "-wgcna-timing.csv")), row.names = FALSE)
    }
  }
  bpstop(bp)

  bam.pair <- readRDS(file.path(d.out, "bam-pair-tests.rds"))$results
  inla.pair <- readRDS(file.path(d.out, "inla-pair-tests.rds"))$results
  key <- c("feature1", "feature2")
  pair.compare <- merge(bam.pair, inla.pair, by = key,
    suffixes = c("_bam", "_inla"), sort = FALSE)
  write.csv(pair.compare, file.path(d.out, "pair-comparison.csv"), row.names = FALSE)
  bam.W <- readRDS(file.path(d.out, "bam-wgcna.rds"))
  inla.W <- readRDS(file.path(d.out, "inla-wgcna.rds"))
  modules <- merge(bam.W$modules, inla.W$modules, by = "feature_id",
    suffixes = c("_bam", "_inla"), sort = FALSE)
  write.csv(modules, file.path(d.out, "module-comparison.csv"), row.names = FALSE)
  original.edges <- fread(file.path(root, edge.files[[dataset]]))[
    direction == "positive" & component_id == "component_001"]
  fwrite(original.edges, file.path(d.out, "original-primary-component-edges.csv"))
  timing <- data.frame(dataset = dataset, genes = length(ids), spots = ncol(Y),
    pairs = nrow(pairs), q = ncol(basis$B), engine = c("bam", "inla"),
    setup_seconds = c(bam.setup.seconds, inla.setup.seconds),
    fit_seconds = c(bam.seconds, inla.seconds),
    compact_seconds = c(compact.seconds, 0),
    fit_seconds_includes_marginal = c(dataset == "celltype", TRUE),
    stringsAsFactors = FALSE)
  write.csv(timing, file.path(d.out, "fit-timings.csv"), row.names = FALSE)
  writeLines(capture.output(sessionInfo()), file.path(d.out, "session-info.txt"))
  rm(Y, bam.fit, inla.fit, inla.model, fits,
    bam.pair, inla.pair, pair.compare, bam.W, inla.W)
  gc()
}
