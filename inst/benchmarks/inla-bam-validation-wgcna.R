#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1", RCPP_PARALLEL_NUM_THREADS = "1")
lib <- normalizePath(Sys.getenv("MGCVST_VALIDATION_LIBRARY",
  "artifacts/inla-bam-validation/library"),
  winslash = "/", mustWork = TRUE)
.libPaths(c(lib, .libPaths()))
Sys.setenv(R_LIBS_USER = paste(.libPaths(), collapse = .Platform$path.sep))
library(mgcvST)
library(mgcv)
library(BiocParallel)
library(CppMatrix)
library(WGCNA)
library(mclust)

repo <- normalizePath(Sys.getenv("MGCVST_VALIDATION_REPO",
  getwd()), winslash = "/", mustWork = TRUE)
old <- normalizePath(Sys.getenv("MGCVST_HISTORICAL_ROOT",
  "C:/Users/yxy1234/Downloads/magicST"), winslash = "/", mustWork = TRUE)
out <- Sys.getenv("MGCVST_WGCNA_OUTPUT",
  file.path(repo, "artifacts/inla-bam-validation/wgcna"))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
only <- Sys.getenv("MGCVST_WGCNA_ONLY", "")
workers <- 4L
flat <- list(prior = "flat", param = numeric(), initial = 0)
para <- list(networkType = "signed", power = 6, TOMType = "signed",
             hclustMethod = "average", minClusterSize = 20L, deepSplit = 1L)

timed <- function(expr) {
  gc(FALSE)
  t0 <- proc.time()[["elapsed"]]
  value <- force(expr)
  list(value = value, seconds = proc.time()[["elapsed"]] - t0)
}

bam_compact <- function(fits, ids) {
  states <- lapply(fits, function(fit) {
    L <- mgcvST:::.gam_training_lpmatrix(fit)
    geometry <- mgcvST:::.mgcvst_model_geometry(fit, L)
    nuisance <- mgcvST:::.mgcvst_nuisance_state(
      fit, geometry, list(L = L, frozen = TRUE)
    )
    W <- rkhs_extract_working_model(fit)
    list(W = W, geometry = geometry, nuisance = nuisance)
  })
  geometry <- states[[1L]]$geometry
  geometry$nuisance_columns <- states[[1L]]$nuisance$columns
  geometry$nuisance_design <- states[[1L]]$nuisance$design
  geometry$nuisance_projection <- "conditional_Vp_block"
  lambda <- stats::setNames(vapply(states, function(x) {
    j <- x$geometry$target[["global"]]
    x$geometry$sp[x$geometry$smooth[[j]]$sp_index]
  }, numeric(1L)), ids)
  structure(list(
    feature_id = ids,
    working_error = do.call(cbind, lapply(states, function(x) x$W$working_error)),
    working_variance = do.call(cbind, lapply(states, function(x) x$W$working_variance)),
    dispersion = stats::setNames(vapply(states, function(x) x$W$dispersion,
                                         numeric(1L)), ids),
    lambda = lambda,
    component_lambda = matrix(lambda, ncol = 1L, dimnames = list(ids, "global")),
    smoothing_parameters = do.call(rbind, lapply(states, function(x) x$geometry$sp)),
    nuisance_covariance = stats::setNames(
      lapply(states, function(x) x$nuisance$covariance), ids
    ),
    geometry = geometry, row_id = geometry$row_id,
    score_components = geometry$score_components, model_setting = "global",
    test_engine = "single_model",
    diagnostics = data.frame(index = seq_along(ids), feature_id = ids,
      converged = vapply(fits, function(x) isTRUE(x$converged), logical(1L)),
      error_message = NA_character_)
  ), class = c("mgcvST_model_fit", "mgcvST_fit", "mgcvST"))
}

common.file <- file.path(old,
  "research/results/wgcna_module_replicates/common.rds")
common <- readRDS(common.file)
basis <- common$basis
basis$component <- basis$score.component <- "global"
D <- common$covariates
T <- common$T
n <- nrow(T)
q <- ncol(T)
p <- 90L
if (n != 2125L || q != 298L || common$genes != p ||
    common$amplitude != 0.6 || common$mean_count != 12 ||
    common$nb_size != 15) {
  stop("The historical WGCNA geometry does not match the frozen design.")
}
groups <- rep(1:3, each = 30L)
ids <- paste0("g", sprintf("%03d", seq_len(p)))
arms <- data.frame(
  arm = c("strong_original", "strong_repeat", "moderate", "independent"),
  rho = c(0.6, 0.6, 0.3, 0),
  first_seed = 20260907:20260910,
  first_path = c(
    "research/results/wgcna_module_smoke/simulation.rds",
    "research/results/wgcna_module_validation/strong_repeat/simulation.rds",
    "research/results/wgcna_module_validation/moderate/simulation.rds",
    "research/results/wgcna_module_validation/independent/simulation.rds"
  )
)

inla_model <- inlaST.set(expression_label ~ offset(offset0), D, basis,
                         family = mgcv::nb(), score_backend = "sparse")
bam_data <- D
bam_data$expression_label <- 0
bam_formula <- expression_label ~ offset(offset0) +
  s(x, y, bs = "spde", xt = basis)
bam_G <- mgcv::bam(bam_formula, data = bam_data, family = mgcv::nb(),
                   method = "fREML", discrete = TRUE, nthreads = 1L,
                   fit = FALSE)
response_index <- attr(bam_G$terms, "response")
family_raw <- serialize(bam_G$family, NULL)
BP <- BiocParallel::SnowParam(workers = workers, type = "SOCK",
  progressbar = FALSE, stop.on.error = TRUE)
BP <- BiocParallel::bpstart(BP)
on.exit(BiocParallel::bpstop(BP), add = TRUE)
worker_library <- BiocParallel::bplapply(seq_len(workers), function(i) {
  list(path = find.package("mgcvST"),
       version = as.character(utils::packageVersion("mgcvST")),
       has_inla_worker = exists(".inlast_fit_feature",
                                envir = asNamespace("mgcvST"), inherits = FALSE))
}, BPPARAM = BP)
if (!all(vapply(worker_library, function(x) {
  identical(normalizePath(x$path, winslash = "/"), file.path(lib, "mgcvST")) &&
    identical(x$version, as.character(utils::packageVersion("mgcvST"))) &&
    isTRUE(x$has_inla_worker)
}, logical(1L)))) {
  stop("SOCK workers did not load the isolated validation package.")
}

for (a in seq_len(nrow(arms))) {
  for (r in seq_len(10L)) {
    seed <- if (r == 1L) arms$first_seed[a] else
      as.integer(202700000L + a * 10000L + r)
    id <- paste0(arms$arm[a], "_r", sprintf("%04d", r))
    if (nzchar(only) && id != only) next
    result.file <- file.path(out, paste0(id, ".rds"))
    run.dir <- file.path(out, id)
    dir.create(run.dir, recursive = TRUE, showWarnings = FALSE)
    chunks <- split(seq_len(p), rep(seq_len(workers), length.out = p))
    bam.checkpoint <- file.path(run.dir, "bam-fits.rds")
    bam.compact.checkpoint <- file.path(run.dir, "bam-compact.rds")
    if (file.exists(result.file)) {
      if (file.exists(bam.checkpoint) &&
          !file.exists(bam.compact.checkpoint)) {
        bam_fit <- readRDS(bam.checkpoint)
        bam_models <- unlist(
          bam_fit$value, recursive = FALSE, use.names = FALSE
        )
        names(bam_models) <- ids[unlist(chunks, use.names = FALSE)]
        bam_models <- bam_models[ids]
        bam_comp <- timed(bam_compact(bam_models, ids))
        saveRDS(bam_comp, bam.compact.checkpoint, compress = TRUE)
      }
      if (file.exists(bam.checkpoint) && !file.remove(bam.checkpoint)) {
        stop("Could not retire completed bam checkpoint for ", id, ".")
      }
      next
    }

    input.file <- if (r == 1L) {
      file.path(old, arms$first_path[a])
    } else {
      file.path(old, "research/results/wgcna_module_replicates/runs",
                id, "simulation.rds")
    }
    d <- readRDS(input.file)
    Y <- d$Y
    Z <- d$Z
    R0 <- arms$rho[a] * outer(groups, groups, "==")
    diag(R0) <- 1
    dimnames(R0) <- list(ids, ids)
    if (!identical(dim(Y), c(p, n)) || !identical(dim(Z), c(q, p)) ||
        !isTRUE(all.equal(d$R0, R0, tolerance = 0))) {
      stop("Cached simulation does not match ", id, ".")
    }
    if (!is.null(d$seed) && d$seed != seed) {
      stop("Cached simulation seed does not match ", id, ".")
    }
    rownames(Y) <- ids
    input.checkpoint <- file.path(run.dir, "input.rds")
    if (!file.exists(input.checkpoint)) {
      saveRDS(list(
        Y = Y, Z = Z, R0 = R0, seed = seed, rho = arms$rho[a],
        source = normalizePath(input.file, winslash = "/"),
        source_md5 = unname(tools::md5sum(input.file))
      ), input.checkpoint, compress = TRUE)
    }

    inla.checkpoint <- file.path(run.dir, "inla-fit.rds")
    if (file.exists(inla.checkpoint)) {
      inla_fit <- readRDS(inla.checkpoint)
    } else {
      BiocParallel::bpRNGseed(BP) <- seed
      inla_fit <- timed(inlaST.estimate(
        Y, inla_model, feature_id = ids, BPPARAM = BP, chunk_size = 1L,
        retain_smooth = FALSE, diagnostics = TRUE, score_backend = "sparse",
        control = list(precision_prior = flat, nb_size_prior = flat)
      ))
      saveRDS(inla_fit, inla.checkpoint, compress = TRUE)
    }
    if (!is.null(inla_fit$value$failures) &&
        nrow(inla_fit$value$failures)) {
      stop("INLA feature failures in ", id, ": ", nrow(inla_fit$value$failures))
    }
    if (!all(inla_fit$value$diagnostics$converged)) {
      stop("INLA feature convergence failure in ", id, ".")
    }
    inla.net.checkpoint <- file.path(run.dir, "inla-wgcna.rds")
    if (file.exists(inla.net.checkpoint)) {
      inla_net <- readRDS(inla.net.checkpoint)
    } else {
      inla_net <- timed(mgcvST.wgcna(
        inla_fit$value, indices = ids, group = "global", wgcna.para = para
      ))
      saveRDS(inla_net, inla.net.checkpoint, compress = TRUE)
    }

    if (file.exists(bam.checkpoint)) {
      bam_fit <- readRDS(bam.checkpoint)
    } else {
      bam_fit <- timed(BiocParallel::bplapply(chunks, function(ii) {
        z <- vector("list", length(ii))
        for (k in seq_along(ii)) {
          j <- ii[k]
          G <- bam_G
          G$y <- Y[j, ]
          G$mf[[response_index]] <- Y[j, ]
          G$family <- unserialize(family_raw)
          z[[k]] <- mgcv::bam(G = G, method = "fREML", discrete = TRUE,
                              nthreads = 1L,
                              control = mgcv::gam.control(maxit = 100L))
        }
        names(z) <- ids[ii]
        z
      }, BPPARAM = BP))
    }
    bam_models <- unlist(
      bam_fit$value, recursive = FALSE, use.names = FALSE
    )
    names(bam_models) <- ids[unlist(chunks, use.names = FALSE)]
    bam_models <- bam_models[ids]
    bam_identity <- vapply(seq_len(p), function(j) {
      identical(as.numeric(bam_models[[j]]$y), as.numeric(Y[j, ]))
    }, logical(1L))
    if (!all(bam_identity)) {
      stop("bam feature identity mismatch in ", id, ": ",
           paste(ids[!bam_identity], collapse = ", "), ".")
    }
    bam_converged <- vapply(
      bam_models, function(x) isTRUE(x$converged), logical(1L)
    )
    if (!all(bam_converged)) {
      stop("bam feature convergence failure in ", id, ": ",
           paste(ids[!bam_converged], collapse = ", "), ".")
    }
    if (file.exists(bam.compact.checkpoint)) {
      bam_comp <- readRDS(bam.compact.checkpoint)
    } else {
      bam_comp <- timed(bam_compact(bam_models, ids))
      saveRDS(bam_comp, bam.compact.checkpoint, compress = TRUE)
    }
    bam_net <- timed(mgcvST.wgcna(
      bam_comp$value, indices = ids, group = "global", wgcna.para = para
    ))

    Rz <- stats::cov2cor(CppMatrix::matrixMultiply(
      Z, Z, transA = TRUE
    ) / q)
    base <- list(Truth = R0, Latent = Rz)
    labels <- data.frame(feature_id = ids, truth = if (arms$rho[a] > 0)
      groups else NA_integer_)
    metrics <- list()
    networks <- list()
    for (nm in names(base)) {
      R <- base[[nm]]
      adj <- WGCNA::adjacency.fromSimilarity(R, type = "signed", power = 6)
      TOM <- WGCNA::TOMsimilarity(adj, TOMType = "signed", verbose = 0)
      H <- stats::hclust(stats::as.dist(1 - TOM), method = "average")
      cls <- dynamicTreeCut::cutreeDynamic(
        H, distM = 1 - TOM, minClusterSize = 20L, deepSplit = 1L, verbose = 0
      )
      labels[[nm]] <- cls
      networks[[nm]] <- list(correlation = R, TOM = TOM, labels = cls)
    }
    labels$bam <- bam_net$value$networks$selected$labels
    labels$inla <- inla_net$value$networks$selected$labels
    networks$bam <- bam_net$value$networks$selected
    networks$inla <- inla_net$value$networks$selected

    for (nm in names(networks)) {
      R <- networks[[nm]]$correlation
      cls <- networks[[nm]]$labels
      ii <- upper.tri(R)
      within <- ii & outer(groups, groups, "==")
      outside <- ii & outer(groups, groups, "!=")
      sizes <- table(cls[cls > 0L])
      metrics[[nm]] <- data.frame(
        id = id, arm = arms$arm[a], replicate = r, seed = seed,
        rho = arms$rho[a], estimator = nm, spots = n, q = q, genes = p,
        modules = length(sizes), grey = sum(cls == 0L),
        module_sizes = paste(as.integer(sizes), collapse = ";"),
        largest_module = if (length(sizes)) max(sizes) else 0L,
        truth_ARI = if (arms$rho[a] > 0)
          mclust::adjustedRandIndex(groups, cls) else NA_real_,
        within_R = if (arms$rho[a] > 0) mean(R[within]) else NA_real_,
        between_R = if (arms$rho[a] > 0) mean(R[outside]) else NA_real_,
        offdiag_mean = mean(R[ii]), offdiag_sd = stats::sd(R[ii]),
        R_error = norm(R - R0, "F") / sqrt(p)
      )
    }
    metrics <- do.call(rbind, metrics)
    metrics$between_estimator_ARI <- NA_real_
    metrics$between_estimator_ARI[metrics$estimator %in% c("bam", "inla")] <-
      mclust::adjustedRandIndex(labels$bam, labels$inla)
    timing <- data.frame(
      id = id, arm = arms$arm[a], replicate = r, seed = seed,
      estimator = c("bam", "inla"),
      fit_seconds = c(bam_fit$seconds + bam_comp$seconds, inla_fit$seconds),
      score_seconds = c(bam_net$value$timing$score_seconds,
                        inla_net$value$timing$score_seconds),
      network_seconds = c(bam_net$value$timing$network_seconds,
                          inla_net$value$timing$network_seconds),
      wgcna_total_seconds = c(bam_net$seconds, inla_net$seconds)
    )
    diagnostics <- list(
      bam = data.frame(feature_id = ids,
        converged = bam_converged,
        theta = vapply(bam_models, function(x) x$family$getTheta(TRUE), numeric(1L)),
        sp = vapply(bam_models, function(x) x$sp[1L], numeric(1L))),
      inla = inla_fit$value$diagnostics
    )
    saveRDS(list(
      id = id, arm = arms$arm[a], replicate = r, seed = seed,
      rho = arms$rho[a], input_file = normalizePath(input.file, winslash = "/"),
      input_md5 = unname(tools::md5sum(input.file)),
      package_version = as.character(utils::packageVersion("mgcvST")),
      package_source_commit = Sys.getenv(
        "MGCVST_PACKAGE_SOURCE_COMMIT", NA_character_
      ),
      checkout_commit = system2(
        "git", c("-C", repo, "rev-parse", "HEAD"), stdout = TRUE
      ),
      priors = list(precision_prior = flat, nb_size_prior = flat),
      wgcna_parameters = para, metrics = metrics, labels = labels,
      timing = timing, diagnostics = diagnostics, networks = networks
    ), result.file, compress = TRUE)
    if (file.exists(bam.checkpoint) && !file.remove(bam.checkpoint)) {
      stop("Could not retire completed bam checkpoint for ", id, ".")
    }
    rm(d, Y, Z, inla_fit, inla_net, bam_fit, bam_models, bam_comp, bam_net,
       networks)
    gc(FALSE)
  }
}

files <- list.files(out, pattern = "_r[0-9]{4}\\.rds$", full.names = TRUE)
results <- lapply(files, readRDS)
for (i in seq_along(results)) {
  results[[i]]$package_source_commit <-
    "f2ba73886abe7382b53e8ae9ff82148ebbdccaea"
  if (is.null(results[[i]]$checkout_commit)) {
    results[[i]]$checkout_commit <- results[[i]]$package_commit
  }
  results[[i]]$package_commit <- NULL
  saveRDS(results[[i]], files[i], compress = TRUE)
}
metrics <- do.call(rbind, lapply(results, function(x) x$metrics))
timing <- do.call(rbind, lapply(results, function(x) x$timing))
labels <- do.call(rbind, lapply(results, function(x) {
  cbind(data.frame(id = x$id, arm = x$arm, replicate = x$replicate,
                   seed = x$seed), x$labels)
}))
write.csv(metrics, file.path(out, "metrics.csv"), row.names = FALSE)
write.csv(timing, file.path(out, "timing.csv"), row.names = FALSE)
write.csv(labels, file.path(out, "labels.csv"), row.names = FALSE)
manifest <- do.call(rbind, lapply(results, function(x) data.frame(
  id = x$id, arm = x$arm, replicate = x$replicate, seed = x$seed,
  rho = x$rho, input_file = x$input_file, input_md5 = x$input_md5,
  package_version = x$package_version,
  package_source_commit = x$package_source_commit,
  checkout_commit = x$checkout_commit
)))
summary <- do.call(rbind, lapply(split(
  metrics[metrics$estimator %in% c("bam", "inla"), ],
  list(metrics$arm[metrics$estimator %in% c("bam", "inla")],
       metrics$estimator[metrics$estimator %in% c("bam", "inla")])
), function(z) data.frame(
  arm = z$arm[1L], estimator = z$estimator[1L],
  valid_datasets = nrow(z), modules_min = min(z$modules),
  modules_max = max(z$modules), grey_max = max(z$grey),
  module_positive_replicates = sum(z$modules > 0L),
  truth_ARI_median = if (all(is.na(z$truth_ARI))) NA_real_ else
    median(z$truth_ARI, na.rm = TRUE),
  truth_ARI_min = if (all(is.na(z$truth_ARI))) NA_real_ else
    min(z$truth_ARI, na.rm = TRUE),
  truth_ARI_max = if (all(is.na(z$truth_ARI))) NA_real_ else
    max(z$truth_ARI, na.rm = TRUE),
  between_estimator_ARI_min = min(z$between_estimator_ARI),
  between_estimator_ARI_median = median(z$between_estimator_ARI)
)))
convergence <- do.call(rbind, lapply(results, function(x) data.frame(
  id = x$id, arm = x$arm,
  bam_converged = sum(x$diagnostics$bam$converged),
  inla_converged = sum(x$diagnostics$inla$converged)
)))
write.csv(manifest, file.path(out, "manifest.csv"), row.names = FALSE)
write.csv(summary, file.path(out, "summary.csv"), row.names = FALSE)
write.csv(convergence, file.path(out, "convergence.csv"), row.names = FALSE)
saveRDS(list(metrics = metrics, timing = timing, labels = labels,
             manifest = manifest, summary = summary,
             convergence = convergence, files = basename(files)),
        file.path(out, "results.rds"))
