#!/usr/bin/env Rscript

# Reproducible estimator and score benchmark for the INLA and bam backends.
# Example:
# Rscript inst/benchmarks/inla-estimator.R --n=600 --mesh-side=8 \
#   --features=2 --repeats=2 --families=gaussian,negative_binomial
# Use --artifact=none to omit the large optional serialized fit snapshot.

options(stringsAsFactors = FALSE)
Sys.setenv(
  OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1",
  BLIS_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1",
  RCPP_PARALLEL_NUM_THREADS = "1"
)

parse_args <- function(x) {
  defaults <- list(
    n = 600L, mesh_side = 8L, features = 2L, repeats = 1L,
    families = "gaussian,negative_binomial", seed = 9102L,
    input = "",
    artifact = "",
    output = file.path("inst", "benchmarks", "inla-estimator-results.csv"),
    report = file.path("inst", "benchmarks", "inla-estimator-results.md")
  )
  for (arg in x) {
    z <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1L]]
    z[1L] <- gsub("-", "_", z[1L], fixed = TRUE)
    if (length(z) != 2L || !(z[1L] %in% names(defaults))) {
      stop("Unknown argument: ", arg)
    }
    defaults[[z[1L]]] <- z[2L]
  }
  for (name in c("n", "mesh_side", "features", "repeats", "seed")) {
    defaults[[name]] <- as.integer(defaults[[name]])
  }
  defaults$families <- strsplit(defaults$families, ",", fixed = TRUE)[[1L]]
  if (defaults$n < 50L || defaults$mesh_side < 4L ||
      defaults$features < 2L || defaults$repeats < 1L) {
    stop("Require n >= 50, mesh-side >= 4, features >= 2, and repeats >= 1.")
  }
  if (!all(defaults$families %in% c("gaussian", "poisson", "negative_binomial"))) {
    stop("families must use gaussian, poisson, or negative_binomial.")
  }
  defaults
}

elapsed <- function(expr) {
  t0 <- proc.time()[["elapsed"]]
  value <- force(expr)
  list(value = value, seconds = proc.time()[["elapsed"]] - t0)
}

family_object <- function(name) switch(
  name,
  gaussian = stats::gaussian(),
  poisson = stats::poisson(),
  negative_binomial = mgcv::nb()
)

simulate_y <- function(name, eta, x, features) {
  out <- vapply(seq_len(features), function(j) {
    eta_j <- eta + 0.08 * (j - 1L) * cos(2 * pi * x)
    switch(
      name,
      gaussian = eta_j + stats::rnorm(length(eta_j), sd = 0.45),
      poisson = stats::rpois(length(eta_j), exp(eta_j)),
      negative_binomial = stats::rnbinom(length(eta_j), mu = exp(eta_j), size = 4)
    )
  }, numeric(length(eta)))
  t(out)
}

bam_compact <- function(fits, feature_id) {
  states <- lapply(fits, function(fit) {
    L <- mgcvST:::.gam_training_lpmatrix(fit)
    geometry <- mgcvST:::.mgcvst_model_geometry(fit, L)
    nuisance <- mgcvST:::.mgcvst_nuisance_state(
      fit, geometry, list(L = L, frozen = TRUE)
    )
    if (is.null(nuisance)) stop("bam fit did not provide conditional nuisance covariance.")
    W <- rkhs_extract_working_model(fit)
    list(fit = fit, W = W, geometry = geometry, nuisance = nuisance)
  })
  geometry <- states[[1L]]$geometry
  geometry$nuisance_columns <- states[[1L]]$nuisance$columns
  geometry$nuisance_design <- states[[1L]]$nuisance$design
  geometry$nuisance_projection <- "conditional_Vp_block"
  diagnostics <- data.frame(
    index = seq_along(fits), feature_id = feature_id,
    converged = vapply(fits, function(x) isTRUE(x$converged), logical(1L)),
    error_message = NA_character_
  )
  structure(list(
    feature_id = feature_id,
    working_error = do.call(cbind, lapply(states, function(x) x$W$working_error)),
    working_variance = do.call(cbind, lapply(states, function(x) x$W$working_variance)),
    dispersion = stats::setNames(
      vapply(states, function(x) x$W$dispersion, numeric(1L)), feature_id
    ),
    lambda = stats::setNames(
      vapply(states, function(x) {
        j <- x$geometry$target[["global"]]
        x$geometry$sp[x$geometry$smooth[[j]]$sp_index]
      }, numeric(1L)), feature_id
    ),
    component_lambda = matrix(vapply(states, function(x) {
      j <- x$geometry$target[["global"]]
      x$geometry$sp[x$geometry$smooth[[j]]$sp_index]
    }, numeric(1L)), ncol = 1L, dimnames = list(feature_id, "global")),
    smoothing_parameters = do.call(rbind, lapply(states, function(x) x$geometry$sp)),
    family_parameters = stats::setNames(
      lapply(states, function(x) x$W$family_parameters), feature_id
    ),
    nuisance_covariance = stats::setNames(
      lapply(states, function(x) x$nuisance$covariance), feature_id
    ),
    geometry = geometry, row_id = geometry$row_id,
    score_components = geometry$score_components,
    model_setting = "global", diagnostics = diagnostics,
    test_engine = "single_model",
    timing = list(backend = "bam", workers = 1L)
  ), class = c("mgcvST_model_fit", "mgcvST_fit", "mgcvST"))
}

spatial_mean_inla <- function(fit) {
  max(abs(unlist(lapply(fit$score_components, function(component) {
    j <- fit$geometry$target[[component]]
    B <- fit$geometry$smooth[[j]]$B
    apply(fit$smooth_coefficients[[component]], 1L, function(beta) {
      mean(as.numeric(B %*% beta))
    })
  }))))
}

spatial_mean_bam <- function(fits) {
  max(abs(vapply(fits, function(fit) {
    j <- which(vapply(fit$smooth, inherits, logical(1L), "spde.smooth"))
    columns <- seq.int(fit$smooth[[j]]$first.para, fit$smooth[[j]]$last.para)
    L <- mgcvST:::.gam_training_lpmatrix(fit)
    mean(as.numeric(L[, columns, drop = FALSE] %*%
                      stats::coef(fit)[columns]))
  }, numeric(1L))))
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
suppressPackageStartupMessages({
  library(mgcvST)
  library(mgcv)
})
if (!requireNamespace("INLA", quietly = TRUE)) stop("Install INLA before benchmarking.")
if (!requireNamespace("geometry", quietly = TRUE)) stop("Install geometry before benchmarking.")

set.seed(args$seed)
input <- NULL
input_name <- "synthetic"
provenance <- "seeded synthetic spatial fields"
if (nzchar(args$input)) {
  input <- readRDS(args$input)
  if (!is.list(input) || is.null(input$data) || is.null(input$mesh) ||
      is.null(input$Y)) {
    stop("--input RDS must contain data, mesh, and Y.")
  }
  data <- as.data.frame(input$data)
  if (!all(c("x", "y", "offset0") %in% names(data)) ||
      any(!is.finite(as.matrix(data[c("x", "y", "offset0")]))) ) {
    stop("input$data must contain finite x, y, and offset0 columns.")
  }
  mesh <- input$mesh
  observed_Y <- as.matrix(input$Y)
  storage.mode(observed_Y) <- "double"
  if (ncol(observed_Y) != nrow(data) || any(!is.finite(observed_Y)) ||
      any(observed_Y < 0) || any(observed_Y != round(observed_Y))) {
    stop("input$Y must be a finite non-negative count feature-by-observation matrix.")
  }
  args$n <- nrow(data)
  args$features <- nrow(observed_Y)
  args$families <- "negative_binomial"
  input_name <- normalizePath(args$input, winslash = "/", mustWork = TRUE)
  if (!is.null(input$provenance)) {
    provenance <- paste(as.character(unlist(input$provenance)), collapse = "; ")
  }
  model_formula <- response ~ offset(offset0)
  kappa <- if (is.null(input$kappa)) 0.7 else as.numeric(input$kappa)
} else {
  vertices <- as.matrix(expand.grid(
    x = seq(0, 1, length.out = args$mesh_side),
    y = seq(0, 1, length.out = args$mesh_side)
  ))
  mesh <- list(loc = vertices, graph = list(tv = geometry::delaunayn(vertices)))
  data <- data.frame(
    x = runif(args$n, 0.01, 0.99), y = runif(args$n, 0.01, 0.99),
    z = runif(args$n, -1, 1), offset0 = runif(args$n, -0.15, 0.15)
  )
  observed_Y <- NULL
  model_formula <- response ~ z + offset(offset0)
  kappa <- 0.7
}
basis_time <- elapsed(spde_basis(
  mesh, as.matrix(data[c("x", "y")]), kappa = kappa,
  project_intercept = TRUE
))
basis <- basis_time$value
basis$component <- "global"
basis$score.component <- "global"
feature_id <- if (!is.null(observed_Y) && !is.null(rownames(observed_Y))) {
  rownames(observed_Y)
} else paste0("feature_", seq_len(args$features))
pairs <- t(utils::combn(feature_id, 2L))
rows <- list()

for (family_name in args$families) for (iteration in seq_len(args$repeats)) {
  set.seed(args$seed + 100L * match(family_name, args$families) + iteration)
  if (is.null(observed_Y)) {
    eta <- 0.55 + 0.35 * data$z + data$offset0 +
      0.25 * sin(2 * pi * data$x) - 0.2 * cos(2 * pi * data$y)
    Y <- simulate_y(family_name, eta, data$x, args$features)
  } else {
    Y <- observed_Y
  }
  dimnames(Y) <- list(feature_id, NULL)
  family <- family_object(family_name)

  inla_setup <- elapsed(inlaST.set(
    model_formula, data, basis, family = family
  ))
  historical_normal <- list(prior = "normal", param = c(0, 1 / 9), initial = 0)
  inla_fit <- elapsed(inlaST.estimate(
    Y, inla_setup$value, retain_smooth = TRUE,
    BPPARAM = BiocParallel::SerialParam(),
    control = list(
      precision_prior = historical_normal,
      nb_size_prior = historical_normal
    ), marginal_args = list(method = "liu")
  ))
  inla_score <- elapsed(mgcvST.test(
    inla_fit$value, pairs = pairs, calibration = "liu",
    BPPARAM = BiocParallel::SerialParam()
  ))
  inla_engine_seconds <- sum(
    inla_fit$value$diagnostics$fit_seconds, na.rm = TRUE
  )
  inla_marginal_seconds <- inla_fit$value$timing$marginal_elapsed
  inla_compact_seconds <- max(
    0, inla_fit$seconds - inla_engine_seconds - inla_marginal_seconds
  )

  data$response <- Y[1L, ]
  bam_formula <- if (is.null(observed_Y)) {
    response ~ z + offset(offset0) + s(x, y, bs = "spde", xt = basis)
  } else {
    response ~ offset(offset0) + s(x, y, bs = "spde", xt = basis)
  }
  bam_setup <- elapsed(mgcv::bam(
    bam_formula, data = data, family = family, method = "fREML",
    discrete = TRUE, nthreads = 1L, fit = FALSE
  ))
  response_index <- attr(bam_setup$value$terms, "response")
  bam_family_raw <- serialize(bam_setup$value$family, NULL)
  bam_fit <- elapsed(lapply(seq_len(args$features), function(j) {
    G <- bam_setup$value
    G$y <- as.numeric(Y[j, ])
    G$mf[[response_index]] <- as.numeric(Y[j, ])
    G$family <- unserialize(bam_family_raw)
    mgcv::bam(
      G = G, method = "fREML", discrete = TRUE, nthreads = 1L
    )
  }))
  bam_marginal <- elapsed(lapply(bam_fit$value, function(fit) {
    mgcvST:::taps_score_test(
      fit, test.component = 1L, method = "liu", n_threads = 1L
    )
  }))
  bam_reduce <- elapsed(bam_compact(bam_fit$value, feature_id))
  bam_score <- elapsed(mgcvST.test(
    bam_reduce$value, pairs = pairs, calibration = "liu",
    BPPARAM = BiocParallel::SerialParam()
  ))

  common <- data.frame(
    timestamp_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    input = input_name, provenance = provenance,
    family = family_name, replicate = iteration, observations = args$n,
    mesh_vertices = basis$raw_dimension, kappa = kappa,
    features = args$features, threads = 1L,
    basis_seconds = basis_time$seconds
  )
  rows[[length(rows) + 1L]] <- cbind(common, data.frame(
    backend = "inlaST", setup_seconds = inla_setup$seconds,
    fit_seconds = inla_engine_seconds,
    marginal_seconds = inla_marginal_seconds,
    estimator_seconds = inla_fit$seconds,
    compact_seconds = inla_compact_seconds,
    score_seconds = inla_score$seconds,
    total_seconds = inla_setup$seconds + inla_fit$seconds + inla_score$seconds,
    seconds_per_feature = inla_fit$seconds / args$features,
    converged = sum(inla_fit$value$diagnostics$converged),
    max_abs_spatial_observation_mean = spatial_mean_inla(inla_fit$value),
    score_p_value = min(inla_score$value$results$p_two_sided, na.rm = TRUE)
  ))
  rows[[length(rows) + 1L]] <- cbind(common, data.frame(
    backend = "bam_fREML_discrete", setup_seconds = bam_setup$seconds,
    fit_seconds = bam_fit$seconds, marginal_seconds = bam_marginal$seconds,
    estimator_seconds = bam_fit$seconds + bam_marginal$seconds +
      bam_reduce$seconds,
    compact_seconds = bam_reduce$seconds,
    score_seconds = bam_score$seconds,
    total_seconds = bam_setup$seconds + bam_fit$seconds + bam_marginal$seconds +
      bam_reduce$seconds + bam_score$seconds,
    seconds_per_feature = (bam_fit$seconds + bam_marginal$seconds +
      bam_reduce$seconds) / args$features,
    converged = sum(vapply(bam_fit$value, function(x) isTRUE(x$converged), logical(1L))),
    max_abs_spatial_observation_mean = spatial_mean_bam(bam_fit$value),
    score_p_value = min(bam_score$value$results$p_two_sided, na.rm = TRUE)
  ))

  if (!is.null(input)) {
    feature_output <- rbind(
      data.frame(
        backend = "inlaST", feature_id = feature_id,
        lambda = as.numeric(inla_fit$value$lambda),
        dispersion = as.numeric(inla_fit$value$dispersion),
        marginal_p_value = inla_fit$value$diagnostics$marginal_p_value,
        converged = inla_fit$value$diagnostics$converged
      ),
      data.frame(
        backend = "bam_fREML_discrete", feature_id = feature_id,
        lambda = vapply(bam_fit$value, function(x) as.numeric(x$sp[1L]), numeric(1L)),
        dispersion = vapply(bam_fit$value, function(x) as.numeric(x$sig2), numeric(1L)),
        marginal_p_value = vapply(
          bam_marginal$value, `[[`, numeric(1L), "smooth.pvalue"
        ),
        converged = vapply(bam_fit$value, function(x) isTRUE(x$converged), logical(1L))
      )
    )
    pair_output <- rbind(
      cbind(backend = "inlaST", inla_score$value$results),
      cbind(backend = "bam_fREML_discrete", bam_score$value$results)
    )
    feature_path <- sub("\\.csv$", "-features.csv", args$output, ignore.case = TRUE)
    pair_path <- sub("\\.csv$", "-pairs.csv", args$output, ignore.case = TRUE)
    utils::write.csv(feature_output, feature_path, row.names = FALSE)
    utils::write.csv(pair_output, pair_path, row.names = FALSE)
    if (!identical(args$artifact, "none")) {
      artifact_path <- args$artifact
      if (!nzchar(artifact_path)) {
        artifact_path <- sub("\\.csv$", "-fits.rds", args$output, ignore.case = TRUE)
      }
      dir.create(dirname(artifact_path), recursive = TRUE, showWarnings = FALSE)
      saveRDS(list(
        inla_fit = inla_fit$value, bam_fits = bam_fit$value,
        data = data, Y = Y, labels = input$labels, genes = input$genes,
        inla_score = inla_score$value, bam_score = bam_score$value,
        basis = basis, kappa = kappa, provenance = input$provenance
      ), artifact_path)
    }
  }
}

result <- do.call(rbind, rows)
dir.create(dirname(args$output), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(result, args$output, row.names = FALSE)

ratio <- aggregate(
  cbind(estimator_seconds, total_seconds) ~ family + backend,
  result, median
)
lines <- c(
  "# INLA estimator benchmark",
  "",
  paste0("Generated: ", format(Sys.time(), tz = "UTC", usetz = TRUE)),
  paste0("R: ", R.version.string),
  paste0("mgcvST: ", as.character(utils::packageVersion("mgcvST"))),
  paste0("INLA: ", as.character(utils::packageVersion("INLA"))),
  paste0("mgcv: ", as.character(utils::packageVersion("mgcv"))),
  paste0("Input: ", input_name),
  paste0("Fixed kappa: ", format(kappa, digits = 16)),
  paste0("Machine: ", Sys.info()[["sysname"]], " ", Sys.info()[["release"]],
         ", ", parallel::detectCores(logical = FALSE), " physical cores"),
  "Threads used by each estimator: 1",
  "The shared SPDE basis construction is reported separately and excluded from both totals.",
  "estimator_seconds covers fit, compact working-model construction, and mandatory Liu marginal tests for both backends; fit_seconds and marginal_seconds show the measured substeps.",
  "Both backends use the same data, mesh, kappa, observation-mean constraint, likelihood family, pair universe, and one thread. INLA uses EB with log-positive-parameter N(0, 3^2) priors; bam uses fREML, so this is a runtime comparison rather than parameter equality.",
  "The first INLA fit includes the external INLA process startup cost.",
  "",
  "```",
  paste(capture.output(print(ratio, row.names = FALSE)), collapse = "\n"),
  "```",
  "",
  "These are runtime measurements for the recorded problem sizes, not an asymptotic claim."
)
writeLines(lines, args$report)
print(result)
cat("\nWrote ", args$output, " and ", args$report, "\n", sep = "")
