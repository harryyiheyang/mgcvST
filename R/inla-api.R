# Build the observation-centred projected basis used by the existing score
# engine and the raw sparse SPDE matrices used by INLA.
.inlast_prepare_basis <- function(basis, coordinates) {
  .spde_basis_validate(basis)
  if (is.null(basis$kappa)) {
    stop("inlaST.set() currently requires every SPDE basis to have fixed kappa.")
  }
  mesh <- list(xy = basis$mesh_vertices, tv = basis$mesh_triangles)
  coordinates <- as.matrix(coordinates)
  storage.mode(coordinates) <- "double"
  if (ncol(coordinates) != 2L || !nrow(coordinates) || any(!is.finite(coordinates))) {
    stop("The model coordinates must be a finite two-column matrix.")
  }
  loc <- sweep(coordinates, 2L, basis$transform$center, "-") /
    basis$transform$scale
  A <- .spde_basis_project(mesh, loc)
  fem <- .spde_basis_fem(mesh)
  Q <- basis$kappa^4 * fem$M0 + 2 * basis$kappa^2 * fem$M1 + fem$M2
  Q <- Matrix::forceSymmetric(Q)

  # This is the observation mean, not a mesh-node sum-to-zero constraint.
  g <- as.numeric(Matrix::crossprod(A, rep(1 / nrow(A), nrow(A))))
  qg <- qr(matrix(g, ncol = 1L))
  if (qg$rank != 1L || length(g) < 2L) {
    stop("The observation mean constraint has invalid rank.")
  }
  Z <- qr.Q(qg, complete = TRUE)[, -1L, drop = FALSE]
  projected_penalty <- list(
    crossprod(Z, as.matrix(fem$M0 %*% Z)),
    2 * crossprod(Z, as.matrix(fem$M1 %*% Z)),
    crossprod(Z, as.matrix(fem$M2 %*% Z))
  )
  projected_penalty <- lapply(projected_penalty, function(x) (x + t(x)) / 2)
  projected_Q <- basis$kappa^4 * projected_penalty[[1L]] +
    basis$kappa^2 * projected_penalty[[2L]] + projected_penalty[[3L]]

  constrained <- basis
  constrained$coordinates <- coordinates
  constrained$coordinate_keys <- .spde_coordinate_keys(coordinates)
  constrained$B <- as.matrix(A %*% Z)
  constrained$penalty <- projected_penalty
  constrained$Q <- (projected_Q + t(projected_Q)) / 2
  constrained$projection <- Z
  constrained$projection_rank <- 1L
  constrained$project_intercept <- TRUE
  constrained$raw_dimension <- ncol(A)
  # Cached PCs refer to the old projection, so rebuild them if requested later.
  constrained$pc_values <- constrained$pc_vectors <- constrained$pc_cumulative <- NULL
  constrained$pc_training_basis <- constrained$pc_mesh_projection <- NULL
  constrained$pc_cached_dimension <- NULL

  list(
    basis = constrained,
    raw = list(
      A = methods::as(A, "dgCMatrix"), Q = methods::as(Q, "dsCMatrix"),
      constraint = g, projection = Z
    )
  )
}

.inlast_family <- function(family) {
  label <- tolower(family$family)
  link <- tolower(family$link)
  if (identical(label, "gaussian") && identical(link, "identity")) return("gaussian")
  if (identical(label, "poisson") && identical(link, "log")) return("poisson")
  if (grepl("^negative binomial", label) && identical(link, "log")) {
    return("negative_binomial")
  }
  stop("inlaST.set() supports Gaussian(identity), Poisson(log), and negative-binomial(log) families.")
}

# Convert frozen mgcv geometry to an INLA latent-model specification. Keeping
# this conversion here makes the public object independently serializable.
.inlast_model_spec <- function(model, raw_component) {
  geometry <- model$geometry
  n <- nrow(model$L)
  all_smooth <- unique(unlist(lapply(geometry$smooth, `[[`, "columns"),
                                    use.names = FALSE))
  parametric <- setdiff(seq_len(ncol(model$L)), all_smooth)
  tested <- sort(unique(unlist(lapply(
    geometry$smooth[geometry$target], `[[`, "columns"), use.names = FALSE
  )))
  nuisance <- setdiff(seq_len(ncol(model$L)), tested)
  geometry$nuisance_columns <- nuisance
  geometry$nuisance_design <- model$L[, nuisance, drop = FALSE]
  geometry$nuisance_projection <- "conditional_INLA_block"
  model$geometry <- geometry

  random <- vector("list", length(geometry$smooth))
  for (j in seq_along(geometry$smooth)) {
    sm <- geometry$smooth[[j]]
    if (sm$fixed || length(sm$penalties) != 1L || length(sm$sp_index) != 1L) {
      stop("inlaST.set() currently requires one fitted penalty for every smoother; unsupported smooth: '",
           sm$label, "'.")
    }
    component <- names(geometry$target)[match(j, unname(geometry$target))]
    target <- length(component) == 1L && !is.na(component)
    if (target) {
      raw <- raw_component[[component]]
      random[[j]] <- list(
        name = component, A = raw$A, Q = raw$Q, kind = "spde",
        target = TRUE, constraint = raw$constraint,
        projection = raw$projection, geometry_index = j,
        sp_index = sm$sp_index, rankdef = 0L
      )
    } else {
      nuisance_Q <- Matrix::forceSymmetric(
        Matrix::Matrix(sm$penalties[[1L]], sparse = TRUE)
      )
      rankdef <- ncol(nuisance_Q) - as.integer(Matrix::rankMatrix(nuisance_Q))
      random[[j]] <- list(
        name = paste0("nuisance_", j),
        A = methods::as(Matrix::Matrix(sm$B, sparse = TRUE), "dgCMatrix"),
        Q = nuisance_Q,
        kind = "nuisance", target = FALSE, constraint = NULL,
        projection = NULL, geometry_index = j, sp_index = sm$sp_index,
        rankdef = rankdef
      )
    }
  }

  nuisance_map <- vector("list", length(nuisance))
  for (k in seq_along(nuisance)) {
    column <- nuisance[k]
    fixed_index <- match(column, parametric)
    if (!is.na(fixed_index)) {
      nuisance_map[[k]] <- list(source = "fixed", block = NA_integer_,
                                index = fixed_index, full_column = column)
    } else {
      block <- which(vapply(geometry$smooth, function(sm) column %in% sm$columns,
                           logical(1L)))
      if (length(block) != 1L) stop("Could not map a nuisance coefficient into INLA geometry.")
      nuisance_map[[k]] <- list(
        source = "random", block = block,
        index = match(column, geometry$smooth[[block]]$columns), full_column = column
      )
    }
  }
  random_offset <- cumsum(c(ncol(geometry$X),
                            vapply(random, function(x) ncol(x$A), integer(1L))))
  nuisance_index <- vapply(nuisance_map, function(x) {
    if (identical(x$source, "fixed")) x$index else
      random_offset[x$block] + x$index
  }, integer(1L))
  combined_design <- do.call(cbind, c(list(geometry$X), lapply(random, `[[`, "A")))
  if (length(nuisance_index) != ncol(geometry$nuisance_design) ||
      !isTRUE(all.equal(as.matrix(combined_design[, nuisance_index, drop = FALSE]),
                        geometry$nuisance_design, tolerance = 1e-10))) {
    stop("The INLA nuisance coefficient map does not match the score geometry.")
  }

  fixed_names <- colnames(geometry$X)
  if (!ncol(geometry$X)) fixed_names <- character()
  spec <- list(
    n = n,
    family = .inlast_family(model$G$family),
    fixed = list(X = geometry$X, names = fixed_names),
    random = random,
    nuisance_design = geometry$nuisance_design,
    nuisance_map = nuisance_map,
    nuisance_index = nuisance_index,
    geometry_sp_length = length(geometry$sp),
    sp_names = names(model$G$sp),
    offset = model$offset,
    mean_constraint = "observation",
    mean_constraint_active = TRUE
  )
  if (spec$family == "negative_binomial" &&
      isTRUE(model$G$family$n.theta == 0)) {
    spec$nb_size_fixed <- as.numeric(model$G$family$getTheta(trans = TRUE))
  }
  list(model = model, spec = spec)
}

#' Prepare a sparse INLA estimator for mgcvST
#'
#' Constructs the same fixed-kappa SPDE score geometry as [model.set()] while
#' retaining raw sparse mesh precision matrices for INLA. Each spatial field is
#' always constrained to have zero mean at the observed locations.
#'
#' @details Complete formulas and frozen `G` designs use the shared design
#' preparation of [mgcvST.set()]. Alternatively, supply SPDE component(s)
#' separately through `basis`, as in [model.set()]. Formula/basis setup enforces
#' the observation mean-zero constraint. A supplied `G` must already have a
#' compatible constrained SPDE design; otherwise rebuild it from formula/data.
#' Supported families are Gaussian with identity link, Poisson
#' with log link, and negative binomial with log link. Full fixed-kappa SPDE
#' bases and nuisance smoothers with one penalty are supported. Spatial
#' mean-zero constraints cannot be disabled. Observation weights and
#' additional cross-penalties are not supported.
#'
#' @inheritParams model.set
#' @param G Optional frozen `gam.prefit` design, supplied instead of formula,
#'   data and basis. Its full SPDE terms must satisfy the observation constraint.
#' @param control Named engine controls stored in the model and inherited by
#'   [inlaST.estimate()]. See that function for priors and numerical controls.
#' @param score_backend Default score implementation inherited by
#'   [inlaST.estimate()]: `"auto"`, `"dense"`, or `"sparse"`.
#' @param precision_scale Scale on which the spatial log-precision prior is
#'   defined. `"raw"` retains the original FEM precision parameterization.
#'   `"observation"` normalizes each constrained spatial field to unit mean
#'   marginal variance at the observed locations when its internal precision
#'   equals one. A user-supplied proper prior then applies to this standardized
#'   precision and generally changes the prior on the original FEM multiplier.
#'   The default flat log prior is invariant to this constant log-scale shift.
#' @return An `inlaST_model` for [inlaST.estimate()].
#' @export
inlaST.set <- function(
    formula = NULL, data = NULL, basis = NULL, family = mgcv::nb(),
    setting = c("global", "global_local"), coordinates = c("x", "y"),
    precision_scale = c("raw", "observation"), G = NULL, control = list(),
    score_backend = c("auto", "dense", "sparse"), ...) {
  t0 <- proc.time()[["elapsed"]]
  setting_supplied <- !missing(setting)
  family_supplied <- !missing(family)
  # Accept mgcvST.set(formula, data, family) while preserving the historical
  # third-position prepared-basis shorthand.
  if (inherits(basis, "family") || inherits(basis, "extended.family")) {
    if (family_supplied) stop("The family was supplied twice.")
    family <- basis
    family_supplied <- TRUE
    basis <- NULL
  }
  requested_setting <- match.arg(setting)
  precision_scale <- match.arg(precision_scale)
  score_backend <- match.arg(score_backend)
  control <- .inlast_control(control)
  dots <- list(...)
  forbidden <- intersect(names(dots),
                         c("weights", "prior.weights", "subset", "na.action",
                           "paraPen", "H", "method", "fit"))
  if (length(forbidden)) {
    stop("inlaST.set() does not support these setup arguments: ",
         paste(forbidden, collapse = ", "), ".")
  }
  if (!is.null(G)) {
    if (!is.null(formula) || !is.null(data) || !is.null(basis) ||
        family_supplied || length(dots)) {
      stop("Supply G alone as the design, or formula, data, family and setup arguments.")
    }
    base <- .mgcvst_set_prepare(G = G, .allow_poisson = TRUE)
    prepared <- .inlast_prepare_frozen_components(base)
  } else if (is.null(basis)) {
    data <- as.data.frame(data)
    if (is.null(data) || !nrow(data)) stop("data must be one shared non-empty data frame.")
    complete <- .inlast_prepare_formula(formula, data)
    base <- .mgcvst_set_prepare(complete$formula, data, family,
                                ..., .allow_poisson = TRUE)
    prepared <- complete$prepared[base$components]
  } else {
    components <- if (requested_setting == "global") "global" else c("global", "local")
    supplied <- if (requested_setting == "global") list(global = basis) else basis
    if (!is.list(supplied) || !all(components %in% names(supplied))) {
      stop("basis must supply the spatial component(s) selected by setting.")
    }
    data <- as.data.frame(data)
    if (!all(coordinates %in% names(data))) {
      stop("Both coordinate columns must be present in data.")
    }
    xy <- as.matrix(data[, coordinates, drop = FALSE])
    prepared <- lapply(supplied[components], .inlast_prepare_basis,
                       coordinates = xy)
    constrained_basis <- lapply(prepared, `[[`, "basis")
    if (requested_setting == "global") constrained_basis <- constrained_basis$global
    base <- model.set(
      formula = formula, data = data, basis = constrained_basis,
      family = family, setting = requested_setting, coordinates = coordinates, ...
    )
    frozen <- .mgcvst_freeze_geometry(base$G)
    base$L <- frozen$L
    base$geometry <- frozen$geometry
    base$shared_design <- TRUE
    base$timing <- list(setup_seconds = frozen$elapsed,
                        elapsed = frozen$elapsed)
  }
  if (setting_supplied && !identical(base$setting, requested_setting)) {
    stop("setting does not match the global/local SPDE components in the design.")
  }
  if (!is.null(base$G$w) && any(base$G$w != 1)) {
    stop("inlaST.set() currently requires unit observation weights.")
  }
  built <- .inlast_model_spec(base, lapply(prepared, `[[`, "raw"))
  model <- built$model
  model$inla_spec <- built$spec
  spatial_scale <- if (precision_scale == "observation") {
    vapply(prepared, .inlast_observation_precision_scale, numeric(1L))
  } else stats::setNames(rep(1, length(prepared)), names(prepared))
  for (j in seq_along(model$inla_spec$random)) {
    block <- model$inla_spec$random[[j]]
    model$inla_spec$random[[j]]$precision_scale <- if (isTRUE(block$target)) {
      unname(spatial_scale[block$name])
    } else 1
  }
  model$precision_scale <- precision_scale
  model$inla_spec$precision_scale_mode <- precision_scale
  model$inla_spec$spatial_precision_scale <- spatial_scale
  model$mean_constraint <- "observation"
  model$mean_constraint_active <- TRUE
  model$inla_control <- .inlast_family_control(model, control)
  model$score_backend <- score_backend
  if (score_backend == "sparse") {
    capability <- .inlast_sparse_score_capability(model)
    if (!capability$eligible) {
      stop("score_backend = 'sparse' is unavailable: ", capability$reason, ".")
    }
  }
  model$timing$elapsed <- proc.time()[["elapsed"]] - t0
  class(model) <- c("inlaST_model", "mgcvST_model")
  model
}

.inlast_validate_estimate <- function(Y, model, feature_id, BPPARAM, chunk_size,
                                      offset, control, retain_smooth, diagnostics) {
  if (!inherits(model, "inlaST_model") || is.null(model$inla_spec)) {
    stop("model must be returned by inlaST.set().")
  }
  if (!isTRUE(model$mean_constraint_active) ||
      !identical(model$mean_constraint, "observation") ||
      !isTRUE(model$inla_spec$mean_constraint_active)) {
    stop("The required observation mean-zero constraint is not active.")
  }
  Y <- as.matrix(Y)
  storage.mode(Y) <- "double"
  if (length(dim(Y)) != 2L || !nrow(Y) || !ncol(Y) || any(!is.finite(Y))) {
    stop("Y must be a non-empty finite numeric feature-by-observation matrix.")
  }
  if (ncol(Y) != model$inla_spec$n) {
    stop("ncol(Y) must equal the number of observations in model.")
  }
  if (model$inla_spec$family %in% c("poisson", "negative_binomial") &&
      (any(Y < 0) || any(Y != round(Y)))) {
    stop("Count responses must be non-negative integers.")
  }
  if (is.null(feature_id)) feature_id <- rownames(Y)
  if (is.null(feature_id)) feature_id <- as.character(seq_len(nrow(Y)))
  feature_id <- as.character(feature_id)
  if (length(feature_id) != nrow(Y) || anyNA(feature_id) ||
      any(!nzchar(feature_id)) || anyDuplicated(feature_id)) {
    stop("feature_id must contain one unique non-empty identifier per feature.")
  }
  if (!inherits(BPPARAM, "BiocParallelParam")) {
    stop("BPPARAM must inherit from 'BiocParallelParam'.")
  }
  if (is.null(chunk_size)) {
    chunk_size <- ceiling(nrow(Y) / max(1L, min(nrow(Y), BiocParallel::bpworkers(BPPARAM))))
  }
  if (!is.numeric(chunk_size) || length(chunk_size) != 1L ||
      !is.finite(chunk_size) || chunk_size < 1 || chunk_size != as.integer(chunk_size)) {
    stop("chunk_size must be one positive integer.")
  }
  if (!is.list(control)) stop("control must be a list of INLA controls.")
  for (name in c("retain_smooth", "diagnostics")) {
    value <- get(name)
    if (!is.logical(value) || length(value) != 1L || is.na(value)) {
      stop(name, " must be TRUE or FALSE.")
    }
  }
  if (!is.null(offset)) {
    valid <- is.numeric(offset) && all(is.finite(offset)) &&
      ((is.null(dim(offset)) && length(offset) == ncol(Y)) ||
       (is.matrix(offset) && identical(dim(offset), dim(Y))))
    if (!valid) {
      stop("offset must be finite numeric: an observation-length vector or a matrix matching Y.")
    }
  }
  list(Y = Y, feature_id = feature_id, chunk_size = as.integer(chunk_size))
}

.inlast_fit_chunk <- function(index, Y, spec, base_offset, extra_offset,
                              control, diagnostics) {
  lapply(index, function(j) {
    total_offset <- base_offset
    if (!is.null(extra_offset)) {
      total_offset <- total_offset + if (is.matrix(extra_offset)) extra_offset[j, ] else extra_offset
    }
    tryCatch(
      .inlast_fit_feature(spec, Y[j, ], offset = total_offset, control = control,
                          diagnostics = diagnostics),
      error = function(e) e
    )
  })
}

.inlast_as_gam <- function(result, model, y, offset) {
  fit <- model$G
  response <- attr(fit$terms, "response")
  fit$y <- as.numeric(y)
  fit$mf[[response]] <- as.numeric(y)
  fit$model <- fit$mf
  fit$linear.predictors <- as.numeric(result$eta)
  fit$fitted.values <- as.numeric(result$mu)
  fit$offset <- as.numeric(offset)
  fit$prior.weights <- rep(1, length(y))
  fit$sig2 <- as.numeric(result$dispersion)
  fit$sp <- as.numeric(result$smoothing_parameters)
  names(fit$sp) <- model$inla_spec$sp_names
  coefficient <- stats::setNames(numeric(ncol(model$L)), colnames(model$L))
  smooth_columns <- unique(unlist(lapply(model$geometry$smooth, `[[`, "columns"),
                                        use.names = FALSE))
  parametric <- setdiff(seq_len(ncol(model$L)), smooth_columns)
  coefficient[parametric] <- result$fixed_mode
  for (j in seq_along(model$geometry$smooth)) {
    columns <- model$geometry$smooth[[j]]$columns
    block <- model$inla_spec$random[[j]]
    coefficient[columns] <- if (isTRUE(block$target)) {
      result$coefficients[[block$name]]
    } else {
      result$random_mode[[j]]
    }
  }
  fit$coefficients <- coefficient
  if (model$inla_spec$family == "negative_binomial") {
    fit$family$putTheta(log(result$family_parameters[1L]))
  }
  fit$.taps_score_X <- model$L
  class(fit) <- c("gam", "glm", "lm")
  fit
}

#' Estimate mgcvST working models with sparse INLA
#'
#' Fits one latent Gaussian model per feature with INLA and returns the compact
#' working-model contract consumed by [mgcvST.test()]. No GAM fit is used.
#'
#' @param Y Numeric feature-by-observation matrix.
#' @param model An object returned by [inlaST.set()].
#' @param feature_id Unique feature identifiers.
#' @param BPPARAM A `BiocParallelParam` controlling feature-level parallelism.
#' @param chunk_size Positive number of features per task.
#' @param offset Optional shared observation offset or matrix matching `Y`.
#' @param control Named INLA engine overrides for controls saved by
#'   [inlaST.set()]. Omitted entries inherit model settings. Each supplied
#'   prior list replaces the whole prior; numerical `control.inla` entries
#'   merge by name. Explicit `NULL` resets optional fixed parameters, except
#'   NB size fixed by the model family. The supported approximation is
#'   `int_strategy = "eb"`, `latent_strategy = "gaussian"`. `num_threads`
#'   defaults to one. Optional positive `fixed_precision`,
#'   `gaussian_precision`, and `nb_size` fix latent precision multipliers,
#'   inverse Gaussian residual variance, and NB size, respectively.
#'   `fixed_precision` values always refer to the original FEM multiplier,
#'   including when the model uses observation-scale precision priors. Prior
#'   lists `precision_prior`, `gaussian_precision_prior` and `nb_size_prior`
#'   contain `prior`, `param`, and logarithmic `initial` values (default zero).
#'   Spatial precision and NB size default to `prior = "flat"` with no
#'   parameters on the log scale. Gaussian observation precision defaults to
#'   `prior = "normal", param = c(0, 1/9)`, encoding
#'   `log(parameter) ~ N(0, 3^2)`. Custom normal, registered scalar INLA priors,
#'   and INLA expression/table priors are supported. Normal parameters are
#'   mean and precision. A flat log-hyperparameter objective corresponds to
#'   density proportional to `1/parameter` on its positive scale.
#'   `control.inla` accepts supported numerical tuning, with Gaussian latent
#'   strategy and EB integration enforced. Unknown controls are rejected.
#' @param retain_smooth Retain estimated score-component coefficients.
#' @param diagnostics Retain per-feature INLA diagnostics and compute the
#'   expected-Fisher nuisance covariance for comparison with native `Vp`.
#'   When `FALSE`, that extra solve is skipped and its per-feature entries
#'   in `expected_nuisance_covariance` are `NULL`.
#' @param marginal_test Optional marginal-score callback. `NULL` uses the
#'   package-local TAPS score test.
#' @param marginal_args Named arguments passed to the marginal score test.
#' @param retain_marginal Retain the frozen state needed by [mgcvST.marginal()].
#' @param score_backend Score-state implementation. `"auto"` uses the sparse
#'   mesh backend for one global fixed-kappa SPDE with fixed-effect nuisance
#'   terms only, and otherwise uses `"dense"`. `"sparse"` requires that
#'   supported structure and errors before fitting when it is unavailable.
#'   When omitted, the setting saved by [inlaST.set()] is used.
#' @details Hyperparameter priors remain part of INLA's empirical-Bayes
#' estimates; these are not mgcv REML estimates. `mgcv::nb(theta = value)`
#' fixes NB size, and a conflicting `control$nb_size` is rejected. The
#' working model uses conditional latent estimates and expected Fisher
#' variances. Its nuisance `Vp` block is extracted directly from INLA's
#' constrained conditional Gaussian posterior precision, rather than rebuilt
#' from the expected Fisher matrix. The existing score calibration is applied to these inputs;
#' this does not establish finite-sample calibration after hyperparameter
#' estimation. Every spatial component's constraint residual and observed
#' spatial mean are retained in the result.
#' The score uses the SPDE covariance conditioned on observation mean zero.
#' Centering this kernel again has no effect. The native nuisance `Vp` with
#' expected working variances does not in general guarantee `P 1 = 0`, so an
#' unconstrained raw kernel cannot be substituted using centering invariance.
#' Flat hyperpriors need not yield proper hyperparameter posteriors. They are
#' supported only as empirical-Bayes optimization objectives in the current
#' single-configuration engine. A returned finite precision or zero optimizer
#' status does not establish that the maximum is interior; zero spatial
#' variance and the Poisson limit of the NB model require boundary checks.
#' @return An `inlaST_fit` that is also an `mgcvST_model_fit`.
#' @export
inlaST.estimate <- function(
    Y, model, feature_id = rownames(Y),
    BPPARAM = BiocParallel::SerialParam(), chunk_size = NULL,
    offset = NULL, control = list(), retain_smooth = FALSE,
    diagnostics = FALSE, marginal_test = NULL, marginal_args = list(),
    retain_marginal = FALSE,
    score_backend = c("auto", "dense", "sparse")) {
  if (!requireNamespace("INLA", quietly = TRUE)) {
    stop("inlaST.estimate() requires the INLA package.")
  }
  checked <- .inlast_validate_estimate(
    Y, model, feature_id, BPPARAM, chunk_size, offset, control,
    retain_smooth, diagnostics
  )
  if (missing(score_backend) && !is.null(model$score_backend)) {
    score_backend <- model$score_backend
  }
  score_backend_requested <- match.arg(score_backend, c("auto", "dense", "sparse"))
  sparse_capability <- .inlast_sparse_score_capability(model)
  if (identical(score_backend_requested, "sparse") &&
      !sparse_capability$eligible) {
    stop("score_backend = 'sparse' is unavailable: ",
         sparse_capability$reason, ".")
  }
  score_backend <- if (identical(score_backend_requested, "auto")) {
    if (sparse_capability$eligible) "sparse" else "dense"
  } else score_backend_requested
  score_sparse <- if (identical(score_backend, "sparse")) {
    .inlast_sparse_score_geometry(model)
  } else NULL
  if (!is.null(marginal_test) && !is.function(marginal_test)) {
    stop("marginal_test must be NULL or a function.")
  }
  if (!is.list(marginal_args) ||
      (length(marginal_args) && (is.null(names(marginal_args)) ||
                                 any(!nzchar(names(marginal_args)))))) {
    stop("marginal_args must be a named list.")
  }
  forbidden_marginal <- intersect(names(marginal_args),
                                  c("fit", "test.component", "n_threads"))
  if (length(forbidden_marginal)) {
    stop("Do not supply these arguments through marginal_args: ",
         paste(forbidden_marginal, collapse = ", "))
  }
  if (!is.logical(retain_marginal) || length(retain_marginal) != 1L ||
      is.na(retain_marginal)) stop("retain_marginal must be TRUE or FALSE.")
  # Validate once in the parent so a misspelled or invalid control cannot turn
  # every feature into an otherwise opaque per-feature failure.
  control <- .inlast_control(.inlast_merge_control(model$inla_control, control))
  control <- .inlast_family_control(model, control)
  Y <- checked$Y
  feature_id <- checked$feature_id
  chunk_size <- checked$chunk_size
  groups <- split(seq_len(nrow(Y)), ceiling(seq_len(nrow(Y)) / chunk_size))
  workers <- max(1L, min(length(groups), BiocParallel::bpworkers(BPPARAM)))
  t0 <- proc.time()[["elapsed"]]
  chunks <- BiocParallel::bplapply(
    groups, .inlast_fit_chunk, Y = Y, spec = model$inla_spec,
    base_offset = model$offset, extra_offset = offset, control = control,
    diagnostics = diagnostics, BPPARAM = BPPARAM
  )
  fit_elapsed <- proc.time()[["elapsed"]] - t0
  fits <- unlist(chunks, recursive = FALSE)

  n <- ncol(Y)
  p <- nrow(Y)
  E <- V <- matrix(NA_real_, n, p, dimnames = list(NULL, feature_id))
  dispersion <- stats::setNames(rep(NA_real_, p), feature_id)
  family_parameters <- stats::setNames(vector("list", p), feature_id)
  n_sp <- model$inla_spec$geometry_sp_length
  smoothing_parameters <- matrix(
    NA_real_, p, n_sp,
    dimnames = list(feature_id, model$inla_spec$sp_names)
  )
  nuisance_covariance <- stats::setNames(vector("list", p), feature_id)
  expected_nuisance_covariance <- stats::setNames(vector("list", p), feature_id)
  diagnostics_table <- data.frame(
    index = seq_len(p), feature_id = feature_id, converged = FALSE,
    marginal_p_value = NA_real_, marginal_requested_method = NA_character_,
    marginal_method = NA_character_, marginal_fallback = NA,
    residual_df = NA_real_, criterion = NA_real_,
    criterion_name = "INLA log marginal likelihood", fit_seconds = NA_real_,
    outer_convergence = NA_character_, error_class = NA_character_,
    error_message = NA_character_, error_call = NA_character_,
    stringsAsFactors = FALSE
  )
  coefficient <- if (retain_smooth) {
    lapply(model$geometry$target, function(j) {
      matrix(NA_real_, p, ncol(model$geometry$smooth[[j]]$B),
             dimnames = list(feature_id, NULL))
    })
  } else NULL
  if (!is.null(coefficient)) names(coefficient) <- names(model$geometry$target)
  marginal_state <- stats::setNames(vector("list", p), feature_id)
  marginal_geometry <- NULL
  marginal_seconds <- numeric(p)
  constraint_residual <- observation_spatial_mean <- matrix(
    NA_real_, p, length(model$inla_spec$random),
    dimnames = list(feature_id, vapply(model$inla_spec$random, `[[`, character(1L), "name"))
  )
  inla_diagnostics <- if (diagnostics) stats::setNames(vector("list", p), feature_id) else NULL

  for (j in seq_len(p)) {
    z <- fits[[j]]
    if (inherits(z, "condition")) {
      diagnostics_table$error_class[j] <- class(z)[1L]
      diagnostics_table$error_message[j] <- conditionMessage(z)
      diagnostics_table$error_call[j] <- paste(deparse(conditionCall(z)), collapse = " ")
      next
    }
    E[, j] <- z$working_error
    V[, j] <- z$working_variance
    dispersion[j] <- z$dispersion
    family_parameters[[j]] <- z$family_parameters
    smoothing_parameters[j, ] <- z$smoothing_parameters
    nuisance_covariance[j] <- list(z$nuisance_covariance)
    expected_nuisance_covariance[j] <- list(z$expected_nuisance_covariance)
    diagnostics_table$converged[j] <- isTRUE(z$converged)
    diagnostics_table$criterion[j] <- z$log_marginal_likelihood
    diagnostics_table$fit_seconds[j] <- z$fit_seconds
    diagnostics_table$outer_convergence[j] <- if (isTRUE(z$converged)) "converged" else "failed"
    constraint_residual[j, ] <- z$constraint_residual
    observation_spatial_mean[j, ] <- z$observation_spatial_mean
    if (diagnostics) {
      inla_diagnostics[[j]] <- z[c(
        "tau", "tau_internal", "precision_scale", "lambda", "mode_status", "mode_status_text",
        "constraint_residual", "constraint_residual_uncorrected",
        "observation_spatial_mean", "estimation"
      )]
    }
    total_offset_j <- model$offset
    if (!is.null(offset)) {
      total_offset_j <- total_offset_j + if (is.matrix(offset)) offset[j, ] else offset
    }
    marginal_t0 <- proc.time()[["elapsed"]]
    faux <- tryCatch(.inlast_as_gam(z, model, Y[j, ], total_offset_j),
                     error = function(e) e)
    marginal <- if (inherits(faux, "condition")) faux else tryCatch(
        .mgcvst_marginal_score(
          faux, marginal_test, marginal_args,
          test_component = unname(model$geometry$target[["global"]])
        ), error = function(e) e)
    marginal_seconds[j] <- proc.time()[["elapsed"]] - marginal_t0
    if (inherits(marginal, "condition")) {
      diagnostics_table$error_class[j] <- class(marginal)[1L]
      diagnostics_table$error_message[j] <- conditionMessage(marginal)
      diagnostics_table$error_call[j] <- paste(deparse(conditionCall(marginal)), collapse = " ")
    } else {
      diagnostics_table$marginal_p_value[j] <- marginal$p_value
      diagnostics_table$marginal_requested_method[j] <- marginal$requested_method
      diagnostics_table$marginal_method[j] <- marginal$method
      diagnostics_table$marginal_fallback[j] <- marginal$fallback
    }
    if (retain_marginal) {
      captured <- if (inherits(faux, "condition")) faux else tryCatch(
        .mgcvst_capture_marginal(
          faux, marginal_geometry,
          test_component = unname(model$geometry$target[["global"]])
        ),
        error = function(e) e)
      if (inherits(captured, "condition")) {
        marginal_state[[j]] <- captured
      } else {
        if (is.null(marginal_geometry)) marginal_geometry <- captured$geometry
        marginal_state[[j]] <- captured$state
      }
    }
    if (!is.null(coefficient)) {
      for (name in names(coefficient)) coefficient[[name]][j, ] <- z$coefficients[[name]]
    }
  }
  target_sp <- vapply(model$geometry$target, function(j) {
    model$geometry$smooth[[j]]$sp_index
  }, integer(1L))
  component_lambda <- smoothing_parameters[, target_sp, drop = FALSE]
  colnames(component_lambda) <- names(target_sp)
  total_offset <- if (is.null(offset)) model$offset else if (is.matrix(offset))
    sweep(offset, 2L, model$offset, "+") else model$offset + offset
  elapsed <- proc.time()[["elapsed"]] - t0
  first_good <- which(!vapply(fits, inherits, logical(1L), what = "condition"))[1L]
  estimation <- if (is.na(first_good)) NULL else fits[[first_good]]$estimation
  if (!is.null(estimation)) estimation$control <- control

  ans <- structure(list(
    feature_id = feature_id,
    working_error = E,
    working_variance = V,
    dispersion = dispersion,
    lambda = component_lambda[, "global"],
    component_lambda = component_lambda,
    smoothing_parameters = smoothing_parameters,
    family_parameters = family_parameters,
    geometry = model$geometry,
    nuisance_covariance = nuisance_covariance,
    expected_nuisance_covariance = expected_nuisance_covariance,
    row_id = model$geometry$row_id,
    offset = total_offset,
    linear_design = model$geometry$X,
    score_components = model$geometry$score_components,
    model_setting = model$setting,
    model = model,
    diagnostics = diagnostics_table,
    timing = list(elapsed = elapsed, fit_elapsed = fit_elapsed,
                  marginal_elapsed = sum(marginal_seconds),
                  compaction_elapsed = max(0, elapsed - fit_elapsed - sum(marginal_seconds)),
                  workers = workers, chunks = length(groups),
                  chunk_size = chunk_size, backend = class(BPPARAM)[1L]),
    smooth_coefficients = coefficient,
    retain_smooth = retain_smooth,
    test_engine = if (length(model$components) == 1L) "single_model" else NULL,
    score_backend = score_backend,
    score_backend_requested = score_backend_requested,
    score_sparse = score_sparse,
    estimator = "INLA",
    estimation = estimation,
    constraint_residual = constraint_residual,
    observation_spatial_mean = observation_spatial_mean,
    inla_diagnostics = inla_diagnostics,
    mean_constraint = "observation",
    mean_constraint_active = TRUE,
    call = match.call()
  ), class = c("inlaST_fit", "mgcvST_model_fit", "mgcvST_fit", "mgcvST"))
  if (retain_marginal) {
    ans$marginal_data <- list(
      version = 1L, geometry = if (is.null(marginal_geometry)) list() else list(marginal_geometry),
      geometry_index = if (is.null(marginal_geometry)) integer(p) else rep.int(1L, p),
      state = marginal_state,
      definition = "frozen_INLA_conditional_marginal_TAPS"
    )
  }
  ans
}
