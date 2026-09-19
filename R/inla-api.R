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

# Soft guard on the width of the nuisance design. The nuisance design is carried
# dense through the score kernel's exact GLS Vp / P construction, so an
# accidentally huge one (a factor with hundreds of levels, a wide interaction)
# is a performance and conditioning hazard rather than a supported model.
.INLAST_MAX_NUISANCE_COLUMNS <- 200L

.inlast_check_nuisance_width <- function(p_x, detail = NULL) {
  p_x <- as.integer(p_x)
  if (length(p_x) == 1L && !is.na(p_x) && p_x > .INLAST_MAX_NUISANCE_COLUMNS) {
    stop("The nuisance design has ", p_x, " columns, above the supported ",
         "maximum of ", .INLAST_MAX_NUISANCE_COLUMNS,
         ". A high-dimensional nuisance design needs a dedicated ",
         "implementation; the INLA path carries the nuisance design densely",
         if (is.null(detail)) "." else paste0(": ", detail, "."))
  }
  invisible(p_x)
}

# One message for both setup paths, so the unsupported structure is stated
# identically wherever a caller meets it.
.INLAST_NUISANCE_SMOOTH_MESSAGE <- paste(
  "nuisance smooths are not yet supported in the INLA path; planned: an",
  "INLA-native smooth (binned rw2) fitted jointly with the spatial field and",
  "profiled out taps-style in the score kernel. Supply parametric covariates",
  "instead."
)

.inlast_reject_nuisance_smooth <- function(labels) {
  labels <- as.character(labels)
  stop("inlaST.set(): ", .INLAST_NUISANCE_SMOOTH_MESSAGE,
       if (length(labels)) paste0(" Offending term(s): ",
                                  paste(labels, collapse = ", "), ".") else "")
}

# Convert frozen mgcv geometry to an INLA latent-model specification. Keeping
# this conversion here makes the public object independently serializable.
#
# The spatial target is the single random block. Every
# remaining (parametric) column is a fixed effect, and is exactly the nuisance
# design the sparse score kernel consumes. A non-target smooth is rejected here.
.inlast_model_spec <- function(model, raw_component) {
  geometry <- model$geometry
  n <- nrow(model$L)
  target_index <- unname(geometry$target)
  nuisance_smooth <- setdiff(seq_along(geometry$smooth), target_index)
  if (length(nuisance_smooth)) {
    .inlast_reject_nuisance_smooth(
      vapply(geometry$smooth[nuisance_smooth], `[[`, character(1L), "label")
    )
  }
  tested <- sort(unique(unlist(lapply(
    geometry$smooth[geometry$target], `[[`, "columns"), use.names = FALSE
  )))
  nuisance <- setdiff(seq_len(ncol(model$L)), tested)
  geometry$nuisance_columns <- nuisance
  geometry$nuisance_design <- model$L[, nuisance, drop = FALSE]
  geometry$nuisance_projection <- "conditional_INLA_block"
  .inlast_check_nuisance_width(length(nuisance))
  model$geometry <- geometry

  random <- vector("list", length(target_index))
  for (k in seq_along(target_index)) {
    j <- target_index[k]
    sm <- geometry$smooth[[j]]
    if (sm$fixed || length(sm$penalties) != 1L || length(sm$sp_index) != 1L) {
      stop("inlaST.set() requires one fitted penalty for the target SPDE smoother; ",
           "unsupported smooth: '", sm$label, "'.")
    }
    component <- names(geometry$target)[k]
    raw <- raw_component[[component]]
    random[[k]] <- list(
      name = component, A = raw$A, Q = raw$Q, kind = "spde",
      target = TRUE, constraint = raw$constraint,
      projection = raw$projection, geometry_index = j,
      sp_index = sm$sp_index, rankdef = 0L
    )
  }

  # With non-target smooths rejected above, the nuisance columns are exactly the
  # parametric columns, so this is `geometry$X` in lpmatrix column order.
  fixed_X <- as.matrix(geometry$nuisance_design)
  fixed_names <- colnames(fixed_X)
  if (!ncol(fixed_X)) fixed_names <- character()
  nuisance_map <- lapply(seq_along(nuisance), function(k) {
    list(source = "fixed", block = NA_integer_, index = k,
         full_column = nuisance[k])
  })
  nuisance_index <- seq_along(nuisance)
  combined_design <- do.call(cbind, c(list(fixed_X), lapply(random, `[[`, "A")))
  if (length(nuisance_index) != ncol(geometry$nuisance_design) ||
      !isTRUE(all.equal(as.matrix(combined_design[, nuisance_index, drop = FALSE]),
                        geometry$nuisance_design, tolerance = 1e-10))) {
    stop("The INLA nuisance coefficient map does not match the score geometry.")
  }

  spec <- list(
    n = n,
    family = .inlast_family(model$G$family),
    fixed = list(X = fixed_X, names = fixed_names),
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

# Poisson twin of a negative-binomial specification. The family is fixed inside
# the spec, so a prescreened feature is fitted by handing the engine this spec
# instead. Everything else -- design, sparse blocks, constraint, offset -- is
# shared by reference, so the twin costs nothing beyond the family switch.
.inlast_poisson_spec <- function(spec) {
  if (!identical(spec$family, "negative_binomial")) {
    stop("The Poisson prescreen applies to negative-binomial models only.")
  }
  spec$family <- "poisson"
  spec$nb_size_fixed <- NULL
  spec
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
#' with log link, and negative binomial with log link. Exactly one full
#' fixed-kappa SPDE basis is the tested spatial target, and it is the only
#' smooth the model may contain. Spatial
#' mean-zero constraints cannot be disabled. Observation weights and
#' additional cross-penalties are not supported.
#'
#' @section Nuisance smooths:
#' The native INLA setup uses one spatial target and parametric nuisance terms.
#' A non-target `s()`, `te()`, `ti()` or `t2()` term is rejected at setup.
#'
#' Parametric covariates are supported: numeric columns, factors and
#' their interactions all enter the nuisance design as fixed effects. A second
#' spatial (`bs = "spde"`) term is rejected, and the total nuisance design is
#' capped at 200 columns.
#'
#' @inheritParams model.set
#' @param G Optional frozen `gam.prefit` design, supplied instead of formula,
#'   data and basis. Its full SPDE terms must satisfy the observation constraint.
#' @param control Named engine controls stored in the model and inherited by
#'   [inlaST.estimate()]. See that function for priors and numerical controls.
#' @param mesh Optional mesh selecting the native sparse setup: an
#'   [spde_mesh()], an `fm_mesh_2d`, or an `fm_mesh_3d`. When supplied, the
#'   single spatial term is built directly from `mesh`, `kappa` and
#'   `coordinates`; `formula` then contains only the response, an optional
#'   `offset()` and parametric terms, and no dense observation-by-coefficient
#'   basis is ever formed. Any smooth term in `formula` is rejected; see the
#'   nuisance-smooth section. Mesh dimension (2 or 3) is detected
#'   automatically.
#' @param kappa Fixed positive spatial scale required by `mesh`. With
#'   `alpha = 2` the Matern smoothness is `nu = 1` in 2D and `nu = 1/2` in 3D,
#'   and the practical range is `sqrt(8 * nu) / kappa`.
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
    setting = "global", coordinates = c("x", "y"),
    precision_scale = c("raw", "observation"), G = NULL, control = list(),
    ..., mesh = NULL, kappa = NULL) {
  t0 <- proc.time()[["elapsed"]]
  if (!identical(setting, "global")) {
    stop("setting must be \"global\". The second \"local\" geographic process ",
         "(setting = \"global_local\") was removed from mgcvST; supply one ",
         "spatial SPDE term.")
  }
  family_supplied <- !missing(family)
  # Accept the formula/data/family calling convention used by mgcvST.set(),
  # together with the prepared-basis shorthand.
  if (inherits(basis, "family") || inherits(basis, "extended.family")) {
    if (family_supplied) stop("The family was supplied twice.")
    family <- basis
    family_supplied <- TRUE
    basis <- NULL
  }
  # Native sparse setup: the spatial term comes from (mesh, kappa, coordinates)
  # and no dense observation-by-coefficient basis is ever formed.
  if (!is.null(mesh)) {
    if (!is.null(basis) || !is.null(G)) {
      stop("Supply mesh with formula and data alone; basis and G belong to the legacy setup.")
    }
    if (length(list(...))) {
      stop("The native mesh setup does not accept extra mgcv setup arguments: ",
           paste(names(list(...)), collapse = ", "), ".")
    }
    native <- .inlast_set_native(
      formula = formula, data = data, family = family, mesh = mesh,
      kappa = kappa, coordinates = coordinates,
      setting = setting, precision_scale = match.arg(precision_scale),
      control = .inlast_control(control)
    )
    native$timing <- list(setup_seconds = proc.time()[["elapsed"]] - t0,
                          elapsed = proc.time()[["elapsed"]] - t0)
    return(native)
  }
  requested_setting <- setting
  precision_scale <- match.arg(precision_scale)
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
    components <- "global"
    supplied <- list(global = basis)
    data <- as.data.frame(data)
    if (!all(coordinates %in% names(data))) {
      stop("Both coordinate columns must be present in data.")
    }
    xy <- as.matrix(data[, coordinates, drop = FALSE])
    prepared <- lapply(supplied[components], .inlast_prepare_basis,
                       coordinates = xy)
    constrained_basis <- lapply(prepared, `[[`, "basis")$global
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
  # INLA is sparse-only: there is no dense score and no fallback. Reject a
  # structure the sparse kernel cannot carry here, at setup, where the reason is
  # still attributable to the design the caller just supplied.
  capability <- .inlast_sparse_score_capability(model)
  if (!capability$eligible) {
    stop("inlaST.set() builds sparse-score models only, and this design is ",
         "not eligible: ", capability$reason, ". The dense INLA score was ",
         "removed; rebuild the design with a single global fixed-kappa SPDE ",
         "target and no extra random blocks.")
  }
  model$score_backend <- "sparse"
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

# Cluster workers start with their own default library stack: a package
# installed into a non-default library (a scratch or site library) is invisible
# to them, and -- worse -- a DIFFERENT build of the same package sitting in the
# default library would be picked up silently. Sending the manager's .libPaths()
# with every task fixes both. The wrapper lives in baseenv() so that a worker
# which cannot yet load mgcvST is still able to deserialize it; only after the
# library stack is set does it reach into the namespace for the real worker
# function. Every serialized argument (sparse Matrix blocks, plain lists and
# vectors) is likewise free of mgcvST classes.
.inlast_chunk_task <- function() {
  task <- function(index, Y, spec, base_offset, extra_offset, control,
                   diagnostics, libpaths, poisson = NULL) {
    if (length(libpaths)) .libPaths(unique(c(libpaths, .libPaths())))
    # A worker can arrive with mgcvST already loaded from another library:
    # deserialising exported globals (e.g. testthat's topLevelEnvironment
    # option under SnowParam's exportglobals) loads the namespace by name
    # BEFORE this body runs, resolving against the worker's default paths.
    # A stale build then shadows the manager's; reload from the right path.
    expected <- find.package("mgcvST", lib.loc = .libPaths(), quiet = TRUE)
    if ("mgcvST" %in% loadedNamespaces() && length(expected)) {
      current <- getNamespaceInfo(asNamespace("mgcvST"), "path")
      same <- identical(normalizePath(current, winslash = "/", mustWork = FALSE),
                        normalizePath(expected[[1L]], winslash = "/", mustWork = FALSE))
      if (!same) try(unloadNamespace("mgcvST"), silent = TRUE)
    }
    fun <- get(".inlast_fit_chunk", envir = asNamespace("mgcvST"))
    fun(index, Y, spec, base_offset, extra_offset, control, diagnostics,
        poisson)
  }
  environment(task) <- baseenv()
  task
}

# `poisson` is the full-length routing vector from the Poisson prescreen,
# indexed by the absolute feature index, or NULL when the screen is off. The
# Poisson twin of the spec is built once per chunk; it is a shallow list copy,
# so the sparse blocks are shared rather than duplicated.
.inlast_fit_chunk <- function(index, Y, spec, base_offset, extra_offset,
                              control, diagnostics, poisson = NULL) {
  route <- if (is.null(poisson)) rep(FALSE, length(index)) else
    as.logical(poisson[index])
  spec_poisson <- if (any(route)) .inlast_poisson_spec(spec) else NULL
  lapply(seq_along(index), function(k) {
    j <- index[k]
    total_offset <- base_offset
    if (!is.null(extra_offset)) {
      total_offset <- total_offset + if (is.matrix(extra_offset)) extra_offset[j, ] else extra_offset
    }
    feature_spec <- if (isTRUE(route[k])) spec_poisson else spec
    tryCatch(
      .inlast_fit_feature(feature_spec, Y[j, ], offset = total_offset,
                          control = control, diagnostics = diagnostics),
      error = function(e) e
    )
  })
}

#' Estimate mgcvST working models with sparse INLA
#'
#' Fits one latent Gaussian model per feature with INLA and returns the compact
#' working-model contract consumed by [mgcvST.test()]. No GAM fit is used.
#'
#' @param Y Numeric feature-by-observation matrix.
#' @param model An object returned by [inlaST.set()].
#' @param feature_id Unique feature identifiers.
#' @param BPPARAM A `BiocParallelParam` distributing feature chunks over
#'   workers. Each worker fits its own features with INLA using
#'   `control$num_threads` (default one), so workers multiply rather than share
#'   threads; downstream marginal work uses OpenMP `threads` in the manager.
#'   `SerialParam()` keeps everything in one process.
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
#'   `prior = "flat", param = numeric(), initial = 0`, a flat prior on the
#'   internal log precision `-log(variance)`. On the variance scale this has
#'   density proportional to `1/variance`; it is a Gaussian dispersion prior,
#'   distinct from the negative-binomial size prior. Custom normal, registered
#'   scalar INLA priors, and INLA expression/table priors remain supported.
#'   Normal parameters are mean and precision. A flat log-hyperparameter
#'   objective corresponds to density proportional to `1/parameter` on its
#'   positive scale.
#'   `control.inla` accepts supported numerical tuning, with Gaussian latent
#'   strategy and EB integration enforced. `poisson_screen_phi` (default `1.1`)
#'   is the Poisson prescreen threshold: with a negative-binomial family, each
#'   feature first gets an offset-and-covariate-only Poisson GLM, and a feature
#'   whose Pearson dispersion `phi = sum((y - mu)^2 / mu) / (n - p)` is at most
#'   the threshold is fitted with the Poisson family instead. Set it to `0` to
#'   disable the screen (`NULL` restores the default); other families ignore it.
#'   The per-feature
#'   `phi` and the family actually used are reported in the diagnostics as
#'   `prescreen_phi` and `family_used`. Unknown controls are rejected.
#' @param retain_smooth Retain estimated score-component coefficients.
#' @param diagnostics Retain per-feature INLA diagnostics and compute the
#'   expected-Fisher nuisance covariance for comparison with native `Vp`.
#'   When `FALSE`, that extra solve is skipped and its per-feature entries
#'   in `expected_nuisance_covariance` are `NULL`.
#' @param retain_marginal Retain the frozen state needed by [mgcvST.marginal()].
#' @details Hyperparameter priors remain part of INLA's empirical-Bayes
#' estimates; these are not mgcv REML estimates. `mgcv::nb(theta = value)`
#' fixes NB size, and a conflicting `control$nb_size` is rejected. The
#' working model uses conditional latent estimates and expected Fisher
#' variances. Sparse downstream scores rebuild the nuisance covariance from
#' the same expected Fisher matrix and use exact Liu trace moments;
#' this does not establish finite-sample calibration after hyperparameter
#' estimation. Every spatial component's constraint residual and observed
#' spatial mean are retained in the result.
#' The score uses the SPDE covariance conditioned on observation mean zero.
#' The sparse score uses a matching expected-curvature nuisance adjustment.
#' The INLA path is sparse-only: the sparse kernel is the sole score and
#' marginal implementation, there is no dense score, no custom marginal
#' callback and no fallback. A design the sparse capability gate rejects makes
#' [inlaST.set()] error at setup, and a model that somehow reaches this
#' function without that geometry errors here.
#' Flat hyperpriors need not yield proper hyperparameter posteriors. They are
#' supported only as empirical-Bayes optimization objectives in the current
#' single-configuration engine. A returned finite precision or zero optimizer
#' status does not establish that the maximum is interior; zero spatial
#' variance and the Poisson limit of the NB model require boundary checks.
#' @param threads OpenMP threads for sparse downstream marginal construction.
#' @return An `inlaST_fit` that is also an `mgcvST_model_fit`.
#' @export
inlaST.estimate <- function(
    Y, model, feature_id = rownames(Y),
    BPPARAM = BiocParallel::SerialParam(), chunk_size = NULL,
    offset = NULL, control = list(), retain_smooth = FALSE,
    diagnostics = FALSE, retain_marginal = FALSE, threads = 1L) {
  if (!requireNamespace("INLA", quietly = TRUE)) {
    stop("inlaST.estimate() requires the INLA package.")
  }
  if (length(threads) != 1L || !is.numeric(threads) || !is.finite(threads) ||
      threads < 1 || threads > .Machine$integer.max ||
      threads != as.integer(threads)) stop("threads must be a positive integer.")
  checked <- .inlast_validate_estimate(
    Y, model, feature_id, BPPARAM, chunk_size, offset, control,
    retain_smooth, diagnostics
  )
  # Sparse-only: the geometry builder is the capability gate. It stops with the
  # reason when the model cannot be scored; there is no dense alternative.
  score_sparse <- .inlast_sparse_score_geometry(model)
  if (!is.logical(retain_marginal) || length(retain_marginal) != 1L ||
      is.na(retain_marginal)) stop("retain_marginal must be TRUE or FALSE.")
  # Validate once in the parent so a misspelled or invalid control cannot turn
  # every feature into an otherwise opaque per-feature failure.
  control <- .inlast_control(.inlast_merge_control(model$inla_control, control))
  control <- .inlast_family_control(model, control)
  Y <- checked$Y
  feature_id <- checked$feature_id
  chunk_size <- checked$chunk_size
  # Cheap per-feature overdispersion prescreen (R/family-prescreen.R): a
  # near-Poisson feature drifts to a flat NB size likelihood, so fit it with the
  # Poisson family instead. The family is fixed inside the spec, so a routed
  # feature is fitted against a Poisson twin of the same spec
  # (.inlast_poisson_spec()); the twin is rebuilt once per chunk in the worker
  # rather than shipped, because it shares every sparse block with the original.
  screen_threshold <- .mgcvst_prescreen_threshold(
    control[["poisson_screen_phi", exact = TRUE]]
  )
  prescreen <- .mgcvst_prescreen_route(
    Y, X = model$inla_spec$fixed$X,
    offset = if (is.null(offset)) model$offset else
      if (is.matrix(offset)) sweep(offset, 2L, model$offset, "+") else
        model$offset + offset,
    threshold = screen_threshold,
    active = identical(model$inla_spec$family, "negative_binomial")
  )
  # Validate the twin once in the parent so a spec the switch cannot handle
  # fails here rather than inside every affected worker. Chunk composition is
  # deliberately left untouched: a routed feature is switched onto the twin
  # inside its chunk, so load balancing and the unscreened code path are the
  # same as before.
  if (any(prescreen$poisson)) .inlast_poisson_spec(model$inla_spec)
  family_used <- rep(model$inla_spec$family, nrow(Y))
  family_used[prescreen$poisson] <- "poisson"
  groups <- split(seq_len(nrow(Y)), ceiling(seq_len(nrow(Y)) / chunk_size))
  workers <- max(1L, min(length(groups), BiocParallel::bpworkers(BPPARAM)))
  # Feature chunks are independent latent Gaussian models. The spec is plain
  # sparse data, INLA's external binary uses per-process working directories,
  # and control$num_threads stays at its default of one inside each worker, so
  # chunk-level BiocParallel parallelism is safe.
  t0 <- proc.time()[["elapsed"]]
  chunks <- BiocParallel::bplapply(
    groups, .inlast_chunk_task(), Y = Y, spec = model$inla_spec,
    base_offset = model$offset, extra_offset = offset, control = control,
    diagnostics = diagnostics, libpaths = .libPaths(),
    poisson = if (any(prescreen$poisson)) prescreen$poisson else NULL,
    BPPARAM = BPPARAM
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
    prescreen_phi = prescreen$phi, family_used = family_used,
    stringsAsFactors = FALSE
  )
  coefficient <- if (retain_smooth) {
    lapply(model$geometry$target, function(j) {
      matrix(NA_real_, p, .inlast_target_width(model, j),
             dimnames = list(feature_id, NULL))
    })
  } else NULL
  if (!is.null(coefficient)) names(coefficient) <- names(model$geometry$target)
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
                  marginal_elapsed = 0,
                  compaction_elapsed = max(0, elapsed - fit_elapsed),
                  workers = workers, chunks = length(groups),
                  chunk_size = chunk_size, backend = class(BPPARAM)[1L]),
    smooth_coefficients = coefficient,
    retain_smooth = retain_smooth,
    test_engine = "single_model",
    score_backend = "sparse",
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
  ans$marginal_data <- list(version = 2L, definition = "sparse_INLA_marginal_TAPS_Liu")
  valid <- which(diagnostics_table$converged)
  if (length(valid)) {
    marginal_t0 <- proc.time()[["elapsed"]]
    marginal <- .inlast_marginal(ans, valid, chunk_size = min(16L, chunk_size), threads = threads)
    ans$marginal_data$result <- marginal
    ans$diagnostics$marginal_p_value[valid] <- marginal$p_value
    ans$diagnostics$marginal_requested_method[valid] <- "liu"
    ans$diagnostics$marginal_method[valid] <- "liu"
    ans$diagnostics$marginal_fallback[valid] <- FALSE
    ans$diagnostics$error_message[valid] <- marginal$error_message
    ans$timing$marginal_elapsed <- proc.time()[["elapsed"]] - marginal_t0
  }
  if (!retain_marginal) ans$marginal_data <- NULL
  ans
}
