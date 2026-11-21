# Build the observation-centred projected basis used by the existing score
# engine and the raw sparse SPDE matrices used by INLA. The basis keeps the
# kappa_internal fixed by spde_basis() from its own observation coordinates.
.inlast_prepare_basis <- function(basis, coordinates) {
  .spde_basis_validate(basis)
  mesh <- list(xy = basis$mesh_vertices, tv = basis$mesh_triangles)
  coordinates <- as.matrix(coordinates)
  storage.mode(coordinates) <- "double"
  if (ncol(coordinates) != 2L || !nrow(coordinates) || any(!is.finite(coordinates))) {
    stop("The model coordinates must be a finite two-column matrix.")
  }
  loc <- sweep(coordinates, 2L, basis$transform$center, "-") /
    basis$transform$scale
  A <- .spde_basis_project(mesh, loc)
  Q <- .spde_fem_precision(.spde_basis_fem(mesh), basis$kappa_internal)

  # This is the observation mean, not a mesh-node sum-to-zero constraint.
  g <- as.numeric(Matrix::crossprod(A, rep(1 / nrow(A), nrow(A))))
  qg <- qr(matrix(g, ncol = 1L))
  if (qg$rank != 1L || length(g) < 2L) {
    stop("The observation mean constraint has invalid rank.")
  }
  Z <- qr.Q(qg, complete = TRUE)[, -1L, drop = FALSE]
  projected_Q <- crossprod(Z, as.matrix(Q %*% Z))

  constrained <- basis
  constrained$coordinates <- coordinates
  constrained$coordinate_keys <- .spde_coordinate_keys(coordinates)
  constrained$B <- as.matrix(A %*% Z)
  constrained$Q <- (projected_Q + t(projected_Q)) / 2
  constrained$projection <- Z
  constrained$projection_rank <- 1L
  constrained$project_intercept <- TRUE
  constrained$raw_dimension <- ncol(A)

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
  if (identical(label, "poisson") && identical(link, "log")) return("poisson")
  if (grepl("^negative binomial", label) && identical(link, "log")) {
    return("negative_binomial")
  }
  stop("inlaST.set() supports Poisson(log) and negative-binomial(log) families.")
}

# Soft guard on the width of the nuisance design. The nuisance design is carried
# dense through the score kernel's exact GLS Vp / P construction, so an
# accidentally huge one (a factor with hundreds of levels, a wide interaction)
# is a performance and conditioning hazard rather than a supported model.
.INLAST_MAX_NUISANCE_COLUMNS <- 1000L

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

# One message for both setup paths: the Stage 1 null is a GLM.
.inlast_reject_nuisance_smooth <- function(labels) {
  stop("The null of inlaST.set() is a GLM: enter covariates as linear terms ",
       "or factors, not smooth terms (",
       paste(as.character(labels), collapse = ", "), ").", call. = FALSE)
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
  geometry$nuisance_projection <- "expected_Fisher_penalized_Vp"
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
  # parametric columns, so this is `geometry$X` in the mgcv setup column order.
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

#' Set the model of the spdeST pipeline
#'
#' Fixes the design, the family and the fixed-kappa SPDE field of the model,
#' with its raw sparse projector and mesh precision. The spatial field is
#' always constrained to have zero mean at the observed locations.
#'
#' @details A complete formula (with its `s(..., bs = "spde")` term) or a
#' frozen `G` design (`mgcv::gam(..., fit = FALSE)`) is parsed by mgcv.
#' Alternatively, supply the SPDE component separately through `basis`, as in
#' [model.set()], or a `mesh` for the native sparse setup. Formula/basis setup enforces
#' the observation mean-zero constraint. A supplied `G` must already have a
#' compatible constrained SPDE design; otherwise rebuild it from formula/data.
#' Supported families are Poisson and negative binomial with log link; the
#' NB size is estimated, so a family with a fixed theta is refused. Exactly one full
#' fixed-kappa SPDE basis is the tested spatial target, and it is the only
#' smooth in a frozen design. Native mesh setup also supports any number of
#' admissible iid nuisance blocks. Spatial mean-zero constraints cannot be
#' disabled. Observation weights and additional cross-penalties are not
#' supported.
#'
#' @section Nuisance smooths:
#' The native INLA setup supports one spatial target and any number of
#' categorical `s(group, bs = "re")` terms or `s(..., bs = "gp")` terms whose
#' penalty is full rank after the observation intercept is projected out.
#' Each GP design and penalty are transformed together and whitened to an iid
#' block. Every iid block retains its estimated precision penalty in the small
#' nuisance `Vp` block alongside the unpenalized fixed effects. Numeric
#' random-effect groups, `by=` nuisance terms, rank-deficient GP penalties and
#' other smooth bases are outside the native INLA model.
#'
#' Parametric covariates are supported: numeric columns, factors and
#' their interactions all enter the nuisance design as fixed effects. A second
#' spatial (`bs = "spde"`) term is rejected, and the total nuisance design is
#' capped at 1000 columns.
#'
#' @inheritParams model.set
#' @param G Optional frozen `gam.prefit` design, supplied instead of formula,
#'   data and basis. Its full SPDE terms must satisfy the observation constraint.
#' @param control `list(spde = )`, a named list of [spde.estimate()] controls
#'   stored in the model and inherited by [spdeST.estimate()].
#' @param mesh Optional mesh selecting the native sparse setup: an
#'   [spde_mesh()], an `fm_mesh_2d`, or an `fm_mesh_3d`. When supplied, the
#'   single spatial term is built directly from `mesh`, `kappa` and
#'   `coordinates`; `formula` then contains only the response, an optional
#'   `offset()`, parametric terms, and admissible nuisance terms. The spatial
#'   projector remains sparse; the bounded nuisance design is dense. See the
#'   nuisance-smooth section. Mesh dimension (2 or 3) is detected automatically.
#' @param kappa Unit-scale SPDE kappa for the native `mesh` setup; the default
#'   is `0.05`. It is fixed and never estimated, so all features share one
#'   Gaussian-process kernel shape and differ only in variance. The unit
#'   length `L` is the largest per-axis span, `max - min`, of the observation
#'   coordinates named by `coordinates`, and `kappa` is the SPDE scale in
#'   coordinates divided by `L`. The package converts it to
#'   `kappa_internal = kappa * s / L` on the mesh, where `s` is the
#'   [spde_mesh()] coordinate scale; for an `fm_mesh_2d` or `fm_mesh_3d` in
#'   raw coordinates, `s = 1` and `kappa_internal = kappa / L`. Larger `kappa`
#'   gives a more local field. With `alpha = 2` the Matern smoothness is
#'   `nu = 1` in 2D and `nu = 1/2` in 3D, and the practical range in unit
#'   lengths is `sqrt(8 * nu) / kappa`. The default `0.05` therefore gives
#'   about 57 unit lengths in 2D and 40 in 3D, a very smooth global field.
#'   `NULL` is an error. The basis, complete-formula and `G` setups take the
#'   unit-scale kappa from [spde_basis()] and do not accept this argument.
#' @param retain_smooth,diagnostics Settings of [spdeST.estimate()] fixed with
#'   the model: retain the estimated score-component coefficients; retain the
#'   per-feature record of the fit.
#' @return An `inlaST_model` for [spdeST.scan()]. It stores `kappa_unit`,
#'   `unit_length`, `coordinate_span` and `kappa_internal`, and the
#'   estimation settings in `estimation`.
#' @export
inlaST.set <- function(
    formula = NULL, data = NULL, basis = NULL, family = mgcv::nb(),
    setting = "global", coordinates = c("x", "y"), G = NULL, control = list(),
    ..., mesh = NULL, kappa = 0.05, retain_smooth = FALSE, diagnostics = FALSE) {
  estimation <- list(retain_smooth = .spdest_flag(retain_smooth, "retain_smooth"),
                     diagnostics = .spdest_flag(diagnostics, "diagnostics"))
  model <- .inlast_set_model(
    formula = formula, data = data, basis = basis, family = family,
    setting = setting, coordinates = coordinates, G = G, control = control, ...,
    mesh = mesh, kappa = kappa,
    .family_supplied = !missing(family), .kappa_supplied = !missing(kappa)
  )
  model$estimation <- estimation
  model
}

.inlast_set_model <- function(
    formula = NULL, data = NULL, basis = NULL, family = mgcv::nb(),
    setting = "global", coordinates = c("x", "y"), G = NULL, control = list(),
    ..., mesh = NULL, kappa = 0.05, .family_supplied = FALSE,
    .kappa_supplied = FALSE) {
  t0 <- proc.time()[["elapsed"]]
  kappa_supplied <- .kappa_supplied
  kappa <- .spde_kappa_check(kappa)
  if (!identical(setting, "global")) {
    stop("setting must be \"global\". The second \"local\" geographic process ",
         "(setting = \"global_local\") was removed from spdeST; supply one ",
         "spatial SPDE term.")
  }
  family_supplied <- .family_supplied
  # Accept the formula, data, family calling convention,
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
      setting = setting, control = .inlast_control(control)
    )
    native$timing <- list(setup_seconds = proc.time()[["elapsed"]] - t0,
                          elapsed = proc.time()[["elapsed"]] - t0)
    return(native)
  }
  if (kappa_supplied) {
    stop("kappa applies only to the native mesh setup. The basis, ",
         "complete-formula and G setups take the unit-scale kappa from ",
         "spde_basis().")
  }
  requested_setting <- setting
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
  }
  if (!is.null(base$G$w) && any(base$G$w != 1)) {
    stop("inlaST.set() currently requires unit observation weights.")
  }
  built <- .inlast_model_spec(base, lapply(prepared, `[[`, "raw"))
  model <- built$model
  model$inla_spec <- built$spec
  model$mean_constraint <- "observation"
  model$mean_constraint_active <- TRUE
  .spde_estimate_check(model$inla_spec)
  model$inla_control <- control
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
  if (!is.matrix(Y)) Y <- as.matrix(Y)
  if (length(dim(Y)) != 2L || !nrow(Y) || !ncol(Y)) {
    stop("Y must be a non-empty finite numeric feature-by-observation matrix.")
  }
  if (ncol(Y) != model$inla_spec$n) {
    stop("ncol(Y) must equal the number of observations in model.")
  }
  # Checked in blocks of rows: no temporary of the size of Y, no copy of a
  # double or integer matrix.
  Y <- .mgcvst_check_response_matrix(
    Y, counts = model$inla_spec$family %in% c("poisson", "negative_binomial"))
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
  chunk_size <- .mgcvst_check_chunk_size(chunk_size)
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
  list(Y = Y, feature_id = feature_id, chunk_size = chunk_size)
}

.inlast_nuisance_mode <- function(fit, spec) {
  map <- spec$nuisance_map
  width <- ncol(spec$nuisance_design)
  if (!is.list(map) || length(map) != width) {
    stop("The INLA nuisance coefficient map is incomplete.")
  }
  value <- vapply(map, function(z) {
    if (identical(z$source, "fixed")) {
      if (length(z$index) != 1L || is.na(z$index) || z$index < 1L ||
          z$index > length(fit$fixed_mode)) {
        stop("The INLA fixed nuisance coefficient map is invalid.")
      }
      return(as.numeric(fit$fixed_mode[z$index]))
    }
    if (identical(z$source, "random")) {
      if (length(z$block) != 1L || is.na(z$block) || z$block < 1L ||
          z$block > length(fit$random_mode) || length(z$index) != 1L ||
          is.na(z$index) || z$index < 1L ||
          z$index > length(fit$random_mode[[z$block]])) {
        stop("The INLA random nuisance coefficient map is invalid.")
      }
      return(as.numeric(fit$random_mode[[z$block]][z$index]))
    }
    stop("The INLA nuisance coefficient map has an unknown source.")
  }, numeric(1L))
  names(value) <- colnames(spec$nuisance_design)
  value
}
