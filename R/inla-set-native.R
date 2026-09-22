# Native sparse setup for the INLA path.
#
# The builder assembles the INLA latent-model specification directly from
# (mesh, coordinates, kappa). It retains the sparse projector `A`, mesh
# precision `Q`, constraint vector `g`, and the small dense nuisance design.
# The mesh dimension is detected from the mesh object, so 2D triangulations and
# 3D tetrahedralisations use the same model contract.
#
# Any number of supported nuisance smooths are carried as iid blocks alongside
# parametric covariates. Categorical `bs = "re"` terms use indicator designs;
# full-rank `bs = "gp"` terms are centred and whitened before entering INLA.
# Their penalized coefficients share the small Vp block with the unpenalized
# fixed coefficients.

# ---------------------------------------------------------------------------
# Statistical contract of the native model
# ---------------------------------------------------------------------------
#
# Latent field.  `u` is the value of the Matern SPDE field at the `q` mesh
# nodes.  Its prior is `u ~ N(0, (tau * Q(kappa))^-1)` with
#
#     Q(kappa) = kappa^4 * M0 + 2 * kappa^2 * M1 + M2,
#
# the standard Lindgren-Rue-Lindstrom finite-element approximation of the
# operator `(kappa^2 - Laplacian)^(alpha/2)` at `alpha = 2`.  `M0` is the lumped
# mass matrix, `M1` the stiffness matrix and `M2 = M1 M0^-1 M1`.  `kappa` is
# fixed by the caller; the single precision multiplier `tau` is the smoothing
# parameter estimated by INLA (`generic0` with `Cmatrix = Q`).
#
# Smoothness and range.  With `alpha = 2` the Matern smoothness is
# `nu = alpha - d/2`, so `nu = 1` on a 2D mesh and `nu = 1/2` on a 3D mesh
# (the 3D field is exponential-correlated, i.e. rougher for the same `kappa`).
# The usual empirical range at correlation 0.1 is `sqrt(8 * nu) / kappa`, hence
# `2 * sqrt(2) / kappa` in 2D and `2 / kappa` in 3D.  Because `nu` differs by
# dimension, the same `kappa` does NOT mean the same field in 2D and 3D; the
# per-dimension `nu` and `range` are recorded on the returned model.
#
# Observation-mean constraint.  The linear predictor is
#
#     eta = offset + X beta + A u,
#
# where `A` is the sparse barycentric projector from mesh nodes to observation
# locations.  Each row of `A` sums to one (barycentric weights), so the constant
# mesh vector `1` maps to the constant observation vector `1`: the intercept
# already in `X` lies exactly in `col(A)`.  Left alone, the intercept and the
# constant component of `u` are confounded, and the confounding is not benign,
# because the prior does penalise that direction: for `u = c * 1`,
#
#     u' (tau Q) u = tau c^2 kappa^4 * 1' M0 1,
#
# which is finite and proportional to `kappa^4`.  The split between the
# intercept and the field mean would therefore be decided by `kappa`, not by the
# data.  The native path removes the direction exactly, with the SAME constraint
# the legacy path uses:
#
#     g' u = 0,     g = crossprod(A, rep(1/n, n)) = colMeans(A).
#
# In words: the fitted smooth is column-centred over the observed locations,
# `mean(A u) = 0`.  This is the observation-mean constraint, deliberately NOT
# the integral (`1' M0 u = 0`) constraint: it is the one that makes `beta0` the
# mean of the linear predictor at the data, and it is what the score test's
# null covariance conditions on.  It is passed to INLA as a hard linear
# constraint (`extraconstr = list(A = matrix(g, 1), e = 0)`) and the `generic0`
# `rankdef` is incremented by one so that `tau` is normalised on the
# `q - 1`-dimensional constrained support.
#
# Implied design.  The legacy path realised the constraint by an explicit
# orthonormal null-space basis `Z` (`q` by `q - 1`, `g' Z = 0`) and then formed
# the dense `B = A Z`.  `B` has exactly zero column means, which is the sense in
# which "column-centre the smooth design" describes the constraint.  The native
# path keeps `(A, Q, g)` and lets INLA impose `g' u = 0`; the two
# parameterisations give the same fitted field `A u` and the same marginal
# likelihood, and the sparse score kernel consumes `(A, Q, g)` directly, so `Z`
# is never needed.  `projection` is consequently `NULL` on native blocks, and
# `result$coefficients` are the `q` mesh-node values (satisfying `g' u = 0`)
# rather than the `q - 1` projected coordinates.

# Normalise the supported mesh containers to (dimension, vertices, cells,
# coordinate transform).  `spde_mesh` carries the unit-width rescaling used when
# the mesh was built, so observation coordinates must be mapped through it.
.inlast_native_mesh <- function(mesh) {
  if (inherits(mesh, "spde_mesh")) {
    core <- mesh$mesh
    transform <- mesh$transform
  } else {
    core <- mesh
    transform <- list(center = NULL, scale = 1)
  }
  if (inherits(core, "fm_mesh_3d")) {
    loc <- as.matrix(core$loc[, 1:3, drop = FALSE])
    cells <- as.matrix(core$graph$tv)
    dim <- 3L
    if (ncol(cells) != 4L) {
      stop("A 3D SPDE mesh must supply four-column tetrahedra.")
    }
  } else {
    flat <- .spde_basis_mesh(core)
    loc <- flat$xy
    cells <- flat$tv
    dim <- 2L
    if (is.null(transform$center)) transform <- flat$transform
  }
  storage.mode(loc) <- "double"
  storage.mode(cells) <- "integer"
  if (any(!is.finite(loc))) stop("The mesh vertices must be finite.")
  if (anyNA(cells) || any(cells < 1L) || any(cells > nrow(loc))) {
    stop("The mesh cells must be valid one-based vertex indices.")
  }
  center <- transform$center
  if (is.null(center)) center <- rep(0, dim)
  scale <- transform$scale
  if (is.null(scale)) scale <- 1
  if (length(center) != dim || !all(is.finite(center)) ||
      length(scale) != 1L || !is.finite(scale) || scale <= 0) {
    stop("The mesh coordinate transform is incompatible with its dimension.")
  }
  list(dim = dim, q = nrow(loc), loc = loc, cells = cells, core = core,
       center = as.numeric(center), scale = as.numeric(scale))
}

# Finite-element matrices.  In 2D the package's own linear-triangle assembler is
# used, so the native and legacy precisions agree to the last bit.  In 3D the
# same quantities come from fmesher (`c0` lumped mass, `g1` stiffness,
# `g2 = g1 c0^-1 g1`), which is the reference implementation for tetrahedra.
.inlast_native_fem <- function(m) {
  if (m$dim == 2L) {
    fem <- .spde_basis_fem(list(xy = m$loc, tv = m$cells))
    return(list(M0 = fem$M0, M1 = fem$M1, M2 = fem$M2, source = "mgcvST_fem2d"))
  }
  if (!requireNamespace("fmesher", quietly = TRUE)) {
    stop("A 3D native SPDE model requires the fmesher package.")
  }
  fem <- fmesher::fm_fem(m$core, order = 2L)
  if (is.null(fem$c0) || is.null(fem$g1) || is.null(fem$g2)) {
    stop("fmesher::fm_fem() did not return c0, g1 and g2 for the 3D mesh.")
  }
  list(M0 = Matrix::forceSymmetric(fem$c0),
       M1 = Matrix::forceSymmetric(fem$g1),
       M2 = Matrix::forceSymmetric(fem$g2), source = "fmesher_fm_fem")
}

# Sparse observation-to-node projector.  Rows sum to one in both dimensions.
.inlast_native_projector <- function(m, loc) {
  if (m$dim == 2L) {
    return(.spde_basis_project(list(xy = m$loc, tv = m$cells), loc))
  }
  if (!requireNamespace("fmesher", quietly = TRUE)) {
    stop("A 3D native SPDE model requires the fmesher package.")
  }
  A <- fmesher::fm_basis(m$core, loc = loc)
  if (is.list(A) && !is.null(A$A)) A <- A$A
  A <- methods::as(methods::as(A, "generalMatrix"), "CsparseMatrix")
  if (nrow(A) != nrow(loc) || ncol(A) != m$q) {
    stop("The 3D mesh projector has incompatible dimensions.")
  }
  rows <- Matrix::rowSums(A)
  if (any(!is.finite(rows)) || max(abs(rows - 1)) > 1e-8) {
    stop("Some observation locations fall outside the supplied 3D mesh.")
  }
  A
}

# Convert a raw mgcv GP smooth into an observation-centred iid block.  If
# `S = R'R`, the coefficient change b = R gamma gives penalty b'b and design
# `B R^-1`.  The QR complement is the same observation-intercept projection
# used by the SPDE bridge, applied here before checking that the remaining
# penalty is full rank.
.inlast_native_gp_bridge <- function(spec, data, label) {
  smooths <- mgcv::smoothCon(
    spec, data, absorb.cons = FALSE, scale.penalty = TRUE
  )
  lapply(smooths, function(sm) {
    if (length(sm$S) != 1L) {
      stop("The INLA nuisance term '", label,
           "' must supply one full-rank penalty; only full-rank-penalty ",
           "smooths are supported on the INLA path.")
    }
    B <- as.matrix(sm$X)
    S <- as.matrix(sm$S[[1L]])
    g <- colMeans(B)
    qg <- qr(matrix(g, ncol = 1L))
    if (qg$rank != 1L || ncol(B) < 2L) {
      stop("The INLA nuisance term '", label,
           "' does not have an intercept direction to project out.")
    }
    K <- qr.Q(qg, complete = TRUE)[, -1L, drop = FALSE]
    B <- B %*% K
    S <- crossprod(K, S %*% K)
    S <- (S + t(S)) / 2
    if (qr(S)$rank != ncol(S)) {
      stop("The INLA nuisance term '", label, "' retains an unpenalized ",
           "null-space direction after intercept projection; only ",
           "full-rank-penalty smooths are supported on the INLA path.")
    }
    R <- chol(S)
    Z <- t(backsolve(R, t(B), transpose = TRUE))
    colnames(Z) <- paste0(sm$label, ".", seq_len(ncol(Z)))
    list(
      type = "gp", name = sm$label, label = label, Z = Z,
      projection = K, centred_penalty = S, penalty_cholesky = R,
      null_space_dim = 0L
    )
  })
}

# The spatial term is supplied by mesh/kappa/coordinates. Separate every
# supported nuisance smooth from the parametric formula while preserving its
# offset and retaining every smooth variable in the common na.fail frame.
.inlast_native_design <- function(formula, data) {
  if (!inherits(formula, "formula") || length(formula) != 3L ||
      !is.symbol(formula[[2L]])) {
    stop("Supply a two-sided formula with a single response name.")
  }
  data <- as.data.frame(data)
  if (!nrow(data)) stop("data must be one shared non-empty data frame.")
  response <- as.character(formula[[2L]])
  if (response %in% all.vars(formula[[3L]])) {
    stop("The response cannot also be a covariate or offset.")
  }
  # The response is a per-feature bridge supplied later by inlaST.estimate();
  # mirror the legacy setup and install a zero placeholder when absent.
  if (!response %in% names(data)) data[[response]] <- numeric(nrow(data))
  terms <- stats::terms(formula, data = data)
  labels <- attr(terms, "term.labels")
  smooth <- grepl("^(s|te|ti|t2)\\(", labels)
  split <- mgcv::interpret.gam(formula)
  mf <- stats::model.frame(
    split$fake.formula, data, na.action = stats::na.fail
  )
  nuisance <- list()
  if (any(smooth)) {
    spatial <- grepl("bs\\s*=\\s*[\"']spde", labels[smooth])
    if (any(spatial)) {
      stop("The native INLA path takes its single spatial term from mesh and ",
           "kappa; remove the spatial smooth term(s) from the formula: ",
           paste(labels[smooth][spatial], collapse = ", "), ".")
    }
    smooth_labels <- labels[smooth]
    specs <- split$smooth.spec
    if (length(specs) != length(smooth_labels)) {
      stop("The native INLA path could not align the nuisance smooth terms.")
    }
    for (j in seq_along(specs)) {
      expr <- str2lang(smooth_labels[j])
      args <- as.list(expr)[-1L]
      arg_names <- names(args)
      if (is.null(arg_names)) arg_names <- rep("", length(args))
      variables <- which(!nzchar(arg_names))
      bs <- which(arg_names == "bs")
      if (!identical(expr[[1L]], as.name("s")) || length(bs) != 1L ||
          !is.character(args[[bs]]) || length(args[[bs]]) != 1L) {
        stop("The native INLA nuisance term must be either ",
             "s(group, bs = 're') or s(..., bs = 'gp'). Offending term: ",
             smooth_labels[j], ".")
      }
      basis <- as.character(args[[bs]])
      if (identical(basis, "re")) {
        if (length(variables) != 1L || !is.symbol(args[[variables]])) {
          stop("The native INLA random-effect term must be ",
               "s(group, bs = 're') with one grouping column.")
        }
        if (any(arg_names == "by")) {
          stop("The native INLA random-effect term does not support by=.")
        }
        group_name <- as.character(args[[variables]])
        nuisance[[length(nuisance) + 1L]] <- list(
          type = "re", name = group_name, data_name = group_name,
          label = smooth_labels[j]
        )
      } else if (identical(basis, "gp")) {
        if (any(arg_names == "by")) {
          stop("The native INLA GP term does not support by=.")
        }
        if (isTRUE(specs[[j]]$fixed)) {
          stop("The INLA nuisance term '", smooth_labels[j],
               "' must supply one full-rank penalty; only ",
               "full-rank-penalty smooths are supported on the INLA path.")
        }
        nuisance <- c(
          nuisance,
          .inlast_native_gp_bridge(specs[[j]], mf, smooth_labels[j])
        )
      } else {
        stop("The native INLA nuisance term must use bs = 're' or bs = 'gp'. ",
             "Offending term: ", smooth_labels[j], ".")
      }
    }
    internal_names <- make.unique(
      vapply(nuisance, `[[`, character(1L), "name"), sep = "."
    )
    for (j in seq_along(nuisance)) nuisance[[j]]$name <- internal_names[j]
  }
  y <- as.numeric(stats::model.response(mf))
  offset <- stats::model.offset(mf)
  offset <- if (is.null(offset)) numeric(nrow(mf)) else as.numeric(offset)
  for (j in seq_along(nuisance)) {
    if (!identical(nuisance[[j]]$type, "re")) next
    data_name <- nuisance[[j]]$data_name
    group <- mf[[data_name]]
    if (!(is.factor(group) || is.character(group))) {
      stop("The native INLA random-effect grouping variable '", data_name,
           "' must be a factor or character column.")
    }
    cls <- droplevels(factor(group))
    Z <- Matrix::sparseMatrix(i = seq_len(nrow(mf)), j = as.integer(cls),
                              x = 1, dims = c(nrow(mf), nlevels(cls)))
    colnames(Z) <- paste0(data_name, ":", levels(cls))
    nuisance[[j]]$levels <- levels(cls)
    nuisance[[j]]$Z <- Z
  }
  X <- stats::model.matrix(stats::terms(split$pf), mf)
  X <- as.matrix(X)
  storage.mode(X) <- "double"
  if (nrow(X) != nrow(mf) || any(!is.finite(X))) {
    stop("The parametric design must be finite with one row per observation.")
  }
  if (any(!is.finite(offset)) || any(!is.finite(y))) {
    stop("The response and offset must be finite.")
  }
  list(y = y, offset = offset, X = X, response = response,
       nuisance = nuisance,
       row_id = as.character(rownames(mf)))
}

# Observation coordinates, mapped through the mesh's own coordinate transform.
.inlast_native_coordinates <- function(data, coordinates, m) {
  if (!is.character(coordinates) || length(coordinates) != m$dim ||
      anyDuplicated(coordinates) || !all(coordinates %in% names(data))) {
    stop("coordinates must name ", m$dim,
         " distinct columns of data for this mesh dimension.")
  }
  loc <- as.matrix(data[, coordinates, drop = FALSE])
  storage.mode(loc) <- "double"
  if (any(!is.finite(loc))) {
    stop("The observation coordinates must be finite.")
  }
  sweep(loc, 2L, m$center, "-") / m$scale
}

# Matern smoothness and practical range for the detected mesh dimension.
.inlast_native_matern <- function(kappa, dim, alpha = 2) {
  nu <- alpha - dim / 2
  if (nu <= 0) stop("alpha = 2 does not give a valid Matern field in dimension ", dim, ".")
  list(alpha = alpha, dim = dim, kappa = kappa, nu = nu,
       range = sqrt(8 * nu) / kappa,
       range_definition = "sqrt(8 nu)/kappa, correlation ~0.1")
}

# Width of a target block's retained coefficient vector.  Legacy blocks store
# the dense projected basis; native blocks store only its column count.
.inlast_target_width <- function(model, j) {
  sm <- model$geometry$smooth[[j]]
  if (!is.null(sm$B)) return(ncol(sm$B))
  width <- sm[["coef_dim", exact = TRUE]]
  if (is.null(width) || length(width) != 1L || !is.finite(width)) {
    stop("The target smooth has no retained coefficient dimension.")
  }
  as.integer(width)
}

# Assemble the native model.  The returned object carries the very same
# `inla_spec` contract that `.inlast_model_spec()` produces for legacy models,
# so `inlaST.estimate()`, the sparse score kernel and the Liu marginal all run
# unchanged. The spatial projector stays sparse. Only the bounded nuisance
# design (fixed columns and nuisance iid blocks) is carried densely.
.inlast_set_native <- function(formula, data, family, mesh, kappa, coordinates,
                               setting, precision_scale, control) {
  if (!identical(setting, "global")) {
    stop("The native INLA path supports setting = 'global' with one spatial field.")
  }
  if (!identical(precision_scale, "raw")) {
    stop("The native INLA path supports precision_scale = 'raw' only; the ",
         "observation-scale normaliser needs the dense projected basis.")
  }
  if (length(kappa) != 1L || !is.numeric(kappa) || !is.finite(kappa) || kappa <= 0) {
    stop("kappa must be one positive finite number.")
  }
  kappa <- as.numeric(kappa)
  m <- .inlast_native_mesh(mesh)
  design <- .inlast_native_design(formula, data)
  .inlast_check_nuisance_width(ncol(design$X) +
    sum(vapply(design$nuisance, function(z) ncol(z$Z), integer(1L))))
  loc <- .inlast_native_coordinates(as.data.frame(data), coordinates, m)
  n <- length(design$y)
  if (nrow(loc) != n) stop("The coordinates and the model frame disagree in length.")

  A <- .inlast_native_projector(m, loc)
  A <- methods::as(methods::as(A, "generalMatrix"), "CsparseMatrix")
  fem <- .inlast_native_fem(m)
  Q <- kappa^4 * fem$M0 + 2 * kappa^2 * fem$M1 + fem$M2
  Q <- methods::as(Matrix::forceSymmetric(Q), "symmetricMatrix")
  if (!inherits(Q, "dsCMatrix")) Q <- methods::as(Q, "dsCMatrix")
  if (nrow(Q) != m$q || ncol(A) != m$q || any(!is.finite(Q@x)) ||
      any(!is.finite(A@x))) {
    stop("The assembled SPDE precision and projector are incompatible.")
  }

  # The observation-mean constraint; see the contract note above.
  g <- as.numeric(Matrix::crossprod(A, rep(1 / n, n)))
  if (length(g) != m$q || any(!is.finite(g)) || !any(abs(g) > 0)) {
    stop("The observation mean constraint has invalid rank.")
  }

  label <- paste0("s(", paste(coordinates, collapse = ","), ")")
  X <- design$X
  p_x <- ncol(X)
  fixed_names <- colnames(X)
  if (!p_x) fixed_names <- character()
  nuisance_map <- lapply(seq_len(p_x), function(k) {
    list(source = "fixed", block = NA_integer_, index = k, full_column = k)
  })
  random <- list(list(
    name = "global", A = A, Q = Q, kind = "spde", target = TRUE,
    constraint = g, projection = NULL, geometry_index = 1L,
    sp_index = 1L, rankdef = 0L, precision_scale = 1
  ))
  U <- X
  nuisance_index <- seq_len(p_x)
  sp_names <- label
  full_start <- p_x + m$q
  for (k in seq_along(design$nuisance)) {
    block <- design$nuisance[[k]]
    width <- ncol(block$Z)
    j <- length(random) + 1L
    random[[j]] <- list(
      name = block$name, A = block$Z, Q = Matrix::Diagonal(width),
      kind = "nuisance", subtype = "iid", target = FALSE,
      constraint = NULL, projection = NULL, geometry_index = NA_integer_,
      sp_index = j, rankdef = 0L, precision_scale = 1,
      nuisance_type = block$type, levels = block$levels
    )
    U <- cbind(U, as.matrix(block$Z))
    block_index <- full_start + seq_len(width)
    nuisance_index <- c(nuisance_index, block_index)
    nuisance_map <- c(nuisance_map, lapply(seq_len(width), function(k) {
      list(source = "random", block = j, index = k,
           full_column = block_index[k])
    }))
    sp_names <- c(sp_names, block$label)
    full_start <- full_start + width
  }
  spec <- list(
    n = n,
    family = .inlast_family(family),
    fixed = list(X = X, names = fixed_names),
    random = random,
    nuisance_design = U,
    nuisance_map = nuisance_map,
    nuisance_index = nuisance_index,
    geometry_sp_length = length(random),
    sp_names = sp_names,
    offset = design$offset,
    mean_constraint = "observation",
    mean_constraint_active = TRUE
  )
  if (spec$family == "negative_binomial" && isTRUE(family$n.theta == 0)) {
    spec$nb_size_fixed <- as.numeric(family$getTheta(trans = TRUE))
  }
  if (!p_x) {
    stop("The native INLA path requires at least one parametric column; keep the intercept.")
  }

  geometry <- list(
    X = X,
    smooth = list(list(
      label = label, B = NULL, coef_dim = m$q, penalties = list(),
      sp_index = 1L, fixed = FALSE, score_component = "global",
      columns = p_x + seq_len(m$q)
    )),
    target = stats::setNames(1L, "global"),
    score_components = "global",
    offset = design$offset,
    row_id = design$row_id,
    sp = rep(NA_real_, length(random)),
    nuisance_columns = nuisance_index,
    nuisance_design = U,
    nuisance_projection = "expected_Fisher_penalized_Vp"
  )

  model <- list(
    geometry = geometry,
    setting = "global",
    components = "global",
    response = design$response,
    y = design$y,
    offset = design$offset,
    kappa = stats::setNames(kappa, "global"),
    family = family,
    inla_spec = spec,
    spde = c(.inlast_native_matern(kappa, m$dim),
             list(mesh_vertices = m$q, mesh_cells = nrow(m$cells),
                  fem_source = fem$source,
                  constraint = "observation mean: colMeans(A) %*% u = 0")),
    mesh = mesh,
    precision_scale = "raw",
    mean_constraint = "observation",
    mean_constraint_active = TRUE,
    native = TRUE,
    timing = list(setup_seconds = NA_real_, elapsed = NA_real_)
  )
  model$inla_spec$precision_scale_mode <- "raw"
  model$inla_spec$spatial_precision_scale <- stats::setNames(1, "global")
  model$inla_control <- .inlast_family_control(model, control)
  model$score_backend <- "sparse"
  class(model) <- c("inlaST_native_model", "inlaST_model", "mgcvST_model")
  capability <- .inlast_sparse_score_capability(model)
  if (!capability$eligible) {
    stop("The native INLA model is not sparse-score eligible: ", capability$reason, ".")
  }
  model
}
