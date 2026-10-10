# Extract the triangulation needed by the dependency-free fitting path.
.spde_basis_mesh <- function(mesh) {
  if (inherits(mesh, "spde_mesh")) {
    core <- mesh$mesh
    transform <- mesh$transform
  } else if (inherits(mesh, "fm_mesh_2d")) {
    core <- mesh
    transform <- list(center = c(0, 0), scale = 1)
  } else if (is.list(mesh) && !is.null(mesh$loc) &&
             !is.null(mesh$graph$tv)) {
    core <- mesh
    transform <- list(center = c(0, 0), scale = 1)
  } else {
    stop("mesh must contain loc and graph$tv, or inherit from spde_mesh or fm_mesh_2d.")
  }
  xy <- as.matrix(core$loc[, 1:2, drop = FALSE])
  tv <- as.matrix(core$graph$tv)
  storage.mode(xy) <- "double"
  storage.mode(tv) <- "integer"
  if (ncol(xy) != 2L || any(!is.finite(xy))) {
    stop("mesh vertex coordinates must be a finite two-column matrix.")
  }
  if (ncol(tv) != 3L || anyNA(tv) || any(tv < 1L) || any(tv > nrow(xy))) {
    stop("mesh triangles must be a valid three-column vertex-index matrix.")
  }
  list(xy = xy, tv = tv, transform = transform)
}

# Barycentric projector for construction and new-coordinate interpolation.
# In geometry's 2D default, tsearchn already delegates to quadtree tsearch.
.spde_basis_project <- function(mesh, loc) {
  hit <- geometry::tsearchn(mesh$xy, mesh$tv, loc)
  idx <- as.integer(hit$idx)
  if (anyNA(idx)) {
    stop("Some basis locations are outside the supplied SPDE mesh.")
  }
  w <- as.matrix(hit$p)
  Matrix::sparseMatrix(
    i = rep(seq_len(nrow(loc)), each = 3L),
    j = as.vector(t(mesh$tv[idx, , drop = FALSE])),
    x = as.vector(t(w)),
    dims = c(nrow(loc), nrow(mesh$xy))
  )
}

# Linear-triangle FEM matrices. The element formulas are adapted from the
# MIT-licensed INLA inla.barrier.fem() implementation (Lindgren et al.).
.spde_basis_fem <- function(mesh) {
  n <- nrow(mesh$xy)
  nt <- nrow(mesh$tv)
  grad <- rbind(c(-1, -1), c(1, 0), c(0, 1))
  mass <- numeric(n)
  ii <- integer(9L * nt)
  jj <- integer(9L * nt)
  xx <- numeric(9L * nt)
  for (tri in seq_len(nt)) {
    px <- mesh$tv[tri, ]
    z <- mesh$xy[px, , drop = FALSE]
    T <- cbind(z[2L, ] - z[1L, ], z[3L, ] - z[1L, ])
    area <- abs(det(T)) / 2
    mass[px] <- mass[px] + area / 3
    local <- area * grad %*% solve(crossprod(T)) %*% t(grad)
    at <- (tri - 1L) * 9L + seq_len(9L)
    ii[at] <- rep(px, times = 3L)
    jj[at] <- rep(px, each = 3L)
    xx[at] <- as.vector(local)
  }
  M0 <- Matrix::Diagonal(n, mass)
  M1 <- Matrix::sparseMatrix(i = ii, j = jj, x = xx, dims = c(n, n))
  M1 <- Matrix::forceSymmetric(M1)
  M2 <- M1 %*% Matrix::Diagonal(n, 1 / mass) %*% M1
  list(M0 = M0, M1 = M1, M2 = Matrix::forceSymmetric(M2))
}

# Proper Matern SPDE precision (alpha = 2) on the mesh scale.
.spde_fem_precision <- function(fem, kappa_internal) {
  Matrix::forceSymmetric(
    kappa_internal^4 * fem$M0 + 2 * kappa_internal^2 * fem$M1 + fem$M2
  )
}

.spde_kappa_check <- function(kappa) {
  if (is.null(kappa)) {
    stop("kappa = NULL is not supported: mgcvST never estimates kappa. ",
         "Supply one positive unit-scale kappa; the default is 0.05.",
         call. = FALSE)
  }
  if (!is.numeric(kappa) || length(kappa) != 1L || !is.finite(kappa) ||
      kappa <= 0) {
    stop("kappa must be one positive finite unit-scale number ",
         "(the default is 0.05).", call. = FALSE)
  }
  as.numeric(kappa)
}

# Convert a unit-scale kappa to the mesh scale. The unit length L is the
# largest per-axis span of the observation coordinates `loc`, in raw
# coordinate units. Mesh coordinates equal raw coordinates divided by
# `mesh_scale` (the spde_mesh() transform; 1 for a raw mesh), so
# kappa_internal = kappa_unit * mesh_scale / L.
.spde_kappa_scale <- function(kappa, loc, mesh_scale = 1) {
  kappa <- .spde_kappa_check(kappa)
  loc <- as.matrix(loc)
  span <- apply(loc, 2L, function(z) max(z) - min(z))
  names(span) <- colnames(loc)
  L <- max(span)
  if (!is.finite(L) || L <= 0) {
    stop("The observation coordinates must span a positive length on at ",
         "least one axis to define the unit length for kappa.", call. = FALSE)
  }
  if (length(mesh_scale) != 1L || !is.finite(mesh_scale) || mesh_scale <= 0) {
    stop("The mesh coordinate scale must be one positive finite number.")
  }
  list(kappa_unit = kappa, unit_length = L, coordinate_span = span,
       kappa_internal = kappa * mesh_scale / L)
}

# The unit-scale kappa record carried by bases, smooths, models and fits.
.spde_kappa_fields <- function(x) {
  list(kappa_unit = x$kappa_unit, unit_length = x$unit_length,
       coordinate_span = x$coordinate_span,
       kappa_internal = x$kappa_internal)
}

.spde_basis_validate <- function(x, loc = NULL) {
  if (!inherits(x, "mgcvST_spde_basis")) {
    stop("xt must be an object returned by spde_basis().")
  }
  if (is.null(x$B) || is.null(x$Q) || is.null(x$coordinates)) {
    stop("xt is missing the prepared SPDE basis matrices.")
  }
  k <- x$kappa_internal
  if (is.null(x$kappa_unit) || length(k) != 1L || !is.finite(k) || k <= 0) {
    stop("xt has no fixed unit-scale kappa; rebuild it with spde_basis(). ",
         "mgcvST never estimates kappa.")
  }
  if (!is.null(loc) && (!is.numeric(loc) || ncol(as.matrix(loc)) != 2L ||
      any(!is.finite(loc)))) {
    stop("The smooth coordinates must be a finite two-column numeric matrix.")
  }
  invisible(x)
}

# Retained coordinate metadata for compatibility with saved basis objects.
.spde_coordinate_keys <- function(loc) {
  x <- loc[, 1L]; y <- loc[, 2L]
  x[x == 0] <- 0; y[y == 0] <- 0
  paste(sprintf("%a", x), sprintf("%a", y), sep = ":")
}

# Evaluate a fixed prepared basis without rebuilding mesh/FEM/precision.
.spde_basis_at <- function(basis, loc, timing = NULL,
                           stage = "prediction") {
  t0 <- proc.time()[["elapsed"]]
  if (!is.null(timing)) on.exit({
    key <- paste0(stage, "_seconds")
    timing[[key]] <- timing[[key]] + proc.time()[["elapsed"]] - t0
    key <- paste0(stage, "_calls")
    timing[[key]] <- timing[[key]] + 1L
  })
  .spde_basis_validate(basis)
  loc <- as.matrix(loc)
  if (!is.numeric(loc) || length(dim(loc)) != 2L || ncol(loc) != 2L ||
      any(!is.finite(loc))) stop("loc must be a finite two-column numeric matrix.")

  transform <- basis$transform
  if (length(transform$center) != 2L || any(!is.finite(transform$center)) ||
      length(transform$scale) != 1L || !is.finite(transform$scale) || transform$scale <= 0) {
    stop("The basis has an invalid saved coordinate transform.")
  }
  mesh <- list(xy = basis$mesh_vertices, tv = basis$mesh_triangles)
  P <- basis$projection
  if (is.null(mesh$xy) || is.null(mesh$tv) || is.null(P) || nrow(P) != nrow(mesh$xy)) {
    stop("The basis does not contain aligned saved mesh and projection matrices.")
  }
  if (!nrow(loc)) return(matrix(numeric(), 0L, ncol(P)))
  scaled <- sweep(loc, 2L, transform$center, "-") / transform$scale
  A <- .spde_basis_project(mesh, scaled)
  as.matrix(A %*% P)
}

.spde_basis_warn_k <- function(object, basis) {
  if (!is.null(object$bs.dim) && object$bs.dim >= 0L) {
    warning(
      "k is ignored for bs = '", basis,
      "'; supply all basis controls through xt.", call. = FALSE
    )
  }
}

# The package carries one global spatial score process.
.spde_basis_component <- function(object) {
  component <- object$xt$component
  if (is.null(component)) component <- "global"
  if (!is.character(component) || length(component) != 1L ||
      is.na(component) || !identical(component, "global")) {
    stop("xt$component must be 'global'. The second 'local' geographic ",
         "process was removed from mgcvST.")
  }
  score.component <- object$xt$score.component
  if (!is.null(score.component) && !identical(score.component, "global")) {
    stop("xt$score.component must be NULL or 'global'. The second 'local' ",
         "geographic process was removed from mgcvST.")
  }
  list(component = component, score.component = score.component)
}

#' Prepare an SPDE basis for mgcvST fitting
#'
#' Constructs the observation projector and the fixed-kappa SPDE precision
#' once. The returned object is self-contained: fitting with `bs = "spde"`
#' requires neither INLA, fmesher nor sf.
#' Basis construction and prediction use geometry for barycentric interpolation
#' on the saved mesh. No FEM or precision is recomputed. Outside-mesh
#' coordinates cause an error. Every evaluation uses the supplied coordinates,
#' including training, subset and reordered rows.
#'
#' @section Unit-scale kappa:
#' `kappa` is a unit-scale value. The unit length `L` is the largest
#' per-axis span, `max - min`, of `loc`, the observation coordinates supplied
#' to the model. The SPDE has scale `kappa` in coordinates divided by `L`, so
#' the same `kappa` gives the same field shape whether the coordinates are
#' recorded in millimetres or micrometres. The package converts it
#' internally. Mesh coordinates are the original coordinates divided by the
#' [spde_mesh()] scale `s` (`s = 1` for an `fm_mesh_2d` or a list mesh), and
#' the precision uses `kappa_internal = kappa * s / L`. The returned object
#' stores `kappa_unit`, `unit_length` (`L`), `coordinate_span` (the span per
#' axis) and `kappa_internal`.
#'
#' @param mesh A [spde_mesh()], `fm_mesh_2d`, or list containing `loc` and
#'   `graph$tv`.
#' @param loc Observation coordinates in the original coordinate system. Their
#'   bounding box defines the unit length for `kappa`.
#' @param kappa Unit-scale SPDE kappa; the default is `0.05`. It is fixed and
#'   never estimated, so every feature fitted with this basis shares one
#'   Gaussian-process kernel shape and differs only in variance. Larger
#'   `kappa` gives a more local field. With `alpha = 2`, the practical range
#'   in unit lengths is `sqrt(8 * nu) / kappa`, where `nu = 1` in 2D and
#'   `nu = 1/2` in 3D. The default `0.05` therefore gives a practical range of
#'   about 57 unit lengths in 2D and 40 in 3D, a very smooth global field.
#'   `NULL` is an error.
#' @param project_intercept Whether to project the intercept from the mesh
#'   coefficient space.
#' @return A self-contained object to pass directly as `xt`.
#' @export
spde_basis <- function(mesh, loc, kappa = 0.05, project_intercept = TRUE) {
  x <- .spde_basis_mesh(mesh)
  loc.raw <- .spde_xy(loc)
  loc.scaled <- sweep(loc.raw, 2L, x$transform$center, "-") /
    x$transform$scale
  kappa <- .spde_kappa_scale(kappa, loc.raw, x$transform$scale)
  if (!is.logical(project_intercept) || length(project_intercept) != 1L ||
      is.na(project_intercept)) {
    stop("project_intercept must be TRUE or FALSE.")
  }

  A <- .spde_basis_project(x, loc.scaled)
  fem <- .spde_basis_fem(x)
  m.raw <- ncol(A)
  if (project_intercept) {
    g <- as.matrix(Matrix::crossprod(A, matrix(1, nrow(A), 1L))) / nrow(A)
    fitqr <- qr(g)
    if (fitqr$rank < 1L || fitqr$rank >= m.raw) {
      stop("The intercept projection does not have a valid coefficient-space rank.")
    }
    Z <- qr.Q(fitqr, complete = TRUE)[,
      (fitqr$rank + 1L):m.raw, drop = FALSE]
  } else {
    fitqr <- list(rank = 0L)
    Z <- diag(m.raw)
  }
  B <- CppMatrix::matrixMultiply(as.matrix(A), Z)
  Q <- .spde_fem_precision(fem, kappa$kappa_internal)
  Q <- CppMatrix::matrixMultiply(t(Z), CppMatrix::matrixMultiply(Q, Z))
  Q <- (Q + t(Q)) / 2

  out <- c(list(B = B, Q = Q, coordinates = loc.raw), kappa, list(
    transform = x$transform, mesh_vertices = x$xy,
    mesh_triangles = x$tv, projection = Z,
    projection_rank = fitqr$rank, project_intercept = project_intercept,
    raw_dimension = m.raw
  ))
  class(out) <- "mgcvST_spde_basis"
  out$coordinate_keys <- .spde_coordinate_keys(loc.raw)
  out
}

#' @rdname spde_basis
#' @param x An object returned by [spde_basis()].
#' @param ... Unused.
#' @return `x`, invisibly.
#' @export
print.mgcvST_spde_basis <- function(x, ...) {
  cat("mgcvST SPDE basis\n")
  cat("  observations:", nrow(x$B), "\n")
  cat("  coefficients:", ncol(x$B), "\n")
  cat("  kappa (unit scale, fixed):", format(x$kappa_unit), "\n")
  cat("  unit length L:", format(x$unit_length), "\n")
  cat("  kappa (internal mesh scale):", format(x$kappa_internal), "\n")
  invisible(x)
}
