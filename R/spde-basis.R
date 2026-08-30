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

# Barycentric projector. This is intentionally used only by spde_basis().
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

.spde_basis_validate <- function(x, loc, pc = FALSE) {
  if (!inherits(x, "mgcvST_spde_basis")) {
    stop("xt must be an object returned by spde_basis().")
  }
  if (is.null(x$B) || is.null(x$penalty) || is.null(x$coordinates)) {
    stop("xt is missing the prepared SPDE basis matrices.")
  }
  if (!identical(dim(loc), dim(x$coordinates)) ||
      !isTRUE(all.equal(unname(loc), unname(x$coordinates), tolerance = 1e-10))) {
    stop("The smooth coordinates must match the rows used by spde_basis().")
  }
  if (pc) {
    if (is.null(x$kappa)) {
      stop("bs = 'spdePC' requires a basis constructed with fixed kappa.")
    }
    cutoff <- x$pc_cutoff
    if (length(cutoff) != 1L || !is.finite(cutoff) ||
        cutoff <= 0 || cutoff > 1) {
      stop("xt$pc_cutoff must be one finite cumulative contribution in (0, 1].")
    }
    if (is.null(x$pc_values) || is.null(x$pc_vectors) ||
        is.null(x$pc_cumulative)) {
      stop("xt does not contain the prepared SPDE principal components.")
    }
  }
  invisible(x)
}

.spde_basis_warn_k <- function(object, basis) {
  if (!is.null(object$bs.dim) && object$bs.dim >= 0L) {
    warning(
      "k is ignored for bs = '", basis,
      "'; supply all basis controls through xt.", call. = FALSE
    )
  }
}

.spde_basis_component <- function(object) {
  component <- object$xt$component
  if (is.null(component)) component <- "global"
  if (!is.character(component) || length(component) != 1L ||
      is.na(component) || !(component %in% c("global", "local"))) {
    stop("xt$component must be 'global' or 'local'.")
  }
  score.component <- object$xt$score.component
  if (!is.null(score.component) &&
      (!is.character(score.component) || length(score.component) != 1L ||
       is.na(score.component) ||
       !(score.component %in% c("global", "local")))) {
    stop("xt$score.component must be NULL, 'global', or 'local'.")
  }
  list(component = component, score.component = score.component)
}

#' Prepare an SPDE basis for mgcvST fitting
#'
#' Constructs the observation projector and finite-element penalties once.
#' The returned object is self-contained: fitting with `bs = "spde"` or
#' `bs = "spdePC"` requires neither INLA, fmesher, sf, nor geometry.
#'
#' Coordinates are transformed using the scale stored by [spde_mesh()]. Thus,
#' `kappa` is dimensionless on a unit-width map. For a map scale `L`, its
#' raw-coordinate value is `kappa / L`.
#'
#' @param mesh A [spde_mesh()], `fm_mesh_2d`, or list containing `loc` and
#'   `graph$tv`.
#' @param loc Observation coordinates in the original coordinate system.
#' @param kappa Positive fixed dimensionless spatial scale. The default is
#'   `0.1`. Set explicitly to `NULL` to prepare the three penalties needed for
#'   joint estimation of kappa and tau by `bs = "spde"`.
#' @param pc_cutoff Cumulative covariance contribution retained by
#'   `bs = "spdePC"`. The default is `0.999`.
#' @param project_intercept Whether to project the intercept from the mesh
#'   coefficient space.
#' @return A self-contained object to pass directly as `xt`.
#' @export
spde_basis <- function(mesh, loc, kappa = 0.1, pc_cutoff = 0.999,
                       project_intercept = TRUE) {
  x <- .spde_basis_mesh(mesh)
  loc.raw <- .spde_xy(loc)
  loc.scaled <- sweep(loc.raw, 2L, x$transform$center, "-") /
    x$transform$scale
  if (!is.null(kappa)) {
    kappa <- as.numeric(kappa)
    if (length(kappa) != 1L || !is.finite(kappa) || kappa <= 0) {
      stop("kappa must be NULL or one positive finite number.")
    }
  }
  pc_cutoff <- as.numeric(pc_cutoff)
  if (length(pc_cutoff) != 1L || !is.finite(pc_cutoff) ||
      pc_cutoff <= 0 || pc_cutoff > 1) {
    stop("pc_cutoff must be one finite cumulative contribution in (0, 1].")
  }
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
  penalty <- list(
    CppMatrix::matrixMultiply(t(Z), CppMatrix::matrixMultiply(fem$M0, Z)),
    2 * CppMatrix::matrixMultiply(t(Z), CppMatrix::matrixMultiply(fem$M1, Z)),
    CppMatrix::matrixMultiply(t(Z), CppMatrix::matrixMultiply(fem$M2, Z))
  )
  penalty <- lapply(penalty, function(S) (S + t(S)) / 2)

  Q <- NULL
  pc.values <- NULL
  pc.vectors <- NULL
  pc.cumulative <- NULL
  if (!is.null(kappa)) {
    Q <- kappa^4 * penalty[[1L]] + kappa^2 * penalty[[2L]] + penalty[[3L]]
    Q <- (Q + t(Q)) / 2
    G <- CppMatrix::matrixSolve(Q, diag(ncol(Q)))
    E <- CppMatrix::matrixEigen((G + t(G)) / 2)
    ord <- order(E$values, decreasing = TRUE)
    pc.values <- as.numeric(E$values[ord])
    pc.vectors <- as.matrix(E$vectors[, ord, drop = FALSE])
    tol <- sqrt(.Machine$double.eps) * max(1, max(abs(pc.values)))
    if (min(pc.values) < -tol) {
      stop("The projected SPDE covariance is not positive semidefinite.")
    }
    pc.values <- pmax(pc.values, 0)
    if (any(pc.values <= 0) || !is.finite(sum(pc.values))) {
      stop("The projected SPDE covariance has invalid eigenvalues.")
    }
    pc.cumulative <- cumsum(pc.values) / sum(pc.values)
  }

  out <- list(
    B = B, penalty = penalty, Q = Q, coordinates = loc.raw,
    kappa = kappa, pc_cutoff = pc_cutoff, pc_values = pc.values,
    pc_vectors = pc.vectors, pc_cumulative = pc.cumulative,
    transform = x$transform, mesh_vertices = x$xy,
    mesh_triangles = x$tv, projection = Z,
    projection_rank = fitqr$rank, project_intercept = project_intercept,
    raw_dimension = m.raw
  )
  class(out) <- "mgcvST_spde_basis"
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
  cat("  kappa:", if (is.null(x$kappa)) "estimated" else x$kappa, "\n")
  cat("  PC cumulative contribution:", x$pc_cutoff, "\n")
  invisible(x)
}
