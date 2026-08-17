# Summarize minimum-angle and area quality for an fmesher mesh.
.spde_mesh_quality <- function(mesh) {
  tv <- mesh$graph$tv
  p1 <- mesh$loc[tv[, 1], 1:2, drop = FALSE]
  p2 <- mesh$loc[tv[, 2], 1:2, drop = FALSE]
  p3 <- mesh$loc[tv[, 3], 1:2, drop = FALSE]
  a <- sqrt(rowSums((p2 - p3)^2))
  b <- sqrt(rowSums((p1 - p3)^2))
  c <- sqrt(rowSums((p1 - p2)^2))
  A <- acos(pmax(-1, pmin(1, (b^2 + c^2 - a^2) / (2 * b * c))))
  B <- acos(pmax(-1, pmin(1, (a^2 + c^2 - b^2) / (2 * a * c))))
  C <- pi - A - B
  ta <- abs((p2[, 1] - p1[, 1]) * (p3[, 2] - p1[, 2]) -
            (p2[, 2] - p1[, 2]) * (p3[, 1] - p1[, 1])) / 2
  c(min.angle = min(c(A, B, C)) * 180 / pi,
    min.area = min(ta), median.area = stats::median(ta),
    min.area.ratio = min(ta) / stats::median(ta))
}

# Convert an fmesher segment object to sf line geometries.
.spde_mesh_lines <- function(b) {
  sf::st_sfc(lapply(seq_len(nrow(b$idx)), function(j) {
    sf::st_linestring(rbind(
      b$loc[b$idx[j, 1], 1:2], b$loc[b$idx[j, 2], 1:2]
    ))
  }))
}

#' Construct a topology-aware finite element mesh
#'
#' Builds a `fmesher` triangulation independently of how the boundary was
#' obtained. If `target` is supplied, the edge length is adapted until a valid
#' mesh close to the requested triangle count is found.
#'
#' @param boundary A `spde_boundary`, `sf`, or `sfc` polygon.
#' @param loc Optional observation coordinates that every accepted mesh must
#'   cover.
#' @param target Requested approximate number of triangles.
#' @param max_edge Optional fixed maximum edge length in original coordinates.
#' @param boundary_tol Optional boundary simplification tolerance.
#' @param cutoff Optional minimum point separation.
#' @param boundary_relax Non-negative boundary expansion on the internally
#'   normalized scale. Defaults to `0.03` for a learned boundary and `0.01`
#'   otherwise.
#' @param engine Either `"fmesher"` or a regular `"grid"` starting lattice.
#' @param min_angle Requested minimum triangle angle in degrees.
#' @param min_area_ratio Minimum accepted triangle area divided by median area.
#' @param max_iter Maximum edge-length adaptation iterations.
#' @return A `spde_mesh` object containing the mesh, transformations, and
#'   diagnostics.
#' @export
spde_mesh <- function(boundary, loc = NULL, target = 100L, max_edge = NULL,
                      boundary_tol = NULL, cutoff = NULL,
                      boundary_relax = NULL,
                      engine = c("fmesher", "grid"),
                      min_angle = 21, min_area_ratio = 0.01,
                      max_iter = 10L) {
  learned <- inherits(boundary, "spde_boundary") &&
    identical(boundary$method, "learned")
  if (is.null(boundary_relax)) boundary_relax <- if (learned) 0.03 else 0.01
  domain <- .spde_domain(boundary)
  engine <- match.arg(engine)
  if (length(target) != 1L || !is.finite(target) || target < 4) {
    stop("target must be one finite number greater than or equal to 4.")
  }
  target <- as.integer(round(target))
  if (!is.null(max_edge) &&
      (length(max_edge) != 1L || !is.finite(max_edge) || max_edge <= 0)) {
    stop("max_edge must be NULL or one positive finite number.")
  }
  if (!is.null(boundary_tol) &&
      (length(boundary_tol) != 1L || !is.finite(boundary_tol) ||
       boundary_tol < 0)) {
    stop("boundary_tol must be NULL or one non-negative finite number.")
  }
  if (!is.null(cutoff) &&
      (length(cutoff) != 1L || !is.finite(cutoff) || cutoff < 0)) {
    stop("cutoff must be NULL or one non-negative finite number.")
  }
  if (length(boundary_relax) != 1L || !is.finite(boundary_relax) ||
      boundary_relax < 0) {
    stop("boundary_relax must be one non-negative finite number.")
  }
  if (length(min_area_ratio) != 1L || !is.finite(min_area_ratio) ||
      min_area_ratio < 0 || min_area_ratio >= 1) {
    stop("min_area_ratio must be in [0, 1).")
  }
  if (length(max_iter) != 1L || !is.finite(max_iter) || max_iter < 1) {
    stop("max_iter must be one positive integer.")
  }
  if (!is.null(loc)) loc <- .spde_xy(loc)

  crs <- sf::st_crs(domain)
  g0 <- sf::st_cast(sf::st_geometry(domain), "POLYGON")
  g0 <- sf::st_sfc(lapply(g0, function(x) {
    sf::st_polygon(lapply(x, .spde_close_ring))
  }), crs = sf::NA_crs_)
  valid <- sf::st_is_valid(g0, reason = TRUE)
  if (!all(valid == "Valid Geometry")) {
    stop("boundary is not a valid planar polygon: ",
         paste(unique(valid[valid != "Valid Geometry"]), collapse = "; "))
  }

  bb <- sf::st_bbox(g0)
  center <- c(mean(bb[c("xmin", "xmax")]), mean(bb[c("ymin", "ymax")]))
  scale <- max(bb[["xmax"]] - bb[["xmin"]],
               bb[["ymax"]] - bb[["ymin"]])
  if (!is.finite(scale) || scale <= 0) stop("boundary must have positive width and height.")
  g <- sf::st_sfc(lapply(g0, function(x) {
    sf::st_polygon(lapply(x, function(z) sweep(z, 2, center, "-") / scale))
  }), crs = sf::NA_crs_)
  expected.euler <- sum(vapply(g, function(x) 2 - length(x), numeric(1)))

  relax.used <- boundary_relax
  for (j.relax in seq_len(12L)) {
    g.mesh <- if (relax.used > 0) {
      sf::st_buffer(g, dist = relax.used, nQuadSegs = 2,
                    joinStyle = "MITRE", mitreLimit = 2)
    } else {
      g
    }
    mesh.euler <- sum(vapply(g.mesh, function(x) 2 - length(x), numeric(1)))
    if (mesh.euler == expected.euler) break
    relax.used <- relax.used / 2
  }
  if (mesh.euler != expected.euler) {
    relax.used <- 0
    g.mesh <- g
    mesh.euler <- expected.euler
  }

  loc.normalized <- if (is.null(loc)) NULL else
    sweep(loc, 2, center, "-") / scale
  b0 <- fmesher::fm_as_segm(g)
  area <- as.numeric(sum(sf::st_area(g.mesh)))
  h <- if (is.null(max_edge)) {
    sqrt(4 * area / (sqrt(3) * target))
  } else {
    max_edge / scale
  }

  best <- NULL
  best.error <- Inf
  history <- vector("list", as.integer(max_iter))
  for (j in seq_len(as.integer(max_iter))) {
    tol <- if (is.null(boundary_tol)) h / 5 else boundary_tol / scale
    cut <- if (is.null(cutoff)) h / 4 else cutoff / scale
    for (k in seq_len(12L)) {
      g1 <- sf::st_simplify(g.mesh, dTolerance = tol,
                            preserveTopology = TRUE)
      b1 <- fmesher::fm_as_segm(g1)
      if (engine == "fmesher") {
        mesh <- fmesher::fm_mesh_2d(
          boundary = b1, max.edge = h, cutoff = cut,
          min.angle = min_angle
        )
      } else {
        grid <- fmesher::fm_hexagon_lattice(g1, edge_len = h)
        mesh <- fmesher::fm_rcdt_2d_inla(
          loc = grid, boundary = b1, extend = FALSE, cutoff = cut,
          refine = list(max.edge = Inf, min.angle = min_angle)
        )
      }
      E <- Matrix::nnzero(mesh$graph$vv) / 2
      euler <- mesh$n - E + nrow(mesh$graph$tv)
      zero.rows <- 0L
      if (!is.null(loc.normalized)) {
        A <- INLA::inla.spde.make.A(mesh, loc = loc.normalized)
        zero.rows <- sum(Matrix::rowSums(A) == 0)
      }
      if (euler == expected.euler && zero.rows == 0L) break
      tol <- tol / 2
      cut <- cut / 2
    }

    nt <- nrow(mesh$graph$tv)
    err <- abs(nt - target) / target
    quality <- .spde_mesh_quality(mesh)
    quality.ok <- quality[["min.area.ratio"]] >= min_area_ratio
    history[[j]] <- c(
      iteration = j, max.edge.normalized = h,
      boundary.tol.normalized = tol, cutoff.normalized = cut,
      triangles = nt, relative.error = err,
      euler.characteristic = euler,
      expected.euler.characteristic = expected.euler,
      location.zero.rows = zero.rows,
      minimum.angle = quality[["min.angle"]],
      minimum.area.ratio = quality[["min.area.ratio"]],
      quality.ok = quality.ok
    )
    valid.candidate <- euler == expected.euler && zero.rows == 0L && quality.ok
    if (valid.candidate && err < best.error) {
      best <- list(mesh = mesh, boundary = b1, domain = g1,
                   max.edge = h, boundary.tol = tol, cutoff = cut,
                   location.zero.rows = zero.rows, quality = quality)
      best.error <- err
    }
    if (valid.candidate && (!is.null(max_edge) || err == 0)) break
    h <- h * min(2, max(0.5, sqrt(nt / target)))
  }

  history <- do.call(rbind, history[seq_len(j)])
  if (is.null(best)) {
    stop("No mesh preserved both boundary topology and all supplied locations.")
  }
  mesh <- best$mesh
  b1 <- best$boundary
  p0 <- sf::st_as_sf(data.frame(x = b0$loc[, 1], y = b0$loc[, 2]),
                     coords = c("x", "y"))
  p1 <- sf::st_as_sf(data.frame(x = b1$loc[, 1], y = b1$loc[, 2]),
                     coords = c("x", "y"))
  l0 <- sf::st_union(.spde_mesh_lines(b0))
  l1 <- sf::st_union(.spde_mesh_lines(b1))
  d01 <- max(as.numeric(sf::st_distance(p0, l1)))
  d10 <- max(as.numeric(sf::st_distance(p1, l0)))
  E <- Matrix::nnzero(mesh$graph$vv) / 2

  out <- list(
    mesh = mesh,
    boundary.input = boundary,
    boundary.raw = b0,
    boundary = b1,
    diagnostics = c(
      target.triangles = target,
      achieved.triangles = nrow(mesh$graph$tv),
      relative.triangle.error = best.error,
      location.zero.rows = best$location.zero.rows,
      expected.euler.characteristic = expected.euler,
      boundary.relax.requested.normalized = boundary_relax,
      boundary.relax.used.normalized = relax.used,
      boundary.relax.used.original = relax.used * scale,
      minimum.angle = best$quality[["min.angle"]],
      minimum.area = best$quality[["min.area"]],
      median.area = best$quality[["median.area"]],
      minimum.area.ratio = best$quality[["min.area.ratio"]],
      candidates.evaluated = nrow(history),
      coordinate.scale = scale,
      max.edge.normalized = best$max.edge,
      max.edge.original = best$max.edge * scale,
      boundary.tol.normalized = best$boundary.tol,
      boundary.tol.original = best$boundary.tol * scale,
      observed.boundary.error.normalized = max(d01, d10),
      observed.boundary.error.original = max(d01, d10) * scale,
      boundary.edges.raw = nrow(b0$idx),
      boundary.edges.used = nrow(b1$idx),
      vertices = mesh$n,
      triangles = nrow(mesh$graph$tv),
      euler.characteristic = mesh$n - E + nrow(mesh$graph$tv)
    ),
    transform = list(center = center, scale = scale, crs = crs),
    engine = engine,
    domain.original = sf::st_set_crs(g0, crs),
    domain.normalized = best$domain,
    history = history
  )
  class(out) <- "spde_mesh"
  out
}

#' Inspect an SPDE mesh
#'
#' @param x A `spde_mesh` object.
#' @param ... Unused.
#' @export
print.spde_mesh <- function(x, ...) {
  d <- x$diagnostics
  cat("SPDE mesh\n")
  cat("  engine: ", x$engine, "\n", sep = "")
  cat("  vertices: ", d[["vertices"]], "\n", sep = "")
  cat("  triangles: ", d[["triangles"]], " (target ",
      d[["target.triangles"]], ")\n", sep = "")
  cat("  Euler characteristic: ", d[["euler.characteristic"]], "\n", sep = "")
  invisible(x)
}

#' Plot an SPDE mesh
#'
#' @param x A `spde_mesh` object.
#' @param ... Arguments passed to the `fmesher` mesh plot method.
#' @export
plot.spde_mesh <- function(x, ...) {
  plot(x$mesh, asp = 1, ...)
  invisible(x)
}
