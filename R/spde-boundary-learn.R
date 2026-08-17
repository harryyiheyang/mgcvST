# Return the unique undirected edges in a triangle index matrix.
.spde_triangle_edges <- function(T) {
  E <- rbind(T[, c(1, 2), drop = FALSE],
             T[, c(1, 3), drop = FALSE],
             T[, c(2, 3), drop = FALSE])
  E <- t(apply(E, 1, sort))
  unique(E)
}

# Label connected components of an undirected graph by union-find.
.spde_graph_components <- function(n, E) {
  parent <- seq_len(n)
  root <- function(i) {
    while (parent[i] != i) {
      parent[i] <<- parent[parent[i]]
      i <- parent[i]
    }
    i
  }
  if (nrow(E)) {
    for (j in seq_len(nrow(E))) {
      a <- root(E[j, 1])
      b <- root(E[j, 2])
      if (a != b) parent[b] <- a
    }
  }
  z <- integer(n)
  for (i in seq_len(n)) z[i] <- root(i)
  match(z, unique(z))
}

# Summarize incident edge lengths at each point.
.spde_local_scale <- function(n, E, d) {
  z <- vector("list", n)
  for (j in seq_len(nrow(E))) {
    z[[E[j, 1]]] <- c(z[[E[j, 1]]], d[j])
    z[[E[j, 2]]] <- c(z[[E[j, 2]]], d[j])
  }
  out <- vapply(z, function(x) if (length(x)) stats::median(x) else NA_real_,
                numeric(1))
  if (anyNA(out)) out[is.na(out)] <- stats::median(d)
  out
}

# Compute triangle circumradii from point and triangle matrices.
.spde_circumradius <- function(P, T) {
  p1 <- P[T[, 1], , drop = FALSE]
  p2 <- P[T[, 2], , drop = FALSE]
  p3 <- P[T[, 3], , drop = FALSE]
  a <- sqrt(rowSums((p2 - p3)^2))
  b <- sqrt(rowSums((p1 - p3)^2))
  c <- sqrt(rowSums((p1 - p2)^2))
  area <- abs((p2[, 1] - p1[, 1]) * (p3[, 2] - p1[, 2]) -
              (p2[, 2] - p1[, 2]) * (p3[, 1] - p1[, 1])) / 2
  a * b * c / (4 * area)
}

# Normalize locations and construct their local Delaunay geometry.
.spde_prepare_points <- function(loc) {
  center <- colMeans(loc)
  scale <- max(diff(range(loc[, 1])), diff(range(loc[, 2])))
  if (!is.finite(scale) || scale <= 0) stop("loc must span a two-dimensional region.")
  z <- sweep(loc, 2, center, "-") / scale
  if (qr(sweep(z, 2, colMeans(z), "-"))$rank < 2L) {
    stop("loc must not be collinear.")
  }
  T <- geometry::delaunayn(z)
  E <- .spde_triangle_edges(T)
  d <- sqrt(rowSums((z[E[, 1], , drop = FALSE] -
                     z[E[, 2], , drop = FALSE])^2))
  h <- .spde_local_scale(nrow(z), E, d)
  score <- d / sqrt(h[E[, 1]] * h[E[, 2]])
  list(loc = loc, z = z, center = center, scale = scale,
       triangles = T, edges = E, local.scale = h, edge.score = score)
}

# Retain the largest polygon and at most the requested largest holes.
.spde_limit_holes <- function(domain, max_holes) {
  domain <- sf::st_cast(domain, "POLYGON")
  ar <- as.numeric(sf::st_area(domain))
  domain <- domain[which.max(ar)]
  ring <- domain[[1]]
  if (length(ring) <= max_holes + 1L) return(domain)
  if (max_holes == 0L) return(sf::st_sfc(sf::st_polygon(ring[1])))
  ha <- vapply(ring[-1], function(x) {
    as.numeric(sf::st_area(sf::st_sfc(sf::st_polygon(list(x)))))
  }, numeric(1))
  use <- order(ha, decreasing = TRUE)[seq_len(max_holes)]
  sf::st_sfc(sf::st_polygon(c(ring[1], ring[-1][use])))
}

# Construct one normalized candidate boundary and its diagnostics.
.spde_learn_candidate <- function(P, max_holes, group_size, seed,
                                  link_scale, alpha_scale, boundary_pad) {
  set.seed(seed)
  n <- nrow(P$z)
  K <- min(n - 1L, max(3L, ceiling(n / group_size)))
  t0 <- proc.time()[["elapsed"]]
  km <- stats::kmeans(P$z, centers = K, iter.max = 40L,
                      nstart = 1L, algorithm = "Lloyd")
  time.kmeans <- proc.time()[["elapsed"]] - t0

  E <- P$edges[P$edge.score <= link_scale, , drop = FALSE]
  Ein <- E[km$cluster[E[, 1]] == km$cluster[E[, 2]], , drop = FALSE]
  comp <- .spde_graph_components(n, Ein)
  grp <- as.integer(interaction(km$cluster, comp, drop = TRUE))

  ng <- tabulate(grp)
  small <- which(ng < max(3L, floor(group_size / 3)))
  for (g in small) {
    ii <- which(grp == g)
    hit <- E[E[, 1] %in% ii | E[, 2] %in% ii, , drop = FALSE]
    if (!nrow(hit)) next
    jj <- unique(c(hit[, 1], hit[, 2]))
    cand <- grp[jj]
    cand <- cand[cand != g]
    if (!length(cand)) next
    to <- as.integer(names(which.max(table(cand))))
    grp[ii] <- to
  }
  grp <- match(grp, unique(grp))
  ng <- tabulate(grp)
  G <- length(ng)
  C <- rowsum(P$z, grp) / ng

  support <- matrix(FALSE, G, G)
  ga <- grp[E[, 1]]
  gb <- grp[E[, 2]]
  ii <- ga != gb
  support[cbind(ga[ii], gb[ii])] <- TRUE
  support[cbind(gb[ii], ga[ii])] <- TRUE

  Tc <- geometry::delaunayn(C)
  Ec <- .spde_triangle_edges(Tc)
  dc <- sqrt(rowSums((C[Ec[, 1], , drop = FALSE] -
                      C[Ec[, 2], , drop = FALSE])^2))
  hc <- .spde_local_scale(G, Ec, dc)
  R <- .spde_circumradius(C, Tc)
  hs <- apply(matrix(hc[Tc], ncol = 3), 1, stats::median)
  keep <- support[cbind(Tc[, 1], Tc[, 2])] &
    support[cbind(Tc[, 1], Tc[, 3])] &
    support[cbind(Tc[, 2], Tc[, 3])] & R / hs <= alpha_scale
  Tc <- Tc[keep, , drop = FALSE]
  if (!nrow(Tc)) stop("No supported triangles remained; increase link_scale or alpha_scale.")

  tri <- sf::st_sfc(lapply(seq_len(nrow(Tc)), function(j) {
    x <- C[Tc[j, ], , drop = FALSE]
    sf::st_polygon(list(rbind(x, x[1, ])))
  }))
  hull <- vector("list", G)
  for (g in seq_len(G)) {
    x <- P$z[grp == g, , drop = FALSE]
    if (nrow(x) >= 3L) {
      h <- grDevices::chull(x)
      x <- x[c(h, h[1]), , drop = FALSE]
    } else {
      r <- stats::median(P$local.scale[grp == g]) / 3
      x <- rbind(C[g, ] + c(-r, -r), C[g, ] + c(r, -r),
                 C[g, ] + c(r, r), C[g, ] + c(-r, r),
                 C[g, ] + c(-r, -r))
    }
    hull[[g]] <- sf::st_polygon(list(x))
  }

  t1 <- proc.time()[["elapsed"]]
  domain <- sf::st_make_valid(sf::st_union(c(tri, sf::st_sfc(hull))))
  pad <- boundary_pad * stats::median(P$local.scale)
  domain <- sf::st_buffer(domain, dist = pad, nQuadSegs = 1)
  domain <- sf::st_simplify(domain, dTolerance = pad / 2,
                            preserveTopology = TRUE)
  if (!all(sf::st_geometry_type(domain) %in% c("POLYGON", "MULTIPOLYGON"))) {
    domain <- sf::st_collection_extract(domain, "POLYGON")
  }
  domain <- .spde_limit_holes(domain, max_holes)
  pts <- sf::st_as_sf(data.frame(x = P$z[, 1], y = P$z[, 2]),
                      coords = c("x", "y"))
  coverage <- mean(lengths(sf::st_covered_by(pts, domain)) > 0L)
  time.geometry <- proc.time()[["elapsed"]] - t1

  list(domain = domain, groups = grp, centers = C, triangles = Tc,
       seed = seed, n.groups = G, holes = length(domain[[1]]) - 1L,
       coverage = coverage, kmeans.seconds = time.kmeans,
       geometry.seconds = time.geometry,
       total.seconds = proc.time()[["elapsed"]] - t0)
}

# Compute polygon intersection over union for candidate consensus.
.spde_domain_iou <- function(a, b) {
  x <- sum(as.numeric(sf::st_area(sf::st_intersection(a, b))))
  y <- sum(as.numeric(sf::st_area(sf::st_union(a, b))))
  x / y
}

#' Learn a spatial boundary from observation coordinates
#'
#' Uses local K-means landmarks, graph-connected group correction, supported
#' Delaunay triangles, and cross-seed consensus to construct a polygon with at
#' most `max_holes` holes, or exactly `n_holes` holes when that argument is
#' supplied. The default group size grows as `ceiling(log(n))`.
#'
#' @param loc Observation coordinates or an `sf` point object.
#' @param max_holes Maximum number of retained holes.
#' @param n_holes Optional exact number of retained holes specified from
#'   scientific knowledge. When supplied, candidate boundaries with a different
#'   number of holes are rejected.
#' @param group_size Target observations per initial K-means group.
#' @param seeds Integer random seeds used to generate candidate boundaries.
#' @param link_scale Maximum locally standardized point-edge length.
#' @param alpha_scale Maximum locally standardized triangle circumradius.
#' @param boundary_pad Boundary expansion in units of local point spacing.
#' @param crs Optional coordinate reference system understood by [sf::st_crs()].
#' @return A `spde_boundary` object containing the selected domain, candidate
#'   diagnostics, and candidate geometries.
#' @export
spde_boundary_learn <- function(loc, max_holes = 0L, n_holes = NULL,
                                group_size = NULL, seeds = 1:10,
                                link_scale = 2.25, alpha_scale = 1.15,
                                boundary_pad = 0.35,
                                crs = sf::NA_crs_) {
  loc <- .spde_xy(loc)
  n <- nrow(loc)
  max_holes <- as.integer(max_holes)
  if (length(max_holes) != 1L || is.na(max_holes) || max_holes < 0L) {
    stop("max_holes must be one non-negative integer.")
  }
  if (!is.null(n_holes)) {
    n_holes <- as.integer(n_holes)
    if (length(n_holes) != 1L || is.na(n_holes) || n_holes < 0L) {
      stop("n_holes must be NULL or one non-negative integer.")
    }
    max_holes <- n_holes
  }
  if (is.null(group_size)) group_size <- ceiling(log(n))
  group_size <- as.integer(group_size)
  if (length(group_size) != 1L || is.na(group_size) || group_size < 3L ||
      group_size >= n) {
    stop("group_size must be an integer from 3 to n - 1.")
  }
  seeds <- as.integer(seeds)
  if (!length(seeds) || anyNA(seeds)) stop("seeds must contain finite integers.")
  pars <- c(link_scale = link_scale, alpha_scale = alpha_scale,
            boundary_pad = boundary_pad)
  if (any(!is.finite(pars)) || link_scale <= 0 || alpha_scale <= 0 ||
      boundary_pad < 0) {
    stop("link_scale and alpha_scale must be positive; boundary_pad must be non-negative.")
  }

  P <- .spde_prepare_points(loc)
  M <- vector("list", length(seeds))
  for (j in seq_along(seeds)) {
    M[[j]] <- .spde_learn_candidate(
      P, max_holes = max_holes, group_size = group_size, seed = seeds[j],
      link_scale = link_scale, alpha_scale = alpha_scale,
      boundary_pad = boundary_pad
    )
  }

  B <- length(M)
  J <- diag(1, B)
  if (B > 1L) {
    for (i in seq_len(B - 1L)) {
      for (j in (i + 1L):B) {
        J[i, j] <- J[j, i] <- .spde_domain_iou(M[[i]]$domain, M[[j]]$domain)
      }
    }
  }
  hole.ok <- rep(TRUE, B)
  if (!is.null(n_holes)) {
    hole.ok <- vapply(M, function(x) x$holes == n_holes, logical(1))
    if (!any(hole.ok)) {
      stop("No candidate boundary had exactly ", n_holes,
           " hole(s); adjust boundary-learning parameters or seeds.")
    }
  }
  eligible <- which(hole.ok &
                    vapply(M, function(x) x$coverage == 1, logical(1)))
  if (!length(eligible)) {
    eligible <- which(hole.ok &
                      vapply(M, function(x) x$coverage >= 0.995, logical(1)))
  }
  if (!length(eligible)) {
    suffix <- if (is.null(n_holes)) "." else
      paste0(" while preserving exactly ", n_holes, " hole(s).")
    stop("No candidate boundary covered at least 99.5% of locations", suffix)
  }
  score <- rep(-Inf, B)
  score[eligible] <- rowMeans(J[eligible, eligible, drop = FALSE])
  selected <- which.max(score)
  best <- M[[selected]]

  rings <- lapply(best$domain[[1]], function(x) {
    sweep(x * P$scale, 2, P$center, "+")
  })
  domain <- .spde_make_polygon(rings, crs = sf::st_crs(crs))
  tab <- data.frame(
    seed = seeds,
    groups = vapply(M, `[[`, integer(1), "n.groups"),
    holes = vapply(M, `[[`, integer(1), "holes"),
    coverage = vapply(M, `[[`, numeric(1), "coverage"),
    consensus = score,
    kmeans_seconds = vapply(M, `[[`, numeric(1), "kmeans.seconds"),
    geometry_seconds = vapply(M, `[[`, numeric(1), "geometry.seconds"),
    total_seconds = vapply(M, `[[`, numeric(1), "total.seconds"),
    selected = seq_len(B) == selected
  )

  .new_spde_boundary(
    domain, method = "learned",
    diagnostics = list(
      holes = best$holes, requested_holes = n_holes,
      coverage = best$coverage,
      group_size = group_size, selected_seed = best$seed,
      consensus = score[selected], candidates = tab,
      transform = list(center = P$center, scale = P$scale)
    ),
    candidates = lapply(M, `[[`, "domain"), loc = loc
  )
}
