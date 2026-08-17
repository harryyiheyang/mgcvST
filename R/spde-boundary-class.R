# Normalize matrix or sf locations to finite x-y coordinates.
.spde_xy <- function(loc) {
  if (inherits(loc, "sf")) {
    loc <- sf::st_coordinates(loc)[, 1:2, drop = FALSE]
  } else {
    loc <- as.matrix(loc)
  }
  if (!is.numeric(loc) || ncol(loc) != 2L || nrow(loc) < 3L ||
      any(!is.finite(loc))) {
    stop("loc must contain at least three finite two-dimensional coordinates.")
  }
  colnames(loc) <- c("x", "y")
  loc
}

# Close and validate one polygon ring after removing adjacent duplicates.
.spde_close_ring <- function(x) {
  x <- as.matrix(x)
  if (!is.numeric(x) || ncol(x) != 2L || nrow(x) < 3L ||
      any(!is.finite(x))) {
    stop("Each boundary ring must contain at least three finite points.")
  }
  keep <- c(TRUE, rowSums(x[-1, , drop = FALSE] !=
                          x[-nrow(x), , drop = FALSE]) > 0)
  x <- x[keep, , drop = FALSE]
  if (nrow(x) < 3L) stop("A boundary ring collapsed below three unique points.")
  if (!all(x[1, ] == x[nrow(x), ])) x <- rbind(x, x[1, ])
  x
}

# Build and validate an sfc polygon from an ordered list of rings.
.spde_make_polygon <- function(rings, crs = sf::NA_crs_) {
  rings <- lapply(rings, .spde_close_ring)
  domain <- sf::st_sfc(sf::st_polygon(rings), crs = crs)
  valid <- sf::st_is_valid(domain, reason = TRUE)
  if (!all(valid == "Valid Geometry")) {
    stop("The boundary does not define a valid polygon: ",
         paste(unique(valid), collapse = "; "))
  }
  domain
}

# Store a validated domain and its construction diagnostics.
.new_spde_boundary <- function(domain, method, diagnostics = list(),
                               candidates = NULL, loc = NULL) {
  if (!inherits(domain, c("sf", "sfc"))) {
    stop("domain must be an sf or sfc polygon object.")
  }
  out <- list(
    domain = sf::st_geometry(domain),
    method = method,
    diagnostics = diagnostics,
    candidates = candidates,
    loc = loc
  )
  class(out) <- "spde_boundary"
  out
}

# Extract polygon geometry from any supported boundary input.
.spde_domain <- function(x) {
  if (inherits(x, "spde_boundary")) return(x$domain)
  if (inherits(x, c("sf", "sfc"))) return(sf::st_geometry(x))
  stop("boundary must be a spde_boundary, sf, or sfc polygon object.")
}

#' Inspect an SPDE boundary
#'
#' @param x A `spde_boundary` object.
#' @param ... Unused.
#' @export
print.spde_boundary <- function(x, ...) {
  d <- x$diagnostics
  holes <- if (!is.null(d$holes)) d$holes else
    sum(vapply(x$domain, function(g) length(g) - 1L, integer(1)))
  cat("SPDE boundary\n")
  cat("  method: ", x$method, "\n", sep = "")
  cat("  holes: ", holes, "\n", sep = "")
  if (!is.null(d$coverage)) {
    cat("  location coverage: ", format(d$coverage), "\n", sep = "")
  }
  invisible(x)
}

#' Plot an SPDE boundary
#'
#' @param x A `spde_boundary` object.
#' @param ... Arguments passed to [plot()].
#' @export
plot.spde_boundary <- function(x, ...) {
  plot(x$domain, ...)
  if (!is.null(x$loc)) graphics::points(x$loc, pch = 16, cex = 0.35)
  invisible(x)
}
