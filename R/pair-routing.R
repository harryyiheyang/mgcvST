# Choice of the pair route. Both routes need only the dense score state
# (a_g, H_g) of every gene, so the choice depends on the score dimension q, the
# number of pairs P, the threads and the memory, not on the estimator.
#
# The exact route evaluates, for every pair, the four trace moments of the
# product of the two q x q curvature matrices (cost cubic in q) and a k x k
# singular value problem on the shared basis. The PCAlearning route replaces
# the trace moments by a contraction of low-rank coefficients and the k x k
# problem by the Cholesky factors of the compressed matrices, so that its pair
# cost does not depend on q.

# Seconds per pair and thread of the exact route at score dimension q.
.mgcvst_exact_pair_seconds <- function(q) 5e-3 * (q / 298)^3

# Seconds per pair and thread of the PCAlearning route with k leading singular
# values and rank r: the k x k singular value problem dominates (0.47 ms at
# k = 50), the contraction of the degree-2 monomials costs r^2 operations.
.mgcvst_pca_pair_seconds <- function(k, r) {
  3.8e-9 * k^3 + 2.5e-8 * r^2 + 1e-6
}

# Budget of the exact route: the longest pair phase that `auto` accepts, and
# the largest score-state store written to disk.
.mgcvst_route_limits <- function() {
  list(exact_seconds = 2 * 3600, store_bytes = 64 * 1024^3, resident_fraction = 0.3)
}

.mgcvst_format_duration <- function(seconds) {
  if (!is.finite(seconds)) return("an unknown time")
  if (seconds < 90) return(sprintf("%.0f s", seconds))
  if (seconds < 5400) return(sprintf("%.0f min", seconds / 60))
  if (seconds < 3 * 86400) return(sprintf("%.1f h", seconds / 3600))
  sprintf("%.1f days", seconds / 86400)
}

# Resolve `moments` to "exact" or "pcalearning". `auto` takes the exact route
# when its pair phase fits the time budget and its resident pair bases and
# score-state store fit the memory and disk budgets, or when there are no more
# genes than the PCAlearning rank; otherwise PCAlearning.
.mgcvst_route_resolve <- function(moments, q, n_used, n_pairs, threads, rank,
                                  k_exact, k_pca, available_memory) {
  k_exact <- as.integer(min(k_exact, q))
  k_pca <- as.integer(min(k_pca, q))
  seconds_exact <- n_pairs * .mgcvst_exact_pair_seconds(q) / threads
  seconds_pca <- n_pairs * .mgcvst_pca_pair_seconds(k_pca, rank) / threads
  resident <- 8 * q * k_exact * n_used
  store <- 8 * (q^2 + q) * n_used
  limits <- .mgcvst_route_limits()
  memory_ok <- !is.finite(available_memory) ||
    resident <= limits$resident_fraction * available_memory
  disk_ok <- store <= limits$store_bytes
  time_ok <- seconds_exact <= limits$exact_seconds
  reason <- NULL
  chosen <- moments
  if (identical(moments, "auto")) {
    if (n_used <= rank) {
      chosen <- "exact"
      reason <- sprintf("%d genes do not exceed the PCAlearning rank %d",
                        n_used, rank)
    } else if (time_ok && memory_ok && disk_ok) {
      chosen <- "exact"
      reason <- sprintf("the exact pair phase takes an estimated %s, within %s",
                        .mgcvst_format_duration(seconds_exact),
                        .mgcvst_format_duration(limits$exact_seconds))
    } else {
      chosen <- "pcalearning"
      reason <- if (!time_ok) {
        sprintf("the exact pair phase would take an estimated %s, above %s",
                .mgcvst_format_duration(seconds_exact),
                .mgcvst_format_duration(limits$exact_seconds))
      } else if (!memory_ok) {
        sprintf("the exact pair bases need %.1f GB of the %.1f GB available",
                resident / 1024^3, available_memory / 1024^3)
      } else {
        sprintf("the exact score-state store needs %.1f GB on disk",
                store / 1024^3)
      }
    }
  }
  list(moments = chosen, reason = reason, seconds = if (identical(chosen, "exact"))
         seconds_exact else seconds_pca,
       seconds_exact = seconds_exact, seconds_pcalearning = seconds_pca,
       k = if (identical(chosen, "exact")) k_exact else k_pca, q = q,
       n_used = n_used, n_pairs = n_pairs)
}

# The route of a checkpoint directory, written when the route is chosen so that
# a resumed run follows it whatever the threads or the memory of the new
# session. NULL when the directory holds none.
.mgcvst_route_stored <- function(checkpoint_dir) {
  file <- file.path(checkpoint_dir, "route.rds")
  if (!file.exists(file)) return(NULL)
  z <- tryCatch(readRDS(file), error = function(e) NULL)
  if (is.list(z) && z$moments %in% c("exact", "pcalearning")) z else NULL
}

.mgcvst_route_save <- function(checkpoint_dir, route) {
  if (is.null(route)) return(invisible(NULL))
  file <- file.path(checkpoint_dir, "route.rds")
  if (file.exists(file)) return(invisible(NULL))
  tmp <- tempfile("route-", tmpdir = checkpoint_dir, fileext = ".tmp")
  on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
  saveRDS(route[c("moments", "k", "q")], tmp, compress = FALSE)
  if (!file.rename(tmp, file) && !file.exists(file)) {
    stop("Could not commit the route record.")
  }
  invisible(NULL)
}
