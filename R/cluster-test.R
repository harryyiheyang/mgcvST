#' Test clusters of features against single features
#'
#' Tests every cluster of features against every feature of `features` given
#' the cluster: the score vectors of the members of a cluster are summed,
#' `a_S = sum_g a_g`, and the score of the cluster with feature `c` is
#' `U = a_S' a_c`. When the spatial field of `c` is unrelated to the fields of
#' the members, `a_c` has mean zero and covariance `M_c` (the state matrix of
#' `c`) whatever the dependence among the members, so given `a_S` the score is
#' a linear form of `a_c` with mean zero and variance `V = a_S' M_c a_S`. The
#' saddlepoint approximation of a linear form of a Gaussian vector is its
#' normal tail, and the two-sided p-value of `z = U / sqrt(V)` is
#' `2 * pnorm(-abs(z))`. The test holds for clusters chosen from the data as
#' long as the choice did not use the tested features, for example modules of
#' [spdeST.wgcna()] built from the RNA features alone and tested against
#' proteins or chromatin features.
#'
#' The estimates are those of the pair test pooled over the members. The
#' information of the cluster with `c` is `T = tr(M_S M_c)`, with `M_S` the
#' summed state matrices of the members, and `rho = U / T` is the latent
#' spatial correlation of `c` with the members, averaged with weights
#' `tr(M_g M_c)`. Within a cluster the same average over the pairs of
#' different members is `rho_cluster = (a_S' a_S - sum_g a_g' a_g) /
#' (tr(M_S^2) - sum_g tr(M_g^2))`. When the members share one spatial factor,
#' `rho_cluster` estimates the squared correlation of a member with that
#' factor, and `k = rho / sqrt(rho_cluster)` the correlation of `c` with the
#' factor of the cluster; `k` is `NA` when `rho_cluster` is not positive.
#'
#' A feature is left out of a cluster that holds it, and `exclude` leaves out
#' further members, for example the RNA feature of the gene of a chromatin
#' window, so that a window is tested against the other genes of its cluster.
#'
#' A group of features is combined against every cluster: its score and
#' information are the sums over its features, `rho = sum U / sum T`, the
#' within-group estimate `rho_group` is formed as `rho_cluster`, and
#' `k = rho / sqrt(rho_group * rho_cluster)` estimates the correlation of the
#' factor of the group with that of the cluster. The p-value of the group is
#' the Cauchy combination of the two-sided p-values of its features, which
#' holds whatever their dependence.
#'
#' The adjustment is applied to the rows with a p-value of `results` as one
#' family and to the rows of `groups` as another, in log space, as in
#' [spdeST.test()].
#'
#' @param fit Result of [spdeST.estimate()].
#' @param clusters Named list of feature IDs or one-based indices, one entry
#'   per cluster, each with at least two available features.
#' @param features Feature IDs or one-based indices of available features,
#'   each tested against every cluster.
#' @param groups `NULL`, or a named list of subsets of `features`, each with at
#'   least two features.
#' @param exclude `NULL`, or a two-column matrix or data frame of feature IDs
#'   or indices: when the feature of the first column is tested against a
#'   cluster that holds the feature of the second column, that member is left
#'   out of the cluster.
#' @param q.value Discovery threshold on the adjusted q-value, in `(0, 1]`.
#' @param adjust Multiple-testing adjustment: `"BY"` (the default), `"BH"`,
#'   `"Sidak"` or `"none"`.
#' @param threads Positive number of OpenMP threads; `NULL` uses one.
#' @param verbose Logical; report progress.
#' @return An object of class `spdeST_clustertest` with
#'   * `results`: one row per cluster and feature, with `cluster`, `feature`,
#'     `members` (the members used), `score` (`U`), `variance` (`V`),
#'     `information` (`T`), `z`, `log_p_two_sided`, `log_q`, `rho` and `k`;
#'     a row without a member left has no p-value;
#'   * `groups`: one row per group and cluster (`NULL` without groups), with
#'     `group`, `cluster`, `features`, `score`, `information`, `rho`,
#'     `rho_group`, `rho_cluster`, `k`, `log_p_two_sided` and `log_q`;
#'   * `clusters`: one row per cluster with `cluster`, `members` and
#'     `rho_cluster`;
#'   * `feature_id`, `threshold`, `basis` (the score dimension `q`), `timing`
#'     and `call`.
#' @export
spdeST.clustertest <- function(fit, clusters, features, groups = NULL,
                               exclude = NULL, q.value = 0.05,
                               adjust = c("BY", "BH", "Sidak", "none"),
                               threads = NULL, verbose = FALSE) {
  started <- proc.time()[["elapsed"]]
  if (!inherits(fit, "spdeST_fit")) {
    stop("fit must be returned by spdeST.estimate().")
  }
  .mgcvst_check_fit_format(fit)
  adjust <- match.arg(adjust)
  q.value <- as.numeric(q.value)
  if (length(q.value) != 1L || !is.finite(q.value) || q.value <= 0 ||
      q.value > 1) {
    stop("q.value must be one finite value in (0, 1].")
  }
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("verbose must be TRUE or FALSE.")
  }
  if (is.null(threads)) threads <- 1L
  if (!is.numeric(threads) || length(threads) != 1L || !is.finite(threads) ||
      threads < 1 || threads != floor(threads)) {
    stop("threads must be one positive integer.")
  }
  threads <- as.integer(threads)
  id <- fit$feature_id
  available <- .mgcvst_feature_available(fit)
  cl <- .mgcvst_cluster_sets(clusters, id, available, "clusters")
  ft <- .mgcvst_cluster_index(features, id, "features")
  if (anyDuplicated(ft)) stop("features contains duplicated features.")
  if (!all(available[ft])) {
    stop("features without a fitted spatial model: ",
         paste(id[ft[!available[ft]]], collapse = ", "), ".")
  }
  gr <- if (is.null(groups)) NULL else .mgcvst_cluster_sets(groups, id, available, "groups")
  if (!is.null(gr) && !all(unlist(gr) %in% ft)) {
    stop("Every feature of groups must be in features.")
  }
  ex <- if (is.null(exclude)) matrix(integer(), 0L, 2L) else {
    x <- as.matrix(exclude)
    if (length(dim(x)) != 2L || ncol(x) != 2L) {
      stop("exclude must be a two-column matrix or data frame.")
    }
    matrix(.mgcvst_cluster_index(as.vector(x), id, "exclude"), ncol = 2L)
  }
  ex <- unique(rbind(ex[ex[, 1L] %in% ft, , drop = FALSE], cbind(ft, ft)))
  old_threads <- .mgcvst_thread_limit()
  on.exit(.mgcvst_thread_restore(old_threads), add = TRUE)

  fit <- .spdest_sparse_prepare(fit)
  basis <- .spdest_check_basis(fit, .spdest_sparse_observation_basis(fit))
  q <- as.integer(basis$rank)
  members <- sort(unique(unlist(cl)))
  keep <- intersect(members, ex[, 2L])
  builder <- .mgcvst_state_builder(fit, unique(c(members, ft)), threads, basis)
  batches <- function(ids) {
    probe <- .mgcvst_memory_probe()
    room <- if (is.finite(probe$available)) 0.3 * probe$available else 512 * 1024^2
    if (room < builder$work) {
      stop("Insufficient available memory for one score-state preparation. ",
           "Increase the job memory allocation.")
    }
    size <- as.integer(max(1L, min(32L, floor(room / builder$work))))
    split(ids, ceiling(seq_along(ids) / size))
  }
  build <- function(ids) {
    s <- builder$build(ids)
    bad <- vapply(s, function(z) !is.null(z$error), logical(1L))
    if (any(bad)) {
      stop("No score state for ", paste(id[ids[bad]], collapse = ", "), ": ",
           s[[which(bad)[1L]]]$error)
    }
    s
  }

  # The summed score vectors and state matrices of the clusters, and the
  # states of the members that some feature leaves out.
  K <- length(cl)
  need <- 8 * q^2 * (K + length(keep) + if (is.null(gr)) 0 else length(gr))
  probe <- .mgcvst_memory_probe()
  if (is.finite(probe$available) && need > 0.5 * probe$available) {
    stop(sprintf(paste0("The summed states of the clusters and groups and the states of the ",
                        "members left out need about %.1f GB, above half of the %.1f GB ",
                        "available."), need / 1024^3, probe$available / 1024^3))
  }
  aS <- matrix(0, q, K)
  MS <- replicate(K, matrix(0, q, q), simplify = FALSE)
  self <- matrix(0, K, 2L)
  kept <- list()
  done <- 0L
  for (ids in batches(members)) {
    s <- build(ids)
    for (j in seq_along(ids)) {
      a <- s[[j]]$a
      M <- s[[j]]$M
      for (k in which(vapply(cl, function(z) ids[j] %in% z, logical(1L)))) {
        aS[, k] <- aS[, k] + a
        MS[[k]] <- MS[[k]] + M
        self[k, ] <- self[k, ] + c(sum(a * a), sum(M * t(M)))
      }
      if (ids[j] %in% keep) kept[[as.character(ids[j])]] <- s[[j]]
    }
    done <- done + length(ids)
    if (verbose) message("Summed the score states of ", done, " of ", length(members),
                         " cluster members.")
  }
  rho_cluster <- vapply(seq_len(K), function(k) {
    (sum(aS[, k]^2) - self[k, 1L]) / (sum(MS[[k]] * t(MS[[k]])) - self[k, 2L])
  }, numeric(1L))
  ae <- if (length(kept)) do.call(cbind, lapply(kept, `[[`, "a")) else matrix(0, q, 0L)
  Me <- lapply(kept, `[[`, "M")
  at <- match(as.character(keep), names(kept))

  # Every feature against every cluster, with the members it leaves out.
  G <- if (is.null(gr)) 0L else length(gr)
  aG <- matrix(0, q, G)
  MG <- replicate(G, matrix(0, q, q), simplify = FALSE)
  gself <- matrix(0, G, 2L)
  out <- vector("list", length(ft))
  for (ids in batches(ft)) {
    s <- build(ids)
    a <- do.call(cbind, lapply(s, `[[`, "a"))
    M <- lapply(s, `[[`, "M")
    x <- mgcvst_cluster_moments_cpp(aS, MS, a, M, threads)
    used <- matrix(lengths(cl), K, length(ids))
    rows <- list()
    for (j in seq_along(ids)) {
      drop <- ex[ex[, 1L] == ids[j], 2L]
      for (k in seq_len(K)) {
        g <- intersect(drop, cl[[k]])
        if (length(g)) {
          used[k, j] <- used[k, j] - length(g)
          rows[[length(rows) + 1L]] <- list(j = j, k = k, g = at[match(g, keep)])
        }
      }
      for (h in which(vapply(gr, function(z) ids[j] %in% z, logical(1L)))) {
        aG[, h] <- aG[, h] + a[, j]
        MG[[h]] <- MG[[h]] + M[[j]]
        gself[h, ] <- gself[h, ] + c(sum(a[, j]^2), sum(M[[j]] * t(M[[j]])))
      }
    }
    if (length(rows)) {
      g <- lapply(rows, `[[`, "g")
      d <- mgcvst_cluster_exclude_cpp(
        aS, a, M, ae, Me,
        vapply(rows, `[[`, integer(1L), "j"), vapply(rows, `[[`, integer(1L), "k"),
        c(0L, cumsum(lengths(g))), as.integer(unlist(g)), threads)
      cell <- cbind(vapply(rows, `[[`, integer(1L), "k"), vapply(rows, `[[`, integer(1L), "j"))
      x$U[cell] <- x$U[cell] - d$dU
      x$V[cell] <- x$V[cell] + d$dV
      x$T[cell] <- x$T[cell] - d$dT
    }
    out[match(ids, ft)] <- lapply(seq_along(ids), function(j) data.frame(
      cluster = names(cl), feature = id[ids[j]], members = used[, j],
      score = x$U[, j], variance = x$V[, j], information = x$T[, j],
      stringsAsFactors = FALSE))
    if (verbose) message("Tested ", sum(lengths(out) > 0L), " of ", length(ft),
                         " features against ", K, " clusters.")
  }
  res <- do.call(rbind, out)
  ok <- res$members > 0L & is.finite(res$variance) & res$variance > 0
  res$z <- ifelse(ok, res$score / sqrt(res$variance), NA_real_)
  res$log_p_two_sided <- log(2) + stats::pnorm(-abs(res$z), log.p = TRUE)
  res$log_q <- NA_real_
  if (any(ok)) res$log_q[ok] <- .mgcvst_log_adjust(res$log_p_two_sided[ok], adjust)$log_q
  res$rho <- ifelse(res$members > 0L, res$score / res$information, NA_real_)
  rc <- rho_cluster[match(res$cluster, names(cl))]
  res$k <- ifelse(rc > 0, res$rho / sqrt(pmax(rc, 0)), NA_real_)
  rownames(res) <- NULL

  grp <- NULL
  if (G) {
    rho_group <- vapply(seq_len(G), function(h) {
      (sum(aG[, h]^2) - gself[h, 1L]) / (sum(MG[[h]] * t(MG[[h]])) - gself[h, 2L])
    }, numeric(1L))
    grp <- do.call(rbind, lapply(seq_len(G), function(h) do.call(rbind, lapply(seq_len(K), function(k) {
      z <- res[res$cluster == names(cl)[k] & res$feature %in% id[gr[[h]]] & ok, ]
      data.frame(group = names(gr)[h], cluster = names(cl)[k], features = nrow(z),
                 score = sum(z$score), information = sum(z$information),
                 log_p_two_sided = .mgcvst_cauchy_log_p(z$log_p_two_sided),
                 stringsAsFactors = FALSE)
    }))))
    grp$rho <- grp$score / grp$information
    grp$rho_group <- rho_group[match(grp$group, names(gr))]
    grp$rho_cluster <- rho_cluster[match(grp$cluster, names(cl))]
    both <- grp$rho_group > 0 & grp$rho_cluster > 0
    grp$k <- ifelse(both, grp$rho / sqrt(pmax(grp$rho_group * grp$rho_cluster, 0)), NA_real_)
    grp$log_q <- NA_real_
    good <- is.finite(grp$log_p_two_sided)
    if (any(good)) grp$log_q[good] <- .mgcvst_log_adjust(grp$log_p_two_sided[good], adjust)$log_q
    grp <- grp[, c("group", "cluster", "features", "score", "information", "rho",
                   "rho_group", "rho_cluster", "k", "log_p_two_sided", "log_q")]
  }
  structure(list(
    results = res, groups = grp,
    clusters = data.frame(cluster = names(cl), members = lengths(cl),
                          rho_cluster = rho_cluster, stringsAsFactors = FALSE),
    feature_id = id, threshold = list(q_value = q.value, adjust = adjust),
    basis = list(kind = basis$kind, q = q),
    timing = list(elapsed = proc.time()[["elapsed"]] - started, threads = threads),
    call = match.call()), class = "spdeST_clustertest")
}

#' Print a cluster test
#'
#' @param x An `spdeST_clustertest` object.
#' @param ... Unused.
#' @return `x`, invisibly.
#' @export
print.spdeST_clustertest <- function(x, ...) {
  cut <- log(x$threshold$q_value)
  r <- x$results
  cat("spdeST cluster tests\n")
  cat("  clusters:", nrow(x$clusters), " features:", length(unique(r$feature)), "\n")
  cat("  adjustment:", x$threshold$adjust, "\n")
  cat("  discoveries at q <=", format(x$threshold$q_value), ":",
      sum(r$log_q <= cut, na.rm = TRUE), "of", sum(!is.na(r$log_q)), "\n")
  if (!is.null(x$groups)) {
    cat("  group discoveries:", sum(x$groups$log_q <= cut, na.rm = TRUE), "of",
        sum(!is.na(x$groups$log_q)), "\n")
  }
  invisible(x)
}

# Feature IDs or one-based indices to indices of `feature_id`.
.mgcvst_cluster_index <- function(x, feature_id, what) {
  if (is.character(x)) {
    i <- match(x, feature_id)
    if (anyNA(i)) {
      stop("Unknown feature ID in ", what, ": ",
           paste(unique(x[is.na(i)]), collapse = ", "), ".")
    }
    return(i)
  }
  if (!is.numeric(x) || any(!is.finite(x)) || any(x != floor(x)) ||
      any(x < 1 | x > length(feature_id))) {
    stop(what, " must hold feature IDs or one-based feature indices.")
  }
  as.integer(x)
}

# A named list of feature sets to indices, each with at least two distinct
# available features.
.mgcvst_cluster_sets <- function(x, feature_id, available, what) {
  if (!is.list(x) || !length(x) || is.null(names(x)) || anyNA(names(x)) ||
      any(!nzchar(names(x))) || anyDuplicated(names(x))) {
    stop(what, " must be a list with unique non-empty names.")
  }
  lapply(stats::setNames(names(x), names(x)), function(nm) {
    i <- .mgcvst_cluster_index(x[[nm]], feature_id, what)
    if (length(i) < 2L || anyDuplicated(i)) {
      stop("Entry '", nm, "' of ", what, " must hold at least two distinct features.")
    }
    if (!all(available[i])) {
      stop("Entry '", nm, "' of ", what, " holds features without a fitted spatial model: ",
           paste(feature_id[i[!available[i]]], collapse = ", "), ".")
    }
    i
  })
}

# Natural log of the Cauchy combination of two-sided p-values given in log
# space, with equal weights: T = mean(tan((1/2 - p) pi)) and p = 1/2 -
# atan(T) / pi. A term with p below 1e-16 is 1 / (p pi) and these terms are
# summed in log space, so that p-values below the double range keep their
# order.
.mgcvst_cauchy_log_p <- function(log_p) {
  log_p <- log_p[!is.na(log_p)]
  n <- length(log_p)
  if (!n) return(NA_real_)
  if (any(log_p == -Inf)) return(-Inf)
  small <- log_p < log(1e-16)
  rest <- sum(tan((0.5 - exp(log_p[!small])) * pi))
  big <- -Inf
  if (any(small)) {
    lt <- -log(pi) - log_p[small]
    big <- max(lt) + log(sum(exp(lt - max(lt))))
  }
  if (big > 700) return(-log(pi) - big + log(n))
  stats::pcauchy((exp(big) + rest) / n, lower.tail = FALSE, log.p = TRUE)
}
