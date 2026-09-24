# Normalize named WGCNA overrides without changing unspecified defaults.
.mgcvst_wgcna_parameters <- function(wgcna.para) {
  out <- list(networkType = "signed", power = 6, TOMType = "signed",
              hclustMethod = "average", minClusterSize = 20L, deepSplit = 1L)
  if (is.null(wgcna.para)) return(out)
  if (!is.list(wgcna.para)) stop("wgcna.para must be NULL or a named list.")
  if (!length(wgcna.para)) return(out)
  nm <- names(wgcna.para)
  if (is.null(nm) || anyNA(nm) || any(!nzchar(nm)) || anyDuplicated(nm)) {
    stop("Every wgcna.para entry must have a unique non-empty name.")
  }
  if (any(!nm %in% names(out))) {
    stop("Unknown wgcna.para setting(s): ", paste(setdiff(nm, names(out)), collapse = ", "),
         ". Supported settings: ", paste(names(out), collapse = ", "), ".")
  }
  for (x in nm) out[x] <- wgcna.para[x]
  choices <- list(networkType = c("unsigned", "signed", "signed hybrid"),
                  TOMType = c("unsigned", "signed", "signed Nowick", "unsigned 2",
                              "signed 2", "signed Nowick 2"),
                  hclustMethod = c("ward.D", "ward.D2", "single", "complete",
                                   "average", "mcquitty", "median", "centroid"))
  for (x in names(choices)) {
    if (!is.character(out[[x]]) || length(out[[x]]) != 1L ||
        is.na(out[[x]]) || !out[[x]] %in% choices[[x]]) {
      stop("wgcna.para$", x, " must be one of: ", paste(choices[[x]], collapse = ", "), ".")
    }
  }
  for (x in c("power", "minClusterSize", "deepSplit")) {
    z <- out[[x]]
    if (!is.numeric(z) || length(z) != 1L || !is.finite(z)) {
      stop("wgcna.para$", x, " must be one finite number.")
    }
  }
  if (out$power <= 0) stop("wgcna.para$power must be positive.")
  if (out$minClusterSize < 2 || out$minClusterSize != floor(out$minClusterSize)) {
    stop("wgcna.para$minClusterSize must be an integer of at least 2.")
  }
  if (out$deepSplit < 0 || out$deepSplit > 4 || out$deepSplit != floor(out$deepSplit)) {
    stop("wgcna.para$deepSplit must be an integer from 0 to 4.")
  }
  out$minClusterSize <- as.integer(out$minClusterSize)
  out$deepSplit <- as.integer(out$deepSplit)
  out
}

# Resolve explicit gene blocks, preserving each requested order.
.mgcvst_wgcna_indices <- function(indices, feature_id) {
  blocks <- if (is.list(indices)) indices else list(selected = indices)
  if (!length(blocks) || is.null(names(blocks)) || anyNA(names(blocks)) ||
      any(!nzchar(names(blocks))) || anyDuplicated(names(blocks))) {
    stop("A list of indices must contain uniquely named, non-empty gene blocks.")
  }
  for (nm in names(blocks)) {
    x <- blocks[[nm]]
    if (is.character(x)) {
      if (anyNA(x) || any(!nzchar(x))) stop("Invalid feature ID in block '", nm, "'.")
      i <- match(x, feature_id)
      if (anyNA(i)) stop("Unknown feature ID in block '", nm, "': ",
                         paste(x[is.na(i)], collapse = ", "), ".")
    } else if (is.numeric(x) && is.null(dim(x))) {
      if (any(!is.finite(x)) || any(x != floor(x)) ||
          any(x < 1 | x > length(feature_id))) {
        stop("Numeric indices in block '", nm, "' must be valid integer feature positions.")
      }
      i <- as.integer(x)
    } else {
      stop("indices must contain feature IDs or integer feature positions, not component numbers or pair rows.")
    }
    if (length(i) < 2L || anyDuplicated(i)) {
      stop("Block '", nm, "' must contain at least two distinct features.")
    }
    blocks[[nm]] <- i
  }
  blocks
}

# Score-only native kernels for both mgcv backends. No per-gene R loop and no
# eigendecomposition fallback: every feature's score vector `a` is built by
# the C++ dense batch kernel with score_only = TRUE (no H, no pair-calibration
# matrices), in batches of 256 genes.
.mgcvst_wgcna_scores <- function(fit, used, verbose, threads = 1L) {
  group <- "global"
  legacy <- is.null(fit$geometry$smooth)
  native <- NULL
  T0 <- field_scale <- NULL
  if (legacy) {
    T0 <- .mgcvst_legacy_shared_score_factor(fit$geometry)
    field_scale <- .mgcvst_field_scale(fit)
    width <- stats::setNames(ncol(T0), "global")
  } else {
    fit$.mgcvst_fixed_factors <- .mgcvst_model_fixed_factors(fit)
    native <- .mgcvst_model_dense_preparation(fit, used)
    if (is.null(native)) {
      stop("mgcvST.wgcna() requires a single marked-SPDE model with nuisance covariance.")
    }
    width <- native$width
  }
  A <- matrix(NA_real_, unname(width), length(used),
             dimnames = list(NULL, fit$feature_id[used]))
  batch_size <- 256L
  first <- 1L
  while (first <= length(used)) {
    ids <- used[first:min(length(used), first + batch_size - 1L)]
    cols <- first:min(length(used), first + batch_size - 1L)
    if (legacy) {
      z <- mgcvst_dense_score_batch_cpp(
        T0, fit$working_variance[, ids, drop = FALSE],
        fit$working_error[, ids, drop = FALSE], field_scale[ids],
        fit$geometry$X, list(), threads, score_only = TRUE
      )
    } else {
      phi <- fit$dispersion[ids]
      sp <- fit$smoothing_parameters[ids, , drop = FALSE]
      bad <- !is.finite(phi) | phi <= 0 |
        rowSums(!is.finite(sp) | sp <= 0) > 0L
      z <- mgcvst_dense_score_batch_cpp(
        native$T0, fit$working_variance[, ids, drop = FALSE],
        fit$working_error[, ids, drop = FALSE],
        phi / sp[, native$sp_index], native$X,
        fit$nuisance_covariance[ids], threads, score_only = TRUE
      )
      for (k in seq_along(ids)) {
        if (bad[k]) z[[k]] <- list(error =
          "The feature has invalid dispersion or smoothing parameters.")
      }
    }
    for (k in seq_along(ids)) {
      if (!is.null(z[[k]]$error)) {
        stop("Feature '", fit$feature_id[ids[k]], "' failed WGCNA score ",
             "construction: ", z[[k]]$error)
      }
      A[, cols[k]] <- z[[k]]$a
    }
    if (verbose && (min(length(used), first + batch_size - 1L) %% 256L == 0L ||
        first + batch_size - 1L >= length(used))) {
      message("Constructed scores for ", min(length(used), first + batch_size - 1L),
              " of ", length(used), " features.")
    }
    first <- first + length(ids)
  }
  if (any(!is.finite(A))) {
    stop("The selected score coordinates contain non-finite values.")
  }
  list(A = A, group = group, width = width[group],
       feature_id = fit$feature_id[used])
}
# --------------------------------------------------------------------------
# Shared WGCNA front end and back end.
#
# mgcvST.wgcna() and inlaST.wgcna() differ in exactly one step: how the score
# matrix A is built. Everything before it (argument validation, gene blocks,
# working-model validity, optional-package check) and everything after it
# (similarity -> correlation -> adjacency -> TOM -> tree -> dynamic cut ->
# colours) is identical and lives here once.
# --------------------------------------------------------------------------

# Validate the arguments shared by both entry points and resolve gene blocks.
.mgcvst_wgcna_prepare <- function(fit, indices, wgcna.para, verbose, threads,
                                  caller) {
  if (!inherits(fit, "mgcvST_fit")) {
    stop(caller, "() requires a compact fit returned by mgcvST.estimate() or ",
         "inlaST.estimate().")
  }
  if (missing(indices) || is.null(indices)) {
    stop("indices must explicitly select the genes to analyze.")
  }
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("verbose must be TRUE or FALSE.")
  }
  threads <- as.integer(threads)
  if (length(threads) != 1L || is.na(threads) || threads < 1L) {
    stop("threads must be one positive integer.")
  }
  para <- .mgcvst_wgcna_parameters(wgcna.para)
  ids <- fit$feature_id
  if (!is.character(ids) || anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) {
    stop("The fit must have unique, non-empty feature IDs.")
  }
  blocks <- .mgcvst_wgcna_indices(indices, ids)
  used <- unique(unlist(blocks, use.names = FALSE))
  geometry <- fit$geometry
  if (is.null(geometry)) stop("The fit does not retain score geometry.")
  # One global spatial score process supplies the WGCNA coordinates.
  available <- if (!is.null(geometry$smooth)) names(geometry$target) else "global"
  if (!identical(available, "global")) {
    stop("The fit must carry exactly one score component named 'global'.")
  }
  if (.mgcvst_inla_downstream(fit)) {
    if (!is.matrix(fit$score_a) || ncol(fit$score_a) != length(ids) ||
        length(fit$dispersion) != length(ids)) {
      stop("The compact fit dimensions are incompatible with feature_id.")
    }
    valid <- is.finite(fit$dispersion[used]) & fit$dispersion[used] > 0 &
      colSums(!is.finite(fit$score_a[, used, drop = FALSE])) == 0L
  } else {
    E <- fit$working_error
    V <- fit$working_variance
    if (length(dim(E)) != 2L || length(dim(V)) != 2L ||
        !identical(dim(E), dim(V)) || ncol(E) != length(ids) ||
        length(fit$dispersion) != length(ids)) {
      stop("The compact fit dimensions are incompatible with feature_id.")
    }
    valid <- is.finite(fit$dispersion[used]) & fit$dispersion[used] > 0 &
      colSums(!is.finite(E[, used, drop = FALSE])) == 0L &
      colSums(!is.finite(V[, used, drop = FALSE]) | V[, used, drop = FALSE] <= 0) == 0L
  }
  if (any(!valid)) {
    stop("Selected features lack valid working models: ",
         paste(ids[used[!valid]], collapse = ", "), ".")
  }
  if (!requireNamespace("WGCNA", quietly = TRUE) ||
      !requireNamespace("dynamicTreeCut", quietly = TRUE) ||
      !requireNamespace("fastcluster", quietly = TRUE)) {
    stop("Install packages 'WGCNA', 'dynamicTreeCut', and 'fastcluster' to use ",
         caller, "().")
  }
  list(para = para, ids = ids, blocks = blocks, used = used, threads = threads)
}

# Similarity, adjacency, TOM, tree and module labels for every requested block.
# `score` is whatever the backend-specific score builder returned.
.mgcvst_wgcna_networks <- function(score, blocks, ids, para, verbose) {
  networks <- modules <- vector("list", length(blocks))
  names(networks) <- names(modules) <- names(blocks)
  for (nm in names(blocks)) {
    id <- ids[blocks[[nm]]]
    A <- score$A[, match(id, score$feature_id), drop = FALSE]
    normalization <- if (is.null(score$normalization)) nrow(A) else
      score$normalization
    S <- .magic_mm(A, A, transA = TRUE) / normalization
    dimnames(S) <- list(id, id)
    if (any(!is.finite(S)) || any(diag(S) <= 0)) {
      stop("Block '", nm, "' has an invalid score covariance.")
    }
    R <- stats::cov2cor(S)
    adj <- WGCNA::adjacency.fromSimilarity(R, type = para$networkType,
                                           power = para$power)
    TOM <- WGCNA::TOMsimilarity(adj, TOMType = para$TOMType, verbose = 0)
    dimnames(adj) <- dimnames(TOM) <- list(id, id)
    H <- fastcluster::hclust(stats::as.dist(1 - TOM), method = para$hclustMethod)
    status <- if (length(id) < para$minClusterSize) "below_minClusterSize" else
      "evaluated"
    labels <- if (status == "below_minClusterSize") integer(length(id)) else
      as.integer(dynamicTreeCut::cutreeDynamic(H, distM = 1 - TOM,
        minClusterSize = para$minClusterSize, deepSplit = para$deepSplit,
        verbose = 0))
    names(labels) <- id
    tab <- data.frame(feature_id = id, module = unname(labels),
                      color = WGCNA::labels2colors(labels),
                      stringsAsFactors = FALSE)
    networks[[nm]] <- list(feature_id = id, covariance = S, correlation = R,
      adjacency = adj, TOM = TOM, tree = H, labels = labels, modules = tab,
      q = nrow(A), status = status)
    modules[[nm]] <- data.frame(component = nm, tab, stringsAsFactors = FALSE)
    if (verbose) {
      message("Block '", nm, "': ", length(id), " genes, ",
              length(unique(labels[labels > 0L])), " modules, ",
              sum(labels == 0L), " grey; ", status, ".")
    }
  }
  tab <- do.call(rbind, modules)
  rownames(tab) <- NULL
  list(modules = tab, networks = networks)
}

# Assemble the public result. Shared so both entry points return one shape.
.mgcvst_wgcna_result <- function(score, prepared, verbose, started,
                                 score_seconds, call) {
  t0 <- proc.time()[["elapsed"]]
  z <- .mgcvst_wgcna_networks(score, prepared$blocks, prepared$ids,
                              prepared$para, verbose)
  structure(list(
    modules = z$modules, networks = z$networks, score = score,
    settings = list(
      group = "global",
      indices = lapply(prepared$blocks, function(i) prepared$ids[i]),
      wgcna.para = prepared$para
    ),
    timing = list(score_seconds = score_seconds,
                  network_seconds = proc.time()[["elapsed"]] - t0,
                  total_seconds = proc.time()[["elapsed"]] - started),
    call = call), class = "mgcvST_wgcna")
}

#' Identify co-expression modules within explicitly selected gene blocks
#'
#' Constructs the original score covariance `crossprod(A) / nrow(A)` from an
#' [mgcvST.estimate()] or [inlaST.estimate()] fit, converts it to correlation,
#' and performs WGCNA. `indices` selects genes, for example the genes in a
#' connected component
#' identified after [mgcvST.test()]. A named list analyzes several gene blocks
#' separately. It never pools blocks, selects significant pairs, or reruns
#' estimation or testing. All pairs within each selected block enter its COV.
#'
#' Shared SPDE factors are constructed once. All backends compute only the score
#' vectors needed here; pair-test calibration matrices are unnecessary. The
#' result retains these vectors and the matrices needed for downstream analysis
#' without modifying the supplied fit. No centering, ridge, or covariance
#' projection is applied.
#'
#' The score coordinates of the single spatial component are divided by their
#' coordinate count, as in the score covariance definition.
#'
#' An [inlaST.estimate()] fit is not accepted here; use [inlaST.wgcna()] for
#' sparse INLA fits.
#'
#' @param fitmgcvST A compact fit returned by [mgcvST.estimate()].
#' @param indices Required gene IDs, integer positions in `fitmgcvST$feature_id`,
#'   or a named list of such vectors. Each block must contain at least two distinct
#'   genes. A vector defines the block named `selected`. Gene order is preserved.
#' @param wgcna.para `NULL`, an empty list, or named partial overrides of:
#'   `networkType = "signed"`, `power = 6`, `TOMType = "signed"`,
#'   `hclustMethod = "average"`, `minClusterSize = 20L`, `deepSplit = 1L`.
#'   Unspecified settings keep these defaults; unknown names are errors.
#'   Blocks smaller than `minClusterSize` retain their matrices and receive
#'   grey (zero) labels, without changing the requested module size.
#' @param verbose Whether to display compact progress messages.
#' @param threads Positive OpenMP thread count controlling score construction
#'   for both backends.
#' @return An `mgcvST_wgcna` object with `modules` (component, feature ID, integer
#'   module, color), named `networks` (feature IDs, covariance, correlation,
#'   adjacency, TOM, tree, labels, modules, coordinate count `q`, and status),
#'   `score` (aligned `A`, group names and widths), `settings`, and `timing`.
#'   A zero module label means unassigned (grey). Module labels are local to
#'   each input block. No biological annotations are inferred automatically.
#' @seealso [inlaST.wgcna()] for sparse INLA fits.
#' @examples
#' \dontrun{
#' W <- mgcvST.wgcna(fit, indices = genes)
#' W <- mgcvST.wgcna(fit, indices = list(block1 = genes1, block2 = genes2),
#'                   wgcna.para = list(power = 6))
#' W$modules
#' W$networks$block1$covariance
#' }
#' @export
mgcvST.wgcna <- function(fitmgcvST, indices,
                         wgcna.para = NULL, verbose = FALSE, threads = 1L) {
  started <- proc.time()[["elapsed"]]
  call <- match.call()
  if (.mgcvst_inla_downstream(fitmgcvST)) {
    stop("mgcvST.wgcna() does not accept inlaST.estimate() fits; use inlaST.wgcna().")
  }
  prepared <- .mgcvst_wgcna_prepare(fitmgcvST, indices, wgcna.para, verbose,
                                    threads, "mgcvST.wgcna")
  t0 <- proc.time()[["elapsed"]]
  score <- .mgcvst_wgcna_scores(fitmgcvST, prepared$used, verbose,
                                threads = prepared$threads)
  score_seconds <- proc.time()[["elapsed"]] - t0
  .mgcvst_wgcna_result(score, prepared, verbose, started, score_seconds, call)
}

#' Identify co-expression modules from a sparse INLA fit
#'
#' The sparse-kernel sibling of [mgcvST.wgcna()]. It takes an
#' [inlaST.estimate()] fit and builds the gene-by-gene similarity
#' \eqn{S_{ij} = (R'a_i)'(R'a_j)}, where `R` is the observation-kernel
#' coordinate basis used by [inlaST.test()] (the same coverage constant), and
#' \eqn{a_i} are the sparse INLA score vectors. It then runs exactly the same
#' WGCNA splitting as [mgcvST.wgcna()]: `WGCNA::adjacency.fromSimilarity()`,
#' `WGCNA::TOMsimilarity()`, `fastcluster::hclust()`,
#' `dynamicTreeCut::cutreeDynamic()` and `WGCNA::labels2colors()`. The
#' similarity, the normaliser and the downstream code are shared with
#' [mgcvST.wgcna()], not re-derived.
#'
#' The score vectors are produced by projecting the sparse INLA score vectors
#' onto `R`; `threads` has no effect on this projection. BiocParallel is not
#' used, because the sparse INLA downstream runs one OpenMP layer in the
#' manager process.
#'
#' The normaliser is the basis rank `r` (the number of retained projection
#' coordinates), not `nrow(A)`.
#'
#' @inheritParams mgcvST.wgcna
#' @param fitmgcvST A compact fit returned by [inlaST.estimate()].
#' @return An `mgcvST_wgcna` object, identical in shape to the
#'   [mgcvST.wgcna()] result.
#' @seealso [mgcvST.wgcna()] for `mgcvST.estimate()` fits.
#' @examples
#' \dontrun{
#' fit <- inlaST.estimate(Y, model)
#' W <- inlaST.wgcna(fit, indices = genes, threads = 4L)
#' W$modules
#' }
#' @export
inlaST.wgcna <- function(fitmgcvST, indices,
                         wgcna.para = NULL, verbose = FALSE, threads = 1L) {
  started <- proc.time()[["elapsed"]]
  call <- match.call()
  if (!.mgcvst_inla_downstream(fitmgcvST)) {
    stop("inlaST.wgcna() requires a fit returned by inlaST.estimate(); use ",
         "mgcvST.wgcna() for an mgcvST.estimate()/model.set() fit.")
  }
  .mgcvst_inla_require_sparse(fitmgcvST)
  prepared <- .mgcvst_wgcna_prepare(fitmgcvST, indices, wgcna.para, verbose,
                                    threads, "inlaST.wgcna")
  t0 <- proc.time()[["elapsed"]]
  score <- .mgcvst_inla_wgcna_scores(fitmgcvST, prepared$used,
                                     prepared$threads, verbose)
  score_seconds <- proc.time()[["elapsed"]] - t0
  .mgcvst_wgcna_result(score, prepared, verbose, started, score_seconds, call)
}

#' @rdname mgcvST.wgcna
#' @param x An `mgcvST_wgcna` result.
#' @param ... Unused.
#' @export
print.mgcvST_wgcna <- function(x, ...) {
  cat("mgcvST WGCNA modules\n")
  cat("  score groups:", paste(x$settings$group, collapse = ", "), "\n")
  for (nm in names(x$networks)) {
    z <- x$networks[[nm]]
    cat(" ", nm, ":", length(z$feature_id), "genes;",
        length(unique(z$labels[z$labels > 0L])), "modules;",
        sum(z$labels == 0L), "grey;", z$status, "\n")
  }
  cat("  elapsed seconds:", format(x$timing$total_seconds), "\n")
  invisible(x)
}
