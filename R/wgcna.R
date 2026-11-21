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

# --------------------------------------------------------------------------
# Shared WGCNA front end and back end.
#
# Argument validation, gene blocks, working-model validity and the
# optional-package check; then similarity -> correlation -> adjacency -> TOM
# -> tree -> dynamic cut -> colours.
# --------------------------------------------------------------------------

# Validate the arguments shared by both entry points and resolve gene blocks.
.mgcvst_wgcna_prepare <- function(fit, indices, wgcna.para, verbose, threads,
                                  caller) {
  if (!inherits(fit, "spdeST_scan")) {
    stop(caller, "() requires a fit returned by spdeST.estimate().")
  }
  .mgcvst_check_fit_format(fit)
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
  if (!identical(names(geometry$target), "global")) {
    stop("The fit must carry exactly one score component named 'global'.")
  }
  if (!is.matrix(fit$score_a) || ncol(fit$score_a) != length(ids) ||
      length(fit$dispersion) != length(ids)) {
    stop("The compact fit dimensions are incompatible with feature_id.")
  }
  valid <- is.finite(fit$dispersion[used]) & fit$dispersion[used] > 0 &
    colSums(!is.finite(fit$score_a[, used, drop = FALSE])) == 0L
  if (any(!valid)) {
    unfitted <- if (is.null(fit$diagnostics$spatial_fitted)) rep(FALSE, length(used)) else
      !fit$diagnostics$spatial_fitted[used]
    if (any(!valid & unfitted)) {
      stop("Selected features have no spatial fit; add one with spdeST.estimate(): ",
           paste(ids[used[!valid & unfitted]], collapse = ", "), ".")
    }
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
    call = call), class = "spdeST_wgcna")
}

#' Identify co-expression modules from a fit
#'
#' Builds the gene-by-gene similarity \eqn{S_{ij} = (R'a_i)'(R'a_j)} of the
#' genes of each selected block, where `R` is the observation-kernel
#' coordinate basis used by [spdeST.test()] (all `q - 1` directions of the
#' constrained field) and \eqn{a_i} are the score vectors stored on the fit,
#' divided by the basis rank `r`. It then runs
#' `WGCNA::adjacency.fromSimilarity()`, `WGCNA::TOMsimilarity()`,
#' `fastcluster::hclust()`, `dynamicTreeCut::cutreeDynamic()` and
#' `WGCNA::labels2colors()`. Blocks are selected explicitly; no gene is
#' selected by estimation or testing, and all pairs within a block enter its
#' covariance. No centering, ridge, or covariance projection is applied, and
#' the supplied fit is not modified.
#'
#' @param fit A fit returned by [spdeST.estimate()].
#' @param indices Required gene IDs, integer positions in `fit$feature_id`,
#'   or a named list of such vectors. Each block must contain at least two distinct
#'   genes. A vector defines the block named `selected`. Gene order is preserved.
#' @param wgcna.para `NULL`, an empty list, or named partial overrides of:
#'   `networkType = "signed"`, `power = 6`, `TOMType = "signed"`,
#'   `hclustMethod = "average"`, `minClusterSize = 20L`, `deepSplit = 1L`.
#'   Unspecified settings keep these defaults; unknown names are errors.
#'   Blocks smaller than `minClusterSize` retain their matrices and receive
#'   grey (zero) labels, without changing the requested module size.
#' @param verbose Whether to display compact progress messages.
#' @param threads Positive OpenMP thread count controlling score construction.
#' @return An `spdeST_wgcna` object with `modules` (component, feature ID, integer
#'   module, color), named `networks` (feature IDs, covariance, correlation,
#'   adjacency, TOM, tree, labels, modules, coordinate count `q`, and status),
#'   `score` (aligned `A`, group names and widths), `settings`, and `timing`.
#'   A zero module label means unassigned (grey). Module labels are local to
#'   each input block. No biological annotations are inferred automatically.
#' @examples
#' \dontrun{
#' scan <- spdeST.scan(Y, model)
#' fit <- spdeST.estimate(scan, Y, features = "all")
#' W <- spdeST.wgcna(fit, indices = genes, threads = 4L)
#' W$modules
#' }
#' @export
spdeST.wgcna <- function(fit, indices,
                         wgcna.para = NULL, verbose = FALSE, threads = 1L) {
  started <- proc.time()[["elapsed"]]
  call <- match.call()
  prepared <- .mgcvst_wgcna_prepare(fit, indices, wgcna.para, verbose,
                                    threads, "spdeST.wgcna")
  t0 <- proc.time()[["elapsed"]]
  score <- .mgcvst_wgcna_scores(fit, prepared$used,
                                     prepared$threads, verbose)
  score_seconds <- proc.time()[["elapsed"]] - t0
  .mgcvst_wgcna_result(score, prepared, verbose, started, score_seconds, call)
}

#' @rdname spdeST.wgcna
#' @param x An `spdeST_wgcna` result.
#' @param ... Unused.
#' @export
print.spdeST_wgcna <- function(x, ...) {
  cat("spdeST WGCNA modules\n")
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
