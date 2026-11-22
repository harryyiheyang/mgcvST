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

# Validate the arguments and resolve the gene blocks.
.mgcvst_wgcna_prepare <- function(test, indices, wgcna.para, verbose) {
  if (!inherits(test, "spdeST_test")) {
    stop("spdeST.wgcna() requires a result of spdeST.test().")
  }
  if (missing(indices) || is.null(indices)) {
    stop("indices must explicitly select the genes to analyze.")
  }
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("verbose must be TRUE or FALSE.")
  }
  if (!all(c("feature_index", "variance", "variance_at_floor") %in% names(test$features))) {
    stop("The test does not hold the field variance of its features; re-run ",
         "spdeST.test().")
  }
  para <- .mgcvst_wgcna_parameters(wgcna.para)
  ids <- test$feature_id
  if (!is.character(ids) || anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) {
    stop("The test must have unique, non-empty feature IDs.")
  }
  blocks <- .mgcvst_wgcna_indices(indices, ids)
  if (!requireNamespace("WGCNA", quietly = TRUE) ||
      !requireNamespace("dynamicTreeCut", quietly = TRUE) ||
      !requireNamespace("fastcluster", quietly = TRUE)) {
    stop("Install packages 'WGCNA', 'dynamicTreeCut', and 'fastcluster' to use ",
         "spdeST.wgcna().")
  }
  list(para = para, ids = ids, blocks = blocks,
       used = unique(unlist(blocks, use.names = FALSE)))
}

# The rows (i, j, score, t1) of the discovered pairs of the test whose two
# features are selected.
.mgcvst_wgcna_pairs <- function(test, used) {
  columns <- c("i", "j", "score", "t1", "log_q")
  q <- test$threshold$q_value
  if (!is.numeric(q) || length(q) != 1L || !is.finite(q)) {
    stop("The test does not hold its discovery threshold.")
  }
  keep <- function(z) {
    z[z$i %in% used & z$j %in% used & !is.na(z$log_q) & z$log_q <= log(q),
      columns[1:4], drop = FALSE]
  }
  if (is.data.frame(test$results)) return(keep(test$results))
  if (!length(test$shards) || !all(file.exists(test$shards))) {
    stop("The test holds neither its pair table nor its result shards.")
  }
  do.call(rbind, lapply(test$shards, function(file) {
    keep(.mgcvst_read_shard(file, columns))
  }))
}

# Correlation, adjacency, TOM, tree and module labels for every requested block.
.mgcvst_wgcna_networks <- function(pairs, features, blocks, ids, para, verbose) {
  networks <- modules <- vector("list", length(blocks))
  names(networks) <- names(modules) <- names(blocks)
  for (nm in names(blocks)) {
    index <- blocks[[nm]]
    id <- ids[index]
    position <- integer(length(ids))
    position[index] <- seq_along(index)
    at <- match(index, features$feature_index)
    z <- mgcvst_he_correlation_cpp(pairs$i, pairs$j, pairs$score, pairs$t1,
                                   position, as.numeric(features$variance[at]),
                                   length(index))
    at_floor <- id[which(features$variance_at_floor[at])]
    R <- z$rho
    dimnames(R) <- dimnames(z$covariance) <- list(id, id)
    adj <- WGCNA::adjacency.fromSimilarity(z$similarity, type = para$networkType,
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
    networks[[nm]] <- list(feature_id = id, correlation = R, pairs = z$pairs,
      limited = z$limited,
      covariance = z$covariance, variance_at_floor = at_floor,
      adjacency = adj, TOM = TOM, tree = H, labels = labels, modules = tab,
      status = status)
    modules[[nm]] <- data.frame(component = nm, tab, stringsAsFactors = FALSE)
    if (verbose) {
      message("Block '", nm, "': ", length(id), " genes, ",
              length(unique(labels[labels > 0L])), " modules, ",
              sum(labels == 0L), " grey; ", z$pairs, " discovered pairs, ",
              z$limited, " limited to [-1, 1]; ", length(at_floor),
              " genes with the variance at its lower limit; ", status, ".")
    }
  }
  tab <- do.call(rbind, modules)
  rownames(tab) <- NULL
  list(modules = tab, networks = networks)
}

#' Identify co-expression modules from the pair test
#'
#' The network is built from the discoveries of [spdeST.test()]. The
#' similarity of two genes whose pair is a discovery (`log_q <= log(q.value)`)
#' is the Haseman-Elston-type estimate of their latent spatial correlation,
#' `rho_ij = U_ij / t1_ij`, with the score `U_ij` and its information `t1_ij`
#' of the pair; the similarity of every other pair is 0, because the ratio of
#' a score and an information that are both near zero estimates nothing; a
#' gene has similarity 1 with itself. No score is recomputed. The estimate of
#' a discovered pair can still lie outside `[-1, 1]`; the adjacency is
#' defined on `[-1, 1]`, so such an estimate enters it as `-1` or `1`, and
#' the number of these pairs is returned. The function then runs
#' `WGCNA::adjacency.fromSimilarity()`, `WGCNA::TOMsimilarity()`,
#' `fastcluster::hclust()`, `dynamicTreeCut::cutreeDynamic()` and
#' `WGCNA::labels2colors()`. Blocks are selected explicitly.
#'
#' @param test A result of [spdeST.test()] over the pairs of the selected
#'   genes; `keep = "discoveries"` holds all that is read.
#' @param indices Required gene IDs, integer positions in `test$feature_id`,
#'   or a named list of such vectors. Each block must contain at least two distinct
#'   genes. A vector defines the block named `selected`. Gene order is preserved.
#' @param wgcna.para `NULL`, an empty list, or named partial overrides of:
#'   `networkType = "signed"`, `power = 6`, `TOMType = "signed"`,
#'   `hclustMethod = "average"`, `minClusterSize = 20L`, `deepSplit = 1L`.
#'   Unspecified settings keep these defaults; unknown names are errors.
#'   Blocks smaller than `minClusterSize` retain their matrices and receive
#'   grey (zero) labels, without changing the requested module size.
#' @param verbose Whether to display compact progress messages.
#' @return An `spdeST_wgcna` object with `modules` (component, feature ID, integer
#'   module, color), named `networks` (feature IDs, `correlation` (the
#'   estimates `U / t1` of the discovered pairs, not limited, and 0 for the
#'   other pairs), `pairs` (the number of discovered pairs of the block),
#'   `limited` (the number of them outside
#'   `[-1, 1]`), `covariance` (the estimates of the covariance of the latent
#'   fields, `U / t1` times the root of the product of the two field
#'   variances, with the variances on the diagonal), `variance_at_floor` (the
#'   genes whose spatial fit ended on the lower limit of the variance: their
#'   correlations are computed at that limit and are not estimates),
#'   adjacency, TOM, tree, labels, modules and status),
#'   `settings`, and `timing`.
#'   A zero module label means unassigned (grey). Module labels are local to
#'   each input block. No biological annotations are inferred automatically.
#' @examples
#' \dontrun{
#' scan <- spdeST.scan(Y, model)
#' fit <- spdeST.estimate(scan, Y, features = "all")
#' test <- spdeST.test(fit, moments = "exact", keep = "discoveries")
#' W <- spdeST.wgcna(test, indices = genes)
#' W$modules
#' }
#' @export
spdeST.wgcna <- function(test, indices, wgcna.para = NULL, verbose = FALSE) {
  started <- proc.time()[["elapsed"]]
  call <- match.call()
  prepared <- .mgcvst_wgcna_prepare(test, indices, wgcna.para, verbose)
  pairs <- .mgcvst_wgcna_pairs(test, prepared$used)
  pair_seconds <- proc.time()[["elapsed"]] - started
  t0 <- proc.time()[["elapsed"]]
  z <- .mgcvst_wgcna_networks(pairs, test$features, prepared$blocks, prepared$ids,
                              prepared$para, verbose)
  structure(list(
    modules = z$modules, networks = z$networks,
    settings = list(
      indices = lapply(prepared$blocks, function(i) prepared$ids[i]),
      wgcna.para = prepared$para
    ),
    timing = list(pair_seconds = pair_seconds,
                  network_seconds = proc.time()[["elapsed"]] - t0,
                  total_seconds = proc.time()[["elapsed"]] - started),
    call = call), class = "spdeST_wgcna")
}

#' @rdname spdeST.wgcna
#' @param x An `spdeST_wgcna` result.
#' @param ... Unused.
#' @export
print.spdeST_wgcna <- function(x, ...) {
  cat("spdeST WGCNA modules\n")
  for (nm in names(x$networks)) {
    z <- x$networks[[nm]]
    cat(" ", nm, ":", length(z$feature_id), "genes;",
        length(unique(z$labels[z$labels > 0L])), "modules;",
        sum(z$labels == 0L), "grey;", z$status, "\n")
  }
  cat("  elapsed seconds:", format(x$timing$total_seconds), "\n")
  invisible(x)
}
