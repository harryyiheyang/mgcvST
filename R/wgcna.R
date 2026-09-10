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

# Score vectors only: one Woodbury solve for the fixed design and error.
.mgcvst_wgcna_score_vector <- function(T, F, variance, X, error) {
  D <- 1 / variance
  DT <- D * T
  K <- diag(ncol(T)) + .magic_mm(T, DT, transA = TRUE)
  Y <- cbind(X, error)
  DY <- D * Y
  U <- DY - .magic_mm(DT, .magic_solve(K, .magic_mm(T, DY, transA = TRUE)))
  Pe <- U[, ncol(U), drop = FALSE]
  if (ncol(X)) {
    VX <- U[, seq_len(ncol(X)), drop = FALSE]
    H <- CppMatrix::matrixGeneralizedInverse(.magic_mm(X, VX, transA = TRUE))
    Pe <- Pe - .magic_mm(VX, .magic_mm(H, .magic_mm(X, Pe, transA = TRUE)))
  }
  as.numeric(.magic_mm(F, Pe, transA = TRUE))
}

# Reuse shared SPDE factors; retain all nuisance smoothers in each fitted V.
.mgcvst_wgcna_scores <- function(fit, used, group, verbose) {
  legacy <- is.null(fit$geometry$smooth)
  if (legacy) {
    Q <- fit$geometry$Q
    if (isTRUE(fit$geometry$score_precision_psd)) {
      Q <- as.matrix(Q)
      E <- CppMatrix::matrixEigen((Q + t(Q)) / 2)
      d <- as.numeric(E$values)
      tol <- sqrt(.Machine$double.eps) * max(1, max(abs(d)))
      if (min(d) < -tol) {
        stop("The shared score precision is not positive semidefinite.")
      }
      keep <- d > tol
      if (!any(keep)) {
        stop("The shared score precision has no positive eigenvalues.")
      }
      V <- as.matrix(E$vectors[, keep, drop = FALSE])
      base <- .magic_mm(
        fit$geometry$B,
        .magic_mm(sweep(V, 2L, 1 / sqrt(d[keep]), "*"), V, transB = TRUE)
      )
    } else {
      base <- .mgcvst_spde_factor(fit$geometry$B, Q, 1)
    }
    scale <- .mgcvst_field_scale(fit)
  } else {
    fit$.mgcvst_fixed_factors <- .mgcvst_model_fixed_factors(fit)
  }
  A <- NULL
  width <- NULL
  rows <- NULL
  for (k in seq_along(used)) {
    i <- used[k]
    if (legacy) {
      T <- F <- sqrt(scale[i]) * base
      a <- .mgcvst_wgcna_score_vector(
        T, F, fit$working_variance[, i], fit$geometry$X,
        fit$working_error[, i]
      )
      widths <- stats::setNames(ncol(F), "global")
    } else if (identical(fit$score_backend, "sparse")) {
      state <- .mgcvst_model_sparse_score_state(fit, i, score_only = TRUE)
      a <- state$a
      widths <- state$width
    } else {
      z <- .mgcvst_model_operator(fit, i)
      target <- z$target
      F <- do.call(cbind, target)
      Pe <- .mgcvst_model_apply_P(z$operator, fit$working_error[, i])
      a <- as.numeric(.magic_mm(F, matrix(Pe, ncol = 1L), transA = TRUE))
      widths <- vapply(target, ncol, integer(1L))
    }
    if (is.null(width)) {
      width <- widths
      if (is.null(names(width)) || anyNA(names(width)) ||
          any(!nzchar(names(width))) || anyDuplicated(names(width))) {
        stop("The fit does not have valid named score components.")
      }
      end <- cumsum(width)
      start <- end - width + 1L
      rows <- unlist(lapply(group, function(x) seq.int(start[[x]], end[[x]])),
                     use.names = FALSE)
      A <- matrix(NA_real_, length(rows), length(used),
                  dimnames = list(NULL, fit$feature_id[used]))
    } else if (!identical(widths, width)) {
      stop("Selected features do not share aligned score-coordinate groups.")
    }
    A[, k] <- a[rows]
    if (verbose && (k %% 100L == 0L || k == length(used))) {
      message("Constructed scores for ", k, " of ", length(used), " features.")
    }
  }
  if (any(!is.finite(A))) {
    stop("The selected score coordinates contain non-finite values.")
  }
  list(A = A, group = group, width = width[group],
       feature_id = fit$feature_id[used])
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
#' `group` selects spatial score components, not connected gene blocks. A single
#' available component is used automatically. With multiple components, select
#' their names explicitly; selected coordinate groups are concatenated and
#' divided by their total coordinate count, as in the score covariance definition.
#'
#' @param fitmgcvST A compact fit returned by [mgcvST.estimate()] or
#'   [inlaST.estimate()].
#' @param indices Required gene IDs, integer positions in `fitmgcvST$feature_id`,
#'   or a named list of such vectors. Each block must contain at least two distinct
#'   genes. A vector defines the block named `selected`. Gene order is preserved.
#' @param group Score-component names. `NULL` uses the sole available component;
#'   it is an error when the fit contains multiple score components.
#' @param wgcna.para `NULL`, an empty list, or named partial overrides of:
#'   `networkType = "signed"`, `power = 6`, `TOMType = "signed"`,
#'   `hclustMethod = "average"`, `minClusterSize = 20L`, `deepSplit = 1L`.
#'   Unspecified settings keep these defaults; unknown names are errors.
#'   Blocks smaller than `minClusterSize` retain their matrices and receive
#'   grey (zero) labels, without changing the requested module size.
#' @param verbose Whether to display compact progress messages.
#' @return An `mgcvST_wgcna` object with `modules` (component, feature ID, integer
#'   module, color), named `networks` (feature IDs, covariance, correlation,
#'   adjacency, TOM, tree, labels, modules, coordinate count `q`, and status),
#'   `score` (aligned `A`, group names and widths), `settings`, and `timing`.
#'   A zero module label means unassigned (grey). Module labels are local to
#'   each input block. No biological annotations are inferred automatically.
#' @examples
#' \dontrun{
#' W <- mgcvST.wgcna(fit, indices = genes)
#' W <- mgcvST.wgcna(fit, indices = list(block1 = genes1, block2 = genes2),
#'                   wgcna.para = list(power = 6))
#' W$modules
#' W$networks$block1$covariance
#' }
#' @export
mgcvST.wgcna <- function(fitmgcvST, indices, group = NULL,
                         wgcna.para = NULL, verbose = FALSE) {
  started <- proc.time()[["elapsed"]]
  if (!inherits(fitmgcvST, "mgcvST_fit")) {
    stop("fitmgcvST must be returned by mgcvST.estimate() or inlaST.estimate().")
  }
  if (missing(indices) || is.null(indices)) stop("indices must explicitly select the genes to analyze.")
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) stop("verbose must be TRUE or FALSE.")
  para <- .mgcvst_wgcna_parameters(wgcna.para)
  ids <- fitmgcvST$feature_id
  if (!is.character(ids) || anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) {
    stop("The fit must have unique, non-empty feature IDs.")
  }
  blocks <- .mgcvst_wgcna_indices(indices, ids)
  used <- unique(unlist(blocks, use.names = FALSE))
  geometry <- fitmgcvST$geometry
  if (is.null(geometry)) stop("The fit does not retain score geometry.")
  available <- if (!is.null(geometry$smooth)) names(geometry$target) else "global"
  if (!length(available) || anyNA(available) || any(!nzchar(available)) || anyDuplicated(available)) {
    stop("The fit does not have valid named score components.")
  }
  if (is.null(group)) {
    if (length(available) != 1L) stop("group must explicitly select score components: ", paste(available, collapse = ", "), ".")
    group <- available
  }
  if (!is.character(group) || !length(group) || anyNA(group) ||
      anyDuplicated(group) || any(!group %in% available)) {
    stop("group must contain distinct available score-component names: ", paste(available, collapse = ", "), ".")
  }
  E <- fitmgcvST$working_error
  V <- fitmgcvST$working_variance
  if (length(dim(E)) != 2L || length(dim(V)) != 2L ||
      !identical(dim(E), dim(V)) || ncol(E) != length(ids) ||
      length(fitmgcvST$dispersion) != length(ids)) stop("The compact fit dimensions are incompatible with feature_id.")
  valid <- is.finite(fitmgcvST$dispersion[used]) & fitmgcvST$dispersion[used] > 0 &
    colSums(!is.finite(E[, used, drop = FALSE])) == 0L &
    colSums(!is.finite(V[, used, drop = FALSE]) | V[, used, drop = FALSE] <= 0) == 0L
  if (any(!valid)) stop("Selected features lack valid working models: ", paste(ids[used[!valid]], collapse = ", "), ".")
  if (!requireNamespace("WGCNA", quietly = TRUE) || !requireNamespace("dynamicTreeCut", quietly = TRUE) ||
      !requireNamespace("fastcluster", quietly = TRUE)) {
    stop("Install packages 'WGCNA', 'dynamicTreeCut', and 'fastcluster' to use mgcvST.wgcna().")
  }
  t0 <- proc.time()[["elapsed"]]
  score <- .mgcvst_wgcna_scores(fitmgcvST, used, group, verbose)
  score_seconds <- proc.time()[["elapsed"]] - t0
  t0 <- proc.time()[["elapsed"]]
  networks <- modules <- vector("list", length(blocks))
  names(networks) <- names(modules) <- names(blocks)
  for (nm in names(blocks)) {
    id <- ids[blocks[[nm]]]
    A <- score$A[, match(id, score$feature_id), drop = FALSE]
    S <- .magic_mm(A, A, transA = TRUE) / nrow(A)
    dimnames(S) <- list(id, id)
    if (any(!is.finite(S)) || any(diag(S) <= 0)) stop("Block '", nm, "' has an invalid score covariance.")
    R <- stats::cov2cor(S)
    adj <- WGCNA::adjacency.fromSimilarity(R, type = para$networkType, power = para$power)
    TOM <- WGCNA::TOMsimilarity(adj, TOMType = para$TOMType, verbose = 0)
    dimnames(adj) <- dimnames(TOM) <- list(id, id)
    H <- fastcluster::hclust(stats::as.dist(1 - TOM), method = para$hclustMethod)
    status <- if (length(id) < para$minClusterSize) "below_minClusterSize" else "evaluated"
    labels <- if (status == "below_minClusterSize") integer(length(id)) else
      as.integer(dynamicTreeCut::cutreeDynamic(H, distM = 1 - TOM,
        minClusterSize = para$minClusterSize, deepSplit = para$deepSplit, verbose = 0))
    names(labels) <- id
    tab <- data.frame(feature_id = id, module = unname(labels),
                      color = WGCNA::labels2colors(labels), stringsAsFactors = FALSE)
    networks[[nm]] <- list(feature_id = id, covariance = S, correlation = R,
      adjacency = adj, TOM = TOM, tree = H, labels = labels, modules = tab,
      q = nrow(A), status = status)
    modules[[nm]] <- data.frame(component = nm, tab, stringsAsFactors = FALSE)
    if (verbose) message("Block '", nm, "': ", length(id), " genes, ",
      length(unique(labels[labels > 0L])), " modules, ", sum(labels == 0L), " grey; ", status, ".")
  }
  tab <- do.call(rbind, modules)
  rownames(tab) <- NULL
  structure(list(modules = tab, networks = networks, score = score,
    settings = list(group = group, indices = lapply(blocks, function(i) ids[i]), wgcna.para = para),
    timing = list(score_seconds = score_seconds,
      network_seconds = proc.time()[["elapsed"]] - t0,
      total_seconds = proc.time()[["elapsed"]] - started), call = match.call()),
    class = "mgcvST_wgcna")
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
