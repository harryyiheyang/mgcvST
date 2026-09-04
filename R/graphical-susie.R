# Fit one nodewise SuSiE regression from Graph-SuSiE sufficient statistics.
.mgcvst_graphical_susie_node <- function(
    g, S, W, n, L, max_iter, susie_tol, residual_floor,
    coverage, min_abs_corr, model_init = NULL) {
  p <- ncol(S)
  idx <- setdiff(seq_len(p), g)
  yty <- n * S[g, g]
  observed_variance <- yty / (n - 1)
  if (!is.finite(observed_variance) || observed_variance <= 0) {
    stop("The Graph-SuSiE pseudo-response variance must be positive.")
  }
  negative_residual <- FALSE
  fit <- tryCatch(
    susieR::susie_ss(
      XtX = n * W[idx, idx, drop = FALSE],
      Xty = n * S[idx, g], yty = yty, n = n, L = L,
      X_colmeans = rep(0, length(idx)), y_mean = 0,
      standardize = TRUE, model_init = model_init,
      estimate_residual_variance = TRUE,
      residual_variance_lowerbound = residual_floor * observed_variance,
      estimate_prior_variance = TRUE,
      coverage = coverage, min_abs_corr = min_abs_corr,
      max_iter = max_iter, tol = susie_tol, check_input = FALSE,
      verbose = FALSE
    ),
    error = function(e) {
      if (!grepl(
        "est_residual_variance() failed: the estimated value is negative",
        conditionMessage(e), fixed = TRUE
      )) {
        stop(e)
      }
      negative_residual <<- TRUE
      NULL
    }
  )
  if (negative_residual) {
    return(list(
      node = g, predictors = idx, support = integer(), cs = list(),
      beta = numeric(length(idx)), pip = numeric(length(idx)),
      model_init = model_init, niter = 0L, converged = TRUE,
      residual_variance = NA_real_, negative_residual = TRUE
    ))
  }
  cs <- fit$sets$cs
  keep <- if (length(cs)) {
    sort(unique(unlist(cs, use.names = FALSE)))
  } else {
    integer()
  }
  beta <- susieR::susie_get_posterior_mean(fit)
  if (length(keep)) beta[-keep] <- 0 else beta[] <- 0
  init <- fit[c("alpha", "mu", "mu2", "V", "sigma2")]
  init <- structure(init, class = "susie")
  list(
    node = g, predictors = idx, support = idx[keep], cs = cs,
    beta = beta, pip = as.numeric(fit$pip), model_init = init,
    niter = fit$niter, converged = isTRUE(fit$converged),
    residual_variance = as.numeric(fit$sigma2), negative_residual = FALSE
  )
}

# Convert one feature's component widths to score-coordinate indexes.
.mgcvst_graphical_susie_group_index <- function(width, group) {
  if (is.null(names(width)) || any(!nzchar(names(width)))) {
    stop("The mgcvST score components do not have valid names.")
  }
  end <- cumsum(width)
  start <- end - width + 1L
  unlist(lapply(group, function(x) seq.int(start[[x]], end[[x]])),
         use.names = FALSE)
}

# Construct aligned score states while reusing fixed-kappa SPDE factors.
.mgcvst_graphical_susie_score_states <- function(fit, feature_index) {
  fit$.mgcvst_fixed_factors <- .mgcvst_model_fixed_factors(fit)
  lapply(feature_index, function(j) .mgcvst_model_score_state(fit, j))
}

#' Estimate a Graph-SuSiE network from an mgcvST fit
#'
#' Constructs score coordinates for explicitly selected spatial score
#' components, runs the first nodewise SuSiE scan in parallel, and uses the
#' purified credible-set AND support in one fixed-support `glasso` fit. The
#' returned covariance and precision from that fit are retained as `W0` and
#' `Theta0`. Subsequent sweeps visit all nodes in a newly randomized order and
#' update the full nodewise SuSiE regressions with warm `model_init` objects.
#' No further `glasso` refit is performed.
#'
#' `group` is deliberately required. A fit constructed with
#' `model.set(setting = "global_local")` contains separate `global` and
#' `local` score-coordinate groups; they are never combined by default.
#' Network edges are defined by purified SuSiE credible-set membership in both
#' nodewise directions.
#'
#' A negative residual-variance estimate defines a zero-edge result for that
#' node, and the remaining nodes continue. During serial sweeps, a proposed
#' node update with a non-positive conditional variance is not inserted into
#' the covariance matrix; that node retains its preceding state and the sweep
#' continues.
#'
#' @param fitmgcvST An `mgcvST_model_fit` returned by
#'   `mgcvST.estimate(Y, model.set(...))`.
#' @param group Non-empty character vector naming the score-component groups
#'   used to construct the network covariance. This argument has no default.
#' @param features Optional character feature IDs or integer feature indexes
#'   defining the network nodes. `NULL` uses every feature with a complete
#'   compact working model.
#' @param L Positive maximum number of SuSiE single effects per node. The
#'   effective value is at most the number of other network nodes.
#' @param max_iter Positive maximum number of IBSS iterations in each node fit.
#' @param max_sweeps Positive maximum number of randomized serial sweeps after
#'   the parallel initialization.
#' @param susie_tol Positive convergence tolerance passed to `susie_ss()`.
#' @param residual_floor Non-negative lower bound for the residual variance,
#'   expressed as a fraction of the observed pseudo-response variance.
#' @param coverage Credible-set coverage passed to `susie_ss()`.
#' @param min_abs_corr Minimum absolute within-set correlation used for SuSiE
#'   credible-set purification.
#' @param workers Positive number of workers used when `BPPARAM` is `NULL`.
#'   The default is 10 and is capped by the number of nodes.
#' @param BPPARAM Optional `BiocParallelParam` for the first scan. `NULL`
#'   selects a SOCK backend on Windows and a multicore backend elsewhere.
#' @param seed Integer seed controlling randomized sweep order.
#' @param glasso_tol Positive glasso convergence multiplier. It is used by the
#'   one fixed-support `glasso` initialization and by the subsequent coordinate
#'   sweeps with the glasso stopping rule: the mean absolute parameter change
#'   must be smaller than `glasso_tol` times the mean absolute off-diagonal
#'   score covariance.
#' @param glasso_maxit Positive maximum iterations for the fixed-support
#'   `glasso` initialization.
#' @param retain_model_init Logical; retain the minimal warm-start objects for
#'   all final node fits.
#' @param verbose Logical; print compact progress messages.
#' @return An object of class `mgcvST_network` containing the selected groups
#'   and features, score covariance, final covariance and precision,
#'   purified-CS AND adjacency, directed support and credible sets, the initial
#'   fixed-support `W0` and `Theta0`, convergence trace, timings, and settings.
#' @export
graphical_susie <- function(
    fitmgcvST, group, features = NULL, L = 5L, max_iter = 100L,
    max_sweeps = 10L, susie_tol = 1e-3,
    residual_floor = 0.1, coverage = 0.95, min_abs_corr = 0.5,
    workers = 10L, BPPARAM = NULL, seed = 1L,
    glasso_tol = 1e-4, glasso_maxit = 10000L,
    retain_model_init = TRUE, verbose = FALSE) {
  call <- match.call()
  if (!inherits(fitmgcvST, "mgcvST_model_fit")) {
    stop(
      "fitmgcvST must be returned by mgcvST.estimate(Y, model.set(...))."
    )
  }
  if (missing(group)) {
    stop("group must explicitly name at least one mgcvST score component.")
  }
  group <- as.character(group)
  if (!length(group) || anyNA(group) || any(!nzchar(group)) ||
      anyDuplicated(group)) {
    stop("group must contain unique non-empty score-component names.")
  }
  available <- as.character(fitmgcvST$score_components)
  unknown <- setdiff(group, available)
  if (length(unknown)) {
    stop(
      "Unknown score component(s): ", paste(unknown, collapse = ", "),
      ". Available components are: ", paste(available, collapse = ", "), "."
    )
  }
  feature_available <- .mgcvst_feature_available(fitmgcvST)
  if (is.null(features)) {
    feature_index <- which(feature_available)
  } else if (is.numeric(features)) {
    feature_index <- as.integer(features)
    if (anyNA(feature_index) || any(feature_index < 1L) ||
        any(feature_index > length(fitmgcvST$feature_id)) ||
        anyDuplicated(feature_index)) {
      stop("Numeric features must contain unique valid feature indexes.")
    }
  } else {
    features <- as.character(features)
    if (!length(features) || anyNA(features) || any(!nzchar(features)) ||
        anyDuplicated(features)) {
      stop("Character features must contain unique non-empty feature IDs.")
    }
    feature_index <- match(features, fitmgcvST$feature_id)
    if (anyNA(feature_index)) {
      stop("Every requested feature must occur in fitmgcvST$feature_id.")
    }
  }
  if (length(feature_index) < 2L) {
    stop("At least two features with complete working models are required.")
  }
  if (any(!feature_available[feature_index])) {
    stop("Every network feature must have a complete working model.")
  }
  feature_id <- fitmgcvST$feature_id[feature_index]
  p <- length(feature_index)
  n <- nrow(fitmgcvST$working_error)
  L <- as.integer(L)
  max_iter <- as.integer(max_iter)
  max_sweeps <- as.integer(max_sweeps)
  workers <- as.integer(workers)
  glasso_maxit <- as.integer(glasso_maxit)
  if (length(L) != 1L || is.na(L) || L < 1L) stop("L must be positive.")
  if (length(max_iter) != 1L || is.na(max_iter) || max_iter < 1L) {
    stop("max_iter must be positive.")
  }
  if (length(max_sweeps) != 1L || is.na(max_sweeps) || max_sweeps < 1L) {
    stop("max_sweeps must be positive.")
  }
  if (length(workers) != 1L || is.na(workers) || workers < 1L) {
    stop("workers must be positive.")
  }
  if (length(glasso_maxit) != 1L || is.na(glasso_maxit) ||
      glasso_maxit < 1L) {
    stop("glasso_maxit must be positive.")
  }
  scalar <- c(
    susie_tol = susie_tol, residual_floor = residual_floor, coverage = coverage,
    min_abs_corr = min_abs_corr, glasso_tol = glasso_tol
  )
  if (any(!is.finite(scalar))) stop("Numerical controls must be finite.")
  if (susie_tol <= 0 || residual_floor < 0 || coverage <= 0 ||
      coverage > 1 || min_abs_corr < 0 ||
      min_abs_corr > 1 || glasso_tol <= 0) {
    stop("One or more numerical controls are outside their valid range.")
  }
  if (!is.logical(retain_model_init) || length(retain_model_init) != 1L ||
      is.na(retain_model_init)) {
    stop("retain_model_init must be TRUE or FALSE.")
  }
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("verbose must be TRUE or FALSE.")
  }
  if (!is.null(BPPARAM) && !inherits(BPPARAM, "BiocParallelParam")) {
    stop("BPPARAM must be NULL or inherit from 'BiocParallelParam'.")
  }
  seed <- as.integer(seed)
  if (length(seed) != 1L || is.na(seed)) stop("seed must be one integer.")
  L_effective <- min(L, p - 1L)

  t0 <- proc.time()[["elapsed"]]
  state <- .mgcvst_graphical_susie_score_states(fitmgcvST, feature_index)
  score <- vector("list", p)
  coordinate_index <- NULL
  width0 <- NULL
  for (j in seq_len(p)) {
    if (is.null(width0)) {
      width0 <- state[[j]]$width
      coordinate_index <- .mgcvst_graphical_susie_group_index(width0, group)
    } else if (!identical(state[[j]]$width, width0)) {
      stop("Selected features do not share the same score-coordinate groups.")
    }
    score[[j]] <- state[[j]]$a[coordinate_index]
  }
  A_score <- do.call(cbind, score)
  colnames(A_score) <- feature_id
  S <- .magic_mm(A_score, A_score, transA = TRUE) / n
  S <- (S + t(S)) / 2
  dimnames(S) <- list(feature_id, feature_id)
  if (any(!is.finite(S)) || any(diag(S) <= 0)) {
    stop("The selected score-coordinate covariance is invalid.")
  }
  score_seconds <- proc.time()[["elapsed"]] - t0

  if (is.null(BPPARAM)) {
    nworkers <- min(workers, p)
    if (nworkers == 1L) {
      BPPARAM <- BiocParallel::SerialParam()
    } else if (.Platform$OS.type == "windows") {
      BPPARAM <- BiocParallel::SnowParam(
        workers = nworkers, type = "SOCK", stop.on.error = TRUE,
        progressbar = FALSE, RNGseed = seed
      )
    } else {
      BPPARAM <- BiocParallel::MulticoreParam(
        workers = nworkers, stop.on.error = TRUE,
        progressbar = FALSE, RNGseed = seed
      )
    }
  }
  nworkers <- min(p, BiocParallel::bpworkers(BPPARAM))
  if (verbose) message("Running the parallel Graph-SuSiE initialization.")
  t0 <- proc.time()[["elapsed"]]
  R0 <- BiocParallel::bplapply(
    seq_len(p), .mgcvst_graphical_susie_node,
    S = S, W = S, n = n, L = L_effective, max_iter = max_iter,
    susie_tol = susie_tol, residual_floor = residual_floor,
    coverage = coverage, min_abs_corr = min_abs_corr,
    model_init = NULL, BPPARAM = BPPARAM
  )
  parallel_seconds <- proc.time()[["elapsed"]] - t0
  D <- matrix(FALSE, p, p, dimnames = list(feature_id, feature_id))
  CS <- vector("list", p)
  names(CS) <- feature_id
  model_init <- vector("list", p)
  names(model_init) <- feature_id
  pip <- matrix(0, p, p, dimnames = list(feature_id, feature_id))
  for (z in R0) {
    if (length(z$support)) D[z$node, z$support] <- TRUE
    CS[[z$node]] <- z$cs
    model_init[[z$node]] <- z$model_init
    pip[z$node, z$predictors] <- z$pip
  }
  adjacency0 <- D & t(D)
  diag(adjacency0) <- FALSE

  zero <- which(!adjacency0 & upper.tri(adjacency0), arr.ind = TRUE)
  if (!nrow(zero)) zero <- NULL
  if (verbose) message("Running the fixed-support glasso initialization.")
  t0 <- proc.time()[["elapsed"]]
  fit0 <- glasso::glasso(
    s = S, rho = 0, zero = zero, thr = glasso_tol,
    maxit = glasso_maxit, penalize.diagonal = FALSE, trace = FALSE
  )
  glasso_seconds <- proc.time()[["elapsed"]] - t0
  if (fit0$errflag != 0) stop("The fixed-support glasso initialization failed.")
  W0 <- fit0$w
  Theta0 <- fit0$wi
  dimnames(W0) <- dimnames(Theta0) <- list(feature_id, feature_id)
  W <- W0
  Theta <- Theta0
  off_diagonal <- row(S) != col(S)
  convergence_scale <- mean(abs(S[off_diagonal]))
  convergence_threshold <- glasso_tol * convergence_scale

  trace <- data.frame(
    sweep = 0L, covariance_change = NA_real_,
    average_absolute_change = NA_real_,
    convergence_threshold = convergence_threshold,
    directed_edges = sum(D), edges = sum(adjacency0) / 2L,
    mean_ibss_iterations = mean(vapply(R0, `[[`, numeric(1L), "niter")),
    converged_node_fraction = mean(vapply(
      R0, `[[`, logical(1L), "converged"
    )),
    negative_residual_nodes = sum(vapply(
      R0, `[[`, logical(1L), "negative_residual"
    )),
    invalid_covariance_nodes = 0L,
    seconds = parallel_seconds + glasso_seconds
  )
  converged <- FALSE
  set.seed(seed)
  sweep_seconds <- 0
  node_order <- vector("list", max_sweeps)
  for (sweep in seq_len(max_sweeps)) {
    W_old <- W
    order <- sample.int(p)
    node_order[[sweep]] <- order
    niter <- numeric(p)
    node_converged <- logical(p)
    negative_residual <- logical(p)
    invalid_covariance <- logical(p)
    t0 <- proc.time()[["elapsed"]]
    for (g in order) {
      D_old <- D[g, ]
      z <- .mgcvst_graphical_susie_node(
        g, S, W, n, L_effective, max_iter, susie_tol,
        residual_floor, coverage, min_abs_corr, model_init[[g]]
      )
      D[g, ] <- FALSE
      if (length(z$support)) D[g, z$support] <- TRUE
      keep <- which(D[g, ] & D[, g])
      pos <- match(keep, z$predictors, nomatch = 0L)
      pos <- pos[pos > 0L]
      beta <- z$beta
      if (length(pos)) beta[-pos] <- 0 else beta[] <- 0
      w12 <- as.numeric(.magic_mm(
        W[z$predictors, z$predictors, drop = FALSE],
        matrix(beta, ncol = 1L)
      ))
      conditional_variance <- W[g, g] - as.numeric(crossprod(beta, w12))
      if (!is.finite(conditional_variance) || conditional_variance <= 0) {
        D[g, ] <- D_old
        niter[g] <- z$niter
        node_converged[g] <- TRUE
        negative_residual[g] <- z$negative_residual
        invalid_covariance[g] <- TRUE
        next
      }
      W[z$predictors, g] <- w12
      W[g, z$predictors] <- w12
      CS[[g]] <- z$cs
      model_init[[g]] <- z$model_init
      pip[g, ] <- 0
      pip[g, z$predictors] <- z$pip
      niter[g] <- z$niter
      node_converged[g] <- z$converged
      negative_residual[g] <- z$negative_residual
    }
    elapsed <- proc.time()[["elapsed"]] - t0
    sweep_seconds <- sweep_seconds + elapsed
    covariance_change <- norm(W - W_old, "F") / norm(W_old, "F")
    average_absolute_change <- mean(abs(W - W_old))
    adjacency <- D & t(D)
    diag(adjacency) <- FALSE
    Theta <- .magic_solve(W, diag(p))
    dimnames(Theta) <- list(feature_id, feature_id)
    trace <- rbind(
      trace,
      data.frame(
        sweep = sweep, covariance_change = covariance_change,
        average_absolute_change = average_absolute_change,
        convergence_threshold = convergence_threshold,
        directed_edges = sum(D), edges = sum(adjacency) / 2L,
        mean_ibss_iterations = mean(niter),
        converged_node_fraction = mean(node_converged),
        negative_residual_nodes = sum(negative_residual),
        invalid_covariance_nodes = sum(invalid_covariance),
        seconds = elapsed
      )
    )
    if (verbose) {
      message(
        "Completed sweep ", sweep, "; average absolute parameter change = ",
        format(average_absolute_change, digits = 4L), "."
      )
    }
    if (all(node_converged) &&
        average_absolute_change < convergence_threshold) {
      converged <- TRUE
      break
    }
  }
  adjacency <- D & t(D)
  diag(adjacency) <- FALSE
  if (!retain_model_init) model_init <- NULL
  structure(
    list(
      feature_id = feature_id, feature_index = feature_index,
      group = group, score_group_width = width0[group],
      score_covariance = S, covariance = W, precision = Theta,
      adjacency = adjacency, directed_support = D,
      credible_sets = CS, pip = pip,
      W0 = W0, Theta0 = Theta0, initial_adjacency = adjacency0,
      model_init = model_init, node_order = node_order[seq_len(sweep)],
      converged = converged, sweeps = sweep, trace = trace,
      settings = list(
        L = L, effective_L = L_effective, max_iter = max_iter,
        max_sweeps = max_sweeps, susie_tol = susie_tol,
        residual_floor = residual_floor,
        coverage = coverage, min_abs_corr = min_abs_corr,
        workers = nworkers, seed = seed,
        glasso_tol = glasso_tol, glasso_maxit = glasso_maxit,
        coordinate_convergence = paste0(
          "mean(abs(W-W_old)) < glasso_tol * ",
          "mean(abs(offdiag(score_covariance)))"
        ),
        symmetrization = "purified-CS AND"
      ),
      timing = list(
        score_coordinates = score_seconds,
        parallel_initialization = parallel_seconds,
        glasso_initialization = glasso_seconds,
        susie_sweeps = sweep_seconds,
        total = score_seconds + parallel_seconds + glasso_seconds +
          sweep_seconds
      ),
      call = call
    ),
    class = c("mgcvST_network", "mgcvST")
  )
}

#' Print a Graph-SuSiE network
#'
#' @param x An object returned by [graphical_susie()].
#' @param ... Unused.
#' @return `x`, invisibly.
#' @export
print.mgcvST_network <- function(x, ...) {
  cat("mgcvST Graph-SuSiE network\n")
  cat("  score group(s):", paste(x$group, collapse = ", "), "\n")
  cat("  features:", length(x$feature_id), "\n")
  cat("  edges:", sum(x$adjacency) / 2L, "\n")
  cat("  L:", x$settings$effective_L, "\n")
  cat("  randomized sweeps:", x$sweeps, "\n")
  cat("  converged:", x$converged, "\n")
  cat("  elapsed seconds:", format(x$timing$total), "\n")
  invisible(x)
}
