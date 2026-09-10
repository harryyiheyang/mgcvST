# Internal INLA estimation engine.  The public adapter lives in the API file;
# keeping the engine independent of mgcv objects also makes its scaling and
# constraint conventions directly testable.

.inlast_merge_control <- function(base = list(), override = list()) {
  validate <- function(x, label) {
    if (is.null(x)) x <- list()
    if (!is.list(x)) stop(label, " must be a list.")
    nm <- names(x)
    if (length(x) && (is.null(nm) || anyNA(nm) || any(!nzchar(nm)) ||
                      anyDuplicated(nm))) {
      stop(label, " must have unique, non-empty names.")
    }
    x
  }
  base <- validate(base, "base control")
  override <- validate(override, "override control")
  # Prior fields are always replaced as one object. Native tuning fields can
  # be overlaid individually, while explicit NULL values remain explicit.
  out <- base
  for (nm in names(override)) {
    if (identical(nm, "control.inla") &&
        is.list(out[[nm]]) && is.list(override[[nm]])) {
      out[[nm]] <- .inlast_merge_control(out[[nm]], override[[nm]])
    } else {
      out[nm] <- list(override[[nm]])
    }
  }
  out
}

.inlast_control <- function(control = list()) {
  if (is.null(control)) control <- list()
  if (!is.list(control)) stop("control must be a list.")
  defaults <- list(
    int_strategy = "eb",
    latent_strategy = "gaussian",
    fixed_effect_precision = 0,
    precision_prior = list(
      prior = "flat", param = numeric(), initial = 0
    ),
    gaussian_precision = NULL,
    gaussian_precision_prior = list(
      prior = "normal", param = c(0, 1 / 9), initial = 0
    ),
    nb_size = NULL,
    nb_size_prior = list(
      prior = "flat", param = numeric(), initial = 0
    ),
    fixed_precision = NULL,
    control.inla = list(),
    num_threads = 1L,
    verbose = FALSE,
    keep_fit = FALSE
  )
  unknown <- setdiff(names(control), names(defaults))
  if (length(unknown)) {
    stop("Unknown INLA control field", if (length(unknown) > 1L) "s" else "",
         ": ", paste(unknown, collapse = ", "), ".")
  }
  # The public API validates once before parallel dispatch and the worker
  # validates again.  Retain explicit NULL controls so the second pass is
  # idempotent and cannot trigger partial `$` matching (for example nb_size
  # accidentally matching nb_size_prior).
  out <- .inlast_merge_control(defaults, control)
  if (!identical(out$int_strategy, "eb") ||
      !identical(out$latent_strategy, "gaussian")) {
    stop("The current engine requires int_strategy = 'eb' and latent_strategy = 'gaussian'.")
  }
  out$fixed_effect_precision <- as.numeric(out$fixed_effect_precision)
  if (length(out$fixed_effect_precision) != 1L ||
      !is.finite(out$fixed_effect_precision) ||
      out$fixed_effect_precision < 0) {
    stop("control$fixed_effect_precision must be one non-negative number.")
  }
  out$num_threads <- as.integer(out$num_threads)
  if (length(out$num_threads) != 1L || is.na(out$num_threads) ||
      out$num_threads < 1L) {
    stop("control$num_threads must be one positive integer.")
  }
  for (nm in c("verbose", "keep_fit")) {
    if (!is.logical(out[[nm]]) || length(out[[nm]]) != 1L ||
        is.na(out[[nm]])) stop("control$", nm, " must be TRUE or FALSE.")
  }
  for (nm in c("gaussian_precision", "nb_size")) {
    if (!is.null(out[[nm]])) {
      out[[nm]] <- as.numeric(out[[nm]])
      if (length(out[[nm]]) != 1L || !is.finite(out[[nm]]) ||
          out[[nm]] <= 0) stop("control$", nm, " must be positive when supplied.")
    }
  }
  # Validate prior structure in the parent process before feature-level
  # parallelism, including priors for parameters fixed in a particular model.
  for (nm in c("precision_prior", "gaussian_precision_prior",
               "nb_size_prior")) {
    validated <- .inlast_prior_hyper(out[[nm]])
    out[[nm]] <- validated[c("prior", "param", "initial")]
  }
  native <- out$control.inla
  if (!is.list(native) || (length(native) &&
      (is.null(names(native)) || anyNA(names(native)) ||
       any(!nzchar(names(native))) || anyDuplicated(names(native))))) {
    stop("control$control.inla must be a uniquely named list.")
  }
  allowed_native <- c(
    "strategy", "int.strategy",
    "h", "dz", "diff.logdens",
    "tolerance", "tolerance.f", "tolerance.g", "tolerance.x",
    "tolerance.step", "restart", "cutoff",
    "adapt.hessian.max.trials", "adapt.hessian.scale", "adaptive.max",
    "step.len", "stencil", "diagonal",
    "numint.maxfeval", "numint.relerr", "numint.abserr", "cmin",
    "step.factor", "global.node.factor", "global.node.degree",
    "stupid.search.max.iter", "stupid.search.factor",
    "constr.marginal.diagonal"
  )
  unsupported <- setdiff(names(native), allowed_native)
  if (length(unsupported)) {
    stop("Unsupported control$control.inla field",
         if (length(unsupported) > 1L) "s" else "", ": ",
         paste(unsupported, collapse = ", "), ".")
  }
  numeric_native <- setdiff(names(native), c("strategy", "int.strategy"))
  for (nm in numeric_native) {
    value <- native[[nm]]
    if (is.null(value)) next
    if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
        (!is.finite(value) && !identical(nm, "cmin"))) {
      stop("control$control.inla$", nm,
           " must be one finite numeric value.")
    }
  }
  positive_native <- intersect(
    names(native),
    c(
      "h", "dz", "diff.logdens",
      "tolerance", "tolerance.f", "tolerance.g", "tolerance.x",
      "tolerance.step", "adapt.hessian.scale",
      "numint.maxfeval", "numint.relerr", "numint.abserr",
      "global.node.factor", "stupid.search.factor",
      "constr.marginal.diagonal"
    )
  )
  for (nm in positive_native) {
    value <- native[[nm]]
    if (!is.null(value) && value <= 0) {
      stop("control$control.inla$", nm, " must be positive.")
    }
  }
  native_strategy <- native[["strategy", exact = TRUE]]
  native_integration <- native[["int.strategy", exact = TRUE]]
  if (!is.null(native_strategy) &&
      !identical(native_strategy, out$latent_strategy)) {
    stop("control$control.inla$strategy must be 'gaussian'.")
  }
  if (!is.null(native_integration) &&
      !identical(native_integration, out$int_strategy)) {
    stop("control$control.inla$int.strategy must be 'eb'.")
  }
  forced_native <- .inlast_merge_control(
    native,
    list(strategy = out$latent_strategy, int.strategy = out$int_strategy)
  )
  if (!requireNamespace("INLA", quietly = TRUE)) {
    stop("Validating control$control.inla requires the INLA package.")
  }
  builder <- get("control.inla", envir = asNamespace("INLA"))
  tryCatch(
    do.call(builder, forced_native),
    error = function(e) stop(
      "Invalid control$control.inla: ", conditionMessage(e), call. = FALSE
    )
  )
  out$control.inla <- native
  out
}

.inlast_prior_hyper <- function(x, fixed = NULL) {
  if (!is.list(x) || is.null(names(x)) || anyNA(names(x)) ||
      any(!nzchar(names(x))) || anyDuplicated(names(x))) {
    stop("An INLA hyperparameter prior must be a uniquely named list.")
  }
  unknown <- setdiff(names(x), c("prior", "param", "initial"))
  if (length(unknown)) {
    if ("fixed" %in% unknown) {
      stop("Do not put fixed in a prior object; use fixed_precision, ",
           "gaussian_precision, or nb_size.")
    }
    stop("Unknown INLA hyperparameter prior field",
         if (length(unknown) > 1L) "s" else "", ": ",
         paste(unknown, collapse = ", "), ".")
  }
  if (!is.character(x$prior) || length(x$prior) != 1L ||
      is.na(x$prior) || !nzchar(trimws(x$prior))) {
    stop("An INLA hyperparameter prior requires one non-empty prior name.")
  }
  initial <- if (is.null(x$initial)) 0 else as.numeric(x$initial)
  if (length(initial) != 1L || !is.finite(initial)) {
    stop("An INLA hyperparameter initial value must be finite.")
  }
  prior_text <- trimws(x$prior)
  prior_key <- tolower(prior_text)
  custom <- if (startsWith(prior_key, "expression:")) {
    "expression:"
  } else if (startsWith(prior_key, "table:")) {
    "table:"
  } else if (startsWith(prior_key, "rprior:")) {
    "rprior:"
  } else NULL
  if (is.null(x$param)) {
    if (identical(prior_key, "flat") ||
        (!is.null(custom) && custom %in% c("expression:", "table:"))) {
      param <- numeric()
    } else {
      stop("INLA prior '", prior_text, "' requires an explicit param vector.")
    }
  } else {
    param <- as.numeric(x$param)
  }
  if (any(!is.finite(param))) {
    stop("INLA hyperparameter prior parameters must be finite.")
  }

  if (!is.null(custom)) {
    if (identical(custom, "rprior:")) {
      stop("INLA rprior definitions are not supported by this parallel backend.")
    }
    definition <- substring(prior_text, nchar(custom) + 1L)
    if (!nzchar(trimws(definition))) {
      stop("An INLA ", custom, " prior requires a non-empty definition.")
    }
    if (length(param)) {
      stop("INLA ", custom,
           " priors encode parameters in the definition; param must be numeric().")
    }
    canonical <- prior_text
  } else {
    if (!requireNamespace("INLA", quietly = TRUE)) {
      stop("Validating an INLA hyperparameter prior requires the INLA package.")
    }
    registry <- INLA::inla.models()$prior
    registered <- names(registry)
    hit <- match(prior_key, tolower(registered))
    if (is.na(hit)) {
      stop("Unknown INLA hyperparameter prior '", prior_text, "'.")
    }
    canonical <- registered[hit]
    if (canonical %in% c("none", "invalid")) {
      stop("INLA prior '", canonical,
           "' is not an estimable hyperparameter prior.")
    }
    nparameters <- as.integer(registry[[canonical]]$nparameters)
    if (length(nparameters) != 1L || is.na(nparameters) ||
        nparameters < 0L) {
      stop("INLA prior '", canonical,
           "' has a variable parameter layout that this scalar backend ",
           "cannot validate.")
    }
    if (length(param) != nparameters) {
      stop("INLA prior '", canonical, "' requires ", nparameters,
           " parameter", if (nparameters == 1L) "" else "s", ".")
    }
  }
  if (prior_key %in% c("normal", "gaussian") && param[2L] <= 0) {
    stop("A normal INLA hyperparameter prior requires positive precision.")
  }
  if (prior_key %in% c("loggamma", "gamma") && any(param <= 0)) {
    stop("A log-Gamma INLA hyperparameter prior requires positive shape and rate.")
  }
  if (identical(prior_key, "flat") && length(param)) {
    stop("An INLA flat hyperparameter prior requires param = numeric().")
  }
  ans <- list(
    prior = canonical, param = param, initial = initial, fixed = FALSE
  )
  if (!is.null(fixed)) {
    fixed <- as.numeric(fixed)
    if (length(fixed) != 1L || !is.finite(fixed) || fixed <= 0) {
      stop("A fixed precision must be one positive finite number.")
    }
    ans$initial <- log(fixed)
    ans$fixed <- TRUE
  }
  ans
}

.inlast_prior_metadata <- function(ctl, family, fixed_precision) {
  describe <- function(prior, fixed = FALSE, value = NULL) {
    if (fixed) return(list(type = "fixed", value = value))
    # ctl has already been normalized by .inlast_control().
    key <- tolower(prior$prior)
    if (identical(key, "flat")) {
      list(
        type = "improper_flat_log_hyperparameter",
        prior = "flat", param = numeric(), initial = prior$initial,
        proper = FALSE
      )
    } else if (key %in% c("normal", "gaussian")) {
      mean <- prior$param[1L]
      sd <- sqrt(1 / prior$param[2L])
      list(
        type = "log_gaussian", prior = prior$prior,
        param = prior$param, initial = prior$initial, proper = TRUE,
        statement = sprintf(
          "log(parameter) ~ N(%s, %s^2)",
          format(mean, digits = 15L), format(sd, digits = 15L)
        )
      )
    } else {
      list(
        type = "INLA_hyperprior", prior = prior$prior,
        param = prior$param, initial = prior$initial,
        proper = if (key %in% c("loggamma", "gamma") ||
                     startsWith(key, "pc")) TRUE else NA
      )
    }
  }
  latent <- describe(
    ctl$precision_prior, fixed = !is.null(fixed_precision),
    value = fixed_precision
  )
  observation <- switch(
    family,
    gaussian = describe(
      ctl$gaussian_precision_prior,
      fixed = !is.null(ctl$gaussian_precision),
      value = ctl$gaussian_precision
    ),
    negative_binomial = describe(
      ctl$nb_size_prior, fixed = !is.null(ctl$nb_size), value = ctl$nb_size
    ),
    poisson = list(type = "none")
  )
  active <- c(list(latent), list(observation))
  improper <- any(vapply(active, function(x) identical(x$proper, FALSE), logical(1L)))
  list(latent_precision = latent, observation = observation,
       any_improper = improper)
}

.inlast_validate_spec <- function(spec, y, offset) {
  if (!is.list(spec) || is.null(spec$fixed) || is.null(spec$random)) {
    stop("spec must contain fixed and random model blocks.")
  }
  X <- as.matrix(spec$fixed$X)
  storage.mode(X) <- "double"
  n <- length(y)
  if (length(dim(X)) != 2L || nrow(X) != n || any(!is.finite(X))) {
    stop("spec$fixed$X must be a finite matrix with one row per observation.")
  }
  if (ncol(X) && qr(X)$rank < ncol(X)) {
    stop("spec$fixed$X must have full column rank.")
  }
  xnames <- spec$fixed$names
  if (is.null(xnames)) xnames <- colnames(X)
  if (is.null(xnames)) {
    xnames <- if (ncol(X)) paste0("fixed", seq_len(ncol(X))) else character()
  }
  if (length(xnames) != ncol(X) || anyNA(xnames) ||
      any(!nzchar(xnames)) || anyDuplicated(xnames)) {
    stop("spec$fixed$names must uniquely name every fixed-effect column.")
  }
  if (!is.list(spec$random) || !length(spec$random)) {
    stop("spec$random must contain at least one latent block.")
  }
  rnames <- vapply(spec$random, function(z) {
    if (is.null(z$name)) "" else as.character(z$name)[1L]
  }, character(1L))
  if (any(!nzchar(rnames)) || anyNA(rnames) || anyDuplicated(rnames)) {
    stop("Every random block must have a unique non-empty name.")
  }
  random <- vector("list", length(spec$random))
  constraints <- vector("list", length(spec$random))
  for (j in seq_along(spec$random)) {
    z <- spec$random[[j]]
    A <- methods::as(z$A, "CsparseMatrix")
    Q <- methods::as(z$Q, "CsparseMatrix")
    if (nrow(A) != n || !ncol(A) || nrow(Q) != ncol(A) ||
        ncol(Q) != ncol(A) || any(!is.finite(A@x)) ||
        any(!is.finite(Q@x))) {
      stop("Random block '", rnames[j], "' has incompatible A and Q matrices.")
    }
    if (!isTRUE(Matrix::isSymmetric(Q, tol = 1e-10))) {
      stop("Random block '", rnames[j], "' requires a symmetric Q matrix.")
    }
    kind <- if (is.null(z$kind)) "nuisance" else as.character(z$kind)[1L]
    if (!(kind %in% c("spde", "nuisance"))) {
      stop("Random block kind must be 'spde' or 'nuisance'.")
    }
    positive_definite <- tryCatch({
      suppressWarnings(Matrix::Cholesky(Matrix::forceSymmetric(Q), LDL = FALSE))
      TRUE
    }, error = function(e) FALSE)
    rankdef <- 0L
    if (!positive_definite && identical(kind, "spde")) {
      stop("Random block '", rnames[j],
           "' requires a positive-definite generic0 Q matrix.")
    }
    if (!positive_definite) {
      rankdef <- ncol(Q) - as.integer(Matrix::rankMatrix(Q)[1L])
      if (rankdef < 1L || rankdef >= ncol(Q)) {
        stop("Nuisance block '", rnames[j],
             "' must have a nonzero positive-semidefinite penalty.")
      }
    }
    if (!is.null(z$rankdef)) {
      supplied_rankdef <- as.integer(z$rankdef)
      if (length(supplied_rankdef) != 1L || is.na(supplied_rankdef) ||
          supplied_rankdef != rankdef) {
        stop("Random block '", rnames[j], "' has an incorrect rankdef.")
      }
    }
    constraint <- NULL
    if (identical(kind, "spde")) {
      # This is deliberately recomputed, rather than merely trusting adapter
      # metadata.  Thus every spatial fit uses the observation mean constraint.
      constraint <- as.numeric(Matrix::crossprod(A, rep.int(1 / n, n)))
      if (!any(abs(constraint) > 0) || any(!is.finite(constraint))) {
        stop("SPDE block '", rnames[j], "' has an invalid observation constraint.")
      }
      if (!is.null(z$constraint) && !isTRUE(all.equal(
          as.numeric(z$constraint), constraint, tolerance = 1e-10,
          check.attributes = FALSE))) {
        stop("SPDE block '", rnames[j],
             "' constraint does not equal crossprod(A, 1 / n).")
      }
    } else if (!is.null(z$constraint)) {
      constraint <- as.numeric(z$constraint)
      if (length(constraint) != ncol(A) || any(!is.finite(constraint)) ||
          !any(abs(constraint) > 0)) {
        stop("Random block '", rnames[j], "' has an invalid constraint.")
      }
    }
    z$A <- A
    z$Q <- Matrix::forceSymmetric(Q)
    precision_scale <- z[["precision_scale", exact = TRUE]]
    if (is.null(precision_scale)) precision_scale <- 1
    precision_scale <- as.numeric(precision_scale)
    if (length(precision_scale) != 1L || !is.finite(precision_scale) ||
        precision_scale <= 0) {
      stop("Random block '", rnames[j],
           "' requires one positive finite precision_scale.")
    }
    z$precision_scale <- precision_scale
    z$name <- rnames[j]
    z$kind <- kind
    z$rankdef <- rankdef
    z$target <- isTRUE(z$target)
    random[[j]] <- z
    constraints[j] <- list(constraint)
  }
  family <- tolower(as.character(spec$family)[1L])
  if (family %in% c("nb", "negative binomial", "negative_binomial")) {
    family <- "negative_binomial"
  }
  if (!(family %in% c("gaussian", "poisson", "negative_binomial"))) {
    stop("spec$family must be 'gaussian', 'poisson', or 'negative_binomial'.")
  }
  if (family %in% c("poisson", "negative_binomial") &&
      (any(y < 0) || any(abs(y - round(y)) > 1e-8))) {
    stop("Poisson and negative-binomial responses must be non-negative counts.")
  }
  list(
    n = n, X = X, xnames = xnames, random = random,
    constraints = constraints, family = family, offset = offset
  )
}

.inlast_hyper_mode <- function(fit, pattern) {
  # fit$mode$theta is the joint EB hyperparameter mode on INLA's internal
  # (log) scale.  summary.hyperpar contains transformed marginal modes, which
  # need not be the parameters used for the conditional latent fit.
  theta <- fit$mode$theta
  tags <- fit$mode$theta.tags
  if (is.null(tags) || length(tags) != length(theta)) tags <- names(theta)
  if (is.null(tags)) return(NA_real_)
  hit <- grep(pattern, tags, ignore.case = TRUE)
  if (length(hit) != 1L) return(NA_real_)
  value <- exp(as.numeric(theta[hit]))
  if (length(value) != 1L || !is.finite(value) || value <= 0) NA_real_ else value
}

.inlast_latent_mode_block <- function(fit, tag, expected_length) {
  contents <- fit$misc$configs$contents
  mode <- fit$mode$x
  if (is.null(contents) || is.null(contents$tag) ||
      is.null(contents$start) || is.null(contents$length) || is.null(mode)) {
    stop("INLA did not retain the joint latent-mode index.")
  }
  hit <- which(contents$tag == tag)
  if (length(hit) != 1L || contents$length[hit] != expected_length) {
    stop("INLA joint latent mode has an invalid block for '", tag, "'.")
  }
  index <- contents$start[hit] + seq_len(expected_length) - 1L
  if (any(index < 1L) || any(index > length(mode))) {
    stop("INLA joint latent mode index is out of bounds for '", tag, "'.")
  }
  value <- as.numeric(mode[index])
  if (length(value) != expected_length || any(!is.finite(value))) {
    stop("INLA joint latent mode is non-finite for '", tag, "'.")
  }
  value
}

.inlast_expected_covariance <- function(X, random, tau, working_variance,
                                        constraints, nuisance_index,
                                        fixed_effect_precision = 0) {
  p_fixed <- ncol(X)
  designs <- c(list(Matrix::Matrix(X, sparse = TRUE)),
               lapply(random, `[[`, "A"))
  design <- do.call(cbind, designs)
  inv_var <- 1 / as.numeric(working_variance)
  weighted <- design * sqrt(inv_var)
  penalties <- c(
    list(Matrix::Diagonal(p_fixed, fixed_effect_precision)),
    lapply(seq_along(random), function(j) tau[j] * random[[j]]$Q)
  )
  H <- Matrix::forceSymmetric(Matrix::crossprod(weighted) +
                              Matrix::bdiag(penalties))
  starts <- p_fixed + c(0L, cumsum(vapply(random, function(z) ncol(z$A),
                                          integer(1L))))
  active <- which(vapply(constraints, Negate(is.null), logical(1L)))
  C <- Matrix::Matrix(0, nrow = length(active), ncol = ncol(design), sparse = TRUE)
  if (length(active)) {
    for (k in seq_along(active)) {
      j <- active[k]
      cols <- starts[j] + seq_len(ncol(random[[j]]$A))
      C[k, cols] <- constraints[[j]]
    }
  }
  nuisance_index <- as.integer(nuisance_index)
  if (!length(nuisance_index)) return(matrix(numeric(), 0L, 0L))
  if (anyNA(nuisance_index) || any(nuisance_index < 1L) ||
      any(nuisance_index > ncol(design)) || anyDuplicated(nuisance_index)) {
    stop("spec$nuisance_index is invalid for the combined coefficient vector.")
  }
  rhs <- Matrix::sparseMatrix(
    i = nuisance_index, j = seq_along(nuisance_index), x = 1,
    dims = c(ncol(design), length(nuisance_index))
  )
  if (length(active)) rhs <- cbind(rhs, Matrix::t(C))
  solved <- Matrix::solve(H, rhs)
  Vn <- solved[, seq_along(nuisance_index), drop = FALSE]
  if (length(active)) {
    HC <- solved[, length(nuisance_index) + seq_len(nrow(C)), drop = FALSE]
    middle <- as.matrix(C %*% HC)
    Vn <- Vn - HC %*% solve(middle, as.matrix(C %*% Vn))
  }
  V <- as.matrix(Vn[nuisance_index, , drop = FALSE])
  (V + t(V)) / 2
}

# Fit one response using fixed-kappa generic0 SPDE blocks.  A generic0 precision
# is tau * Q.  Consequently the mgcv smoothing parameter on the package's
# working scale is lambda = phi * tau (phi is one for negative binomial).
.inlast_fit_feature <- function(spec, y, offset = NULL, control = list(),
                                diagnostics = FALSE) {
  if (!requireNamespace("INLA", quietly = TRUE)) {
    stop("The INLA package is required by inlaST.estimate().")
  }
  y <- as.numeric(y)
  if (!length(y) || any(!is.finite(y))) {
    stop("y must be a non-empty finite numeric vector.")
  }
  if (is.null(offset)) offset <- numeric(length(y))
  offset <- as.numeric(offset)
  if (length(offset) != length(y) || any(!is.finite(offset))) {
    stop("offset must be a finite vector with length(y) entries.")
  }
  ctl <- .inlast_control(control)
  z <- .inlast_validate_spec(spec, y, offset)

  fixed_data <- data.frame(.inlast_offset = offset, check.names = FALSE)
  fixed_internal <- if (ncol(z$X)) {
    paste0(".inlast_x", seq_len(ncol(z$X)))
  } else character()
  for (j in seq_len(ncol(z$X))) fixed_data[[fixed_internal[j]]] <- z$X[, j]
  random_internal <- paste0(".inlast_r", seq_along(z$random))
  effects <- list(fixed_data)
  Astack <- list(1)
  for (j in seq_along(z$random)) {
    effects[[j + 1L]] <- seq_len(ncol(z$random[[j]]$A))
    names(effects)[j + 1L] <- random_internal[j]
    Astack[[j + 1L]] <- z$random[[j]]$A
  }
  stack <- INLA::inla.stack(
    data = list(y = y), A = Astack, effects = effects, tag = "estimation"
  )

  fenv <- new.env(parent = asNamespace("INLA"))
  rhs <- if (length(fixed_internal)) {
    paste(c("-1", fixed_internal, "offset(.inlast_offset)"), collapse = " + ")
  } else {
    "-1 + offset(.inlast_offset)"
  }
  fixed_precision <- ctl$fixed_precision
  if (!is.null(fixed_precision)) {
    fixed_precision <- as.numeric(fixed_precision)
    if (length(fixed_precision) == 1L) {
      fixed_precision <- rep(fixed_precision, length(z$random))
    }
    if (length(fixed_precision) != length(z$random) ||
        any(!is.finite(fixed_precision)) || any(fixed_precision <= 0)) {
      stop("control$fixed_precision must be positive and have one value per random block.")
    }
  }
  for (j in seq_along(z$random)) {
    qname <- paste0(".inlast_Q", j)
    hname <- paste0(".inlast_H", j)
    assign(qname, z$random[[j]]$precision_scale * z$random[[j]]$Q,
           envir = fenv)
    fixed_j <- if (is.null(fixed_precision)) NULL else
      fixed_precision[j] / z$random[[j]]$precision_scale
    assign(hname, list(prec = .inlast_prior_hyper(
      ctl$precision_prior, fixed = fixed_j
    )), envir = fenv)
    constraint_text <- ""
    if (!is.null(z$constraints[[j]])) {
      cname <- paste0(".inlast_C", j)
      assign(cname, list(
        A = matrix(z$constraints[[j]], nrow = 1L), e = 0
      ), envir = fenv)
      constraint_text <- paste0(", extraconstr=", cname)
    }
    # Supplying rankdef overrides INLA's automatic constraint adjustment.
    # generic0's tau normalizer must use the dimension of the constrained
    # support: rank(Q) minus the number of independent exact constraints.
    inla_rankdef <- z$random[[j]]$rankdef +
      as.integer(!is.null(z$constraints[[j]]))
    if (inla_rankdef >= ncol(z$random[[j]]$Q)) {
      stop("Random block '", z$random[[j]]$name,
           "' has no proper constrained support.")
    }
    rhs <- paste0(
      rhs, " + f(", random_internal[j],
      ", model='generic0', Cmatrix=", qname,
      ", rankdef=", inla_rankdef,
      ", constr=FALSE, hyper=", hname,
      constraint_text, ")"
    )
  }
  formula <- stats::as.formula(paste("y ~", rhs), env = fenv)

  control_family <- if (z$family == "gaussian") {
    list(hyper = list(prec = .inlast_prior_hyper(
      ctl$gaussian_precision_prior, fixed = ctl$gaussian_precision
    )))
  } else if (z$family == "negative_binomial") {
    list(hyper = list(theta = .inlast_prior_hyper(
      ctl$nb_size_prior, fixed = ctl$nb_size
    )))
  } else {
    list()
  }
  inla_family <- switch(
    z$family, gaussian = "gaussian", poisson = "poisson",
    negative_binomial = "nbinomial"
  )
  t0 <- proc.time()[["elapsed"]]
  native_control <- .inlast_merge_control(
    ctl$control.inla,
    list(strategy = ctl$latent_strategy, int.strategy = ctl$int_strategy)
  )
  fit <- INLA::inla(
    formula, family = inla_family, data = INLA::inla.stack.data(stack),
    control.predictor = list(A = INLA::inla.stack.A(stack), compute = TRUE),
    control.fixed = list(
      mean = 0, prec = ctl$fixed_effect_precision,
      mean.intercept = 0, prec.intercept = ctl$fixed_effect_precision
    ),
    control.family = control_family,
    control.inla = native_control,
    control.compute = list(config = TRUE),
    num.threads = ctl$num_threads, verbose = ctl$verbose
  )
  fit_seconds <- proc.time()[["elapsed"]] - t0

  beta <- if (ncol(z$X)) vapply(
    fixed_internal,
    function(tag) .inlast_latent_mode_block(fit, tag, 1L),
    numeric(1L)
  ) else numeric()
  names(beta) <- z$xnames
  random_mode <- vector("list", length(z$random))
  names(random_mode) <- vapply(z$random, `[[`, character(1L), "name")
  for (j in seq_along(z$random)) {
    random_mode[[j]] <- .inlast_latent_mode_block(
      fit, random_internal[j], ncol(z$random[[j]]$A)
    )
  }
  constraint_residual_uncorrected <- vapply(seq_along(z$random), function(j) {
    if (is.null(z$constraints[[j]])) return(NA_real_)
    sum(z$constraints[[j]] * random_mode[[j]])
  }, numeric(1L))
  names(constraint_residual_uncorrected) <- names(random_mode)
  active_constraint <- !is.na(constraint_residual_uncorrected)
  if (any(active_constraint)) {
    scale <- vapply(which(active_constraint), function(j) {
      1 + sum(abs(z$constraints[[j]]) * abs(random_mode[[j]]))
    }, numeric(1L))
    if (any(abs(constraint_residual_uncorrected[active_constraint]) >
            1e-6 * scale)) {
      stop("INLA returned a latent mode that violates an active spatial constraint.")
    }
    # INLA's constraint solve is accurate to its numerical tolerance.  Remove
    # that final roundoff component so all returned fields obey g'u = 0 to
    # machine precision, matching the projected mgcvST representation.
    for (j in which(active_constraint)) {
      g <- z$constraints[[j]]
      random_mode[[j]] <- random_mode[[j]] -
        g * sum(g * random_mode[[j]]) / sum(g * g)
    }
  }
  eta <- offset
  if (length(beta)) eta <- eta + as.numeric(z$X %*% beta)
  for (j in seq_along(z$random)) {
    eta <- eta + as.numeric(z$random[[j]]$A %*% random_mode[[j]])
  }

  tau_internal <- numeric(length(z$random))
  for (j in seq_along(z$random)) {
    tau_internal[j] <- if (is.null(fixed_precision)) {
      .inlast_hyper_mode(fit, paste0(
        "^log precision for ", random_internal[j], "$"
      ))
    } else fixed_precision[j] / z$random[[j]]$precision_scale
  }
  if (any(!is.finite(tau_internal)) || any(tau_internal <= 0)) {
    stop("INLA did not return positive latent precision modes.")
  }
  precision_scale <- vapply(z$random, `[[`, numeric(1L), "precision_scale")
  tau <- precision_scale * tau_internal
  names(tau_internal) <- names(random_mode)
  names(precision_scale) <- names(random_mode)
  names(tau) <- names(random_mode)

  if (z$family == "gaussian") {
    observation_precision <- if (is.null(ctl$gaussian_precision)) {
      .inlast_hyper_mode(
        fit, "^log precision for the gaussian observations$"
      )
    } else ctl$gaussian_precision
    if (!is.finite(observation_precision) || observation_precision <= 0) {
      stop("INLA did not return a positive Gaussian observation precision mode.")
    }
    dispersion <- 1 / observation_precision
    family_parameters <- numeric()
    mu <- eta
    working_response <- y
    working_variance <- rep.int(dispersion, length(y))
  } else if (z$family == "negative_binomial") {
    nb_size <- if (is.null(ctl$nb_size)) {
      .inlast_hyper_mode(
        fit, "^log size for the nbinomial observations"
      )
    } else ctl$nb_size
    if (!is.finite(nb_size) || nb_size <= 0) {
      stop("INLA did not return a positive negative-binomial size mode.")
    }
    dispersion <- 1
    family_parameters <- nb_size
    mu <- exp(eta)
    if (any(!is.finite(mu)) || any(mu <= 0)) {
      stop("The INLA negative-binomial conditional mode produced invalid means.")
    }
    working_response <- eta + (y - mu) / mu
    working_variance <- 1 / mu + 1 / nb_size
  } else {
    dispersion <- 1
    family_parameters <- numeric()
    mu <- exp(eta)
    if (any(!is.finite(mu)) || any(mu <= 0)) {
      stop("The INLA Poisson conditional mode produced invalid means.")
    }
    working_response <- eta + (y - mu) / mu
    working_variance <- 1 / mu
  }
  working_error <- working_response - offset
  if (any(!is.finite(working_error)) ||
      any(!is.finite(working_variance)) || any(working_variance <= 0)) {
    stop("The final INLA conditional mode produced an invalid IRLS system.")
  }

  lambda <- dispersion * tau
  # Exact indexing is required here because a public spec also contains
  # `sp_names`; `$sp` would partially match that character vector.
  smoothing_parameters <- spec[["sp", exact = TRUE]]
  if (is.null(smoothing_parameters)) {
    smoothing_parameters <- rep(NA_real_, max(c(0L, unlist(lapply(
      z$random, function(x) if (is.null(x$sp_index)) integer() else x$sp_index
    )))))
  } else {
    smoothing_parameters <- as.numeric(smoothing_parameters)
  }
  for (j in seq_along(z$random)) {
    index <- z$random[[j]]$sp_index
    if (is.null(index) || !length(index)) next
    if (length(index) != 1L || is.na(index) || index < 1L ||
        index > length(smoothing_parameters)) {
      stop("Each random block can map to at most one valid smoothing-parameter index.")
    }
    smoothing_parameters[index] <- lambda[j]
  }

  target <- which(vapply(z$random, `[[`, logical(1L), "target"))
  coefficients <- lapply(target, function(j) {
    projection <- z$random[[j]]$projection
    if (is.null(projection)) return(random_mode[[j]])
    projection <- as.matrix(projection)
    if (nrow(projection) != length(random_mode[[j]])) {
      stop("Target projection has incompatible dimensions.")
    }
    as.numeric(Matrix::crossprod(projection, random_mode[[j]]))
  })
  names(coefficients) <- names(random_mode)[target]

  # This second covariance solve is for diagnostics, never for the score.
  expected_nuisance_covariance <- if (diagnostics) {
    nuisance_index <- spec$nuisance_index
    if (is.null(nuisance_index)) nuisance_index <- seq_len(ncol(z$X))
    .inlast_expected_covariance(
      z$X, z$random, tau, working_variance, z$constraints, nuisance_index,
      fixed_effect_precision = ctl$fixed_effect_precision
    )
  } else NULL
  posterior <- .inlast_posterior_vp(fit, spec)
  nuisance_covariance <- posterior$nuisance_covariance
  if (!is.null(spec$nuisance_design)) {
    nuisance_design <- as.matrix(spec$nuisance_design)
    if (nrow(nuisance_design) != length(y) ||
        ncol(nuisance_design) != nrow(nuisance_covariance)) {
      stop("spec$nuisance_design is incompatible with spec$nuisance_index.")
    }
  }
  constraint_residual <- vapply(seq_along(z$random), function(j) {
    if (is.null(z$constraints[[j]])) return(NA_real_)
    sum(z$constraints[[j]] * random_mode[[j]])
  }, numeric(1L))
  names(constraint_residual) <- names(random_mode)
  observation_spatial_mean <- vapply(seq_along(z$random), function(j) {
    if (!identical(z$random[[j]]$kind, "spde")) return(NA_real_)
    mean(as.numeric(z$random[[j]]$A %*% random_mode[[j]]))
  }, numeric(1L))
  names(observation_spatial_mean) <- names(random_mode)
  mode_status <- fit$mode$mode.status
  mode_status_code <- if (is.null(mode_status)) NA_integer_ else
    suppressWarnings(as.integer(mode_status[1L]))
  mode_ok <- !is.null(mode_status) &&
    !is.na(mode_status_code) && identical(mode_status_code, 0L)
  theta_internal <- as.numeric(fit$mode$theta)
  converged <- mode_ok &&
    all(is.finite(theta_internal)) &&
    !any(!is.finite(c(beta, unlist(random_mode, use.names = FALSE))))
  log_marginal_likelihood <- if (!is.null(fit$mlik) && length(fit$mlik)) {
    as.numeric(fit$mlik[1L, 1L])
  } else NA_real_
  prior_metadata <- .inlast_prior_metadata(
    ctl, z$family, fixed_precision = fixed_precision
  )
  prior_semantics <- if (isTRUE(prior_metadata$any_improper)) {
    paste(
      "empirical Bayes conditional mode; at least one estimated",
      "log-hyperparameter uses an improper flat prior"
    )
  } else {
    paste(
      "empirical Bayes conditional mode with the explicitly recorded INLA",
      "hyperpriors for estimated positive parameters"
    )
  }
  hyper_mode_diagnostics <- list(
    optimizer_status = mode_status_code,
    optimizer_status_text = if (is.null(mode_status)) NA_character_ else
      paste(mode_status, collapse = "; "),
    theta_internal = theta_internal,
    theta_tags = fit$mode$theta.tags,
    finite = all(is.finite(theta_internal)),
    optimizer_trace_available = !is.null(fit$misc$opt.trace),
    covariance_eigenvalues = as.numeric(fit$misc$cov.intern.eigenvalues),
    warnings = fit$misc$warnings,
    boundary_check = if (isTRUE(prior_metadata$any_improper)) paste(
      "INLA reports optimizer status but no explicit unbounded-mode flag;",
      "flat-prior fits require multi-start sensitivity checks"
    ) else paste(
      "INLA reports optimizer status but no explicit unbounded-mode flag;",
      "boundary-sensitive fits require profile or multi-start checks"
    )
  )
  ans <- list(
    working_error = working_error,
    working_variance = working_variance,
    dispersion = dispersion,
    family = z$family,
    family_parameters = as.numeric(family_parameters),
    tau = tau,
    tau_internal = tau_internal,
    precision_scale = precision_scale,
    lambda = lambda,
    smoothing_parameters = smoothing_parameters,
    coefficients = coefficients,
    nuisance_covariance = nuisance_covariance,
    expected_nuisance_covariance = expected_nuisance_covariance,
    fixed_mode = beta,
    fixed_mean = beta,
    random_mode = random_mode,
    random_mean = random_mode,
    eta = eta,
    mu = mu,
    constraint_residual_uncorrected = constraint_residual_uncorrected,
    constraint_residual = constraint_residual,
    observation_spatial_mean = observation_spatial_mean,
    converged = converged,
    mode_status = mode_status_code,
    mode_status_text = if (is.null(mode_status)) NA_character_ else
      paste(mode_status, collapse = "; "),
    log_marginal_likelihood = log_marginal_likelihood,
    fit_seconds = fit_seconds,
    estimation = list(
      engine = "INLA generic0",
      hyperparameter_integration = ctl$int_strategy,
      latent_approximation = ctl$latent_strategy,
      prior_semantics = prior_semantics,
      hyperpriors = prior_metadata,
      hyper_mode_diagnostics = hyper_mode_diagnostics,
      fixed_kappa = TRUE,
      lambda_scaling = "lambda = dispersion * tau",
      nuisance_covariance = paste(
        "INLA conditional Gaussian posterior block at the empirical-Bayes",
        "configuration"
      ),
      posterior_covariance_diagnostics = posterior$diagnostics
    )
  )
  if (ctl$keep_fit) ans$inla <- fit
  ans
}
