# Extract only frozen-fit data, not a GAM/refit recipe or pair-score state.
.mgcvst_marginal_geometry <- function(fit, X = NULL, test_component = 1L) {
  if (is.null(X)) {
    X <- as.matrix(stats::predict(fit, newdata = fit$model, type = "lpmatrix"))
  }
  smooth <- lapply(fit$smooth, function(s) {
    ans <- s[c("first.para", "last.para", "first.sp", "last.sp", "S", "rank",
               "null.space.dim", "fixed", "label")]
    ans$getA <- if (is.null(s$getA)) NULL else TRUE
    ans
  })
  offset <- fit$offset
  if (is.null(offset)) offset <- rep(0, nrow(X))
  list(X = X, smooth = smooth, offset = as.numeric(offset),
       row_id = .mgcvst_row_id(fit, nrow(X)), test_component = test_component)
}

.mgcvst_capture_marginal <- function(fit, geometry = NULL, test_component = 1L) {
  .working_family_id(fit$family$family)
  if (is.null(geometry)) {
    geometry <- .mgcvst_marginal_geometry(fit, test_component = test_component)
  } else {
    current <- .mgcvst_marginal_geometry(fit, X = geometry$X,
                                        test_component = test_component)
    if (!identical(current, geometry) || nrow(geometry$X) != length(fit$y)) {
      stop("Marginal retention requires unchanged smooth geometry, rows and offsets within a chunk.")
    }
  }
  state <- fit[c("linear.predictors", "y", "prior.weights", "sig2",
                 "coefficients", "sp")]
  state$family_raw <- serialize(fit$family, NULL)
  list(state = state, geometry = geometry)
}

.mgcvst_collect_marginal <- function(chunks, p, feature_id) {
  state <- vector("list", p)
  geometry <- list()
  geometry_index <- integer(p)
  for (chunk in chunks) {
    g <- chunk$marginal_geometry
    if (!is.null(g)) {
      equal <- which(vapply(geometry, identical, logical(1L), y = g))
      if (!length(equal)) {
        geometry[[length(geometry) + 1L]] <- g
        equal <- length(geometry)
      }
      geometry_index[chunk$index] <- equal[1L]
    }
    state[chunk$index] <- chunk$marginal_state
  }
  names(state) <- feature_id
  list(version = 1L, geometry = geometry, geometry_index = geometry_index,
       state = state, definition = "frozen_fit_conditional_marginal_TAPS")
}

# Port of the standard/extended-family branches of extract_pseudo_response()
# and gammfast_working() in mgcv.taps; not the pairwise Fisher-IRLS extractor.
.mgcvst_marginal_working <- function(fit) {
  family <- fit$family
  extended <- inherits(family, "extended.family") && is.function(family$Dd)
  eta <- if (extended) as.numeric(fit$linear.predictors) else fit$linear.predictors
  y <- fit$y
  prior_weights <- if (extended && is.null(fit$prior.weights)) rep(1, length(y)) else
    as.numeric(fit$prior.weights)
  if (extended) {
    mu <- family$linkinv(eta)
    theta <- family$getTheta()
    Dval <- family$Dd(y, mu, theta, prior_weights, level = 0)
    mu_eta <- family$mu.eta(eta)
    D_eta <- Dval$Dmu * mu_eta
    w <- Dval$EDmu2 * mu_eta^2 / 2
    z <- eta - D_eta / (2 * w)
  } else {
    mu <- family$linkinv(eta)
    mu_eta <- family$mu.eta(eta)
    variance <- family$variance(mu)
    z <- eta + (y - mu) / mu_eta
    w <- prior_weights * mu_eta^2 / variance
  }
  if (length(z) != length(y) || length(w) != length(y) ||
      any(!is.finite(z)) || any(!is.finite(w)) || any(w < 0) || !any(w > 0)) {
    stop("The family produced an invalid diagonal PIRLS working system.")
  }
  z <- as.numeric(z)
  w <- as.numeric(w)
  phi0 <- fit$sig2
  if (is.null(phi0) || length(phi0) != 1L || !is.numeric(phi0) ||
      !is.finite(phi0) || phi0 <= 0) phi0 <- 1.0
  phi0 <- as.numeric(phi0)
  list(pseudo_response = z, V_phi = phi0 / w, phi0 = phi0,
       valid_idx = if (extended) is.finite(z) & is.finite(w) & w > 1e-12 else NULL)
}

# Preserve the upstream square-root and marginal Liu arithmetic, including
# the original eigensolver, symmetrization order and spectral cutoff.
.mgcvst_marginal_matrixsqrt <- function(A) {
  fit <- CppMatrix::matrixEigen(t(A) / 2 + A / 2)
  d <- c(fit$value)
  d[d <= 0] <- 0
  d <- sqrt(d)
  list(w = CppMatrix::matrixMultiply(fit$vector, t(fit$vector) * d))
}

.mgcvst_marginal_liu <- function(q, moments) {
  c1 <- moments[1L]; c2 <- moments[2L]
  c3 <- moments[3L]; c4 <- moments[4L]
  s1 <- c3 / (c2^(3/2))
  s2 <- c4 / c2^2
  muQ <- c1
  sigmaQ <- sqrt(2 * c2)
  tstar <- (q - muQ) / sigmaQ
  if (s1^2 > s2) {
    a <- 1 / (s1 - sqrt(s1^2 - s2))
    delta <- s1 * a^3 - a^2
    l <- a^2 - 2 * delta
  } else {
    delta <- 0
    l <- 1 / s2
    a <- sqrt(l)
  }
  muX <- l + delta
  sigmaX <- sqrt(2) * a
  stats::pchisq(tstar * sigmaX + muX, df = l, ncp = delta, lower.tail = FALSE)
}

.mgcvst_marginal_moments <- function(lambda) {
  c(sum(lambda), sum(lambda^2), sum(lambda^3), sum(lambda^4))
}

.mgcvst_marginal_davies <- function(z, fallback, max_eps, max_iter) {
  d <- tryCatch(CompQuadForm::davies(q = z$statistic, lambda = z$lambda,
                                    lim = max_iter, acc = max_eps),
                error = function(e) e)
  failed <- inherits(d, "condition") || !is.finite(d$Qq) ||
    d$Qq <= 0 || d$Qq > 1 || d$ifault != 0L
  reason <- if (!failed) NA_character_ else if (inherits(d, "condition")) {
    conditionMessage(d)
  } else {
    paste0("Davies numerical failure: Qq=", d$Qq, ", ifault=", d$ifault)
  }
  ifault <- if (inherits(d, "condition")) NA_integer_ else as.integer(d$ifault)
  if (failed && fallback == "liu") {
    p <- .mgcvst_marginal_liu(z$statistic, .mgcvst_marginal_moments(z$lambda))
    list(p_value = p, method_used = "liu", fallback_used = TRUE,
         fallback_reason = reason, davies_ifault = ifault, error_message = NA_character_)
  } else {
    list(p_value = if (failed) NA_real_ else as.numeric(d$Qq),
         method_used = "davies", fallback_used = FALSE,
         fallback_reason = NA_character_, davies_ifault = ifault,
         error_message = reason)
  }
}

.mgcvst_marginal_chunk <- function(payload, geometry, calibration, fallback,
                                  null.tol, max_eps, max_iter, n_threads = 1L) {
  .mgcvst_thread_limit()
  if (!identical(n_threads, 1L)) stop("Marginal workers require n_threads = 1.")
  lapply(seq_along(payload$index), function(k) {
    tryCatch({
      state <- payload$state[[k]]
      if (inherits(state, "condition")) stop(state)
      if (is.null(state)) stop("No retained marginal state for this feature.")
      z <- .mgcvst_marginal_spectrum(state, geometry[[payload$geometry_index[k]]], null.tol)
      if (!length(z$lambda) || !is.finite(z$statistic) || any(!is.finite(z$lambda))) {
        stop("Marginal TAPS produced an empty or non-finite quadratic-form spectrum.")
      }
      if (calibration == "davies") {
        z$calibration <- .mgcvst_marginal_davies(z, fallback, max_eps, max_iter)
        z$lambda <- NULL
      }
      z
    }, error = function(e) e)
  })
}

#' Marginal TAPS evaluation of already estimated features
#'
#' Evaluates frozen-fit conditional marginal TAPS, ported from mgcv.taps,
#' without refitting. This is not the pairwise covariance test. Estimate with
#' `retain_marginal = TRUE` to retain the needed response, coefficients and
#' family state. Pairwise working quantities alone are not sufficient.
#' Shared design and penalties are stored once; feature spectra are computed
#' through BiocParallel with one numerical thread per worker. For Liu, a
#' separate marginal-only OpenMP kernel sums spectral powers across features;
#' R's scalar noncentral chi-square tail calculation is retained. No nested
#' OpenMP runs inside BiocParallel workers. Davies failures remain NA unless
#' `fallback = "liu"` is explicitly requested; the actual method is reported.
#'
#' @param fitmgcvST Result of [mgcvST.estimate()] with `retain_marginal = TRUE`.
#' @param features Unique feature IDs or one-based indices; NULL selects all.
#' @param calibration Either `"davies"` or `"liu"`.
#' @param fallback Either `"none"` (default) or explicit Davies-to-Liu fallback.
#' @param BPPARAM BiocParallel configuration; defaults to serial execution.
#' @param chunk_size Number of features per interruptible task/block.
#' @param threads OpenMP threads for the parent-process marginal Liu kernel.
#' @param null.tol Upstream penalty-column null-space threshold.
#' @param max_eps,max_iter Davies accuracy and integration limit.
#' @return A data.frame in requested feature order with statistic, p-value,
#' requested/actual method, fallback and numerical failure information. No
#' filtering or multiple-testing correction is applied implicitly.
#' @export
mgcvST.marginal <- function(fitmgcvST, features = NULL,
    calibration = c("davies", "liu"), fallback = c("none", "liu"),
    BPPARAM = BiocParallel::SerialParam(), chunk_size = 100L, threads = 1L,
    null.tol = 1e-10, max_eps = 1e-8, max_iter = 1e5) {
  calibration <- match.arg(calibration)
  fallback <- match.arg(fallback)
  if (!inherits(fitmgcvST, "mgcvST_fit") || is.null(fitmgcvST$marginal_data)) {
    stop("Estimate with retain_marginal = TRUE before calling mgcvST.marginal().")
  }
  if (!inherits(BPPARAM, "BiocParallelParam")) stop("BPPARAM must inherit from 'BiocParallelParam'.")
  for (name in c("chunk_size", "threads", "max_iter")) {
    value <- get(name)
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
        value < 1 || value > .Machine$integer.max || value != as.integer(value)) {
      stop(name, " must be a positive integer.")
    }
  }
  for (name in c("null.tol", "max_eps")) {
    value <- get(name)
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) || value <= 0) {
      stop(name, " must be positive and finite.")
    }
  }
  if (calibration == "davies" && !requireNamespace("CompQuadForm", quietly = TRUE)) {
    stop("calibration = 'davies' requires CompQuadForm.")
  }
  ids <- fitmgcvST$feature_id
  if (is.null(features)) features <- seq_along(ids)
  index <- if (is.character(features)) match(features, ids) else features
  if (!is.numeric(index) || !length(index) || anyNA(index) ||
      any(!is.finite(index)) || any(index < 1 | index > length(ids)) ||
      any(index != as.integer(index)) || anyDuplicated(index)) {
    stop("features must contain unique known feature IDs or valid indices.")
  }
  index <- as.integer(index)
  data <- fitmgcvST$marginal_data
  .mgcvst_thread_limit()
  chunks <- split(seq_along(index), ceiling(seq_along(index) / chunk_size))
  payload <- lapply(chunks, function(rows) {
    i <- index[rows]
    list(index = i, state = data$state[i], geometry_index = data$geometry_index[i])
  })
  evaluated <- BiocParallel::bplapply(payload, .mgcvst_marginal_chunk,
    geometry = data$geometry, calibration = calibration, fallback = fallback,
    null.tol = null.tol, max_eps = max_eps, max_iter = max_iter,
    n_threads = 1L, BPPARAM = BPPARAM)
  evaluated <- unlist(evaluated, recursive = FALSE)
  result <- data.frame(feature_id = ids[index], statistic = NA_real_, p_value = NA_real_,
    method_requested = calibration, method_used = calibration, fallback_used = FALSE,
    fallback_reason = NA_character_, davies_ifault = NA_integer_, error_message = NA_character_,
    stringsAsFactors = FALSE)
  good <- which(!vapply(evaluated, inherits, logical(1L), what = "condition"))
  for (k in setdiff(seq_along(index), good)) {
    result$error_message[k] <- conditionMessage(evaluated[[k]])
  }
  if (calibration == "liu" && length(good)) {
    blocks <- split(good, ceiling(seq_along(good) / chunk_size))
    for (block in blocks) {
      # Powers deliberately use R arithmetic; the OpenMP kernel only sums.
      powers <- lapply(evaluated[block], function(z) {
        x <- z$lambda
        cbind(x, x^2, x^3, x^4)
      })
      moments <- mgcvst_marginal_liu_moments_cpp(powers, as.integer(threads))
      for (j in seq_along(block)) {
        k <- block[j]
        p <- tryCatch(.mgcvst_marginal_liu(evaluated[[k]]$statistic, moments[j, ]),
                      error = function(e) e)
        if (inherits(p, "condition")) result$error_message[k] <- conditionMessage(p) else
          result$p_value[k] <- p
      }
    }
  }
  for (k in good) {
    result$statistic[k] <- evaluated[[k]]$statistic
    if (calibration == "davies") {
      for (name in names(evaluated[[k]]$calibration)) {
        result[k, name] <- evaluated[[k]]$calibration[[name]]
      }
    }
    if (!is.finite(result$p_value[k]) || result$p_value[k] < 0 || result$p_value[k] > 1) {
      result$p_value[k] <- NA_real_
      if (is.na(result$error_message[k])) result$error_message[k] <- "Invalid marginal p-value."
    }
  }
  result
}
