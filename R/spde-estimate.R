#' Spatial fits of count features with one SPDE field
#'
#' `spde.estimate()` fits, feature by feature, the log-link model
#' `y ~ NB(mu, size)`, `log(mu) = offset + X beta + A u`, with the field
#' `u ~ N(0, (tau Q)^-1)` on the constraint `sum(constraint * u) = 0` and a
#' flat `beta`. The precision `tau` and the overdispersion `alpha = 1 / size`
#' maximise the Laplace approximation of the marginal likelihood (the REML
#' criterion of [mgcv::gam()]); the fit is compiled, with a sparse Cholesky
#' factor of `A'WA + tau Q` whose symbolic analysis is done once.
#'
#' The overdispersion is searched on `alpha >= 0`, where the Poisson model is
#' the finite end point `alpha = 0`. The Poisson fit is found first. The
#' negative-binomial maximum is then searched jointly in
#' `(log(tau), log(alpha))`: line searches over each in turn bring the point
#' near the maximum, and Newton steps in both, with the gradient and the
#' Hessian of the criterion by central differences, finish it. The search
#' starts at the best `alpha` at the Poisson `tau` and, when that gains
#' nothing, at the `alpha` of the null model (`size_start`). The
#' negative-binomial fit is returned when its criterion exceeds the Poisson
#' one; a field that absorbs all the overdispersion therefore gives the
#' Poisson fit, not a size that runs to infinity. At every `(tau, alpha)` the
#' fixed effects and the field are solved again to their joint mode.
#'
#' Every stopping rule is a relative change of the objective,
#' `|f - f_old| / (|f| + 0.1)`, as [stats::glm.fit()] measures its deviance;
#' for a Newton iteration the change is the predicted gain of its step. The
#' gradient of the criterion at the returned point is reported, so that each
#' fit can be checked.
#'
#' @param Y Feature-by-observation count matrix (dense or sparse).
#' @param X Observation-by-covariate design of the fixed effects.
#' @param A Sparse observation-by-vertex projector of the field.
#' @param Q Sparse precision of the field at `tau = 1`.
#' @param offset `NULL`, a length-`n` vector, or a feature-by-observation
#'   matrix.
#' @param poisson Optional logical vector: `TRUE` fits the Poisson model only.
#' @param size_start Optional sizes of the null model, one per feature: the
#'   first point of the line search of `alpha = 1 / size`, and the start of
#'   the second search when the first ends on the Poisson model.
#' @param constraint Vector `c` of the constraint `c'u = 0`; the default is
#'   the mean of the field over the observations, `crossprod(A, 1 / n)`.
#' @param feature_id Feature identifiers.
#' @param threads OpenMP threads (one feature per thread).
#' @param chunk_size Features per kernel call; bounds the memory of a call.
#' @param tol Relative change of the criterion at which the line searches
#'   and the joint Newton iteration over `(log(tau), log(alpha))` stop.
#' @param tol_inner Relative predicted gain of the penalised log-likelihood at
#'   which the Newton iteration for the mode stops. The criterion has a
#'   log-determinant that is not stationary at the mode, so its error is of
#'   the order of the square root of this gain: `1e-13` keeps it below `tol`
#'   (on 80 simulated features of 1500 observations at most 1.5e-6, against
#'   5e-4 with `1e-10`).
#' @param maxit_inner,maxit_halving,maxit_line,maxit_outer Limits of the
#'   Newton iteration for the mode, of a step halving, of the evaluations of a
#'   line search, and of the line-search cycles and of the joint Newton steps.
#' @param fd_step Step of the central differences of the criterion in
#'   `log(tau)` and `log(alpha)`.
#' @param step_max Longest joint Newton step, in `log(tau)` or `log(alpha)`.
#' @param log_tau_min,log_tau_max Search interval of `log(tau)`.
#' @param alpha_min,alpha_max Search interval of a positive `alpha`
#'   (`alpha = 0` is the Poisson fit).
#' @param eta_max Largest linear predictor a step may reach.
#' @param rank_tol Relative eigenvalue below which a direction of the fixed
#'   effects is aliased (a factor level whose fitted means have vanished).
#' @return A data frame with one row per feature: `feature_id`, `family`,
#'   `tau`, `size`, `laml` (the criterion, up to a constant), `laml_poisson`,
#'   `converged`, `grad_log_tau`, `grad_log_alpha` (the gradient of the
#'   criterion at the returned point; `NA` where a difference would leave the
#'   search interval, and for `alpha` of a Poisson fit), `evaluations`,
#'   `newton`, `outer`, `error`. A feature whose negative-binomial search
#'   fails everywhere has an error; it is not returned as a Poisson fit. The attributes
#'   `"coefficients"` and `"field"` hold `beta` and `u` by feature.
#' @keywords internal
spde.estimate <- function(Y, X, A, Q, offset = NULL, poisson = NULL,
                          size_start = NULL, constraint = NULL, feature_id = rownames(Y),
                          threads = 1L, chunk_size = 64L, tol = 1e-8,
                          tol_inner = 1e-13, maxit_inner = 50L,
                          maxit_halving = 30L, maxit_line = 100L,
                          maxit_outer = 25L, log_tau_min = -12,
                          log_tau_max = 12, alpha_min = 1e-6, alpha_max = 1e4,
                          eta_max = 700, rank_tol = 1e-11, fd_step = 0.01,
                          step_max = 2) {
  X <- as.matrix(X)
  storage.mode(X) <- "double"
  n <- nrow(X)
  p <- nrow(Y)
  if (ncol(Y) != n) stop("Y must have one column per row of X.")
  A <- methods::as(methods::as(methods::as(A, "CsparseMatrix"), "generalMatrix"),
                   "dMatrix")
  Q <- methods::as(methods::as(methods::as(Matrix::forceSymmetric(
    methods::as(Q, "CsparseMatrix")), "generalMatrix"), "CsparseMatrix"), "dMatrix")
  if (nrow(A) != n || ncol(Q) != ncol(A)) stop("A and Q must match X and each other.")
  if (is.null(constraint)) constraint <- as.numeric(Matrix::crossprod(A, rep.int(1 / n, n)))
  constraint <- as.numeric(constraint)
  if (length(constraint) != ncol(A) || any(!is.finite(constraint)) ||
      !any(constraint != 0)) {
    stop("constraint must be a finite non-zero vector with ncol(A) entries.")
  }
  if (is.null(feature_id)) feature_id <- as.character(seq_len(p))
  if (!is.null(poisson) && length(poisson) != p) stop("poisson needs one flag per feature.")
  if (!is.null(size_start) && length(size_start) != p) stop("size_start needs one value per feature.")
  if (!is.null(offset) && !is.matrix(offset) && length(offset) != n) {
    stop("offset must be NULL, a length-n vector or a feature-by-observation matrix.")
  }
  out <- data.frame(
    feature_id = feature_id, family = NA_character_, tau = NA_real_,
    size = NA_real_, laml = NA_real_, laml_poisson = NA_real_,
    converged = FALSE, grad_log_tau = NA_real_, grad_log_alpha = NA_real_,
    evaluations = NA_integer_, newton = NA_integer_,
    outer = NA_integer_, error = NA_character_, stringsAsFactors = FALSE
  )
  beta <- matrix(NA_real_, ncol(X), p, dimnames = list(colnames(X), feature_id))
  field <- matrix(NA_real_, ncol(A), p, dimnames = list(NULL, feature_id))
  for (rows in split(seq_len(p), ceiling(seq_len(p) / chunk_size))) {
    Yt <- methods::as(methods::as(methods::as(
      Matrix::t(Y[rows, , drop = FALSE]), "CsparseMatrix"), "generalMatrix"),
      "dMatrix")
    O <- if (is.null(offset)) matrix(0, n, 1L) else if (is.matrix(offset))
      t(offset[rows, , drop = FALSE]) else matrix(as.numeric(offset), n, 1L)
    storage.mode(O) <- "double"
    route <- if (is.null(poisson)) rep(FALSE, length(rows)) else as.logical(poisson[rows])
    fit <- mgcvst_spde_estimate_cpp(
      Yt, X, A, Q, constraint, O, route,
      if (is.null(size_start)) rep(NA_real_, length(rows)) else 1 / as.numeric(size_start[rows]),
      tol = tol, tol_inner = tol_inner,
      maxit_inner = as.integer(maxit_inner),
      maxit_halving = as.integer(maxit_halving),
      maxit_line = as.integer(maxit_line), maxit_outer = as.integer(maxit_outer),
      log_tau_min = log_tau_min, log_tau_max = log_tau_max,
      alpha_min = alpha_min, alpha_max = alpha_max, eta_max = eta_max,
      rank_tol = rank_tol, fd_step = fd_step, step_max = step_max,
      threads = as.integer(threads)
    )
    beta[, rows] <- fit$beta
    field[, rows] <- fit$u
    out$family[rows] <- ifelse(is.na(fit$poisson), NA_character_,
                               ifelse(fit$poisson, "poisson", "negative_binomial"))
    out$tau[rows] <- fit$tau
    out$size[rows] <- ifelse(fit$alpha > 0, 1 / fit$alpha, Inf)
    out$laml[rows] <- fit$laml
    out$laml_poisson[rows] <- fit$laml_poisson
    out$converged[rows] <- fit$converged
    out$grad_log_tau[rows] <- fit$grad_log_tau
    out$grad_log_alpha[rows] <- fit$grad_log_alpha
    out$evaluations[rows] <- fit$evaluations
    out$newton[rows] <- fit$newton
    out$outer[rows] <- fit$outer
    out$error[rows] <- fit$error
  }
  attr(out, "coefficients") <- beta
  attr(out, "field") <- field
  out
}

#' The criterion of `spde.estimate()` at given hyperparameters
#'
#' The Laplace marginal likelihood of one feature at the pairs
#' `(log_tau[k], alpha[k])`, each at its own mode of the fixed effects and the
#' field; `alpha = 0` is the Poisson model. It is the objective
#' [spde.estimate()] maximises, for checking a fit on a grid.
#'
#' @param y Counts of one feature.
#' @inheritParams spde.estimate
#' @param offset `NULL` or a length-`n` vector.
#' @param log_tau,alpha Points of evaluation, recycled to a common length.
#' @return The criterion at each point, up to the constant of
#'   [spde.estimate()]; `NA` where the mode is not found.
#' @keywords internal
spde.laml <- function(y, X, A, Q, log_tau, alpha, offset = NULL, constraint = NULL,
                      tol_inner = 1e-13, maxit_inner = 50L, maxit_halving = 30L,
                      eta_max = 700, rank_tol = 1e-11) {
  X <- as.matrix(X)
  storage.mode(X) <- "double"
  n <- nrow(X)
  A <- methods::as(methods::as(methods::as(A, "CsparseMatrix"), "generalMatrix"),
                   "dMatrix")
  Q <- methods::as(methods::as(methods::as(Matrix::forceSymmetric(
    methods::as(Q, "CsparseMatrix")), "generalMatrix"), "CsparseMatrix"), "dMatrix")
  if (is.null(constraint)) constraint <- as.numeric(Matrix::crossprod(A, rep.int(1 / n, n)))
  if (is.null(offset)) offset <- numeric(n)
  k <- max(length(log_tau), length(alpha))
  mgcvst_spde_laml_cpp(
    as.numeric(y), X, A, Q, as.numeric(constraint), as.numeric(offset),
    rep_len(as.numeric(log_tau), k), rep_len(as.numeric(alpha), k),
    tol_inner = tol_inner, maxit_inner = as.integer(maxit_inner),
    maxit_halving = as.integer(maxit_halving), eta_max = eta_max, rank_tol = rank_tol)
}

# The models of engine = "spde": a count family, fixed effects and one SPDE
# field, the precision and the size both estimated. Another model is fitted
# with control = list(engine = "inla").
.spde_estimate_check <- function(spec, control) {
  other <- " Set control = list(engine = 'inla') for this model."
  if (!(spec$family %in% c("poisson", "negative_binomial"))) {
    stop("engine = 'spde' fits Poisson and negative-binomial models only.", other)
  }
  if (length(spec$random) != 1L || !identical(spec$random[[1L]]$kind, "spde")) {
    stop("engine = 'spde' requires fixed effects and one SPDE field, without other random blocks.",
         other)
  }
  if (!is.null(control$nb_size) || !is.null(control$fixed_precision) ||
      !is.null(spec$nb_size_fixed)) {
    stop("engine = 'spde' estimates the precision and the size; neither can be fixed.", other)
  }
  invisible(TRUE)
}

# Stage 2 fits of one chunk of features by spde.estimate(), as the records of
# .inlast_fit_feature(): one list per feature, or a condition for a feature
# whose fit failed. The model is the count model with fixed effects and one
# SPDE field; tau multiplies the raw Q of the block.
.spde_estimate_chunk <- function(index, Y, spec, base_offset, extra_offset,
                                 control, poisson = NULL, threads = 1L,
                                 null = NULL) {
  .spde_estimate_check(spec, control)
  Yk <- Y[index, , drop = FALSE]
  n <- ncol(Yk)
  offset <- if (is.matrix(extra_offset)) {
    sweep(extra_offset[index, , drop = FALSE], 2L, base_offset, "+")
  } else if (is.null(extra_offset)) base_offset else base_offset + extra_offset
  z <- .inlast_validate_spec(spec, Yk[1L, ], numeric(n))
  block <- z$random[[1L]]
  route <- if (identical(z$family, "poisson")) rep(TRUE, length(index)) else
    if (is.null(poisson)) rep(FALSE, length(index)) else as.logical(poisson[index])
  # The size of the Stage 1 null model of each feature: the start of the
  # search, and kept beside the size of the spatial fit.
  null_size <- vapply(seq_along(index), function(k) {
    state <- if (is.null(null)) NULL else null[[index[k]]]
    size <- if (is.list(state)) state$family_parameters else NULL
    if (length(size) == 1L && is.finite(size) && size > 0) as.numeric(size) else NA_real_
  }, numeric(1L))
  t0 <- proc.time()[["elapsed"]]
  fit <- do.call(spde.estimate, c(list(
    Yk, z$X, block$A, block$Q, offset = offset, poisson = route,
    size_start = null_size,
    constraint = z$constraints[[1L]], feature_id = as.character(seq_along(index)),
    threads = threads, chunk_size = length(index)), control$spde))
  seconds <- (proc.time()[["elapsed"]] - t0) / length(index)
  beta <- attr(fit, "coefficients")
  field <- attr(fit, "field")
  sp <- spec[["sp", exact = TRUE]]
  sp <- if (is.null(sp)) rep(NA_real_, max(0L, block$sp_index)) else as.numeric(sp)
  lapply(seq_along(index), function(k) {
    if (!is.na(fit$error[k])) return(simpleError(fit$error[k]))
    u <- field[, k]
    b <- stats::setNames(beta[, k], z$xnames)
    b[is.na(b)] <- 0
    feature_offset <- if (is.matrix(offset)) offset[k, ] else offset
    eta <- feature_offset + as.numeric(block$A %*% u)
    if (length(b)) eta <- eta + as.numeric(z$X %*% b)
    mu <- exp(eta)
    nb <- identical(fit$family[k], "negative_binomial")
    tau <- stats::setNames(fit$tau[k], block$name)
    smoothing <- sp
    if (length(block$sp_index)) smoothing[block$sp_index] <- fit$tau[k]
    random_mode <- stats::setNames(list(u), block$name)
    projection <- block$projection
    coefficients <- stats::setNames(list(if (is.null(projection)) u else
      as.numeric(Matrix::crossprod(as.matrix(projection), u))), block$name)
    residual <- stats::setNames(sum(z$constraints[[1L]] * u), block$name)
    list(
      working_error = eta + (Yk[k, ] - mu) / mu - feature_offset,
      working_variance = 1 / mu + if (nb) 1 / fit$size[k] else 0,
      dispersion = 1, family = fit$family[k],
      family_parameters = if (nb) fit$size[k] else numeric(),
      tau = tau, tau_internal = tau / block$precision_scale,
      precision_scale = stats::setNames(block$precision_scale, block$name),
      lambda = tau, smoothing_parameters = smoothing,
      coefficients = coefficients, fixed_mode = b, fixed_mean = b,
      random_mode = random_mode, random_mean = random_mode, eta = eta, mu = mu,
      constraint_residual_uncorrected = residual, constraint_residual = residual,
      observation_spatial_mean = stats::setNames(mean(eta - feature_offset -
        if (length(b)) as.numeric(z$X %*% b) else 0), block$name),
      converged = isTRUE(fit$converged[k]), mode_status = NA_integer_,
      mode_status_text = NA_character_,
      log_marginal_likelihood = fit$laml[k], fit_seconds = seconds,
      null_size = null_size[k],
      gradient = c(log_tau = fit$grad_log_tau[k], log_alpha = fit$grad_log_alpha[k]),
      estimation = list(
        engine = "spde.estimate",
        criterion = "Laplace marginal likelihood (REML) in (tau, 1 / size)",
        laml_poisson = fit$laml_poisson[k], null_size = null_size[k], evaluations = fit$evaluations[k],
        newton = fit$newton[k], outer = fit$outer[k], fixed_kappa = TRUE,
        lambda_scaling = "lambda = tau"
      )
    )
  })
}
