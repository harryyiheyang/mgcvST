# spde.scan: the Stage 1 null fits of the INLA branch without INLA.
#
# The Stage 1 null model of the INLA branch has no random block when the
# nuisance is fixed effects only: y_j ~ NB(mu, size) or Poisson with
# log(mu) = X beta_j + offset_j. Under the flat priors of the branch its INLA
# conditional mode is the joint maximum-likelihood estimate up to the Laplace
# correction of the size, so every feature is fitted by one batched
# maximum-likelihood kernel (mgcvst_nb_null_batch_cpp, src/nb_null.cpp: mgcv's
# nested scheme on the overdispersion alpha = 1/size, with Poisson as the end
# point alpha = 0 decided exactly by the overdispersion score). The working error
# and variance it returns are those of the INLA null fit, so the Stage 1 score
# statistic and its spectrum are computed by the same sparse kernel.

# A null spec qualifies when it has no random block and a count family;
# options(mgcvST.scan_null = FALSE) restores the per-feature INLA null fits.
.inlast_scan_eligible <- function(null_spec) {
  family <- null_spec$family
  isTRUE(getOption("mgcvST.scan_null", TRUE)) &&
    !length(null_spec$random) &&
    length(family) == 1L && family %in% c("negative_binomial", "poisson") &&
    !is.null(null_spec$fixed$X) && ncol(as.matrix(null_spec$fixed$X)) >= 1L
}

# Null fits of the features `rows` of Y (features x observations), in the
# record format of .inlast_fit_feature() for a spec without random blocks.
.inlast_scan_fits <- function(rows, Y, null_spec, base_offset, extra_offset,
                              control, poisson = NULL, threads = 1L) {
  t0 <- proc.time()[["elapsed"]]
  X <- as.matrix(null_spec$fixed$X)
  storage.mode(X) <- "double"
  n <- nrow(X)
  k <- length(rows)
  O <- matrix(as.numeric(base_offset), n, k)
  if (!is.null(extra_offset)) {
    O <- O + if (is.matrix(extra_offset)) t(extra_offset[rows, , drop = FALSE]) else
      as.numeric(extra_offset)
  }
  Yt <- methods::as(methods::as(methods::as(
    Matrix::t(Y[rows, , drop = FALSE]), "CsparseMatrix"), "generalMatrix"),
    "dMatrix")
  route <- if (identical(null_spec$family, "poisson")) rep(TRUE, k) else
    if (is.null(poisson)) rep(FALSE, k) else as.logical(poisson[rows])
  size <- control[["nb_size", exact = TRUE]]
  z <- mgcvst_nb_null_batch_cpp(
    Yt, X, O, route, alpha_fixed = if (is.null(size)) -1 else 1 / size,
    threads = as.integer(threads)
  )
  seconds <- (proc.time()[["elapsed"]] - t0) / max(1L, k)
  fixed_names <- colnames(null_spec$fixed$X)
  if (is.null(fixed_names)) fixed_names <- paste0("x", seq_len(ncol(X)))
  n_sp <- length(null_spec[["sp", exact = TRUE]])
  empty <- stats::setNames(numeric(), character())
  lapply(seq_len(k), function(i) {
    if (!is.na(z$error[i])) {
      return(simpleError(paste("spde.scan null fit:", z$error[i])))
    }
    beta <- stats::setNames(z$beta[, i], fixed_names)
    eta <- as.numeric(X %*% beta) + O[, i]
    mu <- exp(eta)
    y <- as.numeric(Y[rows[i], ])
    # alpha = 0 is the Poisson end point of the NB2 family.
    nb <- z$alpha[i] > 0
    family <- if (nb) "negative_binomial" else "poisson"
    theta <- z$size[i]
    list(
      working_error = eta + (y - mu) / mu - O[, i],
      working_variance = if (nb) 1 / mu + 1 / theta else 1 / mu,
      dispersion = 1,
      family = family,
      family_parameters = if (nb) theta else numeric(),
      tau = empty, tau_internal = empty, precision_scale = empty, lambda = empty,
      smoothing_parameters = rep(NA_real_, n_sp),
      coefficients = list(),
      fixed_mode = beta, fixed_mean = beta,
      random_mode = list(), random_mean = list(),
      eta = eta, mu = mu,
      constraint_residual_uncorrected = empty, constraint_residual = empty,
      observation_spatial_mean = empty,
      converged = isTRUE(z$converged[i]),
      mode_status = if (isTRUE(z$converged[i])) 0L else 1L,
      mode_status_text = if (isTRUE(z$converged[i])) "spde.scan converged" else
        "spde.scan overdispersion iteration limit",
      log_marginal_likelihood = NA_real_,
      log_likelihood = z$loglik[i],
      fit_seconds = seconds,
      estimation = list(
        engine = "spde.scan batched maximum likelihood",
        size_method = if (nb) "nested: IRLS for beta, Newton on alpha = 1/size" else
          "Poisson IRLS (overdispersion score at alpha = 0 not positive)",
        overdispersion_score = z$score0[i],
        outer_iterations = z$outer[i], inner_iterations = z$inner[i],
        prior_semantics = paste(
          "joint maximum likelihood; equals the flat-prior INLA conditional",
          "mode up to the Laplace correction of the size"
        )
      )
    )
  })
}

# Exact Poisson routing for the scan: a feature is Poisson when the
# overdispersion score of its Poisson fit, sum((y - mu)^2 - y) / 2, is not
# positive, which is the KKT condition of the NB2 maximum-likelihood end point
# alpha = 0. Returns the routing flags, the scores and the Pearson dispersions.
.inlast_scan_route <- function(Y, X, offset, threads = 1L) {
  X <- as.matrix(X)
  storage.mode(X) <- "double"
  n <- nrow(X)
  p <- nrow(Y)
  O <- if (is.null(offset)) matrix(0, n, 1L) else if (is.matrix(offset))
    t(offset) else matrix(as.numeric(offset), n, 1L)
  Yt <- methods::as(methods::as(methods::as(
    Matrix::t(Y), "CsparseMatrix"), "generalMatrix"), "dMatrix")
  z <- mgcvst_nb_null_batch_cpp(Yt, X, O, rep(TRUE, p),
                                threads = as.integer(threads))
  list(phi = z$phi, score = z$score0,
       poisson = is.finite(z$score0) & z$score0 <= 0)
}
