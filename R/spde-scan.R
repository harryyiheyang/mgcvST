#' Stage 1 score scan of an SPDE field against a GLM null
#'
#' `spde.scan()` tests, for every feature, the variance component of an SPDE
#' field against a fixed-effects GLM null. It is Stage 1 of [spdeST.scan()],
#' which calls it with the raw field of the model, and it can be called on
#' its own.
#'
#' Each feature is fitted under the null `log(mu) = X beta + offset` by
#' maximum likelihood: NB2 with overdispersion `alpha = 1/size`, or Poisson
#' (`alpha = 0`) when the overdispersion score
#' `sum((y - mu)^2 - y) / 2` of the Poisson fit is not positive or the feature
#' is routed to Poisson. With working weights `W = mu / (1 + alpha mu)`,
#' working residual `(y - mu) / mu`,
#' `P0 = W - W X (X'WX)^{-1} X'W` and `B B' = Q^{-1}`,
#' the statistic is `T = |B'A'P0 r|^2` and its null law is the mixture
#' `sum(lambda_k chi^2_1)` over the positive eigenvalues of `B'A'P0AB`,
#' calibrated by the saddlepoint approximation. `P0` removes the span of
#' `X`, the constant included, so the field enters as it is: no intercept
#' projection, `Q` stays sparse. All numerical work is compiled
#' (`src/nb_null.cpp`, `src/spde_scan.cpp`, `src/spa_quadform.h`). The null
#' design `X` is a model matrix, so factors enter as 0-1 columns.
#'
#' @param Y Feature-by-observation count matrix (dense or sparse).
#' @param X Null design, observations by columns; it should hold the constant.
#' @param A Sparse observation-by-mesh projector of the field.
#' @param Q Sparse mesh precision of the field. Its scale does not change the
#'   p-values.
#' @param offset `NULL`, a length-n vector, or a feature-by-observation matrix.
#' @param poisson `NULL` or one logical per feature: `TRUE` fits Poisson.
#' @param size `NULL` (estimate the NB size) or a fixed positive size.
#' @param feature_id Feature identifiers.
#' @param threads OpenMP threads of the compiled kernels.
#' @param chunk_size Features per kernel call; bounds the memory of a call.
#' @param block_size Columns per block when the score curvature is built.
#' @param alpha_max Upper bound of the overdispersion (size >= 1 / alpha_max).
#' @param tol Relative change at which the IRLS (deviance) and the size search
#'   (profile log-likelihood) stop, `|l - l_old| / (|l| + 0.1) < tol`, as in
#'   [stats::glm.fit()] (whose default it is).
#' @param maxit_outer,maxit_inner Iteration limits of the size search and of
#'   the IRLS, both 25 as in [stats::glm.control()] and `MASS::glm.nb()` (a
#'   factor level with zero counts needs about 20 IRLS iterations; a very
#'   sparse feature with a large overdispersion needs about 12 size-search
#'   iterations, since the size and coefficient updates alternate and
#'   converge linearly). A null that does not converge gets no p-value.
#' @param alpha_min Floor of the moment start of the overdispersion.
#' @param step_bound,step_fallback A size-search step stays within a factor
#'   `step_bound` of the current overdispersion; where the curvature is not
#'   negative the step is a `step_fallback`-fold increase or decrease.
#' @param maxit_halving,eta_max An IRLS step is halved at most `maxit_halving`
#'   times to keep the linear predictor below `eta_max` and the deviance from
#'   rising.
#' @param rank_tol Relative pivot below which an IRLS step treats a direction
#'   of the design as aliased and does not move along it (the step is the
#'   minimum-norm least-squares solution on that rank); the default is that
#'   of [stats::glm.fit()] (`min(1e-7, epsilon / 1000)`).
#' @param spectrum_tol Eigenvalues below `spectrum_tol` times the largest are
#'   dropped from the mixture.
#' @param spa_centre_tol,spa_root_tol,spa_root_maxit Saddlepoint controls:
#'   within `spa_centre_tol` (relative) of the mixture mean the limit at the
#'   mean is used; the saddlepoint equation is solved to relative accuracy
#'   `spa_root_tol` in at most `spa_root_maxit` steps.
#' @return A data frame with one row per feature: `feature_id`, `family`
#'   (`"negative_binomial"` or `"poisson"`), `size`, `overdispersion_score`,
#'   `log_likelihood`, `converged`, `statistic`, `p_value`, `variance` (the
#'   moment estimate of the field variance `1 / tau` in the unit of `Q`:
#'   with the mixture spectrum `lambda`, `(T - sum(lambda)) / sum(lambda^2)`;
#'   it is not limited to positive values), `variance_se` (its standard error
#'   under the null, `sqrt(2 / sum(lambda^2))`), `error`. The attribute `"coefficients"` holds the
#'   `ncol(X)`-by-feature null coefficients.
#' @keywords internal
spde.scan <- function(Y, X, A, Q, offset = NULL, poisson = NULL, size = NULL,
                      feature_id = rownames(Y), threads = 1L,
                      chunk_size = 256L, block_size = 32L, alpha_max = 1e4,
                      tol = 1e-8,
                      maxit_outer = 25L, maxit_inner = 25L,
                      alpha_min = 1e-12, step_bound = 10, step_fallback = 2,
                      maxit_halving = 30L, eta_max = 700, rank_tol = 1e-11,
                      spectrum_tol = 1e-12, spa_centre_tol = 1e-6,
                      spa_root_tol = 1e-12, spa_root_maxit = 100L) {
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
  if (is.null(feature_id)) feature_id <- as.character(seq_len(p))
  if (!is.null(poisson) && length(poisson) != p) stop("poisson needs one flag per feature.")
  if (!is.null(size) && !(length(size) == 1L && is.finite(size) && size > 0)) {
    stop("size must be NULL or one positive number.")
  }
  if (!is.null(offset) && !is.matrix(offset) && length(offset) != n) {
    stop("offset must be NULL, a length-n vector or a feature-by-observation matrix.")
  }
  out <- data.frame(
    feature_id = feature_id, family = NA_character_, size = NA_real_,
    overdispersion_score = NA_real_, log_likelihood = NA_real_,
    converged = FALSE, statistic = NA_real_, p_value = NA_real_,
    variance = NA_real_, variance_se = NA_real_,
    error = NA_character_,
    stringsAsFactors = FALSE
  )
  beta <- matrix(NA_real_, ncol(X), p, dimnames = list(colnames(X), feature_id))
  for (rows in split(seq_len(p), ceiling(seq_len(p) / chunk_size))) {
    Yt <- methods::as(methods::as(methods::as(
      Matrix::t(Y[rows, , drop = FALSE]), "CsparseMatrix"), "generalMatrix"),
      "dMatrix")
    O <- if (is.null(offset)) matrix(0, n, 1L) else if (is.matrix(offset))
      t(offset[rows, , drop = FALSE]) else matrix(as.numeric(offset), n, 1L)
    storage.mode(O) <- "double"
    route <- if (is.null(poisson)) rep(FALSE, length(rows)) else as.logical(poisson[rows])
    fit <- mgcvst_nb_null_batch_cpp(
      Yt, X, O, route, alpha_max = alpha_max,
      alpha_fixed = if (is.null(size)) -1 else 1 / size, tol = tol,
      maxit_outer = as.integer(maxit_outer), maxit_inner = as.integer(maxit_inner),
      alpha_min = alpha_min, step_bound = step_bound, step_fallback = step_fallback,
      maxit_halving = as.integer(maxit_halving), eta_max = eta_max,
      rank_tol = rank_tol,
      threads = as.integer(threads)
    )
    ok <- is.na(fit$error) & fit$converged
    alpha <- ifelse(ok, fit$alpha, NA_real_)
    score <- mgcvst_spde_score_cpp(
      A, Q, X, Yt, O, fit$beta, alpha, threads = as.integer(threads),
      block_size = as.integer(block_size), spectrum_tol = spectrum_tol,
      spa_centre_tol = spa_centre_tol, spa_root_tol = spa_root_tol,
      spa_root_maxit = as.integer(spa_root_maxit)
    )
    beta[, rows] <- fit$beta
    out$family[rows] <- ifelse(fit$alpha > 0, "negative_binomial", "poisson")
    out$size[rows] <- fit$size
    out$overdispersion_score[rows] <- fit$score0
    out$log_likelihood[rows] <- fit$loglik
    out$converged[rows] <- ok
    out$error[rows] <- ifelse(!is.na(fit$error), fit$error,
      ifelse(!fit$converged, "The null size search did not converge.", score$error))
    out$statistic[rows] <- ifelse(ok, score$statistic, NA_real_)
    out$p_value[rows] <- ifelse(ok, score$p_value, NA_real_)
    out$variance[rows] <- ifelse(ok, score$variance, NA_real_)
    out$variance_se[rows] <- ifelse(ok, score$variance_se, NA_real_)
  }
  attr(out, "coefficients") <- beta
  out
}

