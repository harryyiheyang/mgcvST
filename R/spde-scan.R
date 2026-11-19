#' Stage 1 score scan of an SPDE field against a GLM null
#'
#' `spde.scan()` tests, for every feature, the variance component of an SPDE
#' field against a fixed-effects GLM null. It is Stage 1 of both
#' [mgcvST.estimate()] and [inlaST.estimate()], which call it with the raw
#' field of their basis, and it can be called on its own.
#'
#' Each feature is fitted under the null `log(mu) = X beta + offset` by
#' maximum likelihood: NB2 with overdispersion `alpha = 1/size`, or Poisson
#' (`alpha = 0`) when the overdispersion score
#' `sum((y - mu)^2 - y) / 2` of the Poisson fit is not positive or the feature
#' is routed to Poisson. With working weights `W = mu / (1 + alpha mu)`,
#' working residual `(y - mu) / mu`,
#' `P0 = W - W X (X'WX)^{-1} X'W` and `B B' = Q^{-1}`,
#' the statistic is `T = |B'A'P0 r|^2` and its null law is the mixture
#' `sum(lambda_k chi^2_1)` over the positive eigenvalues of `B'A'P0AB`
#' (Davies, then the saddlepoint when Davies fails). `P0` removes the span of
#' `X`, the constant included, so the field enters as it is: no intercept
#' projection, `Q` stays sparse. All numerical work is compiled
#' (`src/nb_null.cpp`, `src/spde_scan.cpp`).
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
#' @param tol Relative convergence tolerance of the IRLS.
#' @param maxit_outer,maxit_inner Iteration limits of the size search and of
#'   the IRLS.
#' @param spectrum_tol Eigenvalues below `spectrum_tol` times the largest are
#'   dropped from the mixture.
#' @param max_eps,max_iter Accuracy and iteration limit of Davies.
#' @return A data frame with one row per feature: `feature_id`, `family`
#'   (`"negative_binomial"` or `"poisson"`), `size`, `overdispersion_score`,
#'   `log_likelihood`, `converged`, `statistic`, `p_value`, `method`,
#'   `fallback`, `error`. The attribute `"coefficients"` holds the
#'   `ncol(X)`-by-feature null coefficients.
#' @keywords internal
spde.scan <- function(Y, X, A, Q, offset = NULL, poisson = NULL, size = NULL,
                      feature_id = rownames(Y), threads = 1L,
                      chunk_size = 256L, block_size = 32L, alpha_max = 1e4,
                      tol = 1e-10, maxit_outer = 100L, maxit_inner = 100L,
                      spectrum_tol = 1e-12, max_eps = 1e-8, max_iter = 1e5) {
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
    method = NA_character_, fallback = NA, error = NA_character_,
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
      threads = as.integer(threads)
    )
    ok <- is.na(fit$error) & fit$converged
    alpha <- ifelse(ok, fit$alpha, NA_real_)
    score <- mgcvst_spde_score_cpp(
      A, Q, X, Yt, O, fit$beta, alpha, threads = as.integer(threads),
      block_size = as.integer(block_size), spectrum_tol = spectrum_tol
    )
    beta[, rows] <- fit$beta
    out$family[rows] <- ifelse(fit$alpha > 0, "negative_binomial", "poisson")
    out$size[rows] <- fit$size
    out$overdispersion_score[rows] <- fit$score0
    out$log_likelihood[rows] <- fit$loglik
    out$converged[rows] <- ok
    out$error[rows] <- ifelse(!is.na(fit$error), fit$error,
      ifelse(!fit$converged, "The null size search did not converge.", score$error))
    for (j in which(ok & is.na(score$error))) {
      z <- list(statistic = score$statistic[j], lambda = score$lambda[[j]])
      row <- rows[j]
      out$statistic[row] <- z$statistic
      if (!length(z$lambda) || !is.finite(z$statistic)) {
        out$error[row] <- "The score has no finite positive mixture spectrum."
        next
      }
      d <- .mgcvst_marginal_davies(z, max_eps, max_iter)
      out$p_value[row] <- d$p_value
      out$method[row] <- d$method_used
      out$fallback[row] <- d$fallback_used
      if (!is.finite(d$p_value) || d$p_value < 0 || d$p_value > 1) {
        out$p_value[row] <- NA_real_
        out$error[row] <- "Invalid score p-value."
      }
    }
  }
  attr(out, "coefficients") <- beta
  out
}

# The raw field of an spde basis at its coordinates: the sparse projector A and
# the sparse FEM precision Q, with no intercept projection.
.spde_scan_field <- function(basis) {
  .spde_basis_validate(basis)
  mesh <- list(xy = basis$mesh_vertices, tv = basis$mesh_triangles)
  coordinates <- as.matrix(basis$coordinates)
  loc <- sweep(coordinates, 2L, basis$transform$center, "-") / basis$transform$scale
  list(A = .spde_basis_project(mesh, loc),
       Q = .spde_fem_precision(.spde_basis_fem(mesh), basis$kappa_internal))
}
