# Remove the tested spatial block while preserving the nuisance GAM setup.
.mgcvst_null_score_setup <- function(G, target_index, null_spec = NULL) {
  target <- G$smooth[[target_index]]
  target_columns <- seq.int(target$first.para, target$last.para)
  keep_columns <- setdiff(seq_len(ncol(G$X)), target_columns)
  X <- attr(G, "training_X")
  if (is.null(X)) X <- as.matrix(G$X)
  X0 <- as.matrix(X[, keep_columns, drop = FALSE])
  if (is.null(null_spec)) {
    data <- as.data.frame(X0)
    if (ncol(data)) names(data) <- paste0(".mgcvST_null_X", seq_len(ncol(data)))
    data$.mgcvST_null_response <- numeric(nrow(X0))
    rhs <- paste(names(data)[-ncol(data)], collapse = " + ")
    formula <- stats::as.formula(paste(
      ".mgcvST_null_response ~ 0", if (nzchar(rhs)) paste("+", rhs) else "+ 1"
    ))
    null_spec <- list(formula = formula, data = data, offset = G$offset,
                      response = ".mgcvST_null_response")
  }
  colnames(X) <- G$term.names
  if (!is.null(null_spec$X0)) X0 <- null_spec$X0
  list(
    spec = null_spec, X = X, X0 = X0,
    single_parametric = ncol(X0) == 1L &&
      !length(mgcv::interpret.gam(null_spec$formula)$smooth.spec),
    smooth = G$smooth, target_S = G$S[[target$first.sp]],
    term.names = G$term.names, target_index = target_index,
    target_columns = target_columns,
    keep_columns = keep_columns
  )
}

# BAM's one-column QR update drops the matrix dimension in mgcv 1.9-4.
.mgcvst_fit_null <- function(setup, data, family, offset, control, gam_args) {
  if (setup$single_parametric) {
    if (!is.null(gam_args$rho) && !isTRUE(gam_args$rho == 0)) {
      stop("A one-coefficient parametric null cannot preserve BAM's nonzero rho in GAM.")
    }
    gam_args[c("rho", "AR.start", "chunk.size", "gc.level", "use.chol",
               "cluster", "samfrac", "discrete", "nthreads", "coef")] <- NULL
  }
  args <- c(list(formula = setup$spec$formula, data = data, family = family,
                 offset = offset, control = control), gam_args)
  if (setup$single_parametric) {
    do.call(mgcv::gam, c(args, list(method = "REML", discrete = FALSE)))
  } else {
    do.call(mgcv::bam, c(args, list(method = "fREML", discrete = TRUE,
                                  nthreads = 1L)))
  }
}

# Evaluate the spatial score from a pure null GAM and the prepared target block.
.mgcvst_null_score_spectrum <- function(null_fit, setup) {
  working <- .mgcvst_marginal_working(null_fit)
  X0 <- as.matrix(setup$X0)
  B <- setup$X[, setup$target_columns, drop = FALSE]
  S <- as.matrix(setup$target_S)
  theta <- CppMatrix::matrixGeneralizedInverse(S / norm(S, "f"))
  F <- CppMatrix::matrixMultiply(
    B, .mgcvst_marginal_matrixsqrt(theta)$w
  )
  V0 <- working$V_phi
  W0 <- X0 / V0
  W0_Vp <- CppMatrix::matrixMultiply(W0, null_fit$Vp)
  P0_apply <- function(v) {
    v_matrix <- is.matrix(v)
    v <- if (v_matrix) v else matrix(v, ncol = 1L)
    Dv <- v / V0
    adjustment <- CppMatrix::matrixMultiply(
      W0_Vp,
      CppMatrix::matrixMultiply(X0, Dv, transA = TRUE)
    )
    out <- Dv - adjustment
    if (v_matrix) out else as.vector(out)
  }
  offset <- if (is.null(null_fit$offset)) numeric(nrow(X0)) else null_fit$offset
  error <- working$pseudo_response - offset
  a <- CppMatrix::matrixMultiply(F, P0_apply(error), transA = TRUE)
  M <- CppMatrix::matrixMultiply(F, P0_apply(F), transA = TRUE)
  lambda <- eigen(M, symmetric = TRUE, only.values = TRUE)$values
  list(
    statistic = sum(a^2), lambda = lambda,
    smooth.term = setup$smooth[[setup$target_index]]$label
  )
}

.mgcvst_null_score_test <- function(null_fit, setup, method = "davies",
                                    max_eps = 1e-8, max_iter = 1e5) {
  method <- match.arg(method, "davies")
  z <- .mgcvst_null_score_spectrum(null_fit, setup)
  result <- .mgcvst_marginal_davies(z, max_eps, max_iter)
  out <- data.frame(smooth.term = z$smooth.term, smooth.pvalue = result$p_value,
                    method = result$method_used)
  attr(out, "marginal_spectrum") <- list(
    statistic = z$statistic, lambda = z$lambda, smooth.term = z$smooth.term
  )
  out
}

# Working pseudo-response and variance of a null GAM (standard and extended
# families), adapted from extract_pseudo_response() and gammfast_working() of
# mgcv.taps (inst/COPYRIGHTS).
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

# Preserve the upstream square-root arithmetic, including the original
# eigensolver, symmetrization order and spectral cutoff.
.mgcvst_marginal_matrixsqrt <- function(A) {
  fit <- CppMatrix::matrixEigen(t(A) / 2 + A / 2)
  d <- c(fit$value)
  d[d <= 0] <- 0
  d <- sqrt(d)
  list(w = CppMatrix::matrixMultiply(fit$vector, t(fit$vector) * d))
}

# Saddlepoint upper tail of Q = sum(lambda_i chi^2_1) (Kuonen 1999, Biometrika
# 86:929, Barndorff-Nielsen form). K(t) = -1/2 sum log(1 - 2 lambda t); the
# relative error stays bounded in the extreme right tail (Chen and Lumley 2019,
# CSDA 139:75), so it is the Stage 1 fallback when Davies fails there.
# At the mean (t = 0) the formula is 0/0 and its limit
# 1/2 - rho3 / (6 sqrt(2 pi)), rho3 = kappa3 / kappa2^(3/2), is used.
.mgcvst_marginal_saddlepoint <- function(q, lambda) {
  lambda <- lambda[lambda > 0]
  if (!length(lambda) || !is.finite(q)) return(NA_real_)
  if (q <= 0) return(1)
  Kp <- function(t) sum(lambda / (1 - 2 * lambda * t)) - q
  mean_q <- sum(lambda)
  centre <- 0.5 - 8 * sum(lambda^3) / (2 * sum(lambda^2))^1.5 / (6 * sqrt(2 * pi))
  if (abs(q - mean_q) < 1e-6 * mean_q) return(centre)
  t_hat <- if (q > mean_q) {
    stats::uniroot(Kp, c(0, (1 - 1e-12) / (2 * max(lambda))), tol = 1e-15)$root
  } else {
    lo <- -1 / (2 * max(lambda))
    while (Kp(lo) > 0) lo <- 2 * lo
    stats::uniroot(Kp, c(lo, 0), tol = 1e-15)$root
  }
  K <- -0.5 * sum(log1p(-2 * lambda * t_hat))
  K2 <- sum(2 * lambda^2 / (1 - 2 * lambda * t_hat)^2)
  w <- sign(t_hat) * sqrt(max(2 * (t_hat * q - K), 0))
  v <- t_hat * sqrt(K2)
  if (!is.finite(w) || !is.finite(v) || w == 0 || v / w <= 0) return(centre)
  stats::pnorm(w + log(v / w) / w, lower.tail = FALSE)
}

# Stage 1 calibration: Davies, then the saddlepoint when Davies fails. Davies
# fails if it errors, returns a missing or non-finite p-value, or returns
# Qq <= 0 or Qq > 1. ifault is reported but does not decide the method.
.mgcvst_marginal_davies <- function(z, max_eps, max_iter) {
  d <- tryCatch(CompQuadForm::davies(q = z$statistic, lambda = z$lambda,
                                    lim = max_iter, acc = max_eps),
                error = function(e) e)
  error <- inherits(d, "condition")
  failed <- error || length(d$Qq) != 1L || !is.finite(d$Qq) ||
    d$Qq <= 0 || d$Qq > 1
  ifault <- if (error || length(d$ifault) != 1L) NA_integer_ else as.integer(d$ifault)
  if (!failed) {
    return(list(p_value = as.numeric(d$Qq), method_used = "davies",
                fallback_used = FALSE, fallback_reason = NA_character_,
                davies_ifault = ifault, error_message = NA_character_))
  }
  reason <- if (error) conditionMessage(d) else
    paste0("Davies numerical failure: Qq=", paste(d$Qq, collapse = ","),
           ", ifault=", paste(d$ifault, collapse = ","))
  list(p_value = .mgcvst_marginal_saddlepoint(z$statistic, z$lambda),
       method_used = "saddlepoint", fallback_used = TRUE,
       fallback_reason = reason, davies_ifault = ifault,
       error_message = NA_character_)
}
