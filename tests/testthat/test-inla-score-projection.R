.projection_score_fixture <- function(vp_multiplier = 1) {
  set.seed(9021)
  n <- 18L
  g <- rep(1 / n, n)
  Z <- qr.Q(qr(matrix(g, ncol = 1L)), complete = TRUE)[, -1L, drop = FALSE]
  D <- seq(0.7, 1.4, length.out = n)
  e <- sin(seq_len(n) / 2) + seq_len(n) / 50
  X <- matrix(1, n, 1L, dimnames = list(NULL, "(Intercept)"))
  base <- mgcvST:::.rkhs_score_operator_factor(
    Z, D, matrix(numeric(), n, 0L), field_scale = 1, B = NULL, Q = NULL
  )
  WX <- mgcvST:::.mgcvst_model_vsolve(base, X)
  Vp_gls <- solve(crossprod(X, WX))
  Vp <- vp_multiplier * Vp_gls
  fit <- list(
    score_backend = "sparse",
    score_sparse = list(
      A = Matrix::Diagonal(n), Q = Matrix::Diagonal(n), constraint = g,
      projection = Z, coefficient_factor = Z, target = "global", sp_index = 1L
    ),
    dispersion = stats::setNames(1, "gene"),
    smoothing_parameters = matrix(1, 1L, 1L,
      dimnames = list("gene", "global")),
    working_error = matrix(e, ncol = 1L, dimnames = list(NULL, "gene")),
    working_variance = matrix(D, ncol = 1L, dimnames = list(NULL, "gene")),
    nuisance_covariance = list(gene = Vp),
    geometry = list(
      nuisance_design = X, X = X,
      target = c(global = 1L),
      smooth = list(list(
        fixed = FALSE, score_component = "global", B = Z,
        penalties = list(diag(ncol(Z))), sp_index = 1L
      ))
    )
  )
  list(fit = fit, Z = Z, D = D, e = e, X = X, Vp_gls = Vp_gls)
}

.projection_vp_operator <- function(F, D, X, Vp) {
  vsolve <- mgcvST:::.rkhs_score_operator_factor(
    F, D, matrix(numeric(), nrow(F), 0L),
    field_scale = 1, B = NULL, Q = NULL
  )
  WN <- mgcvST:::.mgcvst_model_vsolve(vsolve, X)
  structure(list(
    n = nrow(F), vsolve = vsolve, LN = X, VpN = Vp,
    WN = WN, WN_VpN = WN %*% Vp
  ), class = "mgcvst_vp_score_operator")
}

.projection_factor_state <- function(F, D, X, Vp, e, covariance_factor = F) {
  # P is fixed by the fitted null covariance.  The candidate score factor can
  # be raw or centred without changing that nuisance operator.
  operator <- .projection_vp_operator(covariance_factor, D, X, Vp)
  Pe <- mgcvST:::.mgcvst_model_apply_P(operator, e)
  PF <- mgcvST:::.mgcvst_model_apply_P(operator, F)
  list(a = as.numeric(crossprod(F, Pe)), M = crossprod(F, PF),
       P1 = mgcvST:::.mgcvst_model_apply_P(operator, rep(1, nrow(F))))
}

test_that("conditioned score removes one mean direction, not spatial signal", {
  f <- .projection_score_fixture()
  expect_lt(max(abs(colMeans(f$Z))), 1e-14)
  state <- mgcvST:::.mgcvst_model_score_state(f$fit, 1L)
  expect_gt(sum(abs(state$a)), 1e-6)
  expect_gt(sum(abs(state$M)), 1e-6)
  expect_gt(max(eigen(state$M, symmetric = TRUE, only.values = TRUE)$values), 1e-6)
})

test_that("conditioned sparse and dense production states agree", {
  f <- .projection_score_fixture(0.973)
  sparse <- mgcvST:::.mgcvst_model_score_state(f$fit, 1L)
  dense_fit <- f$fit
  dense_fit$score_backend <- "dense"
  dense_fit$score_sparse <- NULL
  dense_fit$.mgcvst_fixed_factors <- mgcvST:::.mgcvst_model_fixed_factors(dense_fit)
  dense <- mgcvST:::.mgcvst_model_score_state(dense_fit, 1L)
  expect_equal(sparse$a, dense$a, tolerance = 2e-12)
  expect_equal(sparse$M, dense$M, tolerance = 2e-12)
})

test_that("centering an already conditioned factor is an algebraic identity", {
  f <- .projection_score_fixture(0.973)
  C <- diag(nrow(f$Z)) - matrix(1 / nrow(f$Z), nrow(f$Z), nrow(f$Z))
  original <- .projection_factor_state(
    f$Z, f$D, f$X, 0.973 * f$Vp_gls, f$e
  )
  centered <- .projection_factor_state(
    C %*% f$Z, f$D, f$X, 0.973 * f$Vp_gls, f$e
  )
  expect_equal(C %*% f$Z, f$Z, tolerance = 1e-14)
  expect_equal(centered$a, original$a, tolerance = 1e-12)
  expect_equal(centered$M, original$M, tolerance = 1e-12)
})

test_that("raw centering is invariant only for exact GLS nuisance covariance", {
  f <- .projection_score_fixture()
  n <- nrow(f$Z)
  C <- diag(n) - matrix(1 / n, n, n)
  constant_loadings <- seq(0.15, 0.45, length.out = ncol(f$Z))
  raw <- f$Z + tcrossprod(rep(1, n), constant_loadings)
  centered <- C %*% raw

  perturbed_vp <- 0.973 * f$Vp_gls
  raw_perturbed <- .projection_factor_state(
    raw, f$D, f$X, perturbed_vp, f$e, covariance_factor = f$Z
  )
  centered_perturbed <- .projection_factor_state(
    centered, f$D, f$X, perturbed_vp, f$e, covariance_factor = f$Z
  )
  expect_gt(max(abs(raw_perturbed$P1)), 1e-6)
  expect_gt(max(abs(raw_perturbed$a - centered_perturbed$a)), 1e-6)
  expect_gt(max(abs(raw_perturbed$M - centered_perturbed$M)), 1e-6)

  raw_gls <- .projection_factor_state(
    raw, f$D, f$X, f$Vp_gls, f$e, covariance_factor = f$Z
  )
  centered_gls <- .projection_factor_state(
    centered, f$D, f$X, f$Vp_gls, f$e, covariance_factor = f$Z
  )
  expect_lt(max(abs(raw_gls$P1)), 1e-12)
  expect_equal(raw_gls$a, centered_gls$a, tolerance = 1e-11)
  expect_equal(raw_gls$M, centered_gls$M, tolerance = 1e-11)
})
