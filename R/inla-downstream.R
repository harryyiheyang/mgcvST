#' Test cross-feature spatial covariance from an INLA fit
#'
#' Pairwise covariance test for fits returned by [inlaST.estimate()]. The
#' sparse INLA score state stored on the fit is reused; no model is refitted.
#' For sparse INLA fits, both pair methods use a constrained observation-kernel
#' basis that retains at least 0.995 of its eigenvalue sum. This is a
#' pairwise-test approximation; the fitted sparse field is unchanged. Basis and
#' pair-stage timings of the Liu methods are stored in `timing$inla_projection`.
#'
#' @param fitinlaST An object returned by [inlaST.estimate()].
#' @param pairwise_method `"score_liu"` tests the squared cross-gene score with Liu
#'   moment matching. `"conditional_cauchy"` evaluates the two
#'   conditional-normal directions and combines their p-values by the
#'   equal-weight Cauchy rule; see Details. With `pairs = NULL`, both the
#'   conditional method and `score_liu` with `liu_approximation = "exact"`
#'   test every available gene pair.
#' @param liu_approximation Trace evaluation for `pairwise_method = "score_liu"`.
#'   `"exact"` computes the four Liu trace moments in the reduced
#'   observation-kernel coordinates. `"pca_learning"` projects every score
#'   covariance `H_j` onto a rank-`rank` orthonormal basis learned from
#'   stratified training genes; see Details.
#' @param rank Number of basis matrices retained by
#'   `liu_approximation = "pca_learning"`.
#' @param n_per_cell Training genes drawn per stratification cell by
#'   `liu_approximation = "pca_learning"`.
#' @param seed Non-negative integer seed for PCAlearning training-gene sampling;
#'   the caller's random-number state is restored.
#' @param checkpoint_dir Optional checkpoint directory. The conditional method
#'   saves per-gene variance rows; the Liu methods save reusable score states
#'   and pair batches. With `NULL`, temporary storage is removed on exit.
#' @param resume Reuse compatible completed checkpoint entries.
#' @param conditional_precision Precision used for conditional variance
#'   multiplication: `"double"` or `"float32"`. Scores and p-values remain in
#'   double precision.
#' @inheritParams mgcvST.test
#' @param calibration Only `"liu"` is supported for INLA fits.
#' @param chunk_size For conditional pairs, the number of genes materialized
#'   together (1 to 64). For Liu pairs, the maximum number of tested pairs per
#'   native batch; the adaptive memory budget may reduce distinct features.
#' @param BPPARAM Compatibility argument; only `SerialParam()` is accepted.
#' @param threads Positive number of OpenMP threads for sparse INLA feature
#'   preparation, reduced materialization and pair batches. `NULL` uses one thread.
#' @details With `liu_approximation = "pca_learning"`, each gene receives the
#'   variance scales `sigma_g2 = dispersion / lambda` and
#'   `sigma_e2 = 1 + mean(mu) / theta` for negative-binomial genes (1 for
#'   Poisson genes), with `mu` recovered from the working variance. Genes are
#'   stratified into 10 quantile bins of `log(sigma_g2)` crossed with one
#'   Poisson bin and 9 quantile bins of `log(sigma_e2)`, and `n_per_cell` genes
#'   are drawn per cell, with the quota of sparse cells reallocated
#'   proportionally to cell size. The training matrices
#'   `tau_j H_j`, `tau_j = sigma_e2 / sigma_g2`, define an orthonormal basis `B`
#'   through the eigen decomposition of their Gram matrix. Every tested gene is
#'   summarized by its basis coefficients `c_j`, and its residual
#'   `e_j^2 = ||H_j||_F^2 - ||c_j||^2` is reported. Pair traces
#'   `tr((H_i H_j)^s)`, `s = 1, ..., 4`, are those of the projected matrices.
#'   Liu p-values are computed on the log scale and returned in the result
#'   columns `log_p_two_sided`, `log_p_positive` and `log_p_negative`. With
#'   `method = "BY"`, the step-up adjustment is applied to these log p-values,
#'   so decisions remain defined when p-values underflow; `information` and
#'   `effective_rank` are not computed on this path. The result element
#'   `pca_learning` stores the training genes, the Gram eigenvalues and
#'   rotation defining `B`, the coefficients `c_j`, the per-gene table with
#'   `e2_relative = e_j^2 / ||H_j||_F^2`, and the stage timings, including the
#'   trace-table precomputation.
#'
#'   With `pairwise_method = "conditional_cauchy"`, let `S_ij = a_i' a_j` and
#'   `v_(i|j) = a_j' M_i a_j`. Under independent Gaussian null scores,
#'   `S_ij | a_j` is normal with variance `v_(i|j)`, so each directional
#'   two-sided normal p-value is exactly uniform. The directions are combined
#'   using `T = (tan((0.5-p_(i|j))*pi) + tan((0.5-p_(j|i))*pi))/2` and the
#'   standard Cauchy upper tail. Because `T` cannot exceed its larger component,
#'   `p_ij >= min(p_(i|j), p_(j|i))` and the null rejection probability is at
#'   most `2*alpha` under any dependence.
#'
#'   For `X = 1/p_(i|j)` and `Y = 1/p_(j|i)`, both directional tails satisfy
#'   `Pr(X > x) = Pr(Y > x) = 1/x`; in the far tail, `1/p_ij` approaches
#'   `(X+Y)/2`. If both directions become extreme together, their normal
#'   z-scores obey `z_2 = R*z_1`, where `R = sqrt(v_(i|j)/v_(j|i))`.
#'   Under a continuous, nondegenerate distribution of `R`, unequal extremes
#'   occur together only when `abs(R-1)` is of order `1/log(1/alpha)`.
#'   Thus the combined tail approaches the nominal tail as `alpha` tends to
#'   zero under this condition.
#'
#'   The conditional results contain `signed_score`, `statistic`
#'   (`signed_score^2`), the combined `p_two_sided`, the directional
#'   `p_1_given_2` (variance `v_(1|2)`) and `p_2_given_1`, and their natural-log
#'   versions. Multiple testing follows `FDR`, `method` and `q.value` as for
#'   Liu pairs (`p_adjusted`, `log_p_adjusted`, `discovered`); `method = "BY"`
#'   is applied to the log p-values, so tails below the floating-point range
#'   keep their ordering and decisions. `highlight` is unavailable for this
#'   method.
#' @return With `pairwise_method = "score_liu"` and `liu_approximation =
#'   "exact"` (the default), a compact list with integer `i`, `j` and double
#'   `score`, `mlog10p` pairs (materialized in `$result` for an explicit
#'   `pairs` block, or as Parquet `$shards` with `pairs = NULL` or a
#'   checkpointed explicit block), plus `$feature_id`, `$failed` (genes whose
#'   fp16 state could not be built) and `$bh` (BH results); see
#'   `.mgcvst_inla_fp16_run()`. Otherwise, an `mgcvST_test` object.
#' @export
inlaST.test <- function(
    fitinlaST,
    pairwise_method = c("score_liu", "conditional_cauchy"),
    liu_approximation = c("exact", "pca_learning"),
    rank = 10L, n_per_cell = 3L, seed = 1L,
    checkpoint_dir = NULL, resume = TRUE,
    conditional_precision = c("double", "float32"),
    q.value = 0.05, FDR = TRUE, method = "BH",
    pairs = NULL, highlight = NULL, calibration = "liu",
    chunk_size = NULL, threads = NULL, verbose = FALSE,
    BPPARAM = BiocParallel::SerialParam(), ...) {
  pairwise_method <- match.arg(pairwise_method)
  liu_approximation <- match.arg(liu_approximation)
  conditional_precision <- match.arg(conditional_precision)
  if (pairwise_method == "score_liu" && liu_approximation == "exact") {
    if (!is.null(highlight)) {
      stop("highlight is unavailable for the compact fp16 score_liu result; ",
           "filter its i/j columns (or Parquet shards) directly.")
    }
    if (!identical(calibration, "liu")) {
      stop("The exact fp16 score_liu path uses calibration = 'liu' only.")
    }
    if (length(list(...))) stop("Unused arguments in ... for score_liu pairs.")
    .mgcvst_inla_serial_backend(BPPARAM)
    return(.mgcvst_inla_fp16_run(
      fitinlaST, pairs = pairs, checkpoint_dir = checkpoint_dir, resume = resume,
      threads = if (is.null(threads)) 1L else threads,
      chunk_size = if (is.null(chunk_size)) 4000000L else chunk_size,
      verbose = verbose, q.value = q.value, FDR = FDR, method = method
    ))
  }
  if (pairwise_method == "conditional_cauchy") {
    if (liu_approximation != "exact") {
      stop("liu_approximation requires pairwise_method = 'score_liu'.")
    }
    if (!is.null(highlight)) {
      stop("highlight is unavailable for conditional pairwise results.")
    }
    if (!identical(calibration, "liu")) {
      stop("Conditional pairs do not use a non-Liu calibration argument.")
    }
    if (length(list(...))) stop("Unused arguments in ... for conditional pairs.")
    .mgcvst_inla_serial_backend(BPPARAM)
    return(.mgcvst_conditional_test(
      fitinlaST, pairs, q.value, FDR, method, threads, chunk_size,
      checkpoint_dir, resume, conditional_precision, match.call()
    ))
  }
  if (!identical(conditional_precision, "double")) {
    stop("conditional_precision requires pairwise_method = 'conditional_cauchy'.")
  }
  # Only liu_approximation = "pca_learning" reaches this point; "exact" is
  # handled above for both pairs = NULL and an explicit pair block.
  engine <- .mgcvst_test_engine(fitinlaST)
  if (!.mgcvst_inla_downstream(fitinlaST)) {
    stop("liu_approximation = 'pca_learning' requires a sparse INLA fit.")
  }
  engine(
    fitmgcvST = fitinlaST, q.value = q.value, FDR = FDR, method = method,
    BPPARAM = BPPARAM, ..., pairs = pairs, highlight = highlight,
    calibration = calibration, chunk_size = chunk_size,
    threads = threads, verbose = verbose,
    checkpoint_dir = checkpoint_dir, resume = resume,
    liu_approximation = liu_approximation, rank = rank,
    n_per_cell = n_per_cell, seed = seed
  )
}
