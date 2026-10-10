#' Test cross-feature spatial covariance from an INLA fit
#'
#' Pairwise covariance test for fits returned by [inlaST.estimate()]. The
#' sparse INLA score state stored on the fit is reused; no model is refitted.
#' Each gene pair is tested by the squared cross-gene score with Liu moment
#' matching, the only pair test for INLA fits. The score uses a constrained
#' observation-kernel basis that retains at least 0.995 of its eigenvalue sum.
#' This is a pairwise-test approximation; the fitted sparse field is unchanged.
#' Basis and pair-stage timings are stored in `timing$inla_projection`.
#'
#' The result has the compact shape described in [mgcvST.test()]: integer
#' feature indices `i < j`, the signed `score`, the natural-log two-sided,
#' positive and negative p-values, the adjusted two-sided `log_q`, the
#' calibration's `remainder_kind` and a `status` code, streamed to Parquet
#' shards. With `pairs = NULL` every pair of the available genes is tested, one
#' block of left genes at a time. The multiple-testing adjustment is applied
#' once to the two-sided family, in log space, and positive and negative
#' discoveries are split by the sign of the score.
#'
#' @inheritParams mgcvST.test
#' @param fitinlaST An object returned by [inlaST.estimate()].
#' @param rank Number of PCAlearning basis matrices.
#' @param n_per_cell Training genes drawn per PCAlearning stratification cell.
#' @param seed Non-negative integer seed for PCAlearning training-gene sampling;
#'   the caller's random-number state is restored.
#' @details The four Liu trace moments are evaluated by PCAlearning: every score
#'   covariance `H_j` is projected onto a rank-`rank` orthonormal basis learned
#'   from stratified training genes. Each gene receives the
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
#'   The result element `pca_learning` stores the training genes, the Gram
#'   eigenvalues and rotation defining `B`, the coefficients `c_j`, the
#'   per-gene table with `e2_relative = e_j^2 / ||H_j||_F^2`, and the stage
#'   timings, including the trace-table precomputation.
#' @return An `mgcvST_test` object; see [mgcvST.test()].
#' @export
inlaST.test <- function(
    fitinlaST, pairs = NULL, q.value = 0.05,
    adjust = c("BY", "BH", "Sidak", "none"),
    rank = .mgcvst_pca_defaults$rank,
    n_per_cell = .mgcvst_pca_defaults$n_per_cell,
    seed = .mgcvst_pca_defaults$seed,
    threads = NULL, chunk_size = NULL, checkpoint_dir = NULL,
    resume = TRUE, verbose = FALSE) {
  adjust <- match.arg(adjust)
  .mgcvst_test_run(
    fitinlaST, "pcalearning", pairs, q.value, adjust, threads, chunk_size,
    checkpoint_dir, resume, verbose, rank = rank, n_per_cell = n_per_cell,
    seed = seed, call = match.call()
  )
}
