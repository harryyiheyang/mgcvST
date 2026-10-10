#' Test cross-feature spatial covariance from an INLA fit
#'
#' Pairwise covariance test for fits returned by [inlaST.estimate()]. The
#' sparse INLA score state stored on the fit is reused; no model is refitted.
#' Each gene pair is tested by the signed cross-gene score, calibrated by the
#' saddlepoint approximation described in [mgcvST.test()], which also
#' describes the exact and PCAlearning routes and the choice between them. The
#' score uses the full-rank constrained observation-kernel basis (all `q - 1`
#' directions of the constrained field, ordered by eigenvalue), the same basis
#' that [inlaST.wgcna()] uses; no eigenvalue coverage truncation is applied.
#' The exact route reconstructs each gene's reduced curvature in that basis
#' and takes `q` small enough for its cubic pair cost; the PCAlearning route
#' serves large `q`. The fitted sparse field is unchanged. The dimension of the
#' field, the number of basis directions and the basis time are stored in
#' `timing$inla_projection`, the pair-stage timings in `timing$pair_pipeline`
#' (exact route) or `timing$pcalearning`, and the chosen route in
#' `timing$route`.
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
#' @param fitinlaST An object returned by [inlaST.estimate()] (version 0.0.1.9032
#'   or later). Only the features with a spatial model are tested; a feature
#'   that step 2 of the estimation did not fit is unavailable, and `pairs = NULL`
#'   covers the available features.
#' @details The PCAlearning route projects every score covariance `H_j` onto
#'   a rank-`rank` orthonormal basis learned from stratified training genes.
#'   Each gene receives the variance scales `sigma_g2 = dispersion / lambda`
#'   and `sigma_e2 = 1 + mean(mu) / theta` for negative-binomial genes (1 for
#'   Poisson genes), with `mean(mu)` stored at estimation as `mu_bar`. Genes are
#'   stratified into 10 quantile bins of `log(sigma_g2)` crossed with one
#'   Poisson bin and 9 quantile bins of `log(sigma_e2)`, and `n_per_cell` genes
#'   are drawn per cell, with the quota of sparse cells reallocated
#'   proportionally to cell size. The training matrices
#'   `tau_j H_j`, `tau_j = sigma_e2 / sigma_g2`, define an orthonormal basis `B`
#'   through the eigen decomposition of their Gram matrix, which is summed in
#'   a fixed block order so that the result does not depend on the number of
#'   threads. Every tested gene is summarized by its basis coefficients `c_j`,
#'   and its residual `e_j^2 = ||H_j||_F^2 - ||c_j||^2` is reported. The pair
#'   traces `tr(H_i H_j) = c_i' c_j` and `tr((H_i H_j)^2)` are those of the
#'   projected matrices, which the level-2 trace table of `B` gives by
#'   contraction of the degree-2 monomials of the coefficients. The leading
#'   singular values of each pair come from the shared basis `V`
#'   (the `k` leading eigenvectors of the sum of the training matrices, each
#'   divided by its largest absolute entry): with `R_j` the upper Cholesky
#'   factor of `V' H_j V` divided by that entry, they are the singular values
#'   of `R_i R_j'`. The result element
#'   `pca_learning` stores the training genes, the Gram eigenvalues and
#'   rotation defining `B`, the coefficients `c_j`, the per-gene table with
#'   `e2_relative = e_j^2 / ||H_j||_F^2`, and the stage timings, including the
#'   trace-table precomputation.
#' @return An `mgcvST_test` object; see [mgcvST.test()].
#' @export
inlaST.test <- function(
    fitinlaST, pairs = NULL, q.value = 0.05,
    adjust = c("BY", "BH", "Sidak", "none"),
    threads = NULL, chunk_size = NULL, checkpoint_dir = NULL,
    resume = TRUE, verbose = FALSE,
    moments = c("auto", "exact", "pcalearning"),
    rank = .mgcvst_pca_defaults$rank,
    n_per_cell = .mgcvst_pca_defaults$n_per_cell,
    seed = .mgcvst_pca_defaults$seed, k = NULL) {
  adjust <- match.arg(adjust)
  moments <- match.arg(moments)
  .mgcvst_test_run(
    fitinlaST, "inla", pairs, q.value, adjust, threads, chunk_size,
    checkpoint_dir, resume, verbose, moments = moments, rank = rank,
    n_per_cell = n_per_cell, seed = seed, k = k, call = match.call()
  )
}
