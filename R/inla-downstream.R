#' The pair test behind spdeST.test()
#'
#' The method of [spdeST.test()], which documents the arguments and the
#' result. The sparse score state stored on the fit is reused; no model is
#' refitted. The exact route reconstructs each feature's reduced curvature in
#' the observation-kernel basis; its pair cost is cubic in `q`, and the
#' PCAlearning route does not depend on `q` per pair. The dimension of the
#' field, the number of basis directions and the basis time are stored in
#' `timing$inla_projection`, the pair-stage timings in `timing$pair_pipeline`
#' (exact route) or `timing$pcalearning`, and the chosen route in
#' `timing$route`.
#'
#' @inheritParams spdeST.test
#' @param fitinlaST A fit with spatial models. Only the features with a
#'   spatial model are tested; `pairs = NULL` covers the available features.
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
#' @return An `mgcvST_test` object; see [spdeST.test()].
#' @keywords internal
inlaST.test <- function(
    fitinlaST, pairs = NULL, q.value = 0.05,
    adjust = c("BY", "BH", "Sidak", "none"),
    threads = NULL, chunk_size = NULL, checkpoint_dir = NULL,
    resume = TRUE, verbose = FALSE, moments,
    rank = .mgcvst_pca_defaults$rank,
    n_per_cell = .mgcvst_pca_defaults$n_per_cell,
    seed = .mgcvst_pca_defaults$seed, k = NULL, screen = FALSE) {
  if (missing(moments)) .mgcvst_moments_missing()
  adjust <- match.arg(adjust)
  .mgcvst_test_run(
    fitinlaST, pairs, q.value, adjust, threads, chunk_size,
    checkpoint_dir, resume, verbose, moments = moments, rank = rank,
    n_per_cell = n_per_cell, seed = seed, k = k, screen = screen,
    call = match.call()
  )
}
