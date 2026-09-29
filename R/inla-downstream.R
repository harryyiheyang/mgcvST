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
#' @param fitinlaST An object returned by [inlaST.estimate()].
#' @param approximate_test `TRUE` (the default) evaluates the four Liu trace
#'   moments by PCAlearning: every score covariance `H_j` is projected onto a
#'   rank-`rank` orthonormal basis learned from stratified training genes; see
#'   Details. `FALSE` computes the exact Liu trace moments in the reduced
#'   observation-kernel coordinates with the fp16 score-state backend.
#' @param rank Number of PCAlearning basis matrices. Unused when
#'   `approximate_test = FALSE`.
#' @param n_per_cell Training genes drawn per PCAlearning stratification cell.
#'   Unused when `approximate_test = FALSE`.
#' @param seed Non-negative integer seed for PCAlearning training-gene sampling;
#'   the caller's random-number state is restored. Unused when
#'   `approximate_test = FALSE`.
#' @param checkpoint_dir Optional checkpoint directory for reusable score
#'   states and pair batches. With `NULL`, temporary storage is removed on exit.
#' @param resume Reuse compatible completed checkpoint entries.
#' @inheritParams mgcvST.test
#' @param method Multiple-testing adjustment passed to [stats::p.adjust()],
#'   such as `"BH"` or `"BY"`.
#' @param chunk_size Maximum number of tested pairs per native batch; the
#'   adaptive memory budget may reduce distinct features.
#' @param BPPARAM Compatibility argument; only `SerialParam()` is accepted.
#' @param threads Positive number of OpenMP threads for sparse INLA feature
#'   preparation, reduced materialization and pair batches. `NULL` uses one thread.
#' @details With `approximate_test = TRUE`, each gene receives the
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
#'   With `approximate_test = FALSE` and `pairs = NULL`, every available gene
#'   pair is tested.
#' @return With `approximate_test = TRUE`, an `mgcvST_test` object. With
#'   `approximate_test = FALSE`, a compact list with integer `i`, `j` and
#'   double `score`, `mlog10p` pairs (materialized in `$result` for an explicit
#'   `pairs` block, or as Parquet `$shards` with `pairs = NULL` or a
#'   checkpointed explicit block), plus `$feature_id`, `$failed` (genes whose
#'   fp16 state could not be built) and `$bh` (BH results); see
#'   `.mgcvst_inla_fp16_run()`.
#' @export
inlaST.test <- function(
    fitinlaST,
    approximate_test = TRUE,
    rank = 10L, n_per_cell = 3L, seed = 1L,
    checkpoint_dir = NULL, resume = TRUE,
    q.value = 0.05, FDR = TRUE, method = "BH",
    pairs = NULL, highlight = NULL,
    chunk_size = NULL, threads = NULL, verbose = FALSE,
    BPPARAM = BiocParallel::SerialParam(), ...) {
  if (!is.logical(approximate_test) || length(approximate_test) != 1L ||
      is.na(approximate_test)) {
    stop("approximate_test must be TRUE or FALSE.")
  }
  if (!approximate_test) {
    if (!is.null(highlight)) {
      stop("highlight is unavailable for the compact exact Liu result; ",
           "filter its i/j columns (or Parquet shards) directly.")
    }
    unused <- list(...)
    if (length(unused)) {
      stop("Unused arguments in ...: ", paste(names(unused), collapse = ", "))
    }
    .mgcvst_inla_serial_backend(BPPARAM)
    return(.mgcvst_inla_fp16_run(
      fitinlaST, pairs = pairs, checkpoint_dir = checkpoint_dir, resume = resume,
      threads = if (is.null(threads)) 1L else threads,
      chunk_size = if (is.null(chunk_size)) 4000000L else chunk_size,
      verbose = verbose, q.value = q.value, FDR = FDR, method = method
    ))
  }
  engine <- .mgcvst_test_engine(fitinlaST)
  if (!.mgcvst_inla_downstream(fitinlaST)) {
    stop("approximate_test = TRUE requires a sparse INLA fit.")
  }
  engine(
    fitmgcvST = fitinlaST, q.value = q.value, FDR = FDR, method = method,
    BPPARAM = BPPARAM, ..., pairs = pairs, highlight = highlight,
    calibration = "liu", chunk_size = chunk_size,
    threads = threads, verbose = verbose,
    checkpoint_dir = checkpoint_dir, resume = resume,
    liu_approximation = "pca_learning", rank = rank,
    n_per_cell = n_per_cell, seed = seed
  )
}
