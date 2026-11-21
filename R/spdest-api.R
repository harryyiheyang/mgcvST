# spdeST user interface.
#
#   model <- spdeST.set(...)                         # every fixed setting
#   scan  <- spdeST.scan(Y, model, ...)              # Stage 1
#   fit   <- spdeST.estimate(scan, Y, ...)           # Stage 2, chosen features
#   res   <- spdeST.test(fit, ...)                   # pair test

.spdest_flag <- function(value, name) {
  if (!is.logical(value) || length(value) != 1L || is.na(value)) {
    stop(name, " must be TRUE or FALSE.")
  }
  value
}

# Only the arguments the caller supplied, so that an implementation keeps its
# own defaults for the others.
.spdest_supplied <- function(...) {
  args <- list(...)
  args[!vapply(args, is.null, logical(1L))]
}

# spde.scan() controls supplied by name; the others keep its defaults.
.spdest_scan_control <- function(scan_control) {
  if (!is.list(scan_control)) stop("scan_control must be a named list.")
  if (!length(scan_control)) return(list())
  allowed <- c("tol", "maxit_outer", "maxit_inner", "alpha_max", "alpha_min",
               "step_bound", "step_fallback", "maxit_halving", "eta_max", "rank_tol",
               "spectrum_tol", "spa_centre_tol", "spa_root_tol",
               "spa_root_maxit", "block_size")
  bad <- setdiff(names(scan_control), allowed)
  if (is.null(names(scan_control)) || any(!nzchar(names(scan_control))) || length(bad)) {
    stop("scan_control takes named spde.scan() controls: ",
         paste(allowed, collapse = ", "), ".")
  }
  scan_control
}

#' Stage 1 spatial score scan
#'
#' Tests every feature for a spatial SPDE variance component against its null
#' model. The model from [spdeST.set()] carries every fixed setting. Stage 1
#' fits GLM nulls by maximum likelihood (NB2, or Poisson when the overdispersion score
#' `sum((y - mu)^2 - y) / 2` is not positive) and the score test of the raw
#' sparse SPDE field with the null projection `P0`, calibrated by the
#' saddlepoint approximation ([spde.scan()]). The null is a GLM: covariates
#' enter as linear terms or factors (0-1 columns of the model matrix); a null
#' with smooth terms or random effects is refused.
#'
#' @param Y Feature-by-observation response matrix.
#' @param model Result of [spdeST.set()].
#' @param offset `NULL`, a length-n vector, or a feature-by-observation matrix
#'   added to the model offset. It is kept for [spdeST.estimate()].
#' @param feature_id Feature identifiers.
#' @param adjust,q.value Multiple-testing adjustment of the Stage 1 p-values
#'   (`"BY"`, `"BH"`, `"Sidak"` or `"none"`) and the q-value that defines the
#'   discoveries of [spdeST.estimate()].
#' @param BPPARAM BiocParallel backend; `NULL` is serial.
#' @param chunk_size Features per task; `NULL` chooses it from `BPPARAM`.
#' @param threads OpenMP threads per worker.
#' @param checkpoint_dir,resume Chunk checkpoints and their reuse.
#' @param scan_control Named list of [spde.scan()] controls, for example
#'   `list(tol = 1e-8, maxit_outer = 25)`.
#' @return An `spdeST_scan` object: the Stage 1 fit with
#'   `diagnostics$marginal_p_value`, `diagnostics$marginal_q_value`,
#'   `diagnostics$family_used` and `stage1`, carrying the model for
#'   [spdeST.estimate()].
#' @export
spdeST.scan <- function(Y, model, offset = NULL, feature_id = rownames(Y),
                        adjust = c("BY", "BH", "Sidak", "none"), q.value = 0.05,
                        BPPARAM = NULL, chunk_size = NULL, threads = 1L,
                        checkpoint_dir = NULL, resume = TRUE,
                        scan_control = list()) {
  if (!inherits(model, "spdeST_model")) stop("model must be returned by spdeST.set().")
  settings <- utils::modifyList(list(retain_smooth = FALSE, diagnostics = FALSE),
                                if (is.null(model$estimation)) list() else model$estimation)
  args <- c(
    list(Y = Y, model = model, feature_id = feature_id, offset = offset,
         retain_smooth = settings$retain_smooth, diagnostics = settings$diagnostics,
         threads = threads, adjust = match.arg(adjust), q.value = q.value,
         resume = resume, scan_control = scan_control),
    .spdest_supplied(BPPARAM = BPPARAM, chunk_size = chunk_size,
                     checkpoint_dir = checkpoint_dir)
  )
  do.call(.spdest_scan, args)
}

#' Stage 2 spatial estimation of the selected features
#'
#' Fits the spatial model of the features chosen by `features` by
#' [spde.estimate()]: the Laplace marginal likelihood (REML) of each feature
#' is maximised in the precision of the field and the NB size, with the
#' settings fixed in [spdeST.set()] and the offset given to [spdeST.scan()].
#' Stage 1 is not recomputed. The fit keeps the score vectors and the
#' per-feature materials that [spdeST.test()] needs.
#'
#' @param scan Result of [spdeST.scan()], or of an earlier `spdeST.estimate()`
#'   to add features.
#' @param Y The response matrix given to [spdeST.scan()].
#' @param features `"discoveries"` (Stage 1 q-value at most `q.value`),
#'   `"all"`, `"none"`, feature IDs, one-based indices, or one logical per
#'   feature.
#' @param adjust,q.value `NULL` keeps those of [spdeST.scan()].
#' @inheritParams spdeST.scan
#' @return An `spdeST_fit` object for [spdeST.test()].
#' @export
spdeST.estimate <- function(scan, Y, features = "discoveries", adjust = NULL,
                            q.value = NULL, BPPARAM = NULL, chunk_size = NULL,
                            threads = 1L, checkpoint_dir = NULL, resume = TRUE) {
  if (!inherits(scan, "spdeST_scan")) {
    stop("scan must be returned by spdeST.scan().")
  }
  args <- c(
    list(scan, Y = Y, features = features, threads = threads, resume = resume),
    .spdest_supplied(adjust = adjust, q.value = q.value, BPPARAM = BPPARAM,
                     chunk_size = chunk_size, checkpoint_dir = checkpoint_dir)
  )
  do.call(.spdest_estimate, args)
}

#' Pairwise spatial covariance test
#'
#' Tests the cross-feature spatial covariance of every requested feature pair
#' from the compact summaries in a fit returned by [spdeST.estimate()]; no
#' model is refitted. Each pair is tested by the signed cross-feature score
#' `U = sum(s_i x_i y_i)`, where `s_i` are the singular values of the product
#' of the square roots of the two score-covariance matrices. The tail of `U` is
#' a Lugannani-Rice saddlepoint approximation computed in log space on the
#' `k` leading singular values of a basis shared by all features, plus a
#' remainder that matches the remaining power sums of the spectrum. The score
#' uses the full-rank constrained observation-kernel basis (all `q - 1`
#' directions of the constrained field).
#'
#' The argument `moments` selects the route that supplies the shared basis and
#' the remainder; it has no default. With `moments = "exact"`, the shared
#' basis holds the `k = 20` leading eigenvectors of the summed, normalized
#' score covariances, every pair needs the four exact trace moments
#' `tr((H_i H_j)^s)`, `s = 1, ..., 4`, and the remainder is two moment-matched
#' nodes (one node or a Gaussian term when the moments do not allow two nodes).
#' With `moments = "pcalearning"`, the score covariances are projected onto a
#' rank-`rank` basis learned from training genes (see Details of
#' [spdeST.test()]), the pair traces `tr(H_i H_j)` and `tr((H_i H_j)^2)` come from
#' a contraction of the projected coefficients, the shared basis holds the
#' `k = 80` leading eigenvectors of the training genes, and the remainder is
#' one node. The exact route costs time cubic in the score dimension `q` for
#' every pair, and the PCAlearning route does not depend on `q` per pair; the
#' PCAlearning route needs more genes than `rank`. `verbose = TRUE` prints the
#' route and `k`. A checkpoint directory records its route, and a resumed call
#' with another `moments` stops with a message.
#'
#' With `pairs = NULL`, every pair of the available features is tested; the
#' pairs are generated and scored in blocks and never held as one matrix. A
#' feature whose fit failed is not available and is reported in `$failed`.
#' With explicit `pairs`, a pair that contains an unavailable feature is
#' returned with status 3 and missing p-values.
#'
#' The result is compact. Each pair is one row of the integer feature indices
#' `i < j` (positions in `$feature_id`), the signed `score` `U`, its
#' information `t1 = tr(M_i M_j)` (exact on the exact route, `c_i' c_j` on the
#' PCAlearning route), the natural-log
#' two-sided, positive and negative p-values `log_p_two_sided`,
#' `log_p_positive` and `log_p_negative`, the adjusted two-sided
#' log q-value `log_q`, the integer `remainder_kind` of the calibration (0:
#' no remainder; 1: one node; 2: two nodes; 3: Gaussian term) and the integer
#' `status`:
#'
#' | status | meaning | p-values |
#' |---|---|---|
#' | 0 | evaluated | the saddlepoint p-values |
#' | 1 | trace moments non-finite or non-positive | 1 |
#' | 2 | invalid p-value | 1 |
#' | 3 | a feature of the pair has no usable score state | missing, not adjusted |
#' | 4 | screened (`screen = TRUE`) | lower bounds |
#'
#' `U / t1` is the Haseman-Elston-type estimate of the latent spatial
#' correlation of the pair, with null standard error `1 / sqrt(t1)`.
#' `$features` holds, for every feature of a
#' tested pair, `trace_M = tr(M_g)` (the expected `a_g' a_g`) and
#' `trace_M2 = tr(M_g^2)`; the PCAlearning route also keeps the coefficients
#' `c_g` in `$pca_learning$coefficients`.
#'
#' A pair with status 1 or 2 has p = 1 (`log_p_two_sided`, `log_p_positive` and
#' `log_p_negative` are all 0), the conservative choice, and stays in the
#' adjustment family. Rows are written as Parquet shards while the
#' pairs are evaluated, and `$results` is the same table sorted by `(i, j)`
#' when it fits the memory guard (56 bytes per pair, 20% of available memory),
#' and `NULL` otherwise. Without a `checkpoint_dir` the shards are temporary:
#' they are deleted once `$results` is built and `$shards` is empty, and when
#' `$results` is `NULL` they stay in the session's temporary directory and are
#' listed in `$shards`. With a `checkpoint_dir`, `$shards` lists the final
#' files in it. Use a `checkpoint_dir` for runs of many millions of pairs, so
#' that the shards go to a disk of the right size.
#'
#' The adjustment is applied once, to the two-sided family, in log space
#' by the native kernel, so p-values below the double range keep their
#' ordering. `"BY"` is the Benjamini-Yekutieli procedure under arbitrary
#' dependence (`c(m) = sum(1 / seq_len(m))`, as `stats::p.adjust(, "BY")`),
#' `"BH"` the Benjamini-Hochberg step-up, `"Sidak"` the single-step Sidak
#' correction and `"none"` leaves the p-values unadjusted. A pair is a
#' discovery when `log_q <= log(q.value)`; discoveries with a positive score
#' are positive and those with a negative score negative. The adjustment needs
#' the two-sided log p-values and the adjusted values in memory (24 bytes per
#' pair, 40% of available memory); when that does not hold it is skipped, with
#' a warning, and `log_q` is `NA`.
#'
#' @param fit Result of [spdeST.estimate()].
#' @param pairs `NULL` (the default) for every pair of available features, or a
#'   two-column matrix or data frame of feature IDs or one-based indices.
#' @param q.value Discovery threshold on the adjusted q-value, in `(0, 1]`.
#' @param adjust Multiple-testing adjustment of the two-sided family: `"BY"`
#'   (the default), `"BH"`, `"Sidak"` or `"none"`.
#' @param threads Positive number of OpenMP threads for score-state
#'   preparation and the pair kernels. `NULL` uses one.
#' @param chunk_size Maximum number of pairs evaluated per native block and
#'   written per shard. `NULL` uses 10,000 for the exact route and 1,000,000
#'   for the PCAlearning route.
#' @param checkpoint_dir Optional checkpoint directory. Reusable feature score
#'   states and the raw pair shards are saved there, and a repeated call
#'   resumes completed shards. Pair results are keyed by the algorithm
#'   contract; a directory holding pair results from another contract is
#'   refused. With `NULL`, temporary storage is used.
#' @param resume Reuse compatible completed checkpoint entries.
#' @param verbose Logical; report progress and the route.
#' @param moments `"exact"` or `"pcalearning"`, the route of the pair
#'   calibration. Required: there is no default, and a call without it stops.
#' @param rank Number of PCAlearning basis matrices (default 100). Used by the
#'   PCAlearning route.
#' @param n_per_cell Training genes drawn per PCAlearning stratification cell
#'   (default 3).
#' @param seed Non-negative integer seed for PCAlearning training-gene
#'   sampling; the caller's random-number state is restored.
#' @param k Number of leading singular values of the shared basis. `NULL`
#'   uses 20 on the exact route and 80 on the PCAlearning route; a value above
#'   the score dimension `q` is reduced to `q`, and `k = q` reproduces the
#'   full-spectrum saddlepoint of each pair.
#' @param screen Logical; `TRUE` skips the pairs that cannot be discoveries
#'   (both routes). Given one of the two score vectors the pair score is
#'   Gaussian with a variance of mean `t_1`, so for `z = |U| / sqrt(t_1) >= 2`
#'   the two-sided p-value is at least the Gaussian tail `2 * pnorm(-z)`. On
#'   the exact route `t_1 = tr(H_i H_j)`; on the PCAlearning route it is the
#'   projected trace `c_i' c_j`, which the variance of the pair score in that
#'   route's saddlepoint is never below. With `theta` the largest p-value that
#'   `adjust` can reject at `q.value` over the tested pairs (BY:
#'   `q.value / sum(1 / seq_len(m))`), a pair with `2 * pnorm(-z) > 2 * theta`
#'   is not evaluated: it has status 4, its p-values are the lower bound
#'   `2 * pnorm(-max(z, 2))`, and it stays in the adjustment family. The
#'   discoveries, and every `log_q` at or below `log(q.value)`, are those of
#'   `screen = FALSE`; the `log_q` of other pairs are lower bounds. The screen
#'   is off when the bound would be needed below `z = 2` (for example
#'   `adjust = "BH"` with `q.value = 0.05`). `timing$screen` records the
#'   screening value and the number of screened pairs.
#' @details The dimension of the field, the number of basis directions and the
#'   basis time are stored in `timing$score_projection`, the pair-stage timings
#'   in `timing$pair_pipeline` (exact route) or `timing$pcalearning`.
#'
#'   The PCAlearning route projects every score covariance `H_j` onto
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
#' @return An object of class `spdeST_test` with `results`, `shards`,
#'   `feature_id`, `failed` (features without a usable score state and the
#'   reason), `threshold`, `discoveries`, `adjustment`, `timing`
#'   (`timing$route` holds the route, `k` and the score dimension `q`),
#'   `calibration`, `moments` (the route), `contract` and `call`. A PCAlearning
#'   run additionally returns `pca_learning`.
#' @export
spdeST.test <- function(fit, pairs = NULL, q.value = 0.05,
                        adjust = c("BY", "BH", "Sidak", "none"),
                        threads = NULL, chunk_size = NULL, checkpoint_dir = NULL,
                        resume = TRUE, verbose = FALSE, moments,
                        rank = .mgcvst_pca_defaults$rank,
                        n_per_cell = .mgcvst_pca_defaults$n_per_cell,
                        seed = .mgcvst_pca_defaults$seed, k = NULL,
                        screen = FALSE) {
  if (!inherits(fit, "spdeST_fit")) {
    stop("fit must be returned by spdeST.estimate().")
  }
  if (missing(moments)) .mgcvst_moments_missing()
  adjust <- match.arg(adjust)
  .mgcvst_test_run(
    fit, pairs, q.value, adjust, threads, chunk_size,
    checkpoint_dir, resume, verbose, moments = moments, rank = rank,
    n_per_cell = n_per_cell, seed = seed, k = k, screen = screen,
    call = match.call()
  )
}
