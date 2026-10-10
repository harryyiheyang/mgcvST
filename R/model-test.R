# The current mgcv model has one marked SPDE and conditional nuisance covariance.
.mgcvst_model_dense_preparation <- function(fit, features) {
  geometry <- fit$geometry
  if (identical(fit$score_backend, "sparse") || length(geometry$target) != 1L ||
      is.null(geometry$nuisance_design) ||
      length(fit$nuisance_covariance) < max(features) ||
      any(vapply(fit$nuisance_covariance[features], is.null, logical(1L)))) {
    return(NULL)
  }
  j <- unname(geometry$target[[1L]])
  s <- geometry$smooth[[j]]
  T0 <- fit$.mgcvst_fixed_factors[[j]]
  if (s$fixed || length(s$sp_index) != 1L || !is.matrix(T0)) return(NULL)
  list(T0 = T0, X = geometry$nuisance_design, sp_index = s$sp_index,
       width = stats::setNames(ncol(T0), names(geometry$target)))
}

.mgcvst_adjust_choices <- c("BY", "BH", "Sidak", "none")
.mgcvst_moments_choices <- c("exact", "pcalearning")

.mgcvst_moments_missing <- function() {
  stop("moments must be given: \"exact\" (the four exact trace moments of every ",
       "pair, k = 20 leading singular values) or \"pcalearning\" (low-rank trace ",
       "moments, k = 80). There is no default.", call. = FALSE)
}

# Shared orchestration of mgcvST.test() and inlaST.test(): validate the
# arguments, resolve the pair universe, stream the pairs of the route the user
# chose ("exact" moments or "pcalearning") to shards, adjust the two-sided
# family once and assemble the result. `entry` is "mgcv" or "inla".
.mgcvst_test_run <- function(fit, entry, pairs, q.value, adjust, threads,
                             chunk_size, checkpoint_dir, resume, verbose,
                             moments, rank = NULL, n_per_cell = NULL,
                             seed = NULL, k = NULL, call = NULL) {
  q.value <- as.numeric(q.value)
  if (length(q.value) != 1L || !is.finite(q.value) ||
      q.value <= 0 || q.value > 1) {
    stop("q.value must be one finite value in (0, 1].")
  }
  if (!is.character(adjust) || length(adjust) != 1L ||
      !(adjust %in% .mgcvst_adjust_choices)) {
    stop("adjust must be one of ",
         paste0("\"", .mgcvst_adjust_choices, "\"", collapse = ", "), ".")
  }
  if (!is.character(moments) || length(moments) != 1L ||
      !(moments %in% .mgcvst_moments_choices)) {
    stop("moments must be one of ",
         paste0("\"", .mgcvst_moments_choices, "\"", collapse = ", "), ".")
  }
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("verbose must be TRUE or FALSE.")
  }
  if (!is.logical(resume) || length(resume) != 1L || is.na(resume)) {
    stop("resume must be TRUE or FALSE.")
  }
  if (!is.null(checkpoint_dir) && (!is.character(checkpoint_dir) ||
      length(checkpoint_dir) != 1L || is.na(checkpoint_dir) ||
      !nzchar(checkpoint_dir))) {
    stop("checkpoint_dir must be NULL or one directory name.")
  }
  if (is.null(threads)) threads <- 1L
  if (!is.numeric(threads) || length(threads) != 1L || !is.finite(threads) ||
      threads < 1 || threads != floor(threads) ||
      threads > .Machine$integer.max) {
    stop("threads must be one positive integer.")
  }
  threads <- as.integer(threads)
  if (!is.null(chunk_size) && (!is.numeric(chunk_size) ||
      length(chunk_size) != 1L || !is.finite(chunk_size) || chunk_size < 1 ||
      chunk_size != floor(chunk_size))) {
    stop("chunk_size must be one positive integer.")
  }
  if (identical(entry, "inla")) {
    if (!.mgcvst_inla_downstream(fit)) {
      stop("inlaST.test() requires a fit returned by inlaST.estimate().")
    }
    .mgcvst_inla_require_sparse(fit)
  } else {
    if (.mgcvst_inla_downstream(fit)) {
      stop("mgcvST.test() does not accept inlaST.estimate() fits; use inlaST.test().")
    }
    if (!inherits(fit, "mgcvST_model_fit")) {
      stop("mgcvST.test() requires a fit returned by mgcvST.estimate().")
    }
  }
  .mgcvst_check_fit_format(fit)
  if (!is.null(fit$diagnostics$spatial_fitted) && !any(fit$diagnostics$spatial_fitted)) {
    stop("The fit has no spatial model; add spatial models with ",
         if (.mgcvst_inla_downstream(fit)) "inlaST" else "mgcvST",
         ".estimate_spatial().", call. = FALSE)
  }
  # Validated before any basis is built.
  pca <- .mgcvst_pca_check_args(rank, n_per_cell, seed, k)
  # A resumed run follows the route of its checkpoint directory.
  .mgcvst_route_check(checkpoint_dir, resume, moments)
  if (identical(entry, "mgcv") && is.null(fit$geometry)) {
    stop("fitmgcvST has no feature geometry to test.")
  }
  .mgcvst_thread_limit()

  index <- if (is.null(pairs)) NULL else .mgcvst_pair_index(pairs, fit$feature_id)
  available <- .mgcvst_feature_available(fit)
  if (is.null(index) && sum(available) < 2L) {
    stop("At least two available features are required to test all pairs.")
  }
  # Features that estimation did not select for a spatial model are not
  # failures: with pairs = NULL the test covers the features that have one.
  not_selected <- if (is.null(fit$diagnostics$spatial_selected)) {
    rep(FALSE, length(available))
  } else !available & !fit$diagnostics$spatial_selected
  extra <- NULL
  unavailable <- which(!available & !not_selected)
  if (!is.null(index)) {
    ok <- available[index[, 1L]] & available[index[, 2L]]
    if (!all(ok)) {
      extra <- .mgcvst_pairs_frame(index[!ok, 1L], index[!ok, 2L],
        status = .mgcvst_pair_status[["feature"]])
      unavailable <- sort(unique(as.vector(index[!ok, ])))
      unavailable <- unavailable[!available[unavailable]]
      index <- index[ok, , drop = FALSE]
    } else {
      unavailable <- integer()
    }
  }
  n_tested <- if (is.null(index)) {
    sum(available) * (sum(available) - 1) / 2
  } else nrow(index)
  n_used <- if (is.null(index)) sum(available) else
    length(unique(as.vector(index)))
  if (verbose && !identical(adjust, "none")) {
    guard <- .mgcvst_pair_memory_guard(
      n_tested + if (is.null(extra)) 0 else nrow(extra),
      .mgcvst_adjust_bytes_per_pair, 0.4)
    if (!guard$ok) {
      message("The adjustment of ", format(n_tested, big.mark = ","),
              " pairs needs about ", format(guard$need / 1024^3, digits = 3),
              " GiB of the ", format(guard$available / 1024^3, digits = 3),
              " GiB available and may be skipped.")
    }
  }

  basis_elapsed <- 0
  basis_rank <- NA_integer_
  basis_kind <- NA_character_
  test_started <- proc.time()[["elapsed"]]
  routed <- NULL
  route <- NULL
  if (n_tested > 0) {
    basis <- NULL
    if (identical(entry, "inla")) {
      fit <- .inlast_sparse_prepare(fit)
      t_basis <- proc.time()[["elapsed"]]
      basis <- .inlast_check_basis(fit, .inlast_sparse_observation_basis(fit))
      basis_elapsed <- proc.time()[["elapsed"]] - t_basis
      basis_rank <- basis$rank
      basis_kind <- basis$kind
    }
    q <- .mgcvst_state_width(fit, basis)
    k_route <- if (is.null(pca$k)) {
      if (identical(moments, "exact")) .mgcvst_exact_defaults$k else .mgcvst_pca_defaults$k
    } else pca$k
    route <- list(moments = moments, k = as.integer(min(k_route, q)), q = q)
    if (verbose) {
      message("Pair test: ", if (identical(moments, "exact"))
                paste0("exact moments (k = ", route$k, ")") else
                paste0("PCAlearning (rank ", pca$rank, ", k = ", route$k, ")"),
              " on q = ", q, ", ", format(n_tested, big.mark = ","), " pairs and ",
              threads, " thread", if (threads > 1L) "s", ".")
    }
    if (is.null(chunk_size)) {
      chunk_size <- if (identical(moments, "pcalearning")) 1000000L else 10000L
    }
    chunk_size <- min(chunk_size, .Machine$integer.max)
    if (identical(route$moments, "pcalearning")) {
      routed <- .mgcvst_pair_pcalearning(
        fit, index, threads, chunk_size, verbose, basis = basis,
        rank = pca$rank, n_per_cell = pca$n_per_cell, seed = pca$seed,
        k = route$k, checkpoint_dir = checkpoint_dir, resume = resume,
        route = route
      )
    } else {
      routed <- .mgcvst_pair_pipeline(
        fit, index, threads, chunk_size, verbose,
        checkpoint_dir = checkpoint_dir, resume = resume, k = route$k,
        basis = basis, route = route
      )
    }
  } else {
    root <- tempfile("mgcvst-pairs-")
    dir.create(root)
    routed <- list(pair_dir = root, shards = character(), rows = integer(),
                   n_pairs = 0, temporary = TRUE,
                   failed = data.frame(feature_id = character(),
                                       error = character()),
                   elapsed = 0, metadata = list(
                     preparation_elapsed = 0, chunks = 0L,
                     contract = .mgcvst_contract("exact")))
  }
  finalized <- .mgcvst_pairs_finalize(
    routed$pair_dir, routed$shards, routed$rows, extra, adjust, q.value,
    temporary = routed$temporary, verbose = verbose
  )

  failed <- routed$failed
  if (length(unavailable)) {
    reason <- fit$diagnostics$error_message[unavailable]
    reason <- ifelse(is.na(reason), "The feature has no usable fit.", reason)
    reason[not_selected[unavailable]] <-
      "The feature has no spatial fit: it was not selected in step 2."
    failed <- rbind(failed, data.frame(
      feature_id = fit$feature_id[unavailable], error = reason,
      stringsAsFactors = FALSE))
  }
  failed <- failed[!duplicated(failed$feature_id), , drop = FALSE]
  rownames(failed) <- NULL

  summary_elapsed <- routed$metadata$preparation_elapsed
  elapsed <- summary_elapsed + routed$elapsed
  timing <- list(
    elapsed = elapsed, summary_elapsed = summary_elapsed,
    pair_elapsed = elapsed - summary_elapsed,
    workers = threads, chunks = routed$metadata$chunks,
    backend = "C++ OpenMP", preparation_backend = "C++ OpenMP",
    preparation_threads = threads,
    test_wall_elapsed = proc.time()[["elapsed"]] - test_started
  )
  if (identical(entry, "inla") && !is.null(route)) {
    # The observation basis of a sparse INLA fit: q is the dimension of the
    # field, r the number of basis directions (the score dimension).
    timing$inla_projection <- list(
      q = ncol(fit$score_sparse$Q), r = basis_rank, basis_kind = basis_kind,
      basis = "constrained_observation_kernel_A_Qg_inverse_At",
      basis_elapsed = basis_elapsed
    )
  }
  pca_learning <- NULL
  if (!is.null(route)) timing$route <- route
  if (identical(route$moments, "pcalearning")) {
    projection <- routed$metadata
    pca_learning <- projection$pca_learning
    projection$pca_learning <- NULL
    timing$pcalearning <- projection
  } else {
    timing$pair_pipeline <- routed$metadata
  }
  ans <- structure(
    list(
      results = finalized$results,
      shards = finalized$shards,
      feature_id = fit$feature_id,
      failed = failed,
      threshold = finalized$threshold,
      discoveries = finalized$discoveries,
      adjustment = finalized$adjustment,
      pair_contract = if (is.null(pairs)) "all_available_pairs" else
        "explicit_tested_pair_universe",
      test_definition = "single_global_cross_gene_covariance_at_independence",
      timing = timing,
      calibration = "saddlepoint",
      moments = if (is.null(route)) NA_character_ else route$moments,
      contract = routed$metadata$contract,
      checkpoint_dir = if (isTRUE(routed$temporary)) NULL else
        normalizePath(checkpoint_dir, winslash = "/", mustWork = FALSE),
      call = call
    ),
    class = "mgcvST_test"
  )
  if (!is.null(pca_learning)) ans$pca_learning <- pca_learning
  ans
}

#' Test cross-feature spatial covariance
#'
#' Tests the cross-feature spatial covariance of every requested gene pair from
#' the fixed compact summaries in a fit returned by [mgcvST.estimate()]; no GAM
#' is refitted. Each pair is tested by the signed cross-gene score
#' `U = sum(s_i x_i y_i)`, where `s_i` are the singular values of the product
#' of the square roots of the two score-covariance matrices. The tail of `U` is
#' a Lugannani-Rice saddlepoint approximation computed in log space on the
#' `k` leading singular values of a basis shared by all genes, plus a
#' remainder that matches the remaining power sums of the spectrum. Sparse
#' INLA fits from [inlaST.estimate()] are tested with [inlaST.test()], which
#' takes the same arguments.
#'
#' The argument `moments` selects the route that supplies the shared basis and
#' the remainder; it has no default. With `moments = "exact"`, the shared basis
#' holds the `k = 20` leading eigenvectors of the summed, normalized score
#' covariances, every pair needs the four exact trace moments
#' `tr((H_i H_j)^s)`, `s = 1, ..., 4`, and the remainder is two moment-matched
#' nodes (one node or a Gaussian term when the moments do not allow two nodes).
#' With `moments = "pcalearning"`, the score covariances are projected onto a
#' rank-`rank` basis learned from training genes (see [inlaST.test()] for the
#' construction), the pair traces `tr(H_i H_j)` and `tr((H_i H_j)^2)` come from
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
#' `i < j` (positions in `$feature_id`), the signed `score`, the natural-log
#' two-sided, positive and negative p-values `log_p_two_sided`,
#' `log_p_positive` and `log_p_negative`, the adjusted two-sided
#' log q-value `log_q`, the integer `remainder_kind` of the calibration (0:
#' no remainder; 1: one node; 2: two nodes; 3: Gaussian term) and the integer
#' `status` (0: evaluated; 1: trace moments non-finite or non-positive; 2:
#' invalid p-value; 3: a feature of the pair has no usable score state). A
#' pair with a status other than 0 has missing log
#' p-values and is not adjusted. Rows are written as Parquet shards while the
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
#' @param fitmgcvST A fit returned by [mgcvST.estimate()].
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
#'   Both routes serve both estimators.
#' @param rank Number of PCAlearning basis matrices (default 30). Used by the
#'   PCAlearning route.
#' @param n_per_cell Training genes drawn per PCAlearning stratification cell
#'   (default 3).
#' @param seed Non-negative integer seed for PCAlearning training-gene
#'   sampling; the caller's random-number state is restored.
#' @param k Number of leading singular values of the shared basis. `NULL`
#'   uses 20 on the exact route and 80 on the PCAlearning route; a value above
#'   the score dimension `q` is reduced to `q`, and `k = q` reproduces the
#'   full-spectrum saddlepoint of each pair.
#' @return An object of class `mgcvST_test` with `results`, `shards`,
#'   `feature_id`, `failed` (features without a usable score state and the
#'   reason), `threshold`, `discoveries`, `adjustment`, `timing`
#'   (`timing$route` holds the route, `k` and the score dimension `q`),
#'   `calibration`, `moments` (the route), `contract` and `call`. A PCAlearning
#'   run additionally returns `pca_learning`.
#' @export
mgcvST.test <- function(
    fitmgcvST, pairs = NULL, q.value = 0.05,
    adjust = c("BY", "BH", "Sidak", "none"),
    threads = NULL, chunk_size = NULL, checkpoint_dir = NULL,
    resume = TRUE, verbose = FALSE, moments,
    rank = .mgcvst_pca_defaults$rank,
    n_per_cell = .mgcvst_pca_defaults$n_per_cell,
    seed = .mgcvst_pca_defaults$seed, k = NULL) {
  if (missing(moments)) .mgcvst_moments_missing()
  adjust <- match.arg(adjust)
  .mgcvst_test_run(
    fitmgcvST, "mgcv", pairs, q.value, adjust, threads, chunk_size,
    checkpoint_dir, resume, verbose, moments = moments, rank = rank,
    n_per_cell = n_per_cell, seed = seed, k = k, call = match.call()
  )
}

#' Print covariance-test diagnostics
#'
#' @param x An `mgcvST_test` object.
#' @param ... Unused.
#' @return `x`, invisibly.
#' @export
print.mgcvST_test <- function(x, ...) {
  d <- x$discoveries
  cat("mgcvST quadratic-form covariance tests\n")
  cat("  tested pair universe:", format(d$pairs_requested, big.mark = ","), "\n")
  cat("  pairs with p-value:", format(d$pairs_with_p_value, big.mark = ","), "\n")
  cat("  adjustment:", x$threshold$adjust,
      if (!isTRUE(x$adjustment$computed)) "(skipped)" else "", "\n")
  cat("  discoveries at q <=", format(x$threshold$q_value), ":",
      format(d$pairs_discovered, big.mark = ","), "(positive",
      format(d$pairs_discovered_positive, big.mark = ","), ", negative",
      format(d$pairs_discovered_negative, big.mark = ","), ")\n")
  cat("  features without a score state:", nrow(x$failed), "\n")
  cat("  result shards:", length(x$shards), "\n")
  invisible(x)
}
