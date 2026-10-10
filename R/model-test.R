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

# Shared orchestration of mgcvST.test() and inlaST.test(): validate the
# arguments, resolve the pair universe, stream the pairs of the route
# ("exact" moments for mgcv fits, "pcalearning" for sparse INLA fits) to
# shards, adjust the two-sided family once and assemble the result.
.mgcvst_test_run <- function(fit, route, pairs, q.value, adjust, threads,
                             chunk_size, checkpoint_dir, resume, verbose,
                             rank = NULL, n_per_cell = NULL, seed = NULL,
                             call = NULL) {
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
  threads <- as.integer(threads)
  if (length(threads) != 1L || is.na(threads) || threads < 1L) {
    stop("threads must be one positive integer.")
  }
  if (is.null(chunk_size)) {
    chunk_size <- if (identical(route, "pcalearning")) 1000000L else 10000L
  }
  if (!is.numeric(chunk_size) || length(chunk_size) != 1L ||
      !is.finite(chunk_size) || chunk_size < 1 || chunk_size != floor(chunk_size)) {
    stop("chunk_size must be one positive integer.")
  }
  chunk_size <- min(chunk_size, .Machine$integer.max)
  if (identical(route, "pcalearning")) {
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
    if (is.null(fit$geometry)) {
      stop("fitmgcvST has no feature geometry to test.")
    }
  }
  .mgcvst_thread_limit()

  index <- if (is.null(pairs)) NULL else .mgcvst_pair_index(pairs, fit$feature_id)
  available <- .mgcvst_feature_available(fit)
  if (is.null(index) && sum(available) < 2L) {
    stop("At least two available features are required to test all pairs.")
  }
  extra <- NULL
  unavailable <- which(!available)
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

  inla_basis_elapsed <- 0
  test_started <- proc.time()[["elapsed"]]
  routed <- NULL
  if (n_tested > 0) {
    if (identical(route, "pcalearning")) {
      fit <- .inlast_sparse_prepare(fit)
      t_basis <- proc.time()[["elapsed"]]
      basis <- .inlast_sparse_observation_basis(fit)
      inla_basis_elapsed <- proc.time()[["elapsed"]] - t_basis
      routed <- .mgcvst_inla_test_pairs(
        fit, index, threads, chunk_size, verbose, basis = basis, rank = rank,
        n_per_cell = n_per_cell, seed = seed, checkpoint_dir = checkpoint_dir,
        resume = resume
      )
    } else {
      routed <- .mgcvst_pair_pipeline(
        fit, index, threads, chunk_size, verbose,
        checkpoint_dir = checkpoint_dir, resume = resume
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
                     contract = .mgcvst_contract(
                       if (identical(route, "pcalearning")) "pcalearning" else "exact")))
  }
  finalized <- .mgcvst_pairs_finalize(
    routed$pair_dir, routed$shards, routed$rows, extra, adjust, q.value,
    temporary = routed$temporary, verbose = verbose
  )

  failed <- routed$failed
  if (length(unavailable)) {
    message_fit <- fit$diagnostics$error_message[unavailable]
    failed <- rbind(failed, data.frame(
      feature_id = fit$feature_id[unavailable],
      error = ifelse(is.na(message_fit), "The feature has no usable fit.",
                     message_fit), stringsAsFactors = FALSE))
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
    preparation_threads = threads
  )
  pca_learning <- NULL
  if (identical(route, "pcalearning")) {
    projection <- routed$metadata
    pca_learning <- projection$pca_learning
    projection$pca_learning <- NULL
    projection$basis_elapsed <- inla_basis_elapsed
    projection$test_wall_elapsed <- proc.time()[["elapsed"]] - test_started
    timing$inla_projection <- projection
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
      calibration = "liu",
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
#' is refitted. Each pair is tested by the squared signed cross-gene score,
#' calibrated by Liu moment matching of exact trace moments of the two
#' score-covariance matrices. Sparse INLA fits from [inlaST.estimate()] are
#' tested with [inlaST.test()], which takes the same arguments.
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
#' Liu moment matching without remainder) and the integer `status` (0:
#' evaluated; 1: trace moments non-finite or non-positive; 2: invalid
#' p-value; 3: a feature of the pair has no usable score state). Rows are
#' written as Parquet shards while the pairs are evaluated; `$shards` lists the
#' files, and `$results` is the same table sorted by `(i, j)` when it fits the
#' memory guard (56 bytes per pair, 20% of available memory), and `NULL`
#' otherwise.
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
#'   written per shard. `NULL` uses 10,000 for [mgcvST.test()] and 1,000,000
#'   for [inlaST.test()].
#' @param checkpoint_dir Optional checkpoint directory. Reusable feature score
#'   states and the raw pair shards are saved there, and a repeated call
#'   resumes completed shards. Pair results are keyed by the algorithm
#'   contract; a directory holding pair results from another contract is
#'   refused. With `NULL`, temporary storage is used.
#' @param resume Reuse compatible completed checkpoint entries.
#' @param verbose Logical; report progress.
#' @return An object of class `mgcvST_test` with `results`, `shards`,
#'   `feature_id`, `failed` (features without a usable score state and the
#'   reason), `threshold`, `discoveries`, `adjustment`, `timing`,
#'   `calibration`, `contract` and `call`. [inlaST.test()] additionally
#'   returns `pca_learning`.
#' @export
mgcvST.test <- function(
    fitmgcvST, pairs = NULL, q.value = 0.05,
    adjust = c("BY", "BH", "Sidak", "none"),
    threads = NULL, chunk_size = NULL, checkpoint_dir = NULL,
    resume = TRUE, verbose = FALSE) {
  adjust <- match.arg(adjust)
  .mgcvst_test_run(
    fitmgcvST, "exact", pairs, q.value, adjust, threads, chunk_size,
    checkpoint_dir, resume, verbose, call = match.call()
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
