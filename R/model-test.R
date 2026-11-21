.mgcvst_adjust_choices <- c("BY", "BH", "Sidak", "none")
.mgcvst_moments_choices <- c("exact", "pcalearning")

.mgcvst_moments_missing <- function() {
  stop("moments must be given: \"exact\" (the four exact trace moments and ",
       "k = 20 leading singular values) or \"pcalearning\" (low-rank trace ",
       "moments, k = 80). ",
       "There is no default.", call. = FALSE)
}

# The pair test: validate the arguments, resolve the pair universe, stream the
# pairs of the route the user chose ("exact" moments or "pcalearning") to
# shards, adjust the two-sided family once and assemble the result.
.mgcvst_test_run <- function(fit, pairs, q.value, adjust, threads,
                             chunk_size, checkpoint_dir, resume, verbose,
                             moments, rank = NULL, n_per_cell = NULL,
                             seed = NULL, k = NULL, screen = FALSE,
                             call = NULL) {
  if (!is.logical(screen) || length(screen) != 1L || is.na(screen)) {
    stop("screen must be TRUE or FALSE.")
  }
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
  if (!inherits(fit, "inlaST_fit")) {
    stop("The test requires a fit returned by spdeST.estimate().")
  }
  .mgcvst_check_fit_format(fit)
  if (!is.null(fit$diagnostics$spatial_fitted) && !any(fit$diagnostics$spatial_fitted)) {
    stop("The fit has no spatial model; add spatial models with spdeST.estimate().",
         call. = FALSE)
  }
  # Validated before any basis is built.
  pca <- .mgcvst_pca_check_args(rank, n_per_cell, seed, k)
  # A resumed run follows the route of its checkpoint directory.
  .mgcvst_route_check(checkpoint_dir, resume, moments)
  old_threads <- .mgcvst_thread_limit()
  on.exit(.mgcvst_thread_restore(old_threads), add = TRUE)

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
  screen_record <- NULL
  if (n_tested > 0) {
    fit <- .inlast_sparse_prepare(fit)
    t_basis <- proc.time()[["elapsed"]]
    basis <- .inlast_check_basis(fit, .inlast_sparse_observation_basis(fit))
    basis_elapsed <- proc.time()[["elapsed"]] - t_basis
    basis_rank <- basis$rank
    basis_kind <- basis$kind
    q <- as.integer(basis$rank)
    k_route <- if (!is.null(pca$k)) pca$k else if (identical(moments, "pcalearning")) {
      .mgcvst_pca_defaults$k
    } else .mgcvst_exact_defaults$k
    route <- list(moments = moments, k = as.integer(min(k_route, q)), q = q)
    screening <- if (screen) .mgcvst_screen_z(adjust, q.value, n_tested) else
      list(z = 0, theta = NA_real_, tau = NA_real_)
    if (screen && verbose) {
      message(if (screening$z > 0) sprintf(
        "Screen: pairs with |U| / sqrt(t_1) < %.3f (p > %.3g) are not evaluated.",
        screening$z, screening$tau) else paste0(
        "Screen off: the rejection threshold of adjust = \"", adjust,
        "\" at this q.value is too large for the Gaussian bound."))
    }
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
        route = route, screen_z = screening$z
      )
    } else {
      routed <- .mgcvst_pair_pipeline(
        fit, index, threads, chunk_size, verbose,
        checkpoint_dir = checkpoint_dir, resume = resume, k = route$k,
        basis = basis, route = route, screen_z = screening$z
      )
    }
    if (screen) screen_record <- c(screening, list(
      screened = routed$metadata$screened))
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
    backend = "C++ OpenMP", simd = mgcvst_simd_cpp(),
    preparation_backend = "C++ OpenMP",
    preparation_threads = threads,
    test_wall_elapsed = proc.time()[["elapsed"]] - test_started
  )
  if (!is.null(route)) {
    # The observation basis of the fit: q is the dimension of the
    # field, r the number of basis directions (the score dimension).
    timing$inla_projection <- list(
      q = ncol(fit$score_sparse$Q), r = basis_rank, basis_kind = basis_kind,
      basis = "constrained_observation_kernel_A_Qg_inverse_At",
      basis_elapsed = basis_elapsed
    )
  }
  pca_learning <- NULL
  if (!is.null(route)) timing$route <- route
  if (!is.null(screen_record)) timing$screen <- screen_record
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
  # Per-feature traces of the score curvature M_g of every feature in a tested
  # pair: tr(M_g) (the expected a_g' a_g) and tr(M_g^2). With the per-pair t1
  # they are the inputs of the Haseman-Elston-type estimates U / t1.
  if (!is.null(routed$features)) {
    ans$features <- data.frame(
      feature_id = fit$feature_id[routed$features$feature_index],
      routed$features, stringsAsFactors = FALSE)
  }
  ans
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
