# Resolve the score engine recorded by the fitted object.
.mgcvst_test_engine <- function(fit) {
  engine <- fit$test_engine
  if (is.null(engine)) {
    stop("fitmgcvST does not record a score-test engine.")
  }
  engine <- as.character(engine)
  registry <- c(
    spde = ".mgcvst_test_spde",
    single_model = ".mgcvst_test_model_single"
  )
  target <- unname(registry[engine])
  if (length(target) != 1L || is.na(target)) {
    stop("fitmgcvST contains an unknown score-test engine: ", engine)
  }
  get(target, envir = environment(.mgcvst_test_engine), inherits = TRUE)
}

# Evaluate a chunk with one already-selected model score function.
.mgcvst_model_test_chunk <- function(payload, fit, pair_function,
                                      calibration) {
  .mgcvst_thread_limit()
  out <- vector("list", length(payload$rows))
  for (k in seq_along(payload$rows)) {
    i <- payload$pairs[k, 1L]
    j <- payload$pairs[k, 2L]
    z <- tryCatch(
      pair_function(fit, i, j, calibration = calibration),
      error = function(e) e
    )
    if (inherits(z, "condition")) {
      out[[k]] <- data.frame(
        pair_index = payload$rows[k],
        signed_score = NA_real_, information = NA_real_,
        effective_rank = NA_real_, p_two_sided = NA_real_,
        p_positive = NA_real_, p_negative = NA_real_,
        error_message = conditionMessage(z),
        stringsAsFactors = FALSE
      )
      next
    }
    out[[k]] <- data.frame(
      pair_index = payload$rows[k], signed_score = z$score,
      information = z$information,
      effective_rank = z$effective_rank,
      p_two_sided = z$p_two_sided, p_positive = z$p_positive,
      p_negative = z$p_negative,
      error_message = NA_character_,
      stringsAsFactors = FALSE
    )
  }
  out
}

# Shared orchestration for model.set() score engines.
.mgcvst_test_model <- function(
    fitmgcvST, pair_function, test_definition,
    q.value = 0.05, FDR = TRUE, method = "BH",
    BPPARAM = BiocParallel::SerialParam(), ...,
    pairs = NULL, highlight = NULL,
    calibration = c("liu", "davies"),
    chunk_size = NULL,
    threads = NULL, verbose = FALSE) {
  if (!inherits(fitmgcvST, "mgcvST_model_fit")) {
    stop("The model score engine requires a fit from mgcvST.estimate(Y, model).")
  }
  if (is.null(fitmgcvST$geometry)) {
    stop("fitmgcvST has no feature geometry to test.")
  }
  q.value <- as.numeric(q.value)
  if (length(q.value) != 1L || !is.finite(q.value) ||
      q.value <= 0 || q.value > 1) {
    stop("q.value must be one finite value in (0, 1].")
  }
  if (!is.logical(FDR) || length(FDR) != 1L || is.na(FDR)) {
    stop("FDR must be TRUE or FALSE.")
  }
  if (!is.character(method) || length(method) != 1L || is.na(method) ||
      !(method %in% stats::p.adjust.methods)) {
    stop("method must be one of stats::p.adjust.methods.")
  }
  if (!inherits(BPPARAM, "BiocParallelParam")) {
    stop("BPPARAM must inherit from 'BiocParallelParam'.")
  }
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("verbose must be TRUE or FALSE.")
  }
  unused <- list(...)
  if (length(unused)) {
    stop("Unused arguments in ...: ", paste(names(unused), collapse = ", "))
  }
  calibration <- match.arg(calibration)
  if (calibration == "davies" &&
      !requireNamespace("CompQuadForm", quietly = TRUE)) {
    stop("calibration = 'davies' requires the optional CompQuadForm package.")
  }
  if (is.null(threads)) threads <- BiocParallel::bpworkers(BPPARAM)
  threads <- as.integer(threads)
  if (length(threads) != 1L || is.na(threads) || threads < 1L) {
    stop("threads must be one positive integer.")
  }
  .mgcvst_thread_limit()

  index <- .mgcvst_pair_index(pairs, fitmgcvST$feature_id)
  highlight_index <- matrix(integer(), nrow = 0L, ncol = 2L)
  if (!is.null(highlight)) {
    highlight_index <- .mgcvst_pair_index(highlight, fitmgcvST$feature_id)
  }
  n_feature <- length(fitmgcvST$feature_id)
  key <- (index[, 1L] - 1L) * n_feature + index[, 2L]
  highlight_key <- if (nrow(highlight_index)) {
    (highlight_index[, 1L] - 1L) * n_feature + highlight_index[, 2L]
  } else {
    numeric()
  }
  extra <- which(!(highlight_key %in% key))
  if (length(extra)) {
    index <- rbind(index, highlight_index[extra, , drop = FALSE])
    key <- c(key, highlight_key[extra])
  }
  highlighted <- key %in% highlight_key

  available <- .mgcvst_feature_available(fitmgcvST)
  i1 <- index[, 1L]
  i2 <- index[, 2L]
  pair_available <- available[i1] & available[i2]
  result <- data.frame(
    pair_index = seq_len(nrow(index)),
    feature1 = fitmgcvST$feature_id[i1],
    feature2 = fitmgcvST$feature_id[i2],
    signed_score = NA_real_,
    statistic = NA_real_, information = NA_real_,
    effective_rank = NA_real_, p_two_sided = NA_real_,
    p_positive = NA_real_, p_negative = NA_real_,
    p_adjusted = NA_real_, p_positive_adjusted = NA_real_,
    p_negative_adjusted = NA_real_, discovered = FALSE,
    discovered_positive = FALSE, discovered_negative = FALSE,
    highlighted = highlighted, retained = highlighted,
    error_message = NA_character_,
    stringsAsFactors = FALSE
  )
  unavailable_rows <- which(!pair_available)
  for (k in unavailable_rows) {
    missing_feature <- c(i1[k], i2[k])[!available[c(i1[k], i2[k])]]
    result$error_message[k] <- paste(
      paste0(
        fitmgcvST$feature_id[missing_feature], ": ",
        fitmgcvST$diagnostics$error_message[missing_feature]
      ),
      collapse = " | "
    )
  }

  tested_rows <- which(pair_available)
  workers <- if (length(tested_rows)) {
    max(1L, min(length(tested_rows), BiocParallel::bpworkers(BPPARAM)))
  } else {
    0L
  }
  if (is.null(chunk_size)) {
    chunk_size <- if (workers > 0L) ceiling(length(tested_rows) / workers) else 1L
  }
  chunk_size <- as.integer(chunk_size)
  if (length(chunk_size) != 1L || is.na(chunk_size) || chunk_size < 1L) {
    stop("chunk_size must be one positive integer.")
  }
  chunks <- list()
  elapsed <- 0
  if (length(tested_rows)) {
    chunks <- split(tested_rows, ceiling(seq_along(tested_rows) / chunk_size))
    payload <- lapply(chunks, function(rows) list(
      rows = rows, pairs = index[rows, , drop = FALSE]
    ))
    worker_bundle <- .mgcvst_worker_bundle()
    test_chunk <- get(".mgcvst_model_test_chunk", envir = worker_bundle,
                      inherits = FALSE)
    pair_worker <- get(deparse(substitute(pair_function)), envir = worker_bundle,
                       inherits = FALSE)
    t0 <- proc.time()[["elapsed"]]
    evaluated <- BiocParallel::bplapply(
      payload, test_chunk, fit = fitmgcvST, pair_function = pair_worker,
      calibration = calibration, BPPARAM = BPPARAM
    )
    elapsed <- proc.time()[["elapsed"]] - t0
    evaluated <- unlist(evaluated, recursive = FALSE)
    for (z in evaluated) {
      target <- z$pair_index
      names <- intersect(names(z), names(result))
      result[target, names] <- z[1L, names, drop = FALSE]
      result$statistic[target] <- result$signed_score[target]^2
    }
  }

  valid <- is.finite(result$p_two_sided) &
    result$p_two_sided >= 0 & result$p_two_sided <= 1
  if (FDR) {
    result$p_adjusted[valid] <- stats::p.adjust(result$p_two_sided[valid], method)
    result$p_positive_adjusted[valid] <- stats::p.adjust(
      result$p_positive[valid], method
    )
    result$p_negative_adjusted[valid] <- stats::p.adjust(
      result$p_negative[valid], method
    )
  } else {
    result$p_adjusted[valid] <- result$p_two_sided[valid]
    result$p_positive_adjusted[valid] <- result$p_positive[valid]
    result$p_negative_adjusted[valid] <- result$p_negative[valid]
  }
  result$discovered <- valid & result$p_adjusted <= q.value
  result$discovered_positive <- valid & result$p_positive_adjusted <= q.value
  result$discovered_negative <- valid & result$p_negative_adjusted <= q.value
  result$retained <- result$highlighted | result$discovered
  raw_threshold <- if (any(result$discovered)) {
    max(result$p_two_sided[result$discovered])
  } else {
    NA_real_
  }
  structure(
    list(
      results = result,
      threshold = list(
        q_value = q.value, FDR = FDR,
        adjustment_method = if (FDR) method else "none",
        raw_p_threshold = raw_threshold
      ),
      discoveries = list(
        pairs_requested = nrow(index),
        pairs_tested = length(tested_rows), pairs_with_p_value = sum(valid),
        pairs_discovered = sum(result$discovered),
        pairs_discovered_positive = sum(result$discovered_positive),
        pairs_discovered_negative = sum(result$discovered_negative),
        pairs_highlighted = sum(result$highlighted),
        pairs_retained = sum(result$retained)
      ),
      pair_contract = paste0(
        "explicit_tested_pair_universe_with_FDR_discoveries_",
        "union_force_retained_highlights"
      ),
      test_definition = test_definition,
      timing = list(
        elapsed = elapsed, summary_elapsed = 0, pair_elapsed = elapsed,
        workers = workers, chunks = length(chunks),
        backend = class(BPPARAM)[1L]
      ),
      calibration = calibration,
      call = match.call()
    ),
    class = "mgcvST_test"
  )
}

# Single marked-SPDE entry point.
.mgcvst_test_model_single <- function(...) {
  .mgcvst_test_model(
    ..., pair_function = .mgcvst_model_pair_single,
    test_definition = "single_global_cross_gene_covariance_at_independence"
  )
}

#' Test cross-feature spatial covariance
#'
#' Dispatches a compact fit to its registered score engine. Standard one-SPDE
#' fits use the SPDE score path. One-component fits constructed from
#' [model.set()] use the model score path.
#'
#' @inheritParams .mgcvst_test_spde
#' @export
mgcvST.test <- function(
    fitmgcvST, q.value = 0.05, FDR = TRUE, method = "BH",
    BPPARAM = BiocParallel::SerialParam(), ...,
    pairs = NULL, highlight = NULL,
    calibration = c("liu", "davies"),
    chunk_size = NULL,
    threads = NULL, verbose = FALSE) {
  engine <- .mgcvst_test_engine(fitmgcvST)
  engine(
    fitmgcvST = fitmgcvST, q.value = q.value, FDR = FDR, method = method,
    BPPARAM = BPPARAM, ..., pairs = pairs, highlight = highlight,
    calibration = calibration, chunk_size = chunk_size,
    threads = threads, verbose = verbose
  )
}
