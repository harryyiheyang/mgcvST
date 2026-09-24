# Benchmark: fused C++ Liu pair kernel (mgcvst_pair_liu_cpp) vs the old
# mgcvst_pair_trace_powers_cpp + R Liu path, and .mgcvst_pair_pipeline
# end-to-end (origin/main installed to a temporary library) vs the new
# pipeline in this checkout. Not run under testthat; run interactively with
#   Rscript tests/benchmark/bench-liu-exact.R
# Requires an OLD_MGCVST_LIB environment variable pointing at a library
# directory that already has the pre-Liu-exact mgcvST installed (see
# tests/benchmark/README below), otherwise the pipeline comparison is skipped.

suppressPackageStartupMessages({
  library(mgcvST)
})

cat("== Part 1: fused Liu pair kernel vs old trace-powers + R Liu ==\n")

bench_kernel <- function(q, K, threads = 4L, seed = 1L) {
  set.seed(seed)
  H <- lapply(seq_len(K), function(k) {
    z <- matrix(rnorm(q * q), q, q)
    crossprod(z) + diag(q) * 0.1
  })
  a <- matrix(rnorm(q * K), q, K)
  idx <- which(upper.tri(matrix(0, K, K)), arr.ind = TRUE)
  left <- idx[, 2L]
  right <- idx[, 1L]
  ord <- order(left)
  left <- left[ord]
  right <- right[ord]
  pairs <- cbind(left, right)

  t_old <- system.time({
    old_score <- colSums(a[, left, drop = FALSE] * a[, right, drop = FALSE])
    old_moments <- mgcvST:::mgcvst_pair_trace_powers_cpp(
      H, pairs, maxPower = 4L, threads = threads
    )
    old_liu <- mgcvST:::.liu_squared_score_moments(
      abs(old_score), old_moments[, 1L], old_moments[, 2L],
      old_moments[, 3L], old_moments[, 4L]
    )
  })["elapsed"]

  t_new <- system.time({
    new <- mgcvST:::mgcvst_pair_liu_cpp(H, a, left, right, threads = threads)
  })["elapsed"]

  data.frame(q = q, K = K, pairs = nrow(pairs),
             old_seconds = as.numeric(t_old), new_seconds = as.numeric(t_new),
             speedup = as.numeric(t_old) / as.numeric(t_new))
}

kernel_table <- do.call(rbind, lapply(c(30L, 60L, 100L), function(q) {
  bench_kernel(q, K = 400L, threads = 4L)
}))
print(kernel_table, row.names = FALSE)

cat("\n== Part 2: .mgcvst_pair_pipeline end-to-end, old vs new (~2e5 pairs) ==\n")

old_lib <- Sys.getenv("OLD_MGCVST_LIB", unset = NA)
if (is.na(old_lib) || !dir.exists(file.path(old_lib, "mgcvST"))) {
  cat("OLD_MGCVST_LIB is not set to a library containing an origin/main-",
      "vintage mgcvST install; skipping the end-to-end pipeline comparison.\n",
      "To run it: R CMD INSTALL --library=<dir> <origin/main checkout>, then\n",
      "OLD_MGCVST_LIB=<dir> Rscript tests/benchmark/bench-liu-exact.R\n", sep = "")
} else {
  # A hand-built model fit large enough to carry ~2e5 pairs without needing a
  # real mgcv/GAM estimation pass (mirrors the fixture used in
  # test-wgcna.R's "WGCNA retains legacy compact B-Q-X score semantics").
  make_fit <- function(p, n = 30L, q = 6L, seed = 909L) {
    set.seed(seed)
    B <- matrix(rnorm(n * q), n, q)
    Q <- diag(runif(q, 0.5, 2))
    X <- matrix(1, n, 1L)
    ids <- paste0("g", seq_len(p))
    fit <- list(
      feature_id = ids,
      working_error = matrix(rnorm(n * p), n, p),
      working_variance = matrix(runif(n * p, 0.8, 1.3), n, p),
      dispersion = seq(0.9, 1.2, length.out = p),
      lambda = seq(1, 1.6, length.out = p),
      geometry = list(B = B, Q = Q, X = X),
      test_engine = "spde"
    )
    colnames(fit$working_error) <- colnames(fit$working_variance) <- ids
    fit
  }
  p <- 632L  # choose(632, 2) ~= 1.995e5 pairs
  fit <- make_fit(p)
  pairs <- t(utils::combn(p, 2L))
  cat("pairs:", nrow(pairs), "\n")

  new_env <- new.env()
  new_env$fit <- fit
  new_env$pairs <- pairs
  t_new <- system.time({
    new_result <- mgcvST:::.mgcvst_pair_pipeline(
      new_env$fit, new_env$pairs, seq_len(nrow(new_env$pairs)),
      threads = 4L, chunk_size = 20000L, verbose = FALSE
    )
  })["elapsed"]

  old_script <- tempfile(fileext = ".R")
  writeLines(c(
    sprintf('.libPaths(c(%s, .libPaths()))', shQuote(old_lib)),
    'library(mgcvST, lib.loc = commandArgs(trailingOnly = TRUE)[1])',
    'fit <- readRDS(commandArgs(trailingOnly = TRUE)[2])',
    'pairs <- readRDS(commandArgs(trailingOnly = TRUE)[3])',
    't_old <- system.time({',
    '  old_result <- mgcvST:::.mgcvst_pair_pipeline(',
    '    fit, pairs, seq_len(nrow(pairs)), threads = 4L, chunk_size = 20000L,',
    '    verbose = FALSE',
    '  )',
    '})["elapsed"]',
    'saveRDS(list(elapsed = as.numeric(t_old), result = old_result$result),',
    '        commandArgs(trailingOnly = TRUE)[4])'
  ), old_script)
  fit_path <- tempfile(fileext = ".rds")
  pairs_path <- tempfile(fileext = ".rds")
  out_path <- tempfile(fileext = ".rds")
  saveRDS(fit, fit_path)
  saveRDS(pairs, pairs_path)
  status <- system2(
    file.path(R.home("bin"), "Rscript"),
    c(shQuote(old_script), shQuote(old_lib), shQuote(fit_path),
      shQuote(pairs_path), shQuote(out_path))
  )
  if (status != 0 || !file.exists(out_path)) {
    cat("The old-pipeline subprocess failed; skipping the comparison.\n")
  } else {
    old <- readRDS(out_path)
    pipeline_table <- data.frame(
      pairs = nrow(pairs), old_seconds = old$elapsed,
      new_seconds = as.numeric(t_new),
      speedup = old$elapsed / as.numeric(t_new)
    )
    print(pipeline_table, row.names = FALSE)
    cat("max |score| difference:",
        max(abs(old$result$score - new_result$result$score), na.rm = TRUE), "\n")
  }
}
