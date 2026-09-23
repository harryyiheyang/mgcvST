library(mgcvST)
mgcvST:::.mgcvst_thread_limit()
out <- Sys.getenv('MGCVST_CUR_BENCHMARK', 'artifacts/pair-pipeline-validation/landmark-cur-635')
dir.create(out, recursive = TRUE, showWarnings = FALSE)
paths <- c(INLA = 'artifacts/inla-bam-validation/components/unadjusted/inla-fit.rds',
           BAM = 'artifacts/inla-bam-validation/components/unadjusted/bam-compact-fit.rds')
errors <- timings <- list()
for (engine in names(paths)) {
  fit <- readRDS(paths[[engine]])
  P <- t(utils::combn(seq_along(fit$feature_id), 2L))
  path <- file.path(out, paste0(engine, '-exact'))
  if (dir.exists(path)) stop('Benchmark checkpoint already exists: ', path)
  tm <- system.time(exact <- mgcvST:::.mgcvst_pair_pipeline(
    fit, P, seq_len(nrow(P)), threads = 20L, chunk_size = 10000L,
    verbose = TRUE, checkpoint_dir = path))
  by <- p.adjust(exact$result$p_value, method = 'BY')
  timings[[length(timings) + 1L]] <- data.frame(
    engine = engine, method = 'exact', n_ref = 0L, genes = length(fit$feature_id),
    pairs = nrow(P), total_seconds = tm[['elapsed']], builds = exact$metadata$builds,
    tail_pairs = 0L, BY_rejections = sum(by <= .05, na.rm = TRUE),
    BY_disagreements = 0L, invalid_p = sum(!is.finite(exact$result$p_value)))
  saveRDS(exact, file.path(out, paste0(engine, '-exact.rds')))
  for (method in c('random', 'score', 'hyper')) {
    for (R in c(50L, 100L, 200L)) {
      path <- file.path(out, paste(engine, method, R, sep = '-'))
      if (dir.exists(path)) stop('Benchmark checkpoint already exists: ', path)
      tm <- system.time(ans <- mgcvST:::.mgcvst_pair_approximate(
        fit, P, seq_len(nrow(P)), threads = 20L, chunk_size = 10000L,
        verbose = TRUE, n_ref = R, ref_method = method, ref_seed = 20260923L,
        checkpoint_dir = path, tail_recheck = 1e-3, diagnostic_pairs = 10000L))
      tab <- ans$metadata$diagnostics$error
      tab$engine <- engine
      tab$method <- method
      tab$n_ref <- R
      errors[[length(errors) + 1L]] <- tab
      adjusted <- p.adjust(ans$result$p_value, method = 'BY')
      timings[[length(timings) + 1L]] <- data.frame(
        engine = engine, method = method, n_ref = R, genes = length(fit$feature_id),
        pairs = nrow(P), total_seconds = tm[['elapsed']], builds = ans$metadata$builds,
        tail_pairs = ans$metadata$tail_pairs,
        BY_rejections = sum(adjusted <= .05, na.rm = TRUE),
        BY_disagreements = sum((adjusted <= .05) != (by <= .05), na.rm = TRUE),
        invalid_p = sum(!is.finite(ans$result$p_value)))
      saveRDS(ans, paste0(path, '.rds'))
      write.csv(do.call(rbind, errors), file.path(out, 'holdout-errors.csv'), row.names = FALSE)
      write.csv(do.call(rbind, timings), file.path(out, 'timings.csv'), row.names = FALSE)
      rm(ans)
    }
  }
  rm(fit, exact)
  gc()
}
writeLines(capture.output(sessionInfo()), file.path(out, 'session-info.txt'))
