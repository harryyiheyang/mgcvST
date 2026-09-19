# Compare R feature orchestration with the OpenMP score preparation kernel.
library(mgcvST)

out <- "artifacts/package-maintenance-20260919"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
F <- readRDS("artifacts/inla-bam-validation/components/unadjusted/bam-compact-fit.rds")
mgcvST:::.mgcvst_thread_limit()
F$.mgcvst_fixed_factors <- mgcvST:::.mgcvst_model_fixed_factors(F)
ids <- seq_len(32L)
stopifnot(all(mgcvST:::.mgcvst_feature_available(F)[ids]))
native <- mgcvST:::.mgcvst_model_dense_preparation(F, ids)
if (is.null(native)) stop("The benchmark fit lacks conditional mgcv score geometry.")
R <- list()
ref <- list()
for (threads in c(0L, 1L, 2L, 4L)) {
  for (rr in seq_len(3L)) {
    t0 <- proc.time()[["elapsed"]]
    if (threads == 0L) {
      ans <- lapply(ids, function(i) {
        z <- mgcvST:::.mgcvst_model_score_state(F, i)
        list(a = z$a, H = unname(z$M), error = NULL)
      })
    } else {
      ans <- mgcvST:::mgcvst_dense_score_batch_cpp(
        native$T0, F$working_variance[, ids, drop = FALSE],
        F$working_error[, ids, drop = FALSE],
        F$dispersion[ids] / F$smoothing_parameters[ids, native$sp_index],
        native$X, F$nuisance_covariance[ids], threads)
    }
    elapsed <- proc.time()[["elapsed"]] - t0
    if (threads == 0L && rr == 1L) ref <- ans
    if (!isTRUE(all.equal(ans, ref, tolerance = 1e-10))) {
      stop("The batched score states differ from the R reference.")
    }
    da <- max(vapply(seq_along(ids), function(i) max(abs(ans[[i]]$a - ref[[i]]$a)), numeric(1L)))
    dH <- max(vapply(seq_along(ids), function(i) max(abs(ans[[i]]$H - ref[[i]]$H)), numeric(1L)))
    R[[length(R) + 1L]] <- data.frame(
      engine = if (threads == 0L) "R_feature_loop" else "C++_feature_loop",
      threads = max(1L, threads), repeat_id = rr,
      observations = nrow(F$working_error), features = length(ids),
      score_dimension = ncol(native$T0), seconds = elapsed,
      max_a_difference = da, max_H_difference = dH)
    write.csv(do.call(rbind, R), file.path(out, "dense-preparation.csv"), row.names = FALSE)
  }
}
print(do.call(rbind, R))
writeLines(capture.output(sessionInfo()), file.path(out, "dense-preparation-session.txt"))
