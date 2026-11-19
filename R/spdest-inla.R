# spdeST implementation by INLA: the methods behind spdeST.scan(),
# spdeST.estimate() and spdeST.test() for a model from inlaST.set().

.spdest_inla_settings <- function(model) {
  defaults <- list(retain_smooth = FALSE, diagnostics = FALSE)
  utils::modifyList(defaults, if (is.null(model$estimation)) list() else model$estimation)
}

.spdest_scan.inlaST_model <- function(model, Y, offset, feature_id, adjust,
                                      q.value, BPPARAM, chunk_size, threads,
                                      checkpoint_dir, resume) {
  s <- .spdest_inla_settings(model)
  args <- c(
    list(Y = Y, model = model, feature_id = feature_id, offset = offset,
         retain_smooth = s$retain_smooth, diagnostics = s$diagnostics,
         threads = threads, adjust = adjust,
         q.value = q.value, resume = resume),
    .spdest_supplied(BPPARAM = BPPARAM, chunk_size = chunk_size,
                     checkpoint_dir = checkpoint_dir)
  )
  do.call(inlaST.estimate, args)
}

.spdest_estimate.inlaST_fit <- function(scan, Y, features, adjust, q.value,
                                        BPPARAM, chunk_size, threads,
                                        checkpoint_dir, resume) {
  args <- c(
    list(scan, Y = Y, features = features, threads = threads, resume = resume),
    .spdest_supplied(adjust = adjust, q.value = q.value, BPPARAM = BPPARAM,
                     chunk_size = chunk_size, checkpoint_dir = checkpoint_dir)
  )
  do.call(inlaST.estimate_spatial, args)
}

.spdest_test.inlaST_fit <- function(fit, ...) inlaST.test(fit, ...)
