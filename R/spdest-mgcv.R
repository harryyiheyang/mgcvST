# spdeST implementation by mgcv: the methods behind spdeST.scan(),
# spdeST.estimate() and spdeST.test() for a model from mgcvST.set().

.spdest_mgcv_settings <- function(model) {
  defaults <- list(control = mgcv::gam.control(nthreads = 1L), bam_args = list(),
                   retain_smooth = FALSE, diagnostics = FALSE)
  utils::modifyList(defaults, if (is.null(model$estimation)) list() else model$estimation)
}

.spdest_scan.mgcvST_model <- function(model, Y, offset, feature_id, adjust,
                                      q.value, BPPARAM, chunk_size, threads,
                                      checkpoint_dir, resume, scan_control) {
  s <- .spdest_mgcv_settings(model)
  args <- c(
    list(Y = Y, G = model, feature_id = feature_id, offset = offset,
         control = s$control, retain_smooth = s$retain_smooth,
         diagnostics = s$diagnostics, adjust = adjust,
         q.value = q.value, resume = resume, scan_control = scan_control),
    .spdest_supplied(BPPARAM = BPPARAM, chunk_size = chunk_size,
                     checkpoint_dir = checkpoint_dir),
    s$bam_args
  )
  do.call(mgcvST.estimate, args)
}

.spdest_estimate.mgcvST_model_fit <- function(scan, Y, features, adjust, q.value,
                                              BPPARAM, chunk_size, threads,
                                              checkpoint_dir, resume) {
  args <- c(
    list(scan, Y = Y, features = features, resume = resume),
    .spdest_supplied(adjust = adjust, q.value = q.value, BPPARAM = BPPARAM,
                     chunk_size = chunk_size, checkpoint_dir = checkpoint_dir)
  )
  do.call(mgcvST.estimate_spatial, args)
}

.spdest_test.mgcvST_model_fit <- function(fit, ...) mgcvST.test(fit, ...)
