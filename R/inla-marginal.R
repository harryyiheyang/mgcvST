.inlast_marginal <- function(fit, features = NULL, calibration = "liu",
                             BPPARAM = BiocParallel::SerialParam(),
                             chunk_size = 16L, threads = 1L) {
  if (!identical(calibration, "liu")) stop("INLA supports only calibration = 'liu'.")
  if (!inherits(BPPARAM, "SerialParam")) {
    stop("INLA downstream calculations use OpenMP threads, not BiocParallel workers.")
  }
  if (is.null(fit$marginal_data)) stop("Estimate with retain_marginal = TRUE first.")
  if (!identical(fit$score_backend, "sparse")) {
    stop("INLA Liu requires a supported sparse score geometry.")
  }
  for (name in c("chunk_size", "threads")) {
    x <- get(name)
    if (length(x) != 1L || !is.numeric(x) || !is.finite(x) ||
        x < 1 || x > .Machine$integer.max || x != as.integer(x)) {
      stop(name, " must be a positive integer.")
    }
  }
  ids <- fit$feature_id
  if (is.null(features)) features <- seq_along(ids)
  i <- if (is.character(features)) match(features, ids) else features
  if (!is.numeric(i) || !length(i) || anyNA(i) || any(!is.finite(i)) ||
      any(i != as.integer(i)) || any(i < 1 | i > length(ids)) || anyDuplicated(i)) {
    stop("features must contain unique known feature IDs or valid indices.")
  }
  cached <- fit$marginal_data$result
  if (is.data.frame(cached) && all(ids[i] %in% cached$feature_id)) {
    out <- cached[match(ids[i], cached$feature_id), , drop = FALSE]
    rownames(out) <- NULL
    return(out)
  }
  fit <- .inlast_sparse_prepare(fit)
  ans <- data.frame(feature_id = ids[i], statistic = NA_real_, p_value = NA_real_,
    method_requested = "liu", method_used = "liu", fallback_used = FALSE,
    fallback_reason = NA_character_, davies_ifault = NA_integer_,
    error_message = NA_character_, stringsAsFactors = FALSE)
  for (rows in split(seq_along(i), ceiling(seq_along(i) / chunk_size))) {
    z <- .inlast_sparse_batch(fit, i[rows], threads, null_target = TRUE)
    for (j in seq_along(rows)) {
      k <- rows[j]
      if (!is.null(z[[j]]$error) && nzchar(z[[j]]$error)) {
        ans$error_message[k] <- z[[j]]$error
        next
      }
      ans$statistic[k] <- z[[j]]$statistic
      ans$p_value[k] <- .mgcvst_marginal_liu(z[[j]]$statistic, z[[j]]$moments)
      if (!is.finite(ans$p_value[k])) ans$error_message[k] <- "Invalid marginal Liu p-value."
    }
  }
  ans
}
