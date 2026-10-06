# Reuse the successful null state only for a native INLA program crash.
.inlast_spatial_fallback <- function(error, null, spec, y, offset) {
  if (!inherits(error, "inlaCrashError") || inherits(null, "condition") ||
      !isTRUE(null$converged)) return(error)
  target <- which(vapply(spec$random, function(z) isTRUE(z$target), logical(1L)))
  if (length(target) != 1L) stop("The INLA spatial fallback requires one target.")
  blocks <- vapply(spec$random, `[[`, character(1L), "name")
  keep <- setdiff(seq_along(blocks), target)
  index <- match(blocks[keep], names(null$random_mode))
  if (anyNA(index)) stop("The INLA null fallback is missing a nuisance block.")
  z <- null
  z$random_mode <- stats::setNames(vector("list", length(blocks)), blocks)
  z$tau <- z$tau_internal <- z$precision_scale <- z$lambda <-
    stats::setNames(numeric(length(blocks)), blocks)
  z$constraint_residual <- z$constraint_residual_uncorrected <-
    z$observation_spatial_mean <- stats::setNames(rep(NA_real_, length(blocks)), blocks)
  if (length(keep)) {
    z$random_mode[keep] <- null$random_mode[index]
    for (name in c("tau", "tau_internal", "precision_scale", "lambda",
                   "constraint_residual", "constraint_residual_uncorrected",
                   "observation_spatial_mean")) z[[name]][keep] <- null[[name]][index]
  }
  block <- spec$random[[target]]
  z$random_mode[[target]] <- numeric(ncol(block$A))
  scale <- block$precision_scale
  if (is.null(scale)) scale <- 1
  z$tau[target] <- 1e8
  z$precision_scale[target] <- scale
  z$tau_internal[target] <- 1e8 / scale
  z$lambda[target] <- z$dispersion * 1e8
  z$constraint_residual[target] <- z$constraint_residual_uncorrected[target] <-
    z$observation_spatial_mean[target] <- 0
  z$random_mean <- z$random_mode
  z$smoothing_parameters <- rep(NA_real_, spec$geometry_sp_length)
  for (j in seq_along(spec$random)) {
    z$smoothing_parameters[spec$random[[j]]$sp_index] <- z$lambda[j]
  }
  width <- if (is.null(block$projection)) ncol(block$A) else ncol(block$projection)
  z$coefficients <- stats::setNames(list(numeric(width)), blocks[target])
  eta <- offset + as.numeric(spec$fixed$X %*% z$fixed_mode)
  for (j in keep) eta <- eta + as.numeric(spec$random[[j]]$A %*% z$random_mode[[j]])
  z$eta <- eta
  if (z$family == "gaussian") {
    z$mu <- eta
    z$working_error <- y - offset
    z$working_variance <- rep(z$dispersion, length(y))
  } else {
    z$mu <- exp(eta)
    z$working_error <- eta + (y - z$mu) / z$mu - offset
    z$working_variance <- if (z$family == "negative_binomial") {
      1 / z$mu + 1 / z$family_parameters
    } else 1 / z$mu
  }
  z$spatial_fallback <- list(
    method = "null_zero_spatial", spatial_precision = 1e8,
    spatial_precision_scale = "original_FEM", native_converged = FALSE,
    error_class = class(error), error_message = conditionMessage(error),
    error_call = paste(deparse(conditionCall(error)), collapse = " ")
  )
  z$mode_status <- NA_integer_
  z$mode_status_text <- "null_zero_spatial_fallback"
  z$log_marginal_likelihood <- NA_real_
  z$fit_seconds <- NA_real_
  z$inla <- NULL
  z$estimation$engine <- "INLA null fallback"
  z$estimation$latent_mode_source <- "successful null fit; spatial effect set to zero"
  z$estimation$spatial_fallback <- z$spatial_fallback
  z
}
