# Retain exactly one shared nuisance design and one small covariance per feature.
.mgcvst_nuisance_state <- function(fit, geometry, cache) {
  family <- .working_family_id(fit$family$family)
  if (is.null(cache$L)) return(list(error = "missing geometry cache"))
  if (!(family %in% c("gaussian", "poisson", "quasipoisson",
                      "negative_binomial"))) {
    return(list(error = "unsupported family"))
  }
  if (fit$rank != length(fit$coefficients)) {
    return(list(error = sprintf(
      "rank-deficient fit (rank %d of %d)", fit$rank, length(fit$coefficients)
    )))
  }
  tested <- sort(unique(unlist(lapply(
    geometry$smooth[geometry$target], `[[`, "columns"), use.names = FALSE
  )))
  nuisance <- setdiff(seq_len(ncol(cache$L)), tested)
  if (!length(nuisance)) return(list(error = "no nuisance columns"))
  if (!all(dim(fit$Vp) == ncol(cache$L))) {
    return(list(error = "Vp dimension mismatch"))
  }
  VpN <- as.matrix(fit$Vp[nuisance, nuisance, drop = FALSE])
  if (any(!is.finite(VpN)) ||
      !isTRUE(isSymmetric(VpN, tol = 100 * .Machine$double.eps))) {
    return(list(error = "non-finite or asymmetric nuisance Vp block"))
  }
  list(columns = nuisance, design = cache$L[, nuisance, drop = FALSE],
       covariance = VpN)
}

.mgcvst_model_sp <- function(fit) {
  sp <- fit$sp
  if (!is.null(fit$full.sp) && length(fit$full.sp) == length(fit$sp)) sp <- fit$full.sp
  as.numeric(sp)
}

# The shared prediction geometry of the model (seeded from mgcvST.set()) with
# the smoothing parameters and offset of this fit.
.mgcvst_cached_model_geometry <- function(fit, cache) {
  geometry <- cache$geometry
  geometry$sp <- .mgcvst_model_sp(fit)
  geometry$offset <- fit$offset
  geometry
}
