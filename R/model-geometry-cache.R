# Only deterministic, supported prediction geometries may share an exact L.
# This guard is response-free; it is temporary and never stored per feature.
.mgcvst_geometry_signature <- function(fit) {
  known <- c("spde.smooth", "spdePC.smooth", "tprs.smooth", "cr.smooth",
             "pspline.smooth", "cp.smooth", "tensor.smooth", "random.effect")
  if (inherits(fit$family, "general.family") || is.null(fit$model) ||
      length(fit$paraPen) || (!is.null(fit$H) && any(fit$H != 0))) return(NULL)
  for (sm in fit$smooth) {
    cl <- class(sm)[1L]
    if (!(cl %in% known)) return(NULL)
    package <- if (cl %in% c("spde.smooth", "spdePC.smooth")) "mgcvST" else "mgcv"
    registered <- utils::getS3method("Predict.matrix", cl, optional = TRUE,
                                    envir = asNamespace("mgcv"))
    expected <- get0(paste0("Predict.matrix.", cl), asNamespace(package))
    if (is.null(expected) || !identical(registered, expected)) return(NULL)
  }
  terms <- fit$terms
  response <- attr(terms, "response")
  if (length(response) != 1L || response < 1L) return(NULL)
  attr(terms, ".Environment") <- NULL
  list(terms = terms, model = fit$model[-response], smooth = fit$smooth,
       nsdf = fit$nsdf, xlevels = fit$xlevels, contrasts = fit$contrasts,
       coefficient_names = names(fit$coefficients), offset = fit$offset)
}

# Retain exactly one shared nuisance design and one small covariance per feature.
.mgcvst_nuisance_state <- function(fit, geometry, cache) {
  family <- .working_family_id(fit$family$family)
  if (is.null(cache) || is.null(cache$L) ||
      (!isTRUE(cache$frozen) && (is.null(cache$signature) ||
       !identical(.mgcvst_geometry_signature(fit), cache$signature)))) {
    return(list(error = "missing geometry cache"))
  }
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

# Reuse the exact object obtained from the formal predictor, never G$X.
.mgcvst_training_design <- function(fit, cache = NULL) {
  if (!is.null(cache) && isTRUE(cache$frozen)) return(cache$L)
  if (!is.null(cache) && !is.null(cache$signature) &&
      identical(.mgcvst_geometry_signature(fit), cache$signature)) {
    return(cache$L)
  }
  .gam_training_lpmatrix(fit)
}

.mgcvst_cached_model_geometry <- function(fit, cache = NULL, L = NULL) {
  if (is.null(cache)) return(.mgcvst_model_geometry(fit, L = L))
  if (isTRUE(cache$frozen)) {
    geometry <- cache$geometry
    geometry$sp <- .mgcvst_model_sp(fit)
    geometry$offset <- fit$offset
    return(geometry)
  }
  signature <- .mgcvst_geometry_signature(fit)
  if (!is.null(signature) && identical(signature, cache$signature)) {
    geometry <- cache$geometry
    geometry$sp <- .mgcvst_model_sp(fit)
    return(geometry)
  }
  if (is.null(L)) L <- .gam_training_lpmatrix(fit)
  geometry <- .mgcvst_model_geometry(fit, L = L)
  if (!is.null(signature) && is.null(cache$geometry)) {
    cache$signature <- signature
    cache$geometry <- geometry
    cache$L <- L
  }
  geometry
}
