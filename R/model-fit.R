# Expand a penalty to the coefficient dimension used by one smooth.
.mgcvst_expand_penalty <- function(S, n_coef, label) {
  S <- as.matrix(S)
  if (nrow(S) != ncol(S)) stop("Penalty for '", label, "' must be square.")
  if (nrow(S) == n_coef) return(S)
  if (!nrow(S) || n_coef %% nrow(S)) {
    stop("Penalty for '", label, "' is incompatible with its basis.")
  }
  kronecker(diag(n_coef %/% nrow(S)), S)
}

# Extract shared smooth geometry and fitted smoothing parameters.
.mgcvst_model_geometry <- function(fit, L = NULL) {
  if (is.null(L)) L <- .gam_training_lpmatrix(fit)
  smooth_columns <- lapply(
    fit$smooth, function(s) seq.int(s$first.para, s$last.para)
  )
  all_smooth <- unique(unlist(smooth_columns, use.names = FALSE))
  parametric <- setdiff(seq_len(ncol(L)), all_smooth)
  smooth <- vector("list", length(fit$smooth))
  sp_value <- fit$sp
  if (!is.null(fit$full.sp) && length(fit$full.sp) == length(fit$sp)) {
    sp_value <- fit$full.sp
  }
  for (j in seq_along(fit$smooth)) {
    s <- fit$smooth[[j]]
    columns <- smooth_columns[[j]]
    fixed <- isTRUE(s$fixed) || is.null(s$S) || !length(s$S) ||
      (!is.null(s$rank) && isTRUE(all(s$rank == 0)))
    sp_index <- if (fixed) integer() else seq.int(s$first.sp, s$last.sp)
    if (!fixed && length(sp_index) != length(s$S)) {
      stop("Smooth '", s$label, "' has inconsistent penalty indexes.")
    }
    penalties <- if (fixed) list() else lapply(
      s$S, .mgcvst_expand_penalty,
      n_coef = length(columns), label = s$label
    )
    smooth[[j]] <- list(
      label = s$label,
      B = as.matrix(L[, columns, drop = FALSE]),
      penalties = penalties,
      sp_index = sp_index,
      fixed = fixed,
      score_component = s$score.component,
      columns = columns
    )
  }
  score_component <- vapply(
    smooth,
    function(s) if (is.null(s$score_component)) "" else s$score_component,
    character(1L)
  )
  marked <- which(nzchar(score_component))
  if (length(marked) != 1L || sum(score_component == "global") != 1L) {
    stop("The model must mark exactly one global SPDE score component. The ",
         "second 'local' geographic process was removed from spdeST.")
  }
  target <- stats::setNames(marked, "global")
  for (name in names(target)) {
    z <- smooth[[target[[name]]]]
    if (z$fixed || length(z$penalties) != 1L) {
      stop("A score SPDE component must have one fitted full-rank penalty.")
    }
  }
  offset <- fit$offset
  if (is.null(offset)) offset <- numeric(nrow(L))
  list(
    X = as.matrix(L[, parametric, drop = FALSE]),
    smooth = smooth,
    target = target,
    score_components = names(target),
    offset = as.numeric(offset),
    row_id = .mgcvst_row_id(fit, nrow(L)),
    sp = as.numeric(sp_value)
  )
}

