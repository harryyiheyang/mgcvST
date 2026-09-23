# Return variance-scale hyperparameters available for landmark selection.
.mgcvst_landmark_hyper <- function(fit, used) {
  if (!is.list(fit) || is.null(fit$feature_id) ||
      !is.numeric(used) || anyNA(used) || any(!is.finite(used)) ||
      any(used != as.integer(used)) || any(used < 1L) ||
      any(used > length(fit$feature_id)) || anyDuplicated(used)) {
    stop("fit and used must identify distinct valid feature indices.")
  }
  used <- as.integer(used)
  if (!length(used)) stop("used must contain at least one feature.")
  if (identical(fit$test_engine, "spde")) {
    spatial_variance <- .mgcvst_field_scale(fit)[used]
    if (any(!is.finite(spatial_variance)) || any(spatial_variance <= 0)) {
      stop("Legacy SPDE landmark scales must be finite and positive.")
    }
    values <- matrix(as.numeric(spatial_variance), ncol = 1L,
                     dimnames = list(fit$feature_id[used], "variance_spatial"))
    return(values)
  }
  phi <- fit$dispersion[used]
  sp <- fit$smoothing_parameters
  if (!is.numeric(phi) || length(phi) != length(used) ||
      any(!is.finite(phi)) || any(phi <= 0) ||
      !is.matrix(sp) || !is.numeric(sp) || nrow(sp) != length(fit$feature_id) ||
      any(!is.finite(sp[used, , drop = FALSE])) ||
      any(sp[used, , drop = FALSE] <= 0)) {
    stop("dispersion and smoothing parameters must be finite and positive.")
  }

  spec <- if (is.list(fit$model)) fit$model$inla_spec else NULL
  if (is.null(spec)) spec <- fit$inla_spec
  if (identical(fit$score_backend, "sparse") &&
      is.list(spec) && is.list(spec$random)) {
    blocks <- spec$random
    target <- which(vapply(blocks, function(z) isTRUE(z$target), logical(1L)))
    nuisance <- which(vapply(blocks, function(z) {
      !isTRUE(z$target) && identical(z$kind, "nuisance") &&
        identical(z$subtype, "iid")
    }, logical(1L)))
    if (length(target) != 1L) {
      stop("Sparse INLA landmark hyperparameters need one target SPDE block.")
    }
    selected <- c(target, nuisance)
    sp_index <- vapply(blocks[selected], function(z) as.integer(z$sp_index),
                       integer(1L))
    if (anyNA(sp_index) || any(sp_index < 1L) || any(sp_index > ncol(sp))) {
      stop("A sparse INLA random block has an invalid smoothing-parameter index.")
    }
    labels <- vapply(blocks[selected], function(z) {
      if (isTRUE(z$target)) "spatial" else {
        nm <- as.character(z$name)[1L]
        if (is.na(nm) || !nzchar(nm)) "iid" else nm
      }
    }, character(1L))
    labels <- make.unique(paste0("variance_", make.names(labels)))
    values <- vapply(seq_along(sp_index), function(j) {
      phi / sp[used, sp_index[j]]
    }, numeric(length(used)))
    if (is.null(dim(values))) values <- matrix(values, ncol = 1L)
    colnames(values) <- labels
    rownames(values) <- fit$feature_id[used]
    return(values)
  }

  if (identical(fit$score_backend, "sparse")) {
    stop("Sparse INLA landmark hyperparameters need the fitted random-block map.")
  }
  # Dense mgcv smoothing parameters are lambda = phi * tau; phi/lambda is
  # the corresponding covariance scale. Map indexes back to smooth labels.
  smooth <- fit$geometry$smooth
  if (!is.list(smooth) || !length(smooth)) {
    stop("Dense mgcv landmark hyperparameters need fitted smooth metadata.")
  }
  sp_index <- integer()
  labels <- character()
  for (j in seq_along(smooth)) {
    z <- smooth[[j]]
    index <- as.integer(z$sp_index)
    if (!length(index)) next
    if (anyNA(index) || any(index < 1L) || any(index > ncol(sp))) {
      stop("A dense mgcv smooth has an invalid smoothing-parameter index.")
    }
    label <- as.character(z$label)[1L]
    if (is.na(label) || !nzchar(label)) {
      stop("A dense mgcv smooth is missing its semantic label.")
    }
    sp_index <- c(sp_index, index)
    labels <- c(labels, if (length(index) == 1L) label else
      paste0(label, "#", seq_along(index)))
  }
  if (!length(sp_index) || anyDuplicated(sp_index)) {
    stop("Dense mgcv smooth metadata has no unique fitted smoothing parameters.")
  }
  labels <- make.unique(paste0("variance_", make.names(labels)))
  values <- sweep(1 / sp[used, sp_index, drop = FALSE], 1L, phi, FUN = "*")
  colnames(values) <- labels
  rownames(values) <- fit$feature_id[used]
  values
}

.mgcvst_landmark_select <- function(fit, used, n_ref,
                                    method = c("random", "score", "hyper"),
                                    seed, scores = NULL) {
  method <- match.arg(method)
  if (!is.list(fit) || is.null(fit$feature_id) ||
      !is.numeric(used) || anyNA(used) || any(!is.finite(used)) ||
      any(used != as.integer(used)) || any(used < 1L) ||
      any(used > length(fit$feature_id)) || anyDuplicated(used)) {
    stop("fit and used must identify distinct valid feature indices.")
  }
  used <- as.integer(used)
  if (!length(used)) stop("used must contain at least one feature.")
  if (length(n_ref) != 1L || !is.numeric(n_ref) || is.na(n_ref) ||
      !is.finite(n_ref) || n_ref != as.integer(n_ref) || n_ref < 1L) {
    stop("n_ref must be one positive integer.")
  }
  if (length(seed) != 1L || !is.numeric(seed) || is.na(seed) ||
      !is.finite(seed) || seed != as.integer(seed)) {
    stop("seed must be one finite integer.")
  }
  n_ref <- min(as.integer(n_ref), length(used))
  if (n_ref == length(used)) return(sort(used))

  old_exists <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (old_exists) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit(if (old_exists) {
    assign(".Random.seed", old_seed, envir = .GlobalEnv)
  } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  set.seed(as.integer(seed))
  if (method == "random") {
    return(sort(used[sample.int(length(used), n_ref)]))
  }

  standardize <- identical(method, "hyper")
  if (standardize) {
    x <- .mgcvst_landmark_hyper(fit, used)
  } else {
    if (!is.matrix(scores) || !is.numeric(scores) ||
        ncol(scores) != length(fit$feature_id) || is.null(colnames(scores)) ||
        anyDuplicated(colnames(scores)) ||
        !setequal(colnames(scores), fit$feature_id)) {
      stop("scores must be a named q-by-feature matrix aligned to fit$feature_id.")
    }
    x <- t(scores[, match(fit$feature_id[used], colnames(scores)), drop = FALSE])
    if (nrow(x) != length(used) || !ncol(x) || any(!is.finite(x))) {
      stop("scores for used features must be non-empty and finite.")
    }
  }
  if (any(!is.finite(x)) || !ncol(x)) {
    stop("landmark clustering features must be non-empty and finite.")
  }
  varying <- apply(x, 2L, function(z) length(unique(z)) > 1L)
  x <- x[, varying, drop = FALSE]
  if (!ncol(x)) {
    stop("landmark features are constant; use method = 'random'.")
  }
  if (standardize) x <- scale(x)
  if (any(!is.finite(x))) stop("standardized landmark features are not finite.")
  if (nrow(unique(as.data.frame(x))) < n_ref) {
    stop("There are fewer distinct feature vectors than requested landmarks.")
  }
  km <- stats::kmeans(x, centers = n_ref, iter.max = 100L, nstart = 5L)
  selected <- integer(n_ref)
  for (k in seq_len(n_ref)) {
    members <- which(km$cluster == k)
    if (!length(members)) stop("k-means returned an empty landmark cluster.")
    delta <- sweep(x[members, , drop = FALSE], 2L, km$centers[k, ], FUN = "-")
    distance <- rowSums(delta * delta)
    nearest <- members[which.min(distance)]
    selected[k] <- used[nearest]
  }
  if (anyDuplicated(selected)) stop("landmark representatives are not distinct.")
  sort(selected)
}
