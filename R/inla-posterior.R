# Recover the symmetric sparse precision stored by one INLA configuration.
.inlast_config_precision <- function(Q, tolerance = 1e-12) {
  Q <- methods::as(Q, "CsparseMatrix")
  if (nrow(Q) != ncol(Q) || any(!is.finite(Q@x))) {
    stop("The INLA configuration precision must be finite and square.")
  }
  if (inherits(Q, "symmetricMatrix")) return(Matrix::forceSymmetric(Q))
  upper <- Matrix::nnzero(Matrix::triu(Q, 1L))
  lower <- Matrix::nnzero(Matrix::tril(Q, -1L))
  if (upper && !lower) return(Matrix::forceSymmetric(Q, uplo = "U"))
  if (lower && !upper) return(Matrix::forceSymmetric(Q, uplo = "L"))
  if (isTRUE(Matrix::isSymmetric(Q, tol = tolerance))) {
    return(Matrix::forceSymmetric(Q))
  }
  stop("The INLA configuration precision has incompatible triangle storage.")
}

.inlast_config <- function(fit) {
  configs <- fit$misc$configs
  if (is.null(configs) || length(configs$config) != 1L ||
      !identical(as.integer(configs$nconfig), 1L)) {
    stop("The INLA posterior covariance requires one empirical-Bayes configuration.")
  }
  config <- configs$config[[1L]]
  if (is.null(config$Q) || is.null(config$mean)) {
    stop("INLA did not retain its conditional Gaussian configuration.")
  }
  Q <- .inlast_config_precision(config$Q)
  if (nrow(Q) != length(config$mean)) {
    stop("The INLA configuration mean and precision dimensions disagree.")
  }
  list(configs = configs, config = config, Q = Q)
}

# Predictor blocks occur first in fit$mode$x but are absent from config$Q;
# mnpred is the combined APredictor and Predictor length.
.inlast_config_tag_index <- function(configs, tag, within, q) {
  contents <- configs$contents
  if (is.null(contents$tag) || is.null(contents$start) ||
      is.null(contents$length) || is.null(configs$mnpred)) {
    stop("The INLA fit lacks latent block indexing metadata.")
  }
  hit <- which(contents$tag == tag)
  within <- as.integer(within)
  if (length(hit) != 1L || length(within) != 1L || is.na(within) ||
      within < 1L || within > contents$length[hit]) {
    stop("Invalid INLA latent index for tag '", tag, "'.")
  }
  index <- contents$start[hit] + within - 1L - as.integer(configs$mnpred)
  if (index < 1L || index > q) {
    stop("Tag '", tag, "' does not map into the configuration precision.")
  }
  index
}

.inlast_nuisance_config_index <- function(z, spec) {
  map <- spec$nuisance_map
  if (!is.list(map)) stop("spec$nuisance_map is required.")
  index <- vapply(map, function(item) {
    source <- as.character(item$source)[1L]
    if (identical(source, "fixed")) {
      tag <- paste0(".inlast_x", as.integer(item$index))
      within <- 1L
    } else if (identical(source, "random")) {
      tag <- paste0(".inlast_r", as.integer(item$block))
      within <- as.integer(item$index)
    } else {
      stop("Unknown nuisance-map source.")
    }
    .inlast_config_tag_index(z$configs, tag, within, nrow(z$Q))
  }, integer(1L))
  if (anyDuplicated(index)) {
    stop("The nuisance map contains duplicated latent indices.")
  }
  index
}

# Selected covariance of the conditional Gaussian approximation at INLA's EB
# configuration. config$Qinv is a sparse selected inverse, so omitted entries
# cannot be treated as zero.
.inlast_posterior_covariance_selected <- function(
    fit, index, labels = NULL, tolerance = 1e-10, config = .inlast_config(fit)) {
  z <- config
  q <- nrow(z$Q)
  index <- as.integer(index)
  if (!length(index)) {
    return(list(
      covariance = matrix(numeric(), 0L, 0L), indices = integer(),
      constraint_covariance_error = 0, latent_dimension = q,
      predictor_dimension_excluded = as.integer(z$configs$mnpred)
    ))
  }
  if (anyNA(index) || any(index < 1L) || any(index > q) ||
      anyDuplicated(index)) {
    stop("The requested INLA covariance indices are invalid.")
  }
  if (!is.null(labels) && length(labels) != length(index)) {
    stop("Posterior covariance labels must align with its indices.")
  }
  E <- Matrix::sparseMatrix(
    i = index, j = seq_along(index), x = 1,
    dims = c(q, length(index))
  )
  constraints <- z$configs$constr
  C <- if (is.null(constraints) || is.null(constraints$nc) ||
           constraints$nc == 0L) {
    Matrix::Matrix(0, 0L, q, sparse = TRUE)
  } else {
    C0 <- as.matrix(constraints$A)
    if (ncol(C0) != q || nrow(C0) != as.integer(constraints$nc) ||
        any(!is.finite(C0))) {
      stop("The INLA configuration constraint is incompatible with its precision.")
    }
    Matrix::Matrix(C0, sparse = TRUE)
  }
  rhs <- if (nrow(C)) cbind(E, Matrix::t(C)) else E
  solved <- Matrix::solve(z$Q, rhs)
  selected <- solved[, seq_along(index), drop = FALSE]
  constraint_error <- 0
  if (nrow(C)) {
    Hinv_Ct <- solved[, length(index) + seq_len(nrow(C)), drop = FALSE]
    middle <- as.matrix(C %*% Hinv_Ct)
    selected <- selected - Hinv_Ct %*%
      solve(middle, as.matrix(C %*% selected))
    constraint_error <- max(abs(as.matrix(C %*% selected)))
  }
  covariance <- as.matrix(selected[index, , drop = FALSE])
  covariance <- (covariance + t(covariance)) / 2
  if (!is.null(labels)) dimnames(covariance) <- list(labels, labels)
  if (any(!is.finite(covariance)) ||
      min(eigen(covariance, symmetric = TRUE, only.values = TRUE)$values) <
        -tolerance * max(1, max(abs(covariance)))) {
    stop("The selected INLA posterior covariance is not positive semidefinite.")
  }
  list(
    covariance = covariance, indices = index,
    constraint_covariance_error = constraint_error,
    latent_dimension = q,
    predictor_dimension_excluded = as.integer(z$configs$mnpred)
  )
}

.inlast_posterior_vp <- function(fit, spec, tolerance = 1e-10) {
  config <- .inlast_config(fit)
  index <- .inlast_nuisance_config_index(config, spec)
  labels <- colnames(spec$nuisance_design)
  if (is.null(labels)) labels <- paste0("nuisance", seq_along(index))
  answer <- .inlast_posterior_covariance_selected(
    fit, index, labels = labels, tolerance = tolerance, config = config
  )
  list(
    nuisance_covariance = answer$covariance,
    indices = answer$indices,
    diagnostics = list(
      constraint_covariance_error = answer$constraint_covariance_error,
      latent_dimension = answer$latent_dimension,
      selected_dimension = length(index),
      predictor_dimension_excluded = answer$predictor_dimension_excluded,
      configurations = 1L,
      conditional_on_hyperparameters = TRUE,
      hyperparameter_integration = "empirical Bayes single configuration"
    )
  )
}
