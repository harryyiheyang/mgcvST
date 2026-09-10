# Extract the INLA conditional Gaussian posterior covariance needed by the
# mgcvST nuisance projection.  This is an experimental helper, not production
# code.  It never constructs predictor- or observation-space covariance.

.inlast_vp_symmetric_precision <- function(Q, tolerance = 1e-12) {
  if (!requireNamespace("Matrix", quietly = TRUE)) stop("Matrix is required.")
  Q <- methods::as(Q, "CsparseMatrix")
  if (nrow(Q) != ncol(Q) || any(!is.finite(Q@x))) {
    stop("The INLA configuration precision must be finite and square.")
  }
  if (inherits(Q, "symmetricMatrix")) return(Matrix::forceSymmetric(Q))
  upper_nnz <- Matrix::nnzero(Matrix::triu(Q, 1L))
  lower_nnz <- Matrix::nnzero(Matrix::tril(Q, -1L))
  if (upper_nnz && !lower_nnz) {
    return(Matrix::forceSymmetric(Q, uplo = "U"))
  }
  if (lower_nnz && !upper_nnz) {
    return(Matrix::forceSymmetric(Q, uplo = "L"))
  }
  if (isTRUE(Matrix::isSymmetric(Q, tol = tolerance))) {
    return(Matrix::forceSymmetric(Q))
  }
  stop("The INLA configuration precision has incompatible triangle storage.")
}

.inlast_vp_configuration <- function(fit) {
  configs <- fit$misc$configs
  if (is.null(configs) || length(configs$config) != 1L ||
      !identical(as.integer(configs$nconfig), 1L)) {
    stop("The posterior-Vp helper requires one empirical-Bayes configuration.")
  }
  config <- configs$config[[1L]]
  if (is.null(config$Q) || is.null(config$mean)) {
    stop("The INLA fit did not retain its conditional Gaussian configuration.")
  }
  Q <- .inlast_vp_symmetric_precision(config$Q)
  if (nrow(Q) != length(config$mean)) {
    stop("The INLA configuration mean and precision dimensions disagree.")
  }
  list(configs = configs, config = config, Q = Q)
}

# Map one INLA-generated tag and within-block position into config$Q.  The
# predictor blocks are present in fit$mode$x and configs$contents, but absent
# from config$Q; configs$mnpred (APredictor plus Predictor) is the offset.
.inlast_vp_tag_index <- function(configs, tag, within, q) {
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

.inlast_vp_nuisance_indices <- function(fit, spec) {
  z <- .inlast_vp_configuration(fit)
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
    .inlast_vp_tag_index(z$configs, tag, within, nrow(z$Q))
  }, integer(1L))
  if (anyDuplicated(index)) stop("The nuisance map contains duplicated latent indices.")
  index
}

# Compute selected columns/rows of the conditional Gaussian covariance.  INLA
# stores a sparse selected inverse in config$Qinv; its missing off-diagonal
# entries are not zero, so selected columns must be obtained from config$Q.
inlast_posterior_covariance_selected <- function(
    fit, index, labels = NULL, tolerance = 1e-10) {
  z <- .inlast_vp_configuration(fit)
  q <- nrow(z$Q)
  index <- as.integer(index)
  if (!length(index) || anyNA(index) || any(index < 1L) ||
      any(index > q) || anyDuplicated(index)) {
    stop("index must contain unique valid latent configuration indices.")
  }
  if (!is.null(labels) && length(labels) != length(index)) {
    stop("labels must align with index.")
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
      stop("The INLA configuration constraint is incompatible with config$Q.")
    }
    Matrix::Matrix(C0, sparse = TRUE)
  }
  rhs <- if (nrow(C)) cbind(E, Matrix::t(C)) else E
  solved <- Matrix::solve(z$Q, rhs)
  selected_columns <- solved[, seq_along(index), drop = FALSE]
  constraint_error <- 0
  if (nrow(C)) {
    Hinv_Ct <- solved[, length(index) + seq_len(nrow(C)), drop = FALSE]
    middle <- as.matrix(C %*% Hinv_Ct)
    selected_columns <- selected_columns - Hinv_Ct %*%
      solve(middle, as.matrix(C %*% selected_columns))
    constraint_error <- max(abs(as.matrix(C %*% selected_columns)))
  }
  covariance <- as.matrix(selected_columns[index, , drop = FALSE])
  covariance <- (covariance + t(covariance)) / 2
  if (!is.null(labels)) dimnames(covariance) <- list(labels, labels)
  if (any(!is.finite(covariance)) ||
      min(eigen(covariance, symmetric = TRUE, only.values = TRUE)$values) <
        -tolerance * max(1, max(abs(covariance)))) {
    stop("The selected INLA conditional covariance is not positive semidefinite.")
  }
  list(
    covariance = covariance,
    indices = index,
    selected_columns = selected_columns,
    constraint = C,
    diagnostics = list(
      constraint_covariance_error = constraint_error,
      latent_dimension = q,
      selected_dimension = length(index),
      predictor_dimension_excluded = as.integer(z$configs$mnpred),
      configurations = length(z$configs$config),
      conditional_on_hyperparameters = TRUE,
      hyperparameter_integration = "empirical Bayes single configuration"
    )
  )
}

# Extract the nuisance block in exactly spec$nuisance_design column order.
inlast_posterior_vp <- function(fit, spec, tolerance = 1e-10) {
  index <- .inlast_vp_nuisance_indices(fit, spec)
  labels <- colnames(spec$nuisance_design)
  if (is.null(labels)) labels <- paste0("nuisance", seq_along(index))
  answer <- inlast_posterior_covariance_selected(
    fit, index, labels = labels, tolerance = tolerance
  )
  answer$nuisance_covariance <- answer$covariance
  answer$nuisance_map <- spec$nuisance_map
  answer
}
