# Adapted from mgcv.taps fb48abb4409a69dd5f68fa0929e0f8cb21a438be.
# Copyright (c) 2024 Yihe Yang; MIT. See inst/COPYRIGHTS.
# Frozen-fit marginal TAPS, NOT the pairwise covariance statistic.
.mgcvst_marginal_spectrum <- function(state, geometry, null.tol = 1e-10) {
  fit <- state
  if (is.null(fit[["family"]])) fit$family <- unserialize(state$family_raw)
  res <- .mgcvst_marginal_working(fit)
  offset <- if (is.null(state$offset)) geometry$offset else state$offset
  pseudo_response <- res$pseudo_response - offset
  V_phi <- res$V_phi
  phi0 <- res$phi0
  beta <- fit$coefficients
  X <- geometry$X
  smooth_terms <- geometry$smooth
  p <- length(smooth_terms)
  test.component <- geometry$test_component
  if (!is.null(res$valid_idx)) {
    idx <- res$valid_idx
    if (sum(idx) == 0) stop("No valid observations for testing")
    pseudo_response <- pseudo_response[idx]
    V_phi <- V_phi[idx]
    X <- X[idx, , drop = FALSE]
  }
  smooth_index_list <- list()
  random_index_list <- list()
  S_list            <- list()

  get_scaled_penalty <- function(s) {
    if (is.null(s$first.sp) || is.null(s$last.sp)) {
      stop("Penalized smooth has no smoothing-parameter index.")
    }
    sp_idx <- seq.int(s$first.sp, s$last.sp)
    if (length(sp_idx) != length(s$S)) {
      stop("Penalized smooth has inconsistent penalty and smoothing-parameter indexes.")
    }
    sp <- fit$sp[sp_idx]
    if (anyNA(sp)) {
      stop("Penalized smooth has missing smoothing-parameter values.")
    }
    Reduce(`+`, Map(function(S_matrix, sp_value) S_matrix * sp_value / phi0,
                    s$S, sp))
  }

  for (i in seq_len(p)) {
    s       <- smooth_terms[[i]]
    indices <- s$first.para:s$last.para
    reported_null_dim <- s$null.space.dim
    is_zero_rank <- !is.null(s$rank) && isTRUE(all(s$rank == 0))
    is_fixed_smooth <- isTRUE(s$fixed) || is.null(s$S) || length(s$S) == 0 || is_zero_rank

    if (i == test.component && is_fixed_smooth) {
      stop("test.component must be a penalized smooth with fx = FALSE.")
    }

    if (is_fixed_smooth) next

    S_matrix <- get_scaled_penalty(s)
    S_norm <- norm(S_matrix, "f")

    if (!is.finite(S_norm) || S_norm <= 0) {
      if (i == test.component) {
        stop("test.component must have a non-zero penalty matrix.")
      }
      next
    }

    if (is.null(s$getA) == 0) {
      detected_null_indices <- if (reported_null_dim > 0) indices[seq_len(reported_null_dim)] else integer(0)
    } else {
      col_norms <- sqrt(colSums(S_matrix^2))
      detected_null_indices <- indices[which(col_norms < null.tol)]
    }

    smooth_indices <- setdiff(indices, detected_null_indices)

    if (i != test.component) {
      smooth_index_list[[i]] <- indices
      random_index_list[[i]] <- smooth_indices
      S_list[[i]]            <- S_matrix
    }

    if (i == test.component) {
      Bj       <- X[, indices]
      Thetaj   <- CppMatrix::matrixGeneralizedInverse(S_matrix / S_norm)
      Gj_apply <- function(v) {
        v_is_matrix <- is.matrix(v)
        v_mat       <- if (v_is_matrix) v else matrix(v, ncol = 1)
        Bt_v        <- CppMatrix::matrixMultiply(Bj, v_mat, transA = TRUE)
        out         <- CppMatrix::matrixMultiply(Bj, CppMatrix::matrixMultiply(Thetaj, Bt_v))
        if (v_is_matrix) out else as.vector(out)
      }
      random_index_list[[i]] <- smooth_indices
    }
  }

  S_list            <- Filter(Negate(is.null), S_list)
  smooth_index_list <- Filter(Negate(is.null), smooth_index_list)
  random_index_list <- Filter(Negate(is.null), random_index_list)
  smooth_index_vec  <- do.call(c, smooth_index_list)
  random_index_vec  <- do.call(c, random_index_list)
  fixed_index_vec   <- setdiff(seq_len(ncol(X)), random_index_vec)

  A        <- as.matrix(X[, fixed_index_vec])
  alpha    <- beta[fixed_index_vec]
  B_extend <- X[, smooth_index_vec]
  S_All    <- as.matrix(Matrix::bdiag(S_list))
  XtX      <- CppMatrix::matrixMultiply(B_extend, B_extend * (1 / V_phi), transA = TRUE)
  C        <- CppMatrix::matrixInverse(XtX + S_All)

  Vinv_apply <- function(v) {
    v_is_matrix <- is.matrix(v)
    v_mat       <- if (v_is_matrix) v else matrix(v, ncol = 1)
    part1       <- v_mat / V_phi
    Bt_v        <- CppMatrix::matrixMultiply(B_extend, part1, transA = TRUE)
    C_Bt_v      <- CppMatrix::matrixMultiply(C, Bt_v)
    part2       <- CppMatrix::matrixMultiply(B_extend, C_Bt_v)
    out         <- part1 - part2 / V_phi
    if (v_is_matrix) out else as.vector(out)
  }

  Vinv_A      <- Vinv_apply(A)
  XtVinvX     <- CppMatrix::matrixMultiply(A, Vinv_A, transA = TRUE)
  XtVinvX_inv <- CppMatrix::matrixGeneralizedInverse(XtVinvX)

  P_apply <- function(v) {
    v_is_matrix  <- is.matrix(v)
    v_mat        <- if (v_is_matrix) v else matrix(v, ncol = 1)
    Vinv_v       <- Vinv_apply(v_mat)
    AVinv_v      <- CppMatrix::matrixMultiply(A, Vinv_v, transA = TRUE)
    solve_middle <- CppMatrix::matrixMultiply(XtVinvX_inv, AVinv_v)
    out          <- Vinv_v - CppMatrix::matrixMultiply(Vinv_A, solve_middle)
    if (v_is_matrix) out else as.vector(out)
  }

  error  <- pseudo_response - CppMatrix::matrixVectorMultiply(A, alpha)
  r      <- P_apply(error)
  Gj_r   <- Gj_apply(r)
  u      <- max(0,sum(r * Gj_r))
  q      <- ncol(Bj)
  eig_theta <- .mgcvst_marginal_matrixsqrt(Thetaj)
  Theta_sqrt <- eig_theta$w
  N      <- P_apply(Bj)
  BtPB   <- CppMatrix::matrixMultiply(Bj, N, transA = TRUE)
  Q_small <- CppMatrix::matrixListProduct(list(Theta_sqrt, BtPB, Theta_sqrt))
  lambda  <- eigen(Q_small, symmetric = TRUE, only.values = TRUE)$values
  lambda  <- lambda[lambda > 1e-15]

  list(statistic = u, lambda = lambda,
       smooth.term = smooth_terms[[test.component]]$label)
}
