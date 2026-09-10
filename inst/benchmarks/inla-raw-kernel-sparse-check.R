# Reproducible numerical checks for inla-raw-kernel-sparse.R.
#
# This script performs no INLA fits and records no timing comparison.  Its
# dense-factor reference uses observation-by-mesh matrices but never constructs
# an observation-by-observation matrix.
.libPaths(c(".test-library", .libPaths()))
suppressPackageStartupMessages(library(mgcvST))

source("inst/benchmarks/inla-raw-kernel-sparse.R")

out <- Sys.getenv(
  "MGCVST_RAW_SPARSE_CHECK_OUTPUT",
  "artifacts/constraint-type1/raw-kernel-sparse-check.csv"
)
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)

relative_error <- function(actual, expected) {
  max(abs(actual - expected)) / max(1, max(abs(expected)))
}

# Factor reference for the same constrained null P and raw tested kernel.  The
# constrained factor is formed by a rank-one update, not an m by m projector.
factor_reference <- function(A, Q, g, e, D, tau, X) {
  A <- as.matrix(A)
  Q <- as.matrix(Q)
  R <- chol(Q)
  Rinv <- backsolve(R, diag(ncol(Q)))
  Fu <- A %*% Rinv / sqrt(tau)
  h <- as.numeric(crossprod(Rinv, g))
  h_unit <- h / sqrt(sum(h^2))
  F0 <- Fu - tcrossprod(as.numeric(Fu %*% h_unit), h_unit)
  Fs <- Fu - qr.fitted(qr(X), Fu)

  Dinv <- 1 / D
  DinvF0 <- Dinv * F0
  woodbury <- diag(ncol(F0)) + crossprod(F0, DinvF0)
  Rw <- chol((woodbury + t(woodbury)) / 2)
  solve_woodbury <- function(B) {
    backsolve(Rw, forwardsolve(t(Rw), B))
  }
  Vsolve <- function(value) {
    value <- if (is.null(dim(value))) matrix(value, ncol = 1L) else value
    DinvY <- Dinv * value
    DinvY - DinvF0 %*% solve_woodbury(crossprod(F0, DinvY))
  }
  WX <- Vsolve(X)
  Vp <- solve(crossprod(X, WX), diag(ncol(X)))
  Papply <- function(value) {
    Wy <- Vsolve(value)
    Wy - WX %*% Vp %*% crossprod(X, Wy)
  }
  Pe <- Papply(e)
  PF <- Papply(Fs)
  M <- crossprod(Fs, PF)
  list(
    a = as.numeric(crossprod(Fs, Pe)),
    M = (M + t(M)) / 2,
    nuisance_covariance = Vp
  )
}

make_working_states <- function(data, feature = 1:2) {
  sx <- as.numeric(scale(data$x))
  sy <- as.numeric(scale(data$y))
  lapply(feature, function(j) {
    list(
      e = sin((0.17 + 0.03 * j) * sx) +
        cos((0.11 + 0.02 * j) * sy) + 0.05 * sx * sy,
      D = 0.35 + exp(0.12 * j * sx - 0.08 * sy),
      tau = c(0.8, 2.3)[j]
    )
  })
}

check_geometry <- function(name, A, Q, g, data, X) {
  states <- make_working_states(data)
  sparse <- reference <- vector("list", 2L)
  rows <- list()
  for (j in 1:2) {
    z <- states[[j]]
    sparse[[j]] <- inlast_raw_kernel_sparse_state(
      A, Q, g, z$e, z$D, z$tau, nuisance_X = X
    )
    reference[[j]] <- factor_reference(A, Q, g, z$e, z$D, z$tau, X)
    quantities <- list(
      a = c(relative_error(sparse[[j]]$a, reference[[j]]$a), 1e-9),
      M = c(relative_error(sparse[[j]]$M, reference[[j]]$M), 1e-9),
      nuisance_covariance = c(relative_error(
        sparse[[j]]$nuisance_covariance,
        reference[[j]]$nuisance_covariance
      ), 1e-9)
    )
    rows[[length(rows) + 1L]] <- do.call(rbind, lapply(
      names(quantities), function(quantity) data.frame(
        scenario = name, feature = j, quantity = quantity,
        relative_error = quantities[[quantity]][1L],
        tolerance = quantities[[quantity]][2L],
        passed = quantities[[quantity]][1L] <= quantities[[quantity]][2L],
        n = nrow(A), m = ncol(A), nuisance_dimension = ncol(X),
        minimum_M_eigenvalue =
          sparse[[j]]$diagnostics$minimum_M_eigenvalue,
        constraint_solve_error =
          sparse[[j]]$diagnostics$constraint_solve_error,
        nuisance_normal_equation_error =
          sparse[[j]]$diagnostics$nuisance_normal_equation_error,
        raw_direction_nuisance_residual_ratio = sparse[[j]]$diagnostics$
          raw_constraint_direction_nuisance_residual_ratio,
        largest_dense_dimension = sparse[[j]]$diagnostics$dense_mesh_dimension,
        observation_square_matrix_constructed =
          sparse[[j]]$diagnostics$observation_square_matrix_constructed,
        stringsAsFactors = FALSE
      )
    ))
  }
  sparse_pair <- inlast_raw_kernel_sparse_pair(sparse[[1]], sparse[[2]])
  reference_score <- as.numeric(crossprod(reference[[1]]$a, reference[[2]]$a))
  reference_moments <- mgcvST:::.rkhs_score_moments(
    reference[[1]]$M, reference[[2]]$M
  )
  pair_quantities <- list(
    signed_score = relative_error(sparse_pair$signed_score, reference_score),
    trace_moments = relative_error(sparse_pair$moments, reference_moments)
  )
  rows[[length(rows) + 1L]] <- do.call(rbind, lapply(
    names(pair_quantities), function(quantity) data.frame(
      scenario = name, feature = NA_integer_, quantity = quantity,
      relative_error = pair_quantities[[quantity]], tolerance = 1e-9,
      passed = pair_quantities[[quantity]] <= 1e-9,
      n = nrow(A), m = ncol(A), nuisance_dimension = ncol(X),
      minimum_M_eigenvalue = min(vapply(
        sparse, function(z) z$diagnostics$minimum_M_eigenvalue, numeric(1L)
      )),
      constraint_solve_error = max(vapply(
        sparse, function(z) z$diagnostics$constraint_solve_error, numeric(1L)
      )),
      nuisance_normal_equation_error = max(vapply(
        sparse,
        function(z) z$diagnostics$nuisance_normal_equation_error,
        numeric(1L)
      )),
      raw_direction_nuisance_residual_ratio = max(vapply(
        sparse,
        function(z) z$diagnostics$
          raw_constraint_direction_nuisance_residual_ratio,
        numeric(1L)
      )),
      largest_dense_dimension = ncol(A),
      observation_square_matrix_constructed = FALSE,
      stringsAsFactors = FALSE
    )
  ))
  list(rows = do.call(rbind, rows), sparse = sparse)
}

# Real 151673 observation/mesh geometry: n=3611, m=625, kappa=0.7.
large_input <- readRDS("artifacts/pathwaylgm-151673/benchmark-input.rds")
large_basis <- spde_basis(
  large_input$mesh, as.matrix(large_input$data[c("x", "y")]),
  kappa = 0.7, project_intercept = TRUE
)
large_model <- inlaST.set(
  response ~ offset(offset0), large_input$data, large_basis,
  family = mgcv::nb()
)
large_raw <- large_model$inla_spec$random[[1L]]
large_check <- check_geometry(
  "pathwaylgm_151673_k07",
  large_raw$A, large_raw$Q, large_raw$constraint,
  large_input$data, large_model$inla_spec$fixed$X
)
rm(large_basis, large_model, large_raw)
gc(FALSE)

# Estimated-null simulation geometry at kappa=6 with intercept+x nuisance.
small_input <- readRDS(
  "artifacts/constraint-type1/pilot/gaussian_pair_k6-dgp.rds"
)
small_basis <- spde_basis(
  small_input$mesh, as.matrix(small_input$data[c("x", "y")]),
  kappa = 6, project_intercept = TRUE
)
small_model <- inlaST.set(
  response ~ x + offset(offset0), small_input$data, small_basis,
  family = gaussian()
)
small_raw <- small_model$inla_spec$random[[1L]]
small_check <- check_geometry(
  "estimated_null_k6_intercept_x",
  small_raw$A, small_raw$Q, small_raw$constraint,
  small_input$data, small_model$inla_spec$fixed$X
)

answer <- rbind(large_check$rows, small_check$rows)
write.csv(answer, out, row.names = FALSE)
if (!all(answer$passed)) stop("A raw sparse/reference equivalence check failed.")
print(answer, row.names = FALSE)
