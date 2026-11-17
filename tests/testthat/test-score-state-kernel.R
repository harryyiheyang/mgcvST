# The dense score state against its definition, P = W - W X K X' W with
# W = (D + F F')^{-1}: a = F' P e and H = F' P F.
test_that("the dense score state equals its direct definition", {
  set.seed(11)
  n <- 120L; q <- 9L; p <- 2L; G <- 3L
  T0 <- matrix(rnorm(n * q), n, q) / sqrt(n)
  X <- cbind(1, rnorm(n))
  V <- matrix(rexp(n * G) + 0.2, n, G)
  E <- matrix(rnorm(n * G), n, G)
  sc <- c(0.5, 3, 40)
  Vp <- lapply(seq_len(G), function(j) crossprod(matrix(rnorm(p * p), p)) / 20)
  direct <- function(j, K) {
    Fm <- sqrt(sc[j]) * T0
    W <- solve(diag(V[, j]) + tcrossprod(Fm))
    P <- if (is.null(K)) W else {
      if (identical(K, "pinv")) K <- solve(crossprod(X, W %*% X))
      W - W %*% X %*% K %*% crossprod(X, W)
    }
    H <- crossprod(Fm, P %*% Fm)
    list(a = drop(crossprod(Fm, P %*% E[, j])), H = (H + t(H)) / 2)
  }
  check <- function(z, K) for (j in seq_len(G)) {
    ref <- direct(j, if (is.list(K)) K[[j]] else K)
    expect_null(z[[j]]$error)
    expect_equal(z[[j]]$a, ref$a, tolerance = 1e-9)
    if (!is.null(z[[j]]$H)) expect_equal(z[[j]]$H, ref$H, tolerance = 1e-9)
  }
  f <- mgcvst_dense_score_batch_cpp
  check(f(T0, V, E, sc, X, Vp, 1L), Vp)
  check(f(T0, V, E, sc, X, list(), 1L), "pinv")
  check(f(T0, V, E, sc, X[, 0L, drop = FALSE], list(), 1L), NULL)
  z <- f(T0, V, E, sc, X, Vp, 1L, TRUE)
  check(z, Vp)
  expect_null(z[[1L]]$H)
  V[1L, 2L] <- 0
  z <- f(T0, V, E, sc, X, Vp, 1L)
  expect_match(z[[2L]]$error, "strictly positive")
  expect_null(z[[1L]]$error)
})

test_that("the thread limit is restored by the process that set it", {
  vars <- .mgcvst_thread_vars
  before <- Sys.getenv(vars, unset = NA_character_, names = TRUE)
  on.exit({
    Sys.unsetenv(vars[is.na(before)])
    if (any(!is.na(before))) do.call(Sys.setenv, as.list(before[!is.na(before)]))
  })
  Sys.setenv(OMP_NUM_THREADS = "3")
  Sys.unsetenv("MKL_NUM_THREADS")
  old <- .mgcvst_thread_limit()
  expect_identical(Sys.getenv("OMP_NUM_THREADS"), "1")
  expect_identical(Sys.getenv("MKL_NUM_THREADS"), "1")
  .mgcvst_thread_restore(old)
  expect_identical(Sys.getenv("OMP_NUM_THREADS"), "3")
  expect_identical(Sys.getenv("MKL_NUM_THREADS", unset = NA_character_), NA_character_)
  old <- .mgcvst_thread_limit()
  old$pid <- -1L
  .mgcvst_thread_restore(old)
  expect_identical(Sys.getenv("OMP_NUM_THREADS"), "1")
})
