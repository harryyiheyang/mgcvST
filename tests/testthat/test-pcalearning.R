.pca_nb_fit <- local({
  cached <- NULL
  function() {
    skip_if_not_installed("INLA")
    skip_if_not_installed("geometry")
    if (!is.null(cached)) return(cached)
    withr::local_seed(1701L)
    n <- 72L
    vertices <- as.matrix(expand.grid(x = seq(0, 1, length.out = 5L),
                                      y = seq(0, 1, length.out = 5L)))
    mesh <- list(loc = vertices, graph = list(tv = geometry::delaunayn(vertices)))
    data <- data.frame(x = runif(n, 0.02, 0.98), y = runif(n, 0.02, 0.98),
                       z = seq(-1, 1, length.out = n), exposure = runif(n, 0.8, 1.3))
    data$offset0 <- log(data$exposure)
    basis <- spde_basis(mesh, as.matrix(data[c("x", "y")]), kappa = 1.2,
                        project_intercept = TRUE)
    G <- 8L
    Y <- t(vapply(seq_len(G), function(g) {
      eta <- 1 + 0.25 * data$z + data$offset0 +
        0.4 * sin(2 * pi * (data$x + g / G)) + 0.3 * cos(2 * pi * data$y * g / 4)
      rnbinom(n, mu = exp(eta), size = if (g %% 3 == 0) 1e4 else 2 + g)
    }, numeric(n)))
    dimnames(Y) <- list(paste0("g", seq_len(G)), NULL)
    model <- inlaST.set(response ~ z + offset(offset0), data, basis,
                        family = mgcv::nb())
    cached <<- inlaST.estimate(Y, model, BPPARAM = BiocParallel::SerialParam())
    cached
  }
})

.pca_unpack <- function(x, d) {
  M <- matrix(0, d, d)
  M[upper.tri(M, diag = TRUE)] <- x
  M[lower.tri(M)] <- t(M)[lower.tri(M)]
  M
}

test_that("PCAlearning strata and draws reproduce the approx-liu-p training set", {
  fx <- readRDS(test_path("fixtures", "pcalearning-sampling.rds"))
  G <- length(fx$cell)
  scales <- data.frame(nb = fx$nb, sigma_g2 = fx$sigma_g2, sigma_e2 = fx$sigma_e2)
  cell <- mgcvST:::.mgcvst_pca_cells(scales, seq_len(G))
  expect_identical(cell, fx$cell)

  # The experiment first drew its 150 test genes from the same RNG stream.
  elig <- which(!fx$bench100)
  N <- tabulate(cell[elig], 100L)
  withr::local_seed(20260924L)
  n2 <- min(150L - sum(N > 0), sum(N >= 2))
  two <- sample(which(N >= 2), n2)
  qT <- pmin(ifelse(seq_len(100L) %in% two, 2L, 1L), N)
  Tset <- integer()
  for (k in which(qT > 0)) {
    cand <- elig[cell[elig] == k]
    Tset <- c(Tset, cand[sample.int(length(cand), qT[k])])
  }
  expect_identical(sort(Tset), fx$T)
  S <- mgcvST:::.mgcvst_pca_draw(cell, setdiff(elig, Tset), 3L)
  expect_identical(sort(S), fx$S)
})

test_that("PCAlearning training sampling reallocates sparse cells and restores RNG", {
  withr::local_seed(3L)
  G <- 400L
  scales <- data.frame(nb = seq_len(G) > 40L, sigma_g2 = exp(rnorm(G)),
                       sigma_e2 = 1 + exp(rnorm(G, -2)))
  scales$sigma_e2[!scales$nb] <- 1
  before <- .Random.seed
  a <- mgcvST:::.mgcvst_pca_training(scales, seq_len(G), 3L, 11L)
  expect_identical(.Random.seed, before)
  b <- mgcvST:::.mgcvst_pca_training(scales, seq_len(G), 3L, 11L)
  expect_identical(a, b)
  expect_length(a$train, 300L)
  expect_false(anyDuplicated(a$train) > 0L)
  N <- tabulate(a$cell, 100L)
  drawn <- tabulate(a$cell[a$train], 100L)
  expect_true(all(drawn <= N))
  expect_true(all(drawn >= pmin(3L, N)))
  small <- mgcvST:::.mgcvst_pca_training(scales, seq_len(50L), 3L, 11L)
  expect_identical(small$train, seq_len(50L))
})

test_that("liu_approximation = 'exact' keeps the exact INLA Liu path", {
  skip_on_cran()
  fit <- .pca_nb_fit()
  pairs <- t(combn(fit$feature_id, 2L))
  exact <- inlaST.test(fit, pairs = pairs, method = "BY", threads = 2L)
  none <- inlaST.test(fit, pairs = pairs, method = "BY", threads = 2L,
                      liu_approximation = "exact")
  expect_identical(none$results, exact$results)
  expect_null(exact$pca_learning)

  prepared <- mgcvST:::.inlast_sparse_prepare(fit)
  basis <- mgcvST:::.inlast_sparse_observation_basis(prepared)
  units <- mgcvST:::.inlast_sparse_units(prepared, seq_along(fit$feature_id), threads = 1L)
  states <- mgcvST:::.inlast_sparse_materialize_reduced(prepared, units, basis, threads = 1L)
  local <- matrix(match(pairs, fit$feature_id), ncol = 2L)
  moments <- mgcvST:::mgcvst_pair_trace_powers_cpp(lapply(states, `[[`, "M"), local, 4L, 1L)
  score <- vapply(seq_len(nrow(local)), function(j) {
    sum(states[[local[j, 1L]]]$a * states[[local[j, 2L]]]$a)
  }, numeric(1L))
  liu <- mgcvST:::.liu_squared_score_moments(abs(score), moments[, 1L], moments[, 2L],
                                             moments[, 3L], moments[, 4L])
  expect_equal(exact$results$signed_score, score, tolerance = 1e-10)
  expect_equal(exact$results$p_two_sided, liu$p_value, tolerance = 1e-10)
  expect_equal(exact$results$p_adjusted, p.adjust(exact$results$p_two_sided, "BY"))
})

test_that("full-rank PCAlearning reproduces exact Liu p-values", {
  skip_on_cran()
  fit <- .pca_nb_fit()
  G <- length(fit$feature_id)
  pairs <- t(combn(fit$feature_id, 2L))
  exact <- inlaST.test(fit, pairs = pairs, method = "BY", threads = 2L)
  withr::local_seed(5L)
  before <- .Random.seed
  pca <- inlaST.test(fit, pairs = pairs, method = "BY", threads = 2L,
                     liu_approximation = "pca_learning", rank = G)
  expect_identical(.Random.seed, before)
  z <- pca$pca_learning
  expect_identical(z$training$feature_id, fit$feature_id)
  expect_identical(dim(z$coefficients), c(G, G))
  # training matrices are stored in float32
  expect_lt(max(abs(z$genes$e2_relative)), 1e-6)
  expect_equal(pca$results$signed_score, exact$results$signed_score, tolerance = 1e-10)
  expect_equal(pca$results$p_two_sided, exact$results$p_two_sided, tolerance = 1e-5)
  expect_equal(pca$results$p_positive, exact$results$p_positive, tolerance = 1e-5)
  expect_equal(pca$results$p_negative, exact$results$p_negative, tolerance = 1e-5)
  expect_equal(exp(pca$results$log_p_two_sided), pca$results$p_two_sided)
  expect_equal(pca$results$p_adjusted, p.adjust(pca$results$p_two_sided, "BY"),
               tolerance = 1e-12)
  expect_equal(pca$results$p_negative_adjusted,
               p.adjust(pca$results$p_negative, "BY"), tolerance = 1e-12)
  expect_identical(pca$results$discovered, exact$results$discovered)
  expect_identical(pca$results$discovered_positive, exact$results$discovered_positive)
  expect_identical(pca$results$discovered_negative, exact$results$discovered_negative)
  expect_true(all(c("sample", "gram", "tables", "pairs") %in% names(z$elapsed)))
  expect_true("total" %in% names(z$table_timing))

  # Gene blocks of 32: more genes than one block are all materialized.
  prepared <- mgcvST:::.inlast_sparse_prepare(fit)
  basis <- mgcvST:::.inlast_sparse_observation_basis(prepared)
  many <- mgcvST:::.mgcvst_pca_materialize(prepared, rep(seq_len(G), 5L), basis, 2L,
                                            pack = TRUE)
  expect_false(anyNA(many$A))
  expect_false(any(vapply(many$packed, is.null, logical(1L))))
  expect_equal(many$A[, 33:40], many$A[, 1:8])
  expect_identical(pca$timing$inla_projection$pair_schedule, "pcalearning_gene_blocks")

  # A pair list (reversed order, subset) uses the (i, j) kernel with equal results.
  sub <- c(5L, 1L, 20L, 13L)
  listed <- inlaST.test(fit, pairs = pairs[sub, 2:1], method = "BY", threads = 2L,
                        liu_approximation = "pca_learning", rank = G)
  expect_identical(listed$timing$inla_projection$pair_schedule, "pcalearning_pair_list")
  expect_equal(listed$results$log_p_two_sided, pca$results$log_p_two_sided[sub],
               tolerance = 1e-12)
  expect_equal(listed$results$log_p_positive, pca$results$log_p_positive[sub],
               tolerance = 1e-12)

  low <- inlaST.test(fit, pairs = pairs, method = "BY", threads = 2L,
                     liu_approximation = "pca_learning", rank = 3L)
  expect_true(all(low$pca_learning$genes$e2_relative > -1e-12))
  expect_equal(low$pca_learning$genes$e2,
               low$pca_learning$genes$fro2 - unname(rowSums(low$pca_learning$coefficients^2)))
})

test_that("log-space BY matches p.adjust and keeps underflowing tails", {
  withr::local_seed(9L)
  lp <- log(runif(500)) * rexp(500, 0.2)
  by <- mgcvST:::.mgcvst_log_by(lp)
  expect_equal(exp(by), p.adjust(exp(lp), "BY"), tolerance = 1e-12)
  expect_identical(by <= log(0.05), p.adjust(exp(lp), "BY") <= 0.05)
  deep <- mgcvST:::.mgcvst_log_by(c(-3000, -2000, -1))
  expect_true(all(is.finite(deep)))
  expect_lt(deep[1L], deep[2L])
})

test_that("rank-10 trace tables reproduce brute-force projected traces", {
  withr::local_seed(17L)
  q <- 6L
  r <- 10L
  L <- q * (q + 1L) / 2L
  B <- qr.Q(qr(matrix(rnorm(L * r), L, r)))
  unpack <- function(x) {
    M <- matrix(0, q, q)
    M[upper.tri(M, diag = TRUE)] <- x
    off <- upper.tri(M)
    M[off] <- M[off] / sqrt(2)
    M[lower.tri(M)] <- t(M)[lower.tri(M)]
    M
  }
  tables <- mgcvST:::mgcvst_pca_tables_cpp(B, q, 2L)
  expect_identical(dim(tables$Tsym4), rep(as.integer(choose(r + 3L, 4L)), 2L))
  C <- matrix(rnorm(6L * r), 6L, r)
  A <- matrix(rnorm(3L * 6L), 3L, 6L)
  i <- c(1L, 2L, 3L, 5L)
  j <- c(2L, 4L, 6L, 1L)
  out <- mgcvST:::mgcvst_pca_pairs_cpp(A, C, tables, i, j, 2L, moments = TRUE)
  H <- lapply(seq_len(nrow(C)), function(k) unpack(B %*% C[k, ]))
  brute <- t(vapply(seq_along(i), function(k) {
    M <- H[[i[k]]] %*% H[[j[k]]]
    M2 <- M %*% M
    c(sum(diag(M)), sum(diag(M2)), sum(diag(M2 %*% M)), sum(diag(M2 %*% M2)))
  }, numeric(4L)))
  expect_equal(unname(out[, c("t1", "t2", "t3", "t4")]), brute, tolerance = 1e-5)
  expect_equal(out[, "U"], colSums(A[, i] * A[, j]), tolerance = 1e-12)
})

test_that("PCAlearning pair kernel reproduces the approx-liu-p rank-10 traces", {
  fx <- readRDS(test_path("fixtures", "pcalearning-approx-liu-p.rds"))
  r <- fx$rank
  d <- choose(r + 1:4 - 1L, 1:4)
  tables <- list(Tsym1 = fx$tables$Tsym1)
  for (s in 2:4) {
    tables[[paste0("Tsym", s)]] <- .pca_unpack(fx$tables[[paste0("Tsym", s)]], d[s])
  }
  P <- fx$pairs
  out <- mgcvST:::mgcvst_pca_pairs_cpp(fx$A, fx$C, tables, P$i, P$j, 2L)
  expect_equal(out[, "U"], P$U, tolerance = 1e-10)
  approx <- as.matrix(P[c("t1_approx", "t2_approx", "t3_approx", "t4_approx")])
  expect_equal(unname(out[, c("t1", "t2", "t3", "t4")]), unname(approx), tolerance = 1e-6)

  exact <- mgcvST:::.liu_squared_score_moments(abs(P$U), P$t1, P$t2, P$t3, P$t4)
  expected <- mgcvST:::.liu_squared_score_moments(abs(P$U), P$t1_approx, P$t2_approx,
                                                  P$t3_approx, P$t4_approx)
  p_approx <- exp(out[, "logp_two_sided"])
  expect_equal(p_approx, expected$p_value, tolerance = 1e-6)
  delta <- -log10(p_approx) + log10(exact$p_value)
  expect_equal(delta, -log10(expected$p_value) + log10(exact$p_value), tolerance = 1e-5)
  # Bounds of the full 11,175-pair test-set Delta (S, r = 10): [-0.0546, 0.4535].
  expect_true(all(delta >= -0.0547 & delta <= 0.4536))
  expect_identical(p.adjust(p_approx, "BY") <= 0.05,
                   p.adjust(exact$p_value, "BY") <= 0.05)
  e2 <- (fx$fro2 - rowSums(fx$C^2)) / fx$fro2
  expect_true(all(e2 > 0 & e2 <= 0.0267))
})

test_that("PCAlearning checkpoints resume to the uninterrupted result", {
  skip_on_cran()
  fit <- .pca_nb_fit()
  G <- length(fit$feature_id)
  local_mocked_bindings(
    .mgcvst_pca_training = function(scales, universe, n_per_cell, seed) {
      list(train = universe[c(1L, 3L, 5L, 7L)], cell = rep(1L, nrow(scales)))
    },
    .package = "mgcvST"
  )
  basis <- mgcvST:::.inlast_sparse_observation_basis(mgcvST:::.inlast_sparse_prepare(fit))
  run <- function(index, dir, rank = 3L, resume = TRUE, n_per_cell = 3L) {
    mgcvST:::.mgcvst_pair_pcalearning(
      fit, index, seq_len(nrow(index)), 2L, 1000L, FALSE, basis, rank = rank,
      n_per_cell = n_per_cell, seed = 1L, checkpoint_dir = dir, resume = resume
    )
  }
  index <- t(combn(G, 2L))
  reference <- run(index, NULL)
  expect_null(reference$metadata$pca_learning$checkpoint$path)

  dir <- tempfile("mgcvst-pca-checkpoint-")
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  # Interrupted run: only genes 1-6 are requested, so genes 2, 4, 6 are projected.
  part <- run(t(combn(6L, 2L)), dir)
  expect_identical(part$metadata$pca_learning$checkpoint$projected_genes, 3L)
  expect_setequal(list.files(dir), c("manifest.rds", "pca-basis.rds",
                                      "pca-projection-000001.rds"))
  same <- function(x) {
    expect_identical(x$result, reference$result)
    expect_identical(x$metadata$pca_learning$coefficients,
                     reference$metadata$pca_learning$coefficients)
    expect_identical(x$metadata$pca_learning$genes, reference$metadata$pca_learning$genes)
  }
  resumed <- run(index, dir)
  z <- resumed$metadata$pca_learning$checkpoint
  expect_true(z$basis_resumed)
  expect_identical(c(z$resumed_genes, z$projected_genes), c(3L, 1L))
  same(resumed)

  # A lost block is recomputed into a new block.
  unlink(file.path(dir, "pca-projection-000001.rds"))
  again <- run(index, dir)
  expect_identical(again$metadata$pca_learning$checkpoint$projected_genes, 3L)
  expect_true(file.exists(file.path(dir, "pca-projection-000003.rds")))
  same(again)
  complete <- run(index, dir)
  expect_identical(complete$metadata$pca_learning$checkpoint$projected_genes, 0L)
  same(complete)

  expect_error(run(index, dir, rank = 2L), "different fit, score basis, rank")
  expect_error(run(index, dir, resume = FALSE), "already exists")

  public <- tempfile("mgcvst-pca-public-")
  on.exit(unlink(public, recursive = TRUE), add = TRUE)
  test <- function() inlaST.test(fit, pairs = t(combn(fit$feature_id, 2L)),
                                 threads = 2L, liu_approximation = "pca_learning",
                                 rank = 3L, checkpoint_dir = public, resume = TRUE)
  first <- test()
  expect_true(file.exists(file.path(public, "pca-basis.rds")))
  expect_identical(test()$results, first$results)
})

test_that("PCAlearning checks rank, n_per_cell and trace-table memory", {
  skip_on_cran()
  fit <- .pca_nb_fit()
  G <- length(fit$feature_id)
  basis <- mgcvST:::.inlast_sparse_observation_basis(mgcvST:::.inlast_sparse_prepare(fit))
  run <- function(rank, n_per_cell = 3L) {
    index <- t(combn(G, 2L))
    mgcvST:::.mgcvst_pair_pcalearning(
      fit, index, seq_len(nrow(index)), 2L, 1000L, FALSE, basis, rank = rank,
      n_per_cell = n_per_cell, seed = 1L
    )
  }
  expect_error(run(G + 1L), sprintf("achievable PCAlearning rank %d (%d training", G, G), fixed = TRUE)
  expect_error(run(0L), "rank must be one positive integer")
  expect_error(run(2.5), "rank must be one positive integer")
  expect_error(run(3L, n_per_cell = 0L), "n_per_cell must be one positive integer")
  # q = 1404, r = 10: level-4 stage, about 5.1 GiB (5.4 GB peak measured).
  expect_equal(mgcvST:::.mgcvst_pca_table_bytes(1404, 10) / 1024^3, 5.09, tolerance = 1e-2)
  local_mocked_bindings(.mgcvst_memory_probe = function(...) list(available = 1e3),
                        .package = "mgcvST")
  expect_error(run(3L), "trace tables for rank = 3 .* use a smaller rank")
})
