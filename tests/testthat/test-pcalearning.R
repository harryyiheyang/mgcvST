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

test_that("inlaST.test() runs the PCAlearning route on request and returns the compact result", {
  skip_on_cran()
  fit <- .pca_nb_fit()
  pairs <- t(combn(fit$feature_id, 2L))
  default <- inlaST.test(fit, pairs = pairs, threads = 2L, rank = 3L, seed = 4L,
                         moments = "pcalearning")
  expect_s3_class(default, "mgcvST_test")
  expect_identical(default$moments, "pcalearning")
  expect_false(is.null(default$pca_learning))
  expect_named(default$results, c("i", "j", "score", "log_p_two_sided",
    "log_p_positive", "log_p_negative", "log_q", "remainder_kind", "status"))
  expect_identical(default$feature_id[default$results$i], pairs[, 1L])
  expect_identical(default$feature_id[default$results$j], pairs[, 2L])
  expect_true(all(default$results$status == 0L))
  # The PCAlearning remainder is one node, or a Gaussian term.
  expect_true(all(default$results$remainder_kind %in% c(0L, 1L, 3L)))
  contract <- default$contract
  expect_identical(contract$calibration_contract, "spa_v1")
  expect_identical(contract$route, "pcalearning")
  expect_identical(contract$remainder_order, 2L)
  expect_identical(contract$k, min(50L, default$timing$route$q))
  expect_match(contract$basis_sha, "^[0-9a-f]{64}$")
  expect_identical(default$timing$pcalearning$pair_schedule, "pcalearning_pair_list")
  all_pairs <- inlaST.test(fit, threads = 2L, rank = 3L, seed = 4L, moments = "pcalearning")
  expect_identical(all_pairs$timing$pcalearning$pair_schedule, "pcalearning_gene_blocks")
  expect_equal(all_pairs$results, default$results, tolerance = 1e-12)
  expect_identical(all_pairs$pca_learning$training, default$pca_learning$training)
  expect_identical(all_pairs$contract$basis_sha, contract$basis_sha)
})

test_that("approximate inlaST.test() builds the observation basis once", {
  skip_on_cran()
  fit <- .pca_nb_fit()
  pairs <- t(combn(fit$feature_id, 2L))
  r <- mgcvST:::.inlast_sparse_observation_basis(
    mgcvST:::.inlast_sparse_prepare(fit)
  )$rank
  calls <- new.env(parent = emptyenv())
  calls$basis <- 0L
  observation_basis <- mgcvST:::.inlast_sparse_observation_basis
  testthat::local_mocked_bindings(
    .inlast_sparse_observation_basis = function(fit) {
      calls$basis <- calls$basis + 1L
      observation_basis(fit)
    },
    .package = "mgcvST")
  for (moments in c("pcalearning", "exact")) {
    calls$basis <- 0L
    out <- inlaST.test(fit, pairs = pairs, threads = 2L, rank = 3L, seed = 4L,
                       moments = moments)
    expect_identical(calls$basis, 1L)
    expect_identical(out$timing$route$q, r)
  }
  expect_identical(out$timing$inla_projection$r, r)
  expect_identical(out$timing$inla_projection$basis_kind, "full_rank")
})

test_that("inlaST.test shares the mgcvST.test arguments and rejects removed ones", {
  skip_on_cran()
  fit <- .pca_nb_fit()
  pair <- matrix(fit$feature_id[1:2], 1L)
  removed <- list(approximate_test = TRUE, pairwise_method = "conditional_cauchy",
                  liu_approximation = "exact", conditional_precision = "float32",
                  calibration = "liu", BPPARAM = BiocParallel::SerialParam(),
                  FDR = TRUE, method = "BY", highlight = pair, cache_bytes = 1e9)
  for (name in names(removed)) {
    expect_error(do.call(inlaST.test, c(list(fit, pairs = pair, rank = 2L), removed[name])),
                 "unused argument")
  }
  for (name in c("approximate_test", "pairwise_method", "liu_approximation",
                 "calibration", "conditional_precision", "BPPARAM", "FDR", "method",
                 "highlight", "cache_bytes", "...")) {
    expect_false(name %in% names(formals(inlaST.test)))
    expect_false(name %in% names(formals(mgcvST.test)))
  }
  shared <- c("pairs", "q.value", "adjust", "threads", "chunk_size",
              "checkpoint_dir", "resume", "verbose")
  # The shared arguments keep the same positions; the route and the PCAlearning
  # controls follow them, so a positional call means the same in both tests.
  controls <- c("moments", "rank", "n_per_cell", "seed", "k")
  expect_identical(names(formals(inlaST.test)), c("fitinlaST", shared, controls))
  expect_identical(names(formals(mgcvST.test)), c("fitmgcvST", shared, controls))
  expect_identical(formals(inlaST.test)[controls], formals(mgcvST.test)[controls])
  expect_equal(unlist(mgcvST:::.mgcvst_pca_defaults),
               c(rank = 20, n_per_cell = 3, seed = 1, k = 50))
  expect_identical(eval(formals(inlaST.test)$rank), 20L)
  expect_identical(eval(formals(inlaST.test)$moments), c("auto", "exact", "pcalearning"))
  expect_null(formals(inlaST.test)$k)
  expect_error(inlaST.test(fit, pairs = pair, adjust = "holm"), "should be one of")
  expect_error(inlaST.test(fit, pairs = pair, moments = "liu"), "should be one of")
  expect_error(inlaST.test(fit, pairs = pair, k = 0L), "k must be NULL or one positive")
  expect_false(exists(".mgcvst_conditional_test", asNamespace("mgcvST"),
                      inherits = FALSE))
  expect_false(exists("mgcvst_liu_logp_cpp", asNamespace("mgcvST"), inherits = FALSE))
})

test_that("full-rank PCAlearning with k = q reproduces the full-spectrum saddlepoint", {
  skip_on_cran()
  fit <- .pca_nb_fit()
  G <- length(fit$feature_id)
  pairs <- t(combn(fit$feature_id, 2L))
  exact <- .pca_exact_reference(fit, pairs)
  withr::local_seed(5L)
  before <- .Random.seed
  pca <- inlaST.test(fit, pairs = pairs, threads = 2L, rank = G, k = 1000L,
                     moments = "pcalearning")
  expect_identical(.Random.seed, before)
  z <- pca$pca_learning
  expect_identical(pca$contract$k, exact$q)
  expect_identical(z$training$feature_id, fit$feature_id)
  expect_identical(dim(z$coefficients), c(G, G))
  # training matrices are stored in float32
  expect_lt(max(abs(z$genes$e2_relative)), 1e-6)
  expect_equal(pca$results$i, exact$i)
  expect_equal(pca$results$j, exact$j)
  expect_equal(pca$results$score, exact$score, tolerance = 1e-10)
  expect_equal(-pca$results$log_p_two_sided / log(10), -exact$log_p / log(10),
               tolerance = 1e-3)
  expect_equal(pca$results$log_q, mgcvST:::.mgcvst_log_adjust(
    pca$results$log_p_two_sided, "BY")$log_q, tolerance = 1e-12)
  expect_equal(exp(pca$results$log_q), p.adjust(exp(pca$results$log_p_two_sided), "BY"),
               tolerance = 1e-12)
  expect_true(all(c("sample", "gram", "tables", "pairs") %in% names(z$elapsed)))
  expect_true("total" %in% names(z$table_timing))

  # Gene blocks of 32: more genes than one block are all materialized.
  prepared <- mgcvST:::.inlast_sparse_prepare(fit)
  basis <- mgcvST:::.inlast_sparse_observation_basis(prepared)
  producer <- mgcvST:::.mgcvst_pca_producer(prepared, seq_len(G), basis, 2L)
  many <- producer$materialize(rep(seq_len(G), 5L), pack = TRUE)
  expect_false(anyNA(many$A))
  expect_false(any(vapply(many$packed, is.null, logical(1L))))
  expect_equal(many$A[, 33:40], many$A[, 1:8])
  expect_identical(pca$timing$pcalearning$pair_schedule, "pcalearning_pair_list")

  # Every pair of the available genes is streamed by gene blocks, equal to the list.
  all_pairs <- inlaST.test(fit, threads = 2L, rank = G, k = 1000L, chunk_size = 5L,
                           moments = "pcalearning")
  expect_identical(all_pairs$timing$pcalearning$pair_schedule,
                   "pcalearning_gene_blocks")
  expect_gt(all_pairs$timing$pcalearning$chunks, 1L)
  expect_equal(all_pairs$results, pca$results, tolerance = 1e-12)

  # A pair list (reversed order, subset) uses the (i, j) kernel with equal results.
  sub <- c(5L, 1L, 20L, 13L)
  listed <- inlaST.test(fit, pairs = pairs[sub, 2:1], adjust = "none", threads = 2L,
                        rank = G, k = 1000L, moments = "pcalearning")
  expect_identical(listed$timing$pcalearning$pair_schedule, "pcalearning_pair_list")
  by_key <- function(x) paste(x$i, x$j)
  at <- match(by_key(pca$results[sub, ]), by_key(listed$results))
  expect_equal(listed$results$log_p_two_sided[at],
               pca$results$log_p_two_sided[sub], tolerance = 1e-12)
  expect_equal(listed$results$log_p_positive[at],
               pca$results$log_p_positive[sub], tolerance = 1e-12)

  low <- inlaST.test(fit, pairs = pairs, threads = 2L, rank = 3L, moments = "pcalearning")
  expect_true(all(low$pca_learning$genes$e2_relative > -1e-12))
  expect_equal(low$pca_learning$genes$e2,
               low$pca_learning$genes$fro2 - unname(rowSums(low$pca_learning$coefficients^2)))
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
  # Only the levels the saddlepoint needs; for an orthonormal basis the level-1
  # table is the identity, so t1 is the inner product of the coefficients.
  expect_false(any(c("Tsym3", "Tsym4") %in% names(tables)))
  expect_equal(tables$Tsym1, diag(r), tolerance = 1e-12)
  d2 <- as.integer(choose(r + 1L, 2L))
  expect_identical(dim(tables$Tsym2), c(d2, d2))
  C <- matrix(rnorm(6L * r), 6L, r)
  K2 <- mgcvST:::mgcvst_pca_monomials_cpp(C)
  expect_identical(dim(K2), c(d2, 6L))
  A <- matrix(rnorm(3L * 6L), 3L, 6L)
  i <- c(1L, 2L, 3L, 5L)
  j <- c(2L, 4L, 6L, 1L)
  out <- mgcvST:::mgcvst_pca_spa_pairs_cpp(A, C, K2, tables$Tsym2, matrix(1, 1L, 6L),
                                            rep(1, 6L), i, j, 2L)
  H <- lapply(seq_len(nrow(C)), function(k) unpack(B %*% C[k, ]))
  brute <- t(vapply(seq_along(i), function(k) {
    M <- H[[i[k]]] %*% H[[j[k]]]
    c(sum(diag(M)), sum(diag(M %*% M)))
  }, numeric(2L)))
  expect_equal(unname(out[, c("t1", "t2")]), brute, tolerance = 1e-5)
  expect_equal(unname(out[, "t1"]), rowSums(C[i, ] * C[j, ]), tolerance = 1e-12)
  expect_equal(out[, "U"], colSums(A[, i] * A[, j]), tolerance = 1e-12)
  # The block kernel gives the same values for every pair of a range of genes.
  block <- mgcvST:::mgcvst_pca_spa_block_cpp(A, C, K2, tables$Tsym2, matrix(1, 1L, 6L),
                                             rep(1, 6L), 1L, 5L, 2L)
  pair_index <- t(utils::combn(6L, 2L))
  pair_index <- pair_index[pair_index[, 1L] <= 5L, , drop = FALSE]
  single <- mgcvST:::mgcvst_pca_spa_pairs_cpp(A, C, K2, tables$Tsym2, matrix(1, 1L, 6L),
                                              rep(1, 6L), pair_index[, 1L],
                                              pair_index[, 2L], 1L)
  expect_equal(block, single, tolerance = 1e-12)
})

test_that("PCAlearning pair kernel reproduces the approx-liu-p rank-10 traces", {
  fx <- readRDS(test_path("fixtures", "pcalearning-approx-liu-p.rds"))
  r <- fx$rank
  d2 <- choose(r + 1L, 2L)
  T2 <- .pca_unpack(fx$tables$Tsym2, d2)
  P <- fx$pairs
  n <- nrow(fx$C)
  K2 <- mgcvST:::mgcvst_pca_monomials_cpp(fx$C)
  out <- mgcvST:::mgcvst_pca_spa_pairs_cpp(fx$A, fx$C, K2, T2, matrix(1, 1L, n),
                                            rep(1, n), P$i, P$j, 2L)
  expect_equal(out[, "U"], P$U, tolerance = 1e-10)
  approx <- as.matrix(P[c("t1_approx", "t2_approx")])
  expect_equal(unname(out[, c("t1", "t2")]), unname(approx), tolerance = 1e-6)
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
      fit, index, 2L, 1000L, FALSE, basis, rank = rank,
      n_per_cell = n_per_cell, seed = 1L, checkpoint_dir = dir, resume = resume
    )
  }
  rows <- function(x) .pca_pairs(x)
  index <- t(combn(G, 2L))
  reference <- run(index, NULL)
  expect_null(reference$metadata$pca_learning$checkpoint$path)

  dir <- tempfile("mgcvst-pca-checkpoint-")
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  # Interrupted run: only genes 1-6 are requested, so genes 2, 4, 6 are projected.
  part <- run(t(combn(6L, 2L)), dir)
  expect_identical(part$metadata$pca_learning$checkpoint$projected_genes, 3L)
  expect_setequal(list.files(dir, "^(manifest|pca-)"),
                  c("manifest.rds", "pca-basis.rds", "pca-projection-000001.rds"))
  expect_length(list.files(dir, "^pairs-"), 1L)
  same <- function(x) {
    expect_identical(rows(x), rows(reference))
    expect_identical(x$metadata$pca_learning$coefficients,
                     reference$metadata$pca_learning$coefficients)
    expect_identical(x$metadata$pca_learning$genes, reference$metadata$pca_learning$genes)
  }
  resumed <- run(index, dir)
  z <- resumed$metadata$pca_learning$checkpoint
  expect_true(z$basis_resumed)
  expect_identical(c(z$resumed_genes, z$projected_genes), c(3L, 1L))
  same(resumed)
  expect_identical(resumed$metadata$resumed_pairs, 0)

  # Completed pair shards are reused by a repeated run, and a lost block is
  # recomputed into a new block.
  repeated <- run(index, dir)
  expect_identical(repeated$metadata$resumed_pairs, 28)
  unlink(file.path(dir, "pca-projection-000001.rds"))
  again <- run(index, dir)
  expect_identical(again$metadata$pca_learning$checkpoint$projected_genes, 3L)
  expect_true(file.exists(file.path(dir, "pca-projection-000003.rds")))
  same(again)
  complete <- run(index, dir)
  expect_identical(complete$metadata$pca_learning$checkpoint$projected_genes, 0L)
  same(complete)

  expect_error(run(index, dir, rank = 2L), "different fit, score basis, route, rank")
  expect_error(run(index, dir, resume = FALSE), "already exists")

  # A manifest of an earlier algorithm contract is refused, not resumed.
  manifest <- file.path(dir, "manifest.rds")
  current <- readRDS(manifest)
  earlier <- current
  earlier$signature$version <- 2L
  earlier$signature$method <- "pca_learning"
  earlier$signature$contract <- NULL
  saveRDS(earlier, manifest)
  expect_error(run(index, dir), "earlier algorithm contract")
  saveRDS(current, manifest)
  expect_identical(rows(run(index, dir)), rows(reference))

  # Pair shards of another algorithm contract are refused as well.
  old <- file.path(dir, "pairs-0123456789")
  dir.create(old)
  saveRDS(list(), file.path(old, "block-0000000001.rds"))
  expect_error(run(index, dir), "different algorithm contract")
  unlink(old, recursive = TRUE)

  public <- tempfile("mgcvst-pca-public-")
  on.exit(unlink(public, recursive = TRUE), add = TRUE)
  test <- function() inlaST.test(fit, pairs = t(combn(fit$feature_id, 2L)),
                                 threads = 2L, rank = 3L, checkpoint_dir = public,
                                 resume = TRUE, moments = "pcalearning")
  first <- test()
  expect_true(file.exists(file.path(public, "pca-basis.rds")))
  expect_identical(readRDS(file.path(public, "route.rds"))$moments, "pcalearning")
  expect_identical(test()$results, first$results)
  expect_identical(first$checkpoint_dir, normalizePath(public, winslash = "/"))
})

test_that("PCAlearning checks rank, n_per_cell, seed, k and trace-table memory", {
  skip_on_cran()
  fit <- .pca_nb_fit()
  G <- length(fit$feature_id)
  basis <- mgcvST:::.inlast_sparse_observation_basis(mgcvST:::.inlast_sparse_prepare(fit))
  run <- function(rank, n_per_cell = 3L, seed = 1L, k = 50L) {
    index <- t(combn(G, 2L))
    mgcvST:::.mgcvst_pair_pcalearning(
      fit, index, 2L, 1000L, FALSE, basis, rank = rank,
      n_per_cell = n_per_cell, seed = seed, k = k
    )
  }
  expect_error(run(G + 1L), sprintf("achievable PCAlearning rank %d (%d training", G, G), fixed = TRUE)
  expect_error(run(0L), "rank must be one positive integer")
  expect_error(run(2.5), "rank must be one positive integer")
  expect_error(run(3L, n_per_cell = 0L), "n_per_cell must be one positive integer")
  expect_error(run(3L, seed = -1L), "seed must be one non-negative integer")
  expect_error(run(3L, seed = 1.5), "seed must be one non-negative integer")
  expect_error(run(3L, k = 0L), "k must be NULL or one positive integer")
  # q = 1404, r = 20: levels 1 and 2 only, about 3.3 GiB.
  expect_equal(mgcvST:::.mgcvst_pca_table_bytes(1404, 20) / 1024^3, 3.30, tolerance = 1e-2)
  local_mocked_bindings(.mgcvst_memory_probe = function(...) list(available = 1e3),
                        .package = "mgcvST")
  expect_error(run(3L), "trace tables for rank = 3 .* use a smaller rank")
})

test_that("a PCAlearning resume follows the stored schedule whatever chunk_size it is given", {
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
  run <- function(index, dir, chunk_size) {
    mgcvST:::.mgcvst_pair_pcalearning(
      fit, index, 2L, chunk_size, FALSE, basis, rank = 3L, n_per_cell = 3L,
      seed = 1L, checkpoint_dir = dir
    )
  }
  for (index in list(NULL, t(combn(G, 2L)))) {
    dir <- tempfile("mgcvst-pca-chunks-")
    on.exit(unlink(dir, recursive = TRUE), add = TRUE)
    first <- run(index, dir, 5L)
    expect_identical(first$metadata$resumed_pairs, 0)
    expect_gt(first$metadata$chunks, 1L)
    for (later in c(1000L, 2L)) {
      again <- run(index, dir, later)
      expect_identical(again$metadata$chunk_size, 5L)
      expect_identical(again$metadata$chunks, first$metadata$chunks)
      expect_identical(again$metadata$resumed_pairs, choose(G, 2))
      expect_identical(.pca_pairs(again), .pca_pairs(first))
    }
  }
  # The stored schedule is checked against the universe it was written for.
  schedule <- file.path(first$metadata$pair_dir, "schedule.rds")
  stored <- mgcvST:::.mgcvst_schedule_load(schedule)
  expect_identical(stored$kind, "pcalearning")
  stored$n <- stored$n + 1L
  unlink(schedule)
  mgcvST:::.mgcvst_schedule_save(schedule, stored)
  expect_error(run(index, dir, 5L), "pair schedule is damaged or incompatible")
})

test_that("PCAlearning controls fail before the basis, and the pair directory follows the basis", {
  skip_on_cran()
  fit <- .pca_nb_fit()
  calls <- 0L
  observation_basis <- mgcvST:::.inlast_sparse_observation_basis
  local_mocked_bindings(
    .inlast_sparse_observation_basis = function(fit) {
      calls <<- calls + 1L
      observation_basis(fit)
    },
    .package = "mgcvST")
  expect_error(inlaST.test(fit, rank = 0L), "rank must be one positive integer")
  expect_error(inlaST.test(fit, n_per_cell = 1.5), "n_per_cell must be one positive")
  expect_error(inlaST.test(fit, seed = -1), "seed must be one non-negative integer")
  expect_error(inlaST.test(fit, k = 2.5), "k must be NULL or one positive integer")
  expect_error(inlaST.test(fit, threads = 1.5), "threads must be one positive integer")
  expect_identical(calls, 0L)

  dir <- tempfile("mgcvst-pca-order-")
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  local_mocked_bindings(.mgcvst_pca_basis = function(...) stop("basis failure"),
                        .package = "mgcvST")
  expect_error(inlaST.test(fit, rank = 3L, checkpoint_dir = dir, moments = "pcalearning"),
               "basis failure")
  expect_length(list.files(dir, "^pairs-"), 0L)
})

test_that("the PCAlearning Gram matrix and basis are bitwise identical for 1 and 4 threads", {
  withr::local_seed(8L)
  n <- 30L
  L <- 5000L
  packed <- lapply(seq_len(n), function(j) writeBin(rnorm(L), raw(), size = 4L))
  tau <- runif(n, 0.5, 2)
  gram1 <- mgcvST:::mgcvst_pca_gram_cpp(packed, tau, 1L, 256L)
  gram4 <- mgcvST:::mgcvst_pca_gram_cpp(packed, tau, 4L, 256L)
  expect_identical(gram1, gram4)
  dense <- vapply(packed, function(x) readBin(x, "numeric", n = L, size = 4L), numeric(L))
  expect_equal(gram1, (tau %o% tau) * crossprod(dense), tolerance = 1e-12)
  rotation <- qr.Q(qr(matrix(rnorm(n * 4L), n, 4L)))
  expect_identical(mgcvST:::mgcvst_pca_basis_cpp(packed, tau, rotation, 1L, 256L),
                   mgcvST:::mgcvst_pca_basis_cpp(packed, tau, rotation, 4L, 256L))
  # The partial sums of fixed row chunks are added in chunk order: a chunk size
  # that does not divide the rows gives the same result for every thread count.
  expect_identical(mgcvST:::mgcvst_pca_gram_cpp(packed, tau, 1L, 999L),
                   mgcvST:::mgcvst_pca_gram_cpp(packed, tau, 3L, 999L))
})

test_that("PCAlearning tests are bitwise identical for 1 and 4 threads", {
  skip_on_cran()
  fit <- .pca_nb_fit()
  pairs <- t(combn(fit$feature_id, 2L))
  one <- inlaST.test(fit, pairs = pairs, threads = 1L, rank = 4L, seed = 3L,
                     moments = "pcalearning")
  four <- inlaST.test(fit, pairs = pairs, threads = 4L, rank = 4L, seed = 3L,
                      moments = "pcalearning")
  expect_identical(one$results, four$results)
  expect_identical(one$pca_learning$coefficients, four$pca_learning$coefficients)
  expect_identical(one$pca_learning$gram_values, four$pca_learning$gram_values)
  expect_identical(one$contract$basis_sha, four$contract$basis_sha)
})

test_that("a zero-matrix gene is reported as unusable and does not disturb the other genes", {
  .pair_pipeline_mocks()
  testthat::local_mocked_bindings(
    .mgcvst_pair_build_batch = function(fit, ids, threads, native) {
      lapply(ids, function(i) {
        list(a = c(i, i + 0.25), width = 2L,
             M = if (i == 3L) matrix(0, 2L, 2L) else diag(c(i + 0.5, i + 1)))
      })
    }, .package = "mgcvST")
  fit <- .pair_pipeline_fit(6L)
  z <- mgcvST:::.mgcvst_pair_pcalearning(fit, NULL, 1L, 100L, FALSE, rank = 2L)
  out <- .pca_pairs(z)
  bad <- out$i == 3L | out$j == 3L
  expect_true(all(out$status[bad] == 3L))
  expect_true(all(is.na(out$score[bad]) & is.na(out$log_p_two_sided[bad])))
  expect_true(all(out$status[!bad] == 0L))
  expect_true(all(is.finite(out$log_p_two_sided[!bad])))
  expect_identical(z$failed$feature_id, "g3")
  expect_match(z$failed$error, "zero or not finite")
  # The other pairs equal those of a run without the zero-matrix gene.
  testthat::local_mocked_bindings(
    .mgcvst_pair_build_batch = function(fit, ids, threads, native) {
      lapply(ids, function(i) list(a = c(i, i + 0.25), M = diag(c(i + 0.5, i + 1)),
                                   width = 2L))
    }, .package = "mgcvST")
  full <- .pca_pairs(mgcvST:::.mgcvst_pair_pcalearning(fit, NULL, 1L, 100L, FALSE, rank = 2L))
  keep <- !bad
  expect_equal(out$log_p_two_sided[keep], full$log_p_two_sided[full$i != 3L & full$j != 3L],
               tolerance = 1e-3)
})

test_that("an mgcv fit takes the PCAlearning route and agrees with the exact route", {
  .pair_pipeline_mocks()
  fit <- .pair_pipeline_fit(6L)
  exact <- mgcvST:::.mgcvst_pair_pipeline(fit, NULL, 1L, 100L, FALSE, k = 2L)
  pca <- mgcvST:::.mgcvst_pair_pcalearning(fit, NULL, 1L, 100L, FALSE, rank = 2L, k = 2L)
  a <- do.call(rbind, lapply(exact$shards, mgcvST:::.mgcvst_read_shard))
  b <- .pca_pairs(pca)
  a <- a[order(a$i, a$j), ]
  expect_identical(a$i, b$i)
  expect_identical(a$j, b$j)
  expect_equal(a$score, b$score, tolerance = 1e-12)
  expect_true(all(a$status == 0L & b$status == 0L))
  expect_equal(-b$log_p_two_sided / log(10), -a$log_p_two_sided / log(10), tolerance = 1e-3)
  expect_identical(pca$metadata$preparation_backend, "model_native")
  expect_identical(pca$metadata$contract$route, "pcalearning")
  expect_identical(exact$metadata$contract$route, "exact")
})
