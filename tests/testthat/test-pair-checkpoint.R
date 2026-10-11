test_that("pair schedules cluster gene blocks and survive budget changes", {
  pairs <- t(utils::combn(8L, 2L))
  path <- tempfile("pair-checkpoint-")
  dir.create(path)
  on.exit(unlink(path, recursive = TRUE), add = TRUE)
  ord <- mgcvST:::.mgcvst_pair_order(pairs, 1:8, 4L, path)
  expect_identical(sort(ord), seq_len(nrow(pairs)))
  expect_false(identical(ord, seq_len(nrow(pairs))))
  expect_identical(mgcvST:::.mgcvst_pair_order(pairs, 1:8, 2L, path), ord)
  tile <- (pairs - 1L) %/% 2L
  keys <- tile[, 1L] * 4L + tile[, 2L]
  expect_true(all(diff(keys[ord]) >= 0))
  file <- file.path(path, "schedule.rds")
  z <- readRDS(file)
  z$order[1:2] <- z$order[2:1]
  saveRDS(z, file)
  expect_error(mgcvST:::.mgcvst_pair_order(pairs, 1:8, 4L, path), "damaged")
})

test_that("the all-pairs schedule visits every pair once within the block bound", {
  for (n in c(2L, 5L, 9L)) {
    for (width in unique(c(1L, 2L, 4L, n))) {
      for (chunk in c(1L, 3L, 10L, 1000L)) {
        blocks <- mgcvST:::.mgcvst_all_pair_blocks(n, width, chunk)
        pairs <- do.call(rbind, lapply(seq_len(nrow(blocks)), function(b) {
          z <- mgcvST:::.mgcvst_block_pairs(blocks[b, ])
          expect_identical(nrow(z), as.integer(blocks[b, "count"]))
          # A block exceeds the bound only through a single left gene.
          expect_true(nrow(z) <= chunk || blocks[b, "l1"] == blocks[b, "l2"])
          z
        }))
        expect_equal(nrow(pairs), n * (n - 1L) / 2L)
        expect_true(all(pairs[, 1L] < pairs[, 2L]))
        expect_false(anyDuplicated(pairs) > 0L)
        expect_true(all(pairs >= 1L & pairs <= n))
      }
    }
  }
  path <- tempfile("pair-schedule-")
  dir.create(path)
  on.exit(unlink(path, recursive = TRUE), add = TRUE)
  fixed <- mgcvST:::.mgcvst_all_pair_schedule(path, 9L, 4L, 5L)
  expect_identical(fixed$width, 2L)
  again <- mgcvST:::.mgcvst_all_pair_schedule(path, 9L, 100L, 1L)
  expect_identical(again$blocks, fixed$blocks)
  expect_error(mgcvST:::.mgcvst_all_pair_schedule(path, 8L, 4L, 5L), "damaged")
})

test_that("pair shards must match the pairs of their schedule", {
  root <- tempfile("pair-shards-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  contract <- mgcvST:::.mgcvst_contract("exact")
  dir <- mgcvST:::.mgcvst_pairs_open(root, list(index = matrix(1:4, 2L)), contract)
  frame <- mgcvST:::.mgcvst_pairs_frame(c(1L, 1L), c(2L, 3L), score = c(1, -1),
    log_p_two_sided = log(c(0.4, 0.2)), log_p_positive = log(c(0.2, 0.9)),
    log_p_negative = log(c(0.8, 0.1)))
  expect_false(mgcvST:::.mgcvst_shard_complete(dir, 1L, c(1L, 1L), c(2L, 3L)))
  mgcvST:::.mgcvst_write_parquet(frame, mgcvST:::.mgcvst_shard_file(dir, 1L))
  expect_true(mgcvST:::.mgcvst_shard_complete(dir, 1L, c(1L, 1L), c(2L, 3L)))
  expect_identical(mgcvST:::.mgcvst_read_shard(
    mgcvST:::.mgcvst_shard_file(dir, 1L)), frame)
  expect_error(mgcvST:::.mgcvst_shard_complete(dir, 1L, c(1L, 2L), c(2L, 3L)),
               "damaged or does not match")
  index <- cbind(c(1L, 1L, 2L), c(2L, 3L, 3L))
  expect_identical(mgcvST:::.mgcvst_shard_window(dir, 1L, index), 2L)
  expect_true(is.na(mgcvST:::.mgcvst_shard_window(dir, 3L, index)))
  expect_error(mgcvST:::.mgcvst_shard_window(dir, 1L, index[c(2L, 1L, 3L), ]),
               "damaged or does not match")
})

test_that("pair results written under another algorithm contract are refused", {
  root <- tempfile("pair-contract-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  contract <- mgcvST:::.mgcvst_contract("exact")
  universe <- list(index = matrix(1:4, 2L))
  dir <- mgcvST:::.mgcvst_pairs_open(root, universe, contract)
  expect_identical(mgcvST:::.mgcvst_pairs_open(root, universe, contract), dir)
  expect_error(mgcvST:::.mgcvst_pairs_open(root, universe, contract, resume = FALSE),
               "already exist")
  # The rank and the basis sha are results of a run, not part of the algorithm:
  # a directory of the same route and kernel is kept whatever they are.
  expect_identical(mgcvST:::.mgcvst_pairs_open(root, universe,
    mgcvST:::.mgcvst_contract("exact", k = 20L, basis_sha = strrep("a", 64L))) != dir, TRUE)
  expect_silent(mgcvST:::.mgcvst_pairs_refuse_stale(root,
    mgcvST:::.mgcvst_contract("exact", k = 3L, basis_sha = strrep("b", 64L))))

  # Another route, kernel version, remainder order, schema or contract string
  # is refused, so that no directory of an earlier algorithm stays silently.
  expect_identical(contract$kernel_version, 3L)
  for (field in c("calibration_contract", "route", "kernel_version", "remainder_order",
                  "schema")) {
    later <- contract
    later[[field]] <- if (is.character(later[[field]])) paste0(later[[field]], "_x") else
      later[[field]] + 1L
    expect_error(mgcvST:::.mgcvst_pairs_open(root, universe, later),
                 "different algorithm contract", info = field)
  }
  expect_error(mgcvST:::.mgcvst_pairs_open(root, universe,
    mgcvST:::.mgcvst_contract("pcalearning")), "different algorithm contract")
  # Directories written by the earlier kernel versions (spa_v1, kernels 1 and 2:
  # kernel 2 had no degenerate-gene rule, so its p-values are not reusable).
  for (kernel in 1:2) {
    stale <- file.path(root, paste0("pairs-kernel", kernel))
    dir.create(stale)
    old_contract <- contract
    old_contract$kernel_version <- kernel
    saveRDS(list(version = 3L, contract = old_contract, universe_sha = "x"),
            file.path(stale, "contract.rds"))
    expect_error(mgcvST:::.mgcvst_pairs_open(root, universe, contract),
                 sprintf("kernel %d; this call writes spa_v1, exact, kernel 3", kernel))
    unlink(stale, recursive = TRUE)
  }
  old <- file.path(root, "pairs-0123")
  dir.create(old)
  saveRDS(list(first = 1L, last = 1L), file.path(old, "block-0000000001.rds"))
  expect_error(mgcvST:::.mgcvst_pairs_open(root, universe, contract),
               "none recorded")
  unlink(old, recursive = TRUE)
  expect_identical(mgcvST:::.mgcvst_pairs_open(root, universe, contract), dir)
})

test_that("checkpointed pair tests resume and refuse a pre-contract directory", {
  skip_on_cran()
  f <- st_fixture(nuisance = TRUE)
  fit <- mgcvST.estimate(f$Y, f$model, diagnostics = FALSE,
                         BPPARAM = BiocParallel::SerialParam(), spatial = "all")
  dir <- tempfile("mgcvst-checkpoint-")
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  first <- mgcvST.test(fit, checkpoint_dir = dir, chunk_size = 2L, moments = "exact")
  expect_identical(first$timing$pair_pipeline$resumed_pairs, 0)
  expect_true(file.exists(file.path(dir, "manifest.rds")))
  again <- mgcvST.test(fit, checkpoint_dir = dir, chunk_size = 2L, moments = "exact")
  expect_identical(again$timing$pair_pipeline$resumed_pairs, 3)
  expect_identical(again$results, first$results)
  expect_identical(again$timing$pair_pipeline$builds, 0L)
  expect_error(mgcvST.test(fit, checkpoint_dir = dir, resume = FALSE, moments = "exact"), "already exists")

  old <- file.path(dir, paste0("pairs-", strrep("0", 64L)))
  dir.create(old)
  saveRDS(list(first = 1L, last = 3L, result = data.frame()),
          file.path(old, "block-0000000001.rds"))
  expect_error(mgcvST.test(fit, checkpoint_dir = dir, moments = "exact"),
               "different algorithm contract")
  unlink(old, recursive = TRUE)
  reused <- mgcvST.test(fit, checkpoint_dir = dir, chunk_size = 2L, moments = "exact")
  expect_identical(reused$timing$pair_pipeline$builds, 0L)
  expect_identical(reused$results, first$results)
})
