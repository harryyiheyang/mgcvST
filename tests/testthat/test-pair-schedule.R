test_that("the pair order clusters gene blocks", {
  pairs <- t(utils::combn(8L, 2L))
  ord <- spdeST:::.mgcvst_pair_order(pairs, 1:8, 4L)
  expect_identical(sort(ord), seq_len(nrow(pairs)))
  expect_false(identical(ord, seq_len(nrow(pairs))))
  tile <- (pairs - 1L) %/% 2L
  keys <- tile[, 1L] * 4L + tile[, 2L]
  expect_true(all(diff(keys[ord]) >= 0))
  expect_identical(spdeST:::.mgcvst_pair_order(pairs, 1:8, 8L), seq_len(nrow(pairs)))
})

test_that("the all-pairs schedule visits every pair once within the block bound", {
  for (n in c(2L, 5L, 9L)) {
    for (width in unique(c(1L, 2L, 4L, n))) {
      for (chunk in c(1L, 3L, 10L, 1000L)) {
        blocks <- spdeST:::.mgcvst_all_pair_blocks(n, width, chunk)
        pairs <- do.call(rbind, lapply(seq_len(nrow(blocks)), function(b) {
          z <- spdeST:::.mgcvst_block_pairs(blocks[b, ])
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
  fixed <- spdeST:::.mgcvst_all_pair_schedule(9L, 4L, 5L)
  expect_identical(fixed$width, 2L)
  expect_identical(fixed$blocks, spdeST:::.mgcvst_all_pair_blocks(9L, 2L, 5L))
  expect_identical(spdeST:::.mgcvst_all_pair_schedule(9L, 100L, 5L)$width, 9L)
})

test_that("the pair test writes no checkpoint: no such argument, and nothing is left in the temporary directory", {
  skip_on_cran()
  expect_false(any(c("checkpoint_dir", "resume") %in% names(formals(spdeST.test))))
  for (name in c(".mgcvst_route_save", ".mgcvst_route_check", ".mgcvst_shard_complete",
                 ".mgcvst_shard_window", ".mgcvst_pairs_refuse_stale",
                 ".mgcvst_schedule_save", ".mgcvst_pair_signature", ".mgcvst_pca_save")) {
    expect_false(exists(name, asNamespace("spdeST"), inherits = FALSE), info = name)
  }
  fit <- .pca_nb_fit()
  old <- setwd(withr::local_tempdir())
  before <- list.files(tempdir(), all.files = TRUE, no.. = TRUE)
  on.exit(setwd(old), add = TRUE)
  for (moments in c("exact", "pcalearning")) for (keep in c("all", "discoveries")) {
    test <- if (identical(moments, "exact")) {
      spdeST.test(fit, chunk_size = 2L, moments = moments, keep = keep)
    } else spdeST.test(fit, chunk_size = 2L, moments = moments, keep = keep, rank = 2L, k = 3L)
    expect_length(test$shards, 0L)
    expect_null(test$checkpoint_dir)
    expect_identical(setdiff(list.files(tempdir(), all.files = TRUE, no.. = TRUE), before),
                     character(), info = paste(moments, keep))
    expect_length(list.files(".", all.files = TRUE, no.. = TRUE), 0L)
  }
})
