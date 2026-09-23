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
})

test_that("pair checkpoints require complete compatible result batches", {
  path <- tempfile("pair-checkpoint-")
  dir.create(path)
  on.exit(unlink(path, recursive = TRUE), add = TRUE)
  z <- data.frame(pair_index = 1:2, score = c(1, 2), information = 1,
    effective_rank = 1, p_value = c(0.4, 0.2), error_message = NA_character_)
  mgcvST:::.mgcvst_pair_checkpoint_write(path, 1L, 2L, z)
  expect_identical(mgcvST:::.mgcvst_pair_checkpoint_read(path, 1L, 1:3)$result, z)
  expect_identical(mgcvST:::.mgcvst_pair_checkpoint_read(
    path, 1L, stats::setNames(1:3, c("a", "a", "b")))$result, z)
  expect_null(mgcvST:::.mgcvst_pair_checkpoint_read(path, 3L, 1:3))
  expect_error(mgcvST:::.mgcvst_pair_checkpoint_read(path, 1L, 3:1), "damaged")
  file <- file.path(path, "block-0000000001.rds")
  corrupted <- readRDS(file)
  corrupted$result$p_value[1L] <- 0.01
  saveRDS(corrupted, file)
  expect_error(mgcvST:::.mgcvst_pair_checkpoint_read(path, 1L, 1:3), "damaged")
})
