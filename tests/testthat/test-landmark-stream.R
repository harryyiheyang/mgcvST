test_that("native landmark queue preserves states and four trace moments", {
  path <- tempfile("mgcvst-landmark-stream-")
  dir.create(path)
  on.exit(unlink(path, recursive = TRUE), add = TRUE)
  M <- list(
    diag(c(1, 2, 3)),
    crossprod(matrix(c(1, 1, 0, 0, 2, 1, 1, 0, 2), 3L)) + diag(3) / 2
  )
  references <- list(
    tcrossprod(matrix(c(1, 0, 2, 1, 2, 1, 0, 1, 1), 3L)) + diag(3) / 3,
    diag(c(2, 1, 4))
  )
  expect_gt(max(abs(M[[1L]] %*% references[[1L]] -
                    references[[1L]] %*% M[[1L]])), 0)
  scores <- list(c(1, 2, 3), c(2, -1, 1))
  state <- file.path(path, paste0("state-", seq_len(3L), ".bin"))
  out1 <- file.path(path, paste0("one-", seq_len(3L), ".bin"))
  out2 <- file.path(path, paste0("two-", seq_len(3L), ".bin"))
  for (i in seq_along(M)) {
    mgcvST:::mgcvst_state_write_cpp(
      state[i], "state-signature", paste0("g", i), scores[[i]], M[[i]],
      c(spatial = 3L)
    )
  }
  mgcvST:::mgcvst_state_write_cpp(
    state[3L], "state-signature", "g3", numeric(),
    matrix(numeric(), 0L, 0L), integer(), "fit failed"
  )
  restored <- mgcvST:::mgcvst_state_read_cpp(state[1L], "state-signature", "g1")
  expect_equal(restored$M, M[[1L]])
  expect_equal(restored$a, scores[[1L]])
  expect_equal(restored$width, c(spatial = 3))
  expect_null(restored$error)
  expect_identical(mgcvST:::mgcvst_state_read_cpp(
    state[3L], "state-signature", "g3"), list(error = "fit failed"))

  ids <- paste0("g", seq_len(3L))
  one <- mgcvST:::mgcvst_landmark_stream_cpp(
    state, ids, "state-signature", references, out1, "summary-signature", 1L
  )
  two <- mgcvST:::mgcvst_landmark_stream_cpp(
    state, ids, "state-signature", references, out2, "summary-signature", 2L
  )
  expect_identical(one$completed, 3L)
  expect_identical(two$completed, 3L)
  expect_identical(two$source_errors, 1L)
  expect_identical(two$source_error_message[3L], "fit failed")
  expect_gt(two$workspace_bytes_estimate, 0)
  for (i in seq_along(M)) {
    first <- mgcvST:::mgcvst_trace_read_cpp(
      out1[i], "summary-signature", ids[i], length(references)
    )
    second <- mgcvST:::mgcvst_trace_read_cpp(
      out2[i], "summary-signature", ids[i], length(references)
    )
    expect_identical(first, second)
    expect_identical(dim(first$cross), c(2L, 4L))
    expect_identical(first$self, NULL)
    expect_equal(first$a, scores[[i]])
    expected <- mgcvST:::mgcvst_pair_trace_powers_cpp(
      c(list(M[[i]]), references), cbind(rep(1L, 2L), 2:3), 4L, 1L
    )
    expect_equal(first$cross, expected, tolerance = 1e-5)
  }
  expect_identical(mgcvST:::mgcvst_trace_read_cpp(
    out1[3L], "summary-signature", "g3", 2L), list(error = "fit failed"))
})

test_that("native shards reject mismatches, corruption, and duplicate output", {
  path <- tempfile("mgcvst-landmark-integrity-")
  dir.create(path)
  on.exit(unlink(path, recursive = TRUE), add = TRUE)
  M <- matrix(c(1, 0.2, 0.2, 2), 2L)
  state <- file.path(path, "state.bin")
  summary <- file.path(path, "summary.bin")
  mgcvST:::mgcvst_state_write_cpp(state, "state-sig", "g", c(1, 2), M, 2L)
  expect_error(mgcvST:::mgcvst_landmark_stream_cpp(
    NA_character_, "g", "state-sig", list(M), summary, "summary-sig", 1L),
    "cannot be NA")
  expect_error(mgcvST:::mgcvst_state_write_cpp(
    state, "state-sig", "g", c(1, 2), M, 2L), "already exists")
  expect_error(mgcvST:::mgcvst_state_read_cpp(state, "wrong", "g"), "signature")
  expect_error(mgcvST:::mgcvst_state_read_cpp(state, "state-sig", "wrong"),
               "feature ID")
  mgcvST:::mgcvst_landmark_stream_cpp(
    state, "g", "state-sig", list(M), summary, "summary-sig", 1L
  )
  expect_error(mgcvST:::mgcvst_trace_read_cpp(summary, "wrong", "g", 1L),
               "signature")
  expect_error(mgcvST:::mgcvst_trace_read_cpp(summary, "summary-sig", "wrong", 1L),
               "feature ID")
  expect_error(mgcvST:::mgcvst_trace_read_cpp(summary, "summary-sig", "g", 2L),
               "dimensions")
  expect_error(mgcvST:::mgcvst_landmark_stream_cpp(
    state, "g", "state-sig", list(M), summary, "summary-sig", 1L),
    "already exists")
  expect_error(mgcvST:::mgcvst_landmark_stream_cpp(
    rep(state, 2L), rep("g", 2L), "state-sig", list(M),
    rep(file.path(path, "duplicate.bin"), 2L), "summary-sig", 2L),
    "distinct")

  bytes <- readBin(state, "raw", n = file.info(state)$size)
  damaged <- file.path(path, "damaged.bin")
  bytes[73L] <- as.raw(bitwXor(as.integer(bytes[73L]), 1L))
  writeBin(bytes, damaged)
  expect_error(mgcvST:::mgcvst_state_read_cpp(damaged, "state-sig", "g"),
               "checksum")
  truncated <- file.path(path, "truncated.bin")
  writeBin(bytes[seq_len(70L)], truncated)
  expect_error(mgcvST:::mgcvst_state_read_cpp(truncated, "state-sig", "g"),
               "short|Truncated|length")
  expect_error(mgcvST:::mgcvst_landmark_stream_cpp(
    damaged, "g", "state-sig", list(M),
    file.path(path, "unwritten.bin"), "summary-sig", 1L),
    "completed 0 of 1.*checksum")
  expect_false(file.exists(file.path(path, "unwritten.bin")))
})
