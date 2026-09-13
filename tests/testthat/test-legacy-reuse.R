test_that("legacy Davies builds each score state once per call", {
  skip_if_not_installed("CompQuadForm")
  f <- st_fixture()
  fit <- mgcvST.estimate(
    f$Y, f$G, diagnostics = FALSE,
    BPPARAM = BiocParallel::SerialParam()
  )
  pairs <- rbind(
    c(1L, 2L), c(1L, 3L), c(2L, 3L)
  )
  scale <- mgcvST:::.mgcvst_field_scale(fit)
  T0 <- mgcvST:::.mgcvst_legacy_shared_score_factor(fit$geometry)
  state <- lapply(seq_along(fit$feature_id), function(i) {
    op <- mgcvST:::.rkhs_score_operator_factor(
      sqrt(scale[i]) * T0, fit$working_variance[, i], fit$geometry$X,
      field_scale = scale[i]
    )
    mgcvST:::rkhs_score_summary(fit$working_error[, i], op)
  })
  reference <- lapply(seq_len(nrow(pairs)), function(i) {
    s1 <- state[[pairs[i, 1L]]]
    s2 <- state[[pairs[i, 2L]]]
    score <- as.numeric(crossprod(s1$a, s2$a))
    cal <- mgcvST:::rkhs_score_calibrate(
      score, s1$H, s2$H, method = "davies"
    )
    c(score = score, information = cal$information,
      effective_rank = cal$effective_rank, p_two_sided = cal$p_two_sided)
  })
  reference <- do.call(rbind, reference)

  original_factor <- mgcvST:::.mgcvst_legacy_shared_score_factor
  original_summary <- mgcvST:::rkhs_score_summary
  count <- new.env(parent = emptyenv())
  count$factor <- 0L
  count$summary <- 0L
  testthat::local_mocked_bindings(
    .mgcvst_legacy_shared_score_factor = function(geometry) {
      count$factor <- count$factor + 1L
      original_factor(geometry)
    },
    rkhs_score_summary = function(error, operator) {
      count$summary <- count$summary + 1L
      original_summary(error, operator)
    },
    .package = "mgcvST"
  )

  observed <- mgcvST.test(
    fit, pairs = pairs, calibration = "davies", chunk_size = 1L,
    BPPARAM = BiocParallel::SerialParam()
  )
  expect_identical(count$factor, 1L)
  expect_identical(count$summary, length(unique(as.vector(pairs))))
  expect_equal(observed$results$signed_score, reference[, "score"],
               tolerance = 1e-10)
  expect_equal(observed$results$information, reference[, "information"],
               tolerance = 1e-10)
  expect_equal(observed$results$effective_rank, reference[, "effective_rank"],
               tolerance = 1e-10)
  expect_equal(observed$results$p_two_sided, reference[, "p_two_sided"],
               tolerance = 1e-10)
})
