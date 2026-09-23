landmark_sparse_fixture <- function() {
  list(
    feature_id = paste0("g", seq_len(6L)),
    score_backend = "sparse",
    dispersion = c(2, 4, 6, 8, 10, 12),
    smoothing_parameters = cbind(spatial = c(1, 2, 3, 4, 5, 6),
                                 slide = c(2, 2, 3, 4, 5, 6)),
    model = list(inla_spec = list(random = list(
      list(name = "global", target = TRUE, kind = "spde", sp_index = 1L),
      list(name = "slide", target = FALSE, kind = "nuisance",
           subtype = "iid", sp_index = 2L)
    )))
  )
}

test_that("random landmark selection is reproducible and restores RNG", {
  fit <- landmark_sparse_fixture()
  set.seed(710L)
  before <- .Random.seed
  a <- mgcvST:::.mgcvst_landmark_select(fit, c(6L, 2L, 4L, 1L), 2L,
                                        "random", seed = 19L)
  expect_identical(.Random.seed, before)
  b <- mgcvST:::.mgcvst_landmark_select(fit, c(6L, 2L, 4L, 1L), 2L,
                                        "random", seed = 19L)
  expect_identical(a, b)
  expect_length(a, 2L)
  expect_true(all(a %in% c(1L, 2L, 4L, 6L)))
  expect_identical(mgcvST:::.mgcvst_landmark_select(
    fit, 1:4, 10L, "random", seed = 19L
  ), 1:4)
})

test_that("score landmarks are real nearest cluster members with aligned columns", {
  fit <- landmark_sparse_fixture()
  scores <- rbind(c(0, 0.1, 0.2, 9, 9.1, 9.2),
                  c(0.2, 0.1, 0, 9.2, 9.1, 9))
  colnames(scores) <- fit$feature_id
  selected <- mgcvST:::.mgcvst_landmark_select(
    fit, 1:6, 2L, "score", seed = 31L, scores = scores[, 6:1]
  )
  expect_length(selected, 2L)
  expect_true(all(selected %in% 1:6))
  expect_equal(sum(selected %in% 1:3), 1L)
  expect_equal(sum(selected %in% 4:6), 1L)

  x <- t(scores)
  set.seed(31L)
  km <- stats::kmeans(x, centers = 2L, iter.max = 100L, nstart = 5L)
  expected <- vapply(seq_len(2L), function(k) {
    members <- which(km$cluster == k)
    d <- sweep(x[members, , drop = FALSE], 2L, km$centers[k, ], "-")
    members[which.min(rowSums(d * d))]
  }, integer(1L))
  expect_identical(selected, sort(expected))
})

test_that("hyper features use only mapped variance scales", {
  fit <- landmark_sparse_fixture()
  h <- mgcvST:::.mgcvst_landmark_hyper(fit, c(3L, 1L, 5L))
  expect_identical(colnames(h), c("variance_spatial", "variance_slide"))
  expect_identical(rownames(h), c("g3", "g1", "g5"))
  expect_equal(unname(h[, 1L]), fit$dispersion[c(3L, 1L, 5L)] /
                 fit$smoothing_parameters[c(3L, 1L, 5L), 1L])
  expect_equal(unname(h[, 2L]), fit$dispersion[c(3L, 1L, 5L)] /
                 fit$smoothing_parameters[c(3L, 1L, 5L), 2L])
  expect_false(any(grepl("range|nugget|theta", colnames(h))))
  hyper_selected <- mgcvST:::.mgcvst_landmark_select(
    fit, 1:6, 2L, "hyper", seed = 17L
  )
  expect_length(hyper_selected, 2L)

  dense <- list(
    feature_id = fit$feature_id,
    dispersion = fit$dispersion,
    smoothing_parameters = cbind(penalty1 = 1:6, penalty2 = 6:1),
    geometry = list(smooth = list(
      list(label = "s(space)", sp_index = 1L),
      list(label = "s(slide)", sp_index = 2L)
    ))
  )
  hd <- mgcvST:::.mgcvst_landmark_hyper(dense, c(2L, 4L))
  expect_identical(colnames(hd), c("variance_s.space.", "variance_s.slide."))
  expect_equal(unname(hd[, 1L]), dense$dispersion[c(2L, 4L)] /
                 dense$smoothing_parameters[c(2L, 4L), 1L])

  dense$model <- "legacy-model-object"
  expect_equal(mgcvST:::.mgcvst_landmark_hyper(dense, 2L)[1L, 1L], 2)
})

test_that("legacy SPDE hyper features use the existing spatial variance scale", {
  fit <- list(
    feature_id = c("g1", "g2", "g3"),
    test_engine = "spde",
    dispersion = c(2, 4, 6),
    lambda = c(4, 2, 3),
    working_error = matrix(0, 4L, 3L),
    working_variance = matrix(1, 4L, 3L)
  )
  h <- mgcvST:::.mgcvst_landmark_hyper(fit, c(3L, 1L))
  expect_identical(colnames(h), "variance_spatial")
  expect_identical(rownames(h), c("g3", "g1"))
  expect_equal(unname(h[, 1L]), c(2, 0.5))
  expect_error(mgcvST:::.mgcvst_landmark_hyper(
    within(fit, lambda[1L] <- 0), 1:3
  ), "positive")
})

test_that("landmark selectors reject unaligned or non-finite inputs", {
  fit <- landmark_sparse_fixture()
  scores <- matrix(1:12, nrow = 2L)
  colnames(scores) <- paste0("x", 1:6)
  expect_error(mgcvST:::.mgcvst_landmark_select(
    fit, 1:6, 2L, "score", seed = 1L, scores = scores
  ), "aligned")
  colnames(scores) <- fit$feature_id
  scores[1L, 1L] <- NA_real_
  expect_error(mgcvST:::.mgcvst_landmark_select(
    fit, 1:6, 2L, "score", seed = 1L, scores = scores
  ), "finite")
  expect_error(mgcvST:::.mgcvst_landmark_select(
    fit, 1:6, 2L, "score", seed = 1L,
    scores = matrix(1, 2L, 6L, dimnames = list(NULL, fit$feature_id))
  ), "constant")
  bad <- fit
  bad$smoothing_parameters[1L, 1L] <- 0
  expect_error(mgcvST:::.mgcvst_landmark_hyper(bad, 1:3), "positive")
  expect_error(mgcvST:::.mgcvst_landmark_select(
    fit, 1:6, 0L, "random", seed = 1L
  ), "positive integer")
})
