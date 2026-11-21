.pca_nb_fit <- local({
  cached <- NULL
  function() {
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
    model <- spdeST.set(response ~ z + offset(offset0), data, basis,
                        family = mgcv::nb())
    cached <<- spdeST:::.spdest_estimate(
    .spdest_scan(Y, model, BPPARAM = BiocParallel::SerialParam()),
    Y, features = "all", BPPARAM = BiocParallel::SerialParam())
    cached
  }
})

.pca_pairs <- function(z) {
  out <- do.call(rbind, lapply(z$shards, spdeST:::.mgcvst_read_shard))
  out[order(out$i, out$j), , drop = FALSE]
}

# Full-spectrum saddlepoint log p-values in the reduced observation-kernel
# coordinates: the reference that a full-rank PCAlearning basis and the exact
# route with k = q reproduce.
.pca_exact_reference <- function(fit, pairs) {
  prepared <- spdeST:::.spdest_sparse_prepare(fit)
  basis <- spdeST:::.spdest_sparse_observation_basis(prepared)
  units <- spdeST:::.spdest_sparse_units(prepared, seq_along(fit$feature_id), threads = 1L)
  states <- spdeST:::.spdest_sparse_materialize_reduced(prepared, units, basis, threads = 1L)
  local <- matrix(match(pairs, fit$feature_id), ncol = 2L)
  roots <- lapply(states, function(z) .spa_sqrtm(z$M))
  score <- vapply(seq_len(nrow(local)), function(j) {
    sum(states[[local[j, 1L]]]$a * states[[local[j, 2L]]]$a)
  }, numeric(1L))
  log_p <- vapply(seq_len(nrow(local)), function(j) {
    s <- svd(roots[[local[j, 1L]]] %*% roots[[local[j, 2L]]])$d
    .spa_ref(abs(score[j]), s[s > 1e-13 * s[1L]])
  }, numeric(1L))
  list(i = local[, 1L], j = local[, 2L], score = score, log_p = log_p,
       q = basis$rank)
}

.pair_pipeline_fit <- function(p = 4L) {
  ids <- paste0("g", seq_len(p))
  list(
    feature_id = ids, score_a = matrix(0, 2L, p),
    dispersion = rep(1, p), lambda = rep(1, p),
    smoothing_parameters = matrix(1, p, 1L),
    score_sparse = list(Q = diag(2)),
    geometry = list(target = c(global = 1L))
  )
}

# The observation basis of the fixture: two score coordinates.
.pair_pipeline_basis <- list(rank = 2L)

# The state of feature i of the fixture.
.pair_pipeline_state <- function(i) list(a = c(i, i + 0.25), M = diag(c(i + 0.5, i + 1)))

# The pipeline builds states through these three steps; the fixture replaces
# the sparse score kernel with deterministic states. `materialize` maps the
# feature indices to their states.
.pair_pipeline_mocks <- function(env = parent.frame(),
                                 materialize = function(ids) lapply(ids, .pair_pipeline_state)) {
  testthat::local_mocked_bindings(
    .spdest_sparse_prepare = function(fit) fit,
    .spdest_sparse_units = function(fit, features, threads = 1L) {
      lapply(features, function(i) list(id = i))
    },
    .spdest_sparse_materialize_reduced = function(fit, units, basis, threads = 1L) {
      materialize(vapply(units, `[[`, integer(1L), "id"))
    },
    .package = "spdeST", .env = env
  )
}
