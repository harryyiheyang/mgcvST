# Sample pair rows independently of their scores, without changing caller RNG.
.mgcvst_landmark_holdout <- function(index, references, n, seed) {
  if (!is.matrix(index) || ncol(index) != 2L ||
      !is.numeric(index) || anyNA(index) ||
      !is.numeric(n) || length(n) != 1L || !is.finite(n) ||
      n < 0 || n != floor(n) ||
      !is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
      seed < 0 || seed != floor(seed)) {
    stop("Holdout pairs, count, and seed must be valid integers.")
  }
  candidates <- which(!(index[, 1L] %in% references) &
                        !(index[, 2L] %in% references))
  if (!length(candidates) || n == 0L) return(integer())
  old_exists <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (old_exists) old <- get(".Random.seed", envir = .GlobalEnv)
  on.exit(if (old_exists) assign(".Random.seed", old, envir = .GlobalEnv) else
    if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv), add = TRUE)
  set.seed(seed)
  sort(candidates[sample.int(length(candidates), min(n, length(candidates)))])
}

# Central chi-square tails must omit ncp = 0 to retain finite log tails.
.mgcvst_liu_logp <- function(liu) {
  if (!is.list(liu)) stop("liu must be a Liu calibration result.")
  n <- length(liu$transformed)
  if (length(liu$df) != n || length(liu$ncp) != n) {
    stop("liu must contain aligned transformed, df, and ncp values.")
  }
  log_tail <- rep(NA_real_, n)
  valid <- is.finite(liu$transformed) & is.finite(liu$df) & liu$df > 0 &
    is.finite(liu$ncp) & liu$ncp >= 0
  central <- valid & liu$ncp == 0
  if (any(central)) {
    log_tail[central] <- stats::pchisq(
      liu$transformed[central], df = liu$df[central],
      lower.tail = FALSE, log.p = TRUE
    )
  }
  noncentral <- which(valid & liu$ncp > 0)
  if (length(noncentral)) {
    log_tail[noncentral] <- stats::pchisq(
      liu$transformed[noncentral], df = liu$df[noncentral],
      ncp = liu$ncp[noncentral], lower.tail = FALSE, log.p = TRUE
    )
  }
  log_tail[which(is.nan(log_tail) | log_tail == Inf | log_tail > 0)] <- NA_real_
  log_tail / log(10)
}

# Read already constructed double matrices; exact checks never refit a gene.
.mgcvst_landmark_exact <- function(store, index, threads, chunk_size,
                                    cache_bytes) {
  if (!inherits(store, "mgcvst_score_store") || store$storage != "double") {
    stop("Landmark checks require original double score states.")
  }
  if (!is.matrix(index) || ncol(index) != 2L || !is.numeric(index) ||
      anyNA(index) || any(index < 1L) ||
      any(index > length(store$feature_ids)) || any(index != floor(index))) {
    stop("index must contain two valid feature indices per row.")
  }
  if (!is.numeric(threads) || length(threads) != 1L ||
      !is.finite(threads) || threads < 1L || threads != floor(threads) ||
      !is.numeric(chunk_size) || length(chunk_size) != 1L ||
      !is.finite(chunk_size) || chunk_size < 1L ||
      chunk_size != floor(chunk_size) ||
      !is.numeric(cache_bytes) || length(cache_bytes) != 1L ||
      !is.finite(cache_bytes) || cache_bytes <= 0) {
    stop("threads, chunk_size, and cache_bytes must be positive finite values.")
  }
  rows <- nrow(index)
  moments <- matrix(NA_real_, rows, 4L)
  score <- rep(NA_real_, rows)
  errors <- rep(NA_character_, rows)
  if (!rows) return(list(score = score, moments = moments,
                         p_value = numeric(), log10_p = numeric(),
                         error_message = errors, elapsed = 0,
                         builds = 0L, chunks = 0L,
                         cache_peak_bytes = 0, store_reads = 0L))
  order_rows <- order(index[, 1L], index[, 2L], seq_len(rows))
  ordered <- index[order_rows, , drop = FALSE]
  cache <- new.env(parent = emptyenv())
  cache$states <- list()
  cache$size <- numeric()
  cache$last <- numeric()
  cache$clock <- 0
  cache$bytes <- 0
  cache$peak <- 0
  cache$reads <- 0L
  elapsed <- proc.time()[["elapsed"]]
  first <- 1L
  chunks <- 0L
  while (first <= rows) {
    last <- first - 1L
    active <- integer()
    active_bytes <- 0
    while (last < rows && last - first + 1L < chunk_size) {
      pair <- ordered[last + 1L, ]
      added <- setdiff(unique(pair), active)
      loaded <- vector("list", length(added))
      sizes <- numeric(length(added))
      for (k in seq_along(added)) {
        key <- as.character(added[k])
        z <- cache$states[[key]]
        if (is.null(z)) {
          z <- .mgcvst_store_read(store, added[k])
          cache$reads <- cache$reads + 1L
        }
        loaded[[k]] <- z
        sizes[k] <- if (key %in% names(cache$size)) cache$size[[key]] else
          as.numeric(object.size(z))
      }
      if (active_bytes + sum(sizes) > cache_bytes) {
        if (last < first) {
          stop("cache_bytes cannot hold both double score states for one pair.")
        }
        rm(loaded, sizes)
        if (exists("z", inherits = FALSE)) rm(z)
        break
      }
      pinned <- as.character(c(active, added))
      for (k in seq_along(added)) {
        key <- as.character(added[k])
        if (is.null(cache$states[[key]])) {
          while (cache$bytes + sizes[k] > cache_bytes) {
            stale <- setdiff(names(cache$last), pinned)
            if (!length(stale)) {
              stop("cache_bytes cannot hold the selected double score states.")
            }
            drop <- stale[which.min(cache$last[stale])]
            cache$bytes <- cache$bytes - cache$size[[drop]]
            cache$states[[drop]] <- NULL
            cache$size <- cache$size[names(cache$size) != drop]
            cache$last <- cache$last[names(cache$last) != drop]
          }
          cache$states[[key]] <- loaded[[k]]
          cache$size[key] <- sizes[k]
          cache$bytes <- cache$bytes + sizes[k]
          cache$peak <- max(cache$peak, cache$bytes)
        }
      }
      active <- c(active, added)
      active_bytes <- active_bytes + sum(sizes)
      cache$clock <- cache$clock + 1
      cache$last[as.character(unique(pair))] <- cache$clock
      last <- last + 1L
      rm(loaded, sizes)
      if (exists("z", inherits = FALSE)) rm(z)
    }
    selected <- order_rows[seq.int(first, last)]
    state <- lapply(active, function(id) cache$states[[as.character(id)]])
    failed <- vapply(state, function(z) !is.null(z$error), logical(1L))
    local <- matrix(match(index[selected, ], active), ncol = 2L)
    good <- !failed[local[, 1L]] & !failed[local[, 2L]]
    for (k in which(!good)) {
      ids <- local[k, ][failed[local[k, ]]]
      errors[selected[k]] <- paste(vapply(state[ids], `[[`, character(1L),
                                           "error"), collapse = " | ")
    }
    if (any(good)) {
      valid <- which(!failed)
      mapped <- matrix(match(local[good, ], valid), ncol = 2L)
      matrices <- lapply(state[valid], `[[`, "M")
      moments[selected[good], ] <- mgcvst_pair_trace_powers_cpp(
        matrices, mapped, maxPower = 4L, threads = threads
      )
      a <- lapply(state[valid], `[[`, "a")
      score[selected[good]] <- vapply(seq_len(nrow(mapped)), function(k)
        sum(a[[mapped[k, 1L]]] * a[[mapped[k, 2L]]]), numeric(1L))
    }
    rm(state)
    if (exists("matrices", inherits = FALSE)) rm(matrices, a)
    chunks <- chunks + 1L
    first <- last + 1L
  }
  valid <- is.finite(score) &
    apply(is.finite(moments) & moments > 0, 1L, all) &
    moments[, 1L] > 1e-10
  p <- logp <- rep(NA_real_, rows)
  errors[!valid & is.na(errors)] <-
    "Exact Liu score or trace moments are non-finite or non-positive."
  if (any(valid)) {
    selected <- which(valid)
    m <- moments[selected, , drop = FALSE]
    A <- m[, 1L]; B <- m[, 2L]; C <- m[, 3L]; D <- m[, 4L]
    c2 <- A^2 + 3 * B
    c3 <- A^3 + 9 * A * B + 15 * C
    c4 <- A^4 + 18 * A^2 * B + 60 * A * C + 24 * B^2 + 105 * D
    s1 <- c3 / c2^(3 / 2)
    s2 <- c4 / c2^2
    calibrated <- is.finite(s1) & is.finite(s1^2) & s1 > 0 &
      is.finite(s2) & is.finite((score[selected]^2 - A) / sqrt(2 * c2))
    errors[selected[!calibrated]] <-
      "Exact Liu calibration has non-finite moments."
    if (any(calibrated)) {
      z <- selected[calibrated]
      m <- m[calibrated, , drop = FALSE]
      liu <- .liu_squared_score_moments(abs(score[z]), m[, 1L], m[, 2L],
                                        m[, 3L], m[, 4L])
      p[z] <- liu$p_value
      logp[z] <- .mgcvst_liu_logp(liu)
      invalid <- !is.finite(p[z]) | p[z] < 0 | p[z] > 1
      errors[z[invalid]] <- "Exact Liu calibration returned an invalid p-value."
    }
  }
  list(score = score, moments = moments, p_value = p, log10_p = logp,
       error_message = errors, elapsed = proc.time()[["elapsed"]] - elapsed,
       builds = 0L, chunks = chunks, cache_peak_bytes = cache$peak,
       store_reads = cache$reads)
}

.mgcvst_landmark_error_table <- function(approximate, exact, logp, exact_logp) {
  values <- cbind(abs(approximate / exact - 1), abs(logp - exact_logp))
  quantities <- c("t1", "t2", "t3", "t4", "abs_log10_p")
  do.call(rbind, lapply(seq_len(ncol(values)), function(k) {
    z <- values[, k]
    good <- is.finite(z)
    q <- if (any(good)) stats::quantile(z[good], c(0.5, 0.95, 0.99, 1),
                                       names = FALSE) else rep(NA_real_, 4L)
    data.frame(quantity = quantities[k], pairs = length(z), invalid = sum(!good),
               median = q[1L], p95 = q[2L], p99 = q[3L], maximum = q[4L])
  }))
}
