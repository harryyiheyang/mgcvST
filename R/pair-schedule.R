# Pair schedules: the order and the blocks in which the pairs are evaluated.

# Reuse two resident gene blocks before loading the next pair of blocks.
.mgcvst_pair_order <- function(index, used, capacity) {
  width <- max(1L, floor(capacity / 2))
  block1 <- (match(index[, 1L], used) - 1L) %/% width
  block2 <- (match(index[, 2L], used) - 1L) %/% width
  if (capacity >= length(used)) seq_len(nrow(index)) else
    order(pmin(block1, block2), pmax(block1, block2), method = "radix")
}

# Blocks of the all-pairs universe of `n` genes (positions 1, ..., n). Genes
# are cut into tiles of `width`; for every tile pair (a, b), a <= b, the left
# genes of tile a are cut into sub-blocks whose partner count stays at or
# below `chunk_size` (a single left gene may exceed it). Each row holds the
# left position range l1:l2, the right position range r1:r2, whether the tile
# is the diagonal one (partners v > u only) and the pair count.
.mgcvst_all_pair_blocks <- function(n, width, chunk_size) {
  starts <- seq.int(1L, n, by = width)
  ends <- pmin(n, starts + width - 1L)
  rows <- list()
  for (a in seq_along(starts)) {
    left <- starts[a]:ends[a]
    for (b in a:length(starts)) {
      diagonal <- a == b
      partners <- if (diagonal) ends[a] - left else
        rep(ends[b] - starts[b] + 1L, length(left))
      cs <- cumsum(as.numeric(partners))
      first <- 1L
      while (first <= length(left)) {
        base <- if (first > 1L) cs[first - 1L] else 0
        stop_at <- max(first, findInterval(base + chunk_size, cs))
        count <- sum(partners[first:stop_at])
        if (count == 0) break
        rows[[length(rows) + 1L]] <- c(l1 = left[first], l2 = left[stop_at],
          r1 = starts[b], r2 = ends[b], diagonal = diagonal, count = count)
        first <- stop_at + 1L
      }
    }
  }
  do.call(rbind, rows)
}

# Gene-position pairs (u < v) of one block of .mgcvst_all_pair_blocks().
.mgcvst_block_pairs <- function(block) {
  left <- block[["l1"]]:block[["l2"]]
  if (block[["diagonal"]] > 0) {
    len <- block[["r2"]] - left
    cbind(i = rep(left, len), j = sequence(len, from = left + 1L))
  } else {
    right <- block[["r1"]]:block[["r2"]]
    cbind(i = rep(left, each = length(right)),
          j = rep(right, times = length(left)))
  }
}

# The block schedule of an all-pairs universe.
.mgcvst_all_pair_schedule <- function(n, capacity, chunk_size) {
  width <- if (capacity >= n) n else max(1L, as.integer(floor(capacity / 2)))
  chunk_size <- as.integer(min(chunk_size, .Machine$integer.max))
  list(width = as.integer(width), chunk_size = chunk_size,
       blocks = .mgcvst_all_pair_blocks(as.integer(n), as.integer(width), chunk_size))
}
