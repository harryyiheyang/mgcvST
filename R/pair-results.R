# Compact pair-test results: one schema, streamed Parquet shards, one
# multiple-testing adjustment.
#
# Every pair of a test is one row of integer feature indices i < j, the signed
# score, the natural-log two-sided, positive and negative p-values, the
# adjusted two-sided log q-value, the kind of remainder used by the
# calibration (0 none, 1 one node, 2 two nodes, 3 Gaussian) and a status code.
# Feature names are looked up from the `feature_id` vector stored beside the
# shards; no per-pair character column exists. Shards are written as the pairs
# are evaluated, so a run never holds a per-pair table of all pairs in memory.

.mgcvst_final_columns <- c(
  "i", "j", "score", "log_p_two_sided", "log_p_positive", "log_p_negative",
  "log_q", "remainder_kind", "status"
)

# status: 0 evaluated; 1 trace moments non-finite or non-positive; 2 invalid
# p-value; 3 a gene of the pair has no usable score state; 4 a gene was fitted
# spatially under spatial = "all" but was not selected by the Stage 1 null score
# test of the fit. A pair with status 1, 2 or 4 has p = 1: two-sided log p = 0,
# and both one-sided log p = 0, and it stays in the adjustment family. A pair
# with status 3 has no p-value and is not adjusted.
.mgcvst_pair_status <- c(ok = 0L, moments = 1L, p_value = 2L, feature = 3L,
                         degenerate = 4L)

# Algorithm contract of the pair p-values. Every checkpoint of pair results is
# keyed by it, and a checkpoint written under another contract is refused. The
# rank `k` and the sha of the shared basis are results of the run, not part of
# the algorithm: a directory is refused when its contract differs in the
# contract string, the route, the remainder order, the kernel version or the
# schema.
# Stage 2 is calibrated by a saddlepoint approximation on k leading singular
# values of the pair spectrum plus a remainder that matches the remaining power
# sums: four moments with two nodes on the exact route (`remainder_order` 4),
# two moments with one node on the PCAlearning route (2). The rank `k` and the
# sha of the shared basis are known only once the basis exists, so the pair
# directory of a run is opened (and the contract written) after the basis is
# built; only the refusal of stale directories happens earlier.
.mgcvst_contract <- function(route, k = NA_integer_, basis_sha = NA_character_) {
  stopifnot(route %in% c("exact", "pcalearning"))
  list(calibration_contract = "spa_v1", route = route, k = as.integer(k),
       remainder_order = if (identical(route, "exact")) 4L else 2L,
       basis_sha = as.character(basis_sha), kernel_version = 3L,
       schema = "compact_v1")
}

# The part of a contract that defines the algorithm.
.mgcvst_contract_algorithm <- function(contract) {
  contract[setdiff(names(contract), c("k", "basis_sha"))]
}

.mgcvst_pairs_frame <- function(i, j, score = NA_real_, log_p_two_sided = NA_real_,
                                log_p_positive = NA_real_,
                                log_p_negative = NA_real_, remainder_kind = 0L,
                                status = 0L) {
  n <- length(i)
  data.frame(
    i = as.integer(i), j = as.integer(j),
    score = rep_len(as.numeric(score), n),
    log_p_two_sided = rep_len(as.numeric(log_p_two_sided), n),
    log_p_positive = rep_len(as.numeric(log_p_positive), n),
    log_p_negative = rep_len(as.numeric(log_p_negative), n),
    remainder_kind = rep_len(as.integer(remainder_kind), n),
    status = rep_len(as.integer(status), n),
    stringsAsFactors = FALSE
  )
}

# Natural-log adjusted p-values of the two-sided family: BY (default), BH,
# Sidak or none, computed in log space by the native kernel.
.mgcvst_log_adjust <- function(log_p, adjust = c("BY", "BH", "Sidak", "none")) {
  adjust <- match.arg(adjust)
  mgcvst_log_adjust_cpp(as.numeric(log_p), adjust)
}

# Bytes of the adjustment vectors per pair: log p, log q and the ranking.
.mgcvst_adjust_bytes_per_pair <- 24
# Bytes of one compact result row held in memory.
.mgcvst_result_bytes_per_pair <- 56

# Memory guard: does `bytes_per_pair * n_pairs` stay below `fraction` of the
# memory available to the process? Unknown memory passes. Memory the run has
# released but R has not yet returned is collected first, so that the probe
# does not count it as used.
.mgcvst_pair_memory_guard <- function(n_pairs, bytes_per_pair, fraction) {
  gc(FALSE)
  available <- .mgcvst_memory_probe()$available
  need <- bytes_per_pair * as.numeric(n_pairs)
  list(ok = !is.finite(available) || need <= fraction * available,
       need = need, available = available)
}

.mgcvst_shard_file <- function(dir, id) {
  file.path(dir, sprintf("shard-%010d.parquet", as.integer(id)))
}

.mgcvst_read_shard <- function(file, columns = NULL) {
  # No memory map: a mapped shard cannot be replaced or deleted on Windows.
  z <- if (is.null(columns)) arrow::read_parquet(file, mmap = FALSE) else
    arrow::read_parquet(file, col_select = arrow::all_of(columns), mmap = FALSE)
  as.data.frame(z)
}

.mgcvst_write_parquet <- function(frame, file) {
  tmp <- tempfile("shard-", tmpdir = dirname(file), fileext = ".tmp")
  on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
  arrow::write_parquet(frame, tmp)
  if (file.exists(file)) unlink(file)
  if (!file.rename(tmp, file)) {
    stop("Could not commit the pair result shard ", basename(file), ".")
  }
  invisible(file)
}

# Refuse pair results written under another algorithm contract. An older
# checkpoint is never resumed silently: its p-values come from another
# calibration.
.mgcvst_pairs_refuse_stale <- function(root, contract) {
  dirs <- list.files(root, "^pairs-", full.names = TRUE)
  dirs <- dirs[dir.exists(dirs)]
  for (dir in dirs) {
    record <- file.path(dir, "contract.rds")
    found <- if (file.exists(record)) {
      tryCatch(readRDS(record)$contract, error = function(e) NULL)
    } else NULL
    if (!identical(.mgcvst_contract_algorithm(found),
                   .mgcvst_contract_algorithm(contract))) {
      stop("The checkpoint directory ", root, " holds pair results written ",
           "under a different algorithm contract (",
           if (is.null(found)) "none recorded" else
             paste0(found$calibration_contract, ", ", found$route, ", kernel ",
                    found$kernel_version),
           "; this call writes ", contract$calibration_contract, ", ",
           contract$route, ", kernel ", contract$kernel_version, "): ", basename(dir),
           ". Use a new checkpoint_dir, or delete the pairs-* directories to ",
           "keep the reusable feature score states.", call. = FALSE)
    }
  }
  invisible(NULL)
}

# Open (or resume) the directory of raw pair shards of one pair universe.
.mgcvst_pairs_open <- function(root, universe, contract, resume = TRUE) {
  if (!dir.exists(root) && !dir.create(root, recursive = TRUE)) {
    stop("Could not create the pair result directory.")
  }
  .mgcvst_pairs_refuse_stale(root, contract)
  # The universe of an explicit pair list holds the full index; it is hashed in
  # bounded pieces rather than serialised whole.
  sha <- digest::digest(list(version = 3L, contract = contract,
                             universe = .mgcvst_pair_input_hash(universe)),
                        algo = "sha256")
  dir <- file.path(root, paste0("pairs-", sha))
  record <- list(version = 3L, contract = contract, universe_sha = sha)
  if (dir.exists(dir)) {
    if (!resume) {
      stop("The pair results in ", dir, " already exist; use resume = TRUE or ",
           "a new checkpoint_dir.")
    }
    if (!identical(tryCatch(readRDS(file.path(dir, "contract.rds")),
                            error = function(e) NULL), record)) {
      stop("The pair results in ", dir, " are damaged or were written for ",
           "another pair universe.")
    }
  } else {
    if (!dir.create(dir)) stop("Could not create the pair result directory.")
    tmp <- tempfile("contract-", tmpdir = dir, fileext = ".tmp")
    saveRDS(record, tmp, compress = FALSE)
    if (!file.rename(tmp, file.path(dir, "contract.rds"))) {
      stop("Could not commit the pair result contract.")
    }
  }
  dir
}

# TRUE when the raw shard `id` exists and holds exactly the pairs
# (`i`, `j`); an existing shard of other pairs is an error.
.mgcvst_shard_complete <- function(dir, id, i, j) {
  file <- .mgcvst_shard_file(dir, id)
  if (!file.exists(file)) return(FALSE)
  z <- .mgcvst_read_shard(file, c("i", "j"))
  if (!identical(z$i, as.integer(i)) || !identical(z$j, as.integer(j))) {
    stop("The pair result shard ", basename(file), " is damaged or does not ",
         "match this pair schedule.")
  }
  TRUE
}

# Number of pairs of an existing raw shard starting a window at `first` of the
# ordered pair matrix `index`; NA when the shard is absent.
.mgcvst_shard_window <- function(dir, id, index) {
  file <- .mgcvst_shard_file(dir, id)
  if (!file.exists(file)) return(NA_integer_)
  z <- .mgcvst_read_shard(file, c("i", "j"))
  n <- nrow(z)
  rows <- seq.int(id, length.out = n)
  if (!n || max(rows) > nrow(index) || !identical(z$i, as.integer(index[rows, 1L])) ||
      !identical(z$j, as.integer(index[rows, 2L]))) {
    stop("The pair result shard ", basename(file), " is damaged or does not ",
         "match this pair schedule.")
  }
  n
}

# Adjust the two-sided family once, write the final shards (the raw columns
# plus log_q) and summarize the discoveries. `shards` are raw shard files in
# evaluation order with `rows` pairs each; `extra` optionally holds pairs that
# were not evaluated (a gene without a usable state), written as one more
# shard. Positive and negative discoveries are the adjusted two-sided
# discoveries split by the sign of the score.
#
# When the compact table fits the memory guard it is built in one pass while
# the pairs are written: the columns are preallocated, filled shard by shard,
# and put in (i, j) order by a single permutation. A temporary run then holds
# the table in memory and deletes its shards and pair directory, so that
# nothing is left in the temporary directory; a run that does not fit keeps
# its final shards, and so does every run with a checkpoint directory.
.mgcvst_pairs_finalize <- function(pair_dir, shards, rows, extra, adjust,
                                   q.value, temporary, materialize = TRUE,
                                   verbose = FALSE) {
  n_extra <- if (is.null(extra)) 0 else nrow(extra)
  n_pairs <- sum(as.numeric(rows)) + n_extra
  out_dir <- file.path(pair_dir, paste0("adjusted-", adjust))
  if (!dir.exists(out_dir) && !dir.create(out_dir)) {
    stop("Could not create the adjusted pair result directory.")
  }
  old <- list.files(out_dir, "^results-.*[.]parquet$", full.names = TRUE)
  if (length(old)) unlink(old)

  adjustment <- list(method = adjust, computed = FALSE, reason = NA_character_,
                     n_adjusted = NA_real_)
  log_q <- NULL
  if (identical(adjust, "none")) {
    adjustment$computed <- TRUE
  } else {
    guard <- .mgcvst_pair_memory_guard(n_pairs, .mgcvst_adjust_bytes_per_pair, 0.4)
    if (!guard$ok) {
      adjustment$reason <- paste0(
        "the log p-value vectors of ", format(n_pairs, big.mark = ","),
        " pairs would exceed the safe memory line (",
        format(guard$need / 1024^3, digits = 3), " GiB needed, ",
        format(guard$available / 1024^3, digits = 3), " GiB available)")
      warning("Multiple-testing adjustment skipped: ", adjustment$reason,
              "; log_q is NA. The raw log p-values are in the shards.",
              call. = FALSE)
    } else {
      lp <- numeric(n_pairs)
      at <- 0
      if (n_extra) {
        lp[seq_len(n_extra)] <- extra$log_p_two_sided
        at <- n_extra
      }
      for (file in shards) {
        z <- .mgcvst_read_shard(file, "log_p_two_sided")[[1L]]
        lp[at + seq_along(z)] <- z
        at <- at + length(z)
      }
      adj <- .mgcvst_log_adjust(lp, adjust)
      log_q <- adj$log_q
      adjustment$computed <- TRUE
      adjustment$n_adjusted <- adj$n
      rm(lp, adj)
    }
  }

  # The in-memory table is decided before the write loop, after the adjustment
  # vectors are released.
  keep <- materialize && n_pairs > 0 && n_pairs <= .Machine$integer.max &&
    .mgcvst_pair_memory_guard(n_pairs, .mgcvst_result_bytes_per_pair, 0.2)$ok
  columns <- NULL
  if (keep) {
    columns <- new.env(parent = emptyenv())
    columns$i <- integer(n_pairs)
    columns$j <- integer(n_pairs)
    for (name in c("score", "log_p_two_sided", "log_p_positive",
                   "log_p_negative", "log_q")) {
      columns[[name]] <- numeric(n_pairs)
    }
    columns$remainder_kind <- integer(n_pairs)
    columns$status <- integer(n_pairs)
  }
  write_files <- !(temporary && keep)

  threshold <- log(q.value)
  counts <- c(tested = 0, with_p = 0, discovered = 0, positive = 0, negative = 0)
  log_p_threshold <- -Inf
  any_discovery <- FALSE
  final <- character()
  at <- 0
  write_part <- function(part, name) {
    k <- nrow(part)
    part$log_q <- if (identical(adjust, "none")) {
      ifelse(is.na(part$log_p_two_sided), NA_real_,
             pmin(part$log_p_two_sided, 0))
    } else if (is.null(log_q)) {
      rep(NA_real_, k)
    } else log_q[at + seq_len(k)]
    if (keep) {
      slot <- at + seq_len(k)
      for (column in .mgcvst_final_columns) columns[[column]][slot] <- part[[column]]
    }
    at <<- at + k
    counts[["tested"]] <<- counts[["tested"]] +
      sum(part$status != .mgcvst_pair_status[["feature"]])
    counts[["with_p"]] <<- counts[["with_p"]] + sum(!is.na(part$log_p_two_sided))
    found <- !is.na(part$log_q) & part$log_q <= threshold
    if (any(found)) {
      any_discovery <<- TRUE
      counts[["discovered"]] <<- counts[["discovered"]] + sum(found)
      counts[["positive"]] <<- counts[["positive"]] +
        sum(found & !is.na(part$score) & part$score > 0)
      counts[["negative"]] <<- counts[["negative"]] +
        sum(found & !is.na(part$score) & part$score < 0)
      log_p_threshold <<- max(log_p_threshold, part$log_p_two_sided[found])
    }
    if (write_files) {
      file <- file.path(out_dir, name)
      .mgcvst_write_parquet(part[, .mgcvst_final_columns, drop = FALSE], file)
      final <<- c(final, file)
    }
    invisible(NULL)
  }
  if (n_extra) write_part(extra, "results-0000000000.parquet")
  for (file in shards) {
    part <- .mgcvst_read_shard(file)
    write_part(part, sub("^shard-", "results-", basename(file)))
    if (temporary) unlink(file)
  }
  rm(log_q)

  results <- NULL
  if (keep) {
    perm <- order(columns$i, columns$j, method = "radix")
    if (is.unsorted(perm)) {
      for (column in .mgcvst_final_columns) {
        columns[[column]] <- columns[[column]][perm]
      }
    }
    rm(perm)
    results <- stats::setNames(mget(.mgcvst_final_columns, envir = columns),
                               .mgcvst_final_columns)
    attr(results, "row.names") <- .set_row_names(as.integer(n_pairs))
    class(results) <- "data.frame"
    rm(columns)
  }
  if (temporary && (keep || !length(final))) {
    parent <- dirname(pair_dir)
    unlink(pair_dir, recursive = TRUE)
    if (startsWith(basename(parent), "mgcvst-pairs-")) {
      unlink(parent, recursive = TRUE)
    }
  }
  discovered <- if (adjustment$computed) {
    list(pairs_discovered = counts[["discovered"]],
         pairs_discovered_positive = counts[["positive"]],
         pairs_discovered_negative = counts[["negative"]])
  } else {
    list(pairs_discovered = NA_real_, pairs_discovered_positive = NA_real_,
         pairs_discovered_negative = NA_real_)
  }
  list(
    results = results, shards = final, n_pairs = n_pairs,
    adjustment = adjustment,
    threshold = list(q_value = q.value, adjust = adjust,
      log_p_threshold = if (any_discovery) log_p_threshold else NA_real_),
    discoveries = c(list(pairs_requested = n_pairs,
                         pairs_tested = counts[["tested"]],
                         pairs_with_p_value = counts[["with_p"]]), discovered)
  )
}
