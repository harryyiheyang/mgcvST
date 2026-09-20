# Small penalized nuisance blocks share the dense expected-information solve
# used for fixed effects. Fixed columns have zero precision; every proper iid
# block repeats its fitted precision over the block's coefficient columns.
.inlast_iid_nuisance_precision <- function(fixed_columns, block_width,
                                           block_precision, features = NULL) {
  fixed_columns <- as.integer(fixed_columns)
  block_width <- as.integer(block_width)
  block_precision <- as.matrix(block_precision)
  storage.mode(block_precision) <- "double"
  if (any(!is.finite(block_precision)) || any(block_precision < 0)) {
    stop("block_precision must contain one non-negative finite row per iid block.")
  }
  if (is.null(features)) features <- ncol(block_precision)
  features <- as.integer(features)
  out <- matrix(0, fixed_columns + sum(block_width), features)
  first <- fixed_columns
  for (j in seq_along(block_width)) {
    rows <- first + seq_len(block_width[j])
    out[rows, ] <- rep(block_precision[j, ], each = block_width[j])
    first <- first + block_width[j]
  }
  out
}
