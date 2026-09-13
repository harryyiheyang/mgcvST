# Full-point native P1 coverage checks; no expression model is fitted.
library(Matrix)
library(INLA)
library(fmesher)
options(warn = 2)
args <- commandArgs(TRUE)
if (length(args) != 1L) stop("Supply the extracted spde3d transfer directory.")
src <- normalizePath(args[1L], winslash = "/", mustWork = TRUE)
out <- "artifacts/inla3d-transfer/full-geometry"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
R1 <- list()
bad <- integer(4L)
for (i in 1:4) {
  key <- c("magic", "langlieb", "mosta", "langlieb")[i]
  name <- if (i == 4L) "aggregated_points.tsv.gz" else "aligned_points.tsv.gz"
  p <- file.path(src, "data", key)
  spec <- jsonlite::read_json(file.path(p, "mesh_contract.json"), simplifyVector = TRUE)
  xyz <- as.matrix(read.csv(file.path(p, "nodes.csv")))
  tv <- as.matrix(read.csv(file.path(p, "tetrahedra.csv")))
  mesh <- fm_mesh_3d(loc = xyz, tv = tv)
  h <- names(read.delim(gzfile(file.path(p, name)), nrows = 0L))
  cls <- rep("NULL", length(h))
  cls[h %in% c("x_aligned", "y_aligned", "z", "total_umi", "bead_count")] <- "numeric"
  t0 <- proc.time()[["elapsed"]]
  d <- read.delim(gzfile(file.path(p, name)), colClasses = cls)
  loc <- sweep(sweep(as.matrix(d[, c("x_aligned", "y_aligned", "z")]),
    2L, spec$source_origin, "-"), 2L, spec$source_units_per_mesh_unit, "/")
  read_seconds <- proc.time()[["elapsed"]] - t0
  t0 <- proc.time()[["elapsed"]]
  A <- inla.spde.make.A(mesh, loc = loc)
  A_seconds <- proc.time()[["elapsed"]] - t0
  re <- max(abs(rowSums(A) - 1))
  ae <- max(abs(as.matrix(A %*% xyz) - loc))
  stopifnot(all(is.finite(A@x)), min(A@x) > -1e-10, max(rowSums(A != 0)) <= 4L)
  b <- which(abs(rowSums(A) - 1) >= 1e-7)
  repaired <- 0L
  repair_error <- 0
  if (length(b)) {
    write.csv(data.frame(source_data_row = head(b, 1000L),
      loc[head(b, 1000L), , drop = FALSE]),
      file.path(out, paste0(key, "-", i, "-uncovered.csv")), row.names = FALSE)
  }
  # Repair only roundoff-sized violations using the supplied tetrahedra.
  for (j in b) {
    best <- -Inf
    for (t in seq_len(nrow(tv))) {
      v <- xyz[tv[t, ], , drop = FALSE]
      w <- solve(t(sweep(v[2:4, , drop = FALSE], 2L, v[1L, ], "-")), loc[j, ] - v[1L, ])
      w <- c(1 - sum(w), w)
      if (min(w) > best) {
        best <- min(w)
        w0 <- w
        t0 <- t
      }
    }
    if (best >= -1e-12) {
      w0 <- pmax(w0, 0)
      w0 <- w0 / sum(w0)
      e <- max(abs(as.numeric(w0 %*% xyz[tv[t0, ], ]) - loc[j, ]))
      if (e > 1e-10) stop("Boundary weight repair exceeds the coordinate tolerance.")
      A[j, tv[t0, ]] <- w0
      repaired <- repaired + 1L
      repair_error <- max(repair_error, e)
    }
  }
  bad[i] <- sum(abs(rowSums(A) - 1) >= 1e-7)
  good <- abs(rowSums(A) - 1) < 1e-7
  covered_error <- max(abs(as.matrix(A[good, ] %*% xyz) - loc[good, ]))
  stopifnot(covered_error < 1e-7)
  R1[[i]] <- data.frame(dataset = key, points = name, n = nrow(d),
    nodes = mesh$n, tetrahedra = nrow(tv), nnz_A = nnzero(A),
    A_MiB = as.numeric(object.size(A)) / 2^20, read_seconds = read_seconds,
    A_seconds = A_seconds, row_error = re, affine_error = ae,
    uncovered_points = length(b), repaired_points = repaired,
    remaining_uncovered_points = bad[i], repair_coordinate_error = repair_error,
    covered_affine_error = covered_error,
    total_umi = if ("total_umi" %in% names(d)) sum(d$total_umi) else NA_real_,
    beads = if ("bead_count" %in% names(d)) sum(d$bead_count) else NA_real_)
  write.csv(do.call(rbind, R1), file.path(out, "coverage.csv"), row.names = FALSE)
  rm(d, loc, A)
  gc()
}
stopifnot(R1[[2L]]$total_umi == R1[[4L]]$total_umi,
  R1[[2L]]$n == R1[[4L]]$beads)
capture.output(sessionInfo(), file = file.path(out, "session.txt"))
if (any(bad)) stop("Full-point coverage failed; see coverage.csv and uncovered-point coordinates.")
