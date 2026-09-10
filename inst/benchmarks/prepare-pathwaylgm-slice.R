# Read the original local slice without changing any PathwayLGM files.
suppressPackageStartupMessages(library(Matrix))
input <- Sys.getenv("MGCVST_DLPFC_INPUT",
                    "D:/PathwayLGM/data/real_dlpfc_go/samples/151673/sample.rds")
out <- Sys.getenv("MGCVST_INLA_OUTPUT", "artifacts/pathwaylgm-151673")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
x <- readRDS(input)
i <- sort(x$analysis_spot_indices)
counts <- x$counts[, i, drop = FALSE]
xy <- as.matrix(x$coords[i, , drop = FALSE])
labels <- x$labels[i]
stopifnot(identical(colnames(counts), rownames(xy)),
          identical(colnames(counts), names(labels)), !anyNA(labels),
          all(counts@x >= 0), all(counts@x == round(counts@x)))
lib <- as.numeric(Matrix::colSums(counts))
stopifnot(all(lib > 0))
detected <- Matrix::rowMeans(counts > 0)
average <- Matrix::rowMeans(counts)
eligible <- which(detected >= .05 & average > .05)
symbols <- as.character(x$genes$gene_symbol[match(rownames(counts), x$genes$gene_id)])
# Fixed named panel plus six expression strata; manual layer labels are unused.
named <- c("MBP", "PLP1", "MOBP", "RELN", "GAD1", "SLC17A7")
selected <- vapply(named, function(s) {
  j <- intersect(which(symbols == s), eligible)
  if (length(j)) j[which.max(average[j])] else NA_integer_
}, integer(1L))
selected <- unique(selected[!is.na(selected)])
pool <- setdiff(eligible[order(average[eligible], rownames(counts)[eligible])], selected)
strata <- unique(round(seq(1, length(pool), length.out = 12L - length(selected))))
selected <- c(selected, pool[strata])
stopifnot(length(selected) == 12L)
Y <- as.matrix(counts[selected, , drop = FALSE])
d <- data.frame(x = xy[, 1L], y = xy[, 2L],
                offset0 = log(lib) - mean(log(lib)), row.names = rownames(xy))
# Same fixed mesh for both estimators, including boundary points outside tissue.
margin <- .025 * max(diff(range(d$x)), diff(range(d$y)))
nodes <- as.matrix(expand.grid(
  x = seq(min(d$x) - margin, max(d$x) + margin, length.out = 25L),
  y = seq(min(d$y) - margin, max(d$y) + margin, length.out = 25L)
))
mesh <- list(loc = nodes, graph = list(tv = geometry::delaunayn(nodes)))
genes <- data.frame(gene_id = rownames(Y), symbol = symbols[selected],
                    mean_count = average[selected], detection_fraction = detected[selected])
provenance <- list(
  sample = x$sample_id, donor = x$donor_id, source_path = normalizePath(input),
  source_md5 = unname(tools::md5sum(input)), source_spots = ncol(x$counts),
  retained_spots = length(i), retained_rule = "analysis_spot_indices (labeled in-tissue)",
  labels = "Original layer_guess_reordered: Layer1-Layer6 and WM; not used in fitting",
  selection = "Six predefined symbols plus six mean-count strata, eligibility detected>=5% and mean>0.05",
  offset = "log full-library UMI minus its observation mean",
  mesh = "25x25 rectangular finite-element grid; same mesh for both estimators",
  constraints = "Each spatial field has mean zero at the actual 3611 model observations"
)
saveRDS(list(data = d, mesh = mesh, Y = Y, labels = labels,
             genes = genes, provenance = provenance), file.path(out, "benchmark-input.rds"))
write.csv(genes, file.path(out, "selected-genes.csv"), row.names = FALSE)
write.csv(data.frame(barcode = rownames(d), d, label = labels),
          file.path(out, "retained-spots.csv"), row.names = FALSE)
writeLines(capture.output(str(provenance)), file.path(out, "input-provenance.txt"))
print(genes)
print(table(labels))
cat("Prepared", nrow(Y), "genes x", ncol(Y), "spots;", nrow(nodes), "mesh vertices\n")
