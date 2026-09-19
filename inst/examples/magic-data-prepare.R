# Assemble the available MAGIC observations in one portable research object.
src <- "artifacts/inla3d-transfer/spde3d_transfer_2026-09-13/data/magic"
out <- "artifacts/datasets/MAGIC"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
d <- read.delim(gzfile(file.path(src, "aligned_points.tsv.gz")), check.names = FALSE)
g <- read.delim(gzfile(file.path(src, "Snap25.tsv.gz")), check.names = FALSE)
s <- read.csv(file.path(src, "slices.csv"))
spec <- jsonlite::read_json(file.path(src, "mesh_contract.json"), simplifyVector = TRUE)
ix <- match(d$point_id, g$point_id)
if (anyNA(ix) || anyDuplicated(d$point_id) || anyDuplicated(g$point_id) ||
    !all(d$total_umi == g$total_umi[ix])) {
  stop("MAGIC coordinate and count identities do not match.")
}
if (anyNA(d) || anyNA(g$count[ix]) || any(d$total_umi <= 0) ||
    any(g$count < 0) || any(g$count != round(g$count))) {
  stop("The available MAGIC observations contain invalid input values.")
}
s <- s[order(s$z), , drop = FALSE]
rownames(s) <- NULL
if (anyDuplicated(s$slice_id) || anyDuplicated(s$z)) stop("The slice-to-z mapping is not unique.")
ix_s <- match(d$slice_id, s$slice_id)
if (anyNA(ix_s) || !all(d$z == s$z[ix_s]) ||
    !all(tabulate(ix_s, nrow(s)) == s$observations)) {
  stop("The MAGIC slice table does not agree with the observations.")
}
s$slice_order <- seq_len(nrow(s))
s$gap_um <- c(NA_real_, diff(s$z))
d$x_mm <- (d$x_aligned - spec$source_origin[1L]) / spec$source_units_per_mesh_unit[1L]
d$y_mm <- (d$y_aligned - spec$source_origin[2L]) / spec$source_units_per_mesh_unit[2L]
d$z_mm <- (d$z - spec$source_origin[3L]) / spec$source_units_per_mesh_unit[3L]
d$slice_order <- s$slice_order[ix_s]
d$exposure <- d$total_umi / 10000
Y <- matrix(as.integer(g$count[ix]), ncol = 1L,
             dimnames = list(d$point_id, "Snap25"))
genes <- data.frame(symbol = "Snap25", ensembl_id = "ENSMUSG00000027273")
extra <- file.path(out, "additional-genes", c("Foxp1.tsv.gz", "Tfap2b.tsv.gz"))
for (p in extra) {
  h <- read.delim(gzfile(p), check.names = FALSE)
  j <- match(d$point_id, h$point_id)
  if (anyNA(j) || anyDuplicated(h$point_id) || nrow(h) != nrow(d) ||
      !all(h$total_umi[j] == d$total_umi) || anyNA(h$count[j]) ||
      any(h$count[j] < 0) || any(h$count[j] != round(h$count[j])) ||
      length(unique(h$gene_symbol)) != 1L || length(unique(h$ensembl_gene_id)) != 1L) {
    stop("The additional MAGIC gene does not match the retained observations: ", p)
  }
  Y <- cbind(Y, matrix(as.integer(h$count[j]), ncol = 1L,
    dimnames = list(d$point_id, unique(h$gene_symbol))))
  genes <- rbind(genes, data.frame(symbol = unique(h$gene_symbol),
    ensembl_id = unique(h$ensembl_gene_id)))
}
nodes <- as.matrix(read.csv(file.path(src, "nodes.csv")))
tv <- as.matrix(read.csv(file.path(src, "tetrahedra.csv")))
storage.mode(tv) <- "integer"
if (min(tv) < 1L || max(tv) > nrow(nodes)) stop("Invalid tetrahedral node index.")
dictionary <- data.frame(
  name = names(d),
  role = c("identifier", "group identifier", "spatial coordinate", "spatial coordinate",
           "spatial coordinate", "sequencing depth", "source index", "spatial coordinate",
           "spatial coordinate", "spatial coordinate", "ordered slice index", "offset exposure"),
  unit = c("", "", "micrometre", "micrometre", "micrometre", "UMI count", "zero-based row",
           "millimetre", "millimetre", "millimetre", "retained-slice rank", "total UMI / 10000"),
  source = c(rep("transfer", 7L), rep("derived without changing observation order", 5L)),
  stringsAsFactors = FALSE)
meta_file <- file.path(out, "additional-genes", "raw-covariates.tsv.gz")
slice_file <- file.path(out, "additional-genes", "raw-slice-metadata.csv")
h <- read.delim(gzfile(meta_file), check.names = FALSE)
hs <- read.csv(slice_file, check.names = FALSE)
j <- match(d$point_id, h$point_id)
js <- match(s$slice_id, hs$slice_id)
if (anyNA(j) || anyNA(js) || anyDuplicated(h$point_id) ||
    anyDuplicated(hs$slice_id) || nrow(h) != nrow(d) || nrow(hs) != nrow(s)) {
  stop("The raw H5AD metadata do not match the retained MAGIC observations.")
}
h <- h[j, , drop = FALSE]
hs <- hs[js, , drop = FALSE]
stopifnot(identical(paste0("sample-", h$sample), d$slice_id),
  all(d$source_row == seq_len(nrow(d)) - 1L),
  all(h$h5ad_source_row == d$source_row - ave(d$source_row, d$slice_id, FUN = min)),
  all(hs$n_obs == s$observations), all(h$section_seq_id == hs$section_seq_id[ix_s]))
new_names <- setdiff(names(h), "point_id")
if (any(new_names %in% names(d))) stop("Raw metadata would overwrite an existing covariate.")
d <- cbind(d, h[, new_names, drop = FALSE])
s <- cbind(s, hs[, setdiff(names(hs), "slice_id"), drop = FALSE])
role <- rep("stored count or QC summary", length(new_names))
role[new_names %in% c("sample", "reg")] <- "section or capture-region identifier"
role[new_names == "h5ad_source_row"] <- "zero-based index within source H5AD section"
role[new_names == "section_seq_id"] <- "author-supplied section ordering key"
role[new_names %in% c("raw_spatial_x", "raw_spatial_y", "Spot_col", "Spot_row")] <- "original spatial or capture-grid coordinate"
unit <- rep("as stored in source H5AD", length(new_names))
unit[grepl("^pct_", new_names)] <- "percent"
unit[new_names %in% c("raw_spatial_x", "raw_spatial_y")] <- "not specified by source metadata"
dictionary <- rbind(dictionary, data.frame(name = new_names, role = role, unit = unit,
  source = "official raw H5AD / bundled section ordering CSV", stringsAsFactors = FALSE))
files <- file.path(src, c("aligned_points.tsv.gz", "Snap25.tsv.gz", "slices.csv",
                          "nodes.csv", "tetrahedra.csv", "mesh_contract.json"))
files <- c(files, extra, meta_file, slice_file)
provenance <- data.frame(
  file = basename(files),
  sha256 = vapply(files, digest::digest, character(1L), algo = "sha256", file = TRUE),
  stringsAsFactors = FALSE)
MAGIC <- list(
  covariates = d, expression = Y, slices = s,
  meshes = list(native3d = list(loc = nodes, tv = tv, contract = spec)),
  genes = genes,
  dictionary = dictionary, provenance = provenance,
  metadata = list(
    dataset = "MAGIC", status = "research snapshot; three available genes",
    paper_url = "https://pmc.ncbi.nlm.nih.gov/articles/PMC11525186/",
    source = "spde3d_transfer_2026-09-13; MAGIC full retained observations",
    alignment = "supplied GEASO rigid aligned coordinates",
    count_scale = "raw UMI counts; expression rows match covariates rows",
    unavailable = c("cell-type labels", "anatomical-region labels", "confirmed donor/batch IDs",
                    "sex or biological replicate IDs", "complete transcriptome matrix"),
    raw_metadata = "All 93 source H5AD files; preserve raw total_counts separately from transfer total_umi",
    covariate_scope = "QC fields are available metadata, not additional terms in the existing fitted model",
    selection = "Snap25 transfer example; Foxp1 and Tfap2b selected from the source paper Fig. 5d",
    retained_subset = "97,830 aligned observations; source paper reports 98,192 spots",
    current_fitted_model = "intercept + 3D SPDE; log(exposure) offset",
    slice_effect = "candidate for exploration; not present in the saved fit"))
saveRDS(MAGIC, file.path(out, "MAGIC.rds"), compress = "xz")
check <- readRDS(file.path(out, "MAGIC.rds"))
stopifnot(identical(check, MAGIC), identical(rownames(check$expression), d$point_id))
write.csv(dictionary, file.path(out, "covariate-dictionary.csv"), row.names = FALSE)
write.csv(s, file.path(out, "slices.csv"), row.names = FALSE)
write.csv(provenance, file.path(out, "source-files.csv"), row.names = FALSE)
