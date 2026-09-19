# MAGIC E18.5 three-dimensional brain data

The reproducible research object is written to
`artifacts/datasets/MAGIC/MAGIC.rds`. It contains 97,830 observations, three
raw-count genes, 39 point-level covariates, metadata for 93 measured sections,
and the 1,962-node tetrahedral mesh used by the native three-dimensional
examples. The three genes are Snap25, Foxp1, and Tfap2b. Foxp1 and Tfap2b were
selected because the source paper displays them as representative regional
markers in Figure 5d.

The source article reports 98,192 spots for the E18.5 brain atlas. The 93 raw
H5AD files in the official Zenodo archive contain exactly 97,830 observations,
and all 97,830 match the transfer data. The reason for the difference of 362
spots between the article and the archived raw H5AD files has not been
resolved. The current data should therefore be described by their observed
size rather than as a documented QC subset of the article total.

## Source files

The official raw archive is
`Mouse_Embryo_Brain_3D_T9_70_50um-raw.zip` from
[Zenodo record 13934668](https://zenodo.org/records/13934668). Its published
MD5 is `d0485e65ba20b960c03739d8851fc68e`. The archive contains 93 H5AD files
and `sample-section-info.csv`. The portable alignment, Snap25 counts, and mesh
are under
`artifacts/inla3d-transfer/spde3d_transfer_2026-09-13/data/magic`.

The H5AD files do not contain cell-type or anatomical-region labels. Their
`reg1`–`reg9` values identify capture regions, and their `sample` values match
the section filenames. Neither field is treated as a biological batch. The
author-supplied `seq_id` is retained as `section_seq_id`, an ordering key whose
more specific interpretation is not documented in the archive or the authors'
analysis notebook.

## Build

Run the commands from the package root. Python requires `h5py`, `numpy`, and
`scipy`. Pass the locations of the downloaded archive and its bundled section
order file explicitly:

```sh
python inst/examples/magic-extract-paper-genes.py /path/to/Mouse_Embryo_Brain_3D_T9_70_50um-raw.zip

python inst/examples/magic-extract-covariates.py /path/to/Mouse_Embryo_Brain_3D_T9_70_50um-raw.zip /path/to/sample-section-info.csv
```

The first command extracts Foxp1 and Tfap2b counts. It stops before writing
outputs if either feature is absent from any section, if point identifiers are
not unique, or if counts are not non-negative integers. The second command
checks all 93 metadata schemas, verifies the global and within-section row
mapping, and extracts the raw point and slice metadata. Both commands default
to the canonical transfer and output paths in this repository; `--transfer`
and `--out` can override those paths.

Assemble and validate the R object after both extraction commands finish:

```sh
Rscript inst/examples/magic-data-prepare.R
```

The R script joins every file by `point_id`, verifies count and exposure
values, keeps raw H5AD `total_counts` separate from transfer `total_umi`, adds
millimetre coordinates and the model exposure, and records SHA-256 hashes for
the source files. It writes `MAGIC.rds`, `covariate-dictionary.csv`,
`slices.csv`, and `source-files.csv` under `artifacts/datasets/MAGIC`.

## Object structure

- `covariates`: 97,830 rows and 39 columns. These include aligned and raw
  coordinates, transfer exposure, section identifiers, and the stored H5AD QC
  summaries.
- `expression`: a 97,830 by 3 integer matrix with Snap25, Foxp1, and Tfap2b
  raw counts. Its row names equal `covariates$point_id`.
- `slices`: 93 rows with the aligned z coordinate, observed spacing, retained
  observation count, author-supplied ordering key, and stored H&E dimensions.
- `meshes$native3d`: mesh nodes, tetrahedra, and the coordinate contract.
- `genes`, `dictionary`, and `provenance`: gene identifiers, field meanings,
  and source-file hashes.
- `metadata`: scope, provenance, unavailable annotations, gene-selection
  rationale, and the model used by the saved three-dimensional fit.

The object makes the available observations and their provenance directly
inspectable. It does not reproduce the complete transcriptome, and the added
QC fields are metadata rather than terms in the existing fitted model.
