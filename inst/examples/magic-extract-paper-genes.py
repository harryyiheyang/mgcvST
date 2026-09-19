import csv
import argparse
import gzip
import io
import os
import zipfile

import h5py
import numpy as np
import scipy.sparse as sp

parser = argparse.ArgumentParser(description="Extract paper marker genes from the official MAGIC raw H5AD archive.")
parser.add_argument("archive", help="Mouse_Embryo_Brain_3D_T9_70_50um-raw.zip")
parser.add_argument("--transfer", default="artifacts/inla3d-transfer/spde3d_transfer_2026-09-13/data/magic")
parser.add_argument("--out", default="artifacts/datasets/MAGIC/additional-genes")
args = parser.parse_args()
archive = args.archive
transfer = args.transfer
out = args.out
genes = {
    "Foxp1": "ENSMUSG00000030067",
    "Tfap2b": "ENSMUSG00000025927",
}

os.makedirs(out, exist_ok=True)
z = zipfile.ZipFile(archive)
members = sorted(x.filename for x in z.infolist() if x.filename.endswith(".h5ad"))
if len(members) != 93:
    raise RuntimeError(f"expected 93 H5AD files, found {len(members)}")
values = {}

for name in members:
    sec = name.split("adata_")[1][:-5]
    with h5py.File(io.BytesIO(z.read(name)), "r") as f:
        ens = [x.decode() for x in f["var"][dict(f["var"].attrs)["_index"]][:]]
        pos = {g: i for i, g in enumerate(ens)}
        n, p = dict(f["X"].attrs)["shape"]
        X = sp.csr_matrix((f["X"]["data"][:], f["X"]["indices"][:],
                           f["X"]["indptr"][:]), shape=(n, p))
        absent = [symbol for symbol, ensembl in genes.items() if ensembl not in pos]
        if absent:
            raise RuntimeError(f"{', '.join(absent)} absent from the feature table in sample-{sec}")
        cols = [pos[ensembl] for ensembl in genes.values()]
        sub = np.asarray(X[:, cols].todense())
        if (not np.all(np.isfinite(sub)) or np.any(sub < 0) or
                not np.all(np.equal(sub, np.round(sub)))):
            raise RuntimeError(f"counts are not non-negative integers in sample-{sec}")
        obs = [x.decode() for x in f["obs"][dict(f["obs"].attrs)["_index"]][:]]
        if len(obs) != n:
            raise RuntimeError(f"expression and observation rows differ in sample-{sec}")
        for i, barcode in enumerate(obs):
            point_id = f"sample-{sec}:{barcode}"
            if point_id in values:
                raise RuntimeError(f"duplicate source point ID: {point_id}")
            values[point_id] = sub[i]

points = []
with gzip.open(os.path.join(transfer, "aligned_points.tsv.gz"), "rt") as f:
    points = list(csv.DictReader(f, delimiter="\t"))

missing = [x["point_id"] for x in points if x["point_id"] not in values]
point_ids = [x["point_id"] for x in points]
if len(point_ids) != len(set(point_ids)):
    raise RuntimeError("transfer point IDs are not unique")
extra = set(values) - set(point_ids)
if missing or extra:
    raise RuntimeError(f"point-ID join failed: missing={len(missing)}, extra={len(extra)}")

for j, (symbol, ensembl) in enumerate(genes.items()):
    path = os.path.join(out, f"{symbol}.tsv.gz")
    with gzip.open(path, "wt", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["point_id", "gene_symbol", "ensembl_gene_id", "count", "total_umi"])
        for x in points:
            w.writerow([x["point_id"], symbol, ensembl,
                        int(values[x["point_id"]][j]), int(x["total_umi"])])

with open(os.path.join(out, "extraction-summary.csv"), "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["gene_symbol", "ensembl_gene_id", "observations", "id_matches",
                "nonzero", "detection_fraction", "mean_count", "max_count"])
    for j, (symbol, ensembl) in enumerate(genes.items()):
        x = np.array([values[p["point_id"]][j] for p in points])
        w.writerow([symbol, ensembl, len(x), len(x), int(np.count_nonzero(x)),
                    float(np.mean(x > 0)), float(np.mean(x)), int(np.max(x))])
