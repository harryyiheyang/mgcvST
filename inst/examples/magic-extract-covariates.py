import csv
import argparse
import gzip
import io
import json
import os
import zipfile

import h5py
import numpy as np

parser = argparse.ArgumentParser(description="Extract point metadata from the official MAGIC raw H5AD archive.")
parser.add_argument("archive", help="Mouse_Embryo_Brain_3D_T9_70_50um-raw.zip")
parser.add_argument("order_file", help="Author-supplied sample-section-info.csv")
parser.add_argument("--transfer", default="artifacts/inla3d-transfer/spde3d_transfer_2026-09-13/data/magic/aligned_points.tsv.gz")
parser.add_argument("--out", default="artifacts/datasets/MAGIC/additional-genes")
args = parser.parse_args()
archive = args.archive
transfer = args.transfer
order_file = args.order_file
out = args.out

def column(group, name):
    x = group[name]
    if isinstance(x, h5py.Group):
        categories = x["categories"][:]
        categories = [v.decode() if isinstance(v, bytes) else str(v) for v in categories]
        codes = x["codes"][:]
        return np.array([categories[i] if i >= 0 else None for i in codes], dtype=object)
    x = x[:]
    if x.dtype.kind in "SO":
        return np.array([v.decode() if isinstance(v, bytes) else str(v) for v in x], dtype=object)
    return x

os.makedirs(out, exist_ok=True)
with gzip.open(transfer, "rt") as f:
    points = list(csv.DictReader(f, delimiter="\t"))

with open(order_file, newline="") as f:
    order_rows = list(csv.DictReader(f))
sequence = {r["ID"]: int(r["seq_id"]) for r in order_rows}
if len(sequence) != len(order_rows):
    raise RuntimeError("order file contains duplicate IDs")

z = zipfile.ZipFile(archive)
members = sorted(x.filename for x in z.infolist() if x.filename.endswith(".h5ad"))
if len(members) != 93:
    raise RuntimeError(f"expected 93 H5AD files, found {len(members)}")
records = {}
schemas = []
slices = []
obs_columns = None

for name in members:
    suffix = name.split("adata_")[1][:-5]
    slice_id = f"sample-{suffix}"
    if slice_id not in sequence:
        raise RuntimeError(f"{slice_id} is absent from the order file")
    with h5py.File(io.BytesIO(z.read(name)), "r") as f:
        cols = list(f["obs"].keys())
        if obs_columns is None:
            obs_columns = cols
        if cols != obs_columns:
            raise RuntimeError(f"obs schema differs in {slice_id}")
        obs = {c: column(f["obs"], c) for c in cols}
        spatial = f["obsm"]["spatial"][:]
        n = len(obs["_index"])
        if spatial.shape != (n, 2):
            raise RuntimeError(f"unexpected spatial shape in {slice_id}: {spatial.shape}")
        for i in range(n):
            barcode = str(obs["_index"][i])
            point_id = f"{slice_id}:{barcode}"
            if point_id in records:
                raise RuntimeError(f"duplicate source point ID: {point_id}")
            records[point_id] = {c: obs[c][i] for c in cols if c != "_index"}
            records[point_id]["raw_spatial_x"] = spatial[i, 0]
            records[point_id]["raw_spatial_y"] = spatial[i, 1]
            records[point_id]["h5ad_source_row"] = i
            records[point_id]["section_seq_id"] = sequence[slice_id]
        spatial_uns = f["uns"]["spatial"]
        slices.append({
            "slice_id": slice_id,
            "section_seq_id": sequence[slice_id],
            "n_obs": n,
            "HE_row": int(f["uns"]["HE_row"][()]),
            "HE_col": int(f["uns"]["HE_col"][()]),
            "uns_spatial_keys": ";".join(spatial_uns.keys()),
        })
        schemas.append({
            "slice_id": slice_id,
            "obs_columns": cols,
            "obsm_keys": list(f["obsm"].keys()),
            "uns_keys": list(f["uns"].keys()),
        })

missing = [p["point_id"] for p in points if p["point_id"] not in records]
point_ids = [p["point_id"] for p in points]
if len(point_ids) != len(set(point_ids)):
    raise RuntimeError("transfer point IDs are not unique")
extra = set(records) - set(point_ids)
global_row_mismatch = [p["point_id"] for i, p in enumerate(points)
                       if int(p["source_row"]) != i]
slice_starts = {}
for p in points:
    slice_starts.setdefault(p["slice_id"], int(p["source_row"]))
local_row_mismatch = [p["point_id"] for p in points
                      if p["point_id"] in records and
                      int(p["source_row"]) - slice_starts[p["slice_id"]] !=
                      records[p["point_id"]]["h5ad_source_row"]]
umi_mismatch = [p["point_id"] for p in points
                if p["point_id"] in records and
                int(p["total_umi"]) != int(records[p["point_id"]]["total_counts"])]
if missing or extra or global_row_mismatch or local_row_mismatch:
    raise RuntimeError(f"join failed: missing={len(missing)}, extra={len(extra)}, global_row={len(global_row_mismatch)}, local_row={len(local_row_mismatch)}")

fields = ["point_id", "h5ad_source_row", "section_seq_id", "raw_spatial_x", "raw_spatial_y"]
fields += [c for c in obs_columns if c != "_index"]
with gzip.open(os.path.join(out, "raw-covariates.tsv.gz"), "wt", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(fields)
    for p in points:
        r = records[p["point_id"]]
        vals = [p["point_id"]] + [r[c] for c in fields[1:]]
        w.writerow(["NA" if x is None else x for x in vals])

with open(os.path.join(out, "raw-slice-metadata.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(slices[0]))
    w.writeheader()
    w.writerows(sorted(slices, key=lambda x: x["section_seq_id"]))

audit = {
    "h5ad_files": len(members),
    "transfer_points": len(points),
    "source_h5ad_points": len(records),
    "exact_point_id_matches": len(points) - len(missing),
    "extra_source_points_not_in_transfer": len(extra),
    "source_row_vs_global_transfer_order_mismatches": len(global_row_mismatch),
    "source_row_minus_slice_start_vs_h5ad_row_mismatches": len(local_row_mismatch),
    "total_umi_vs_total_counts_mismatches": len(umi_mismatch),
    "total_umi_equals_raw_total_counts": len(umi_mismatch) == 0,
    "obs_schema_identical_all_sections": len({tuple(x["obs_columns"]) for x in schemas}) == 1,
    "obsm_schema_identical_all_sections": len({tuple(x["obsm_keys"]) for x in schemas}) == 1,
    "uns_schema_identical_all_sections": len({tuple(x["uns_keys"]) for x in schemas}) == 1,
    "obs_columns": obs_columns,
    "obsm_keys": schemas[0]["obsm_keys"],
    "uns_keys": schemas[0]["uns_keys"],
}
with open(os.path.join(out, "raw-covariates-audit.json"), "w") as f:
    json.dump(audit, f, indent=2)
