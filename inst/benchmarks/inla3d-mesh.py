"""Geometry-driven tetrahedral pilot; run from the package checkout."""
import csv
import json
import time
from pathlib import Path

import gmsh
import numpy as np

out = Path("artifacts/inla3d/mesh")
out.mkdir(parents=True, exist_ok=True)
rows = []
gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 0)
gmsh.option.setNumber("General.NumThreads", 1)
gmsh.option.setNumber("Mesh.MaxNumThreads3D", 1)
gmsh.option.setNumber("Mesh.RandomSeed", 20260913)
for kind in ("uniform", "adaptive"):
    for budget in (1500, 2800):
        start = time.perf_counter()
        scale = 0.22
        best = None
        trials = []
        for attempt in range(12):
            gmsh.clear()
            gmsh.model.add("stacked_slices")
            gmsh.model.occ.addBox(0, 0, 0, 3, 2, 1)
            gmsh.model.occ.synchronize()
            gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
            gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
            gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
            field = gmsh.model.mesh.field.add("MathEval")
            expr = str(scale)
            if kind == "adaptive":
                expr += "*(1-0.55*Exp(-((x-1.0)^2+(y-1.0)^2+(z-0.5)^2)/0.16))"
            gmsh.model.mesh.field.setString(field, "F", expr)
            gmsh.model.mesh.field.setAsBackgroundMesh(field)
            gmsh.model.mesh.generate(3)
            gmsh.model.mesh.optimize("Netgen")
            tags, xyz, _ = gmsh.model.mesh.getNodes()
            xyz = xyz.reshape(-1, 3)
            types, etags, enodes = gmsh.model.mesh.getElements(3)
            if list(types) != [4]:
                raise RuntimeError("Expected only first-order tetrahedra.")
            order = np.argsort(tags)
            tags, xyz = tags[order], xyz[order]
            tv = np.searchsorted(tags, enodes[0]).reshape(-1, 4)
            if not np.array_equal(tags[tv.ravel()], enodes[0]):
                raise RuntimeError("Tetrahedron node mapping failed.")
            n = len(xyz)
            trials.append({"attempt": attempt + 1, "scale": scale, "nodes": n})
            if n <= budget and (best is None or n > len(best[0])):
                quality = np.asarray(gmsh.model.mesh.getElementQualities(etags[0], "minSICN"))
                best = (xyz.copy(), tv.copy(), quality, scale)
            if best is not None and len(best[0]) >= 0.94 * budget:
                break
            scale *= (n / (0.97 * budget)) ** (1 / 3)
        if best is None or len(best[0]) < 0.85 * budget:
            raise RuntimeError(f"Could not attain the node budget for {kind}, {budget}.")
        xyz, tv, quality, scale = best
        if quality.min() <= 0:
            raise RuntimeError("Inverted or degenerate tetrahedron.")
        label = f"{kind}-{budget}"
        np.savetxt(out / f"{label}-vertices.csv", xyz, delimiter=",", header="x,y,z", comments="")
        np.savetxt(out / f"{label}-tetrahedra.csv", tv + 1, fmt="%d", delimiter=",", header="v1,v2,v3,v4", comments="")
        edges = np.concatenate([tv[:, [a, b]] for a in range(4) for b in range(a + 1, 4)])
        edges = np.unique(np.sort(edges, axis=1), axis=0)
        lengths = np.linalg.norm(xyz[edges[:, 0]] - xyz[edges[:, 1]], axis=1)
        mid = xyz[edges].mean(axis=1)
        roi = np.linalg.norm(mid - [1, 1, 0.5], axis=1) < 0.35
        row = dict(mesh=label, kind=kind, budget=budget, nodes=len(xyz), tetrahedra=len(tv),
                   mesh_seconds=time.perf_counter() - start, target_scale=scale,
                   min_quality=float(quality.min()), median_quality=float(np.median(quality)),
                   median_edge=float(np.median(lengths)),
                   roi_median_edge=float(np.median(lengths[roi])),
                   outside_median_edge=float(np.median(lengths[~roi])))
        rows.append(row)
        (out / f"{label}-search.json").write_text(json.dumps(trials, indent=2))
        print(row, flush=True)
gmsh.finalize()
with (out / "manifest.csv").open("w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=rows[0])
    writer.writeheader()
    writer.writerows(rows)
(out / "software.json").write_text(json.dumps({"gmsh": gmsh.__version__, "numpy": np.__version__}, indent=2))
