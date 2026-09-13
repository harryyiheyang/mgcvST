# Real tissue tetrahedral transfer audit

The 2026-09-13 transfer supplies three real tissue meshes and one executable
gene example. MAGIC supports an INLA-only fit on all 97,830 observations with
1,962 mesh nodes. Two real-data fits and ten conditional simulation fits
completed locally with flat spatial log precision and NB log size.
The original archive and its successful script are preserved in
`artifacts/inla3d-transfer/spde3d_transfer_2026-09-13/`.

## Input identity and available evidence

Archive: `spde3d_transfer_2026-09-13.zip`, 193,862,890 bytes, SHA-256
`6fe5cdd6bd95047ded6580a43495ec2818d2a37ff1dcf481f8fad7a80191a7dc`.
All supplied SHA256SUMS entries were verified before execution. The archive
contains 41 files. Its source and alignment descriptions remain the authority
for provenance; this audit does not independently revalidate registration.

| Dataset | Slices | Full points | Mesh nodes | Tetrahedra | Expression supplied |
|---|---:|---:|---:|---:|---|
| MAGIC | 93 | 97,830 | 1,962 | 7,676 | Snap25 and total UMI |
| Langlieb | 98 | 5,429,770 | 4,048 | 16,409 | Total UMI only |
| MOSTA E9.5 | 70 | 646,893 | 4,013 | 17,903 | None |

Langlieb also supplies 918,923 aggregate points and bead-to-bin mappings.
Summed total UMI agrees between full and aggregate tables at 2,778,586,291;
aggregate bead counts sum to 5,429,770. Mapping-level membership was not
independently re-audited in this task.

MAGIC coordinates convert from micrometres to millimetres. Langlieb and MOSTA
use the supplied separate per-axis normalizations; their physical XY:Z
calibration remains unresolved. Langlieb retains historical registration
drift and missing-slice spacing. MOSTA uses the author's aligned coordinates.
Layerwise convex envelopes can bridge holes and separate tissue pieces.

The transfer contains completed mesh nodes and tetrahedra, but no mesh
generation or adaptive refinement program. Therefore it supports reuse and
fitting of the supplied meshes, rather than reproduction of their construction.
Only MAGIC is below 3,000 nodes. The two approximately 4,000-node meshes have
geometry validation here, without an expression fit.

## Original script reproduction

The unchanged `scripts/check_inla_3d.R --fit` passed all three sampled native
3D mesh/A/Q checks and its 5,000-point Snap25 fit on this computer. Local fit
time was 6.14 seconds, count RMSE 2.19206, and count correlation 0.73318.
The supplied run reported 11.78 seconds; these are separate computer runs.
Local versions were R 4.6.1, INLA 26.6.8 and fmesher 0.8.0; the supplied run
used INLA 26.8.7. Full records are in `inst/validation/inla3d-transfer/baseline/`.

The original script uses a standard normal prior on spatial log sigma, a
normal prior on NB log size with mean log(10) and SD 1.5, and intercept
precision 1/9. Its approximation is simplified Laplace with CCD integration.
These are distinct from the requested flat-prior configuration.

## Flat-prior implementation

The new standalone script retains the supplied nodes, tetrahedra, kappa,
coordinate transform, point/count matching, UMI exposure `total_umi/10000`,
NB variant 0 and native volume-integral constraint. It uses a flat intercept
and Gaussian empirical-Bayes approximation with VB disabled, consistent with
the existing package's flat-objective estimator convention.

The native alpha=2 model in three dimensions has nu=1/2 and

```
Q(1) = (kappa^4 C + 2 kappa^2 G + G C^-1 G) / (8 pi kappa)
Q(sigma) = Q(1) / sigma^2.
```

The relative maximum difference between native Q(1) and the assembled FEM
expression was 2.09e-16. The integral constraint is taken directly from the
native SPDE object; it is not replaced with the earlier synthetic study's
observation-mean constraint. Native spatial `theta` is log(sigma), while
generic0 uses log precision = -2 log(sigma), so a constant density on either
internal scale changes only by a constant Jacobian. See the
[official INLA flat-prior definition](https://www.inla.r-inla-download.org/r-inla.org/doc/prior/prior-flat.pdf).

On local INLA 26.6.8, specifying `prior="flat"` directly for native spde2
failed in `inla_parse_ffield` before optimization. A verbose diagnostic
confirmed the parser received the flat prior. A separate attempt to fix the
native spatial hyperparameter also crashed; no native/generic fitted-result
equivalence is claimed. Diagnostics remain local. The completed route uses
generic0 with the native Q and integral constraint, diagonal=0 and rankdef=1.
Spatial log precision and NB log size both explicitly use `prior="flat"`.

The flat objective is used for point estimation with conditional Gaussian
uncertainty. This run does not establish propriety of an integrated posterior,
interior optimality from status=0 alone, or interval coverage. See
`inla-flat-prior.md`. Original-versus-new fit differences also include the
approximation and fixed-effect prior changes, and cannot be attributed solely
to hyperpriors.

## Completed fitting checks

| Fit | Observations | Seconds | NB size | Spatial sigma |
|---|---:|---:|---:|---:|
| Flat, original evenly spaced subset | 5,000 | 3.88 | 17.9408 | 0.4360 |
| Flat, full MAGIC Snap25 | 97,830 | 6.86 | 15.0437 | 0.3847 |
| Ten simulations, median | 97,830 each | 6.28 | 15.1768 | 0.3455 |

All 12 fits had mode status 0, finite positive parameter estimates, and no
reported INLA warnings. Maximum absolute integral constraint error was
5.04e-15. The full real-data count RMSE was 2.30824 and correlation 0.68532;
these describe in-sample fit, rather than SVG significance.

Ten independently seeded NB response vectors were generated from the full
fit's fixed log-rate surface, original exposures and NB size 15.0437. The
median log-rate RMSE was 0.05812 (range 0.05696-0.05981), with correlation
0.96975 (0.96800-0.97099). This conditional recovery experiment fixes one
smoothed surface. It does not draw ten independent latent SPDE fields or
test unbiased spatial variance estimation, mesh optimality, repeated-sample
interval coverage, or a 3D BAM comparison. Raw responses, seeds and estimated
fields are retained in the local RDS checkpoints.

Times describe individual runs with two INLA threads; some checks overlapped
other local work. They exclude data preparation and do not measure peak RAM.

## Full-point coverage and boundary handling

Native P1 interpolation covers all MAGIC and MOSTA points. Langlieb has one
unlocated raw point (data row 5,423,701) and the corresponding aggregate point
(row 917,612). Both have the same supplied coordinates on the upper boundary.
Direct barycentric evaluation in supplied tetrahedron 3,980 gives minimum
weight -1.38e-14. Clamping that roundoff-sized negative weight and renormalizing
reconstructs the supplied coordinate to 2.66e-15 mesh units in the independent
Python diagnostic.

The full-geometry script records native failures, permits a correction only
when the best minimum weight is at least -1e-12, verifies coordinate error
at most 1e-10, and fails for remaining uncovered points. It does not move
observations, alter the supplied mesh, silently remove points, or relax
tolerance for substantially out-of-domain coordinates. See `coverage.csv`
for the actual repaired count and errors. Its row/affine error columns record
the native matrix before repair; `covered_affine_error` is checked after repair.

The interpolation matrices remain sparse, at no more than four weights per
point. Approximate sparse A storage is 4.48 MiB for MAGIC, 248.46 MiB for
Langlieb raw points, 42.04 MiB for its aggregates and 29.54 MiB for MOSTA.
These object sizes are not process peak memory.

## Reproduction and next implementation step

Run from the canonical checkout after extracting the archive:

```powershell
Rscript inst/benchmarks/inla3d-transfer.R artifacts/inla3d-transfer/spde3d_transfer_2026-09-13
Rscript inst/benchmarks/inla3d-transfer-geometry.R artifacts/inla3d-transfer/spde3d_transfer_2026-09-13
```

The fit script stops if its output marker already exists, preserving completed
checkpoints. Original archive hashes and local results identify this run;
bulk data and fitted objects remain under ignored `artifacts/`.

Adaptive mesh comparisons should hold the supplied kappa fixed across node
budgets. Recomputing kappa from each mesh's edge length changes both the
discretization and the stochastic model. Here MAGIC kappa=5.1803213/mm gives
the conventional sqrt(8 nu)/kappa scale of 0.3860764 mm; the supplied h=0.54599
mm is an edge-length statistic. The next geometry work needs the original
mesh-generation program to preserve its envelope and quality rules. Full
Langlieb/MOSTA expression fitting additionally needs gene count matrices.

The public mgcvST API remains 2D, package version 0.0.1.9006. These are standalone
INLA-only research checks; this task does not add a public 3D API.
