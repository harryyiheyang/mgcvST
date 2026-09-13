# Dense mgcv computation reuse

The mgcv backend retains its dense score mathematics and both Liu and Davies calibration. These changes concern the lifetime of computed quantities, not the statistical estimator or calibration.

Model-set pair tests construct each unique endpoint's score state once per call. Legacy/SPDE Davies tests also construct the shared precision factor once per call and each unique endpoint once. The existing legacy Liu lifecycle is unchanged. Pair workers consume temporary packed score states rather than rebuilding feature operators for each pair chunk.

A temporary dense score unit contains the score vector, component widths, and one triangle of the symmetric coefficient-space score matrix. This is still quadratic in the score-coordinate count; it is not a sparse representation or a constant-size sufficient statistic. Reconstruction fills the symmetric triangle and requires no matrix factorization or solve. These units are scoped to one test call and are removed afterward; they are not retained in the fitted object or pair-test result.

With `retain_marginal = TRUE`, the built-in marginal calculation retains the statistic and eigenvalue vector already computed at estimation. Subsequent marginal calls reuse them for the same target component and `null.tol`, then perform the requested Liu or Davies calibration. Changed tolerance or component triggers the original spectrum calculation. Custom callbacks and custom design matrices do not enter this cache. The retained marginal state describes the frozen estimate; its internal fields should not be edited independently of its cached spectrum.

WGCNA already constructs one score vector per unique selected gene within a call and does not construct pair-calibration matrices. There is currently no persistent score-vector cache shared between separate pair-test and WGCNA calls. Repeated calls can therefore reconstruct those vectors. Geometry-only marginal penalty algebra is also unchanged to preserve the existing floating-point evaluation order.

## Sparse INLA factorization locations

`src/inla_sparse.cpp` uses Eigen sparse LLT for the shared field precision `Q`. `apply_B` and `apply_Bt` apply its inverse triangular factors with the stored permutation; no full inverse of `Q` is formed. The prepared factor is reused across features.

For pair scores and WGCNA, the feature-dependent expected working curvature is `H = tau * Q + A' W A`. Sparse LDLT solves supply the field and fixed-effect nuisance adjustments. A reconstruction unit retains this sparse factor and its permutation and diagonal; `stored_ldlt_solve` reuses them without factorizing `H` again. The marginal null-target path bypasses this field-curvature factorization.

There is also a dense LDLT of the small fixed-effect nuisance information `J`; it supplies `Vp`. This is distinct from the large sparse field system. All of these are linear-system operations preceding Liu trace calibration, not eigenvalue decompositions for quadratic-form calibration.

More precisely, `A' W A` is the conditional likelihood Fisher block. `H` is the Woodbury auxiliary system for the null marginal covariance, which retains each feature's spatial variance under the pairwise cross-covariance null. It is not the likelihood Fisher block alone. Under the marginal target-variance null, that field is absent from the null covariance, while its covariance kernel remains the tested direction.

## Validation on Windows, 2026-09-13

Version 0.0.1.9008 installed successfully into the default R 4.6 user library. Installed-package tests passed: legacy reuse (6 assertions), dense cache (11), marginal (37), hot path (27), WGCNA (46), INLA sparse equivalence (51), and INLA OpenMP routing (19). Three hot-path checks requiring an explicitly configured historical source checkout were skipped. This was focused regression validation, not a full R CMD check or a new calibration simulation.

The dense cache tests include two local Snow workers, both package fit entry points, and both Liu and Davies. Serial and parallel result tables agree to tolerance 1e-12. Direct pre-cache score formulas are also checked, along with exact triangle restoration and one feature-state construction per unique endpoint. A parameter-name collision with `bplapply` discovered during validation was fixed before these successful runs.
