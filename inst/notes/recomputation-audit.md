# mgcvST unnecessary recomputation audit

Implementation update: version 0.0.1.9008 addresses call-wide model-set pair reuse, legacy Davies shared-factor/feature reuse, and retained marginal spectrum replay. See [dense-reuse.md](dense-reuse.md) for the implemented lifecycle and remaining cross-call work. The findings below preserve the original audit snapshot.

Scope: independent, read-only audit of the canonical checkout `C:/Users/yxy1234/Downloads/mgcvST`, snapshot `96e5fca46119c878a8aa6f310d83dd1c807b98e4` plus the uncommitted working-tree state present on 2026-09-13. This report evaluates repeated computation only. It does not assess statistical validity, calibration, style, or propose new analysis behavior.

## Executive finding

There are material avoidable repetitions. Pair testing for a model-set fit uses chunk-local feature-state caching for both Liu and Davies, so either calibration repeats a feature summary when that feature occurs in multiple chunks. Only the legacy/SPDE Liu engine has a call-wide feature-summary lifecycle. The legacy/SPDE Davies path additionally refactorizes the same shared precision once per feature occurrence. Marginal evaluation repeats some geometry-only penalty work per feature and a later `mgcvST.marginal()` call reconstructs the entire spectrum already computed during estimation. WGCNA correctly deduplicates selected genes within one call, but it does not reuse score states already constructed by pair testing or an earlier WGCNA call.

## Findings, ordered by priority

### P1 — Model-set Liu and Davies pair testing recompute feature states at every pair-chunk boundary

**Model-set mgcv branch, both Liu and Davies.** Dispatch through `.mgcvst_test_model_single()` / `.mgcvst_test_model()` reaches `.mgcvst_model_test_chunk()` for either calibration. That chunk function creates a fresh cache environment for every task (`R/model-test.R:20-25`). `.mgcvst_model_cached_state()` stores only `a` and `M` in that environment (`R/model-score.R:201-210`). Pair chunks are formed before `bplapply()` (`R/model-test.R:185-208`). Therefore, for feature *g*, `.mgcvst_model_score_state()` runs once in every chunk containing at least one pair incident on *g*, rather than once per `mgcvST.test()` call, regardless of whether calibration is Liu or Davies. Its exact count is

`count(g) = number of pair chunks whose unique endpoint set contains g`,

bounded by 1 and the number of pair chunks (and by the degree of *g*). Each run rebuilds the feature operator, applies `P` to the error and the complete marked factor, and forms `a` and `M` (`R/model-score.R:181-198`). A hub gene distributed through 20 chunks is consequently summarized 20 times, potentially on different workers.

The fixed SPDE penalty factors are already hoisted correctly: `.mgcvst_model_fixed_factors()` is run once in the parent (`R/model-test.R:205-208`) and the resulting factors are sent with `test_fit`. The remaining state cannot be shared merely by worker persistence because task scheduling is not feature-affine and every task deliberately replaces `.mgcvst_state_cache`.

**Reusable scope.** `a_g` and `M_g` depend on that feature's fitted `W`, dispersion/smoothing parameters, nuisance covariance, and working error, but none of those change between pairs in the same test call. They can be constructed once for every unique endpoint before pair calibration, following the call-wide lifecycle already used by the legacy/SPDE Liu path. They remain valid across all chunks/workers for the call and, while the fit is immutable, across repeated test calls with the same score-component definition. Pair calibration still must run once per pair.

**Legacy/SPDE mgcv branch, Davies only.** For non-Liu calibration, payloads contain each chunk's unique endpoints (`R/mgcvst-api.R:1323-1337`) and `.mgcvst_test_chunk()` creates one summary per unique endpoint in that chunk (`R/mgcvst-api.R:1051-1075`). Thus the same count formula applies. There is no cache across chunks or workers. Its Liu path is separate and call-wide.

### P1 — Legacy/SPDE Davies refactorizes a shared `Q` once per feature occurrence

Within `.mgcvst_test_chunk()`, every endpoint summary calls `rkhs_score_operator()` (`R/mgcvst-api.R:1058-1066`). For positive-definite `Q`, that function performs an eigenvalue decomposition for validation, then a Cholesky factorization and triangular inverse (`R/score-operator.R:82-103`). For PSD `Q`, it performs one eigendecomposition and constructs `Q^{-1/2}` (`R/score-operator.R:82-94`). Yet `geometry$B` and `geometry$Q` are shared by every feature.

Consequently, in the PD case there is one shared-`Q` eigen plus one shared-`Q` Cholesky/solve **per unique feature per pair chunk**. In the PSD case there is one shared-`Q` eigen plus inverse-square-root construction at the same frequency. Only multiplication by `sqrt(field_scale[g])`, the Woodbury matrices driven by `W_g`, and projections driven by `W_g` must vary by feature.

The default Liu route shows the eliminable boundary: it factors `Q` once per complete test call (`R/mgcvst-api.R:897-921`), then scales the common `T0` once per used feature (`R/mgcvst-api.R:927-936`). Davies can use the same call-wide base factor without changing its calibration.

### P2 — A repeated downstream marginal call reconstructs a spectrum already computed at estimation

`mgcvST.estimate()` always runs the marginal score test (`R/mgcvst-api.R:626-635`), and each successful feature invokes `.mgcvst_marginal_score()` (`R/mgcvst-api.R:470-478`; model-set equivalent `R/model-fit.R:191-197`). With `retain_marginal = TRUE`, it separately stores frozen state and shared geometry. Later, `mgcvST.marginal()` calls `.mgcvst_marginal_spectrum()` again for each requested feature (`R/marginal-api.R:160-177, 248-276`).

For the common later request using the same test component, tolerance, and Liu calibration, the second call repeats the working-system reconstruction, scaled penalties, matrix inverses/generalized inverses, `P` applications, test-penalty square root, `B'PB`, and `Q_small` eigenvalues from the estimation call. Count: once per feature during estimation plus once per requested feature per downstream call. Geometry is deduplicated, but computed spectra/results are not.

**Reusable scope.** A successfully computed `(statistic, lambda)` is frozen-fit state and can serve Liu and Davies calibration; the expensive spectrum is independent of the selected calibration method. It is reusable across chunks/workers and subsequent calls when `test_component` and `null.tol` are unchanged. If retaining every `lambda` is considered too large, retaining the four Liu moments plus statistic eliminates repeat work only for Liu. Changed tolerance/component requires reconstruction.

### P2 — Marginal test-penalty generalized inverse and square root are often geometry-only but run per gene

For each feature, `.mgcvst_marginal_spectrum()` builds the target scaled penalty, computes `Thetaj <- ginv(S_matrix / ||S_matrix||_F)` (`R/marginal-taps.R:28-41, 57-84`), then eigendecomposes `Thetaj` to obtain its square root (`R/marginal-taps.R:139-145`). With the package's marked score component contract (one fitted penalty; `R/model-fit.R:61-65`),

`S_matrix = (sp_g / phi_g) S0`, so `S_matrix / ||S_matrix||_F = S0 / ||S0||_F`

for positive `sp_g/phi_g`. Therefore `Thetaj` and `Theta_sqrt` are identical for all ordinary-family features sharing the geometry, yet each is reconstructed once per feature per marginal evaluation. They can be cached once per retained geometry/test component and reused across features, chunks, and workers. For extended-family features whose `valid_idx` changes, `Bj` is row-subsetted, but the coefficient-space normalized penalty remains unchanged; `B'PB` and the final `Q_small` eigen still vary.

By contrast, the following marginal work must remain per feature: PIRLS `V_phi`; nuisance penalty combinations when smoothing-parameter ratios differ; `XtX + S_All` inversion; fixed-effect projection; residual score; `B'PB`; and the final `Q_small` eigendecomposition. These all depend directly or indirectly on feature-specific `W`, `phi`, or smoothing parameters.

### P3 — Marginal retention rebuilds and compares metadata once per feature

Inside each fit chunk, the first retained feature creates `.mgcvst_marginal_geometry`; every later feature calls it again with the shared `X`, constructs the smooth metadata list, zeroes the offset, and performs full `identical()` against the saved geometry (`R/marginal-api.R:19-35`; callers `R/model-fit.R:178-188` and `R/mgcvst-api.R:455-466`). Count: one construction per successful retained feature, scoped to its fitting chunk; collection later deduplicates chunk geometries using further full `identical()` comparisons (`R/marginal-api.R:39-53`).

No lpmatrix is regenerated because the existing `X` is passed, so this is lower priority than the numerical repetitions. For frozen `mgcvST.set()` geometry, the metadata is already contractually shared and can be built once in the parent. For non-frozen/custom geometry, per-feature validation may be required, though a signature/hash or once-per-worker validated seed could avoid repeated deep comparisons where the existing geometry signature has already established equality.

### P3 — WGCNA score construction is efficient within one call but has no cross-call/state reuse

`mgcvST.wgcna()` unions all indices across possibly overlapping blocks (`R/wgcna.R:246-247`) and `.mgcvst_wgcna_scores()` constructs each used feature's score vector exactly once (`R/wgcna.R:127-170`). The shared legacy SPDE factor is decomposed once per WGCNA call (`R/wgcna.R:100-123`), and model-set fixed factors are likewise built once per call (`R/wgcna.R:124-125`). Per-block covariance, correlation, adjacency, TOM, and clustering are then necessarily block-specific (`R/wgcna.R:283-305`), except exact duplicate blocks supplied under different names would repeat identical network work.

The avoidable repetition is across APIs/calls: WGCNA rebuilds `a_g` even if the same immutable fit's pair test just built `a_g` (with `M_g` built by both Liu and Davies on the model-set path), and a repeated WGCNA call rebuilds all selected scores. Count: once per unique selected gene per WGCNA invocation. The state is reusable across blocks and already is; it is also mathematically reusable across pair testing and repeated WGCNA calls for the same component selection, provided cache ownership and memory bounds are explicit. WGCNA correctly avoids computing `M_g`, so forcing it through a full pair-summary cache could add work; a reusable score-only state is the appropriate granularity.

## Existing caches that are effective

- Model fitting seeds the formal training lpmatrix and extracted geometry from the first successful feature and distributes it to later chunks (`R/model-fit.R:294-328`; `R/model-geometry-cache.R:52-83`). For frozen `mgcvST.set()` designs, reuse is unconditional. This removes repeated design/geometry construction, while feature-specific smoothing parameters and offsets are refreshed.
- Model-set pair testing factors each marked fixed penalty once in the parent (`R/model-score.R:24-33`; `R/model-test.R:205-208`). Scaling by feature-specific `sqrt(phi/sp)` correctly remains per feature.
- Default Liu legacy/SPDE testing factors shared `Q` once and constructs each unique feature summary once per complete call (`R/mgcvst-api.R:897-956, 1303-1311`). Pair trace powers remain correctly once per pair.
- WGCNA unions overlapping blocks before score construction, so overlapping blocks do not cause repeated per-feature operators within that invocation (`R/wgcna.R:246-247`).

## Recommended implementation order

1. Give the model-set pair path a call-wide unique-feature summary stage for both Liu and Davies, and give legacy/SPDE Davies the same lifecycle. For legacy/SPDE Davies, hoist the shared `Q` factor at the same time. This removes repetition proportional to feature degree and chunk count.
2. Separate marginal spectrum construction from calibration and retain/reuse the estimation-time spectrum (or Liu moments) when `retain_marginal = TRUE`.
3. Cache normalized target-penalty `Thetaj`/`Theta_sqrt` by retained geometry and test component.
4. If profiling shows value after the numerical fixes, replace repeated marginal metadata reconstruction/deep equality checks with the established shared-geometry contract.
5. Consider an optional fit-owned, bounded score-vector cache shared by pair tests and WGCNA only if repeated downstream calls are common; avoid computing `M` for WGCNA alone.
