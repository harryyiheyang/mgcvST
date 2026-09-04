# Fitting and marginal evaluation

Prepare the common design before estimating genes:

```r
model <- mgcvST.set(
  response ~ celltype + batch + offset(log_depth) +
    s(z, k = 5) + s(x, y, bs = "spdePC", xt = basis),
  data = dat, family = mgcv::nb()
)
# Or wrap an externally prepared mgcv::gam(..., fit = FALSE) setup:
model <- mgcvST.set(G = G)
fit <- mgcvST.estimate(Y, model)
pair <- mgcvST.test(fit, pairs = pairs)
model$timing
lapply(model$G$smooth, function(s) if (!is.null(s$timing)) as.list(s$timing))
```

The formula is complete and parsed by mgcv. Gaussian is supported for SCT;
negative binomial is supported for raw counts. `mgcvST.set()` computes the
formal L without fitting a gene. All workers reuse this L and the same prepared
mesh/penalties. There is no initial gene fit to establish the design and no
per-gene lpmatrix reconstruction. Covariates, factor coding, smooths and row
order are common across genes. `source_files` and `worker_init` are unavailable
on this prepared path. Rebuild the setting to change its design.

`mgcvST.estimate(Y, model, offset = O)` adds `O` to the shared formula/setup
offset. `O` is either an observation-length vector or a matrix matching the
feature-by-observation dimensions and ordering of Y. L stays shared; working
quantities, retained marginal offsets and nuisance Vp blocks remain per gene.
The result's `offset` stores the total offsets, in the same orientation as O.

SPDE and SPDE-PC construction and prediction always interpolate at the input
coordinate rows. They keep mesh, transform, Z, Q and PC eigenvectors fixed.
SPDE-PC evaluates `A(L2) (Z V)`; its full score basis is also updated to
`A(L2) Z`. No stale training rows or tolerance-based coordinate matching are
used. Each smooth's timing environment records cumulative basis/prediction
seconds and calls in the current process. `model$timing` records total setup
and formal-lpmatrix seconds. Re-evaluation of PC multiplication can change
floating-point rounding: the fixed upstream p-value checks differed by at
most 4.7e-13 in absolute value.

The `model.set()` and raw-G paths below remain available for compatibility;
the unconditional common-design contract applies to `mgcvST.set()`.

`mgcvST.estimate()` fits each gene and immediately calls
the package-local `taps_score_test()` with its exact cached training lpmatrix.
The original TAPS arithmetic is included in mgcvST; no mgcv.taps installation
or sourced score function is needed.
Marginal testing is mandatory. Wood diagnostics remain optional.

```r
fit <- mgcvST.estimate(Y, G, BPPARAM = bp)
fit <- mgcvST.estimate(Y, G, BPPARAM = bp,
                      marginal_args = list(method = "liu"))
pair <- mgcvST.test(fit, pairs = pairs)
```

The default TAPS method is Davies with Liu fallback for numerical failures.
Direct Liu skips Davies. The registered BiocParallel backend runs gene fits
and their marginal tests; each worker uses one numerical thread. This inline
path does not collect all gene spectra into a parent OpenMP batch. Pairwise
calibration and its OpenMP implementation are unchanged.

`fit$diagnostics` retains `marginal_requested_method`, `marginal_method` and
`marginal_fallback` for each gene. The fallback flag identifies Davies-to-Liu
changes; direct Liu is FALSE. Custom callbacks without method metadata retain
their p-values and leave the unavailable method/fallback fields as NA.

The matrix comes from the formal mgcv predictor, never directly from G$X.
Within a gene, compaction and TAPS share that matrix. Supported deterministic
model geometries also reuse it across genes. Unknown methods and custom worker
initialization obtain the current gene's matrix without sharing across genes.
Full GAM objects and transient design caches are discarded after fitting.

`marginal_test` still accepts an explicitly supplied custom callback, and
`marginal_args` selects calibration for the built-in implementation. The
existing `retain_marginal` / `mgcvST.marginal()` API is retained for explicit
recalibration of stored fits; it is not required by the estimation workflow.
No additional score-test interface is introduced in mgcvST.

For an HPC installation without INLA, sf, fmesher or mgcv.taps, prepare the boundary,
mesh and `spde_basis` on Windows, then transfer the saved basis together with
the aligned data. The HPC process can read the basis, construct `mgcvST.set`,
estimate, test and predict without loading those four packages. Uploading
only a boundary still requires mesh construction; that operation uses sf and
fmesher. Existing-mesh interpolation uses geometry, whose required dependency
chain does not include INLA, sf or fmesher.

Strict-equivalence pair optimizations are independent of these API additions:

* one F-transpose-PF multiplication per model feature state;
* successful a/M states cached lazily within each pair chunk;
* fixed marked-SPDE factors built once per test call, with failure conditions
  replayed inside the existing pair error handler;
* no changes to the existing ordinary-SPDE Liu summary cache, to the pair
  calibration formula, to signed dot products, or to numerical cutoffs.

No warm starts, family reuse across fits, or replacement of predicted
lpmatrices by G$X is introduced. For model.set batches with registered SPDE
methods and supported ordinary mgcv low-rank smoothers, the first successful
fit obtains the exact formal training lpmatrix once. One nuisance design LN is
retained for the batch, while each successful feature retains only its
conditional Vp nuisance block; full Vp and GAM objects are discarded.

The primary model score operator uses only the tested SPDE in Vs and applies
the conditional penalized nuisance adjustment through LN and VpN. It never
forms an n-by-n P. Projection products use CppMatrix matrixMultiply and solves
use CppMatrix matrixSolve. Unknown/custom prediction methods, general.family,
cross-penalty setups and coefficient rank drops retain the previous operator
path. Supplying source_files or worker_init also disables shared prediction
geometry, because worker initialization may change method dispatch. This
path is also skipped when a parallel batch has fewer than two features per
worker, avoiding an extra fit wave merely to establish the shared cache. This
statistical refactor is tested for numerical rather than bitwise
equivalence. Pure caching and duplicate-product changes remain subject to
strict identity. The original diagnostic and result fields remain present;
disabled optional results are NA.

Within one feature state, Vs-inverse times LN and its product with VpN are each
computed once and reused for the working response and tested factors. These
intermediates are released after a and M have been cached for the current pair
chunk; they are not retained across all features.
