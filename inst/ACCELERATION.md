# A: high-throughput fitting and marginal evaluation

The high-throughput defaults now skip optional computations:

```r
fit <- mgcvST.estimate(Y, G, marginal = FALSE, diagnostics = FALSE)
pair <- mgcvST.test(fit, pairs = pairs)
```

To reproduce the former fitted summaries and callback behavior, explicitly use
`marginal = TRUE, diagnostics = TRUE`. Existing marginal_test/marginal_args and
source-first callbacks remain available. This legacy adapter retains the
callback's own calibration behavior, including any upstream fallback; it is
not the new explicit-calibration API described below.

For independent marginal TAPS, save its inputs without testing during fitting:

```r
fit <- mgcvST.estimate(Y, G, retain_marginal = TRUE)
marginal <- mgcvST.marginal(fit, calibration = "liu", threads = 2L)
marginal <- mgcvST.marginal(fit, calibration = "davies", BPPARAM = bp)
# This is an explicit calibration change for numerical failures only:
marginal <- mgcvST.marginal(fit, calibration = "davies", fallback = "liu",
                          BPPARAM = bp)
```

`retain_marginal` saves minimal fitted response, linear predictor, prior weights,
coefficients, smoothing parameters, dispersion and serialized family state;
common lpmatrix and penalty metadata are stored once. It does not retain GAMs,
run summary.gam, construct a TAPS spectrum, or refit anything. It adds O(n*p)
storage and is deliberately opt-in. Old compact objects without this data
cannot reconstruct exact TAPS from pairwise working quantities alone.

The independent engine ports the standard/extended-family, refit=FALSE TAPS
definition from mgcv.taps fb48abb. This is a conditional frozen-fit evaluation,
not a null-refitted Rao test. Its fit-space design (including reduced PC fit
space) is separate from the pairwise test's full score-space geometry. The
first release covers mgcvST-supported families; it does not add Cox, ZIP, qgam
or ordered-categorical support or offer null refitting.

BiocParallel computes independent feature spectra. Every worker uses one
numerical thread. Liu spectral powers are summed by a separate marginal-only
OpenMP kernel in bounded blocks in the parent process, after worker jobs have
returned. The kernel never invokes R callbacks from OpenMP. Chi-square tail
probabilities still use the original scalar R calculation. This is not a
reuse of the pairwise trace-power engine, and it is not a C++ rewrite of TAPS
spectrum construction.

The new Davies path treats nonzero ifault, non-finite/out-of-range probability
or a zero tail as a numerical failure. It reports NA unless the caller
explicitly requests `fallback="liu"`. Actual method, fallback reason and
Davies ifault are returned for every feature. No implicit screening, FDR or
pair-universe change is performed.

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
statistical refactor is tested for numerical rather than bitwise
equivalence. Pure caching and duplicate-product changes remain subject to
strict identity. The original diagnostic and result fields remain present;
disabled optional results are NA.

Within one feature state, Vs-inverse times LN and its product with VpN are each
computed once and reused for the working response and tested factors. These
intermediates are released after a and M have been cached for the current pair
chunk; they are not retained across all features.
