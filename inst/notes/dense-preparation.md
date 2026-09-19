# Parallel feature preparation for covariance tests

The mgcv test paths now construct feature score summaries in a C++ OpenMP
batch before evaluating pairs. The legacy Liu and Davies paths share this
kernel. Current `model.set()` fits use the same kernel with their fitted
conditional nuisance covariance. Older serialized model fits without that
covariance retain the existing R/BiocParallel path.

The shared field factor is formed once. Each OpenMP iteration then computes
one feature's working precision, nuisance adjustment, score vector and score
covariance. The calculation retains the original Woodbury formulas and never
constructs an observation-by-observation covariance or projection matrix.
R handles input and output objects, feature errors and packed-state files.
Preparation uses batches of at most 32 features. Pair calculations begin after
the preparation threads finish; Snow pair workers and preparation threads do
not run as nested parallel layers.

For example, with a previously fitted object and an explicit pair universe:

```r
ans <- mgcvST.test(fit, pairs = pairs, threads = 4L)
```

`threads` controls the native preparation batch and the legacy Liu pair
kernel. Model pair evaluation and Davies calibration continue to use
`BPPARAM`. The default thread count is `bpworkers(BPPARAM)` for mgcv fits.
INLA defaults to one thread and already uses C++ OpenMP for feature preparation
and matrix materialization. Its reconstruction and cache policy are unchanged.

## Measurement

The benchmark uses the first 32 available features of the saved unadjusted
paired-validation mgcv fit: 2,125 observations and 298 score coordinates.
The same fitted object, shared factor and conditional nuisance covariances
are supplied to both implementations. Three consecutive runs were measured
per configuration on the same Windows host using R 4.6.1. These times measure
feature-state computation, including call input/output conversion; they exclude
fitting, the common factor, file storage, and pair testing.

| Preparation | Threads | Median seconds |
|---|---:|---:|
| R feature loop with CppMatrix operations | 1 | 8.79 |
| C++ feature loop | 1 | 8.66 |
| C++ feature loop | 2 | 4.31 |
| C++ feature loop | 4 | 2.31 |

The four-thread batch is 3.81 times faster than the serial R loop for this
input. Single-thread performance is nearly unchanged: parallel feature work
accounts for the improvement. The maximum absolute differences from the R
reference were `9.06e-14` for score vectors and `8.88e-16` for score covariance
entries. Results agreed across thread counts. End-to-end test speedup also
depends on the size of the requested pair universe and its calibration cost.

The flat benchmark is `inst/benchmarks/dense-preparation.R`; it reads the saved
fit at `artifacts/inla-bam-validation/components/unadjusted/bam-compact-fit.rds`.
Per-run measurements are retained in
`inst/validation/dense-preparation.csv`. Full logs and session details are in
`artifacts/package-maintenance-20260919/`.

Package tests compare against direct observation-space precision matrices,
the previous feature-state implementation, and public Liu/Davies results.
They cover conditional nuisance covariance, a zero-column nuisance design,
rank deficiency, reordered features, failed-feature isolation, and multiple
thread counts. Model timing now separates preparation from pair evaluation.
