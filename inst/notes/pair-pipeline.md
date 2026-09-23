# Pairwise Liu execution and reference-trace approximation

## Audit baseline

The audit baseline is GitHub/main `219d25ab9d1dfbcb8f3c86c54b731f3a9630e4f0`
(0.0.1.9014). The common INLA observation-kernel projection retaining 0.995
of the eigenvalue sum remains unchanged. It is computed before the gene
states and is shared by every gene. There is no gene-specific 0.995 cutoff.

| Route in the baseline | Gene preparation | Pair execution | Bottleneck |
|---|---|---|---|
| INLA, `R/inla-test.R:.mgcvst_inla_test_pairs` | One feature on a cache miss; 512 MiB LRU eviction allows reconstruction | Singleton calls to the native pair kernel | Neither feature nor pair OpenMP receives a useful batch; repeated sparse factorization and materialization |
| Legacy mgcv Liu, `R/mgcvst-api.R:.mgcvst_liu_summaries` | Unique features in native batches | Native pair batches | All full double matrices remain resident, without an adaptive memory budget |
| `model.set()` mgcv Liu, `R/model-test.R:.mgcvst_test_model` | Unique states written to temporary packed shards | BiocParallel groups, with serial pair kernels | No pair-level OpenMP batching; temporary shards cannot resume another call |

The 512 MiB setting was a software default, not evidence of an HPC limit.
`src/score_pair.cpp` and `src/inla_sparse.cpp` explicitly disable Eigen inner
parallelism. `src/dense_prepare.cpp` calls Armadillo BLAS/LAPACK inside a
feature OpenMP loop. Environment variables alone cannot reliably change an
already initialized BLAS runtime; the shared thread limiter now also calls
`RhpcBLASctl::blas_set_num_threads(1)`.

Davies calibration retains its existing execution route. The new checkpoint
and adaptive-cache arguments currently apply to Liu calibration.

## Exact pipeline

`R/pair-pipeline.R` now prepares each required gene once, writes an atomic
upper-triangle double shard, and subsequently only reads these shards. A
persistent `checkpoint_dir` validates a fingerprint of the statistical inputs
and shared coordinates. It reuses both completed feature states and completed
pair-result batches. Different pair requests can reuse the feature states.
The caller's pair order is restored after computation.

Persistent fingerprints hash large dense inputs in bounded blocks instead
of serializing the entire fit to a temporary disk file. Runs without a
persistent checkpoint skip the full content fingerprint. Fingerprint version
2 deliberately rejects development checkpoints made by the older scheme.

Feature preparation uses dynamic OpenMP scheduling, with sequential work
inside each feature. Pair batches use the native OpenMP trace kernel. Gene
blocks group pair work so two resident blocks are reused before proceeding
to another block. A small cache causes shard reads, not reconstruction of M.

`cache_bytes = NULL` probes system headroom, Linux cgroup v1/v2 limits and
Slurm allocations. The minimum available signal controls the automatic
budget. Feature workspaces and pair buffers are reserved separately; probes
are refreshed during preparation and at pair-batch boundaries. A supplied
`cache_bytes` is a resident-cache ceiling, not a total process-memory limit.
If no memory source is available, the reported conservative fallback is
512 MiB. This fallback is not a hardware assumption.

All q-by-q covariance storage and native arithmetic remain double in the
exact branch. A 150 GB allocation therefore does not imply that 11,099 full
1500-by-1500 double matrices all fit: those matrices alone require about
186 GiB. Disk shards and block scheduling remain useful in this case.

## Approximate pipeline

`approximate = TRUE` is opt-in and now uses real-gene landmark trace-CUR. The
selection methods are uniform random sampling, k-means on the unnormalized
feature score vectors (`score`), and k-means on standardized fitted covariance
variance scales (`hyper`). Selection is reproducible and restores the caller's
RNG state. The common 0.995 observation-kernel projection is unchanged. Matrix
B learning was explored separately and is not connected to this pipeline.

For each trace power, the method forms the landmark trace matrix W and the
gene-to-landmark trace matrix C. The two GEMM stages constructing C use
float32 arithmetic with trace sums accumulated in double; W is computed in
double. Both C and W are stored in double. Only
reference self traces provide diagonal scaling; non-reference self traces are
not required. W is symmetrized, then inverted through a signed truncated
eigendecomposition: eigenvalues with absolute value below
`ref_tol * max(abs(eigenvalues))` are dropped, while the signs of retained
eigenvalues are preserved. This is needed because higher trace powers need not
be positive semidefinite. The resulting low-rank trace reconstruction and exact
score inner products feed the existing Liu formula; no full G-by-G trace table
is formed.

In score selection, feature units and exact score vectors are written to
shards before selecting landmarks. Sparse units retain the feature H
factorization and are reused to materialize the double curvature matrix once;
they are not rebuilt after score-based selection. Other selection methods can
choose landmarks before building feature states. The non-sparse score
selection path currently constructs operator units serially in R; other dense
native preparation routes remain available.

Missing double M states are materialized once per gene in adaptive R batches
and written as packed upper-triangle double native binary shards. The
implementation then makes one R-to-C++ call for all missing gene trace tasks.
C++ dynamically schedules those
genes across workers, reads each double state shard, performs the two float32
GEMM stages, accumulates the four trace powers in double, and writes native
binary summaries for non-reference genes. Reference summaries remain RDS
files, which the R summary reader prefers when resuming. C++17 is required for
the native file queue. Although trace summaries are streamed by the native
workers, the full double M shards remain on disk for exact holdout diagnostics
and exact tail rechecks; the algorithm is not summary-only streaming. The
temporary state store is removed at the end by default, while `checkpoint_dir`
preserves states, summaries, and exact pair batches for resumption.

By default, up to 10,000 pairs whose two endpoints are non-landmarks are
sampled for diagnostics; set `diagnostic_pairs = 0` to disable this check.
Approximate pairs with p-values below `tail_recheck` (default `1e-3`) are
recomputed exactly before the outer wrapper performs multiple-testing
adjustment; set
`tail_recheck = 0` to disable exact tail checks. This threshold rule does not
provide a strict BY error guarantee. The approximation still reports invalid
moments as failures rather than clipping them.

## CUR precheck and validation status

### Current real-gene CUR precheck (in progress)

A uniform 500-gene sample was drawn from the 11,184-gene fit. It uses the
unchanged common rank r = 1404 projection and separates 340 selector-training
genes from 160 held-out genes; the comparison uses the same 10,000 pair set.
For 500 nearest-neighbor pairs in standardized hyperparameter space, the
relative M difference had median 0.3849487 and 95th percentile 1.548682;
500 random pairs had median 0.5043123 and 95th percentile 2.083156. Among
164 directed near pairs whose two variance parameters each differed by less
than 5%, the median relative M difference
was 0.4015642. These curvature discrepancies show that similar hyperparameter
scales do not ensure similar gene-specific curvature.

The sample has fixed range 0.386 and no fitted observation-level nugget; neither
is a varying landmark feature. The completed R = 100 comparison gives:

| Selection | t2 median relative error | t4 95th percentile | t4 maximum | Absolute log10-p 95th percentile | Absolute log10-p maximum |
|---|---:|---:|---:|---:|---:|
| Score k-means | 0.01176% | 0.13173% | 0.56731% | 0.001350 | 0.02913 |
| Hyperparameter k-means | 0.00904% | 0.09070% | 0.56078% | 0.001169 | 0.04466 |
| Random | 0.00933% | 0.09194% | 1.11416% | 0.001226 | 0.05905 |

These errors precede tail replacement. All 10,000 held-out pairs are common
to the three methods and no held-out gene was eligible as a landmark. All
124,750 approximate pair p-values were finite. For hyper selection, replacing
the held-out approximate p-values below 1e-3 by their double results reduced
absolute log10-p error to a 95th percentile of 0.0001428 and maximum 0.001840.
All four reference trace matrices retained rank 100 at tolerance 1e-6.

The R = 100 approximate tail contains 52,594 pairs (42.16%) for hyper and
random selection, and 52,595 for score selection. Under comparable cost per
double pair, rechecks alone therefore bound pair-compute speedup near
2.37-fold even before landmark and diagnostic work. This is a cost bound
based on the observed sample fraction, not a measured total-run speedup or
an established fraction for all 11,184 genes. R = 50/200 and the exhaustive
double baseline are still being computed.

A separate closest-gene check found that genes 32 and 229 differed by less
than 0.35% in two reported statistics, while their relative M difference was
0.71468753. This is another precheck showing that nearby summary statistics
alone do not establish interchangeable curvature matrices.

### Earlier reference-trace prototype results

The following values are retained from the earlier reference-trace prototype
and its validation snapshots. They are historical development checks, not
measurements of the current real-gene CUR implementation.

The prior 97,830-observation MAGIC three-gene fit had q = 1962 and common rank
r = 1407. All three exact-pipeline p-values agree bit-for-bit with direct
double native Liu evaluation. A fresh run built three states; a resumed run
built zero and reused all three pair results. This is a correctness and
resumption check, not a large-G throughput benchmark.

For the 635-gene historical data, 64 fixed-seed random references were used.
The 10,000 validation pairs exclude every reference gene at both endpoints.
Exact comparisons use the same float32 checkpoint matrices as reconstruction,
so this table isolates reconstruction error from storage error.

| Route | Quantity | Median relative error | 95th percentile | Maximum |
|---|---|---:|---:|---:|
| INLA, common r = 29 | t1 | 0.000645% | 0.007842% | 0.114715% |
| INLA, common r = 29 | t2 | 0.001394% | 0.015129% | 0.282194% |
| INLA, common r = 29 | t3 | 0.002306% | 0.024377% | 0.368986% |
| INLA, common r = 29 | t4 | 0.003166% | 0.032566% | 0.476483% |
| mgcv/BAM | t1 | 0.006741% | 0.051207% | 0.469147% |
| mgcv/BAM | t2 | 0.010502% | 0.083117% | 0.878289% |
| mgcv/BAM | t3 | 0.012501% | 0.100536% | 1.144822% |
| mgcv/BAM | t4 | 0.013757% | 0.111705% | 1.302537% |

Absolute log10-p error (median / 95th percentile / maximum) was
0.000003715 / 0.000077630 / 0.002868027 for reduced INLA, and
0.000066324 / 0.000987213 / 0.049684075 for mgcv/BAM. All 201,295 reconstructed
p-values in each run were finite. These are empirical errors on these data,
not uniform error bounds or a proof of FDR validity for approximate p-values.

The 635-gene INLA rank of 29 is substantially smaller than MAGIC's rank of
1407. Its timing cannot establish high-dimensional 3D speed or accuracy.
On the three-gene MAGIC storage check, the maximum relative float32 trace
error across all four powers was 5.24e-9. This is likewise a small observed
example, not a universal bound on roundoff or tail p-values.

The updated `inst/benchmarks/pair-trace-approximation.R` benchmarks the current
three selectors at R = 50, 100, and 200 with double holdout and tail checks;
the historical values above are preserved in the earlier validation artifacts.
Local validation artifacts are under `artifacts/pair-pipeline-validation/`,
including held-out trace/p-value data, reference-count comparisons and the
separate MAGIC double-versus-float32 check. The current three-selection CUR
error and end-to-end timing comparison on the supplied 3D fit remain in progress.

## Supplied HPC fit and measured parallelism

The supplied `inlaST-estimate.rds` contains 11,184 available genes and 97,818
observations, with q = 1962 and common 0.995 rank r = 1404. The input file is
unchanged. Working error and working variance are retained as observation-by-gene
matrices; their combined approximately 16.3 GiB is useful model state.

The following historical local timings use the supplied fit, optimized native
kernels, and single-threaded inner BLAS. Preparation includes both sparse units and
reduced covariance materialization. Trace timing excludes reading states and
computing the common projection.

| Work | One thread | Twenty threads | Observed speedup |
|---|---:|---:|---:|
| Prepare 20 gene states | 21.09 s | 3.05 s | 6.91x |
| Four trace powers for 32 pairs | 14.68 s | 1.44 s | 10.19x |

All compared scores, matrices and traces were bit-identical across thread
counts. This confirms that useful batches reached OpenMP in that snapshot. It
does not measure end-to-end performance for all gene pairs or the current CUR
implementation. The successful `R CMD check --no-manual` result applies only
to that earlier source snapshot; the current branch has not completed a full
package check.

The current separate class-mean experiment compares means of covariance
matrices with covariance matrices constructed from mean fitted working fields.
It learns arbitrary symmetric matrix directions within the unchanged common
1404-dimensional coordinate system. It is distinct from the reference-trace
CUR implementation above and is not yet an exported basis-learning API.
