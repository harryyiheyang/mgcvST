# Pairwise Liu execution

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

Davies calibration retains its existing execution route. The checkpoint
arguments and the internal adaptive cache apply to Liu calibration.

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

The internal `cache_bytes = NULL` default probes system headroom, Linux cgroup v1/v2 limits and
Slurm allocations. The minimum available signal controls the automatic
budget. Feature workspaces and pair buffers are reserved separately; probes
are refreshed during preparation and at pair-batch boundaries. A supplied
`cache_bytes` (internal) is a resident-cache ceiling, not a total process-memory limit.
If no memory source is available, the reported conservative fallback is
512 MiB. This fallback is not a hardware assumption.

All q-by-q covariance storage and native arithmetic remain double in the
exact branch. A 150 GB allocation therefore does not imply that 11,099 full
1500-by-1500 double matrices all fit: those matrices alone require about
186 GiB. Disk shards and block scheduling remain useful in this case.

## Exact pipeline validation

The prior 97,830-observation MAGIC three-gene fit had q = 1962 and common rank
r = 1407. All three exact-pipeline p-values agree bit-for-bit with direct
double native Liu evaluation. A fresh run built three states; a resumed run
built zero and reused all three pair results. This is a correctness and
resumption check, not a large-G throughput benchmark.

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
does not measure end-to-end performance for all gene pairs. The successful `R CMD check --no-manual` result applies only
to that earlier source snapshot; the current branch has not completed a full
package check.

The approximate Liu trace path for sparse INLA fits is PCAlearning; see
`pca-learning.md`.
