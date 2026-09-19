> Historical design audit at commit `96e5fca`. `graphical_susie()` referenced below was deleted from the package after this audit; only `mgcvST.wgcna()`/`inlaST.wgcna()` remain as sparse-score downstream consumers. The implementation contract was subsequently narrowed to INLA-only, exact Liu, C++ sparse arithmetic and OpenMP. BiocParallel and approximate-calibration proposals below are not the selected implementation. See [the implemented route](inla-openmp.md). mgcv calculation paths are unchanged.

# Sparse INLA downstream parallel design on Windows

## Scope and conclusion

This audit covered a former sparse INLA score path used by marginal tests,
pairwise tests, `mgcvST.wgcna()` and `graphical_susie()`. The historical safe
Windows unit of parallelism is an independent R process created by
`BiocParallel::SnowParam(type = "SOCK")`, with one numerical thread inside
each process. Sparse `A` and `Q` matrices may be serialized to those processes.
A numerical Cholesky factor should be constructed and reused inside one worker;
it should not be treated as a concurrently shared factor across R processes or
OpenMP threads.

Sparse storage does not imply a sparse factor of comparable size. Ordering and
the graph of `tau * Q + A'WA` determine Cholesky fill-in. Resource planning must
measure the factor nonzeros and worker peak memory, in addition to `nnzero(Q)`.
That former route was also not sparse end to end: `score_sparse$coefficient_factor`
is the dense projected `Tbase`, and a full pair state contains dense `M`.

## Current implementation

`inlaST.estimate()` sends feature chunks, the response matrix, and one shared
model specification through `bplapply()`. On Windows, every SOCK worker is a
separate R process. The model specification contains `A` and `Q` as compressed
sparse Matrix objects. They remain sparse during serialization, although each
worker owns its deserialized copy.

For pairwise sparse scores, `.mgcvst_test_model()` sends one compact fit as a
common argument and defaults to about one pair chunk per worker. Each worker
creates a private state cache. For each feature first encountered in that
chunk, `.mgcvst_model_sparse_score_state()` constructs

```
H = tau * Q + A' W A
```

and factors `H` with `Matrix::Cholesky(..., LDL = FALSE, super = FALSE)`. The
factor is captured by a private solver closure only while that score state is
being assembled. The cached downstream state contains dense low-rank `a` and
`M`, rather than the CHOLMOD factor. A feature repeated in the same chunk
therefore reuses `a` and `M`; a feature appearing in different workers is
factored independently. Many small chunks repeat those feature factorizations.
The large-component validation illustrated this cost: its 41-chunk BAM run
cannot be interpreted as an estimator-only speed comparison with a four-chunk
run.

`mgcvST.marginal()` distributes retained per-feature states and shared geometry
in feature chunks. Each worker is forced to one numerical thread. Liu spectral
moments are then summed by a separate OpenMP kernel in the manager process;
OpenMP is not nested inside the BiocParallel workers.

`mgcvST.WGCNA()` obtains sparse scores serially with `score_only = TRUE`. This
still forms and factors sparse `H`, but returns before constructing `KT`, `SKT`
and the dense pair-calibration matrix `M`. WGCNA receives the resulting compact
score matrix. `graphical_susie()` similarly constructs score states before its
parallel nodewise initialization; its `bplapply()` sends the already formed
dense covariance `S`, not sparse factors.

The package C++ OpenMP kernels do not factor sparse matrices. `score_pair.cpp`
works on pointers to dense R matrices, and `marginal_liu.cpp` sums dense
spectral-power columns. Sparse Matrix/CHOLMOD operations occur in ordinary R
code outside those OpenMP regions. `.mgcvst_thread_limit()` sets OpenMP, BLAS,
RcppParallel and data.table thread counts to one in workers. Thus the current
process-parallel path avoids concurrent calls on one CHOLMOD factor and avoids
worker-by-BLAS oversubscription.

The `threads` argument accepted by the model-set pair path is currently
validated but is not forwarded into its sparse factorization or calibration.
Its effective parallelism is the number of BiocParallel workers. This should be
reported accurately rather than described as nested OpenMP acceleration.

## Serialization and factor ownership

`dgCMatrix` and `dsCMatrix` objects have ordinary compressed-column slots and
are suitable transport objects. Matrix documents `CHMfactor` as an R wrapper
around a CHOLMOD factor structure. Even when a factor can be serialized and
successfully reconstructed, that behavior does not make one factor a shared,
concurrent object. SOCK serialization creates an independent copy in each
process; OpenMP threads in one process would instead touch the same object.

The safe contract is therefore:

1. Serialize sparse matrices and plain numeric metadata.
2. Build symbolic and numerical factor state inside the worker that owns it.
3. Reuse that state only sequentially within that worker.
4. Return plain R matrices, vectors and diagnostics to the manager.
5. Do not put a solver closure or CHOLMOD factor into a payload intended for
   concurrent use.

The SuiteSparse build documentation also notes that CHOLMOD performance can
depend on OpenMP and parallel BLAS. That is a build-level property, not evidence
that a single factor object may safely be called concurrently from package
OpenMP code. The conservative package design keeps CHOLMOD calls process-local
and single-threaded when several SOCK workers are active.

## Recommended Windows architecture

Use a persistent SOCK cluster for a complete downstream stage. Partition by
feature ownership rather than by arbitrary pair rows:

- Send the shared sparse geometry once to each persistent worker.
- Assign each feature to one owner worker and build its sparse score state once.
- Return compact `a`, `M`, width, convergence and factor diagnostics. For WGCNA,
  request only `a` and width.
- Run pair calibration on compact states. If pair calculation remains
  distributed, send compact states needed by a coarse pair block; avoid sending
  the full fitted object with every small task.
- Apply BH adjustment once in the manager over the complete set of valid p
  values.
- Keep `OMP_NUM_THREADS`, BLAS threads and INLA threads at one whenever more
  than one SOCK worker is used.

Returning every full dense `M` to one manager is not a viable large-data
default. At score dimension 3,000, one double-precision `M` occupies about 72
MB; 10,000 such matrices require about 720 GB before list and allocator
overhead. Full-state work therefore needs bounded feature/pair blocks with
explicit eviction, disk-backed checkpoints, or an operator representation.
WGCNA should continue using score-only states. The choice among bounded dense
states and a new operator design requires separate numerical implementation
work; sparse `Q` alone does not solve this capacity limit.

For the present implementation, the lowest-risk tuning is a small number of
large chunks, normally one chunk per worker. This preserves the tested
statistical calculation and increases cache reuse without introducing shared
factor state. A future persistent worker cache can key symbolic structure by a
hash of the `A`, `Q`, constraint and nuisance geometry, and numerical state by
feature ID plus `tau` and working-weight identity. Cache invalidation must be
explicit because `W` and `tau` change `H`.

The design should expose these measurements per worker: serialized payload
bytes, `nnzero(A)`, `nnzero(Q)`, `nnzero(H)`, factor `colcount` sum or expanded
factor nonzeros, factorization count, cache hits, peak RSS and elapsed time.
This separates transport cost, fill-in, and repeated factorization.

## Small validation plan

No additional long benchmark is needed. A package test or short benchmark can
use the existing ten-gene real component fixture and compare `SerialParam()`,
`SnowParam(2)` and `SnowParam(4)` with one thread per process:

1. Assert worker PID, package version, library path and all thread-control
   variables.
2. Assert `A`, `Q` and worker-received geometry retain sparse classes and the
   same dimensions, nonzero counts and numeric values.
3. Count sparse factorizations and cache hits. With one coarse chunk per worker,
   each feature should be factored at most once in each worker that needs it.
4. Compare `a`, width, `M`, signed score, information and all 45 pair p values
   against serial results within recorded numerical tolerances.
5. Compare WGCNA score-only `a` and width with the full sparse state, and assert
   that score-only does not construct `M`.
6. Record `nnzero(H)`, factor `colcount`, serialized sizes and peak process RSS.
7. Start and stop the persistent cluster explicitly, then repeat a small call
   to confirm that no stale cache crosses a new model identity.

The completed public-path ten-gene check already supplies the numerical
Serial/Snow2/Snow4 comparison. The additional test should focus on sparse class,
factorization-count and memory assertions; it should not repeat the large
component run whose 7.3-minute burden the user has accepted.

This audit did not execute a new serialization or factor-concurrency experiment.
Statements about the present data flow and cache contents come from source-code
inspection; the proposed ownership and capacity rules are design recommendations
grounded in that inspection and the cited library documentation.

## Primary references

- Bioconductor, [BiocParallel reference manual](https://bioconductor.org/packages/release/bioc/manuals/BiocParallel/man/BiocParallel.pdf): `SnowParam` uses distributed-memory workers and documents worker/task controls and export behavior.
- Bioconductor, [Random Numbers in BiocParallel](https://bioconductor.org/packages/release/bioc/vignettes/BiocParallel/inst/doc/Random_Numbers.html): SOCK workers are independent R processes and results can be invariant to worker/task count under the documented RNG contract.
- Matrix maintainers, [Matrix reference manual](https://stat.ethz.ch/CRAN/web/packages/Matrix/refman/Matrix.html): compressed sparse classes, sparse solves and the `CHMfactor` representation backed by CHOLMOD.
- SuiteSparse maintainers, [SuiteSparse build and threading documentation](https://github.com/DrTimothyAldenDavis/SuiteSparse/blob/dev/README.md): CHOLMOD is the sparse Cholesky component and its performance can depend on the SuiteSparse OpenMP and BLAS build.
