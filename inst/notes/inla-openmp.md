# INLA sparse OpenMP arithmetic

> Historical implementation note. This describes the pre-0.0.1.9014 exact
> full-q downstream route. Current INLA Liu pair testing uses the constrained
> observation-kernel score approximation described in
> `inla-pairwise-projection.md`; its full sparse fit remains unchanged.

The INLA single-global fixed-kappa SPDE downstream path uses the same expected-curvature score statistics and exact Liu moments, implemented with sparse Eigen solves and OpenMP. The test congruence (square root of residual precision, kernel, square root of residual precision) is preserved algebraically; it is not a removable projection. No stochastic trace estimator or new calibration is introduced. The mgcv estimators and numerical kernels are unchanged.

INLA rejects Davies and non-Serial BiocParallel parameters. `threads` controls downstream OpenMP; fitting uses INLA's own `control$num_threads`. The current sparse capability remains one constrained SPDE target with fixed-effect nuisance terms. Unsupported INLA downstream geometries raise an explicit error rather than falling through to an unintended dense or process-parallel implementation.

## Computation and reusable units

Sparse A and Q remain sparse through weighted crossproducts and factorization. The shared Q factor is prepared once and read concurrently through const Eigen operations. Feature-specific H factors are constructed once per feature in a pair-test call. Each OpenMP task owns its working matrices. Dense curvature is constructed in bounded right-hand-side blocks, only when required for exact trace products. No observation-by-observation matrix is constructed and no eigendecomposition is used by the new Liu kernels.

Pairwise reconstruction units contain the score vector, scale, sparse weighted crossproduct, sparse H triangular factor, diagonal and permutation, constraint-solve vector and denominator, and small nuisance quantities. Units contain neither a full curvature matrix M nor observation-length response/weight vectors. Small collections remain in memory; larger collections are temporarily serialized as units, with cleanup on return. M is materialized only for the active pair block and is never written to disk. Reading a unit requires no H refactorization. Repeated pair blocks may materialize M again; the implementation does not claim to eliminate that cost.

The shared Q cache is an in-process external pointer. After serializing and restoring a complete fit in a new R session, it is rebuilt once from Q; it is never treated as a portable pointer. Unit factors themselves are ordinary serializable sparse matrices and vectors.

Marginal calls retain the already computed result for replay. Marginal matrix storage is released after the four moments are computed. WGCNA requests only score vectors and preserves the previous coordinate-rank normalization in its Gram matrix.

## Small verification

Tests compare dense expected-curvature references with the C++ route for scores, Gram matrices, four trace moments and Liu p-values. They also compare one and two OpenMP threads and the materialized states after a serialization round trip. An additional two-gene NB API smoke test uses 30 observations with the flat spatial and NB-size priors. No type-I simulation campaign is repeated for this arithmetic change.

A separate 256-node, two-feature arithmetic fixture stores 61,213 bytes of reconstruction units, versus 1,049,008 bytes for the two dense M objects (about 17 times smaller). Reading and materializing this small fixture were below the elapsed timer resolution; this is not a large-mesh timing claim. Sparse factor fill-in determines the storage and reconstruction cost on a real mesh.
