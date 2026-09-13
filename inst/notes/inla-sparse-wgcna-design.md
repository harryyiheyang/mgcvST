> Historical design audit at commit `96e5fca`. The implementation contract was subsequently narrowed to INLA-only, exact Liu, C++ sparse arithmetic and OpenMP. BiocParallel and approximate-calibration proposals below are not the selected implementation. See [the implemented route](inla-openmp.md). mgcv calculation paths are unchanged.

# Sparse INLA to WGCNA design audit

## Current data path

The current sparse path preserves the observation matrix and precision matrix as sparse objects through the feature score solve. In `R/inla-sparse-score.R:33-35`, `A` and `Q` are stored as `CsparseMatrix` objects. For each feature, `R/inla-sparse-score.R:104-124` forms `K=A'WA`, factors `tau*Q+K` with sparse Cholesky, solves jointly for the nuisance columns and score right-hand side, and returns the score vector `a`. WGCNA requests `score_only=TRUE` at `R/wgcna.R:135-138`, so it does not form the per-feature `KT`, `SKT`, or dense calibration matrix `M` at `R/inla-sparse-score.R:128-135`.

Two qualifications are material. First, the canonical source at this audit constructs a dense projected coefficient factor: `R/inla-sparse-score.R:36-40` converts `Z` and `Q%*%Z` to dense matrices, applies dense Cholesky, and stores dense `Tbase`. Thus the fitted score geometry is sparse in observation space but not fully sparse in coefficient coordinates. The constrained `sqrt(P) Q sqrt(P)` statistic must remain unchanged; only an algebraically equivalent sparse representation may replace this dense arithmetic. Second, `R/wgcna.R:120-121` calls `.mgcvst_model_fixed_factors()` before dispatching to the sparse branch. This is computationally redundant for sparse fits, although `R/model-score.R:24-27` currently returns only an empty list structure and does not factor or densify `Q`.

WGCNA then stores all selected score vectors in a dense coordinate-by-gene matrix at `R/wgcna.R:157-162`. For each requested gene block, `R/wgcna.R:273-290` constructs dense covariance, correlation, adjacency, TOM, and clustering objects. These are gene-by-gene matrices. Their density is part of the stated WGCNA calculation and is distinct from the spatial (n\times q) projection and (q\times q) precision algebra.

The current WGCNA covariance is the uncentered Gram matrix `crossprod(A)/nrow(A)`, followed by `cov2cor`; the documentation explicitly says no centering is applied. An orthogonal change of aligned score coordinates preserves this Gram matrix. If a future analysis instead requests ordinary Pearson correlation after centering across score coordinates, inner products alone are insufficient unless coordinate sums are also retained, and a whitening rotation need not preserve those sums. That alternative must be treated as a change of estimand rather than an implementation optimization.

## Pair-test duplication

The WGCNA path does not evaluate pairs and performs no per-pair decomposition. It constructs each selected feature score once and obtains all within-block similarities by one crossproduct. The expensive pair behavior belongs to `mgcvST.test()`.

For model-based pair tests, `R/model-test.R:20-26` creates a feature-state cache inside each pair chunk. The cache prevents repeated construction within one chunk, but the same feature is recomputed in every chunk or worker that contains one of its pairs. A sparse full state contains a score vector and a dense (q\times q) `M`; at (q=1,962), one double matrix is about 29.4 MiB. Forty cached feature matrices require about 1.15 GiB before R object overhead in each worker. Sending the complete fit to multiple SOCK workers and recreating overlapping caches can multiply both this memory and the feature-specific sparse Cholesky work.

Davies calibration at `R/model-score.R:214-218` consumes two feature matrices for every pair and performs the calibration's spectral work per pair. The optimized Liu path in `R/mgcvst-api.R:897-1047` constructs one feature summary, stores all dense `H` matrices, and evaluates trace powers in bounded C++ pair blocks. `src/score_pair.cpp:56-73` distributes pairs with OpenMP, but each thread creates dense `product=left*right` and, for higher moments, `product2`. At (q=1,962), each temporary product is another approximately 29.4 MiB per active OpenMP thread, in addition to the feature matrices. This path avoids repeated eigendecomposition, but it remains cubic coefficient-space matrix multiplication for each requested pair.

## Exact scaling design under the current INLA-only requirement

The implementation applies only to INLA fits. The mgcv fitting, marginal, pair, and WGCNA paths remain unchanged. The INLA implementation should separate three stages with explicit resource controls:

1. **Feature score stage.** Construct each feature state once in bounded C++ blocks. For WGCNA, retain only `a`. For marginal tests, retain the statistic and four exact Liu trace moments. For pair tests, retain `a` and the exact feature matrix or factor required by the current trace identities. Shared sparse `A` and `Q` remain resident and feature-specific states are written in original feature order.
2. **Pair stage.** Reuse the feature summaries. Process the explicit pair universe in bounded C++ OpenMP blocks and return result rows immediately. INLA Davies is unavailable by design; Liu is the sole calibration. Only one OpenMP layer may be active, with BLAS and INLA numerical threads fixed to one during downstream scoring.
3. **Network stage.** Assemble the exact score Gram matrix for each user-specified block, then run the existing adjacency, TOM, clustering, and tree cut. Independent named blocks may be processed sequentially under a memory budget. One large block remains one exact network job because covariance, adjacency, and TOM each require the complete within-block matrix.

The downstream INLA interface should expose `threads`, `feature_chunk_size`, and `pair_chunk_size`, and should use C++ sparse/OpenMP kernels without BiocParallel. OpenMP's `OMP_NUM_THREADS` controls threads inside a process; `OMP_MAX_ACTIVE_LEVELS=1`, BLAS thread limits of one, and explicit `num_threads` arguments prevent nested oversubscription. Earlier BiocParallel measurements remain execution evidence for the existing public code, but they are not the design for this INLA-only extension.

The applicable runtime reference is the [OpenMP `OMP_NUM_THREADS` specification](https://www.openmp.org/spec-html/5.1/openmpse59.html).

## Memory limits

Let (p) be the number of selected genes, (q_s) the score-coordinate count, and (b) the size of one WGCNA block. The dense score matrix requires approximately `8*q_s*p` bytes. Each retained feature calibration matrix requires approximately `8*q_s^2` bytes. These feature matrices dominate pair-test memory and should be bounded by a cache budget, with completed summaries written to compact immutable shards when all cannot remain resident.

For a WGCNA block, each dense double matrix requires `8*b^2` bytes: 7.45 MiB at 1,000 genes, 190.7 MiB at 5,000 genes, and 762.9 MiB at 10,000 genes. The current returned object retains covariance, correlation, adjacency, and TOM, requiring at least `32*b^2` bytes before temporary matrices, distance objects, names, and clustering state. A conservative execution budget should reserve at least six to eight `b^2` double matrices during TOM construction and reject a requested block before allocation when that estimate exceeds a user-set memory limit.

Exact user-defined blocks are already supported and analyzed separately. They are the scaling boundary because they preserve the current declared semantics. Automatically splitting one requested block, using WGCNA `blockwiseModules`, sparsifying similarities, or thresholding adjacency would remove between-block relationships or alter TOM and module labels. These approximations are excluded from the current implementation.

## Recommended order

The minimal safe changes are: bypass `.mgcvst_model_fixed_factors()` in the INLA WGCNA sparse branch; preserve the constrained `sqrt(P) Q sqrt(P)` score definition while changing only its arithmetic representation; add a C++ score-only feature-block stage whose output order is checked against feature IDs; and add an explicit network memory preflight. Pair-test optimization should move caching from pair chunks to a single feature-summary stage and use bounded C++ OpenMP blocks. Equality tests must compare score vectors, Liu trace moments and p-values, score Gram matrices, adjacency, TOM, and module labels against the current serial INLA calculation. The mgcv implementation is outside this change.

## Marginal Liu trace identities

The frozen marginal TAPS calculation defines `Theta` as the generalized inverse of the scaled tested penalty and `C=B'PB`. Its current mixture matrix is `H=Theta^(1/2) C Theta^(1/2)`. Liu calibration needs only `c_k=sum(lambda^k)` for `k=1,...,4`. By cyclic invariance,

`c_k = trace(H^k) = trace((Theta C)^k)`.

Thus INLA can compute `G=Theta C`, multiply `G` successively, and record `trace(G)`, `trace(G^2)`, `trace(G^3)`, and `trace(G^4)` without constructing a square root or eigendecomposition. An equivalent factor form uses a valid factor `L` of `Theta` and `H=L' C L`, then takes the same four traces. The marginal statistic remains `max(0, r' B Theta B' r)`.

The existing implementation first truncates nonpositive eigenvalues of `Theta` to zero and later retains only mixture eigenvalues greater than `1e-15`. Direct trace powers include any small positive or negative numerical components that this spectrum path discards. The equivalence test must therefore report both the four moment differences and the sums of discarded eigenvalue powers. Production code must reject non-finite or invalid Liu moments explicitly; it must not silently apply absolute values, threshold the score matrix, or clamp moments to make them positive. A small equivalence fixture should compare the statistic, all four moments, and final Liu p-value, while recording the old truncation contribution as the numerical boundary.
