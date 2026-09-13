# Sparse expected-curvature score design

## Current implementation boundary

The public INLA sparse backend is not sparse from end to end. In
`R/inla-sparse-score.R`, `.inlast_sparse_score_geometry()` converts the
constraint projection to a dense matrix, forms the dense projected precision
`Z'QZ`, and stores its dense inverse-Cholesky factor (`lines 35--43`). Each
feature then uses this factor to construct the score coordinates and the dense
calibration matrix (`lines 89--115`). The same feature calculation uses
`fit$nuisance_covariance[[feature]]`, which is the native INLA conditional
covariance, although `inlaST.estimate()` also retains
`expected_nuisance_covariance` (`R/inla-api.R`, lines 562--564 and 681--682).
Thus the current pairwise and WGCNA sparse entry points inherit both the dense
projection and the native/expected curvature mismatch.

The retained marginal route is a separate dense GAM representation. It forms
dense penalty inverses, dense nuisance covariance matrices, and a dense
eigenproblem in `R/marginal-taps.R` (especially lines 68--72, 101--123, and
139--145). It therefore cannot be the large-data INLA marginal implementation.
Pairwise tests call `.mgcvst_model_sparse_score_state()` through
`R/model-score.R`, lines 181--220. WGCNA calls the score-only form of the same
function (`R/wgcna.R`, lines 135--138), but that form still multiplies by the
stored dense coefficient factor before returning.

## Projection-free expected-curvature state

Let `A` be the sparse observation matrix, `Q` the sparse fixed-kappa SPDE
precision, `g'x=0` the active linear constraint (the observation-mean or
volume constraint supplied by the model), `W=D^-1`, and `X` the fixed/nuisance
design. For one feature define

`K=A'WA`, `L=A'WX`, `t=A'We`, and `H=tau*Q+K`.

All applications of the constrained inverse of `H` can use one sparse
Cholesky factor and the rank-one constraint formula already implemented in
`.mgcvst_model_sparse_constrained_solver()`. With `S` denoting this operator,
the expected nuisance information and covariance are

`J=X'WX-L'SL` and `Vp_expected=J^-1`.

Only `J` and `Vp_expected` must be dense; their dimension is the nuisance
coefficient count, rather than the number of observations or mesh nodes. The
feature score state can then retain

`U=L-KSL`, `q=X'We-L'St`, and
`h=t-KSt-U Vp_expected q`.

The expected covariance operator for `h` is

`B(v)=Kv-K S(Kv)-U Vp_expected U'v`.

This expression requires sparse products, constrained sparse solves, and a
small dense nuisance adjustment. It does not require `Q^-1`, `H^-1`, an
observation-by-observation matrix, or a projected precision.

If `m` is the mesh size and `p_x` is the nuisance dimension, each feature uses
one sparse factorization of `H`, `p_x+1` initial right-hand sides, sparse
matrix-vector products, and `O(p_x^3)` dense work for `J`. Memory is governed
by the Cholesky fill of `H` plus `O(m*p_x+p_x^2)`. Applying `B` thereafter uses
sparse products, one constrained `H` solve, and the small nuisance matrices.
No step has `O(n^2)` storage. The cost of a `Q` application is similarly a
sparse triangular solve after one factorization shared across features.

The target covariance operator `C` is the constrained inverse of `Q`. For a
proper fixed-kappa `Q`, factor sparse `Q` once and apply the same rank-one
constraint formula using `g`. If a future precision is intrinsically singular,
use a sparse symmetric-indefinite factorization of the augmented KKT matrix
`[Q g; g' 0]`. Converting a dense projection back to `sparseMatrix` does not
meet this contract because fill and memory have already been incurred.

## Marginal and pairwise tests

The expression `h'C h/tau` is a candidate random-field quadratic for a
projection-free marginal test. It is not yet established as algebraically
equivalent to the existing marginal TAPS statistic. The production design must
first derive the tested variance-component score and null law from the current
TAPS definition, then compare statistics and spectra on small exact fixtures.
If that derivation establishes this quadratic, its null mixture spectrum is
the nonnegative spectrum represented by `C^(1/2) B C^(1/2)/tau`, equivalently
the constrained generalized eigenproblem `Bv=lambda*tau*Qv`. Under the selected
expected-curvature policy, negative computed eigenvalues are set to zero. Their
pre-truncation magnitudes may be retained as internal numerical diagnostics,
without adding them to the public result.

For features 1 and 2, the signed pair score is
`h1'C h2/sqrt(tau1*tau2)`. The Liu moments are
`tr[(C B1 C B2/(tau1*tau2))^k]`, for `k=1,...,4`. These formulas use only
applications of `B1`, `B2`, and constrained sparse solves with `Q`.

There is a computational boundary. Exact traces of all four powers and the
complete Davies spectrum are generally full-rank problems even when `Q`, `K`,
and `H` are sparse. They do not mathematically require storing a dense matrix:
exact traces can be accumulated with blocked basis-vector applications and
sparse solves. That route avoids dense storage but can require `m` operator
applications per power and is therefore expensive. Powers of `C B1 C B2` need
not remain sparse if they are materialized.

The first production implementation should preserve the existing exact Liu
semantics by using blocked exact trace accumulation and expose its runtime
cost. Reproducible Hutchinson or stochastic Lanczos trace estimation is a
separate research option; it must not silently replace `calibration="liu"`.
Likewise, restarted Lanczos with a residual-spectrum approximation is a
research option for Davies and must be labelled approximate. Exact Davies
requires the complete relevant spectrum and may be impractical for a large
mesh even without a dense stored operator. An optional exact small-mesh route
may materialize a mesh-by-mesh matrix after an explicit threshold check. No
route should allocate an observation-by-observation matrix.

## WGCNA

WGCNA does not require explicit whitened score coordinates. For selected
features, solve `Q z_j=h_j` subject to `g'z_j=0` and form the gene Gram matrix

`S_ij=h_i'z_j/[q*sqrt(tau_i*tau_j)]`, where `q` is the constrained score
dimension used by the current `crossprod(A)/nrow(A)` definition.

The sparse `Q` factor is shared across genes, and solves can be performed in
bounded right-hand-side blocks. The gene-by-gene covariance, correlation,
adjacency, and TOM matrices are necessarily dense because WGCNA consumes all
gene pairs. This cost is in the selected gene count and is distinct from an
observation-by-observation or mesh-by-mesh dense matrix. The current internal
contract returning `score$A` should therefore gain a sparse-INLA Gram path;
the public WGCNA result can retain the resulting covariance while recording
that coordinate scores were not materialized.

## Proposed internal contracts

- `.inlast_sparse_score_geometry()` retains only sparse `A`, sparse `Q`, `g`,
  the target name, and the smoothing-parameter index. It does not retain `Z`
  or `coefficient_factor`.
- `.mgcvst_sparse_feature_state(fit, feature)` returns `h`, `tau`, an operator
  for `B`, the small `J`/`Vp_expected`, and solve diagnostics.
- `.mgcvst_sparse_q_solver(score_geometry)` factors `Q` once and returns a
  constrained solve operator shared by marginal, pairwise, and WGCNA calls.
- `.mgcvst_sparse_pair_calibrate(state1, state2, q_solver, calibration,
  control)` preserves the existing calibration contract; any future approximate
  route requires a separately named, validated option.
- `.mgcvst_sparse_marginal_calibrate(state, q_solver, calibration, control)`
  uses the same expected-curvature operators and numerical controls.
- `.mgcvst_sparse_wgcna_gram(states, q_solver, block_size)` forms the required
  gene Gram matrix without score coordinates.

Tests should compare these operators with the existing dense formulas on small
fixtures, including the constraint residual, score, first four moments, and
WGCNA Gram matrix. Memory tests should reject stored dense mesh projections,
dense mesh covariance factors, and any `n`-by-`n` allocation in the sparse
fit and downstream paths.

This note is a source review and implementation design. It does not establish
the marginal-score equivalence or select an approximate calibration method;
both require separate derivation and validation before production changes.
