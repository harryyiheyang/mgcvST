# Sparse raw-kernel score plan

## Scope

This note audits the pairwise covariance score for one fixed-kappa SPDE
component.  The fitted null field always satisfies the observation-weighted
constraint

\[
g^T u=0,\qquad g=A^T\mathbf 1/n,
\]

but the tested kernel is the raw

\[
G=AQ^{-1}A^T/\tau.
\]

The derivation does not identify the conditional coefficient covariance from
the INLA fit with either `G` or the constrained null covariance.  INLA supplies
the empirical-Bayes working variance and precision `tau`; the score operator is
then reconstructed from the expected-Fisher working model.

## Mesh-coordinate derivation

Write `D` for the diagonal working covariance and

\[
C_g=Q^{-1}-Q^{-1}g(g^TQ^{-1}g)^{-1}g^TQ^{-1}.
\]

The fitted-null covariance is `Vc = D + A Cg A' / tau`.  For

\[
S=A^TD^{-1}A,\qquad H=\tau Q+S,
\]

define the constrained inverse action

\[
H_c^{-1}B=H^{-1}B-H^{-1}g
 (g^TH^{-1}g)^{-1}g^TH^{-1}B.
\]

Then Woodbury gives, without constructing `Vc`,

\[
W y=V_c^{-1}y=D^{-1}y-D^{-1}A H_c^{-1}A^TD^{-1}y.
\]

For a nuisance design `X`, set

\[
T=A^TD^{-1}X,\quad U=X^TD^{-1}X,
\]

\[
A^TWA=S-SH_c^{-1}S,\quad
A^TWX=T-SH_c^{-1}T,
\]

\[
X^TWX=U-T^TH_c^{-1}T,\quad
V_N=(X^TWX)^{-1}.
\]

The correct fitted-null residual precision is

\[
P=W-WXV_NX^TW.
\]

For a working error `e`, the mesh-coordinate score and information are

\[
b=A^TPe,\qquad K=A^TPA.
\]

Both follow from the quantities above using only sparse products, constrained
solves, and rank-`p` nuisance corrections.

The raw kernel has a large near-intercept direction on the current FEM mesh.
Because the correct `P` obeys `PX=0`, it is algebraically safe to replace the
score basis by

\[
A_0=A-X(X^TX)^{-1}X^TA.
\]

The implementation never materializes `A0`: the matrix
`E=(X'X)^-1 X'A`, of size `p` by `m`, supplies the required low-rank updates to
`A'WA`, `A'WX`, and `A'We`.  This is the same stabilization used by the frozen
constraint simulation, and it avoids changing the deliberately invalid
keep-nuisance reference variant.

Let `Q = L L'` be an unpermuted sparse Cholesky factorization.  An aligned raw
score factor would be `F=A0 L^-T/sqrt(tau)`, but it is not formed.  Instead,

\[
a=L^{-1}b/\sqrt{\tau},\qquad
M=L^{-1}K L^{-T}/\tau.
\]

Thus for two features the observed signed score is `crossprod(a1,a2)`, and the
first four Liu moments are `trace((M1 M2)^k)`, `k=1,...,4`.  Davies uses the
same two mesh-space `M` matrices.  These are exactly the nonzero spectral
quantities associated with `sqrt(P) G sqrt(P)`; computing an observation-space
matrix square root is unnecessary.

## Storage and dense work that remain

The following objects stay sparse: `A`, `Q`, `D^-1/2 A`, `S`, `H`, and their
Cholesky factors.  Constraint application adds one mesh vector and scalar
Schur complement.  Nuisance removal adds `m` by `p`, `p` by `p`, and rank-`p`
updates, where `p` is normally small.

Exact use of the current calibration API still produces one dense `m` by `m`
matrix `M` per active feature.  Forming `K` requires applying the constrained
`H` solve to the `m` columns of `S`; the result is generally dense even though
the FEM inputs are sparse.  The `p` by `p` nuisance covariance is also dense.
These dense matrices are mesh-space, never observation-space.  Avoiding the
`m` by `m` result would require a new operator-valued or approximate trace
calibrator; it is not a consequence of merely rewriting `sqrt(P)Gsqrt(P)`.

For many Liu pairs, retaining every dense `M` remains the main memory cost.
An exact redesign could expose coefficient-space solve operators to a revised
trace-power backend, but products such as `trace((M1 M2)^k)` still require
substantial mesh-space work.  Randomized trace estimates would change the
calibration numerically and need separate statistical validation.

## Current package call sites

`mgcvST.test()` dispatches through `.mgcvst_model_pair_single()`.  Each feature
is summarized by `.mgcvst_model_score_state()`, which calls
`.mgcvst_model_operator_vp()`.  That path currently constructs projected
observation-by-innovation factors from `geometry$B` and its projected penalty,
then computes `a` and `M` through `.mgcvst_model_apply_P()`.

The older high-throughput path in `.mgcvst_liu_summaries()` likewise constructs
`T0 = B Q^-1/2` and stores dense per-feature `H` matrices before
`mgcvst_pair_trace_powers_cpp()` evaluates pair moments.  The experimental raw
implementation in `inst/benchmarks/inla-raw-kernel-sparse.R` can replace the
single-component state construction conceptually, but it is not wired into
either production path and does not alter the frozen type-I scripts.

The INLA adapter already retains the raw sparse `A`, `Q`, and exact constraint
inside `model$inla_spec$random`.  The fitted object currently exposes projected
score geometry to the common `mgcvST.test()` path.  A production raw-kernel
option would therefore need an explicit fit contract carrying the raw matrices
or a compact immutable raw-score specification; silently substituting raw
geometry into the present projected fields would mix two test definitions.

## Runnable reference

`inlast_raw_kernel_sparse_state()` implements the derivation, and
`inlast_raw_kernel_sparse_pair()` passes its `a` and `M` to the existing
calibrator.  The intended reference check uses the same `raw_A`, `raw_Q`, `g`,
working error, working variance, `tau`, and nuisance design with
`inlast_constraint_reference_states(...)$states$raw_kernel_only`.  Compare
`a`, `M`, signed pair score, and the four trace moments.  The check is numerical
equivalence only; no timing claim should be inferred while the type-I
simulation is sharing the machine.

This check was run on the saved `nb3_pair_k07` pilot geometry with `n=200`,
`m=36`, two independently fitted NB features, and the actual sparse `A`, `Q`
and fitted working states.  Relative errors against the factor-based reference
were at most `1.32e-12` for `a`, `6.26e-12` for `M`, `7.60e-13` for the signed
pair score, and `9.32e-14` for the four trace moments.  The two-sided Liu
p-values differed by `2.14e-13`.  The sparse constrained-solve residual was at
most `1.36e-16`; the minimum `M` eigenvalues were `1.17e-17` and `-1.19e-17`,
the latter being roundoff at the scale of the matrices.  No runtime comparison
was recorded.

A reproducible two-geometry check is in
`inst/benchmarks/inla-raw-kernel-sparse-check.R`; it writes
`artifacts/constraint-type1/raw-kernel-sparse-check.csv`.  It performs no INLA
fit.  The first case uses the actual PathwayLGM 151673 geometry (`n=3611`,
`m=625`, `kappa=0.7`) with deterministic positive working variances, errors,
and two precisions.  Its maximum relative errors were `3.55e-15` for `a`,
`3.47e-16` for `M`, `3.93e-15` for the signed score, and `2.22e-16` for the
trace moments.  The raw constraint direction's residual norm after intercept
projection was `0.121759` of its total norm; despite the dominant intercept
component, the two `M` minima were only `-1.28e-16` and `-6.01e-17`.

The second case uses the simulation geometry (`n=200`, `m=36`, `kappa=6`)
with a nontrivial `intercept + x` nuisance design.  All relative errors were at
most `1.39e-15`, including the nuisance covariance.  Its most negative `M`
eigenvalue was `-3.91e-16`.  The CSV records a largest dense dimension of `m`
and `observation_square_matrix_constructed = FALSE` for every row.  The dense
reference also avoids observation-square matrices, though it deliberately
forms observation-by-mesh factors for an independent comparison.
