# INLA pairwise observation-kernel projection

Version 0.0.1.9014 adds a score-only approximation for sparse INLA Liu pair
tests. The basis is formed lazily from the constrained observation kernel and
retains the smallest number of directions whose cumulative eigenvalue sum is
at least 0.995. The native observation-mean constraint, feature-specific
working state, sparse Woodbury systems, and expected-curvature nuisance
correction remain those of the full sparse fit.

For each feature, the direct probe constructs its score and curvature only in
the retained coordinates. The pair evaluator keeps a bounded reduced-state
cache and records projection metadata and phase timing in
`timing$inla_projection`. The full-r internal route is retained for arithmetic
comparison; it is not a public option.

The current MAGIC 3D check used Snap25, Foxp1, and Tfap2b at 97,830
observations with q = 1,962. The 0.995 basis retained r = 1,407 directions
(coverage 0.9950175559; tail 0.0049824441). Direct reduced scores and moments
agreed with the matching projected full state. Against full-q Liu, the three
log10-p differences were -0.0582411, 0.0381902, and -0.161442. The initial
three-pair public call took 10.92 s, including 6.72 s for basis construction;
the reduced four-trace kernel took 1.43 s and the full-q kernel took 3.86 s.
Thus this small run measures a reduced trace kernel, not an end-to-end speedup.
