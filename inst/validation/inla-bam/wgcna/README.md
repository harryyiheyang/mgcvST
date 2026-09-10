# Paired WGCNA validation

The local validation used 40 paired datasets: 10 independent seeds in each
of the `strong_original`, `strong_repeat`, `moderate`, and `independent`
arms. Each dataset contains 90 negative-binomial features at 2,125 locations
with 298 spatial score coordinates. Current bam fREML/discrete and sparse
INLA used the same response matrix. Both INLA hyperpriors were flat.

The frozen installed package was version 0.0.1.9005 from commit
`f2ba73886abe7382b53e8ae9ff82148ebbdccaea`. The checkout commit recorded in
the manifest is separate because source development continued while the
validation ran.

`summary.csv` reports module recovery by arm and estimator. Both estimators
identified three modules in every non-null dataset. Their median truth
ARI was 1 for both strong arms and 0.869 for the moderate arm. In the
independent arm, bam detected no module in 10 of 10 datasets. INLA detected no
module in 9 of 10 datasets and assigned 25 genes to one module in seed
202740010. Both estimators converged for all 90 features in that dataset; their
correlation matrices differed by at most 0.00914 and their TOM matrices by at
most 0.00127. The WGCNA settings were retained unchanged.

`convergence.csv` reports feature-fit convergence counts. Full covariance,
correlation, TOM, labels, diagnostics, timing, and provenance remain in the
local artifact directory.
