# MAGIC marker and slide-effect exploration

This analysis compares three representative genes on the same 97,830 MAGIC
observations and 93 sections. We use the comparison to examine expression
patterns and the separation of a three-dimensional spatial field from a
section-level effect. The model and its covariance test remain under study.

## Marker choice and data

The source article's Figure 5d displays Foxp1, Nr4a3, Satb2 and Tfap2b as
regional markers. Foxp1 and Tfap2b are present in all 93 locally available
source H5AD feature tables and were extracted by exact point identity.
Snap25 was the original transfer example. The new marker choice is guided
by the published figure. See the [source article](https://doi.org/10.1038/s41588-024-01906-4)
and the [dataset construction note](magic-data.md).

The portable research object is `artifacts/datasets/MAGIC/MAGIC.rds`. It contains
a 97,830 by 3 raw-count matrix, 39 metadata columns, a 93-row section table,
the supplied tetrahedral mesh, and source hashes. QC variables are available
for subsequent analyses; the current fitted models use an intercept, the
three spatial coordinates, and the existing total-UMI exposure. Section identity
defines the added random effect. Raw H5AD count totals are preserved separately
from the transfer exposure because they differ slightly.

The source article reports 98,192 spots. The 93 available raw H5AD files and
the transfer both contain 97,830 spots. The reason for this difference remains
unresolved.

## Fitting specification

For each gene, we fit a negative-binomial model with linear predictor
`log(mu) = log(exposure) + beta + A u + Z b`. The comparison includes the
spatial field alone, an additional iid section effect, and an additional
Ornstein-Uhlenbeck (OU) section effect. The mesh has 1,962 nodes and 7,676
tetrahedra. The spatial kappa, count exposure, and historical integral
constraint are held fixed across these comparisons. These exploratory fits
use standalone INLA; the current package API uses its observation-mean
constraint and supports one target field.

All estimated hyperparameters have explicit flat priors on INLA's internal
scales: log negative-binomial size, log spatial precision, log section
precision, and log OU decay. The starting section precision is 100 and the
starting OU decay is 20 per millimetre. These are optimization starting values.
The integration settings are Gaussian/empirical Bayes, with two INLA threads
and no variational correction. Marginal likelihood is not used for comparison.

The OU process uses the actual z coordinates in millimetres and correlation
`exp(-phi * abs(z_s - z_t))`. Adjacent observed sections are 40 to 240 micrometres
apart. A rank AR1 over 1 to 93 would instead measure decay per retained-section
transition. INLA supports both [AR1](https://www.inla.r-inla-download.org/r-inla.org/doc/latent/ar1.pdf)
and [OU](https://www.inla.r-inla-download.org/r-inla.org/doc/latent/ou.pdf).
An equally spaced 10-micrometre AR1 lattice retaining unobserved intermediate
states is another possible distance-based construction; this pilot used OU.

## Findings from the flat-prior comparison

All nine fits returned mode status zero and no recorded INLA warnings. The
Snap25 spatial-only predictor reproduces the saved baseline within 5.7e-8.
The correlations below compare observed and fitted counts per 10,000 UMI;
they describe in-sample agreement.

| Gene | Spatial-only rate correlation | With iid sections | With OU sections | Correlation between baseline and OU spatial components | OU correlation at 60 micrometres |
|---|---:|---:|---:|---:|---:|
| Snap25 | 0.359 | 0.396 | 0.396 | 0.909 | 0.477 |
| Foxp1 | 0.614 | 0.621 | 0.621 | 0.651 | 0.935 |
| Tfap2b | 0.789 | 0.791 | 0.791 | 0.485 | 0.994 |

Tfap2b shows a localized pattern in the observed sections. For Foxp1 and
Tfap2b, adding a section process changes the spatial component substantially
while changing total fitted expression relatively little. For Tfap2b, the
baseline and OU total fitted rates correlate at 0.990, whereas their spatial
components correlate at 0.485. Thus similar overall fits can assign the
along-z signal differently to the spatial and section components.

Section adjustment also reduces the root mean square of the section-mean
Pearson residuals. For Snap25, the value changes from 0.178 to 0.014 under
either section model. This supports examining section-level structure. It
does not identify that structure as a technical effect: section-specific
biology and a genuine anatomical gradient can contribute to the same pattern.
The very high fitted OU correlation for Tfap2b makes this separation an
especially useful target for simulation and model diagnostics.

## Extending the spatial covariance test

The current sparse test accepts exactly one global SPDE block. Supporting
an additional section process requires the joint latent design `[A, Z]`
and block prior precision `diag(tau_s Q_s, Q_b)`. Its working posterior system
includes the cross block `A' W Z`. The section field must be integrated when
constructing the observation precision and fixed-effect adjustment. The
tested factor remains the constrained spatial **prior** factor; the posterior
Schur complement is a solver operation, not a replacement target kernel.

A further question concerns section effects shared across genes. Under a
null of zero cross-gene spatial covariance, correlated section effects can
still produce a nonzero expected pair score. For score coordinates
`a_j = F_j' P_j e_j`, the section contribution to its expectation is
`tr(F_j' P_j Z Cov(b_j,b_k) Z' P_k F_k)`. Integrating a separate marginal
section covariance for each gene does not generally remove this term.
The pairwise null therefore needs an explicit shared-section covariance
model, or a declared conditional target that projects out section effects.

The next calibration should distinguish independent and cross-gene shared
section effects, and genuine anatomical z variation. Known-parameter Gaussian
checks can first verify the algebra, followed by negative-binomial simulations
with estimated flat hyperparameters on the retained geometry. Agreement with
mgcv, type-I error, and null p-value uniformity remain the validation criteria.
The present real-data fits do not produce validated p-values for this extension.

## Reproduction and displays

Run the data construction steps in `magic-data.md`, then:

```r
source("inst/examples/magic-slide-exploration.R")
source("inst/examples/magic-slide-diagnostics.R")
source("inst/examples/magic-marker-visualization.R")
source("inst/examples/magic-marker-interactive.R")
```

The fit directory stores configurations, per-gene fits, section diagnostics,
and validation summaries. `magic-slide-comparison.csv` records the nine
exploratory fits summarized here. Marker plots show linear expression rates
with zero in white and increasing blue. The observed and fitted interactive
views share a within-gene scale, capped at the 98th percentile for display.
Actual section plots preserve all observations in the same six selected planes.
These remain exploration outputs; no representative gene has been finalized
for the README.
