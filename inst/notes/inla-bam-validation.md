# Local paired validation of bam and flat-prior INLA

## Purpose and source

This validation compares the frequentist mgcv/bam workflow with supplementary
INLA estimation using the author's flat log-hyperparameter priors. The package
source at launch is commit `f2ba738`, version `0.0.1.9005`. A separate local
library under `artifacts/inla-bam-validation/library` isolates the run from
the installed package used by earlier research jobs. Commit `3e5a60d`, version
`0.0.1.9006`, corrects nuisance-penalty rank validation. The corrected nuisance
reruns and global/local simulations use a second frozen installation under
`artifacts/inla-bam-validation/validated-library`; the other lanes retain the
original installation and do not use the affected nuisance block.

The primary estimator is `mgcv::bam(method = "fREML", discrete = TRUE)`.
INLA uses `list(prior = "flat", param = numeric(), initial = 0)` for spatial
log precision and NB log size; Gaussian simulations also explicitly use a
flat Gaussian observation log precision objective. Each comparison holds the response,
mesh, fixed kappa, observation-mean constraints, covariates and offsets fixed
between estimators. INLA retains sparse fitting and its sparse score backend
where supported; models with extra random blocks use the established backend.

## Replication and coverage

Each simulation case has ten distinct, fixed seeds. Both estimators receive
the same generated data for each seed. The per-lane CSV manifests identify
every attempted dataset; results retain failed or non-converged states without
substituting seeds. Reusing a completed checkpoint does not create a new
replicate or a new runtime measurement.

- Inference: Gaussian and NB responses, low counts, independent and correlated
  spatial fields, positive and negative association, nuisance terms, true
  zero-spatial-effect settings. Estimation, marginal tests and pair tests are
  compared using supported Liu and Davies calibration paths. A separate set
  of ten global/local datasets validates estimation and score-network extraction
  for the global, local and combined components.
- Modules: the four original WGCNA arms each use their first ten manifest
  seeds. The original design has 90 genes, 2,125 observations, 298 score
  coordinates, three planted groups of 30 genes, field amplitude 0.6,
  mean count 12 and NB size 15. The two strong arms have correlation 0.6,
  the moderate arm 0.3 and the independent arm zero.
- Scaling: ten fresh paired datasets at each of 2,000, 8,000 and 32,000
  observations, with the existing 15-by-15 mesh, two NB features, fixed
  kappa 6, mean count 0.3, size 2 and field amplitude 0.6.
- Real data: complete Visium-B positive components, 635 genes without cell-type
  adjustment and 491 genes with adjustment, including their 443 shared genes.
  Full component gene sets and pair universes are retained. These are observed
  analyses, not ten independent biological replicates.

## Comparison and interpretation

For each case, compare fitted spatial fields, smoothing and dispersion
estimates, marginal and pair p-values, signed scores, information, discovery
decisions, score correlation, TOM and module assignments as applicable.
Report all ten paired results and their distribution, alongside agreement
with known simulation truth. Separate feature fitting, marginal testing,
score construction, pair testing and WGCNA time; record thread counts and
concurrent workload when interpreting runtimes.

Agreement is assessed along these scientific outputs rather than by requiring
bitwise equality of different estimators. Backend equivalence on the same fit
is a separate numerical check. Differences are first checked for input,
geometry, prior, convergence or implementation mismatches; remaining differences
are reported as estimator differences. mgcv/bam remains the primary analysis,
and INLA provides the supplementary comparison.

Ten simulation replicates provide a paired validation sample for each case.
Calibration and discovery-rate summaries are reported with their denominators;
these results alone do not constitute a precise type-I-error or formal
equivalence study.

## Files and execution

The runnable scripts are in `inst/benchmarks/inla-bam-validation-*.R`.
Raw fits, input snapshots, logs and per-replicate checkpoints remain under
`artifacts/inla-bam-validation/`. Compact tables and simulation manifests are
versioned under `inst/validation/inla-bam/`. The existing research inputs,
outputs, seed manifests and running jobs are preserved.

Run the scripts from the package repository root with `Rscript`. Install the
chosen source in a separate R library first, and set
`MGCVST_VALIDATION_LIBRARY` to its absolute path. The scripts retain completed
checkpoints, so choose a fresh output location when a new fit is intended.

| Script suffix | Invocation / scope |
| --- | --- |
| `inference.R` | `--task=1` through `--task=100`, one fixed manifest row per process |
| `multigroup.R` | `--replicate=1` through `--replicate=10` |
| `scaling.R` | No arguments runs all 30 manifest rows; `--task=1` selects one |
| `wgcna.R` | No arguments runs all four arms, ten seeds each |
| `components.R` | No arguments runs both full observed components |
| `summarize.R` | Collects the completed inference and scaling checkpoints |
| `degeneracy.R` | Audits the recorded low-information BAM pairs without refitting |
| `feature-audit.R` | Checks real-component identities, convergence, parameters and marginal results |
| `network-audit.R` | Builds complete real-component networks and the non-convergence sensitivity comparison |
| `component-summary.R` | Recomputes full real-data pair, marginal, module and parameter summaries |
| `mode2-repeat.R` | Repeats each recorded nonzero-mode-status feature once with the original controls |

The table abbreviates the common filename prefix
`inst/benchmarks/inla-bam-validation-`. The historical WGCNA inputs and Visium-B
research archive remain required for their respective lanes. Their paths can
be supplied through `MGCVST_HISTORICAL_ROOT` and `MGCVST_REAL_DATA_ROOT`.
The inference output root is configurable through `MGCVST_INFERENCE_OUTPUT`;
the ten corrected random-nuisance tasks were directed to
`artifacts/inla-bam-validation/nuisance-correction`. The summary script uses
those corrected tasks and the original other 90 tasks.

## Completed simulation results

The ten inference cases each contain ten independent paired datasets and three
features. All ten corrected random-nuisance fits converge. These reruns use the
same responses, covariates and basis as their original attempts, six of which
failed because a numerically factorable spline penalty was assigned a rank
inconsistent with the adapter. The correction uses the same numerical rank
definition for nuisance-penalty validation and construction, while retaining
checks against zero and indefinite penalties.

Across the inference cases, the pair-test decisions at the raw 0.05 threshold
agree for every comparison with valid results from both estimators. Marginal
decisions also agree closely. The table retains numerical failures in the
attempted denominators.

| Test | Calibration | Attempted | Both valid | Invalid BAM | Invalid INLA | Raw decision agreement |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| Pair | Liu | 300 | 283 | 17 | 0 | 283/283 |
| Pair | Davies | 300 | 283 | 17 | 0 | 283/283 |
| Marginal | Liu | 300 | 300 | 0 | 0 | 298/300 |
| Marginal | Davies | 300 | 239 | 0 | 61 | 236/239 |

The 300 feature or pair comparisons arise from 100 independent datasets, with
three dependent comparisons within each dataset. Pair-score signs agree in
278/283 valid comparisons. BH-adjusted pair decisions agree in 282/283 Liu
comparisons and 283/283 Davies comparisons. Only feature1-feature2 is correlated
in the planted-pair cases; the other pairs remain null. Case- and pair-specific
outputs preserve this distinction when assessing detection and calibration.
The 61 invalid INLA Davies marginal results comprise 59 zero and two small
negative tail probabilities returned by the numerical routine, all with
`ifault = 0`. No fallback or replacement p-value was applied.

The 17 invalid BAM pairs all meet the shared calibration guard
`information <= 1e-10`; their positive information values range from
`1.08e-11` to `6.31e-11`. This guard assigns an `effective_rank = 0` sentinel
before evaluating the mixture spectrum. It does not establish algebraic rank
zero. Fifteen such pairs arise in the zero-spatial-effect cases and two in
low-count cases. Both dense BAM and sparse INLA use the same guard; the
corresponding INLA information values exceed the threshold. These results
remain uncalibrated in the summaries, with no conversion to non-rejections.
The current API applies BH jointly to the valid p-values within each requested
pair-test call. Invalid rows remain in the attempted count but are excluded
from the adjustment denominator. The validation preserves that behavior.

The scaling experiment contains ten independent datasets at each size and two
features per dataset. The table pools the twenty paired feature comparisons
for field agreement and summarizes dataset runtimes by their medians.

| Observations | Minimum field correlation between estimators | Median field RMSE, BAM / INLA | Median total seconds, BAM / INLA | Median paired BAM/INLA runtime ratio |
| ---: | ---: | ---: | ---: | ---: |
| 2,000 | 0.999973 | 0.249508 / 0.249659 | 1.765 / 1.910 | 0.92 |
| 8,000 | 0.999995 | 0.181380 / 0.181377 | 6.380 / 3.060 | 2.07 |
| 32,000 | 0.999777 | 0.118410 / 0.118248 | 16.360 / 7.410 | 2.06 |

Both estimators recover closely matching fields, with similar error against
the simulated latent field. Total time includes model preparation, estimation,
marginal testing, compaction and one pair test; shared data generation and basis
construction are excluded. At 32,000 observations, median times for fitting
plus the required marginal testing and compaction are similar
(BAM 6.85 s; INLA 6.94 s), and the sparse INLA pair stage accounts for
much of the total-time advantage (0.04 s versus 9.25 s). These measurements
describe the two-feature experiment on this computer under concurrent local
workloads.

All ten global/local fits converge for both estimators. The median, across
datasets, of the maximum absolute between-estimator score-correlation difference
is 0.0375 for global, 0.0352 for local and 0.0338 for combined extraction. This
three-feature case checks the component interfaces; the 90-feature simulation
provides the substantive module-recovery assessment. The current global/local
fit does not expose a pair-test engine or separate component marginal tests.
These interfaces are recorded as unavailable and are not counted as validated
pair or marginal comparisons.

The full package check for version `0.0.1.9006` completed with zero errors,
warnings and notes. Tests recorded 766 passing expectations and four skipped
pinned-source equivalence checks requiring an external `MGCVST_BASELINE`.
The passing tests include four nuisance-rank regression checks.

## WGCNA simulation results

All forty paired datasets are complete, with 3,600/3,600 feature fits
converged for each estimator. The versioned tables retain every seed,
per-dataset metric, module label and stage timing.

Each strong-correlation arm (`rho = 0.6`) recovers three modules in every
one of its ten datasets. BAM and INLA partitions agree exactly across all
twenty strong-signal datasets. The median ARI against the three planted
groups is one; the minimum is 0.9665 for both estimators.

At moderate correlation (`rho = 0.3`), both estimators recover three modules
in all ten datasets. Their ARI against truth has median 0.8691 and range
0.6974--0.9665. The minimum between-estimator ARI is 0.9310. Thus the
estimators give similar module structure, while finite-sample recovery and
some module boundaries differ under weaker signal. These are separate
comparisons: agreement between estimators and recovery of simulated truth.

With independent fields (`rho = 0`), BAM assigns no genes to modules in all
ten datasets. INLA gives the same result in nine datasets; in replicate ten
(seed `202740010`) it identifies one 25-gene module and leaves 65 genes
unassigned. This difference is retained as part of the independent-field
assessment. The ten repetitions characterize this validation sample rather
than establish a precise module false-positive rate.
Both fits use the same response matrix, feature order and settings in that
replicate. Maximum between-estimator differences are 0.00914 in correlation
and 0.00127 in TOM, consistent with a clustering boundary sensitive to small
estimator differences. No network parameters or seeds were changed.

## Real-component feature audit

The full 635-gene unadjusted and 491-gene cell-type-adjusted analyses retain
identical feature order and observations between estimators. Their feature
sets overlap at 443 genes. All BAM fits converge; INLA records 634/635 and
490/491 converged fits, respectively. The unadjusted feature
`ENSDARG00000058105` and adjusted feature `ENSDARG00000098458` return INLA
`mode_status = 2`. Both retain finite working states and are therefore available
under the current compact-fit API, whose availability check is based on
finite parameters rather than the convergence flag. Downstream summaries
distinguish finite comparisons from comparisons with both fits converged.

All retained marginal Liu p-values are finite, and both estimators reject
the spatial marginal null for every feature in these selected components.
The corresponding paired denominators with both fits converged are 634 and
490. This describes spatial signal in the selected components; these feature
sets were selected from earlier analyses and do not estimate a genome-wide
false-positive rate.

Median smoothing parameters are 0.006010 (BAM) and 0.006278 (INLA) without
cell-type adjustment, and 0.014106 and 0.014455 with adjustment. Under the
specified flat log-size prior, INLA also returns very large NB size values:
26 unadjusted and 15 adjusted features exceed `1e12`. These are consistent
with numerical estimates toward the Poisson limit. Finite size values and
optimizer status are recorded separately when assessing stability. The NB
working dispersion is one and is distinct from the estimated NB size.

One serial repeat with the same inputs, priors and controls returns
`mode_status = 0` for each previously non-converged feature. The smoothing
parameter changes by 0.81% without adjustment and 0.23% with adjustment.
Maximum absolute changes in the working residual/variance are
`0.000700/0.000904` and `0.004337/0.005181`, respectively. NB size changes more
substantially, from 1,552 to 6,144 and from 158 to 581. Spatial working states
are comparatively stable, while the flat size objective gives less stable
large-size estimates. The repeats remain separate from the original fits;
the difference in optimizer status is not attributed to a specific cause.

Full-component WGCNA identifies six modules without adjustment and three
with adjustment for both estimators. Between-estimator module ARI is 0.9588
and 0.9926, respectively. Unadjusted module sizes are
173, 151, 147, 78, 48 and 38 for BAM, compared with
168, 154, 144, 75, 55 and 39 for INLA. Adjusted sizes are
221, 208 and 62 for BAM, compared with 220, 209 and 62 for INLA.

Removing the single non-converged INLA feature from both estimators gives
between-estimator ARI 0.9573 without adjustment and 0.9926 with adjustment.
The score correlations among retained genes are unchanged by this removal.
The adjusted partitions and the unadjusted INLA partition are unchanged;
the unadjusted BAM partition changes slightly after reclustering
(ARI 0.9783 against its full-component partition).

Full-network construction takes 100.52 s for BAM and 2.01 s for INLA in the
unadjusted component, and 85.90 s and 2.29 s in the adjusted component. Most
of this difference is in score construction; adjacency, TOM and clustering
take about 0.25--0.33 s. These timings start from the compact fitted models
and exclude feature estimation. Mean absolute between-estimator correlation
differences are 0.00486 without adjustment and 0.00546 with adjustment.

For the adjusted component, all 120,295 requested pairs have valid p-values
from both estimators. The median absolute p-value difference is 0.000866;
the maximum is 0.37987, so individual comparisons can differ substantially.
Raw 0.05 decisions agree in 99.12% of pairs and BH-adjusted decisions in
99.10%. BAM discovers 58,091 pairs and INLA 57,955, including 57,481 shared
discoveries. Whole-call pair-test times are 1,303.06 s for BAM and 1,195.06 s
for INLA. Exhaustive calibrated pair tests show a smaller timing difference
than score-only network construction in this component.
Signed scores have Pearson correlation 0.99976, and the discovery-set Jaccard
index is 0.98149. Restricting the comparison to the 119,805 pairs with both
fits converged gives score correlation 0.99976, decision agreement 99.11%
and discovery Jaccard 0.98161.
This subset retains the BH decisions from the original 120,295-valid-pair
call; it is not a new multiple-testing family.

For the unadjusted component, all 201,295 requested pairs have valid p-values
from both estimators. Signed-score Pearson correlation is 0.99976, raw
decisions agree in 99.17% of pairs, and BH-adjusted decisions agree in 99.09%.
The median absolute p-value difference is 0.000664 and the maximum is 0.44240.
BAM discovers 103,733 pairs and INLA 102,555, with 102,230 shared discoveries
among a union of 104,058. Score directions agree in 99.68% of pairs.

The unadjusted whole-call pair times are 5,147.89 s for BAM and 2,004.05 s for
INLA. This run used 41 BAM chunks of at most 5,000 pairs and four INLA chunks
of at most 50,324 pairs, with four workers for each estimator. Chunk-local
state caching and output construction contribute to these observed times.
They describe the executed workflows rather than isolate an estimator-only
speed ratio. The adjusted comparison used four chunks for each estimator.
All pair results were combined before the API applied BH to the full set of
valid p-values in each model. The preserved non-converged feature fits are
finite and remain included in these original full-component results.

Within these component-wise pair universes, all 9,744 recorded positive
edges in the original unadjusted component and all 6,543 in the adjusted
component retain positive discoveries with both current estimators. The
archived edge-preservation tables document the comparison to the previous
analysis alongside the complete current within-component test results.
The saved edge lists retain their original screening and multiplicity
adjustments, including BY for the unadjusted network; current discoveries
use BH within each requested component. Edge retention is interpreted
within these stated testing families.

All 180 prespecified paired simulation datasets and both complete observed
components have finished. The two observed analyses contain 321,590 paired
comparisons in total, and all have valid Liu pair p-values for both estimators.
The package source, compact result tables, manifests and runnable scripts are
versioned together; raw inputs, fit objects and exhaustive pair tables remain
in the documented local artifact directory.
