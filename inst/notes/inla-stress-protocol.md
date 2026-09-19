> Historical protocol. Its Davies calculations and public BiocParallel
> downstream experiments describe the recorded study only. Current INLA
> downstream APIs use the single-global sparse Liu/OpenMP implementation.

# INLA stress and score-calibration protocol

Source at launch: `7e80306`, public package 0.0.1.9006. The preceding transfer
study established fitting and geometry checks only. This protocol adds
statistical null experiments and both requested parallel workloads.

## Frozen null design

Each of eight cases has 500 independent replicate datasets. The 2D public
API cases cross pair/marginal nulls with NB mean counts 0.3/3. Pair nulls have
two independent nonzero spatial fields; marginal nulls have no spatial field.
The design uses 200 observations and 36 mesh vertices, with kappa 5,
observation-mean constraint, average latent variance 0.45^2 for pair nulls,
NB size 2 and fixed heterogeneous offsets. Both estimated hyperpriors are
flat. Public INLA estimation, pair calibration and marginal testing are used.

Four additional standalone 3D cases cross independent/two-gene fitting and
40-gene shared-hyperparameter fitting with means 0.3/3. They use 200
observations, a genuine 64-node tetrahedral volume, fixed kappa 5, a
volume-integral constraint and average spatial variance 0.45^2. Each of the
40 fields is independently generated. The tested pair is fixed in advance
as genes 1 and 2; the other 38 inform the joint hyperparameters. Genes 1/2
have identical response vectors in independent and joint experiments because
each gene uses its own deterministic seed. The 780 pairs in one dataset are
not treated as 780 independent simulations.

The standalone 3D score uses the same bilinear Gaussian-mixture calibration
function as the package. It combines the expected NB working state with the
native conditional Gaussian intercept variance. It is a research extension,
not validation of an existing public 3D score API. This Gaussian EB setup
uses no hyperparameter integration or posterior sampling. The shared
hyperparameter estimates may induce dependence, which this Monte Carlo study
measures rather than assumes away.

Davies and Liu are both reported. Nominal thresholds are 0.05 and 0.01,
with exact binomial intervals on valid replicates and all-attempt bounds
[rejections/500, (rejections+unavailable)/500]. Fits or scores that return no
p-value are retained, including matrices that trigger the implementation's
absolute -1e-10 eigenvalue cutoff. Triggering that cutoff is a software
outcome; interpreting it requires the actual spectrum and its scale.
Unavailable p-values are not replaced with p=1. The research worker applies
no additional spectral repair. The existing package Davies path does zero
negative eigenvalues within its tolerance, whereas Liu uses trace moments.
Convergence status, extreme NB sizes, and covariance/constraint issues remain
separate diagnostics. Five hundred repetitions give a Monte Carlo standard
error about 0.00975 at a true 5% rate; compatibility with 5% does not prove
exact size control or validate extreme tails.

The first ten 2D low-count pair replicates were an execution pilot and belong
to the fixed 500-replicate set. An earlier standalone 3D pilot had a different
random-stream layout and is excluded. No rejection-rate-dependent stopping
or expansion is used. A process failure records its first missing replicate;
subsequent unattempted replicates continue without retrying that failure.

Before any formal 3D job launched, diagnostic output was expanded at the
author's request to retain each original M, its factor F, G, w, working
variances, native Vp and expected-information Vp. This output-only addition
preserves the frozen simulations, priors, fit controls, statistic and
calibration decisions. It permits full-spectrum inspection without fitting
again. A separate ten-replicate spectrum diagnostic compares curvature
conventions; those reruns do not replace formal Monte Carlo results.

## Stress workloads

MAGIC geometry: 97,830 observations and 1,962 nodes, native Q and integral
constraint. Forty independent constrained fields are simulated once with
sigma 0.4, NB size 15, intercept 0.75 and the real UMI exposures. The same
response arrays are retained for every configuration.

Feature-parallel configurations are workers/INLA threads 1/1, 2/1, 4/1,
8/1 and 4/2, each fitting the same ten full-size genes. Compare parameter
estimates and field modes against 1/1, task completion, process-tree peak
memory, batch wall time and throughput. These are ten fit tasks per
configuration, not ten independent batches.

The joint workload has 40 distinct fields, 40 flat intercepts and two shared
flat hyperparameters (spatial log precision and NB log size). It uses sparse
block A/Q matrices and one integral constraint per field. Ten repeated fits
are requested at 5,000 observations per gene and ten at the full 97,830,
each using four INLA threads. The first unmonitored 5,000-point smoke run
is retained separately and is not counted as a monitored stress repetition.
Shared and separate estimators need not give identical estimates; numerical
invariance is evaluated within a fixed model configuration.

The supervisor samples owned process trees every 0.5 seconds, records summed
RSS/private bytes and available system RAM, and terminates only its own jobs
if one tree exceeds 24 GiB private memory or system available RAM falls below
6 GiB. Unattempted stress repetitions after a resource stop are explicitly
labelled. No out-of-memory event is intentionally induced. Other user research
processes are preserved. Timing under their background load is recorded as
loaded-workstation performance, not an isolated hardware benchmark. RSS sums
can count shared pages repeatedly; sampled peaks can miss shorter spikes.

All raw inputs, checkpoints, failures and logs are under
`artifacts/inla-stress-calibration/`. Completed results will be reported
separately from this protocol. Public API and priors are unchanged by these
research scripts.
