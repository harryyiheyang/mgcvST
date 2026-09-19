# INLA estimator for mgcvST

`inlaST.set()` prepares a shared model and `inlaST.estimate()` fits each
feature with INLA. The output preserves the compact working-model inputs
used by `inlaST.test()` and `inlaST.marginal()`. The preferred setup is a
complete mgcv formula containing `s(..., bs="spde", xt=basis)`. The same
adapter also accepts a prepared `G` from `mgcvST.set()` and the historical
separate `basis` argument. These routes share mgcv formula parsing and score
geometry; feature estimation does not call `gam()` or `bam()`.
The current package accepts exactly one global SPDE target. `setting =
"global_local"`, a second SPDE term, and a named global/local basis list were
removed. A native sparse setup may instead supply `mesh`, fixed `kappa`, and
two- or three-dimensional `coordinates`; its formula contains only the
response, offset and parametric covariates.

## Required centering constraint

For the single spatial component, let `A` map the raw finite-element mesh
coefficients `u` to the **actual model observations**. The enforced constraint
is

```
g = crossprod(A, rep(1 / nrow(A), nrow(A)))
crossprod(g, u) = 0
```

Thus `mean(A %*% u) = 0`. This is not generally equivalent to summing mesh
coefficients to zero, or to integrating the field over the mesh. The
constraint cannot be disabled by an estimator control. An unconstrained input
basis is converted to the observation-centered basis during setup.

INLA receives the raw sparse precision matrix and a linear constraint.
The sparse score uses raw coordinates and constrained precision solves.
Legacy setup retains equivalent projected geometry for formula compatibility;
native mesh setup retains only the sparse projector and precision.

## Statistical interpretation

The implementation fixes the spatial range parameter `kappa` and
uses INLA empirical-Bayes hyperparameter estimation. INLA's hyperparameter
priors remain part of that estimate. Spatial log precision and NB log size
default to improper flat objectives. Gaussian observation log precision
defaults to `N(0, 3^2)`. Each can be replaced with a validated INLA prior,
including a normal prior with user-selected mean and precision, log-Gamma, or
an INLA expression/table definition.
Choosing an empirical-Bayes integration
strategy does not turn it into mgcv REML. Parameter and test results are
therefore not promised to equal those from an independently optimized GAM.
The returned estimation metadata records the actual controls and priors.
The generic0 marginal-likelihood diagnostic omits a constant depending on
the fixed precision geometry; do not use it to compare different meshes,
precision matrices or `kappa` values.

After hyperparameter estimation, the conditional field estimate is used to
construct the likelihood's working response and Fisher working variance.
Both latent coefficients and hyperparameters are read from INLA's joint
mode, rather than independently transformed marginal modes. Explicit
constraints also reduce the effective dimension of the precision
normalizer; the constrained SPDE has `m - 1` free field coordinates.
INLA's native constrained conditional posterior covariance block (`Vp`) is
used for nuisance adjustment. The expected-Fisher reconstruction is computed
only with `inlaST.estimate(..., diagnostics = TRUE)` and is not substituted
for that block. Otherwise each `expected_nuisance_covariance` entry is `NULL`,
avoiding an extra precision-matrix construction and solve per feature.
The frequentist score construction is then applied to those fixed inputs.
It does not substitute INLA credible intervals for frequentist calibration.

The score kernel is already conditioned on the observation-mean constraint.
Applying the same centering projection to that conditioned kernel again
leaves it unchanged. This statement concerns the kernel. The nuisance
operator combines native INLA `Vp` with an expected-likelihood working
state, and therefore need not obey `P %*% 1 = 0` numerically. Consequently,
the conditioned kernel must not be identified with the unconditioned raw
kernel or with a separately constructed `C %*% G_raw %*% C`.

For Gaussian observations and identical fixed hyperparameters, the
constrained conditional solution can be checked directly against the
projected penalized least-squares equations. For count observations, the
working-model score calibration remains approximate. Small computational
tests do not establish nominal tail probabilities after fitting
hyperparameters on the same data.

## Estimator controls

The `control` argument is a named list. Unknown control names are rejected.

| Control | Default | Meaning |
| --- | --- | --- |
| `int_strategy` | `"eb"` | The supported hyperparameter integration strategy. |
| `latent_strategy` | `"gaussian"` | The supported conditional latent approximation. |
| `num_threads` | `1L` | Threads per INLA feature fit. |
| `fixed_precision` | `NULL` | Optional fixed positive multiplier of each raw latent precision matrix. |
| `gaussian_precision` | `NULL` | Optional fixed inverse Gaussian residual variance. |
| `nb_size` | `NULL` | Optional fixed negative-binomial size, with variance `mu + mu^2 / size`. |
| `precision_prior` | `list(prior="flat", param=numeric(), initial=0)` | Improper flat objective on log latent precision. |
| `gaussian_precision_prior` | `list(prior="flat", param=numeric(), initial=0)` | Flat objective on log inverse Gaussian residual variance. |
| `nb_size_prior` | `list(prior="flat", param=numeric(), initial=0)` | Improper flat objective on log NB size. |
| `control.inla` | `list()` | Validated numerical INLA tuning such as `tolerance`; Gaussian/EB strategy remains mandatory. |
| `fixed_effect_precision` | `0` | Fixed effects are unpenalized by default. |
| `verbose` | `FALSE` | INLA engine output. |

Prior `initial` values are on INLA's internal logarithmic scale. Fixed values
such as `fixed_precision`, `gaussian_precision` and `nb_size` are on their
natural positive scales. The mgcv-compatible smoothing multiplier is
`lambda = dispersion * precision`, where `dispersion` is the Gaussian
residual variance and is one for Poisson and negative-binomial observations.
The Gaussian observation distribution remains Gaussian. Its residual variance
`sigma^2` is the Gaussian dispersion parameter; INLA parameterizes its inverse
as `precision = 1/sigma^2`. A flat prior on `log(precision) = -log(sigma^2)`
corresponds to density proportional to `1/sigma^2` on the variance scale. This
is distinct from NB size: the NB conditional variance is `mu + mu^2/size`,
and its working dispersion multiplier is one in this interface.

An explicitly requested normal hyperprior takes `(mean, precision)`, so
`1/9` encodes standard deviation 3. A flat prior may be abbreviated as
`list(prior="flat")`; it is
an improper density constant on the internal logarithmic scale. Registered
scalar INLA priors and explicit `expression:` or `table:` definitions are
validated and passed through with their actual parameters recorded in the
fit metadata. Prior objects are replaced whole rather than recursively
merged.

Controls supplied to `inlaST.set(..., control=...)` are stored with the
shared model. Controls supplied later to `inlaST.estimate(..., control=...)`
override matching values; nested `control.inla` tuning values are merged by
name. Native tuning cannot change the required Gaussian latent approximation
or single empirical-Bayes configuration. See the
[flat-prior investigation](inla-flat-prior.md) for boundary limitations.

The public estimator defaults to `BiocParallel::SerialParam()`. Use
`BPPARAM` for feature-level parallelism; keep `num_threads=1` when fitting
several features concurrently. All workers must load the same package
version as the parent process.

## Performance comparison

Use `inst/benchmarks/inla-estimator.R` to compare actual fits on the same
observations, mesh, centering constraint and fixed `kappa`. Report setup,
feature fitting and score calculation separately, including INLA startup
overhead. Use equal thread counts and avoid timing while other workers run.

`bam(method = "fREML", discrete = TRUE)` reduces design evaluation and
crossproduct costs, while this INLA backend retains the raw SPDE precision
sparsity. Either may win for a particular mesh size. Fitting speedup need not
equal total pipeline speedup.

`inlaST.estimate()` always uses the sparse score path, which requires exactly
one global target SPDE and only fixed nuisance terms. It preserves the
conditioned score kernel and native nuisance `Vp`; it does not change the
statistical test. The dense INLA score has been removed: a model with
additional random blocks is rejected by `inlaST.set()` with the capability
reason, and there is no `score_backend` argument any more. Numerical errors on
the sparse path are not silently replaced with a different method.

The sparse path forms `A' D^-1 A` and solves constrained systems with
`tau Q + A' D^-1 A`. It avoids building an observation-by-mesh score factor
per feature or any observation-by-observation covariance. The final score
summary remains a dense mesh-dimensional matrix. Shared model setup and
marginal testing retain dense projected geometry, so increasing mesh size
still costs substantially more than increasing spots at a fixed mesh.

## Implementation references

- [INLA generic0 precision model](https://www.inla.r-inla-download.org/r-inla.org/doc/latent/generic0.pdf):
  the positive precision is parameterized on the log scale.
- [INLA latent-field constraints and rank deficiency](https://www.r-inla.org/learnmore/docs/reference/f.html):
  an explicit `rankdef` must account for the extra constraints.
- [mgcv bam documentation](https://stat.ethz.ch/R-manual/R-devel/library/mgcv/html/bam.html):
  the comparison uses `method="fREML", discrete=TRUE` and a reused setup object.
