# Flat log-hyperparameter objectives and sparse scaling

Spatial log precision and NB log size use the following defaults:

```r
flat <- list(prior = "flat", param = numeric(), initial = 0)
model <- inlaST.set(
  complete_formula, data, family,
  control = list(precision_prior = flat, nb_size_prior = flat)
)
fit <- inlaST.estimate(Y, model)
```

In version 0.0.1.9003, omitting either `precision_prior` or
`nb_size_prior` selects this flat objective. Gaussian observation log
precision still defaults to N(0,9). Any prior can be overridden with a
validated INLA prior object; a prior object supplied to
`inlaST.estimate()` replaces the corresponding object stored by
`inlaST.set()`. Every spatial observation-mean constraint stays active.
`initial` for a flat prior is any finite value on the internal log scale;
it is an optimizer starting value, not a prior center.

INLA's `flat` density is constant on its internal parameter, so flat log
precision means density proportional to `1/tau` on the positive precision
scale. It removes the normal-prior penalty from that optimization objective.
See the [official flat-prior definition](https://www.inla.r-inla-download.org/r-inla.org/doc/prior/prior-flat.pdf).
For a fixed mesh and kappa, multiplying Q by a positive constant shifts log
precision. A flat log prior is unchanged by that shift; the ideal physical
optimizer is therefore unchanged when the remainder of the objective and
starting conditions are transformed consistently. Numerical optimizers can
still stop at different points on flat or multimodal profiles.

This option does not establish posterior propriety. As tau tends to infinity,
the spatial field disappears and the likelihood can approach a positive
limit. A constant prior over the unbounded log-tau tail then has infinite
posterior mass. The current engine uses a single empirical-Bayes Gaussian
configuration, and returns the native conditional covariance at its estimated
hyperparameters. It does not integrate an improper hyperparameter posterior.
An optimizer status of zero and finite theta do not certify an interior
maximum. Zero spatial variance and the NB Poisson limit need explicit
boundary checks; no arbitrary precision clipping has been introduced.

Historical simulation results in this repository must be read with their
recorded prior configuration. In particular, the earlier low-count
500-pair comparison used N(0,9) priors on spatial log precision and NB log
size, before the current flat defaults.

The package uses sparse raw SPDE precision during INLA fitting. From version
0.0.1.9001, `score_backend="auto"` also uses sparse precision solves for a
single global SPDE with fixed nuisance terms. This computes the same score
as `score_backend="dense"`, including native INLA `Vp`, without constructing
a dense observation-by-mesh score factor for each feature. Models with
additional random blocks retain the dense backend. Shared model setup and
marginal testing still use dense constrained geometry, and the score summary
is mesh-dimensional and dense. Consequently the relevant dimensions are
both the number of observations and the number of mesh nodes. The benchmark
reports fitting, preparation, marginal and pair testing separately. Neither
score backend allocates an n-by-n observation covariance matrix.

Here “the same score” means the score for the SPDE kernel already conditioned
on the observation-mean constraint. Reapplying that centering to the
conditioned kernel has no effect. The native-`Vp` nuisance operator combined
with the expected working state need not satisfy `P %*% 1 = 0`; this does
not turn the conditioned kernel into the raw kernel or into an independently
formed `C %*% G_raw %*% C`.

Research entry points are `inst/benchmarks/inla-flat-null.R`,
`inla-flat-boundary-audit.R`, and `inla-large-sparse-benchmark.R`. Raw results,
frozen scripts, protocols and the final comparison are stored locally under
`artifacts/flat-prior-investigation/`.
