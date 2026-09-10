# INLA low-count investigation and posterior covariance contract

The implementation now obtains the nuisance `Vp` block directly from INLA's
conditional Gaussian posterior at its empirical-Bayes configuration. It uses
selected sparse solves of the stored posterior precision and applies all exact
constraints. This replaces the earlier expected-Fisher reconstruction. The
earlier matrix is retained separately as `expected_nuisance_covariance` for
diagnostics, and is not the default score `Vp`.

## Main mechanism

For the investigated constrained FEM geometry, the average observation-space
prior variance at the original precision multiplier `tau=1` is only
`c=0.0016131772`. The spatial field used in the null experiment has average
variance `0.36`, so its true FEM multiplier is `c/0.36=0.0044810478`.
The prior `log(tau) ~ N(0, 9)` is consequently centered on a field whose
median average spatial standard deviation is about `0.0402`, whereas the
experiment's standard deviation is `0.6`. With weak count information, the
prior can push the fitted spatial variance down substantially. The resulting
plug-in null covariance does not adequately represent the actual null fields.

This is an inference supported by matched component swaps, true-parameter
controls, and objective-function checks. Three fixed-size NB profiles agreed
with an explicit constrained Laplace calculation; their objective difference
had slope below `0.001` against log precision. A missing or duplicated
constraint dimension would produce a slope around `0.5` in magnitude. The
Gaussian check agreed with the exact restricted likelihood. No such
normalizer error was found. The same explicit objective reproduces the
strong shrinkage on the unscaled parameter, so it is not evidence of an
INLA optimizer malfunction.

## Independent validation

The first 100 pairs were used for exploration. A subsequent fixed 500-pair
validation used independent L'Ecuyer-CMRG seeds `161000+i`, 200 irregular
observations, 36 mesh nodes, fixed kappa 6, NB mean 0.3 and size 2, and two
independent mean-zero spatial fields with average variance 0.36.

At nominal 5%, using the current conditioned SPDE score kernel:

| Construction | Rejections / 500 | Rate | Exact 95% binomial CI |
| --- | ---: | ---: | ---: |
| Original scale, legacy expected-Fisher Vp | 72 | 14.4% | 11.4%–17.8% |
| Original scale, native Vp + expected working state | 72 | 14.4% | 11.4%–17.8% |
| Fixed true precision and NB size, native Vp + expected working state | 20 | 4.0% | 2.5%–6.1% |
| Observation scale, native Vp + expected working state | 4 | 0.8% | 0.2%–2.0% |

All entries were valid with no Davies fallback. Thus replacing Vp alone did
not resolve the inflation. Observation-scale priors removed inflation in
this experiment but were substantially conservative; this is not a claim
of exact 5% calibration or a universal robustness guarantee.
Native Vp with the observed working state rejected 73/500 on the original
scale and 4/500 on the observation scale.

The fixed-truth run reused all 500 validation responses. Its centered-raw
kernel result with native Vp was 24/500 (4.8%). This supports a major
hyperparameter-estimation contribution without proving that every part of
the plug-in count reference law is exact.

Two additional prespecified matched 200-pair checks used native Vp and the
observation-scale prior. At NB mean 0.1 and zero cross-correlation, the
conditioned score rejected 3/200 (1.5%, exact CI 0.31%–4.32%). At NB mean
0.3 and latent cross-correlation 0.7, its power was 19/200 (9.5%, CI
5.82%–14.44%). Matched bam/fREML results were 4/188 and 32/198 respectively.
The 12 and 2 bam invalid scores imply all-attempt rejection bounds of
2%–8% and 16%–17%. INLA had no invalid results or Davies fallbacks in these
checks. Normalized INLA rejected less often than bam under this alternative,
consistent with a possible power cost. Without a paired raw-versus-scaled
INLA alternative comparison, this difference cannot be attributed to scaling
alone and does not establish superiority of the observation-scale option.

## Explicit observation-scale option

```r
model <- inlaST.set(
  response ~ offset(offset0), data, basis, family = mgcv::nb(),
  precision_scale = "observation"
)
fit <- inlaST.estimate(Y, model, diagnostics = TRUE)
```

This sets `Q_internal=c*Q_original`, with the required
`log(tau_internal) ~ N(0, 9)` prior. At `tau_internal=1`, the constrained
field's average observation marginal variance is one. Returned smoothing
multipliers remain `tau_original=c*tau_internal`, so the common score
geometry and `lambda=dispersion*tau_original` contract retain their units.
`control$fixed_precision` is also expressed in original FEM units.

This is explicitly a different physical prior from placing `N(0,9)` on the
original FEM multiplier. It is optional: `precision_scale="raw"` remains
the default. No PC prior is introduced, and no mean-zero constraint is
removed. Nuisance smooths retain their existing precision scale.

## Covariance consistency and raw kernels

For symmetric positive-semidefinite P and the ordinary orthogonal centering
projector C, `sqrt(P) C G C sqrt(P) = sqrt(P) G sqrt(P)` remains correct
when `P 1=0`. Native posterior Vp combined with separately reconstructed
expected-Fisher working variance need not satisfy that condition exactly.
The exploratory study therefore reported the conditioned kernel, the fully
raw kernel, and the explicitly centered raw factor separately. Fully raw
kernel attempts often produced an indefinite plug-in score matrix under
that mixture of curvatures. They are not silently discarded or represented
as an equivalent no-centering method.
If P is indefinite, the real symmetric positive-semidefinite square root
used in this identity is itself unavailable.

The production score geometry remains the conditioned SPDE geometry. A
posterior-native score representation is investigated separately; it should
not be confused with merely substituting one nuisance covariance block.
The experimental `a=E(xi|y), M=I-Cov(xi|y)` construction had approximately
`1.3e-4` discrepancy from the Gaussian reference in four fixed-parameter
fits. Its raw rank-one coordinate identity held to `6.66e-16`; four cached
NB fits were closer to observed-working than expected-working curvature.
These local checks do not establish count-data type-I calibration.

## Reproducibility correction

The original serial bam experiment and the parallel INLA experiment had
different RNG kinds despite identical integer seeds. The corrected bam
comparison explicitly fixes L'Ecuyer-CMRG and verifies response arrays
against the saved INLA cache. On the original 500 pairs, bam rejected
15/496 (3.02%) with the conditioned kernel and 18/496 (3.63%) with the
centered raw factor. Four very large smoothing estimates gave score
information below the existing package threshold. Treating those outcomes
as unknown gives rejection bounds 3.0%–3.8% and 3.6%–4.4% over all 500
attempts. Thus the corrected matched comparison still does not reproduce
the original INLA inflation.

Research scripts are `inst/benchmarks/inla-lowcount-*.R` and
`inst/benchmarks/inla-posterior-*.R`; raw results and frozen scripts are
stored under `artifacts/lowcount-investigation/` in the development checkout.
