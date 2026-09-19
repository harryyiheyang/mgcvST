> Historical diagnostic record. The standalone Davies comparisons below are
> retained as evidence for the study and are not calls supported by the current
> INLA downstream API, which uses sparse Liu calibration.

# INLA score spectrum diagnostics

Inspecting the saved matrices identifies two sources of small negative
eigenvalues: a mismatch between expected NB working curvature and the native
INLA conditional intercept variance, and smaller numerical differences in
that native variance. Their statistical consequences are evaluated separately
from the existing implementation's absolute eigenvalue cutoff. Both fitted
hyperpriors remain flat on the internal log-precision and log-size scales.

## Matrix identity

The standalone 3D calculation uses

\[
V=D+FF^T,\quad v=V^{-1}\mathbf1,\quad
P=V^{-1}-c vv^T,\quad M=F^TPF,
\]

where \(D_i=1/\mu_i+1/\mathrm{size}\) is the expected NB working
variance and \(c\) is the native conditional INLA intercept variance.
Writing \(G=F^TV^{-1}F\), \(w=F^Tv\), and
\(c_0=(\mathbf1^TV^{-1}\mathbf1)^{-1}\) gives the exact identity

\[
M-M_0=-(c-c_0)ww^T,\qquad M_0=G-c_0ww^T.
\]

Thus a native variance larger than \(c_0\) subtracts an additional
rank-one term from the expected-curvature matrix. Whether it creates a
negative eigenvalue of \(M\) also depends on the spatial factor. In the
fixed pilot geometry, the 200-by-63 observation factor has numerical rank
62, and the constant observation vector lies in its column space to a
relative residual of approximately \(1.1\times10^{-14}\). A volume-integral
constraint on the field does not itself force the sampled spatial factor
to be orthogonal to the observation intercept.

For the NB log link, the observed negative Hessian weight is

\[
W_{o,i}=\frac{\mu_i(1+y_i/\mathrm{size})}
 {(1+\mu_i/\mathrm{size})^2},
\]

whereas the expected weight is \(W_{e,i}=1/D_i\). The diagnostic compares
both resulting Schur-complement intercept variances with the native one.
The supplied Matérn precision has positive fixed kappa and is proper;
conditioning it on the zero volume integral is an additional constraint.

## Ten saved pilot replicates

The following entries come from the complete, unmodified matrix spectra.
These diagnostic reruns remain separate from the formal simulation outputs.

| Replicate / feature | Minimum eigenvalue | Maximum eigenvalue | Minimum / maximum | Native / expected variance | Native / observed variance |
|---|---:|---:|---:|---:|---:|
| 3 / 1 | -5.5691e-6 | 0.64375 | -8.6509e-6 | 1.005502 | 0.9999895 |
| 5 / 1 | -1.4875e-5 | 0.75832 | -1.9615e-5 | 1.010226 | 0.9999850 |
| 3 / 2 | -1.4303e-10 | 0.22993 | -6.2206e-10 | 1.000000729 | 1.000000729 |

For the first two entries, relative matrix asymmetry and eigendecomposition
residuals are about 2e-16 to 6e-16. The negative values are therefore much
larger than the numerical error of that eigendecomposition. Their native
variances closely track observed curvature, supporting curvature mismatch
as the principal cause. The third entry has almost identical expected and
observed curvature and a much smaller native-variance discrepancy. Its
spectral effect is correspondingly smaller. The rank-one identity is
verified to relative error at most 1.2e-15 in the pilot.

The larger two negative components account for only approximately 7.2e-7
and 1.5e-6 of total absolute spectral mass. This small mass motivates a
sensitivity calculation; it does not determine calibration by itself.

Three calculations retain the original results as the reference:

1. Keep the original statistic and use its fixed-working-Gaussian covariance
   \(F^TPVPF\). This sandwich calculation does not remove intercept mean
   leakage or account for estimated nuisance parameters.
2. Use expected curvature consistently and recompute the feature scores,
   pair statistic, and covariance.
3. Use observed curvature consistently, including the corresponding Newton
   working response, and recompute the scores and covariance.

Across the eight pilot replicates with an available native p-value, the
largest absolute Davies p-value changes are 4.60e-6, 0.00309, and 0.00220,
respectively. A separate spectral-projection diagnostic measures the mass
that would be removed; projected matrices are not used to generate p-values.
These ten replicates establish sensitivity magnitude, not type-I control.

## Completed low-count 3D null case

All 500 independent-field low-count datasets have been attempted. The
original implementation returns 313 p-values per calibration method.
Its fixed -1e-10 guard blocks 184 datasets; three more produce an unavailable
degenerate calibration. Across all 1,000 feature states, the most negative
eigenvalue is approximately -5.074e-5, the most negative relative eigenvalue
is -5.756e-5, and the largest negative spectral-mass fraction is 3.285e-6.
These quantities describe the matrices independently of the cutoff.

Four individual fit states have extremely negative fitted NB log size,
nonzero mode status, infinite working variance, and zero original score
matrices. Their expected intercept variance is infinite, so the expression
for an expected-curvature matrix contains 0 times infinity. Those states
are retained with their explicit diagnostic reason.

The post-hoc sensitivity calculations use the same 500 saved datasets and
are additional diagnostics following discovery of the curvature mismatch.
The sandwich and expected calculations each return 496 p-values. At a 5%
threshold, both Davies and Liu reject 12/496, or 2.42%, with an exact 95%
binomial interval of 1.26%–4.19%. Allowing either outcome for all four
unavailable datasets gives a descriptive all-500 rejection range of
2.4%–3.2%. Observed-curvature Davies has the same count; observed-curvature
Liu rejects 11/496. All three routes reject one dataset at the 1% threshold.
The original available subset rejects 7/313 at 5%; its 187 unavailable
datasets prevent that subset from establishing overall calibration.

These results show that the original cutoff substantially reduces output
availability while the examined covariance alternatives give similar,
conservative rejection rates in this low-count case. Other count levels and
shared-hyperparameter fits are reported separately in the full stress and
calibration study. No package API, prior, or production calibration rule is
changed by this diagnostic.

At mean count 3, all 500 independent-field datasets have also completed.
The original hybrid calculation returns only nine p-values per method.
The most negative relative eigenvalue across the 1,000 states is about
-2.865e-4, and the largest negative spectral-mass fraction is about 1.10e-5.
All three additional routes return 500 p-values. At the 5% threshold,
sandwich Davies/Liu reject 21/19 datasets, expected-curvature Davies/Liu
reject 23/20, and observed-curvature Davies/Liu reject 23/22. The resulting
rates of 3.8%–4.6% have exact 95% intervals containing 5%. Every route
rejects two datasets at the 1% threshold. These post-hoc comparisons
support investigating a consistent score construction; they do not
retroactively validate the original hybrid rule.

## Reproduction

The scripts are `inst/benchmarks/inla-score-spectrum.R`,
`inla-score-spectrum-sensitivity.R`, and `inla-score-null-sensitivity.R`.
The formal results retain the unmodified matrices and reconstruction inputs
alongside the calibration outputs. Raw states are under
`artifacts/inla-stress-calibration/score-spectra/` and `null3d-paired/`.
Compact pilot results are under `inst/validation/inla-stress/spectrum-pilot/`.
The existing Davies implementation zeroes tolerated negative eigenvalues
inside its PSD factorization; Liu uses trace moments. The research workers
introduce no additional clipping to make the formal results pass.
