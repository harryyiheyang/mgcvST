# Interpretation of the re-estimated-hyperparameter null simulation

This note records the independent review of the simulation protocol and worker code, followed by the completed results. All seven scenarios have now finished with 500 independent replicates each (5,500 INLA feature fits). Earlier sections retain the detailed reasoning from the sequential review; the final results and the matched bam comparison are summarized below.

These historical INLA results used the original expected-Fisher reconstruction
of nuisance Vp. The current implementation instead extracts INLA's actual
conditional posterior Vp, as requested by the user. Subsequent investigation
of this correction and the prior scale is described in
`inla-lowcount-investigation.md`.

## Completed study

At nominal 5%, the Davies rejection rates were:

| Scenario | Conditioned SPDE kernel | Raw G with correct null P | Raw G 95% binomial CI |
| --- | ---: | ---: | ---: |
| Gaussian pair, kappa 6 | 6.2% | 6.4% | 4.4%–8.9% |
| NB mean 3 pair, kappa 6 | 6.0% | 6.8% | 4.8%–9.4% |
| NB mean 0.3 pair, kappa 6 | 15.2% | 15.4% | 12.3%–18.9% |
| NB mean 3 pair, kappa 0.7 | 5.6% | 5.6% | 3.8%–8.0% |
| Gaussian marginal, kappa 0.7 | 4.0% | 4.0% | 2.5%–6.1% |
| NB mean 3 marginal, kappa 0.7 | 2.8% | 2.8% | 1.5%–4.7% |
| NB mean 0.3 marginal, kappa 0.7 | 4.8% | 4.8% | 3.1%–7.1% |

All entries have 500 valid results and no Davies fallback. Every fitted
spatial field retained observation mean zero (maximum absolute error
`1.23e-16`), with log-hyperparameter `N(0,3^2)` priors and no PC priors.
The minimum working weight was `0.02949`, above the marginal masking
threshold. Confidence intervals containing 5% do not prove exact size
control; the NB mean-3 marginal result is conservative in this experiment.

The low-count pair result clearly fails size control. The first serial bam
experiment used Mersenne-Twister whereas the INLA Snow workers used
L'Ecuyer-CMRG. It was therefore a same-DGP comparison, not the identical
response pairs originally claimed. After explicitly fixing the RNG kind
and checking the first 100 response pairs against the saved INLA cache,
the corrected matched bam experiment rejected 15/496 (3.02%) with the
conditioned kernel and 18/496 (3.63%) with raw G. Four converged replicates
had score information below the package's degeneracy threshold. Their
unknown outcomes give all-attempt bounds 3.0%–3.8% and 3.6%–4.4%.
There were no Davies fallbacks; the production score oracle error was
at most `9.99e-15`. The corrected paired comparison still does not
reproduce the INLA inflation. Full corrected artifacts are under
`artifacts/constraint-type1/mgcv-low-count-matched/`.

## What the pair experiment tests

The pair null is zero cross-feature spatial covariance, not zero spatial variation within either feature. Each replicate generates two independent, nonzero, observation-mean-zero SPDE fields. Thus the null parameter is the association between the two fields, with both marginal spatial variances positive.

For fixed nuisance quantities and fixed hyperparameters, the pair statistic is

\[
U=a_1^\mathsf{T}a_2,
\qquad a_j\sim N(0,H_j),
\]

with independent `a1` and `a2` under the null. If `s_l` are the singular values associated with the two covariance factors, its null law is the symmetric bilinear Gaussian mixture

\[
U\;\overset{d}{=}\;\sum_l s_l Z_lW_l,
\]

or, equivalently, a signed quadratic-form mixture with weights `+s_l/2` and `-s_l/2`. The two-sided Davies calculation is exact for this fixed-spectrum reference law. Liu is a moment approximation to that law.

This argument concerns the usual `sqrt(P) G sqrt(P)` score construction. The simulation does not challenge the corresponding `CGC` algebraic identity. Its purpose is to measure what happens when the quantities entering that identity are estimated from the same observations used in the score.

## What hyperparameter re-estimation changes

Every feature is fitted separately by INLA in every replicate. The latent precision, observation parameter, working residual, working variance, and nuisance covariance therefore depend on that feature's response. Separate fits preserve independence between the two features, but within a feature `a_j` is not independent of its estimated `H_j`.

Consequently, a Davies p-value that treats the realized spectrum as fixed is a plug-in p-value. The unconditional null law of the full procedure is a mixture over response-dependent spectra and fitted working models. The Monte Carlo experiment estimates the unconditional rejection probability of this complete plug-in procedure; fixed-spectrum theory alone does not establish its finite-sample size.

The fitted Gaussian diagnostics were consistent with the intended experiment. The largest fitted observation-mean error was about `1.22e-16`, and the projected versus raw-constrained operator equivalence error was at most about `3.96e-13`. Median fitted Gaussian dispersions were approximately `0.251` and `0.250`, compared with the true value `0.25`. Median fitted latent precisions were approximately `0.00495` and `0.00515`, with appreciable replicate-to-replicate variation.

## Completed Gaussian pair result

At nominal alpha `0.05`, the completed Davies rejection counts were:

| Variant | Rejections / 500 | Rate | Exact 95% binomial CI |
|---|---:|---:|---:|
| projected | 31 | 0.062 | [0.0425, 0.0869] |
| raw constrained | 31 | 0.062 | [0.0425, 0.0869] |
| raw kernel only | 32 | 0.064 | [0.0442, 0.0892] |
| raw full, nuisance recomputed | 29 | 0.058 | [0.0392, 0.0822] |
| raw full, centered nuisance kept | 237 | 0.474 | [0.4295, 0.5188] |

For the three coherent sensitivity variants, the confidence intervals include `0.05`. The net difference between 31 and 32 rejections does not mean that only one replicate was discordant. Comparing projected and raw-kernel-only Davies decisions replicate by replicate gives 27 rejected by both, 5 rejected only by raw, 4 rejected only by projected, and 464 rejected by neither. Thus 9 replicates were discordant and their imbalance was only one. These results are compatible with nominal size at the resolution of 500 simulations; they neither prove equality nor tightly rule out modest inflation.

At alpha `0.01`, projected, raw-constrained, and raw-kernel-only each rejected 5 of 500 replicates. At alpha `0.001`, each rejected 0 of 500. With 500 replicates, the Monte Carlo standard error is about `0.00975` at alpha `0.05`, `0.00445` at alpha `0.01`, and the expected number of rejections at alpha `0.001` is only `0.5`. The `0.001` result therefore provides almost no tail-calibration evidence.

The corresponding Liu rejection rates at alpha `0.05` were `0.060` for projected, `0.056` for raw-kernel-only, and `0.058` for raw-full with nuisance recomputed. Differences between Liu and Davies reflect approximation as well as Monte Carlo variation.

The raw-full variant that kept the nuisance covariance from the centered null was intentionally incoherent. Its rejection rate of `0.474`, together with 5 Davies fallbacks, confirms that it is an invalid control. It should not be used to choose between projected and raw score kernels.

## Completed NB mean-3 pair result

For `nb3_pair_k6` at alpha `0.05`, the conditioned/projected Davies procedure rejected 30 of 500 replicates (`0.060`). The raw-kernel-only Davies procedure rejected 34 of 500 (`0.068`, exact 95% binomial CI `[0.0476, 0.0937]`), and its Liu procedure rejected 31 of 500 (`0.062`). There were no failures or fallbacks.

The paired Davies decisions comprise 22 rejected by both, 12 rejected only by raw, 8 rejected only by conditioned/projected, and 458 rejected by neither. The raw Davies estimate is numerically above `0.05`, but its interval includes `0.05`; this is evidence of limited precision, not proof that size is controlled. The Gaussian and NB mean-3 raw estimates, `0.064` and `0.068`, point in the same direction. Because this comparison was noticed after inspecting the results and the scenarios have distinct data-generating laws, they should not be pooled post hoc or treated as a prespecified combined test of inflation. Additional replicates or a prespecified joint analysis would be needed to distinguish modest systematic inflation from Monte Carlo variation.

## Completed NB mean-0.3 pair result

The low-count case gives a qualitatively different result. At alpha `0.05`, conditioned/projected Davies rejected 76 of 500 replicates (`0.152`, exact 95% binomial CI `[0.1217, 0.1865]`) and raw-kernel-only Davies rejected 77 of 500 (`0.154`, exact 95% CI `[0.1235, 0.1887]`). Liu rejected at rates `0.154` and `0.156`, respectively. All 500 replicates were valid and there were no Davies fallbacks.

The current plug-in score procedure therefore does not control type-I error in this low-count scenario. The paired Davies decisions were 58 rejected by both methods, 19 rejected only by raw, 18 rejected only by conditioned/projected, and 405 rejected by neither. The nearly identical marginal rates and balanced discordances show that this failure cannot be attributed to removing the projection from the score kernel.

The case-specific DGP files give the same true latent precision, `tau = 0.004481`, for the mean-3 and mean-0.3 cases. The fitted distributions changed substantially:

| Diagnostic | NB mean 3 | NB mean 0.3 |
|---|---:|---:|
| median minimum working weight | 0.463 | 0.107 |
| 5th percentile minimum working weight | 0.242 | 0.051 |
| proportion with minimum weight below 0.1 | 0.004 | 0.438 |
| median `tau1 / tau_truth` | 1.30 | 3.40 |
| median `tau2 / tau_truth` | 1.19 | 1.92 |
| proportion `tau1 > 10 * tau_truth` | 0.034 | 0.464 |
| proportion `tau2 > 10 * tau_truth` | 0.010 | 0.320 |
| median fitted NB size, feature 1 / 2 | 1.99 / 2.01 | 1.81 / 1.85 |
| 5th-95th percentile NB size, feature 1 | 1.45-2.97 | 0.57-8.33 |
| 5th-95th percentile NB size, feature 2 | 1.51-2.73 | 0.72-8.64 |

The low-count tau distributions are especially irregular: their 75th-percentile ratios to truth are about `191` and `142`, while their medians remain `3.40` and `1.92`. This indicates a large upper group of fits assigning very high latent precision. In the mean-0.3 case, fitted NB size was below 1 in `23.4%` and `12.8%` of feature fits, and above 5 in `11.8%` and `15.0%`; none of these threshold events occurred for mean 3. NB `dispersion1` and `dispersion2` are identically 1 in these output files by the family parameterization, so the estimated NB-size columns carry the relevant overdispersion variation.

As a descriptive check only, Pearson and Spearman correlations between `-log10(p)` and each individual fitted tau, NB size, or minimum working weight were small in the low-count conditioned result: all absolute values were below about `0.08`. This does not identify or exclude a mechanism. The distortion may depend on nonlinear combinations of fitted quantities, the response-dependent spectrum, or other aspects of the working-model approximation. The diagnostics establish coexistence of weak information, unstable hyperparameter estimates, and inflated rejection; they do not isolate a causal contribution from any one of them.

The full quantiles, threshold rates, Pearson and Spearman correlations for both score kernels, and paired rejection cells are stored in `artifacts/constraint-type1/estimated/nb-low-count-diagnostics.csv`.

The alpha `0.05` paired counts are summarized below.

| Case | Both reject | Raw only | Conditioned/projected only | Neither rejects |
|---|---:|---:|---:|---:|
| Gaussian pair, `kappa = 6` | 27 | 5 | 4 | 464 |
| NB mean-3 pair, `kappa = 6` | 22 | 12 | 8 | 458 |
| NB mean-0.3 pair, `kappa = 6` | 58 | 19 | 18 | 405 |

## Meaning of the operator variants

- `projected` uses the mean-zero fitted null and the projected score kernel.
- `raw_constrained` is an algebraically equivalent representation and serves as a numerical invariant check. Its agreement with `projected` is not independent calibration evidence.
- `raw_kernel_only` retains the constrained fitted null and nuisance projection and changes only the tested kernel. The helper residualizes the raw score factor against the nuisance design by QR. Since the ideal nuisance projector satisfies `PX = 0`, this leaves the mathematical statistic unchanged and prevents amplification of floating-point leakage along the nearly unpenalized intercept direction.
- `raw_full_recompute_nuisance` reconstructs a raw null covariance with precision estimated under the constrained fit, then recomputes nuisance GLS. It is a sensitivity calculation, not a fully re-estimated unconstrained model.
- `raw_full_keep_nuisance` combines a raw null with a nuisance covariance from the centered null. It deliberately preserves that mismatch as a negative control and receives no QR correction.

The completed pair results support the narrower conclusion that, with a coherent constrained null and nuisance projection, there is not yet decisive evidence that using the raw kernel changes type-I error. The point estimates are too imprecise to establish equivalence. They also show that the current procedure fails badly in the mean-0.3 scenario with either kernel, so score-kernel projection does not repair its low-count calibration. These results do not support removing the mean-zero constraint from the fitted model.

## The marginal experiment has a different null law

The marginal null sets a feature's spatial field exactly to zero. Its statistic is a positive quadratic form,

\[
Q=a^\mathsf{T}a,
\]

whose fixed-working-model reference distribution is a weighted chi-square sum, `sum(lambda_l * chi-square_1,l)`. A right-tail Davies calculation is the relevant fixed-spectrum calibration; Liu is again an approximation. This is distinct from the signed, two-sided bilinear distribution used by the pair test.

The current marginal worker fits the spatial alternative and then plugs its final IRLS quantities and estimated hyperparameters into a nuisance-adjusted, spatial-null score calculation. Under the boundary null, those fitted quantities can be correlated with the statistic. Any finite-sample size distortion therefore belongs to the complete alternative-fit plug-in procedure and should not automatically be attributed to whether the score kernel is projected.

A uniform multiplicative scaling of both the statistic and its spectrum, such as a common kernel precision scale, cancels from the p-value. Estimated fitted means, working variances, negative-binomial size, active observations, and nonuniform changes to the spectrum do not generally cancel.

## Interpreting the negative-binomial cases

For a log-link NB2 model with mean `mu` and size `r`, the IRLS working precision is approximately

\[
W=\frac{\mu^2}{\mu+\mu^2/r}=\frac{\mu r}{r+\mu}.
\]

At mean count `0.3`, many responses are zero and the working information can be very small. The Gaussian working approximation is then discrete and skewed, latent precision and NB size can be weakly identified, and the required `N(0, 3^2)` log-hyperparameter prior can have more influence. A response-dependent cutoff for working weights can also change the active sample and effective rank across replicates. These mechanisms can produce bias or numerical Davies fallbacks without implying that the mean-zero projection itself is wrong.

The NB results should therefore be read with all of the following diagnostics:

- attempted and valid replicate counts, with rejection rates never reported only after silently discarding failures;
- Davies failure and fallback counts;
- minimum and distribution of working weights and the active-observation count;
- fitted NB size and latent-precision distributions;
- the contrast between mean-count `3` and `0.3`, which isolates much of the low-count effect;
- the contrast between `kappa = 6` and `kappa = 0.7`, which changes the strength of the extra raw direction.

A valid-only rejection rate can be selected if numerical failure depends on the response. Treating every failure as a nonrejection gives a lower bound, not a validated size estimate. Both should be reported when failures occur.

## Scope of the Monte Carlo evidence

The binomial intervals quantify simulation uncertainty for each fixed design and data-generating scenario. They are not uncertainty intervals for arbitrary real datasets. Variants and calibration methods are evaluated on the same replicates, so differences should be analyzed as paired rejection indicators rather than by comparing marginal confidence-interval overlap. Shared observation designs across cases also do not create a general calibration guarantee.

The current 500-replicate pair results are useful finite-sample checks at alpha `0.05`. They are sufficient to expose the large low-count failure, whose confidence interval lies well above `0.05`, but more replicates are required for precise comparisons between coherent variants and for meaningful validation at alpha `0.01` or `0.001`.
