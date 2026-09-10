# Low-count INLA score inflation: mechanism review

This review concerns the `nb03_pair_k6` null experiment: 200 observations, true NB size 2, mean count 0.3, two independent nonzero mean-zero SPDE fields, and zero cross-feature covariance. The completed estimated-hyperparameter procedure rejected about 15% at nominal 5% with either the conditioned or raw score kernel. The near equality of those two rates rules out score-kernel projection as the explanation for the large inflation.

## Score quantities and where hyperparameters enter

For feature `j`, the reference score has

\[
a_j=F_{s,j}^{\mathsf T}P_j e_j,
\qquad
M_j=F_{s,j}^{\mathsf T}P_jF_{s,j},
\]

and the pair statistic is `U = a1' a2`. With fixed coherent working models, its symmetric null law is determined by the singular values formed from `M1` and `M2`.

Latent precision enters two different places. Scaling only the score factor from `Fs` to `c Fs`, while holding `P` fixed, multiplies `a` by `c` and `M` by `c^2`. The observed score and every signed-mixture weight then receive the same pairwise scale, so the Davies p-value is unchanged. This is a useful invariant check. Precision can change the p-value through the null factor inside

\[
V=D+F_0F_0^{\mathsf T}
\]

and hence through `P`; that effect is not a uniform score rescaling.

For NB2 with log link, expected Fisher working quantities are

\[
D_{\mathrm{exp}}^{-1}=\frac{r\mu}{r+\mu},
\qquad
z_{\mathrm{exp}}=\eta+\frac{y-\mu}{\mu}.
\]

Thus fitted NB size changes `D`, while fitted latent eta changes both `D` through `mu = exp(eta)` and the working residual through `z`. Swapping eta, size, or the residual separately can create a working state that no fitted likelihood produced. Such swaps are useful algebraic diagnostics but must be labeled as operator ablations rather than alternative estimators.

The negative observed Hessian and Newton working response are

\[
W_{\mathrm{obs}}=\frac{r\mu(y+r)}{(r+\mu)^2},
\qquad
z_{\mathrm{obs}}=\eta+
\frac{(y-\mu)(r+\mu)}{\mu(y+r)}.
\]

Their expectations recover the Fisher weight. These formulas provide a posterior-native candidate if INLA's fixed-effect covariance comes from the same observed joint curvature.

## Actual posterior covariance and score calibration

There is a more direct construction when the full spatial posterior covariance is available. Write the constrained field as `u = Z gamma`, let `Qc = Z'QZ` and `Rc = chol(Qc)` with `Qc = Rc'Rc`, and define the prior-standardized innovation

\[
\xi=\sqrt{\tau}R_c\gamma.
\]

For a linear Gaussian working model with flat fixed nuisance effects, standard Gaussian conditioning gives

\[
a=E(\xi\mid y)=\widehat\xi,
\qquad
M=I-\operatorname{Cov}(\xi\mid y).
\]

If `Vp_u` is the spatial block of the full joint INLA posterior covariance, already marginalized over the fixed nuisance block, then

\[
\operatorname{Cov}(\xi\mid y)=
\tau R_c Z^{\mathsf T}V_{p,u}Z R_c^{\mathsf T}.
\]

This absorbs fixed-effect uncertainty through the full joint inverse and avoids separately combining an expected working `D` with a posterior nuisance-only covariance.

For the observation-centered raw score, let `Rraw = chol(Q)` and `H = I - 1 g'`, where `g'1 = 1`. The innovation transform is

\[
T=R_c Z^{\mathsf T}H R_{raw}^{-1},
\qquad
a_{raw}=T^{\mathsf T}a,
\qquad
M_{raw}=T^{\mathsf T}MT.
\]

Because `H` maps into the constraint plane, `ZZ'H = H`. Hence the raw centered factor is also

\[
F_{raw,c}=\tau^{-1/2}AHR_{raw}^{-1}
=\tau^{-1/2}\{AR_{raw}^{-1}-1(g^{\mathsf T}R_{raw}^{-1})\},
\]

a rank-one column-centering operation that needs no explicit `Z`.

The bounded validation in `inst/benchmarks/inla-posterior-score.R` used four fixed-hyperparameter Gaussian fits. Posterior versus exact-P discrepancies in projected `a/M` were at most about `1.30e-4/1.25e-4`, while the no-Z raw factor identity held to `6.66e-16`. In four cached NB fits, posterior `M` was much closer to the observed-working construction (`1.5e-6` to `7.4e-5`) than to expected working in the worst cases (up to about `0.0236`). This validates the geometry as a candidate; it does not validate NB tail calibration.

The user-target nuisance covariance is INLA's own posterior fixed-effect covariance, rather than an expected-Fisher reconstruction. Suppose it is inserted into

\[
B=V^{-1}-V^{-1}X C_{\beta}X^{\mathsf T}V^{-1}.
\]

Unless `C_beta = (X' V^-1 X)^-1` for the same `V`, `B X` is not zero. Moreover, if `a = Fs' B e`, its fixed-working covariance is

\[
\operatorname{Cov}(a)=F_s^{\mathsf T}BVB^{\mathsf T}F_s,
\]

not generally `Fs' B Fs`. The latter equality holds for the exact GLS residual projector. Therefore an actual-posterior-Vp implementation should report all of the following before using the central signed-mixture reference:

- `max(abs(B X))`, so nuisance leakage is visible;
- the difference between naive `Fs' B Fs` and sandwich `Fs' B V B' Fs`;
- actual posterior Vp versus exact GLS covariance for the selected expected or observed working state;
- p-values from naive and sandwich spectra in paired simulations.

Using actual posterior Vp with expected `D` mixes observed/Laplace posterior uncertainty with expected working curvature. Pairing actual Vp with `Wobs` and `zobs` may make the components more coherent, especially under empirical-Bayes hyperparameter integration, but it remains response dependent. The Gaussian signed-mixture law is then a plug-in approximation and requires Monte Carlo calibration.

The current raw-kernel helper QR-residualizes the raw factor against `X`. This is algebraically inert only when the defining `P` exactly annihilates `X`. With actual posterior Vp and `P X != 0`, QR residualization changes `a` and `M`. The unmodified raw factor and the QR-residualized factor must therefore be reported as distinct variants. The conditioned `Z` kernel is the clean primary comparison.

## Recommended paired ablations

The smallest interpretable refit experiment is a two-by-two design on the same cached responses:

| Latent precision | NB size | What changes |
|---|---|---|
| estimated | estimated | observed baseline |
| true fixed | estimated | removes precision estimation while eta adapts |
| estimated | true fixed | removes NB-size estimation while eta adapts |
| true fixed | true fixed | removes both hyperparameter estimation paths |

Each cell should retain its own fitted eta and actual posterior Vp. This design measures the full fitting pathway and allows an interaction between tau and size.

A post-fit operator factorial can then use exact GLS Vp to isolate algebraic pathways:

1. choose tau in the null covariance from the estimated or true value;
2. choose NB size in `D` from the estimated or true value;
3. choose eta and working residual from the estimated-hyperparameter or fixed-true-hyperparameter fit;
4. recompute `D`, `e`, `V`, and exact GLS Vp coherently within each cell.

The score-factor-only tau swap should be added as an invariant: it should leave p-values unchanged when `P` is fixed. Cross-swapping `D` and `e` can follow as a clearly labeled diagnostic if the coherent factorial indicates that the eta pathway matters.

Finally, repeat selected cells with actual posterior Vp under both expected and observed working states, reporting unmodified conditioned/raw factors and both naive and sandwich `M`. This separates the hyperparameter pathway from posterior-covariance/working-curvature compatibility.

## Fixed-true-hyperparameter 100-pair probe

`inst/benchmarks/inla-lowcount-oracle-null.R` replayed the cached original responses for replicates 1 through 100. It fixed latent precision at the saved DGP value `0.0044810478` and NB size at 2, while INLA still fitted the constrained latent field. The positive-hyperparameter configuration retained the required `log(parameter) ~ N(0, 3^2)` entries and used no PC prior.

This probe uses the legacy expected-P raw-kernel construction. It is an experimental mechanism ablation, not the user-target posterior-Vp method.

| Procedure on the same 100 responses | Rejections | Rate | Exact 95% binomial CI |
|---|---:|---:|---:|
| fixed true tau and NB size | 3 | 0.03 | [0.0062, 0.0852] |
| original estimated tau and NB size | 10 | 0.10 | [0.0490, 0.1762] |

The paired rejection cells were: both 1, fixed-only 2, estimated-only 9, and neither 88. All 100 oracle fits and Davies calculations were valid, with no fallback. The maximum fitted observation-mean error was `6.84e-17`, and the maximum constraint residual was `3.88e-17`.

The paired shift from 10 to 3 rejections makes hyperparameter estimation a strong candidate contributor for this construction. With only 100 pairs, the fixed-hyperparameter confidence interval remains wide and includes 5%. This probe neither proves calibrated size nor determines whether tau, NB size, or the changed latent eta is responsible; the factorial ablation is needed for that separation.

## Relation to the old bam result

The earlier bam low-count experiment reported approximately 2.4%-2.6% rejection among 499 valid replicates. It is not yet a matched benchmark for this question. Its Snow-versus-serial RNG streams generated different responses, and its score combined expected-Fisher working quantities with the fitted `bam` covariance block, whose observed-versus-expected interpretation is under audit. It should not be used as evidence that one estimator controls size until the response matching and covariance construction are corrected.
