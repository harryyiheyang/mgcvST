# Intercept centering in the score test

## The score identity relevant here

Let `1` be the observation-space intercept, let `C = I - 1 1' / n`, and let
`G` be any symmetric spatial score kernel. The null residual precision is

```
P = V^{-1} - V^{-1} X (X' V^{-1} X)^{-} X' V^{-1}.
```

When `X` contains the intercept, `P 1 = 0` and, by symmetry, `1' P = 0`.
Hence

```
P C = P,       C P = P,
P (C G C) P = P G P.
```

For the symmetric positive-semidefinite square root this is equivalently

```
sqrt(P) C G C sqrt(P) = sqrt(P) G sqrt(P),
```

because `sqrt(P) 1 = 0`. Thus observation-space centering of the tested
kernel is redundant in the score statistic and its null spectrum. This
identity is exact and does not require uniform locations, equal IRLS weights,
constant row sums of `G`, or a special SPDE boundary geometry. Numerically,
the only qualification is the error with which the fitted operator satisfies
`P 1 = 0`.

The same result can be expressed with a factor `G = F F'`. Replacing `F` by
`CF = F - 1 colMeans(F)` does not change `F' P e` or the nonzero spectrum of
`F' P F`. The nonzero eigenvalues of `F' P F` are also those of
`sqrt(P) G sqrt(P)`. Therefore QR residualization of the raw score factor
against an intercept, as used by the `raw_kernel_only` reference
implementation, is a stable way to calculate the same score. It is not an
additional spatial `Z` projection and does not change the test.

With a general nuisance design, subtracting any component of `F` in
`col(X)` preserves `F' P e` and `F' P F`, because `P X = 0`.

## A different object: conditioning the Gaussian field

The fitted latent field may separately be constrained to have observation
mean zero. If raw mesh coefficients have precision `Q`, `A` maps them to
observations, and `g=A'1/n`, conditioning on `g'u=0` gives

```
Kcond = K - K 1 (1' K 1)^{-1} 1' K,
K = A Q^{-1} A'.
```

In general `Kcond` is not `C K C`. Consequently the earlier condition based
on whether `P K 1` vanishes answers whether the conditional Gaussian
covariance may be replaced by the raw covariance. It does not answer whether
the intercept component `C G C` may be omitted from a tested score kernel.
Applying that condition to the user's score identity was a category error.

This distinction separates two implementation choices:

- The null marginal covariance `V` must describe the fitted null model. For
  a pair test, zero cross-covariance still permits nonzero marginal spatial
  variance in each feature, so its fitted mean-zero field constraint remains
  part of `V` and hence of `P`.
- Once that correct `P` is fixed, centering the tested cross-covariance kernel
  with `C` is redundant. A raw sparse factor may be used, with QR removal of
  its nuisance-column component for numerical stability.

For a single-feature spatial-variance null, the tested field is absent under
the null. Its score kernel may likewise be used as `G` or `C G C` when the
null nuisance design includes an intercept. This statement concerns the
score tangent after nuisance residualization; it does not assert that raw and
conditional Gaussian field distributions are identical.

## Type-I interpretation

Under the Gaussian null, using the correct `P` and matching `H=F'PF`
calibrates the score. The identity predicts identical statistics and spectra
for `F` and `CF`, up to floating-point error. Type-I simulation remains useful
for fitted hyperparameters and non-Gaussian working approximations, but it is
not needed to prove the intercept-centering identity.

Replacing the constrained null marginal covariance by an unconstrained one
is a separate experiment. It changes `P`, may break `P V P = P` for data from
the fitted constrained null, and can change calibration. Results from that
experiment should not be interpreted as evidence against omitting `C` from
the tested kernel while retaining the correct null `P`.
