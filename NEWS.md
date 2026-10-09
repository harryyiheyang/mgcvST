# mgcvST 0.0.1.9029

* The marginal saddlepoint no longer calls Liu: at the mean, where the
  Lugannani-Rice formula is 0/0, it uses the limit 1/2 - rho3 / (6 sqrt(2 pi)),
  rho3 = kappa3 / kappa2^(3/2). A statistic q <= 0 now has saddlepoint p-value
  1 (was NA). Davies p-values are unchanged.

# mgcvST 0.0.1.9028

* The Stage 1 (marginal) Davies calibration now falls back to a saddlepoint
  approximation (Kuonen 1999, Biometrika 86:929) instead of Liu when Davies
  fails, which happens in the extreme upper tail. The saddlepoint has the
  correct exponential tail rate and bounded relative error there (Chen and
  Lumley 2019, CSDA 139:75); single-chi-square moment matching decays too
  fast and understates small p-values. Diagnostics report
  `marginal_method = "saddlepoint"` and `marginal_fallback = TRUE` for these
  features. `mgcvST.marginal(fallback = )` now takes `"none"` or
  `"saddlepoint"`. Davies p-values, Liu calibration when requested, and the
  INLA marginal Liu test are unchanged.

# mgcvST 0.0.1.9027

* Removed internal code left unreachable after the fp16 exact path and
  PCAlearning replaced the double-precision sparse INLA pair pipeline: the
  sparse unit store (`R/unit-store.R`) and the exact-pair chunk-size and memory
  helpers `.mgcvst_inla_pair_chunk_size()` and `.mgcvst_inla_memory_plan()`.
  Fits, scores and p-values are unchanged.
* Tests of the deleted pipeline were removed. The `mgcvST.set()` factor
  interaction test now places its nuisance smooth on a separate covariate,
  because `z + s(z)` is rank deficient and such fits correctly fail without a
  conditional nuisance covariance. The check that approximate `inlaST.test()`
  builds the observation basis once now runs on a fitted INLA model. The dense
  Woodbury check of fitted-model curvature with multiple iid nuisance penalties
  now targets the current compact reconstruction units, which replaced the
  curvature output of the sparse batch kernel.

# mgcvST 0.0.1.9026

* API change: every user-facing `kappa` is now a unit-scale value, and its
  default is `0.05` in both `spde_basis()` (previously `0.1`) and
  `inlaST.set()` (previously `NULL`). The unit length `L` is the largest
  per-axis span, `max - min`, of the observation coordinates supplied to the
  model, and `kappa` is the SPDE scale in coordinates divided by `L`. The
  package converts it internally: the native mesh path uses
  `kappa_internal = kappa / L` on a raw-coordinate mesh, and the mgcv path uses
  `kappa_internal = kappa * s / L`, where `s` is the `spde_mesh()` coordinate
  scale. The same `kappa` therefore gives the same kernel shape whether
  coordinates are recorded in millimetres or micrometres. With `alpha = 2`,
  the practical range is `sqrt(8 * nu) / kappa` unit lengths, about 57 in 2D
  and 40 in 3D at the default.
* Callers that passed a physical or raw-coordinate kappa, for example
  `kappa_mm`, must now pass the unit-scale kappa instead. A script that
  computed `kappa_mm = kappa_unit / L` now passes `kappa_unit` directly; a
  stored physical kappa converts as `kappa_unit = kappa_mm * L`.
* Bases, smooths, models and fits store `kappa_unit`, `unit_length`,
  `coordinate_span` and `kappa_internal`, and their print methods report
  `kappa_unit` and `L`. These fields replace the former `kappa` field. Native
  INLA models also record `spde$range_unit` and `spde$range`.
* mgcvST never estimates kappa. Every SPDE term uses one fixed kappa, so all
  features share one Gaussian-process kernel shape and differ only in variance.
  The joint kappa/tau REML path was removed: `spde_basis(kappa = NULL)` and
  `inlaST.set(kappa = NULL)` are errors, bases no longer carry the three
  projected FEM penalties, and smooths no longer carry `kappa.estimated`.
  Bases saved by earlier versions must be rebuilt with `spde_basis()`.
* `inlaST.set()` accepts `kappa` only with `mesh`; the basis, complete-formula
  and `G` setups take kappa from `spde_basis()`.
* `spde_precision(model, kappa_internal, tau)` renames its `kappa` argument to
  `kappa_internal`, because it takes the mesh-scale value.

# mgcvST 0.0.1.9025

* The default Poisson prescreen threshold is now `1.01` for both mgcv and
  INLA estimation. Explicit thresholds retain their supplied values, and
  `poisson_screen_phi = 0` continues to disable screening. Family targets,
  model fitting and score formulas are unchanged.

# mgcvST 0.0.1.9024

* When INLA's native spatial fit crashes after a successful null fit,
  `inlaST.estimate()` reuses the null fixed and nuisance estimates, sets the
  spatial effect to zero, and assigns spatial precision `1e8` on the original
  FEM scale. The fallback and original error are explicit in diagnostics.
  Existing null p-values and downstream score formulas are retained. Input
  errors and failed null fits retain their existing failure behavior.

# mgcvST 0.0.1.9023

* `rkhs_score_calibrate()` and the fused mgcv Liu pair kernel normalize
  covariance states before calibration. Positive information is no longer
  rejected solely because it is at or below 1e-10 in the input units.
  Squared-score Liu and Davies calibration use the same statistical formulas;
  public scores, information, moments and cumulants retain their original units.
  Davies keeps the existing positive-semidefinite check on the original matrices.
  Pair-result checkpoints use a new calibration signature; existing feature
  states remain reusable and earlier pair results are retained separately.

# mgcvST 0.0.1.9022

* `inlaST.test(approximate_test = FALSE)` scales each gene's reduced curvature
  and projected score before fp16 storage. The equivalent squared-score Liu
  calculation now retains very small or large curvature states, while the
  reported signed score remains on its original scale. The reduced geometry
  and its coverage are unchanged. State and pair checkpoint formats now
  include the scaling contract; earlier checkpoint directories must be replaced.

# mgcvST 0.0.1.9021

* `inlaST.estimate()` now extracts the fixed-effect modes of null models with
  no random blocks. Empty random-block tags are kept empty, avoiding an
  internal latent-mode specification error for intercept-only or covariate
  null models. Nonempty random blocks, offsets, priors, and score calculations
  are unchanged.

# mgcvST 0.0.1.9020

* `mgcvST.estimate()` fits a null model with one parametric coefficient and no
  smooth terms with `mgcv::gam(method = "REML")`. This avoids the one-column
  QR dimension error in `mgcv::bam()` 1.9-4, which made intercept-only null
  scores unavailable. The formula, offsets, family, and score calculation are
  preserved. Other null models and all spatial fits keep their BAM path.

# mgcvST 0.0.1.9019

* `inlaST.test()` now exposes only the squared-score Liu test. The conditional
  Cauchy pair test is removed from the package, together with its R and C++
  implementation, tests, example and HPC benchmark scripts.
* The pair-test selectors `pairwise_method`, `liu_approximation`,
  `calibration` and `conditional_precision` are removed from `inlaST.test()`.
  The new argument `approximate_test` selects the trace evaluation:
  `TRUE` (the default) uses the PCAlearning approximation to the four Liu trace
  moments, and `FALSE` uses the exact fp16 Liu path. Both routes keep their
  numerical implementation and result shapes. `method` still selects the
  multiple-testing adjustment.

# mgcvST 0.0.1.9018

* `inlaST.estimate()` now retains the fitted modes of native iid nuisance
  blocks, such as `s(slide, bs = "re")`, in `nuisance_coefficients` alongside
  the fixed-effect modes. Compact working-state reconstruction therefore
  reproduces the fitted predictor for models with nuisance random effects,
  restoring the intended Stage 2 score covariance for those models.
* `inlaST.test()` now rejects `conditional_precision = "float32"` before
  entering either `score_liu` backend unless
  `pairwise_method = "conditional_cauchy"` is selected explicitly.

# mgcvST 0.0.1.9017

* A `model.set()` feature without a usable conditional nuisance covariance now
  fails at estimation (`mgcvST.estimate()` reports it as
  "nuisance covariance unavailable: <reason>" and it becomes unavailable to
  downstream testing/WGCNA) instead of silently falling back to a per-feature
  eigendecomposition. `.mgcvst_pair_pipeline()`, the davies branch of
  `.mgcvst_test_model()`, and `.mgcvst_model_state_shard()` no longer carry a
  per-feature R-loop fallback for model.set() fits; each now stops with
  "Model score states require the conditional nuisance covariance; re-estimate
  with the current mgcvST.estimate()." when the native dense preparation is
  unavailable. Legacy `spde`-engine fits are unaffected.
* The now-dead model-score chain reachable only through the unused
  `pair_function` parameter of the internal `.mgcvst_test_model()`
  (`.mgcvst_model_pair_single()`, `.mgcvst_model_cached_state()`,
  `.mgcvst_model_score_state()`, `.mgcvst_model_operator()`,
  `.mgcvst_model_operator_legacy()`, `.mgcvst_model_operator_vp()`,
  `.mgcvst_model_apply_P()`, `.mgcvst_model_vsolve()`, and the now-orphaned
  `.mgcvst_full_rank_design()`) is removed.

# mgcvST 0.0.1.9016

* `.mgcvst_pair_pipeline()` and the `calibration = "liu"` branch of
  `.mgcvst_test_model()` build their result columns with preallocated vectors
  written by index instead of quadratic row-wise `data.frame` assignment; the
  numeric results are unchanged.
* Liu-calibrated pairs are evaluated by a single fused, double-precision C++
  kernel (`mgcvst_pair_liu_cpp`): score, the four Liu trace moments, and the
  log-space Liu tail (`log_p_two_sided`, `log_p_positive`, `log_p_negative`)
  in one call per pair chunk, with full-rank `H`. `mgcvST.test()` reports the
  new `log_p_*` columns and adjusts on the log scale, so tails below the
  double range keep their BH/BY decisions.
* `mgcvST.test()`/`mgcvST.wgcna()` no longer accept `inlaST.estimate()` fits;
  use `inlaST.test()`/`inlaST.wgcna()`. The now-unreachable exact fp16
  dispatch inside `.mgcvst_test_model()` and the INLA dispatch inside
  `.mgcvst_wgcna_scores()` were removed. The 0.995 sparse-INLA projection
  coverage is now the single constant `.inlast_projection_coverage`.
* `inlaST.wgcna()` builds its score matrix from the same observation-kernel
  coordinate basis `inlaST.test()` uses, normalized by the basis rank, instead
  of the raw sparse score vectors normalized by `m - 1`.
* `mgcvST.wgcna()`'s mgcv score construction uses a score-only native C++
  kernel (`mgcvst_dense_score_batch_cpp(..., score_only = TRUE)`, batches of
  256 genes) for both the legacy SPDE and model backends; the per-gene R
  fallback and eigendecomposition path were removed.

# mgcvST 0.0.1.9015

* The public test entry points are split. `mgcvST.test()` keeps the mgcv
  interface (exact Liu or Davies calibration, `checkpoint_dir`, `resume`) and
  no longer takes `approximate`, `rank`, `n_per_cell`, `seed`,
  `pairwise_method`, `conditional_precision` or `cache_bytes`.
  `inlaST.test()` calls the internal score engines directly and takes
  `pairwise_method = c("score_liu", "conditional_cauchy")`,
  `liu_approximation = c("exact", "pca_learning")`, `rank`, `n_per_cell`,
  `seed`, `checkpoint_dir`, `resume` and `conditional_precision`. The former
  argument names are not kept as aliases.
* `inlaST.test(liu_approximation = "pca_learning")` projects score
  covariances onto a rank-`rank` basis learned from stratified training genes
  (`n_per_cell`, `seed`), and the four Liu trace moments are obtained by
  contraction with trace tables computed once. Liu log p-values are returned
  in `log_p_two_sided`, `log_p_positive` and `log_p_negative`, and BY
  adjustment is applied on the log scale. The result element `pca_learning`
  stores the training genes, basis rotation, coefficients, per-gene residuals
  and stage timings.
* `inlaST.test(pairwise_method = "conditional_cauchy")` returns the raw
  conditional Cauchy results: `signed_score`, `statistic`, `p_two_sided`, the
  directional `p_1_given_2` and `p_2_given_1`, and their log versions.
  Multiple testing follows `FDR`, `method` and `q.value` as for Liu pairs
  (`p_adjusted`, `log_p_adjusted`, `discovered`); BY is no longer forced, and
  the `S`, `p`, `p_BY` and `BY_reject` columns are replaced.
* The real-gene landmark trace-CUR approximation (`approximate = TRUE` or
  `"landmark"`, `n_ref`, `ref_method`, `ref_seed`, `ref_tol`,
  `diagnostic_pairs`) and its native score-state streaming were removed.
* `cache_bytes` is an internal argument of the score engines; the default
  `NULL` probes available memory.
* Liu pair tests prepare each required gene once and evaluate native pair
  batches. Adaptive caches use available system, cgroup and Slurm memory.
  Gene blocks reduce repeated shard reads. `checkpoint_dir` preserves gene states and completed
  exact pair batches for resumption.
* Feature preparation uses dynamic OpenMP scheduling with single-thread BLAS.
  The shared INLA 0.995 observation-kernel projection is unchanged.
* Persistent checkpoint fingerprints use bounded in-memory blocks instead of
  writing the complete fit to a temporary file. Runs without persistent
  checkpoints skip the content fingerprint.
* See `inst/notes/pair-pipeline.md` for scope and validation. Held-out exact
  comparisons are development checks, not mandatory production work.

# mgcvST 0.0.1.9014

* Sparse INLA Liu pair tests now use a score-only constrained observation-kernel
  basis retaining at least 0.995 of the eigenvalue sum. The full sparse INLA
  fit and its nuisance construction are unchanged.
* A native MAGIC 3D check (97,830 observations, q = 1,962, three pairs) kept
  r = 1,407 directions (coverage 0.9950176). Direct reduced scores agreed
  with the corresponding projected full scores; relative to full-q Liu, the
  three log10-p differences were -0.0582, 0.0382, and -0.1614. The first
  three-pair call took 10.92 s including a 6.72 s basis setup; the reduced
  four-trace kernel took 1.43 s versus 3.86 s for full q. These timings do not
  establish an end-to-end speedup.

# mgcvST 0.0.1.9013

* INLA's variational-Bayes mean and variance correction is always disabled.
  The estimator uses the joint latent mode and expected-Fisher covariance, so
  the posterior-marginal correction is outside the fitted score-test contract.

# mgcvST 0.0.1.9012

* Native INLA models accept any number of categorical random-intercept and
  full-rank Gaussian-process nuisance smooths. GP blocks are projected off the
  intercept and whitened to native iid blocks before fitting and score testing.
* INLA fixed effects now use explicit zero precision, and nuisance covariance
  is reconstructed from the expected Fisher information with every fitted iid
  penalty retained in the small coefficient block.
* INLA latent modes are projected onto the model's supplied constraints
  without rejecting fits at an additional residual threshold. Unprojected
  and projected residuals remain available in diagnostics.

# mgcvST 0.0.1.9011

* Native INLA models support one `s(group, bs = "re")` iid nuisance term.
  Its precision penalty is retained in the small coefficient covariance block
  shared with fixed effects, in both marginal and pairwise score calculations.
* INLA fits disable configuration retention and reconstruct expected-Fisher
  nuisance covariance through sparse solves. Joint latent modes are read from
  the retained mode vector.
* INLA thread settings are preserved unless supplied explicitly. INLA tests
  no longer reset the caller's thread environment.
* Parallel INLA workers receive only their feature chunk, including matching
  offsets and Poisson routing. Temporary fit lists are released before scoring.
* INLA fitting and score batches omit repeated static precision checks and
  full working-matrix scans. Required numerical factorizations retain their
  failure handling.
* Pairwise calibration continues to treat cross-feature iid effects as
  independent. The mgcv estimation and testing paths are unchanged.
