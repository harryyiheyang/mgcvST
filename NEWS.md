# mgcvST 0.0.1.9032

* Estimation and memory release (the estimation side of the saddlepoint
  roadmap). The Stage 2 pair calibration is unchanged: Liu moment matching of
  the four trace moments.
* API change: `mgcvST.estimate()` and `inlaST.estimate()` estimate in two
  steps. Step 1 fits only the null model of every feature and computes its
  Stage 1 null-first p-value (Davies, then the saddlepoint); the p-values are
  adjusted by `adjust` (`"BY"` by default, or `"BH"`, `"Sidak"`, `"none"`)
  into `diagnostics$marginal_q_value`. Step 2 fits the spatial model of the
  features chosen by `spatial` only: `"discoveries"` (default; Stage 1
  q-value at most `q.value`, default 0.05), `"all"`, `"none"`, feature IDs,
  one-based indices or a logical vector. The new arguments are `spatial`,
  `adjust`, `q.value`, `checkpoint_dir` and `resume`; they follow the
  existing ones, so positional calls keep their meaning. A null fit keeps no
  working vector.
* API change: `mgcvST.estimate_spatial()` and `inlaST.estimate_spatial()` add
  spatial models for more features to a step 1 fit and return the updated fit.
  The responses are checked against per-row digests stored at step 1, the
  controls, offset and family routing of the original call are reused, and
  the supplied fit is not changed.
* Features without a spatial model have `spatial_selected` or
  `spatial_fitted` `FALSE` in the diagnostics and missing working quantities,
  and are unavailable to the pair tests and to WGCNA. `pairs = NULL` tests
  every pair among the spatially fitted features, and such features are not
  reported as failures in `$failed`; an explicit pair with an unselected
  feature gets status 3 and is explained in `$failed`. `mgcvST.test()` and
  `inlaST.test()` stop with a clear message for a fit without any spatial
  model.
* Resumable estimation: with `checkpoint_dir` each worker saves its chunk of
  either step when it completes, and a repeated call resumes from the saved
  chunks. A chunk is keyed by its step, its features and the digests of their
  responses; the directory manifest records the estimator, the fit format
  and a signature of the model, offset and controls. A directory written by
  another estimator, by a version before this one, or for another model,
  offset or control is refused. Without `chunk_size` a run with a checkpoint
  directory uses chunks of at most 50 features.
* INLA memory: the workers compute the null fits, the Stage 1 p-values, the
  score vectors `a_j` and the mean of the fitted mean `mu_bar` themselves and
  return compact per-feature results only. The manager no longer receives or
  holds the observation-length working vectors (`working_error`,
  `working_variance`, `eta`, `mu`) of any fit. `mu_bar` is stored on the fit;
  PCAlearning reads it instead of recomputing the working state.
* The INLA observation basis is full rank: all `q - 1` directions of the
  constrained field are kept, ordered by their observation-kernel eigenvalue,
  and the 0.995 eigenvalue-coverage truncation is gone. `inlaST.test()` and
  `inlaST.wgcna()` use the same basis, and the WGCNA normaliser is `q - 1`.
  The fit records the basis kind and rank at estimation (`basis_spec`), the
  tests refuse a different basis, and the basis kind is part of the pair
  checkpoint signature, so pair checkpoints written with the truncated basis
  are not reused. The test metadata reports `basis_kind` instead of
  `target_coverage`.
* Fits estimated before this version: an `inlaST.estimate()` fit lacks
  `mu_bar` and the two-step bookkeeping, and `inlaST.test()`,
  `inlaST.wgcna()` and `inlaST.estimate_spatial()` refuse it with a message
  to re-run the estimation. An `mgcvST.estimate()` fit of the earlier format
  holds every quantity the mgcv tests read and is still accepted.
* Measured estimation memory (`inlaST.estimate()`, negative binomial, `n = 40000`
  observations, `m = 144` mesh nodes): a 9031 manager held 1.27 MB of
  observation-length vectors per null fit and per spatial fit (`4 n 8` bytes
  each, 1.25 MB of the 1.27 MB), that is 35 GB for the MAGIC dimensions
  (`n = 97,818`, 11,184 features). The 9032 manager receives 14.9 KB per null
  fit and 7.0 KB per spatial fit, of which `2 m 8` bytes are the score vector
  and the target coefficients; for MAGIC (`m = 1962`) that is about 0.6 GB for
  all features and about 0.25 GB for 10% of them, in the final fit as well as
  in transit. A worker holds the vectors of at most 16 features at a time. On a
  200-feature case (`n = 40000`, `m = 36`) the end-to-end peak above the
  starting memory fell from 837 MB to 583 MB; the rest is the data, the model
  and the per-chunk copies of `Y`, which are unchanged.
* Examples, tests and the README call the estimators with `spatial = "all"`
  where they need spatial fits for every feature.

# mgcvST 0.0.1.9031

* Stage A of the saddlepoint release (cleanup and API). The Stage 2 pair
  calibration is unchanged: the pair p-values are still Liu moment matching of
  the four trace moments.
* API change: `mgcvST.test()` and `inlaST.test()` share one argument list,
  `(fit, pairs, q.value, adjust, threads, chunk_size, checkpoint_dir, resume,
  verbose)`; `inlaST.test()` adds `rank`, `n_per_cell` and `seed` after
  `verbose`, so that a positional call means the same in both tests. Removed: `BPPARAM`, `calibration`, `approximate_test`,
  `liu_approximation`, `cache_bytes`, `FDR`, `method`, `highlight`, and the
  `...` that reached `cache_bytes`. The fp16 exact INLA Liu path
  (`approximate_test = FALSE`) and the mgcv `calibration = "davies"` path are
  deleted; `inlaST.test()` always uses PCAlearning. `rank`, `n_per_cell` and
  `seed` take their defaults from one constant.
* API change: `pairs = NULL` tests every pair of the available features in
  both tests. Pairs are generated and scored in blocks and streamed to
  shards; no pair matrix or per-pair table is allocated.
* API change: one compact result for both tests. Each pair is a row of the
  integer feature indices `i < j`, `score`, the natural-log two-sided,
  positive and negative p-values, the adjusted two-sided `log_q`,
  `remainder_kind` (0: Liu, no remainder) and `status` (0 evaluated, 1 invalid
  trace moments, 2 invalid p-value, 3 a gene without a usable score state).
  A pair with a status other than 0 has missing log p-values and is not
  adjusted. Rows are written as Parquet shards while the pairs are evaluated.
  `$results` is the same table sorted by `(i, j)` when it fits the memory
  guard (56 bytes per pair, 20% of available memory) and is `NULL` otherwise.
  Without a `checkpoint_dir` the shards are temporary: they are deleted once
  `$results` is built and `$shards` is empty; when `$results` is `NULL` they
  stay in `tempdir()` and `$shards` lists them. With a `checkpoint_dir`,
  `$shards` lists the final files. Feature names are looked up as `$feature_id[i]`. The
  19-column per-pair table with character columns, `information`,
  `effective_rank`, `p_*` columns and the highlight/retained flags is gone;
  genes without a usable state are listed in `$failed`.
* Multiple testing is one adjustment of the two-sided family, computed in log
  space by a native kernel: `adjust = "BY"` (default; Benjamini-Yekutieli
  under arbitrary dependence, `c(m) = sum(1 / i)`, equal to
  `stats::p.adjust(, "BY")`), `"BH"`, `"Sidak"` or `"none"`. p-values down to
  `exp(-5000)` keep their order and decisions; the former non-BY path
  exponentiated before adjusting and underflowed. Positive and negative
  discoveries are the adjusted two-sided discoveries split by the sign of the
  score. When 24 bytes per pair exceed 40% of the available memory, the
  adjustment is skipped with a warning and `log_q` is `NA`.
* Checkpoints are keyed by an algorithm contract (`calibration_contract`,
  route, `k`, remainder order, basis sha, kernel version). A checkpoint
  directory that holds pair results written under another contract, including
  every pair block written before this version, and a PCAlearning manifest of
  an earlier contract, are refused with an error instead of being resumed.
  Feature score states keep their signature and remain reusable; delete the
  `pairs-*` directories to reuse them.
* The two mgcv test engines are one. `mgcvST.estimate(Y, G)` with a raw
  `gam(..., fit = FALSE)` setup is converted by `mgcvST.set(G = G)` and fitted
  by the model path, so the mgcv score always uses the conditional nuisance
  covariance `Vp`; the former fixed pseudo-inverse projection of the raw-`G`
  path is gone. The two projections agree up to the accuracy of `bam`'s `Vp`:
  the signed score of the Mapt-Map1b pair of `MISO_E13` (baseline mesh,
  `kappa = 0.1`) moves from 28.295 to 28.182 (two-sided p-value 9.0e-07 to
  9.9e-07), and a toy fixture moves by up to 1.5e-3 (relative). The raw-`G`
  path no longer returns `fit_basis`, `fit_penalty` or `basis_metadata`, and
  accepts Gaussian, negative-binomial, Poisson and quasi-Poisson families.
* The `spdePC` smooth is removed, with `spde_basis(pc_cutoff = )`, the
  principal-component fields of a basis, the score-precision PSD branch and
  `spdePC_g999`-based examples (the data set keeps its finer mesh under that
  name). The mgcv examples use the full-basis `spde`.
* Removed internals: `.mgcvst_liu_summaries()`, `.mgcvst_test_chunk()`,
  `.mgcvst_test_spde()`, `.liu_squared_score()`,
  `.rkhs_covariance_score_direct()`, the Davies-only helpers and worker bundle
  entries, the `state_store` argument, the float32 score-state storage,
  `mgcvst_pca_pack_cpp()` and the seven `mgcvst_fp16_*_cpp()` exports. The
  shared INLA working-state kernels moved from `src/inla_fp16.cpp` to
  `src/inla_working_state.cpp`. `rkhs_score_calibrate()` and
  `rkhs_covariance_score()` lose `method = "davies"`. Also removed:
  `.mgcvst_capture_marginal()`, `.mgcvst_model_sparse_constrained_solver()`,
  `.mgcvst_inla_test_pairs()` (folded into the test driver), the test-only
  `cache_bytes` argument of the exact pair pipeline, the stale metadata names
  `liu_approximation` and `storage`, and the example scripts
  `inst/examples/null.R` and `inst/examples/alternative.R`, which depended on
  removed functions. `inst/examples/inla_estimator.R` runs on the current API.
* Review fixes after Stage A. A PCAlearning resume follows the block schedule
  (including `chunk_size`) stored in the pair directory, as the exact route
  does, so a resume with another `chunk_size` reuses every completed shard. The
  pair directory and its contract are opened after the basis exists in both
  routes; stale directories are still refused before any work. In the exact
  route a pair with an invalid p-value (status 2) has missing log p-values and
  can no longer count as a discovery. `$results` is built in one pass into
  preallocated columns and ordered by one permutation, and a temporary run
  deletes its shards and pair directory after it. The memory guard collects
  garbage before it probes. The universe hash of an explicit pair list is
  computed in bounded pieces. The single-step Sidak adjustment is exact for
  log p below -700. `threads` must be an integer.

# mgcvST 0.0.1.9030

* Stage 1 (the marginal score test) has one calibration in both the mgcv and
  INLA branches: Davies, then the saddlepoint approximation when Davies fails.
  Davies fails if it errors, returns a missing or non-finite p-value, or
  returns Qq <= 0 or Qq > 1. `ifault` is still reported but no longer
  decides, so a Davies p-value in (0, 1] with `ifault = 1` is kept.
* The INLA marginal test calibrates the eigenvalues of the full-space
  curvature 0.5 (M + M') above 1e-12 times the largest, in place of Liu on
  four trace moments. M is built as before and the run time is unchanged.
  Diagnostics report `marginal_method = "davies"` or `"saddlepoint"` and the
  `marginal_fallback` flag, as in the mgcv branch.
* API change: Liu is removed from Stage 1. `mgcvST.marginal()` loses its
  `calibration`, `fallback` and `threads` arguments, and
  `mgcvST.estimate(marginal_args = list(method = "liu"))` is an error. The
  marginal Liu helper and the marginal-only C++ moment kernel
  `mgcvst_marginal_liu_moments_cpp()` were removed. Stage 2 pair tests,
  including their Liu calibration, PCAlearning and the fp16 exact path, are
  unchanged.

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
