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
