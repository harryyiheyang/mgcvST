# PCAlearning approximate Liu path

## Settings

- Branch `pca-learning-liu` from `32ead29`.
- Entry point: `inlaST.test(fit, approximate = "PCAlearning", rank = 10, n_per_cell = 3, seed = 1)`.
- Package compile flags unchanged (no AVX2/FMA).
- Machine: i7-14700K (20 cores, 28 logical processors), 66 GB RAM, Windows 11, R 4.6.1.
- Threads: kernel verification and audit timings in this worktree used 28 threads; the
  1,762-gene run and the 550-gene reproduction run used 20 threads; experiment
  (`experiments/approx-liu-p`) timings used 20 threads; R-level Liu and BY are single-threaded.

## Per-step cost, experiment versus package

Training set S (300 genes), r = 10, q = 1404, eval set T (11,175 pairs).

| Step | Experiment: form, time (20 threads) | Package: form, time (28 threads) | Equivalence |
|---|---|---|---|
| Projection | full q x q M returned to R and stored (15.8 MB double per gene), R crossprod: 0.00076 s/gene after 0.196 s/gene preparation + materialization | c_j, ||H_j||_F^2 and float32 packed copy formed inside the OpenMP materialization; M not returned: 56 genes 10.01 s versus 9.93 s for materialization alone | c_j max rel 1.2e-13; ||H||^2 max rel 3.0e-14; packed copy byte-identical; a identical |
| Gram / basis | double V (4.4 GB): state read 45.9 s, Gram 1.59 s, eigen 0.01 s, B 0.47 s | float32 packed V (1.10 GB), double accumulation: read + pack 21.3 s, Gram 0.47 s, B 0.20 s | Gram max rel 1.95e-9; eigenvalues 1-10 max rel 1.1e-8 vs experiment; B'B - I 9.3e-14; projector difference 1.1e-6; c_j vs experiment max rel 6.9e-9 |
| Trace tables | r^4 tuple products + cyclic Gram, float32 GEMM: 249 s (REPORT), 137.6 s re-measured; 2.7 GB peak | symmetric coefficient products (N, W) + Gram, float32 GEMM with double accumulation: 43.4 s at 28 threads (reviewer: 80.7 s at 20 threads under load); 5.4 GB peak | max rel vs experiment tables: Tsym1 3.0e-14, Tsym2 1.9e-7, Tsym3 3.4e-7, Tsym4 4.0e-7 |
| Pair contraction + Liu | per-pair dot products 3.45e-8 s/pair + R Liu 4.75e-7 s/pair (1 thread) | block GEMM, Liu inside the kernel: 5.7e-8 s/pair (7.5e-8 s/pair with t1..t4 returned) | t1..t4 vs experiment max rel 1.7e-10, 5.2e-7, 1.9e-6, 8.6e-6; vs double brute force (20 pairs) 4.9e-15, 3.4e-7, 3.6e-7, 2.4e-6; U vs exact 1.4e-12 |
| Liu alone | R `pchisq`: 4.5e-7 s/pair (1 thread) | log-space kernel: 5.0e-7 s/pair (1 thread), 3.9e-8 s/pair (28 threads) | 16,125 exact pairs: p max rel 2.7e-14 vs R; log p max rel 6.2e-15 vs experiment `liu_lp`; one-sided max rel 2.7e-14 |
| BY, 6.25e7 log p | not measured | `p.adjust(exp(lp), "BY")` 7.02 s; log-space step-up 4.65 s | identical decisions |

Audit scripts and log (session scratchpad, not in the package): `testA/t1_liu.R`, `t1b_nc.R`,
`t1c_nc.R`, `t2_real.R`, `t2_real.log`, `t2a_tables_syn.R`, `t3_stream.R`, `t4_pairs_by.R`.

Exact path (`approximate = "none"`): p-values change only in the noncentral far tail, where the
C++ log-space Liu is more accurate than R `pchisq(ncp)` (reviewer: old R relative error 3.5e-3 at
p 1e-14 to 1e-50 for ncp < 80; 0 or off by 3.3 in log p for ncp >= 80). Central-branch agreement
with R: 2.4e-15 relative (p >= 1e-5) to 1.7e-13 (p to 1e-280). All 16,125 exact T and bench100
pairs are in the central branch (ncp = 0).

## 1,762-gene run (20 threads)

1,773-gene fit (1,762 sub-fit genes + 11 S genes); all 1,551,441 pairs of the 1,762 genes.

One-time per-gene cost:

| Stage | Genes | Seconds |
|---|---|---|
| Fit read | 1773 | 2.2 |
| Sparse preparation + observation basis | | 6.67 |
| Sampling | | 3.26 |
| Training materialization + pack | 300 | 57.44 |
| Gram | 300 | 0.53 |
| Eigen | | 0.03 |
| Basis B | | 0.19 |
| Trace tables | | 46.68 |
| Projection (non-training materialization) | 1473 | 314.61 |

Pairwise cost:

| Stage | Units | Seconds |
|---|---|---|
| Pair contraction + Liu | 1,551,441 pairs | 0.18 |
| BY, log space, three directions | 4,654,323 p | 0.30 |
| BY, `p.adjust` reference, three directions | 4,654,323 p | 0.40 |

`inlaST.test` wall 431.93 s; peak process working set 11.0 GB; R heap max used 6.6 GB.

Non-training genes (1473): e2_relative median 0.0022, max 0.072.

## Full 11,184-gene fit (HPC): pending

## Reproduction of experiments/approx-liu-p (S, r = 10)

Training set: the experiment's S (300 genes); 11 S genes absent from `sub-fit.rds` were
added from the full fit.

| Eval | Pairs | max abs dlog10 p vs experiment approx | dlog10 p vs exact q01 | q99 |
|---|---|---|---|---|
| T (550-gene run) | 11175 | 2.6e-7 | -0.0103 | 0.0197 |
| bench100 (550-gene run) | 4950 | 2.5e-7 | -0.0100 | 0.0202 |
| T within 1,762-gene run | 10585 | 2.6e-7 | -0.0102 | 0.0197 |

BY discordance versus experiment approx: 0 in two-sided, positive, and negative for all three rows.

BY versus exact:

| Eval | Direction | PCA rejections | PCA only | Exact only |
|---|---|---|---|---|
| T | two-sided | 6012 | 0 | 0 |
| T | positive | 3779 | 0 | 0 |
| T | negative | 2241 | 0 | 0 |
| bench100 | two-sided | 2351 | 0 | 1 (of 2352 exact) |
| bench100 | positive | 1366 | 0 | 0 |
| bench100 | negative | 986 | 0 | 0 |
