The RDS contains the 17 unique BAM pairs identified by
`inst/validation/inla-bam/bam-pair-degeneracy-cases.csv`: four Gaussian spatial
null, eleven negative-binomial spatial null, and two low-count examples.
The original 9005 validation library reconstructed their saved `a` and `H`
states from `artifacts/inla-bam-validation/inference/*/bam-compact-fit.rds`.
No model was refitted. Scores and information were checked against
`inst/validation/inla-bam/inference-paired-results.csv`.

Each record retains its original missing Liu/Davies result, score, information,
simulation identifiers, both covariance matrices and score vectors. Reference
p-values use the 9022 implementation on the same states divided by each
matrix's maximum absolute entry, with the score divided by the corresponding
geometric mean. This equivalent normalization yields finite results without
changing either statistical formula. The original information lies between
1.07544917434142e-11 and 6.30525702731001e-11.

The artifact SHA256 hashes and extraction/replay scripts are recorded in the
research workspace under
`paper_workspace/05_analysis/null_score_revalidation/score_calibration_scale_bug`.
