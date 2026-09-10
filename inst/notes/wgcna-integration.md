# WGCNA source reconciliation — 2026-09-10

The package release source is the Git checkout at
`C:/Users/yxy1234/Downloads/mgcvST`, with remote
`https://github.com/harryyiheyang/mgcvST.git` and branch `main`.

## Version correspondence before integration

| Location or version | Commit | Content |
| --- | --- | --- |
| GitHub previous version and local Git checkout | `d73591ceb7a79ff2e9d735924ac396c6000180da` | HPC package replacement and SOCK worker library verification |
| GitHub latest version used as the integration base | `480365f74dad2a167fb263c4aea05230d5547165` | Constrained sparse INLA estimation and score workflow |
| Temporary manuscript audit checkout | `f7557135a98467a135ed2a7fc7b96f26a7363056` | Older shared-design implementation; two commits behind the integration base |
| `magicST/mgcvST` research source | No Git metadata | Earlier package snapshot with local WGCNA development; not a complete copy of either recent GitHub version |

The canonical local checkout had no tracked modifications before integration.
Its untracked R history and existing package archive were preserved. The local
branch was fast-forwarded from `d73591c` to `480365f` before adding WGCNA.

## WGCNA content recovered

The local WGCNA source was recovered from `magicST/mgcvST/R/wgcna.R`, together
with its help page and regression checks. That source has SHA256
`0EF851FE7E503F7FFDEEF4528BA3DD4131B215AA8A90AB92E9706D90C9A5E796`, matching
`magicST/hpc_wgcna/source/wgcna.R` and the earlier transition record.

Version `0.0.1.9005` exports `mgcvST.wgcna()` and its print method, declares
WGCNA, dynamicTreeCut and fastcluster as optional dependencies, documents the
workflow in the README, and includes package regression tests. The public
gene-block and parameter interfaces preserve the earlier local implementation.
Each selected block uses uncentered score covariance `crossprod(A) / q` and
retains its requested gene order.

The score extraction is integrated with the current package backend so that
WGCNA uses the same fitted nuisance adjustment and score coordinates as the
current mgcv and INLA workflows. This is a necessary adaptation: the earlier
local WGCNA source reconstructed an older score operator. Saved results from
that research source therefore retain their original implementation provenance.

WGCNA requests score vectors only. The dense path applies the current fitted
operator to the working error; the sparse INLA path uses the existing constrained
sparse solver and returns before constructing the pair-test calibration matrix.
The default full sparse score state used by pair testing remains unchanged.

The INLA and HPC updates remain in place. Existing Graph-SuSiE functionality
also remains available. Historical research scripts, cached fits, simulation
outputs and the frozen HPC WGCNA source remain unchanged. The earlier local
removal of Graph-SuSiE is not propagated as part of this additive release.

## Validation

On Windows with R 4.6.1, the source package passed `R CMD check --no-manual`
with zero errors, warnings or notes. With `NOT_CRAN=true`, the full suite
reported 762 passing expectations and four skips requiring an external
`MGCVST_BASELINE` checkout. The WGCNA-specific tests passed all 46 expectations,
including a real sparse INLA estimate, current conditional-Vp scores, manual
network reconstruction, reordered and overlapping blocks, multiple score
groups, old compact fits, and sparse score-only/full-state equivalence.

An independent synthetic check at 300 score coordinates found identical
score vectors and widths for the sparse full-state and score-only paths.
The default full-state return fields remain unchanged.

## Development workflow

Use the canonical Git checkout for subsequent package changes. Fetch and
reconcile remote work before editing, then validate, commit and push completed
updates in the same task. Use the commit hash to identify a release across
computers; a loose source directory or an installed package is not a Git
version identifier.
