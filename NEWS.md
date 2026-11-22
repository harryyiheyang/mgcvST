# spdeST 0.0.1.9082

* Removed what the package no longer uses: the testthat suite with `testthat`
  and `CompQuadForm` in Suggests, the estimation signature stored in the fit
  of `spdeST.estimate()` (nothing read it), `.mgcvst_store_has()`, the figure
  `man/figures/inla3d-snap25.png`, and the ignore rules of directories that no
  longer exist. The entries of earlier development versions are removed from
  this file; they described mechanisms that have been deleted.
* No function of the package writes a checkpoint. `spdeST.test()` keeps its
  score states and pair shards in `tempdir()` (which follows `TMPDIR`) and
  removes them when the call returns or fails.
* All p-values are computed by saddlepoint approximation.
