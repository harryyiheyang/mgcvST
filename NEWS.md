# spdeST 0.0.1.9084

* `spdeST.wgcna()` builds the network from the discoveries of
  `spdeST.test()`: the similarity of a discovered pair is `U / t1`, that of
  every other pair is 0. The ratio of a score and an information that are
  both near zero estimates nothing, and it was the source of the estimates
  far outside `[-1, 1]`. The requirement that every pair of a block be held
  by the test is removed (`keep = "discoveries"` is enough); the network
  returns `pairs`, the number of discovered pairs of the block.

# spdeST 0.0.1.9083

* Removed the compiled entries that only the deleted tests called
  (`mgcvst_spa_quadform_cpp`, `mgcvst_spa_cpp`, `mgcvst_pair_trace_powers_cpp`)
  and `.mgcvst_pair_input_hash()`, which nothing calls.

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
