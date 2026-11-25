# spdeST 0.0.1.9085

* New `spdeST.clustertest()` tests clusters of features against single
  features. The score vectors of the members of a cluster are summed,
  `a_S = sum_g a_g`, and the score of the cluster with a feature `c` is
  `U = a_S' a_c`. Given `a_S`, the score of a feature whose field is unrelated
  to those of the members is a linear form of `a_c` with variance
  `a_S' M_c a_S`, whatever the dependence among the members, so its p-value
  is the normal tail of `U / sqrt(a_S' M_c a_S)` (the saddlepoint
  approximation of a linear form of a Gaussian vector). The test holds for
  clusters chosen from other features, such as RNA modules of
  `spdeST.wgcna()` tested against proteins or chromatin features. Each row
  also gives the pooled estimate `rho = U / tr(M_S M_c)`, the within-cluster
  estimate `rho_cluster` and `k = rho / sqrt(rho_cluster)`, the correlation of
  the feature with the common factor of the cluster. `exclude` leaves chosen
  members out of a cluster for a feature (for example the RNA of the gene of
  a chromatin window), and `groups` combines sets of features against every
  cluster: pooled estimates and the Cauchy combination of their p-values.
  The products run in two compiled kernels.

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
