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
