# B: arbitrary-coordinate prediction on a fixed mesh

New-coordinate prediction is interpolation, not SPDE construction or fitting.
It uses only the saved raw-to-mesh transform, existing triangles, barycentric
weights and coefficient-space projection. Any outside-mesh point is an error.
The same internal `.spde_basis_at(basis, loc, pc)` implements both predictors.

Full SPDE uses `A_new %*% Z`. PC prediction uses the saved
`pc_mesh_projection = Z %*% V_q`, directly computing `A_new %*% (Z V_q)`.
Neither a dense observation-by-mesh A nor a dense full projected PC intermediate
is formed. No custom C++ spatial index or interpolation loop has been added.

Training PC values intentionally retain the previous `B %*% V_q` arithmetic
through CppMatrix, not the reassociated newdata multiplication. Exact training
coordinates return the cached basis; exact subsets, row reordering and mgcv
prediction blocks return cached rows using exact hexadecimal-double keys.
No tolerance-based coordinate rounding is used to classify new locations.

The reduced PC fit/prediction cache does not replace `score_basis = basis$B` or
`pc_score_Q`: pairwise score geometry stays in the full projected coordinates.

`predict.gam(fit, newdata=..., type="link"/"response"/"lpmatrix")` therefore
works with genuinely new in-mesh coordinates, including one-row prediction
blocks. Offset and linear terms continue to be handled by mgcv itself.

`geometry::tsearchn()` already uses `tsearch(..., bary=TRUE)` with the default
quadtree backend for 2D input (geometry 0.5.2). This release retains the existing
tsearchn call. The benchmark separately measures it against explicit quadtree,
then measures coordinate transform, sparse A creation, A*Z, A*(ZV), and end-to-end
full and PC prediction at 100, 1,000, 10,000 and 100,000 new locations.

Run `inst/benchmarks/projector.R` with `MGCVST_BENCH_OUT` pointing to a result
directory. Stage measurements are separate median wall times, not components
that necessarily add up exactly to the independently measured total. Prediction
totals include coordinate-key lookup and validation. The benchmark uses random
strictly interior barycentric samples from the saved MISO mesh and records
agreement between both search calls and between direct PC projection and the
manual `(A Z) V` reference.

The newdata reference uses numerical tolerance because associativity and
sparse/dense multiplication can change rounding. Training compatibility and
the separate A hot-path tests require strict identical results instead.
