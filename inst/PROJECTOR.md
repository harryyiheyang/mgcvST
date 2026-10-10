# B: arbitrary-coordinate prediction on a fixed mesh

New-coordinate prediction is interpolation, not SPDE construction or fitting.
It uses only the saved raw-to-mesh transform, existing triangles, barycentric
weights and coefficient-space projection. Any outside-mesh point is an error.
The internal `.spde_basis_at(basis, loc)` implements the predictor.

Prediction uses `A_new %*% Z`; no dense observation-by-mesh A is formed. No
custom C++ spatial index or interpolation loop has been added.

Exact training coordinates return the cached basis; exact subsets, row
reordering and mgcv prediction blocks return cached rows using exact
hexadecimal-double keys. No tolerance-based coordinate rounding is used to
classify new locations.

`predict.gam(fit, newdata=..., type="link"/"response"/"lpmatrix")` therefore
works with genuinely new in-mesh coordinates, including one-row prediction
blocks. Offset and linear terms continue to be handled by mgcv itself.

`geometry::tsearchn()` already uses `tsearch(..., bary=TRUE)` with the default
quadtree backend for 2D input (geometry 0.5.2). This release retains the existing
tsearchn call. The benchmark separately measures it against explicit quadtree,
then measures coordinate transform, sparse A creation, A*Z, and end-to-end
prediction at 100, 1,000, 10,000 and 100,000 new locations. The benchmark
script predates the removal of the `spdePC` smooth and still measures a
principal-component prediction that no longer exists.

Run `inst/benchmarks/projector.R` with `MGCVST_BENCH_OUT` pointing to a result
directory. Stage measurements are separate median wall times, not components
that necessarily add up exactly to the independently measured total. Prediction
totals include coordinate-key lookup and validation. The benchmark uses random
strictly interior barycentric samples from the saved MISO mesh and records
agreement between both search calls.

The newdata reference uses numerical tolerance because associativity and
sparse/dense multiplication can change rounding. Training compatibility and
the separate A hot-path tests require strict identical results instead.
