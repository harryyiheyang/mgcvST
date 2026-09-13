# Mesh figure QA

The figure shows how a prespecified spatial size field reallocates tetrahedral
resolution within two node budgets. It is a geometry comparison, with a 2 by 2
panel layout: uniform and adaptive meshes at budgets of 1,500 and 2,800 nodes.
Each panel is the actual intersection of all tetrahedra with z = 0.5 mm;
polygons are obtained by linear edge-plane intersection. No intersecting
tetrahedron is sampled out. The dashed circle marks a prespecified region,
not a boundary fitted to expression. The domain is a synthetic 3 by 2 by 1 mm
box, not an anatomical reconstruction.

R was the saved plotting backend. The figure was drawn and exported in R,
using ggplot2, patchwork, cairo_pdf, svglite and ragg. All panels have the same
physical axes and scales. The four complete panels and assembled PNG were
visually inspected: labels, region outlines and mesh edges are legible and
do not collide. There are no stochastic summaries or uncertainty bars in
this geometry figure. Actual node counts, rather than budget labels, appear
in each panel title.

The final canvas is 183 by 145 mm. PDF and SVG retain editable text. The PNG
is a 300 dpi report preview; this is not a journal submission raster bundle.
The static source audit returned 16 passes, four warnings and no failures.
R successfully executed the plotting script, resolving its syntax warning;
explicit device dimensions resolve the validator's unrecognized-width warning.
The remaining warnings concern optional TIFF and 600 dpi submission output.

The simple PDF Tf scanner incorrectly reports 1 pt for Cairo's unit-font
representation, which applies font sizes through the text transformation
matrix. A transform-aware inspection of the actual PDF characters found
effective sizes of 7, 8, 9, 9.6 and 11 pt. Rotated y-axis labels have an 8 pt
transformation; their character bounding-box heights are not font sizes.
Thus all rendered glyphs meet the 5 pt floor. Source coordinates are retained
in `mesh-section-source.csv`; the raw meshes are reproducible with the meshing
script and its versioned parameter settings.
