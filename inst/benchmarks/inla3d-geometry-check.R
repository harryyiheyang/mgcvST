library(Matrix)
library(fmesher)
options(warn = 2)

xyz <- rbind(c(0, 0, 0), c(1, 0, 0), c(0, 1, 0), c(0, 0, 1))
mesh <- fm_mesh_3d(xyz, matrix(1:4, nrow = 1L))
E <- fm_fem(mesh)
D <- rbind(c(-1, -1, -1), diag(3))
stopifnot(abs(E$ta - 1 / 6) < 1e-12,
          max(abs(diag(E$c0) - 1 / 24)) < 1e-12,
          max(abs(as.matrix(E$g1) - tcrossprod(D) / 6)) < 1e-12)
loc <- rbind(c(0.1, 0.2, 0.3), c(0.25, 0.25, 0.25))
A <- fm_basis(mesh, loc = loc)
stopifnot(max(abs(as.matrix(A) - cbind(1 - rowSums(loc), loc))) < 1e-12)
kappa <- 5
Q <- kappa^4 * E$c0 + 2 * kappa^2 * E$g1 + E$g2
Q1 <- fm_matern_precision(mesh, alpha = 2, rho = 2 / kappa, sigma = 1)
stopifnot(max(abs(Q1 * (8 * pi * kappa) - Q)) / max(abs(Q)) < 1e-12)
cat("PASS: tetrahedron volume, 3D stiffness, barycentric weights, and Matérn scaling.\n")
