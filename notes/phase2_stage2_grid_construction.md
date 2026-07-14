Stage 2 - candidate charge-grid construction.

Candidate positions for the external point charges are generated on constant-distance shells
around the substrate's van der Waals surface, following surface curvature and guaranteeing a
minimum inter-point spacing. The pipeline is: (i) an analytic union-of-spheres signed distance
field, SDF(r) = min_atoms(|r - r_atom| - R_vdW), with Bondi radii, evaluated on a voxel grid;
(ii) constant-offset shells extracted as triangular meshes by marching cubes at additive offsets
d = 2.0, 3.0, 4.0 A (absolute distances, not multiples of atomic radius - contact distances are
absolute, so additive offsets are more physically motivated than the multiplicative scaled-vdW
convention of ESP-fitting grids); (iii) area-proportional sampling of each mesh into a dense,
areally-uniform point cloud; and (iv) a single global Poisson-disk (blue-noise) thinning of the
pooled cloud to a minimum spacing of 1.5 A. Thinning is global rather than per-shell so the
spacing guarantee holds across shells (adjacent shells are ~1 A apart radially, so per-shell
thinning alone leaves near-degenerate cross-shell pairs); each retained point keeps its shell
tag for later per-shell analysis.

For the placeholder reactant geometry this yields 326 candidate positions with a verified global
minimum spacing of 1.500 A. Shell offsets are reproduced to +/- 0.001-0.003 A and each shell is a
closed manifold (Euler characteristic 2). The grid is converged with respect to voxel resolution:
refining from 0.30 A to 0.20 A changes shell areas by <0.2% and the candidate count by <1% (326 vs
329), confirming the grid is a property of the molecular geometry rather than the discretization.
The minimum-spacing constraint both ensures selected charges are physically distinct and
conditions the downstream optimizer (the inter-point Coulomb matrix, which appears as the
quadratic term in the charge-selection problem, is well-conditioned at 1.5 A spacing).
