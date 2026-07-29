# Extended vdW grid (shells 2-8 A) - built with the COMMITTED 04_grid pipeline

## What
Extends the candidate charge grid from 3 shells (2/3/4A) to 7 shells (2-8A) so the shell-ladder can test
K=10 discrete +-1 at incrementally larger distances. Built with the committed 04_grid scripts, only change:
margin 5.0 -> 10.0 in step_2_1 (so the SDF reaches past 8A). Downstream steps unchanged (they take shell
distances as CLI args). Built from the TRUE R/TS/P geometries (01_geometry), union-of-spheres envelope.

## Result
grid_ext_points.tsv: 1224 candidate points. Per shell: 2A 94, 3A 86, 4A 124, 5A 156, 6A 194, 7A 206, 8A 364.
Global 1.5A spacing (Poisson-thinned, pooled across shells). All shells closed surfaces (chi=2).

## IMPORTANT: this is a FRESH consistent grid, not original-331 + new points
Global thinning re-pools ALL shells, so the 2/3/4A counts differ slightly from the original (94/86/124 vs
93/75/163). Use this grid CONSISTENTLY for the whole ladder (don't mix with the old dv_grid.tsv). The
Part-1 jobs already running used the OLD grid's 2/3/4A designs; when comparing, note the grid differs
slightly, OR re-optimise 2/3/4A on this new grid for a fully consistent 2-8A sweep (recommended).

## Steps to finish (on HPC)
1. scp the *_ext.py + grid_ext_points.tsv + run_vpot_grid_ext.sh to 04_grid/ on the HPC.
2. Run vpot:  cd 04_grid && bash run_vpot_grid_ext.sh   (uses committed 02_singlepoints wavefunctions;
   E3 convention: sp_reactant.scfp member + sp_reactant.densities container). -> dv_grid_ext.tsv
   SANITY: the 2/3/4A Dv should be close to the original dv_grid.tsv (same geometry, re-thinned points).
3. Optimise K=10 per shell:  for s in 2.0 3.0 4.0 5.0 6.0 7.0 8.0; do
      python3 05d_shell_ladder/optimise/optimise_shell.py 04_grid/dv_grid_ext.tsv $s 01_geometry/reactant.xyz --K 10; done
4. Relax each (05d_shell_ladder/relax, own dir per job, PBS wrappers shell2p0..8p0 already built).

## Files
step_2_1_sdf_ext.py (margin 10), step_2_2/2_3/2_4_ext.py (load sdf_grid_ext.npz), grid_ext_points.tsv,
grid_final.xyz (viewable), run_vpot_grid_ext.sh
