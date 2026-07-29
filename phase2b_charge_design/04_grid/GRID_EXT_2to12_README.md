# Full vdW grid, shells 2-12 A (11 shells, 2816 points) - committed 04_grid pipeline, margin 16

## Change from the 2-8A version
Rebuilt SDF with margin 16 (was 10) so shells up to 12A are clean closed surfaces (11/12A clipped the
smaller box -> chi>>2; now all chi=2). Re-extracted/sampled/thinned ALL shells 2-12 on the bigger SDF ->
ONE consistent grid. Built from TRUE R/TS/P geometries (01_geometry), union-of-spheres envelope.

## Points per shell
2A:86 3A:88 4A:116 5A:160 6A:189 7A:229 8A:264 9A:315 10A:368 11A:372 12A:629  (total 2816)
Global 1.5A spacing. K=10 feasible at min_dist 4.7 on EVERY shell (2A fits 11 sites, rest more).

## IMPORTANT - supersedes the 2-8A grid
This 2816-pt grid REPLACES the earlier 1224-pt (margin-10) grid for a consistent 2-12A ladder. The
re-thinning shifted 2-8A counts slightly (2A 86 vs 94 etc). Use THIS grid for the whole ladder. The
2-8A relaxations already running (on the 1224-pt grid) become a consistency preview, not the final run.

## Steps (HPC)
1. scp *_ext.py + grid_ext_points.tsv + run_vpot_grid_ext.sh to 04_grid/
2. cd 04_grid && bash run_vpot_grid_ext.sh   -> dv_grid_ext.tsv (2816 pts, all shells, with Dv)
   (vpot driver already fixed: arg5 = base sp_reactant, member sp_reactant.scfp. ~1-2s per wavefn.)
   SANITY: 2/3/4A max|Dv| must still match the original (0.0092/0.0061/0.0043).
3. On Ubuntu PC (highspy venv): optimise K=10 per shell:
   for s in 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0 10.0 11.0 12.0; do
     python3 optimise_shell.py dv_grid_ext.tsv $s reactant.xyz --K 10; done
4. scp new-shell coords/charges to HPC relax/, relax each (own dir, PBS). shells 9-12 need new PBS
   wrappers (copy run_shell8p0.pbs -> run_shell9p0..12p0.pbs, edit -N and shell name).
