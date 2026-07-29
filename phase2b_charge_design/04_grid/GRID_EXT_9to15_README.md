# Shell ladder EXTENSION: shells 9-15 A (run now, in parallel with 2-8 already running)

## What & why
Extends the discrete K=10 +-1 ladder further out to test the "still too close?" hypothesis: do the far
shells pull C1-C6 INTO the window as we ease off (field was holding it open -> 4a), or plateau at the mu=0
floppy control (~5.0 -> too weak, 4b), or hold near bare reactant (~3.1-3.5 -> 4c)? 9-15A resolves it.

## Points (grid_ext_9to15_points.tsv, 3579 pts, idx 9000+)
9A:321 10A:363 11A:433 12A:484 13A:537 14A:559 15A:882. K=10 feasible at min_dist 4.7 on all.
Built with committed 04_grid pipeline, margin-16 SDF (box 42A, all shells chi=2), true R/TS/P geometry.

## Steps
1. HPC: scp grid_ext_9to15_points.tsv + run_vpot_9to15.sh + *_ext.py to 04_grid/
   cd 04_grid && bash run_vpot_9to15.sh    -> dv_grid_9to15.tsv (3579 pts, Dv)
   (vpot driver already fixed: sp_reactant.scfp member + sp_reactant base. Same convention, verified.)
2. Ubuntu PC (highspy venv): optimise K=10 per shell:
   for s in 9.0 10.0 11.0 12.0 13.0 14.0 15.0; do
     python3 optimise_shell.py dv_grid_9to15.tsv $s reactant.xyz --K 10; done
3. scp shell9p0..15p0 coords/charges + relax_pbs/*.pbs to HPC 05d_shell_ladder/relax/
   Relax each (own dir, dry-run, submit) - same flow as 2-8. All 7 in parallel.
4. analyse: analyse_relaxed.py per shell. Read C1-C6 trend 8->15A:
     decreasing toward window = 4a (hypothesis CONFIRMED, field was holding open, distance is the fix)
     increasing toward 5.0     = 4b (too weak, floppy) | flat near 3.1-3.5 = 4c (held near reactant)

## NOTE - DEFERRED unified re-optimisation (user-flagged)
The 2-8A jobs currently running use the OLD 1224-pt (margin-10) grid; 9-15 use the margin-16 grid. These
are physically consistent (same geometry/wavefunction, non-overlapping shells) but thinned separately. Once
TARGET shells are identified from this ladder, RE-OPTIMISE ALL shells on ONE unified grid for the final
clean result. The full 2-15A grid is saved as grid_full_2to15_points.tsv (4713 pts) ready for that
re-vpot + re-optimise. Do NOT mix dv_grid.tsv (old) and dv_grid_9to15.tsv for a single design - they're
separate grids; the unified re-run resolves it.
