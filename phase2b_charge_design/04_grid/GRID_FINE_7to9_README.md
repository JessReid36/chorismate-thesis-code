# Fine shell sweep 7.0-9.0 A in 0.1A steps (21 shells) - completeness around the 8A discrete optimum

## Why
8A was the closest-to-window intact discrete design (C1-C6 3.655, just 0.155A outside). This sweep checks
if any sub-A shell between 7-9A dips INTO the window. Expectation (from mechanism): unlikely to cross, but
makes the discrete negative unimpeachable. Each intact shell then gets a real barrier (see barrier8/).

## Grid
grid_fine_7to9_points.tsv: 1041 pts, 21 shells (7.0,7.1,...,9.0). idx 20000+. Built committed pipeline,
margin-16 SDF, true R/TS/P geometry. NOTE: 0.1A spacing + 1.5A global thinning = uneven per-shell counts
(18-139); all still fit K=10 at min_dist 4.7 (sparsest 8.3A: 18 sites, ~13 fit). Fine for a trend sweep.

## Steps
1. HPC: scp grid_fine_7to9_points.tsv + run_vpot_fine.sh to 04_grid/; bash run_vpot_fine.sh -> dv_grid_fine.tsv
2. Ubuntu PC: for s in 7.0 7.1 ... 9.0; do python3 optimise_shell.py dv_grid_fine.tsv $s reactant.xyz --K 10; done
   (21 designs; use a loop: for i in $(seq 0 20); do s=$(python3 -c "print('%.1f'%(7+0.1*$i))"); ...; done)
3. scp coords/charges + PBS to HPC relax/, relax all 21 in parallel (own dir each).
4. analyse_relaxed per shell; the intact ones -> barrier via barrier8/ protocol.

## Then
Whichever intact shell is closest to / in the window -> compute its real barrier. If NONE lower the barrier
below bare (+17.47) -> discrete +-1 confirmed non-catalytic at all distances -> fractional QP.
