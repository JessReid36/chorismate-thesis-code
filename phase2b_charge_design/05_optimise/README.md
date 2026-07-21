# 05_optimise — Tier-1 certified screen (the novelty)

tier1_sweep.py : certified MILP over the real Dv field (dv_grid.tsv). Two objectives, net-free,
K=1..10, configurable charge menu (any magnitudes/signs), any shell mixed freely.
  A) max-lowering : ceiling; goes unphysical (barrier<0) past K~3 -> frozen-objective validity limit
  B) distributed  : min largest single-charge |contribution|, total held in a data-derived band
Solver: HiGHS (highspy). Certified gap 0.000 on ALL K, both objectives (CBC stalled >90% on the
distributed tail; HiGHS closes it in <1s/solve). Provable gap = the differentiator vs GOCAT's GA.
Key result (bare barrier +17.47, menu +/-1): K=1 distributed / K=2 max-lowering place +1 at the
ether-O valley (3.84 A from O3) = the Arg90 role, with zero enzyme knowledge. See tier1_result_menu_pm1.txt.
NOTE: certified optimality is of the frozen-Dv PROXY, not of catalysis (see ../DESIGN_DECISIONS.md).
Which K / which design actually catalyses is decided at Tier-2 (ranking preservation).
