# Phase 2b directory guide (what each step is, outcome, where it led)

## Committed before this session (00-08)
00_admin, 01_geometry (R/TS/P true geoms), 02_singlepoints (R/TS/P wavefns for vpot),
03_dvpot (original dv_grid.tsv), 04_grid (candidate grid pipeline), 05_optimise (certified MILP core),
06_polarize, 07_stage0 (bare-field implosion finding), 07_validate, 08_relaxed_validation,
relaxation_attempts (the surrogate-method development log). DESIGN_DECISIONS, TIER2_KNOWN_ISSUES, GOCAT
comparison, NAMING_CROSSWALK.

## Added this session (05b, 05c, 05d + grid ext)
### 05b_optimise_realism — certified designs relaxed under field; the "charges alone" tests
- K2-K10 x {a_pot,a_force} preorg-penalty designs, relaxed. OUTCOME: discrete +-1 can't recover near-attack
  geometry (open ~4A or fragment). step3_higherK/ has the K6/K10 results + r2SCAN calibration.
- r2SCAN-3c calibration: FAILED (SIE fragments anion) -> DFT-only. [documented dead-end]
- continuous/QP exploration: linear objective rails to +-1; quadratic QP gives fractional (kept for the
  fractional-optimiser stage). [exploratory, informs next]
- OUTCOME -> led to 05c (try sterics) and 05d (try distance).

### 05c_cradle — neutral steric cage (the "sterics alone" test)
- heuristic 6-bead LJ cage + certified set-cover cage. OUTCOME: FAILS, substrate opens past the cage.
  Both single mechanisms now fail. [dead-end, documented]
- OUTCOME -> led to 05d (distance sweep) and the fractional-optimiser plan.

### 05d_shell_ladder — THE distance experiment (main result)
- K=10 discrete +-1 per vdW shell, 2-15A + 0.1A fine 7-9A. Grid extended via committed 04_grid pipeline.
- controls: bare (no field/bodies) + 2A q=0 (bodies inert check). CONFIRM bodies inert, effects are field.
- OUTCOME: NO shell catalyses. Close ruptures; far pushes open (net-charge imbalance). 8A barrier +35.51
  (anti-catalytic). all_barriers/ = exhaustive barrier scan of every shell.
- OUTCOME -> forces the certified fractional-charge QP (next).

### 04_grid additions — extended grid 2-15A
- step_2_*_ext.py (margin 16), grid_ext_points / dv_grid_ext / dv_grid_fine, run_vpot_*.sh. vpot validated
  (2/3/4A Dv matches original to ~1%).

## Next (not yet built)
- Certified fractional-charge QP (convex, proximity-weighted or RESP regulariser) — the certified negative
  from 05d forces it. Research report ranks the options.
