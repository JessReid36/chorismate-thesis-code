# Phase 2b reconciliation — full arc since commit ad09dfe (GOCAT comparison)

This documents everything built/run since the last commit, INCLUDING dead-ends, so results and rationale
are reproducible. Directory-by-directory below. Each has its own README with purpose / method / outcome / next.

## The narrative arc (why each step happened)
1. **05b_optimise_realism** — certified MILP designs (K2-K10, a_pot/a_force preorg penalties) relaxed under
   field. FINDING: discrete +-1 charges alone can't recover near-attack geometry (open ~4A or fragment).
   Includes the continuous-charge + QP exploration (dead-end: linear objective rails to +-1; quadratic QP
   gives fractional but modest novelty).
2. **r2SCAN-3c calibration** (in 05b/step3_higherK) — tested cheaper functional. FAILED the geometry gate
   (SIE fragments the anion). DEAD-END: DFT-only. Kept as a documented negative (methods-chapter point).
3. **05c_cradle** — neutral steric cage to hold geometry where charges can't. FAILED: frozen LJ walls leak,
   substrate opens past the cage. Both single mechanisms (charges-alone, sterics-alone) now fail.
4. **Research report** — certified non-GA optimiser survey (convex QP/SOCP/MISOCP/Lasserre). Confirms
   novelty; ranks convex-QP-with-physical-regulariser as the path if discrete fails.
5. **05d_shell_ladder** — THE big experiment. Discrete K=10 +-1 per vdW shell, swept 2-15A (+ 0.1A fine
   7-9A), to answer "is there a distance where discrete charges catalyse without rupturing?"
   Grid extended 2->15A via committed 04_grid pipeline (margin 5->16). vpot validated (Dv matches original
   to ~1%). FINDING: NO shell works - close ruptures, far pushes open (net-charge imbalance). 8A barrier
   +35.51 (anti-catalytic, DOUBLE bare +17.47). Negative controls (bare + 2A q=0) confirm bodies inert,
   effects are field. -> forces the fractional QP.
6. **all_barriers** — exhaustive barrier scan of EVERY shell (in progress) to make the negative unimpeachable.

## Current status at reconciliation
- Barrier jobs (all shells 2-15A + fine) RUNNING. collect_barriers.sh assembles the table.
- Next after barriers confirm: build the certified fractional-charge QP (convex, proximity-weighted or RESP
  regulariser) - the certified negative forces it.

## Commit intent
Everything here is [retained-final] EXCEPT clearly-marked [exploratory]/[dead-end] items, which are kept
deliberately as documented negatives (they explain why we did what we did).
