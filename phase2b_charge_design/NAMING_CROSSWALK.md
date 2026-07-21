# Phase 2 -> Phase 2b naming crosswalk (back-track record)

phase2b rebuilds phase2_charge_design on the validated consistent-region geometry.
Use this table to trace any phase2b stage back to its phase2 origin.

| phase2b          | phase2 (old)                        | change |
|------------------|-------------------------------------|--------|
| 00_admin         | 00_admin                            | reuse; paths repoint to 01_geometry |
| 01_geometry      | 01_substrate                        | REPLACED: committed extraction from 18c/18e/18f (was manual NEB-TS swap) |
| 02_singlepoints  | 02_baseline + 03_dvpot(step_1_1b)   | R/TS/P CPCM(eps=4) single points |
| 03_dvpot         | 03_dvpot (step_1_2)                 | reuse engine, new geometry |
| 04_grid          | 04_grid                             | rebuild grid on new R/TS/P envelope |
| 05_optimise      | 05_select                           | reuse optimiser; + whole-path NEB min-max |
| 06_polarize      | 06_polarize                         | reuse; = Route 1 density-relaxed validation |
| 07_validate      | 08_validate + Route 1 diagnostics   | enzyme overlay + D1-D5 + native-BsCM control |
| sandbox          | 05_select experiments               | formulation prototyping |
| (dropped)        | 10_surrogate                        | demoted to corroboration only |

## Superseded in phase2 (NOT for final numbers)
- 01_substrate R/TS/P geometries (stale NEB-TS 2.173/2.547)
- all *_realgeom result dirs (built on the stale geometry)
- frozen-geometry field set (-10.888 / -11.354), polarised (-12.201 / -14.009), K-sweep, design_K2-8.pc
- shell_restricted_sweep.py (retracted)
- +15.771 (frozen) and +15.08 (mixed-region) barriers -> superseded by the 18e-consistent barrier (pending)

## Two-tier update (see DESIGN_DECISIONS.md)
| phase2b            | role   | notes |
|--------------------|--------|-------|
| 05_optimise        | Tier 1 | certified screen (THE novelty): MILP + linear whole-path min-max |
| 08_relaxed_validation | Tier 2 | relaxed-path validation on top candidates; not claimed novel |
07_validate remains the frozen-geometry Route-1 headline + D1-D5 + native-BsCM check.
