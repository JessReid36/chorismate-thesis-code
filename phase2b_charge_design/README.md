# phase2b_charge_design

Rebuild of the Phase 2 electrostatic charge-field catalyst design on the VALIDATED,
consistent-region reaction geometry — replacing the stale-geometry phase2_charge_design.
Methods are inherited unchanged; every number is regenerated. See NAMING_CROSSWALK.md
for the stage-by-stage old->new mapping and the superseded list.

## Why the rebuild
phase2_charge_design is methodologically sound but every result sits on a stale geometry:
01_substrate held the NEB-TS (break/form 2.173/2.547), not the frequency-verified
reduced-region OptTS (2.111/2.526), and R/P used a 1959-atom active region while the TS
used 102 atoms. phase2b re-seeds the identical pipeline on the consistent-region triple:
  TS       -> 05_qmmm/18c_reduced_region (OptTS, -313.30 cm^-1)
  Reactant -> 05_qmmm/18e_reactant_reduced (reduced-region min)
  Product  -> 05_qmmm/18f_product_reduced  (reduced-region min)

## Theory level (fixed)
ORCA 6.0.1, B3LYP-D3BJ / def2-SVP / def2/J RIJCOSX / CPCM eps=4, charge -2, mult 1.

## Pipeline
00_admin        central paths + reacting-atom map + provenance
01_geometry     committed extraction of R/TS/P substrate (24 atoms) from Phase 1
02_singlepoints CPCM(eps=4) R/TS/P single points -> densities
03_dvpot        Dv = V_TS - V_R (Sokalski catalytic field) via orca_vpot
04_grid         candidate charge grid (union-of-spheres SDF -> shells -> sample -> thin)
05_optimise     certified MILP + distributed min-max + whole-path NEB min-max
06_polarize     self-consistent polarization (Route 1 density-relaxed validation)
07_validate     Route 1 headline + D1-D5 diagnostics + native-BsCM positive control
sandbox         optimiser-formulation prototyping (geometry-independent; runs now)

## Status
Scaffold. 01_geometry waits on Phase 1 18e/18f convergence; sandbox proceeds now.

## Design decisions
See DESIGN_DECISIONS.md for the two-tier structure (certified screen + relaxed validation),
the GOCAT-grounded novelty positioning, and the ranking-preservation check.
