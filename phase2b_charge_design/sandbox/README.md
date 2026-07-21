# sandbox — optimiser-formulation prototyping

Geometry-independent: the charge optimiser operates only on a Dv grid + inter-point
distances, so we settle the formulation here BEFORE the final grid exists (while Phase 1
runs). Prototype and compare, on a placeholder/synthetic Dv grid:

1. certified max-lowering MILP   — reproduce; confirm gap 0.000; confirm it concentrates
                                    to K=2 and drives K>=3 to a (nonphysical) negative barrier
2. distributed min-max           — band-constrained reformulation making K=3-8 physical
3. whole-path NEB min-max        — minimise the max per-image contribution across path
                                    images; stays linear so the certificate survives (the
                                    home for the multi-charge question)

Output: a settled, validated optimisation formulation ready to take the real Dv grid.

## Validated (synthetic-field prototypes)
- synth_field.py    — single-image grid + multi-image path generators (illustrative magnitudes)
- formulations.py   — certified_milp, distributed_minmax, wholepath_minmax (+ single-point scorer)
- demo_compare.py   — reproduces the findings; run: `python3 demo_compare.py` (needs python-mip)

Findings (gap 0.0000 throughout):
1. Certified MILP reproduces K=2 ddE ~= -10.7 (real -10.888) and K>=3 barrier going negative.
2. Distributed min-max: max|contrib| falls 4.67 -> 1.22 across K=2-8 (real 4.68 -> 1.18).
3. Single-point is DECEPTIVE for multi-charge (crushes the TS image while another image is the
   true max); whole-path min-max gives the honest, bounded, certified barrier and beats it for K>2.
   -> whole-path is the correct multi-charge objective.

Real-data requirement (whole-path): per-image Dv, i.e. orca_vpot on each NEB-band image density
(from 18_nebts) at CPCM eps=4, not just the two endpoints. New step for 03_dvpot.
