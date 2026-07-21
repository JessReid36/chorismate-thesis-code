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
