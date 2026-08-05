# 06_convex_cutplane — certified convex QP + cutting-plane charge optimiser

Novel optimiser for the external-charge catalyst-design problem. Controlled comparison against
GOCATs genetic algorithm and our prior discrete MILP: the EVALUATION harness is held fixed, only the
charge optimiser changes.

## Controlled variable
oracle/ is a PINNED, UNCHANGED copy of the GOCAT-style evaluator (gocat_screen.py + neb_pass.py). It is
the single fixed measurement. Do not edit it here; if it must change, re-pin and note it. Any barrier
difference vs the MILP/GA is then attributable to the optimiser alone.

## Representation
Abstract fractional point charges on the 331-site vdW grid (inputs/grid_final.xyz), each a single
LJ pseudo-atom behind the existing LJ wall (GOCAT-style charged quasi-atoms), replacing molecular
surrogates. The oracle embeds fractional charges unchanged.

## Inner objective (convex)
minimise  sum_i q_i * Dv_i  +  lambda * ||q||^2   (Dv_i = V_TS - V_R, from inputs/dv_grid.tsv)
subject to box |q_i| <= q_max, and (if system_is_charged) net charge sum_i q_i = 0 + low-monopole cap.
Solve to certified optimality (HiGHS QP; gap 0). This is iteration 0; a_i is recomputed at the relaxed
geometry each outer iteration via the dvpot engine.

## Reaction-agnostic I/O contract
inputs: reactant.xyz, ts.xyz, product.xyz, grid_final.xyz, dv_grid.tsv, and flags {system_is_charged,
q_max, lambda}. No reaction-specific hard-coding in solver/; reacting-atom indices live only in the
oracle integrity checks, not the optimiser.

## Status
Scaffold only. solver/ (convex_solve.py + cutplane_loop.py) to be written next.

## Oracle patch
gocat_screen.py carries one fix over the 05d original: the post-convergence RMSG block reads the gradient from the ash run() return value (energy, gradient) instead of a non-existent qmmm.grad attribute. Geometry optimisation and energies are unchanged; only RMSG reporting was affected.
