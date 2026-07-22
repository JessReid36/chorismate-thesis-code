# Stage 0 - bare-field endpoint relaxation: observed failure modes (FINAL)

Relax reactant and product of the four max-lowering designs (K1-K4) under their BARE, unguarded
external point-charge fields, to test whether the certified fields distort or destroy the
substrate before the Tier-2 batch. 8 jobs = 4 designs x {reactant, product}; native ORCA 6.0.1
%pointcharges + CPCM(eps=4) + Opt (coordsys cartesian); NO geometric guard by design. All jobs
ran to the 72-cycle cap. Raw outputs + per-design results table: results repo
phase2b_charge_design/07_stage0/ (runs/<tag>/job.out, STAGE0_results_table.md).

## Headline
Bare native %pointcharges is NOT viable for endpoint relaxation. 0/8 converged; 6/8 imploded.
Nuclear implosion (register A1) is CONFIRMED for every design carrying the +1 at 3.84 A from the
ether O3 - K2, K3, K4 at BOTH endpoints - a nucleus pulled to 0.01-1.06 A of an external charge
from a >=3.5 A start. SCF converges throughout (scf_fail=0): the failure is GEOMETRIC (nuclear
implosion) + optimiser non-convergence, NOT electronic SCF breakdown.

## 1. Implosion CONFIRMED (K2, K3, K4 - both endpoints)
Min atom-charge distance at the final geometry (bare premise: >=2 A):
  K2_reactant 0.98 A (O3->+1 1.72; O3-C4 1.54)      K2_product 1.05 A (O3->+1 1.05; O3-C4 torn to 5.70)
  K3_reactant 0.97 A (O3->+1 1.81; O3-C4 1.51)      K3_product 1.06 A (O3->+1 1.06; O3-C4 torn to 5.92)
  K4_reactant 1.03 A (multi-charge collapse)         K4_product 0.01 A DIVERGE (MAX grad ~4767; O3->+1 0.75)
Product endpoints implode hardest (O3 unanchored once O3-C4 breaks). K2, which did NOT implode at
the mid-run snapshot, imploded by the 72-cycle cap - i.e. implosion is progressive, not immediate.

## 2. K1 - intact but non-converged (the control)
K1 (lone -1, 8.4 A from O3; no +1) is the ONLY geometrically intact design: no implosion
(min atom-charge 7.6-8.2 A, O3-C4 intact), but it STILL did not converge - MAX grad plateaus at
~0.0022 vs the 3e-4 tol. So even without implosion, bare `coordsys cartesian` cannot tighten the
gradient in the fixed field.

## 3. Root causes
- Implosion (K2/K3/K4): real, unopposed Coulomb attraction of the electron-rich O3 (and, in K4,
  other atoms) to a bare +/-1 point charge; no Pauli wall on a bare charge to stop a nucleus.
- Non-convergence (all 8): a fixed external field exerts a net rigid-body force/torque; a
  curvature-poor Cartesian optimiser tightens it only slowly. Redundant internals are
  contraindicated (translation/rotation-invariant -> blind to that force).

## Implications for Stage 1 (uniform, all designs)
1. Convergence: keep coordsys cartesian + add a computed Hessian (Calc_Hess true, Recalc_Hess 5)
   to tighten the field-induced gradient. Tests on the K1-type plateau.
2. Implosion handling: detect via min atom-charge distance (<2 A) and RETAIN as a ranking failure
   (register G1); do not fabricate a floored barrier on a collapsed geometry.
3. Soft one-sided harmonic floor is NOT available in %geom Opt (rigid constraints only; one-sided
   walls exist only in the %md module, and whether they are honoured by its L-BFGS Minimize is
   unverified) - hence detect-and-fail rather than a native floor.
