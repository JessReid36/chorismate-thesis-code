# Stage 0 — bare-field endpoint relaxation: observed failure modes (PRELIMINARY)

Relax reactant and product of the four max-lowering designs (K1–K4) under their BARE,
unguarded external point-charge fields, to observe whether the certified fields distort or
destroy the geometry before committing to the Tier-2 batch. 8 jobs = 4 designs × {reactant,
product}; native ORCA 6.0.1 %pointcharges + CPCM(eps=4) + Opt (coordsys cartesian); NO
geometric guard by design. Source: runs/<design>_<endpoint>/job.out (snapshot at cycles
~27–45; jobs left running to accumulate data — K1/K2 trajMin and final values to be appended).

## Headline
Nuclear implosion (register A1) is CONFIRMED for the strong stacked-charge designs K3 and K4,
at BOTH endpoints: a nucleus is pulled to well under 1 A of an external charge during
relaxation, from a >=3.5 A start. SCF converges throughout (scf_fail=0 on all 8) — the failure
is GEOMETRIC (nuclear implosion) + optimiser non-convergence, NOT electronic SCF breakdown.
The positive HOMO on K3/K4 is a consequence of the collapsed geometry, not independent spill-out.

## 1. Implosion CONFIRMED (K3, K4) — trajMin = min atom-charge distance reached along the opt
  K4_product  trajMin 0.024 A (atom essentially ON the +1 site); O3-C4 opened to 5.25 A  [runs/K4_product/job.out]
  K4_reactant trajMin 0.872 A; O3-C4 3.93 A                                              [runs/K4_reactant/job.out]
  K3_product  trajMin 0.910 A; O3-C4 opened to 5.71 A                                    [runs/K3_product/job.out]
  K3_reactant trajMin 0.931 A; O3-C4 1.52 A                                              [runs/K3_reactant/job.out]
Product endpoints implode hardest (O3 unanchored once O3-C4 breaks), as predicted. The +1 at
3.84 A from O3 (the Arg90 valley) is the attractor; floor breached by ~1.1-2.0 A.

## 2. Two distinct implosion dynamics
  K4_product DIVERGES: MAX gradient 1.07 -> 1.76 -> 3.26 -> 7.39 Eh/Bohr (~doubling/cycle)
    — a blow-up, not slow convergence.                                                  [runs/K4_product/job.out]
  K3 (both) and K4_reactant IMPLODE-then-SETTLE: MAX gradient SHRINKS (~0.004-0.011 Eh/Bohr)
    walking into a collapsed local minimum — dangerous: may FALSELY report convergence on a
    garbage geometry.                                      [runs/K3_*/job.out, runs/K4_reactant/job.out]

## 3. Non-convergence, no implosion at snapshot (K1, K2)
Plateau above the 3e-4 Eh/Bohr MAX-gradient tolerance, not tightening:
  K1: MAX grad ~0.003-0.018 over ~27-28 cycles          [runs/K1_reactant/job.out, runs/K1_product/job.out]
  K2: MAX grad up to ~0.12 (K2_reactant) over ~26-27 cycles   [runs/K2_reactant/job.out, runs/K2_product/job.out]
Pure `coordsys cartesian` (chosen for frame-lock) on a field-roughened surface plateaus. Root
cause: a fixed external field exerts a net rigid-body force/torque that a curvature-poor
Cartesian optimiser tightens only slowly. K1/K2 trajMin >=3.25 A at last read (no implosion).

## Implications for Stage 1
1. A guard/handling for implosion is MANDATORY: certified bare +/-1 charges DO destroy the
   substrate for K3/K4 under relaxation.
2. `coordsys cartesian` is correct (it must relax the field-induced rigid-body force) but needs
   curvature help (computed Hessian) to converge — redundant internals are contraindicated in a
   fixed field (they cannot relax the overall force/torque).
3. Ranking taxonomy (register G1): K4_product-type divergence and K3/K4-type implosion are
   RETAINED as ranking failures (penalty rank); K1/K2-type non-convergence is a tooling issue to
   fix, not a design failure.
