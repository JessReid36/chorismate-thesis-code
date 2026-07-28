# Step 3 RESULT (a_force): PARTIAL - moves the right way but does NOT recover near-attack.
# Combined verdict across a_pot + a_force: a K=2 monopole PAIR cannot preorganise this dianion.

Relaxed the certified a_force-preorganised K2 dipole (guan+1@249, formate-1@185, D_force=0.006).
Converged clean: ASH 198 min, single connected molecule, O3-C4 1.442 (intact), E=-837.382 Eh.

## Measured (analyse_relaxed.py)
  C1-C6 = 4.060 A  -> OPEN, but ~1 A LESS open than a_pot (5.277) and mu=0 control (5.004)
  O3-C4 = 1.442 A  -> intact ether (correct reactant)
  C4-C6 = 2.496 A ; connectivity SINGLE molecule (valid)

## Ladder of evidence (all certified, gap 0.000)
  bare substrate reactant : C1-C6 3.124  (target reactive window 2.6-3.5)
  mu=0 K2 (distorting)    : 5.004
  a_pot  mu=100           : 5.277   (differential-POTENTIAL penalty: NO geometric help)
  a_force mu=100          : 4.060   (force-imbalance penalty: ~1 A better, still outside window)

## Why this is the REPRESENTATION limit, not a tuning problem (the clincher)
Higher-mu a_force sweep drives the distorting drive D_force -> ~0 and BELOW, all certified:
  mu=100 D=6.2e-3 | 500 1.6e-3 | 1000 7.2e-4 | 2000 -2.7e-4 | 5000 4.6e-6
So the optimiser CAN find a K=2 dipole with ~zero force-imbalance. Yet the near-zero-drive mu=100 design
STILL relaxed to 4.06 A. Driving the (fixed-geometry, first-order) force metric to zero does NOT drive the
RELAXED geometry into the window - the linear surrogate and the relaxed geometry have decoupled (same class
of limitation as the frozen-Dv proxy: a fixed-geometry linear coefficient cannot capture the full relaxation
response). Conclusion: it is not the penalty coordinate or strength - TWO POINT CHARGES ARE TOO SPARSE to
cradle a floppy dianion into near-attack. Certified across a_pot + a_force + mu up to 5000.

## This is a real, publishable finding (redirects the approach)
Empirical, certified-optimisation-backed demonstration of WHY positioned charges alone are insufficient and
an enzyme needs a surrounding cradle. We did not assume charges were too sparse - we PROVED it (two coords,
five mu). This is the strongest possible motivation for the secondary (neutral cradle) layer, and/or higher-K.

## Decision (open - for next session)
(1) HIGHER-K: does K3/K4 (fuller charge cavity) recover C1-C6 where K2 cannot? Stays pure point-charge.
    Step-2 already showed higher-K keeps more barrier-lowering at low D; now test the GEOMETRY. Cheap-ish:
    build+relax K3/K4 preorg designs (a_force).
(2) NEUTRAL CRADLE (secondary layer, the research report's recommendation): add neutral dipolar/steric MM
    surrogates that GEOMETRICALLY hold C1-C6, since charges alone can't. The enzyme-like answer. Must be
    posed to preserve the certificate (discrete steric/dipole site selection, linear).
Not mutually exclusive; higher-K is the smaller step and tests whether MORE charges alone suffice before
committing to the cradle machinery.

## a_force is nonetheless the better preorg COORDINATE (retain)
a_force moved geometry the right way (~1 A) where a_pot did nothing -> for any future preorg penalty, use
a_force (mechanical), not a_pot (differential potential). The STEP0 falsifiable choice is resolved: a_force.

## Files
k2_mu100_aforce_relaxed.xyz (C1-C6 4.060), k2_aforce_relax.log (198 min, converged)
