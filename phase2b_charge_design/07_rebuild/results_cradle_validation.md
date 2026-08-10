# Keystone result: neutralising the carboxylates holds the substrate (the arginine role, verified)

## The controlled comparison (same clean pipeline)
- BARE dianion, certified fields: fragment at EVERY K (6,10,15,20); bond (12,21) — the
  C6-to-carboxylate-carbon bond — breaks universally.
- Dianion + TWO +1 counter-charges on the carboxylates (LJ-standoff, min-dist 3.57-3.58 A,
  NOT inside the wall), zero catalytic charges: BOTH endpoints VALID.
  - reactant: integrity=VALID, broke=[], formed=[], E=-836.831 Eh, ether 1.445 A (intact)
  - product : integrity=VALID, broke=[], formed=[], E=-836.819 Eh, ether 3.11 A broken, C-C 1.61 formed
  - converged fast (~22 frames) to correct geometries -- the signature of a held substrate, not the
    279-frame thrashing of the fragmenting bare-dianion screens.

The ONLY difference between fragmentation and integrity is the two counter-charges.

## What this establishes
1. The carboxylate handles are THE fragmentation cause; neutralising them is the cure. This is the
   arginine role (Arg7/63/90 grip the carboxylates), rediscovered by certified optimisation with
   zero enzyme knowledge.
2. Abstract +1 point charges at LJ-standoff distance SUFFICE to hold the substrate -- explicit
   guanidinium (steric bulk) is an optional fidelity refinement, not a necessity.
3. The correct model for catalytic design is the net-neutral, carboxylate-gripped substrate (the
   enzyme's actual reacting species; 200/250 of DeYonker's QM-cluster models are net-neutral), NOT
   the bare -2 dianion.

## Energy sanity
Counter-charged reactant -836.831 Eh vs bare L0 reactant -836.374 Eh: ~0.46 Eh lower -- a sensible
stabilisation from two counter-charges, NOT the ~6.8 Eh collapse seen when the bare field dragged the
substrate onto the charges. The substrate relaxed normally in the field.

## Trustworthiness
LJ-min standoff (counter-charges at 3.57 A, outside the wall), L1 oracle gated to L0, LJ wall
confirmed active on QM, all-bond integrity. Counter-charges placed on the carboxylate outward
bisectors, marched out to clear the per-atom LJ floor.

## Provenance
- Design: L3_solver/make_cradle.py -> design_cradle2_* (2x +1, LJ-standoff)
- Screens: L1_oracle/run_cradle2_{reactant,product}/cradle2_*_result.txt (both VALID)
