# Result: catalytic ±1 layer re-fragments the RING (cradle holds the carboxylate)

## Verdict
Adding a certified K=4 catalytic layer (±1 charges, gap = 0.000) ON TOP of the validated cradle
(two +1 counter-charges) re-fragments the substrate. Both endpoints DISTORTED. But the fragmentation
MOVED, and where it moved is the finding.

## The diagnostic shift in failure location
- Bare dianion (every K 6-20): bond (12,21) breaks -- the carboxylate-to-RING bond (handle torn off).
- Cradle + K4: (12,21) INTACT -- the cradle held the carboxylate. Instead:
  - reactant: broke (3,4),(17,18); formed (3,23),(15,18)  [E=-836.996 Eh]
  - product : broke (15,17); formed (8,15)                 [E=-836.209 Eh]
  RING / reacting-region bonds, not the carboxylate handle.

## Interpretation
The cradle held the carboxylates (12,21 intact), solving the bare-dianion failure. But the catalytic
±1 charges in the reacting region at ~3.6 A are themselves strong enough to tear bonds THERE. The
fragmentation relocated from carboxylate (fixed) to ring (where the catalytic charges sit).

## What this establishes
±1 point charges at LJ-surface distance are TOO STRONG for the catalytic layer: they catalyse by
brute force that also breaks local bonds, independent of carboxylate neutralisation. The cradle is
NECESSARY (fixed the carboxylate failure) but NOT SUFFICIENT -- the catalytic charge MAGNITUDE is the
remaining problem. ±1 is GOCAT's maximum, not a requirement.

## Next (decided separately)
Reduce catalytic charge magnitude (sweep QMAX_catalytic down) or optimal fractional charges. Open
question: is there a magnitude that stabilises the TS without fragmenting the ring, or does any
barrier-lowering point-charge field fragment? Cheaply testable.

## Trustworthiness
LJ-min standoff (all charges on surface, none inside wall), L1 oracle gated to L0, LJ wall active on
QM, all-bond integrity. Catalytic certified gap 0, net-neutral; cradle validated alone.

## Provenance
- Design: L3_solver/finalise_K_cradle.py -> design_cradleK4_* (2 counter + 4 catalytic, gap 0)
- Screens: L1_oracle/run_cradleK4_{reactant,product}/cradleK4_*_result.txt (both DISTORTED, ring)
