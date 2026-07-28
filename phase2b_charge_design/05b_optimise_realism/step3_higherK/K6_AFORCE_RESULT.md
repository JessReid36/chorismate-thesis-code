# K6_aforce (DFT) RESULT: converged-DISSOCIATED (a_force ruptures the substrate at K6)

Converged (213.9 min) but the substrate FRAGMENTED under the K6/a_force field:
  C1-C6  = 5.612 A
  O3-C4  = 6.004 A   ether bond COMPLETELY BROKEN (reactant ~1.43)
  C4-C6  = 2.485 A
  connectivity: FRAGMENTED (only 8/24 atoms in the main fragment)  -> INVALID endpoint

## Reading
The a_force coordinate at K6 drives DISSOCIATION, not preorganisation. This confirms the mid-run
observation (O3-C4 climbing to ~5.7) and completes the K6 x {coordinate} comparison:
  K6_apot   : C1-C6 4.243, intact ether, SINGLE molecule -> VALID (but still open, outside window)
  K6_aforce : FRAGMENTED (O3-C4 6.0) -> INVALID (ruptured)
So at K6 the WINNING coordinate is a_pot (valid, if open); a_force ruptures. This mirrors the flip we
already recorded: K2 -> a_force better; K6 -> a_pot better; the winning coordinate is K-dependent and must
be determined empirically per K (PROTOCOL NOTE in 05c PLAN).

## Why a_force ruptures at K6 but was fine at K2
K2_aforce (net 0 dipole) relaxed cleanly to 4.06. K6_aforce packs more charges (incl. more negative sites);
the denser/more-negative field near the -2 dianion over-polarises and breaks the ether bond - the same
over-polarisation limit seen in K10_apot (fragmented) and K10_aforce (dissociating). Consistent picture:
beyond a threshold of charge density near the substrate, the point-charge field ruptures it.

## Status
Retained as a valid NEGATIVE result (the design over-polarises). Not usable as a cradle baseline (K6 cradle
uses K6_apot, the valid one). Relaxed geometry + log committed for the full record.
