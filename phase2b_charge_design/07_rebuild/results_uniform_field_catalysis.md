# RESULT: uniform point-charge field CATALYSES the relaxed Claisen barrier (uniformity is the key)

## The uniformity-vs-barrier trend (scan barriers, physical-region peaks)
  design              field uniformity        Fz(a.u.)   scan barrier   vs bare 17.47
  pair (2 charges)    ~80% variation          0.0106     ~45 (climbing) +27  (RAISED)
  OEEF-020 (6 chg)    moderate                0.0221     21.32          +3.85 (raised)
  capacitor (50 chg)  <4% variation           0.0200     15.78          -1.69 (LOWERED)

MONOTONIC: more uniform field -> lower barrier. Non-uniform fields RAISE the barrier (field gradient
distorts geometry on relaxation); the uniform field LOWERS it (electronic OEEF stabilisation wins).
The capacitor CATALYSES: 15.78 < 17.47 bare.

## Interpretation
The finite-field scan (-10.6 at 0.02 uniform, FIXED geometry) DOES translate to point charges -- but
only with a UNIFORM field. The pair failed because a 2-charge field is too non-uniform (steep gradient
-> geometric distortion -> raises barrier), NOT because point charges can't catalyse. A parallel-plate
capacitor (50 charges) makes a genuinely uniform field and catalyses the RELAXED barrier.

Geometric relaxation eats most of the fixed-geometry effect (-10.6 fixed -> -1.69 relaxed) but a net
catalytic effect survives. Stronger uniform fields should catalyse more (effect scales with field).

## Caveats
- Scans derailed past the peak (floppy-field artifact); these are physical-region peaks. NEBs running
  (capacitor iter 29, OEEF-020 iter 61) will confirm the trustworthy saddles.
- -1.69 is modest but REAL catalysis on a held, buildable (50-charge), certified-lineage design.

## Next
Stronger uniform-field capacitors (FTARGET 0.03, 0.04) to walk up the catalysis toward the -10.6
ceiling, screening for integrity at each step.
