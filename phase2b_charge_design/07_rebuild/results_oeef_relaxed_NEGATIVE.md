# CORRECTED RESULT: external point-charge fields do NOT catalyse the RELAXED Claisen barrier

## Correction of a premature claim
An earlier note claimed the capacitor catalyses (scan barrier 15.78 < 17.47 bare). THIS WAS WRONG --
the SCAN DERAILED before reaching the true TS (a known failure mode on this floppy surface: the
constrained scan drops into the product basin, reporting the derailment height 15.78 rather than the
true saddle). The NEB (trustworthy connected-path saddle) refutes it.

## NEB barriers (the trustworthy method) -- all FAR ABOVE bare
  design              field           NEB HEI barrier (converging)   vs bare 17.47
  capacitor (50 chg)  uniform 0.020   ~45 (69->53->47->45, iter 36)  +28  (RAISED)
  OEEF-020 (6 chg)    0.022           ~37-46 (oscillating ~iter 65)  +20-28 (RAISED)
  pair (2 chg)        0.011           was climbing >45 (killed)      +28  (RAISED)

## Conclusion
The strong FIXED-geometry OEEF effect (-10.6 kcal/mol at 0.02 a.u. along -z) does NOT survive geometric
relaxation. On the relaxed path, the external field distorts the substrate in a way that RAISES the
barrier by ~20-28 kcal/mol, REGARDLESS of field uniformity (capacitor highly uniform, still +28). The
"monotonic uniformity->barrier trend" seen in the SCANS was an artifact of the scans derailing at
different points -- not a real trend. External point-charge fields do not catalyse this relaxed
reaction.

## Methodological lesson (important)
On this floppy, field-embedded surface the relaxed SCAN is UNRELIABLE -- it derails before the TS and
underestimates the barrier. Only the NEB (connected-path saddle) gives the true barrier. All scan
"barriers" here (15.78, 21.32) are derailment artifacts and must NOT be trusted. Use NEB only.

## Consequence: pivot to the STERIC CRADLE (the plan's documented fork)
Since electrostatic (field) catalysis does not survive relaxation for this substrate, the path is
GEOMETRIC catalysis: hold the substrate in the near-attack conformation via a neutral steric cradle,
not a field. Burschowsky PNAS 2014 supports this (the enzyme uses a placed charge AND geometry). This
is fork (b) of DECISION_next_optimizer.md.
