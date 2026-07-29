# Negative controls for the shell ladder (field-off / body-only baselines)

## Why (the gap this fills)
Every ladder design has BOTH a field (charges) AND bodies (guanidinium/formate LJ walls). When a substrate
opens/fragments we can't attribute it to field vs bodies. And "bare 3.124" is the EXTRACTED grid geometry,
not a relaxation under our production protocol. Two controls fix this:

## Control A - BARE substrate, no bodies, production protocol (relax_bare.py)
24-atom dianion alone, pure QM, same B3LYP-D3BJ/def2-SVP/CPCM + loose conv, NO frozen atoms, NO bodies,
NO field. Tells us where the substrate settles with ZERO external influence. This is the TRUE baseline
(replaces the extracted 3.124 as the reference the ladder is measured against).
  Run: (dry) python3 relax_bare.py    (real) qsub run_bare_control.pbs
  needs reactant.xyz in the dir.

## Control B - closest shell (2A) bodies present, q=0 (relax_shell.py + zeroed charges)
Same 2A surrogate BODIES in the same positions, all charges 0 -> LJ walls act, field is OFF. Differs from
the real shell2p0 design by EXACTLY one variable (the charges). Isolates the steric-body effect at closest
(worst-case) approach.
  Build the zero-charge file:  python3 make_control.py shell2p0_charges.txt shell2p0ctrl_charges.txt
  and copy coords:             cp shell2p0_coords.xyz shell2p0ctrl_coords.xyz
  Run: python3 relax_shell.py shell2p0ctrl   (dry) ; qsub run_ctrl_shell2.pbs  (real)

## Reading the controls
  A (bare)      -> baseline C1-C6 (call it B0). If B0 >> 3.124, the substrate is intrinsically floppy.
  B (2A q=0)    -> if ~B0, close bodies are sterically INERT -> all charged-design opening is FIELD.
                   if != B0, bodies perturb -> must subtract body effect from charged results.
  Then each charged shell's deviation from A = TOTAL effect; from B = FIELD-only effect.

## Note
Only 2 controls needed (not 14): A = baseline, B = body-effect at CLOSEST shell (worst case). If close
bodies are inert, far bodies certainly are, so the whole ladder is field-attributable. Add a far-shell q=0
control only if B shows the bodies do perturb.
