# r2SCAN-3c calibration variant: cheap-geometry cross-check of the higher-K DFT jobs

## What this is
The SAME four higher-K designs (K6/K10 x a_pot/a_force) run at r2SCAN-3c, IN PARALLEL with the DFT
(B3LYP-D3BJ/def2-SVP) jobs, to test whether the cheaper method reproduces the DFT relaxed geometry on
this dianionic / pericyclic / field-sensitive system. This is the DESIGN_DECISIONS Tier-2a calibration
gate run explicitly (do NOT skip it - r2SCAN-3c has a documented SIE risk for the zwitterionic TS).

## Method (the ONLY difference from relax_higherK.py)
  DFT     : ! B3LYP D3BJ def2-SVP def2/J RIJCOSX CPCM TightSCF
  r2SCAN  : ! r2SCAN-3c CPCM TightSCF          (composite: def2-mTZVPP + D4 + gCP built in)
Everything else identical: coords, charges, per-pair LJ, elstat embedding, frozen surrogates,
QM charge -2, CPCM eps=4, LOOSE conv. So any geometry difference is a PURE method effect.
r2SCAN-3c is ~3-5x cheaper and DFT-quality for geometry; energies COARSE only (final barrier stays DFT).

## Run (own dir per job - scratch-collision discipline, per the K6/K10 DFT lesson)
For each name in {K6_apot,K6_aforce,K10_apot,K10_aforce}:
  mkdir -p run_<name>_r2scan
  cp relax_higherK_r2scan.py <name>_coords.xyz <name>_charges.txt run_<name>_r2scan.pbs run_<name>_r2scan/
  cd run_<name>_r2scan && qsub run_<name>_r2scan.pbs && cd ..

## Compare all 8 (calibration verdict)
python3 compare_methods.py .
Reports per-design DFT-vs-r2SCAN C1-C6, O3-C4, connectivity, and:
  - geometry agreement (C1-C6 within 0.2 A)
  - near-attack WINDOW verdict agreement (both in/out of 2.6-3.5 A)
  - connectivity agreement
GATE: geometry agrees within 0.2 A AND window verdict AND connectivity all match on all 4
  -> r2SCAN-3c usable as the cheap geometry driver for subsequent (higher-K / cradle) rounds.
  -> if it FAILS, keep DFT for geometry; r2SCAN-3c demoted or abandoned per DESIGN_DECISIONS.

## Files
relax_higherK_r2scan.py (method-only variant), run_<name>_r2scan.pbs (x4), compare_methods.py
Reuses the committed <name>_coords.xyz + <name>_charges.txt (same designs as the DFT jobs).
