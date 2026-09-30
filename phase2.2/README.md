# phase2.2 — charge design, rebuilt

Scripts for Phase 2.2. The plan, the rationale and the literature basis are in
`chorismate-thesis-results/phase2.2/PHASE22_PLAN.txt`; read that first.

## Why 2.2 rather than continuing 2b

Phase 1 was audited step by step (see the `consolidated_works` repository). Two
findings force a rebuild rather than a patch:

- The committed `dv_grid_ensemble_mean.tsv` states in its own header that it was
  built from ten frames, all of them below the 20,000 ps equilibration cut and all
  run at the 400 kJ spring. Phase 1's finalised position is that this pilot set
  must not be used for reported results.
- Point charges plus CPCM eps = 4 matches neither published GOCAT setting, and the
  continuum supplies the stabilisation the charges are being optimised to provide.

## The three Hartke-group sources are NOT interchangeable

    [D2018]   Dittner & Hartke, JCTC 2018, 10.1021/acs.jctc.8b00151
              abstract point charges, FIXED path, PM7, Menshutkin SN2
    [DTHESIS] Dittner PhD thesis (gocat_2.pdf)
              Ch.2 surface algebra; Ch.3.6 Coulomb implosion and the static vs
              ADAPTIVE fitness distinction; Ch.7 adaptive scheme on a Diels-Alder
    [B2021]   Behrens & Hartke, Top. Catal. 2021, 10.1007/s11244-021-01486-1
              REAL molecular fragments, MEP recomputed each fitness call,
              GFN-xTB/OPLS-AA, SN2 in acetone and KSI from 1OH0 with waters removed

Mixing them has already caused errors in this project. Name the source and its
setting whenever one is cited.

## Scripts

- `s10_cheap_method_benchmark.py` — is a cheap method usable inside the design
  loop for THIS reaction? Tests XTB2, HF-3c, B97-3c and r2SCAN-3c against the
  committed `barrier_vac` column, plus a field-response test using a +1 charge at
  Burschowsky's 3.2 A from the ether oxygen. Verdict rules fixed before running.
- `s11_smearing_calibration.py` — sets the Gaussian width of the smeared charge so
  that one +1 at 3.2 A reproduces Burschowsky's measured 5.9 kcal/mol of
  differential transition-state stabilisation.

## Conventions

ORCA is not on PATH on this cluster. Invoke it by full path with the library
settings in `s8_invacuo_new.pbs`; jobs cannot run on a login node.
