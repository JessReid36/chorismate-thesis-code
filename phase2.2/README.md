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
- `s11_probe_distance_scan.py` — where does the bare point-charge model break for
  this substrate, and what does one +1 at Burschowsky's 3.2 A actually buy? A
  distance scan with three diagnostics whose expected behaviour is stated in
  advance: log-log slope of the differential stabilisation (dipole-like coupling
  predicts about -2), Loewdin charge on the nearest atom, and the reactant HOMO.
  Calibrates against the DIFFERENTIAL, 5.9 - 0.6 = 5.3 kcal/mol, not 5.9.
  Replaces the smearing-width calibration that was planned: a width cannot be chosen
  before the failure is measured, and [DTHESIS] Ch.3.6 attributes the failure to
  missing REPULSION rather than to the singularity alone, so a Gaussian may treat a
  symptom.
- `s12_build_grid.py` — candidate-site grid on the union of atom-centred vdW
  spheres, per [DTHESIS] Ch.2 Eq.2.16. Pure function of its configuration; every
  parameter is written into the output header. `min_approach` acts as a FLOOR on the
  site radius and is set by `s11`; `min_approach=0` recovers [DTHESIS] exactly.

## Conventions

ORCA is not on PATH on this cluster. Invoke it by full path with the library
settings in `s8_invacuo_new.pbs`; jobs cannot run on a login node.

The login node runs **Python 3.6.8**. No `math.dist`, no walrus operator, no
f-string `=` specifier.

**Do not use `%pal` for small single points.** Measured 2026-09-30 on 24-atom jobs:
ORCA reported 2-4 s of compute while consecutive outputs appeared ~180 s apart, i.e.
98 per cent MPI startup and teardown. Serial is roughly 40x faster in wall time for
jobs this size.
