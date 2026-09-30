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
- `s13_lj_charge_site_test.sh` — can a designed charge site carry a Pauli wall? Tests
  whether `!QMMM` runs without CPCM, whether a hand-edited `ORCAFF.prms` with a
  free-floating site is accepted, and where a substrate oxygen stops under relaxed
  optimisation, with a bare point charge as the control. The criterion is step 7's own:
  2.0–3.2 Å is physical salt-bridge range, below 2.0 Å the wall is too weak, beyond
  4.0 Å the charge cannot act.
- `s14_design_milp.py` — the Tier 1 optimiser. Certified global optimum of the
  linear-response design problem by mixed-integer programming, plus exact enumeration of
  the gauge-degenerate optimal set by no-good cuts. Solved by HiGHS through
  `scipy.optimize.milp`, so no commercial licence is needed: the mean-penalised-by-spread
  objective is LP-representable if spread is measured as a semi-deviation rather than a
  standard deviation (Konno & Yamazaki, Management Science 1991, 37, 519–531). Run
  `selftest` before trusting it; one of its checks is against a closed-form optimum.
- `s15_dv_matrix.py` — the multi-frame difference-potential matrix `s14` consumes.
  Reworked from `phase2b_charge_design/03_dvpot/s3_dv_on_grid_v2.py`. Uses the IN VACUO
  densities, not the QM/MM ones, because the design environment is bare substrate plus
  charges and the QM/MM densities already contain the protein field the design is meant
  to replace. Reads its frame list from the grid header so the two cannot disagree.
- `s15b_regen_density.sh.REJECTED` — a wrong turn, kept with its reasoning. Read it
  before attempting anything with ORCA density files.
- `s12_build_grid.py` — candidate-site grid on the union of atom-centred vdW
  spheres, per [DTHESIS] Ch.2 Eq.2.16. Pure function of its configuration; every
  parameter is written into the output header. `min_approach` acts as a FLOOR on the
  site radius and is set by `s11`; `min_approach=0` recovers [DTHESIS] exactly.

## Conventions

ORCA is not on PATH on this cluster. Invoke it by full path with the library
settings in `s8_invacuo_new.pbs`; jobs cannot run on a login node.

The login node runs **Python 3.6.8**. No `math.dist`, no walrus operator, no
f-string `=` specifier.

### ORCA 6 density files: orca_vpot takes a NAME, not a path

`orca_vpot`'s second argument is the density's name INSIDE the `.densities` container,
not a file on disk. ORCA 6 no longer writes a standalone `.scfp`, and `KeepDens` does
not change that. The manual says so only obliquely, via the note that a mismatched
container basename must be passed as a FIFTH argument.

Verified 2026-09-30:

        orca_vpot sp_24883_R.gbw sp_24883_R.scfp points.xyz out.txt
        -> "Electrostatic potential evaluated in 0.428 sec"

Four routes were tried before the manual was read closely enough: passing `.densities`
directly, `KeepDens` on a fresh single point, `MORead`/`NoIter` regeneration, and
`orca_plot` extraction. None was checked against a one-line test first; the test took
two minutes.

### Three staleness bugs, all of which produced plausible numbers

Recorded because they share a shape and will recur.

1. **A hardcoded frame list.** `s8_invacuo_new.pbs` computed in vacuo references for
   five frames while the ensemble had grown to thirty. `s15_dv_matrix.py` had the same
   defect and now reads its frame list from the grid header instead.
2. **A silent skip.** `s12_build_grid.py` was asked for thirty frames, found geometries
   for twenty-two, enclosed those and reported "22 frames" without complaint. It now
   names what it could not find and records the gap in the grid header.
3. **A stale commit.** The committed `s11_probe_distance_scan.py` lacked the
   charge-magnitude scan that the running copy had, because the repository copy came
   from an earlier tarball than the one on the cluster.

The common lesson: anything that FILTERS must report what it filtered out. A count that
looks reasonable is not evidence that nothing was dropped.

### ORCA job cost on this cluster — measured, after a parsing error

**First, the error that caused three wrong conclusions.** ORCA prints

        TOTAL RUN TIME: 0 days 0 hours 4 minutes 12 seconds 340 msec

so the fields are days, hours, MINUTES, seconds. A parser reading `$4*3600+$6*60+$8`
returns days*3600 + hours*60 + minutes and drops the seconds entirely. "4.0 s" was
really 4 minutes. Every timing conclusion drawn before 2026-09-30 came from that.

        correct:  awk '{t=$4*86400+$6*3600+$8*60+$10}'

Sanity-check any parsed timing against something independent — `qstat -f` reports
`resources_used.cput` and `walltime` — before building on it.

**What actually drives cost: basis and grid, not I/O.** Measured on 24-atom single
points of the chorismate dianion, one core:

        hf3c      14 s    MINIX basis, no DFT grid, 12 SCF cycles
        b97-3c   202 s    mTZVP plus grid integration, 20 SCF cycles

A 14x time difference on 1.7x the cycles, so it is cost PER CYCLE. Not a convergence
pathology, despite the substrate being an electronically unbound dianion.

**Parallelism must be chosen per method, by measured cost.**

        b97-3c   202 s serial  ->  71-89 s on 8 ranks    about 2.5x
        profile: 81 % SCF iterations, 12 % startup

That is the profile where `%pal` pays. A 14 s job is startup-dominated and `%pal`
costs more than it saves. `s10` therefore sets ranks per method: xtb2 and hf3c
serial, b97-3c / r2scan-3c / ref on 8. `s11` is entirely at the production level and
runs parallel throughout. A blanket setting either way is wrong.

**Node-local scratch is real but secondary.** The same input, measured with `time`:

        shared filesystem   real 28.7 s   user 4.9 s
        node-local scratch  real  6.6 s   user 3.6 s

Worth keeping — both scripts run each job in `${TMPDIR:-/tmp}` and copy back only the
`.out` — but it is a smaller effect than basis and grid, and the earlier claim that it
was *the* cause was wrong.
