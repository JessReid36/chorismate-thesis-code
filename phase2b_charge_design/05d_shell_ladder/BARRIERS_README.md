> **SUPERSEDED — do not use this protocol.** The method below computes the barrier as a difference of two
> *independently* field-relaxed endpoint energies (E(TS under field) − E(reactant under field)). That
> endpoint-difference barrier is **invalid** under a fixed external field: it is dominated by the
> substrate's displacement through the fixed potential and does not cancel, giving unreliable/impossible
> values. The retained barrier method is the **NEB-under-field on a connected path** between screened
> endpoints — see `../SHELL_LADDER_NEB_RESULT.md` and `../SHELL_LADDER_NEB_CORRECTION.md`. Any
> endpoint-difference number produced by this script (including the "+35.51 at 8A" quoted below) is not
> reportable and is retained only for historical context. [archived]

---

# ALL-layer barriers: 2-15A coarse + 7.0-9.0 fine (0.1A). Exhaustive discrete-charge barrier scan.

## What
Computes the REAL barrier dE! = E(TS under field) - E(reactant under field) for EVERY shell design,
including fragmented and open ones (user: all layers, even throwaways). Makes the discrete negative
exhaustive: "we computed the actual barrier for all ~32 designs 2-15A" not "we inferred the rest".

## Requires (all on HPC, in 05d_shell_ladder/relax)
- each shell's design: shell<tag>_coords.xyz + shell<tag>_charges.txt  (tag = 2p0..15p0, 7p1..8p9)
- completed reactant ladder job: run_shell<tag>/ladder_<tag>_relax.log  (for E_react)
- ts.xyz (copied by the script from 01_geometry)
- relax_barrier.py, make_ts_any.py

## Run
1. put make_ts_any.py, relax_barrier.py, submit_all_barriers.sh, collect_barriers.sh in relax/
2. bash submit_all_barriers.sh     # builds TS-under-field design + PBS per shell, submits each (own dir)
   (skips any shell whose design coords aren't present; ~32 jobs, run in parallel)
3. when done: bash collect_barriers.sh   # prints the full barrier table

## Interpretation
- barrier < 17.47 = CATALYTIC (would be a real hit - investigate)
- barrier >= 17.47 = anti/non-catalytic
- FRAG note = reactant was fragmented; barrier is broken-vs-broken, NOT physically meaningful (report as
  "fragmented, no valid barrier" - the geometry already disqualifies it; the number is for completeness only)

## Expectation
8A already gave +35.51 (anti-catalytic). Fragmented layers (3-6A) have no meaningful barrier. Open layers
(7-15A) expected >= bare. If ALL >= 17.47 (or fragmented) -> discrete +-1 fails on the BARRIER metric across
the entire 2-15A range -> the certified negative that forces the fractional QP.
