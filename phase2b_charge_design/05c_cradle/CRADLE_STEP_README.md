# 05c cradle: first steric-pocket test (K2/a_force + K6/a_pot baselines)

## What this is
Adds a NEUTRAL 6-bead steric cradle (q=0 LJ walls, per-bead sigma) around the reacting C1-C6 region of a
fixed charge design, testing whether pure geometric confinement pulls the field-opened substrate into the
near-attack window [2.6-3.5 A] where CHARGES ALONE could not (K2_aforce 4.06, K6_apot 4.24).

## Design (adaptive placement, all DOFs, documented walls)
6-bead pocket covering 3 opening DOFs: 2 axial (beyond C1/C6), 2 out-of-plane, 2 in-plane lateral.
Beads grown OUTWARD until clear of substrate AND charge surrogates (margin 2.3 A); per-bead sigma set so
each is a real wall (guard/sigma ~1.1), capped 5.0. Neutral (q=0), frozen, carbon-like eps.
  K2/a_force: FULL 6 beads (2 charges leave the pocket open) - clean cage.
  K6/a_pot  : 5 beads (1 lateral dropped - the 6 charges crowd the pocket; still all 3 DOFs covered).
  => more charges crowd the cradle: an early K-dependence signal.
Steric-verified: bead-substrate 2.31, bead-charge >=2.32, bead-bead >=4.44. Start C1-C6 = 3.124 (compact).

## The test
Relax from the COMPACT reactant with (charges + cage). Does the substrate STAY compact (cage holds it in
the window) or does the field still pry it open past the walls?
  relaxed C1-C6 in [2.6-3.5] + single molecule -> STERIC CRADLE WORKS (pure geometry, zero added charge).
  relaxed C1-C6 still open -> steric confinement alone insufficient -> need dipolar cradle (a)/combos.

## Run (own dir per job - scratch discipline)
For each name in {K2aforce, K6apot}:
  mkdir -p run_<name> && cp relax_cradle.py <name>_cradle_*.txt <name>_cradle_coords.xyz run_<name>_cradle.pbs run_<name>/
  cd run_<name> && qsub run_<name>_cradle.pbs && cd ..
Dry-run first (on a compute node, NOT login): python3 relax_cradle.py <name>
Analyse: python3 <path>/analyse_relaxed.py run_<name>/<name>_cradle_relaxed.xyz
  (analyse_relaxed reads the first 24 atoms = substrate; beads/charges ignored by the C1-C6/connectivity check)

## Files
01_placement/place_cradle_beads.py (adaptive placement), <name>_beads.xyz (the cages)
02_relax/relax_cradle.py (runner: per-bead sigma, Ar beads q=0), <name>_cradle_{coords.xyz,charges.txt,beadsigma.txt}, run_<name>_cradle.pbs
