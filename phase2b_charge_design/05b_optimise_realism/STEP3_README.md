# Step 3 — relaxed validation of K2 mu=100 preorganised design (the arbiter)

PLAN Step 3. Tests whether the Step-2 preorganisation penalty actually keeps C1-C6 compact under
RELAXATION (not just on the frozen proxy). This is the make-or-break check for the whole 05b approach.

## The controlled comparison
- CONTROL (mu=0, distorting): K2 sites 150/249 (+1/-1), D=+0.105. ALREADY relaxed (production_batch,
  this session) -> C1-C6 = 5.00 A (distorted open, away from near-attack).
- TEST (mu=100, preorganised): sites 105/121 (TWO +1, net surrogate +2), D=+0.0006 (~0 distorting drive).
  Same K=2, same substrate. Certified gap 0.000. Steric realizability CONFIRMED (guan-guan 2.85 A,
  guan-substrate 2.20-2.38 A, no clash).

## Decision rule (STEP0 falsifiable arbiter)
- relaxed C1-C6 in [2.6, 3.5] A  -> preorg penalty WORKS; a_pot was the right mode. HEADLINE RESULT:
  a certified, preorganisation-aware design that stays reactive where the un-penalised optimum distorted.
- relaxed C1-C6 stays open (>~4 A) -> a_pot was the WRONG mode; revisit with a_force (force-imbalance).

## Method
Validated MM-molecular-surrogate production method (relaxation_attempts/ash_guarded/k2_mmsurr_spike.py),
adapted for TWO guanidiniums (no formate): 44 atoms (24 QM -2 + 2x10 MM guanidinium +1 each). Explicit
QM-MM LJ (480 pairs), frozen surrogates, HDLC, LOOSE conv. write_xyzfile (not the write_xyz bug).

## Files
k2_mu100_coords.xyz  : 44-atom geometry (guanidinium central C at each site, verified delta 0.0000)
k2_mu100_mmsurr.py   : ASH runner (dry-run: no arg; real: 'run')
run_k2_mu100.pbs     : PBS wrapper
build_coords.py note : coords built by translating validated guanidinium template to sites 105/121

## Caveat carried from Step 2
mu=100 K2 is TWO +1 (no dipole) - an unusual motif, flagged not-trusted. Step 3 is exactly what decides
whether it is a real preorganised catalyst or proxy-gaming. Interpret C1-C6 AND check the substrate stays
a single connected molecule (no dissociation) as in ENDPOINT_VALIDITY_ASSESSMENT.
