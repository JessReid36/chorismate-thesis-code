# 05c_cradle — neutral secondary-layer (cradle) extension: PLAN

## Why this phase exists (the earned mandate)
05b proved, certified across K2/K6/K10 x {a_pot, a_force}, that POINT CHARGES ALONE cannot cradle the
chorismate dianion into the near-attack window (C1-C6 2.6-3.5 A):
  K2_apot 5.28 | K2_aforce 4.06 | K6_apot 4.24 (all valid, intact, OPEN) ; K10 over-polarises/FRAGMENTS.
Two limits bracket the pure-charge design space: too few charges -> geometry stays open (~4-5 A); too many
-> over-polarisation ruptures the substrate. The geometric preorganisation job needs NEUTRAL STRUCTURE
(dipolar cradle + steric walls), not more monopoles. This phase adds that secondary layer.
Source: STEP3_HIGHERK_RESULT.md; enzyme-cradle research report (secondary-layer options a-e).

## What the cradle IS (and is not)
IS: a layer of NEUTRAL surrogates - dipolar (amide/hydroxyl/water mimics) and/or steric (LJ-only walls) -
    placed AROUND the substrate TOGETHER WITH a fixed chosen K-charge design, to hold C1-C6 compact while
    the K-charges do the electrostatic TS-stabilisation work.
IS NOT: more charges (proven insufficient); NOT a replacement for the K-charges (they stay fixed); NOT a
    reconstruction of chorismate mutase (enzyme is a yardstick, not a template).

## Baseline design to wrap (DECISION)
Wrap K2 (the simplest, cleanest control): fixed guanidinium(+1)@249 + formate(-1)@185 (the a_force K2
design, relaxed C1-C6 4.06) OR the max-lowering K2 (150/249). Pick ONE baseline so the cradle's effect is
a clean before/after. K2 chosen over K6/K10 because: simplest, K6 no better geometrically, K10 fragments.
[Baseline charge design to be finalised in 01_placement - default: a_force K2, the best-behaved pair.]

## Secondary-layer options to test (research report a-e), in priority order
(b) STERIC LJ-only wall/cap: neutral LJ sites (q=0) that geometrically forbid C1-C6 > ~3.5 A. Cannot
    distort electronics (no charge). Directly counters the opening. PRIORITY 1 (cleanest, purely geometric).
(a) NEUTRAL DIPOLAR cradle: amide/hydroxyl/water surrogates (net-neutral, partial charges + LJ) cradling
    the carboxylates + ring OH, buffering the K-monopoles' polarising field. PRIORITY 2.
(d) ORIENTED neutral H-bond donors at the breaking/forming region. PRIORITY 3 (accent).
(c) WEAK distributed partial charges (+-0.1-0.5) forming a fuller cavity. PRIORITY 4 (refinement).
(e) COMBINATION (b)+(a)+(d) = the enzyme-like charge+cradle+pocket architecture. THE GOAL.

## Steps
0. cradle_lib: build the neutral surrogate templates (LJ-only wall bead; NMA amide dipole; methanol/phenol
   OH; TIP3P-like water), each with frozen geometry, partial charges (net 0 for dipoles, 0 for LJ walls),
   LJ params. Documented, reusable.
1. 01_placement: choose cradle site positions. START with the pure STERIC cap (b): place neutral LJ beads
   to bound C1-C6 in [2.6,3.5] without adding charge. Placement can be geometric (around the C1-C6 axis) -
   does NOT need the MILP (it is neutral, not a charge-selection problem). If we later want CERTIFIED
   cradle-site selection, pose it as discrete linear (excluded-volume + coverage) - deferred.
2. 02_relax: relax (fixed K-charges + cradle) with the validated MM-surrogate machinery (relax_higherK.py
   generalised). Measure relaxed C1-C6. Same LOOSE conv, frozen surrogates, QM -2.
3. 03_compare: C1-C6 recovery vs the bare-K baseline (K2 ~4-5). Does the cradle pull it into [2.6,3.5]?
   Sweep cradle "strength" (n LJ beads / wall stiffness) -> the interpolation from bare-K to enzyme-like.

## Method
Relaxation: B3LYP-D3BJ/def2-SVP/CPCM(eps=4) (unchanged), OR r2SCAN-3c if the running calibration passes.
Cradle surrogates FROZEN (like the charge surrogates); only the substrate relaxes. QM stays 24 atoms, -2.
Steric-realizability check (surr-surr, surr-substrate >= 2.0 A) MANDATORY before every relax (the K3/K4
clash lesson + the surrogate-sized min_dist discipline).

## Decision / success criterion
Cradle (esp. the steric cap b) that pulls relaxed C1-C6 into [2.6,3.5] while keeping the substrate a single
connected molecule = the fix the whole 05b arc was chasing. Then: does the ideal cradle depend on K (does
K2 need more than K6)? = the "interpolation between bare charges and full enzyme" thesis result.

## Certificate note
The K-charge selection stays certified (05/05b). The cradle is NEUTRAL geometric structure; its placement is
NOT a charge-optimisation and does not touch the certificate. If we want certified cradle-site selection
later, it is a separate discrete-linear coverage problem - deferred, not required for the geometric result.

## PROTOCOL NOTE (added): per-K coordinate testing + K20
For EVERY K, BOTH a_pot AND a_force must be run and relaxed to see which gives the better (most closed,
valid) baseline AT THAT K - the winning coordinate FLIPS with K and cannot be assumed:
  K2: a_force better (4.06 vs a_pot 5.28)
  K6: a_pot better (4.24 valid; a_force dissociates)
  K10: both invalid so far (a_pot fragments, a_force dissociates) - awaiting final states
So the cradle (and any future charge round) wraps the BEST VALID baseline per K, determined empirically,
not a fixed coordinate. ADD K20 to the charge-count ladder (K2/K6/K10/K20) - pushes further toward the
GOCAT working regime (81); K20 x {a_pot,a_force} to be generated with the surrogate-sized script and
relaxed, same as K6/K10, then cradled.
