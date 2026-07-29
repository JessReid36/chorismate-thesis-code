# 05c cradle Option A: CERTIFIED set-cover cage selection

## What (user's insight made rigorous)
Heuristic bead placement (01_placement) is uncertified. A neutral bead's STERIC confinement of the RELAXED
geometry is nonlinear -> not directly certifiable. Option A pushes the geometry nonlinearity OFFLINE into
per-site GUARD LABELS (does a bead here block opening-DOF k?), leaving a pure LINEAR integer SET-COVER that
HiGHS solves to gap 0.000. Cage SELECTION is certified (cover all 6 opening (DOF,side) targets + no-clash,
min beads); Tier-2 relaxation confirms confinement. The solver is NOT capped on bead count (it may pick as
many as needed; it minimises count as the objective).

## Candidate grid = UNION (charge lattice + inner shell)  [corrected design]
build_bead_grid_union.py builds the bead-candidate set as the UNION of:
  (a) the charge dv_grid points  -> beads share the SAME lattice as charges (direct charge-vs-bead compare)
  (b) a dense INNER shell reaching to 1.9 A from the substrate -> the neutral-bead-only close region that
      charges cannot occupy (strong walls).
435 points = 98 charge-grid + 337 inner-shell; 96 strong-wall candidates (guard<3.5). Built from the correct
24-atom reactant (C1=0,C6=12,O3=7,C4=8; C1-C6 3.124).
[An earlier bead-only grid gave a spurious "K10 inpl+ uncoverable" result -> that was a GRID ARTIFACT
 (no candidate point in that direction), NOT a physical limit. The union grid covers it. Do not claim
 "cradle impossible at high K": with an adequate grid all tested K admit a full cage.]

## Certified results on the UNION grid (gap 0.000 all)
  K2/a_force : FULL 6-target cover, 3 beads (multi-purpose)
  K6/a_pot   : FULL 6-target cover, 3 beads
  K10/a_pot  : FULL 6-target cover, 3 beads
=> All three charge counts admit a full certified neutral cage. The set-cover minimises bead COUNT -> 3
multi-purpose beads (each guards 2 targets) vs the heuristic's 6. Beads sit guard ~3.9-5.8 A.
NO "charges preclude the cradle" claim - the earlier such result was an under-resolved-grid artifact.

## On co-adaptation (user's point)
Charges are placed FIRST and FROZEN, then the cradle fits the remaining space -> the cradle is constrained
by (dictated by) the charge set. Any crowding is a consequence of STAGING, not a physical law: joint
charge+bead selection (one linear program, deferred) could trade a charge position for a better cage. The
staged certified cover is the current deliverable; certified JOINT selection is the stronger future version.

## Files
build_bead_grid_union.py, bead_grid_union.tsv (435 pts, shared lattice) ; certified_cradle_cover.py ;
certU_<K>_beads.xyz (the cages). [build_bead_grid.py + bead_grid.tsv = earlier bead-only grid, retained but
superseded by the union grid.]

## Next
Assemble (charges + certified union cage) coords, relax, compare relaxed C1-C6 to (i) bare baseline
(K2af 4.06, K6ap 4.24), (ii) heuristic cage (01_placement, running). HOLD certified relaxations until the
heuristic cradle jobs converge (stronger walls; if they fail, weaker certified cages will too).

## Certificate status
Cage SELECTION certified (linear set-cover, gap 0.000). Confinement is Tier-2 (validated by relaxation).
