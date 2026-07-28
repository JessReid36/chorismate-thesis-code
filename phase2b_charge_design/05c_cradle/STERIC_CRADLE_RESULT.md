# 05c cradle: heuristic steric-cage RESULT — steric confinement alone does NOT recover geometry

## Converged results (both jobs finished)
  K2aforce + 6-bead cage : C1-C6 5.306  O3-C4 1.436 (intact)  SINGLE molecule -> OPEN (worse than bare 4.06)
  K6apot   + 5-bead cage : C1-C6 5.931  O3-C4 1.372  FRAGMENTED (23/24)        -> OPEN + fragmenting

## Finding: neutral steric confinement alone is INSUFFICIENT
The adaptive 6-bead cage (q=0 LJ walls, all 3 opening DOFs covered, strong walls guard/sigma~1.1) did NOT
hold the substrate compact. Both cases opened PAST the cage:
  - K2aforce ended MORE open WITH the cage (5.31) than the bare charge design without it (4.06)
  - K6apot opened to 5.93 and began fragmenting
So the field pries the substrate open along conformational routes the fixed cage does not intercept: the
ring flexes and the reacting carbons slip past/around the frozen walls rather than colliding with them. A
flexible field-driven substrate has more escape routes than a 6-bead frozen cage plugs.

## Why (mechanism)
Frozen LJ beads act only via repulsion when a substrate atom pushes INTO a bead. Beads placed around the
COMPACT geometry do not track the substrate as it distorts; the opening motion is not purely the DOFs the
cage guards, so the carbons evade the walls. Pure steric caging of a flexible substrate is leakier than the
3-DOF analysis implied.

## Combined picture (both single-mechanism approaches now fail)
  - CHARGES ALONE (05b): cannot recover near-attack geometry (K2-K6 open ~4A; K10 over-polarises/fragments)
  - STERIC CRADLE ALONE (05c heuristic): cannot hold it either (opens past the cage / fragments)
Neither pure electrostatics nor pure sterics (as formulated) recovers the reactive conformation.

## Implications for next steps
1. Pure-steric neutral cage is FALSIFIED as a standalone geometry fix (as built - frozen LJ, this grid).
2. Path C (test first): CONTINUOUS fractional charges (GOCAT model) - can a gentle DISTRIBUTED fractional
   field do what harsh +-1 could not? Note: continuous alone with a raw max-lowering objective RE-RAILS to
   +-1 (no reason to pick small q) -> needs a gentleness mechanism (preorg penalty or magnitude-reg) to
   actually shape a soft field. If a gentle GOCAT-scale field holds the geometry -> no cradle needed.
3. If Path C also fails -> the cradle must be DIPOLAR/active (oriented net-neutral dipoles that H-bond/orient
   the substrate), not passive steric walls; OR a MOVING/adaptive cage; OR accept external environment
   cannot fully preorganise this substrate.

## Certified-cover note
The certified set-cover cages (03_certified_cover, union grid) were NOT relaxed - the heuristic (stronger-
wall) cage already failed, so the weaker grid-restricted certified cages would also fail. Certified-cover
machinery retained as method; relaxation deferred/moot pending a working confinement mechanism.
