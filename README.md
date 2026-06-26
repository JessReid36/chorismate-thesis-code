# Chorismate Mutase QM/MM + Charge-Grid Project — Code Archive

Source of truth for runnable code and provenance. Large data
(trajectories, restarts, topologies) lives outside git; see provenance
notes for paths and acceptance rationale.

Directory structure grows as steps are added — folders appear when the
first script for that step is committed, not before.

## Canonical references
- Protein: paper-aligned A/B/C trimer; repaired with 2CHT (scaffold)
  + 1DBF (terminal repair/alignment). J/K/L work is exploratory/provenance
  only — NOT retained final methodology.
- MD reference: accepted 1 ns unrestrained NPT production segment
  (segment_001_accepted). The 17-20 ns branches are NOT canonical (audit flags).
- Grid ligand: cha_a_from_liga.mol2 — PLACEHOLDER reactant reference,
  to be re-run on QM/MM-optimized reactant + NEB-TS geometries.

## Status flags
exploratory | provisional | retained-final | archived
