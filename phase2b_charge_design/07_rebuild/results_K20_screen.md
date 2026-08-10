# Clean-room result: certified K=20 field electrostatically tears the chorismate dianion

## Verdict
The certified-optimal K=20 point-charge design (Layer 3b MILP, gap = 0.000) **fragments the bare
chorismate dianion** on relaxation. Reactant endpoint screened through the rebuilt oracle:
`integrity = DISTORTED`, 13 non-reacting bonds broken (ring + both carboxylates + several C-H),
2 spurious bonds formed. Genuine electrostatic tearing, not an artifact.

## Why this result is trustworthy (unlike prior fragmentation results)
- **Standoff grid (L2):** all charge sites >= per-atom LJ-minimum from every substrate atom
  (min_dist 2.65 A); no charge inside the wall.
- **Verified evaluator (L1):** zero-charge evaluate() reproduces L0 exactly (O3-C4 1.447; anchor 17.47).
- **Correct integrity test:** all-bond change vs reference, ignoring only the two reacting bonds.
- **LJ wall confirmed active on QM atoms:** direct QM/MM gradient test (1 QM atom + 1 +1 MM site at
  2.0 A) gives force -0.953 Eh/Bohr on the QM atom (strongly repulsive, pushed away). ASH elstat
  explicitly includes the QM-MM Lennard-Jones interaction. The fragmentation is NOT a wall-absence
  artifact.

## Mechanism
The wall resists close approach of any single atom but cannot prevent the COLLECTIVE tearing of a
connected molecule by 20 strong (+/-1) charges pulling from all directions (mean substrate
displacement 4.58 A). The design is certified-optimal for pure differential stabilisation (sum q.dV);
that objective, maximised with +/-1 point charges, is destructive to substrate integrity.

## Scientific content
Maximal Sokalski differential stabilisation and substrate integrity are in direct tension for the
bare dianion: the provably-optimal pure-stabilisation field fragments the substrate. This is the
clean, certified demonstration that GOCAT's integrity penalty (Phi_RMSG, at every stationary point)
is NECESSARY, not optional. The force-integrity constraint direction is legitimately motivated, now
on trustworthy machinery.

## Provenance
- Grid: L2_grid/grid_final.tsv (17028 sites, LJ-min standoff, |dV|max 0.0096 a.u.)
- Design: L3_solver/design_K20_* (certified MILP, gap 0, net-neutral, charges 3.2-3.9 A on LJ
  surface; +1 rediscovers the Arg90-like ether-O position)
- Screen: L1_oracle/run_K20_reactant/K20_reactant_result.txt (E -843.146 Eh, DISTORTED)
- Wall diagnostic: L1_oracle/wall_check.out (fx -0.953 Eh/Bohr, wall active on QM)
