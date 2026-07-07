# H++ predicted pKa values (pH 7.0, salinity 0.15 M, internal diel 10, external diel 80)

Input: hpp_input_protein_only_with_TER.pdb | Isoelectric point 5.34 | Total charge at pH 7 = -3

Active-site residues marked. State column = protonation at pH 7 (chain A representative).

| Residue | pKa A | pKa B | pKa C | State @pH7 | Active site |
|---|---|---|---|---|---|
| NT | 8.3 | 8.4 | 7.9 | protonated (+) |  |
| ARG-4 | >12.0 | >12.0 | >12.0 | protonated (+) |  |
| ARG-7 | >12.0 | >12.0 | >12.0 | protonated (+) | **yes** |
| GLU-13 | 2.6 | 3.3 | 3.4 | deprotonated (-) |  |
| ARG-14 | >12.0 | >12.0 | >12.0 | protonated (+) |  |
| ASP-15 | 1.0 | 1.1 | 1.8 | deprotonated (-) |  |
| GLU-17 | 3.0 | 2.3 | 2.5 | deprotonated (-) |  |
| GLU-18 | 4.2 | 4.3 | 4.2 | deprotonated (-) |  |
| GLU-19 | 0.9 | 0.2 | 1.9 | deprotonated (-) |  |
| LYS-23 | >12.0 | >12.0 | 11.5 | protonated (+) |  |
| LYS-25 | 11.4 | 11.8 | >12.0 | protonated (+) |  |
| GLU-29 | 3.5 | 4.2 | 2.9 | deprotonated (-) |  |
| LYS-30 | 11.6 | >12.0 | >12.0 | protonated (+) |  |
| GLU-33 | 4.9 | 4.7 | 4.5 | deprotonated (-) |  |
| GLU-34 | 3.1 | 2.7 | 2.0 | deprotonated (-) |  |
| HIS-36 | 7.1 | 7.1 | 7.0 | protonated (+) |  |
| LYS-38 | 11.3 | 11.4 | >12.0 | protonated (+) |  |
| GLU-40 | 4.3 | 4.4 | 4.1 | deprotonated (-) |  |
| ASP-41 | 2.6 | 2.1 | 1.9 | deprotonated (-) |  |
| ASP-52 | 0.1 | 1.1 | 0.6 | deprotonated (-) |  |
| HIS-54 | 9.1 | 8.4 | 9.3 | protonated (+) |  |
| LYS-60 | 11.3 | 11.9 | 11.2 | protonated (+) | **yes** |
| ARG-63 | >12.0 | >12.0 | >12.0 | protonated (+) | **yes** |
| GLU-64 | 3.8 | 2.9 | 4.6 | deprotonated (-) |  |
| TYR-70 | 10.8 | 10.4 | >12.0 | neutral |  |
| CYS-75 | >12.0 | >12.0 | >12.0 | neutral | **yes** |
| GLU-78 | <0.0 | <0.0 | <0.0 | deprotonated (-) | **yes** |
| ASP-80 | 2.5 | 2.4 | 2.1 | deprotonated (-) |  |
| LYS-86 | >12.0 | 10.9 | 11.6 | protonated (+) |  |
| LYS-87 | 10.7 | 11.8 | 11.1 | protonated (+) |  |
| CYS-88 | >12.0 | >12.0 | >12.0 | neutral |  |
| ARG-90 | >12.0 | >12.0 | >12.0 | protonated (+) | **yes** |
| ASP-98 | 2.6 | 4.7 | 4.1 | deprotonated (-) |  |
| ASP-102 | 3.5 | 3.4 | 3.2 | deprotonated (-) |  |
| ARG-105 | >12.0 | >12.0 | >12.0 | protonated (+) |  |
| HIS-106 | <0.0 | 0.5 | 1.0 | neutral |  |
| TYR-108 | 8.7 | 9.9 | 10.1 | neutral | **yes** |
| GLU-110 | 2.8 | 1.9 | 2.2 | deprotonated (-) |  |
| LYS-111 | 10.8 | 10.9 | 11.0 | protonated (+) |  |
| ARG-116 | >12.0 | >12.0 | >12.0 | protonated (+) |  |
| ASP-118 | 3.2 | 2.2 | 3.4 | deprotonated (-) |  |
| LYS-123 | 10.7 | 10.4 | 10.8 | protonated (+) |  |
| GLU-126 | 4.4 | 4.2 | 4.2 | deprotonated (-) |  |
| CT | 4.2 | 4.2 | 4.4 | deprotonated (-) |  |

## Active-site residues summary (all three chains consistent)

| Residue | pKa (A/B/C) | State @pH7 | Role |
|---|---|---|---|
| Arg7  | >12/>12/>12 | protonated (+) | substrate orientation |
| Arg63'| >12/>12/>12 | protonated (+) | cross-subunit, orientation |
| Arg90 | >12/>12/>12 | protonated (+) | TS electrostatic stabilization |
| Glu78 | <0/<0/<0 | deprotonated (-) | Glu78-Arg90-substrate H-bond network |
| Tyr108| 8.7/9.9/10.1 | neutral | orientation |
| Lys60'| 11.3/11.9/11.2 | protonated (+) | cross-subunit |
| Cys75'| >12/>12/>12 | neutral (SH) | cross-subunit |

*Active-site set = published Agbaglo QM-cluster residues (SI, S3). Arg116 is second shell (~5.5-6.8 A, step 04) and is NOT an active-site/QM residue. Val73 is in the set but non-titratable, so it does not appear in this pKa table.*

## Histidine tautomers (non-active-site) - CONFIRMED from H++ placed hydrogens (step 07b)

| Residue | pKa (A/B/C) | Assignment @pH7 (all 3 chains) | Evidence |
|---|---|---|---|
| His36 | 7.1/7.1/7.0 | HIP (doubly protonated, +) | HD1 + HE2 placed |
| His54 | 9.1/8.4/9.3 | HIP (doubly protonated, +) | HD1 + HE2 placed |
| His106| <0/0.5/1.0 | HIE (neutral, epsilon)     | HE2 only |

His36 note: pKa ~7.1 is borderline, but it resolves to HIP at both pH 7.0 (our
run) and pH 6.5 (H++ default) - fraction protonated ~0.56 and ~0.80 respectively
- so the pH choice does not change any histidine assignment. Recorded in
00_admin/step07b_his_assignment.tsv.
