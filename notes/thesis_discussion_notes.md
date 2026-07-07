# Thesis discussion notes and reference tracker

Running file for talking points, open decisions, and papers to read/cite.
Populate as the project develops. Not part of the runnable workflow — this is
a thinking/writing aid for the Discussion and Introduction chapters.

---

## Open decisions to revisit

| Item | Current state | Decision point | Notes |
|---|---|---|---|
| Chorismate reactant coordinates | Using CP2K tutorial mol2 (LigA–C, 24 atoms) as placeholder | End of system-development phase | Provenance is pedagogical, not primary. Options if replaced: (A) anchor geometry to crystallographic TSA pose in 2CHT; (B) build from PubChem/ChEBI chorismate + align to TSA pose; (C) validate + re-provenance CP2K coords. Reactant charge confirmed = **2−** (dianion). Charges parameterised **AM1-BCC** in step 08 (MOE template charges from the CP2K mol2 replaced; derived off-HPC, see step08a_provenance.txt); **geometry still placeholder**. |
| Grid reference ligand | Built around placeholder chorismate | After QM/MM-optimized reactant + NEB-TS available | Grid must be regenerated on the real optimized reactant and TS geometries before any charge-placement claims. |
| Charge-optimization method | Undecided | Phase 2, later | Candidates: ORCA single-points, point-charge embedding, QM/MM point-charge perturbation. Do not write as finalized. |

---

## Discussion talking points

### 1. The reactant-destabilization vs TS-stabilization debate (CENTRAL)

The literature on chorismate mutase contains a genuine, named controversy about
*where* the catalytic effect comes from, and this project's ΔΔE‡ framework sits
directly on top of it.

- The Pauling paradigm holds that enzymes accelerate reactions by binding the
  transition state (TS) more tightly than the substrate, lowering the activation
  barrier relative to the uncatalyzed reaction.
- This paradigm was **challenged** for chorismate mutase: some computational work
  predicted that the decisive factor is **ground-state (reactant) destabilization**
  — i.e. the enzyme strains/compresses the bound substrate — rather than TS
  stabilization.
- Later crystallographic and QM/MM work argues the opposite: catalysis is driven
  mainly by **electrostatic transition-state stabilization**. The strongest
  experimental support is the Arg90→citrulline variant (see refs): removing one
  cationic active-site charge cripples catalysis (kcat and kcat/Km down by >4
  orders of magnitude) while barely affecting substrate binding — showing the
  positive charge stabilizes the TS, not the ground state.
- Several QM/MM studies land on a **combination**: a *minor* substrate
  destabilization by compression plus a *major* electrostatic TS stabilization.

**Why this matters for the thesis:** the whole point of the ΔΔE‡ approach is to be
specific about whether an added charge lowers the TS energy, raises the reactant
energy, or changes the barrier — and by how much. The project therefore addresses
a real open question rather than a settled one. In the write-up, position the work
against BOTH sides of this debate; do not present TS stabilization as uncontested.

Reminder for the whole write-up: never say "the charge lowers the energy" without
stating whether it lowers the TS, raises/lowers the reactant, or changes the
activation barrier.

### 2. Enzyme scaffold vs reactant come from different sources

2CHT provides the protein and the active-site location (via the bound TSA
transition-state analogue). The chorismate reactant geometry comes from a separate
source (currently the CP2K tutorial files). State this explicitly wherever 2CHT is
introduced, so it is never read as "modelled the inhibitor as the substrate."

### 3. Relationship to "optimal catalytic field" prior work

At least one prior study (Sokalski/DTSS group) computes the *optimal catalytic
field* — where charge/field placement best stabilizes the TS. This is the closest
published analogue to the grid-based charge-placement idea in this project. MUST
read carefully and articulate clearly how this project differs (e.g. explicit
grid-sampled candidate positions around the real optimized reactant/TS geometries
vs their analysis). This is both a precedent and the work to differentiate from.

---

## Papers to read / potential citations

Verify each yourself before citing — pulled from literature search, confirm they
say what they're cited for.

### Structural (primary sources for 2CHT)
- **Chook, Ke, Lipscomb (1993)** *Crystal structures of the monofunctional
  chorismate mutase from B. subtilis and its complex with a transition state
  analog.* PNAS 90(18): 8600–8603. DOI 10.1073/pnas.90.18.8600. PMID 8378335.
  → PRIMARY citation for 2CHT. Structure at 1.9 Å, R ≈ 0.201, 12 monomers/asym unit.
- **Chook, Gray, Ke, Lipscomb (1994)** *The monofunctional chorismate mutase from
  B. subtilis...* J. Mol. Biol. 240(5): 476–500. DOI 10.1006/jmbi.1994.1462.
  PMID 8046752. → Fuller structure determination (TSA + prephenate complexes),
  spells out the pericyclic mechanism and active-site detail. Better for the
  "why it works" mechanism discussion than the terse PNAS communication.

### Computational / mechanistic (core to Phase 2 rationale)
- **Ranaghan, Ridder, Szefczyk, Sokalski, Hermann, Mulholland (2004)**
  *Transition state stabilization and substrate strain in enzyme catalysis:
  ab initio QM/MM modelling of the chorismate mutase reaction.* Org. Biomol.
  Chem. 2(7): 968–980. DOI 10.1039/b313759g. → KEY methodological precedent for
  QM/MM of this exact system. Barrier estimates 7.4–11.0 (MP2) and 12.7–16.1
  (B3LYP) kcal/mol vs experimental ΔH‡ = 12.7 ± 0.4 kcal/mol. Benchmark to
  compare own barrier against.
- **Kienhöfer / Lipscomb group (2014)** *Electrostatic transition state
  stabilization rather than reactant destabilization provides the chemical basis
  for efficient chorismate mutase catalysis.* PNAS. DOI 10.1073/pnas.1408512111.
  → Arg90→citrulline experiment. KEY motivation for charge-placement rationale:
  direct experimental evidence that an active-site positive charge stabilizes the
  TS more than the reactant.
- **Sokalski group (2004)** *Differential Transition-State Stabilization in Enzyme
  Catalysis: Quantum Chemical Analysis... and Prediction of the Optimal Catalytic
  Field.* PMID 15584751. → DTSS analysis, ~−23 kcal/mol TS stabilization
  (MP2/6-31G(d)) relative to bound substrate. Closest analogue to the grid/
  charge-field idea — read carefully, differentiate from.
- **DFT-based QM/MM applied to chorismate mutase.** J. Phys. Chem. B.
  DOI 10.1021/jp036236h. → Tested QM-region choices: substrate only vs substrate
  + Glu78/Arg90 charged side chains. Informs QM-region selection decision.
  Concludes minor substrate destabilization + major electrostatic TS stabilization.

---

## Reference tracking table (populate as write-up proceeds)

| Ref (short) | Full citation / DOI | What it supports in the write-up | Chapter/section | Read? | Verified? |
|---|---|---|---|---|---|
| Chook 1993 | PNAS 90:8600, 10.1073/pnas.90.18.8600 | 2CHT source structure; enzyme is homotrimer, analogue-bound | Methods (inputs) | ☐ | ✅ |
| Chook 1994 | JMB 240:476, 10.1006/jmbi.1994.1462 | Mechanism, active site, TSA + prephenate complexes | Intro / Discussion | ☐ | ☐ |
| Ranaghan 2004 | OBC 2:968, 10.1039/b313759g | QM/MM precedent; barrier benchmark | Methods / Discussion | ☐ | ☐ |
| Kienhöfer/Lipscomb 2014 | PNAS, 10.1073/pnas.1408512111 | TS-stabilization evidence; charge-placement motivation | Intro / Discussion | ☐ | ☐ |
| Sokalski DTSS 2004 | PMID 15584751 | Optimal catalytic field; grid-idea precedent | Discussion / Phase 2 | ☐ | ☐ |
| DFT QM/MM CM | JPCB, 10.1021/jp036236h | QM-region choice (substrate ± Glu78/Arg90) | Methods (QM/MM) | ☐ | ☐ |

---

## Verify-later checks (post-minimisation / post-MD)

- **Step 02 terminal graft joins (A115-A116, B1-B2):** repair closed these to a
  canonical C-N distance by rigid translation of the graft block (moved 0.79 A
  and 1.32 A respectively). This fixes the C-N DISTANCE only, not the bond
  angles/dihedrals. After energy minimisation, confirm these joins relaxed to
  physical backbone geometry. If minimisation cannot relax them, the terminal
  repair needs revisiting.
- **Active-site side-chain conformation:** Step 02 confirmed active-site residues
  are PRESENT (7,57,59,60,63,73,74,75,78,90,108,115). Their conformational
  integrity - especially catalytic Arg90 and Glu78 - is validated only later by
  visual inspection and MD stability, not by the presence check.

- **Origin of removed hydrogens (step 03):** the 601 hydrogens removed were
  genuine, identified by the authoritative PDB element column (not name-guessing).
  They entered via the 1DBF-derived terminal residues (1DBF is 1.30 A and
  hydrogen-bearing). Removal was correct ahead of controlled protonation, and was
  explicitly verified rather than assumed. Note for provenance: the input to
  step 03 is the REPAIRED structure (2CHT+1DBF), not raw 2CHT - raw 2CHT (1.9 A,
  1994) carries no deposited hydrogens.

- **Step 04 catalytic contacts - RESOLVED by measurement:** the earlier
  verify-later item (confirm Arg90/Glu78 contact) is now answered directly.
  Closest heavy-atom approaches, consistent across all three sites:
  Arg90 2.8-3.0 A, Glu78 2.9-3.4 A, Arg7 2.7-2.8 A, Tyr108 2.6-3.0 A, Val73' 3.5 A,
  Thr74' 3.3-3.6 A, Cys75' 3.2-3.4 A (all direct contact). Arg63' is a genuine
  cross-chain contact at 2.9-3.1 A in sites A/B and a swung rotamer at 6.4 A in
  site C (see step 05/05b). Arg116 sits at 5.5-6.8 A - SECOND SHELL, reported for
  provenance but not an active-site/QM residue.
  RESOLVED (active-site set): matched to the published Agbaglo QM-cluster set
  {7,57,59,60,63,73,74,75,78,90,108,115}. Val73 IS in the set (now measured);
  Arg116 is NOT (reclassified second shell); Arg63 IS in the set (the earlier
  ~21 A same-chain reading was an artifact of same-chain-only measurement).
- **Still to check visually (ChimeraX, next step):** overall pose sanity and
  Arg116 second-shell position look reasonable.

## Inter-subunit active site (CONFIRMED - important)

The B. subtilis CM active site is INTER-SUBUNIT: each of the three sites sits at
the interface of two adjacent monomers (Chook 1993 PNAS; Chook 1994 JMB). Residue
origin relative to the substrate's own subunit:
  - same-subunit: Arg7, Arg90, Glu78, Tyr108, Leu115
  - adjacent-subunit ('): Arg63', Lys60', Val73', Thr74', Cys75', Phe57', Ala59'
  - second shell (measured, NOT in the active-site/QM set): Arg116 (~5.5-6.8 A)
Step 04 measured closest contacts across ALL chains and recovered this exactly
(A's site completed by C, B's by A, C's by B - the C3 pattern), zero mismatches.

CORRECTION LOGGED: an earlier same-chain-only measurement wrongly reported Arg63
at ~21 A and dropped it as non-contact. Arg63 IS a key active-site residue,
contributed in trans by the adjacent subunit. The drop was a measurement artifact,
now corrected. This is why cross-chain measurement matters for every ' residue.

Residue roles (BRENDA / FMO QM/MM): Arg63, Arg7, Tyr108 orient the substrate;
the Glu78-Arg90-substrate arrangement controls H-bond strength and TS
stabilisation, with Arg90's positive charge polarising the substrate in the TS.
Arg116 is cited elsewhere as substrate-orienting but sits in the second shell
here (~5.5-6.8 A) and is outside the Agbaglo QM-cluster set, so it is not treated
as an active-site/QM residue.
-> directly supports the DDE-doubledagger charge-placement rationale.

## Verify-later: cha_c Arg63' rotamer

Site C's Arg63' points its guanidinium AWAY from the substrate (nearest approach
6.4 A via CB, vs 2.9-3.1 A via NH1 in sites A/B) - a rotamer flip, matching the
visual ChimeraX observation. OPEN QUESTION: is this inherited from raw 2CHT
(crystallographic heterogeneity, fine to report) or introduced by our repair/
placement (artifact, must fix)? Reproduce the old step-09 raw-TSA-vs-placed check
in the streamlined tree before freezing the complex. Also worth tracking whether
MD equilibration converges site C's Arg63 toward the other two.

## Arg63 site-C rotamer - RESOLVED (crystallographic, proven)

Chain of evidence (all measured, not asserted):
  1. Active site is inter-subunit; step 04 recovered same/cross residue origins
     with zero mismatches.
  2. Step 05: in raw 2CHT, Arg63' is conformationally heterogeneous across the 12
     TSA sites (guanidinium-analogue 2.9-7.3 A; 7/12 swung >4.5 A). The cha_c
     source site (raw L211) was already swung at 7.28 A in the crystal.
  3. Step 05b: each repaired A/B/C Arg63 side chain is identical to raw 2CHT
     (side-chain RMSD 0.000 A, all three chains) - repair/cleanup did NOT move them.
  => site-C swung Arg63' is INHERITED crystallographic flexibility of a solvent-
     exposed inter-subunit arginine, NOT a workflow artifact. Safe to freeze.

Residual (genuine MD-stage check): during MD analysis, note whether site-C Arg63'
converges toward the A/B conformation or stays swung. Either outcome is reportable.

## Protonation (step 07) - pH and dielectric decisions

- **pH 7.0 chosen** (not H++ default 6.5). Justification: physiological pH, and the
  closest-matching 2CHT/chorismate preparation in the literature (Steinmann/
  Claeyssens lineage) protonated at pH 7; general enzyme-MD convention is pH 7.
  Reproduces our own original attempt_3 run (0.15_80_10_pH7). NOTE for write-up:
  this is a field-convention choice, NOT strictly Agbaglo's "default parameters"
  (H++ literal default = 6.5). State honestly in limitations.
- **Dielectrics: internal 10, external 80** (VERIFIED against H++ literature - the
  protein interior is the LOW value ~10-20, solvent/water is 80). Filename order
  0.15_80_10 lists salinity, external, internal. Do not invert these.
- **H++ topology/coordinate files DISCARDED** - H++ uses its own internal force
  field/water model (ff19SB/OPC-style); we use only the protonated geometry and
  rebuild with ff14SB + TIP3P in tleap. This is the correct, paper-faithful choice.
- **Active-site protonation confirmed at pH 7** (from pKa table): Arg7/63/90/116 all
  protonated (+), Glu78 deprotonated (-), Tyr108 neutral. Catalytic electrostatics
  correct. Total protein charge -3; with 3x chorismate (-6) => complex -9 => 9 Na+
  neutralization (matches Agbaglo).
- **Citations to verify:** Gordon 2005 (NAR 33:W368); Anandakrishnan 2012 H++ 3.0
  (NAR 40:W537).
