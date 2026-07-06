# Methods chapter structure (living outline)

Populated section-by-section as each workflow step passes. Heading NAMES and
grouping are fixed; NUMBERING is placeholder (X.*) until the thesis chapter
number and template depth are confirmed.

Status per subsection: reserved | drafted | finalised

---

## Computational Methods (chapter)

### X.1 System preparation

#### X.1.1 Input structures and provenance  [drafted]
Covers: step 01 (download 2CHT/1DBF + CP2K ligands, checksums, manifest) and
step 01b (ligand-TSA registration, folded in as provenance footnote).
Supporting table: Table X - raw input structures.
Draft paragraph: written (see thesis draft / discussion notes).

#### X.1.2 Protein structure repair  [drafted]
Covers: step 02 (A/B/C terminal repair using 1DBF; global superposition +
terminal graft; join closure). J/K/L route NOT included (cut).
Supporting: alignment RMSD, peptide-continuity table.

#### X.1.3 Structure cleanup and validation  [drafted]
Covers: step 03 (altloc resolution, hydrogen removal, heavy-atom validation,
active-site presence, peptide continuity).

#### X.1.4 Ligand placement  [reserved]
Covers: step 04 (chorismate superposed onto A/B/C TSA pose; fit RMSD; contact
audit). This is where the J/K/L->A/B/C superposition is described.

#### X.1.5 Protonation and parameterisation  [drafted]
Reserved: H++ protonation, Antechamber/GAFF ligand parameterisation
(chorismate net charge = -2, -nc -2). ff14SB assignment.

### X.2 Molecular dynamics

#### X.2.1 System solvation and topology  [reserved]
Reserved: tleap build, ff14SB + GAFF + TIP3P, ion placement (record RNG seed),
box definition.

#### X.2.2 Minimisation and equilibration  [reserved]
Reserved: restrained minimisation series, NVT heating, NPT equilibration
(restraint ladder).

#### X.2.3 Production simulation  [reserved]
Reserved: accepted 1 ns unrestrained NPT production segment (canonical MD
reference; record seed; note 17-20 ns branches were audited and not accepted).

### X.3 QM/MM calculations  [reserved]
Reserved: QM-region definition (substrate +/- Glu78/Arg90 - see JPCB precedent),
method/basis, ORCA setup, single-point and optimisation.

### X.4 Reaction-path and TS characterisation  [reserved]
Reserved: constrained product search, NEB-TS, frequency validation of the TS.

---

## Phase 2 - Ligand-centred grid (chapter or major section, TBD)  [reserved]
Reserved: SDF construction, marching-cubes shells, area-weighted sampling,
tangential thinning. Charge-optimization method undecided - do not pre-write.
