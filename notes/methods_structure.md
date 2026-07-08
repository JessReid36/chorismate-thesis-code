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

#### X.1.4 Ligand placement  [drafted]
Covers: step 04 (chorismate superposed onto A/B/C TSA pose by Kabsch; fit RMSD
0.08-0.15 A; inter-subunit catalytic-contact audit vs the published Agbaglo
QM-cluster set {7,57,59,60,63,73,74,75,78,90,108,115}). The J/K/L->A/B/C
superposition is described here. Active site is inter-subunit (C3 pattern: A
completed by C, B by A, C by B); Arg116 is reported as second shell (~5.5-6.8 A),
not an active-site/QM residue. Draft paragraph: written (see thesis draft).

#### X.1.5 Protonation and parameterisation  [drafted]
Protonation (steps 07a/07b): protein-only structure submitted to H++ (pH 7.0,
0.15 M, internal diel 10, external diel 80). Server-stripped chain IDs and the
generic HIS names were restored/assigned downstream: chains + 1-127 numbering
recovered by residue-order mapping onto the 07a input; histidines named
HID/HIE/HIP from the placed ring protons (His36/His54 = HIP, His106 = HIE; none
active-site). Drift is classified on three axes by exact set membership
(backbone/side-chain, flippable/non-flippable, active-site/not): backbone,
active-site, active-site-backbone and non-flippable heavy-atom drift are all
0.000 A - the fold and every catalytic residue are provably untouched - and the
only shifts >1 A (max 2.245 A) are Gln44/Gln101 amide flips (flippable,
non-active-site, benign). 378 peptide links intact incl. step-02 graft joins;
6 CYS reduced (no disulfide). An explicit acceptance rule (no missing/extra
heavy atoms; backbone <=0.05, non-flippable <=0.2, active-site-backbone <=0.05 A;
zero active-site atoms >1 A; no non-flippable or active-site large shift) returns
ACCEPT. Protein trimer net charge -3 at pH 7. Accepted structure:
abc_protonated_hpp_accepted.pdb. Draft paragraph: written (see thesis draft).
Parameterisation (steps 08a/08b): chorismate (residue CHA, 24 atoms, net -2)
GAFF-typed with Antechamber + parmchk2 (0 ATTN, 0 MISSING; types
c,c2,c3,ce,h1,ha,ho,o,oh,os). Charges are AM1-BCC - the standard Antechamber/GAFF
scheme and the charge model implied by the protocol ("Antechamber/GAFF"). AM1-BCC
needs sqm, unavailable on the HPC (missing libopenblas.so.0), so charges were
derived off-HPC on a local Ubuntu AmberTools (08a) and read into the HPC
GAFF-typing run via -c rc (08b); this split is a recorded deviation. The MOE
template charges carried in the CP2K mol2 were replaced. Outputs: cha_gaff.mol2 +
cha.frcmod. Draft paragraph: written (see thesis draft).

### X.2 Molecular dynamics

#### X.2.1 System solvation and topology  [drafted]
Covers: step 09 (09a combine + 09b tleap build/audit; collapses attempt_3's
16/17/18/18b/18c). ff14SB (protein) + GAFF (CHA, AM1-BCC charges from step 08)
+ TIP3P assigned in tleap. 09a combines the accepted protonated trimer (07b)
with the 3 placed chorismates (06) into one pre-tleap complex (6255 atoms =
6183 protein + 72 CHA); atom/charge counts derived from inputs, not hard-coded.
09b writes a dry topology (NATOM 6255, net -9) then solvates: TIP3P box with a
10 A solute buffer (solvatebox -> rectangular, actual 92.9 x 76.7 x 95.4 A),
16,472 waters, neutralised with 9 Na+ (addions Na+ 0; no background electrolyte,
matching the protocol). Final solvated system 55,680 atoms, net 0. addions is
deterministic - there is NO ion-placement RNG seed (the Langevin seed belongs to
X.2.2). Counts and charge are read from the prmtop/inpcrd (the solvated PDB
saturates at 9999 residues), not the PDB. AMBER22 substituted for the protocol's
AMBER18 (identical ff/water/charge choices). tleap warnings classified and
benign: 12 terminal-name conversions (cosmetic), 2 pre-neutralisation dry-charge
notices (expected), 7 solute-internal H-atom close contacts from H++ placement
(relaxed by the first restrained minimisation); 0 missing-parameter/unknown-type,
0 errors. Draft paragraph: written (X.2.1 prose).

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
