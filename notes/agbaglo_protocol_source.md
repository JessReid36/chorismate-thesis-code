# PROTOCOL SOURCE PAPER (critical - was previously confused with Claeyssens)

## Your actual MD/parameterization protocol source:
Agbaglo, Summers, Cheng, DeYonker (2024) "The Influence of Model Building Schemes
and Molecular Dynamics Sampling on QM-cluster Models: The Chorismate Mutase Case
Study." Phys. Chem. Chem. Phys. 26(16): 12467-12482. DOI 10.1039/d3cp06100k.
(NIH author manuscript nihms-1987271; PMC 2025-04-24.)

## Two DIFFERENT reference papers - do not conflate:
- Agbaglo/DeYonker 2024 = PROTOCOL source (AMBER18, ff14SB, GAFF, TIP3P, H++,
  1DBF repair, A/C interface, min/equil/production ladder). Pure AMBER. THIS is
  what the system_development workflow reproduces.
- Claeyssens 2011 = SCIENCE anchor only (TS stabilization electrostatics, Arg90/
  Arg7/Glu78). Uses CHARMM27 - NOT the protocol source. Cite for the science,
  never for the AMBER method.

## Exact protocol from Agbaglo 2024 (page 7), matches our workflow:
- Structure: 2CHT, chain C missing residues replaced from 1DBF.
- Protonation: H++ server, default parameters (ref 76).
- Substrate: native chorismate in pre-reactive conformation (NOT the TSA).
- MD: AMBER18 (ref 77), ff14SB, PBC, 9 A nonbonded cutoff.
- Ligand: Antechamber + GAFF (refs 77, 78 = Wang 2004). CHARGE METHOD NOT
  EXPLICITLY STATED - "Antechamber/GAFF" only. AM1-BCC is the standard default
  and the defensible reading; cite Wang 2004 (GAFF) + Jakalian (AM1-BCC) for the
  methods themselves, note consistency with Agbaglo 2024. Do NOT claim the paper
  said AM1-BCC.
- Solvation: TIP3P (ref 79), cubic 10 A box, 9 Na+ neutralization (ref 80).
- Minimization: protein heavy atoms restrained at kpos = 200 kcal/mol/A^2.
- Equilibration: restraints relaxed over FIVE 20 ps NPT (Langevin) runs at
  300 K / 1 atm; SHAKE on H bonds (ref 81).
- Production: 20 ns free run, 1 ps/frame = 20,000 frames.
- Analysis: cpptraj (ref 82) for protein RMSD. [Our cpptraj is currently broken -
  see software_environment.md]
- Active site: Chain A/C interface (least conformational fluctuation of the three).
  CONFIRMS our inter-subunit finding (Arg63' etc. from adjacent chain).

## IMPORTANT consequence for our workflow:
Original attempt_3 KEPT CP2K USER_CHARGES (charge diff = 0). That was a DEVIATION
from Agbaglo 2024, which parameterizes chorismate fresh with Antechamber/GAFF.
Switching to AM1-BCC now brings us CLOSER to the real protocol paper. Good.

## Benchmark values from Agbaglo 2024 (for our Results comparison):
- Mean DG-doubledagger (250 MD QM-cluster models) = 10.34 +/- 2.62 kcal/mol.
- Mean DG-rxn = -15.38 +/- 3.40 kcal/mol.
- Experimental DG-doubledagger = 15.4 kcal/mol at 25 C (ref 86 = Kast/Hilvert).
- QM level: B3LYP-GD3BJ / small Pople basis / CPCM (they note it underestimates
  the barrier; QM level quality deliberately not the focus).
- F-SAPT residue interaction energies (dianionic chorismate): Arg7 -140.5,
  Arg63 -133.2, Arg90 -110.3, Lys60 -78.1, Tyr108 -15.5, Thr74 +10.9 kcal/mol.
  -> directly supports charge-placement/DDE-doubledagger rationale.
