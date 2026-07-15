# Validating the Model Against the Enzyme: A Fidelity Framework

*Working notes for the Methods/Discussion chapter — how we justify that the computational system is a faithful proxy for Bacillus subtilis chorismate mutase, and what "faithful" is allowed to mean.*

## The problem this section answers

A QM/MM model of an enzyme is not the enzyme. It is a deliberate simplification: one solvated trimer in a periodic box, a fixed protonation state, a classical force field for the environment, a density functional for the reacting core, and a handful of nanoseconds of sampling standing in for a continuous biological process at 310 K. The living enzyme is a fluctuating quantum-mechanical object embedded in a crowded, heterogeneous cellular medium. No amount of computation collapses that gap entirely.

The question is therefore not whether the model *is* the enzyme, but whether it is faithful in the specific respects that bear on the claim being made. The claim in this work is narrow and well-defined: that the enzyme accelerates the chorismate-to-prephenate Claisen rearrangement by stabilising a particular transition state, and that the associated activation barrier can be estimated. Fidelity is judged relative to that claim, not in the abstract. A model can be an excellent proxy for the reaction chemistry and a poor one for, say, long-timescale conformational dynamics or allostery; validity does not transfer across questions. This section sets out the axes along which fidelity is assessed, states which the present pipeline addresses and with what evidence, and is explicit about which remain open.

## Fidelity is convergent, not singular

There is no single test that certifies a model. Confidence is built from independent lines of evidence that would be unlikely to agree by chance if the model were wrong. The axes below are ordered roughly from necessary-but-weak (structural) to strong-but-demanding (independent observables). Each is necessary; none is sufficient alone.

### 1. Structural fidelity — does the form match?

The model is anchored to experiment at the outset: the starting coordinates derive from the crystal structure (PDB 2CHT), so the fold and quaternary arrangement are experimentally determined rather than assumed. The relevant question is whether that form survives the transition from a static, cryo-cooled crystal to a solvated, room-temperature dynamical system without drifting into an unphysical state.

The evidence in this pipeline: the equilibrated trajectory remains close to the crystallographic reference (protein backbone RMSD stabilising rather than diverging), and the catalytically essential active-site contacts persist throughout the production trajectory. Because the active site is inter-subunit, this check was performed across chains rather than within a single subunit — the substrate in one monomer is coordinated by residues from the adjacent monomer, and measuring same-chain contacts alone would misrepresent the site. The persistent first-shell interactions (the guanidinium of Arg90, the carboxylate of Glu78, and Arg7, with cross-subunit contributions) are the same interactions the experimental and prior computational literature identify as responsible for binding and transition-state stabilisation.

This axis establishes that the *form* is retained. It is necessary but weak on its own: a structure can look correct and still catalyse incorrectly, because function depends on energetics the geometry alone does not reveal.

### 2. Functional fidelity — does it do the right chemistry at the right cost?

The enzyme's function is quantitative: it lowers the activation free energy of a specific reaction by a measurable amount. A faithful model must reproduce not just the connectivity of reactant, transition state, and product, but the *cost* of the transformation.

The evidence here is the located transition state itself. The pipeline yields a converged first-order saddle point with a single dominant imaginary vibrational mode (−313 cm⁻¹) corresponding to the concerted C–O bond cleavage and C–C bond formation of the Claisen rearrangement, connecting a clean chorismate reactant to an exothermic prephenate product. The potential-energy barrier for the studied snapshot is 15.3 kcal/mol.

Two points of honesty about this number. First, it is a *potential-energy* barrier for a single conformational snapshot, not a free energy and not an ensemble average; the appropriate experimental comparator is the activation enthalpy ΔH‡ = 12.7 kcal/mol rather than the activation free energy ΔG‡ = 15.4 kcal/mol, since a single static path does not contain the entropic and conformational-averaging contributions that separate the two. Against ΔH‡ the value is approximately 2.6 kcal/mol high — within the range reported in the prior QM/MM literature for this enzyme at comparable levels of theory, and on the high side of it, consistent with a single high-barrier conformer. Second, the near-coincidence between the computed 15.3 and the experimental ΔG‡ of 15.4 is not evidence of accuracy; it compares different physical quantities and should not be presented as agreement.

### 3. The differential is more trustworthy than the absolute

Absolute barriers carry the full systematic error of the method — the choice of functional, basis set, and embedding. Differences between two barriers computed with the *same* protocol allow much of that error to cancel. For chorismate mutase the most informative single quantity is not the absolute enzymatic barrier but the *difference* between the enzyme-catalysed and the uncatalysed solution-phase reaction, because that difference is the enzyme's actual catalytic effect and corresponds to the experimentally measured rate acceleration.

This is the strongest functional validation available for this system, and in the present pipeline it is still to be done: the solution-phase reference reaction has not yet been computed with the same protocol. Reproducing the enzyme-minus-solution barrier lowering — rather than either absolute barrier — would be the most direct demonstration that the model captures the catalytic function rather than merely landing on a plausible energy. This is flagged as a priority for completing the validation argument.

### 4. Internal consistency — are we reading the model correctly?

Distinct from whether the model matches the enzyme is whether the reported number is a genuine property of the model or an artefact of one algorithm. This is established by locating the same transition state through independent methods. Here the restrained adiabatic scan (barrier maximum ≈ 16.6 kcal/mol), the nudged-elastic-band climbing image (≈ 15.2 kcal/mol), and the reduced-active-region saddle optimisation (15.3 kcal/mol) agree to within roughly 1 kcal/mol. The imaginary frequency of the converged saddle (−313 cm⁻¹) is close to the value obtained in independent prior work on the same system.

This agreement does not prove the model matches the enzyme — three methods can consistently read the same flawed model — but it removes method choice as a source of doubt and confirms the barrier is a real feature of the constructed system. It is internal validity, a prerequisite for taking any external comparison seriously.

### 5. Independent observables — the demanding test

The most persuasive validation is reproducing an observable the model was not constructed to fit. For this reaction the natural candidate is the kinetic isotope effect, which is a direct experimental probe of transition-state structure: a model that reproduces the measured KIE is functioning like the enzyme in a way that is difficult to achieve by coincidence, because the KIE reports on the specific geometry and force constants at the saddle. Active-site mutant effects on the barrier are a second such probe. Neither is addressed in the present work, and both are noted as routes to stronger validation than any absolute-barrier comparison can provide.

### 6. Ensemble fidelity — matching the statistical nature of the enzyme

The living enzyme is not one structure but a Boltzmann-weighted ensemble of conformations, and the observed rate reflects an average over that ensemble. A single snapshot's barrier is one draw from a distribution; the 15.3 kcal/mol value should be understood in those terms. Fidelity to the enzyme in a statistical-mechanical sense requires that the distribution of barriers over representative conformations, properly weighted, reproduces the observable rate. Multiple reactant snapshots have been prepared for this purpose, but the ensemble of barriers has not yet been computed. Completing it serves two functions: it reduces the influence of any single unrepresentative conformer (the present snapshot appears to sit at the high end), and it aligns the calculation with the ensemble nature of the real system rather than treating one geometry as definitive.

## What the present pipeline establishes, and what remains

Addressed, with evidence:
- Structural fidelity: crystal-anchored fold, dynamically stable, with persistent inter-subunit active-site contacts.
- A verified transition state: converged saddle, single reaction-coordinate imaginary mode, connecting clean reactant and exothermic product.
- A barrier within the published range for the system, correctly compared against ΔH‡ rather than ΔG‡.
- Internal consistency: three independent TS-location methods in agreement.

Open, and required for a complete fidelity argument:
- The catalytic differential: the same protocol applied to the uncatalysed solution reaction, compared against the experimental rate acceleration. This is the single most valuable outstanding step.
- The conformational ensemble: a distribution of barriers over representative snapshots, weighted to compare against the measured rate.
- An independent observable: a kinetic isotope effect or mutant-effect calculation as a direct, unfitted test of transition-state fidelity.
- Higher-level energetics: correlated single-point corrections, since the density functional used is known to underestimate this barrier and the absolute value should be treated as provisional.

## The epistemic bottom line

Certainty is not available, and claiming it would be a methodological error. What is available is a convergent argument: a model whose form is anchored to experiment and stable under dynamics, whose transition state is verified and method-independent, and whose barrier falls in the expected range, is a defensible proxy *for the specific question of the catalytic mechanism and barrier*. The strength of that claim scales with how many independent axes agree, which is why the outstanding items — above all the enzyme-versus-solution differential and the ensemble — matter: they convert a plausible single-snapshot result into a validated statement about the enzyme's function. The framing throughout should be that the model is faithful enough to support the mechanistic claim being made, with the fidelity of that claim resting on the weight of convergent evidence rather than on any single number matching experiment.
