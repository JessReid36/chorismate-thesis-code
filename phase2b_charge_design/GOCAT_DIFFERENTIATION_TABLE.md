# GOCAT vs this thesis — side-by-side differentiation (the key defensive artifact)

Grounded in the GOCAT primary sources (Dittner 2019 dissertation = gocat_2; Behrens 2024 = gocat_1),
verified by direct quotation, not secondary summary. This table pre-empts the single most dangerous
near-miss (GOCAT) on every axis where the novelty check flagged a threat. Line refs are to the uploaded
dissertation text files.

## The one-glance table

| Axis | GOCAT (Dittner/Hartke; Behrens/Hartke) | This thesis | Differentiator |
|------|----------------------------------------|-------------|----------------|
| **Optimiser** | Evolutionary/genetic algorithm (ogolem suite) | Mixed-integer linear program (HiGHS) | Heuristic vs exact |
| **Guarantee** | NONE. States complete enumeration + deterministic global optimization are "impossible"; uses "unbiased metaheuristical optimization algorithms that avoid complete enumeration" (gocat_2 L52-54) | Provable branch-and-bound optimality gap = 0.000 at every solve | No certificate vs certified global optimum |
| **Objective** | NONLINEAR: NEB/AutoNEB reaction-barrier fitness; Behrens Eq 3.1 weights the largest barrier 10x + sum of the rest (= 11·dE for a single-barrier candidate) | LINEAR: whole-path differential electrostatic potential min-max, Dv = V_TS - V_R (Sokalski catalytic field) | Nonlinear barrier vs linear certifiable proxy |
| **Why not certified** | Nonlinear objective + vast combinatorial space -> global opt "impossible" (their word) -> GA is the only option | Linear objective + discrete sites -> MILP is exactly solvable to gap 0 | The linear objective is what MAKES certification possible |
| **Building blocks** | Point charges only (2018/2020); later REAL molecular / QM-MM embeddings (Behrens 2022) | Fixed optimised ±1 charges PLUS a separately-designed NEUTRAL steric/dipolar cradle | Two independently-optimised layers (electrostatics vs geometry) vs charges-only or real-molecule embedding |
| **Charge model** | q_i in [-1,+1]; min-distance r_min; net-charge constraint Sum q_i = const (e.g. 0 neutrality) (gocat_2 L1341-1345,1755); vdW embedding surface | Discrete ±1 on a vdW-shell grid; net-FREE (optimal at each K); surrogate-sized excluded volume | SHARED realism ingredients - NOT claimed novel; the LINEAR CERTIFIABLE ENCODING is |
| **Geometric-ceiling finding** | None | Certified across K2-K10: point charges alone cannot hold the dianion in near-attack; too few under-constrain, too many over-polarise/rupture | Only a CERTIFIED optimiser can prove "no arrangement of N charges suffices" |
| **System framing** | Bare reacting system, enzyme-free (shared) | Bare chorismate dianion, enzyme-free; CM as benchmark not template | Shared framing; not a differentiator |

## Verbatim primary-source anchors (for the thesis to quote)
- GA / no certificate: "chemical space is vast even for small compounds, which makes complete enumerations
  and deterministic global optimization impossible. Thus, additional ingredients are necessary such as ...
  leveraging unbiased metaheuristical optimization algorithms" (gocat_2 L52-54). "evolutionary algorithms
  (EAs) are harnessed as implemented in our global optimization suite ... ogolem" (gocat_2 L62-64).
- Charge model: q_i in [-1,1], net-charge "qi = const. = 0 (neutrality) are to be conserved" (gocat_2
  L1345); vdW embedding surface, charge pushed back onto surface if inside (gocat_2 L1755-1764).
- Objective: nonlinear barrier fitness (Behrens Eq 3.1, 11·dE weighting; gocat_1 fitness section).

## Independent corroboration (a SECOND group also finds global opt of point charges impractical)
Weymuth & Reiher, "Gradient-Driven Molecule Construction" (Int. J. Quantum Chem. 2014, 114, 838-850;
arXiv:1401.1491): inverse-designs a point-charge "jacket potential" around a fragment, but targets
VANISHING GEOMETRY GRADIENTS (stability), NOT a barrier/TS objective, and states directly: "A global
optimization is extremely time-consuming already for very small search spaces. The global optimization of
a large set of point charges can thus be expected not to be possible in a straightforward fashion." -> a
second independent primary source treating certified/global point-charge optimisation as intractable,
which STRENGTHENS the certification-novelty claim. (Report mis-cited as "Krieg & Reiher"; correct is
Weymuth & Reiher 2014. Different objective (geometry stability) + heuristic/local -> does NOT threaten
Claim 1.)

## How to deploy this in the thesis
1. Claim 1 (certified optimiser): this table + the two verbatim GOCAT quotes + Weymuth&Reiher = the whole
   defence. Novelty is the CERTIFICATE; point-charge design is GOCAT's, the Sokalski objective is
   Sokalski's. GOCAT itself says global opt is "impossible" for their formulation -> the linear objective
   is precisely what lets us certify.
2. Claim 2 (two-layer): the "Building blocks" row - GOCAT is charges-only then real-molecule embeddings;
   ours is charges + an independently-designed NEUTRAL cradle. Foreground the separation + bare substrate.
3. Claim 3 (ceiling): the "Geometric-ceiling finding" row - only a certified optimiser can establish it;
   GOCAT's GA cannot prove non-existence of a sufficient arrangement.
4. Do NOT over-claim the realism ingredients (bounded ±1, net constraint, excluded volume, vdW surface):
   they are SHARED with GOCAT (charge-model row). Novelty is encoding them as linear MILP constraints that
   preserve the certificate - which GOCAT's nonlinear GA cannot do.
