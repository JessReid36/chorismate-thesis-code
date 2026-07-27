# GOCAT primary-source comparison (verified against Dittner 2019 + Behrens 2024 dissertations)

Verified directly from the two GOCAT dissertations (not secondary summaries) to make the novelty
differentiation airtight and to avoid over-claiming shared ingredients as novel.

## What GOCAT ACTUALLY is (Dittner 2019, Behrens 2024)
- OBJECTIVE (nonlinear): Behrens Eq 3.1 - fitness weights the largest barrier 10x plus the sum of all
  other barriers (= 11*dE for a single-barrier candidate), computed from NEB path energies. A max-plus-sum
  over NEB-computed barriers -> NONLINEAR in charge positions/values. Plus penalty terms.
- SEARCH (heuristic, NO certificate): Dittner states complete enumeration + deterministic global
  optimization are "impossible" -> uses evolutionary algorithms (EAs) in the ogolem GA suite. No bound,
  no optimality gap.
- CHARGE MODEL (Dittner Sec 2, Eqs 2.15-2.19): N_Ch partial point charges, q_i in [-1,+1] e, minimum
  inter-charge distance r_min = 1.0 A, and a NET-CHARGE constraint Sum q_i = const (e.g. 0 for neutrality).
  Confined to a vdW surface (or sphere/ellipsoid) of all path frames. Optimized over CONTINUOUS q and
  position by GA. H-bond centres represented "by suitable combinations of partial charges and vdW centres";
  steric via LJ / rare-gas pseudo-potentials.
- PATH: 2018/static = FIXED path ("only small catalytic effects", "unrealistic"); 2020/Behrens = on-the-fly
  NEB/AutoNEB re-optimization per candidate -> the nonlinearity that forces the GA.

## CRITICAL correction to our positioning (do NOT over-claim)
GOCAT ALREADY HAS: bounded charges (+-1 e), minimum-distance/excluded-volume-like constraint, net-charge
(neutrality) constraint, vdW embedding surface, and molecular/LJ/H-bond representations. So "adding realism"
(bounded charges, net-neutrality, excluded volume, vdW surface, dipolar/H-bond groups) is NOT novel vs GOCAT
- they have all of it. We must NOT claim these as our contribution.

## What IS genuinely novel (unchanged, now sharper and primary-source-grounded)
Our contribution is making the SAME physically-motivated charge model a CERTIFIED DISCRETE problem:
  - DISCRETE candidate sites + charge values (vs GOCAT's continuous q and position),
  - LINEAR whole-path differential-potential (Dv = V_TS - V_R) min-max objective (vs GOCAT's NONLINEAR
    max-plus-sum-of-NEB-barriers),
  - solved to a PROVABLE branch-and-bound gap = 0.000 (vs GOCAT's GA with no bound/gap).
"Nobody certifies the search, and nobody has a certifiable multi-image path objective" - CONFIRMED: GOCAT's
own text says deterministic global optimization is "impossible" for their (nonlinear) formulation. We make
it possible by choosing a LINEAR objective + discrete decisions. The realism ingredients are ENCODED AS
LINEAR/MILP CONSTRAINTS THAT PRESERVE THE CERTIFICATE - that encoding (not the ingredients themselves) is
the novelty, precisely because GOCAT's nonlinear objective CANNOT be certified.

## Consequence for the thesis framing
Frame as: "we adopt GOCAT's physically-motivated charge model (bounded +-1, net-charge-constrained,
surface-confined, dipolar/H-bond-capable) but replace their nonlinear GA-searched barrier objective with a
LINEAR whole-path differential-potential objective admitting a certified global optimum." The K2 reactant-
distortion finding is then the concrete demonstration of why the PROXY (linear Dv) must be augmented with a
linear preorganisation term - and why Tier-2 relaxed validation (the GOCAT-like nonlinear step) is retained
but confined to the top few, so the certificate is never contaminated.
