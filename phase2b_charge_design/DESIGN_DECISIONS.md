# phase2b design decisions

Running record of the scientific/structural decisions behind phase2b, so the reasoning
can be reconstructed later. Newest at top.

## Two-tier optimiser structure (certified screen + relaxed validation)
The scientific question is: WHICH combination of external charges recreates the reaction
on the bare substrate, R -> P. "Recreates the reaction" is a WHOLE-PATH property (a design
that only crushes the single-point TS can leave another path image as the true bottleneck --
the source of the K>=3 negative-barrier artefact), so it cannot be scored on Dv = V_TS - V_R
alone. The enzyme charge comparison is an interesting side-check, NOT the target.

This forces a two-tier structure:

  Tier 1  05_optimise      CERTIFIED frozen-path screen (THE NOVELTY).
                           Certified MILP + linear whole-path min-max over a FIXED
                           consistent-region NEB. Provable branch-and-bound gap.
                           Generates + RANKS candidate charge arrangements. Fast, linear.

  Tier 2  08_relaxed_valid RELAXED-path validation (NOT claimed as novel).
                           For the top ~3-5 ranked candidates only: re-optimise the path/TS
                           UNDER the designed charge field to get the real catalytic effect,
                           confirm a connected R->P path, and allow the TS to shift.
                           Standard (GOCAT-like) methods; inherits the Route-1 / surrogate /
                           charge-attenuation machinery (bare point charges collapse the path).

## Where the novelty lives (PhD defensibility)
Entirely in Tier 1: first CERTIFIED (provable-gap) discrete optimiser for electrostatic
catalyst-field design, and specifically the LINEAR whole-path min-max objective that keeps
the certificate across the whole reaction path. Contrast:
  - GOCAT (Dittner & Hartke, JCTC 2018 14:3547; JCP 2020 152:114106; Behrens & Hartke,
    Top Catal 2022 65:281): heuristic genetic algorithm, NONLINEAR summed-barrier objective,
    NO optimality certificate.
  - Sokalski/Beker (JCTC 2020 16:3420): rotamer scan.
  - Head-Gordon (Nat Catal 2018 1:649): heuristic / directed-evolution field optimisation.
Nobody certifies the search, and nobody has a certifiable multi-image path objective.
Tier 2 is deliberately NOT novel -> keeps the contribution unambiguous.

## GOCAT findings that shaped this (why the fork is real)
GOCAT 2018 kept the path FIXED and reported this gave "only small catalytic effects" and was
"unrealistic". GOCAT 2020 added on-the-fly path re-optimisation -> "excursions ... all the way
to drastic mechanistic changes", incl. a field-dependent concerted -> zwitterionic two-step
transition in their Diels-Alder. Consequences for us:
  1. The TS can shift (even mechanistically) under a strong field -> Tier 2 must allow it.
     The chorismate Claisen is pericyclic (same family) with concerted-vs-polar tension ->
     plausibly susceptible.
  2. Relaxing the path per candidate is FEASIBLE (GOCAT did it) but is exactly what makes the
     objective nonlinear -> why GOCAT needs a GA and cannot certify. Doing it for EVERY
     candidate = becoming GOCAT. So we do it only for the top-ranked few (Tier 2).

## Load-bearing check: ranking preservation
Because Tier 1 is a frozen-path proxy and GOCAT warns frozen-path underestimates real catalysis,
the certificate only MATTERS if the screen's RANKING survives Tier-2 relaxation. Demonstrating
that the frozen-path optimum faithfully predicts which candidates recreate the reaction is the
headline validation of the whole approach, not an optional extra.

## Tier-2 method: three-level funnel (cheap geometry, DFT energy)
Tier 2 stays validation-grade but is made cheaper WITHOUT a fully semiempirical step, which
would corrupt ranking preservation for this dianionic/pericyclic/field-sensitive system.
  Tier 1   certified MILP on the B3LYP-D3BJ/def2-SVP Dv grid (unchanged, DFT-quality)
  Tier 2a  r2SCAN-3c relaxed-path optimisation under the designed field (cheap GEOMETRY;
           collapse mitigations on). Geometry + connectivity + TS-shift are the deliverables;
           r2SCAN-3c energies only coarsely shortlist.
  Tier 2b  B3LYP-D3BJ single point (diffuse-augmented basis, e.g. def2-SVPD/ma-def2-SVP) on the
           top 2-3 relaxed paths -> FINAL barrier + ranking. This is the word on ranking.
Rationale: r2SCAN-3c geometry is DFT-quality (on par with M06-2X-D3/TZP); all cheap methods
degrade on barriers, worst for GFN2-xTB on pericyclic (BHPERI) + anions + field response, so
GFN2-xTB is at most a geometry pre-filter, never in the ranking. def2-SVP lacks diffuse
functions for a -2 species -> confirmation uses a diffuse-augmented basis.
Calibration go/no-go (run on first real system, all ~5 shortlisted through r2SCAN-3c AND DFT):
  keep cheap triage      if Spearman rho >= 0.9 AND barrier-change error spread <= 1 kcal/mol
  demote to geometry-only if energies reorder (rho < 0.9) but geometries/TS agree -> DFT SP on all
  abandon cheap level    if relaxed geometry/connectivity disagree (spurious collapse/TS) -> full DFT
The SPREAD (candidate-dependence), not the mean offset, is the decision metric.
Source: phase2b research report "Cheap-Geometry, DFT-Energy Funnel for Tier 2".

## Certified optimality is optimality of the PROXY, not of catalysis (critical framing)
The Tier-1 certificate (branch-and-bound gap 0.000) proves optimality w.r.t. the OBJECTIVE we
write down -- currently "minimise Sum q*Dv", i.e. lower the frozen-single-point barrier.
"Provably best at minimising the Dv-barrier" is NOT "provably best at catalysing the reaction".
Conflating the two would be the most attackable claim in the thesis. The certificate guarantees
a PROXY optimum, not a catalytic optimum.

Three gaps between the proxy and real catalysis:
1. Frozen-density gap: Dv uses a fixed electron density. Real charges repolarise the substrate;
   past ~1-3 charges the frozen estimate diverges (max-lowering barrier -> physically impossible
   negatives -- the giveaway).
2. Frozen-path gap (the big one): catalysis = lowering the whole R->TS->P barrier with a still-
   connected path, and charges can MOVE the TS (GOCAT). Minimising Dv at the fixed TS point can
   "win" by crushing that point while another path point becomes the real bottleneck -> a
   Dv-optimal design can be catalytically WORSE than a Dv-suboptimal one. The certificate cannot
   see this by construction.
3. Over-optimisation gap: because more charges trivially lower the proxy, "optimal" under
   max-lowering actively selects the unphysical, spill-out-prone, path-distorting arrangements --
   it optimises hardest in the direction that departs from catalysis.

Consequence: the "more charges -> lower barrier" trend under max-lowering is monotone BY
CONSTRUCTION (a non-helping charge is placed as ~zero contribution), so it is NOT a chemical
finding. The max-lowering column is a DIAGNOSTIC only: (a) the certified ceiling at each K,
(b) the K where it goes unphysical (~K=4 on the real field) = the frozen-objective validity limit
and the concrete argument for Tier-2.

## Two-level optimality claim (how the thesis must carry it)
- Tier 1 = certified optimality of the SCREEN: not "this is the best catalyst" but "this is
  provably the best candidate under our physically-motivated proxy objective, with a certificate".
  Novel precisely because no prior method certifies its screen (GOCAT's GA cannot even claim its
  proxy-optimum is a proxy-optimum).
- Catalysis is decided at Tier 2 and is NOT claimed certified. Dv-optimal + near-optimal candidates
  are carried into relaxed-path validation (path relaxed, TS free); the real barrier there judges
  catalysis. The certificate's job is to choose WHICH candidates get the expensive validation --
  provably rather than heuristically.
Why we downweight max-lowering: it is the proxy that departs MOST from catalysis (gap 3). The
distributed and whole-path objectives write a proxy that CORRELATES BETTER with catalysis (spread
load -> less spill-out/path distortion; whole-path min-max -> cannot cheat by crushing one point).
Still proxies, but better ones; a certificate over a better proxy is worth more.

## Ranking preservation is THE load-bearing experiment (not a footnote)
The certificate is valuable only if the certified Tier-1 ORDERING predicts the Tier-2 catalytic
ordering. Ranking preservation tests exactly that. Holds -> the certificate provably finds the
candidates that turn out to catalyse. Fails -> we have certified a proxy that does not track the
goal, and the objective must be improved. This is where "optimal at lowering the barrier" is
converted into (or fails to convert into) "optimal at promoting catalysis".

## Contribution framing (one line)
Not "we found the optimal catalyst" but "we built the first CERTIFIED screen for electrostatic
catalyst design, and we test honestly whether certified-optimal-on-the-proxy tracks actually-
catalytic". Keeping these two claims rigorously separate is what makes the thesis defensible.
