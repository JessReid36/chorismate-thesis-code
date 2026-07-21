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
