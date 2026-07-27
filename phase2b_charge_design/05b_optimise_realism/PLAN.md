# 05b — Realism-augmented certified Tier-1 MILP: PLAN

Anchor document so the switch to a realism-augmented certified optimiser does not lose the thread.
References:
  - Research report: "Certified MILP Enzyme Design: Adding Active-Site Realism Without Losing the
    Optimality Gap" (the recommendations this plan implements).
  - phase2b_charge_design/GOCAT_PRIMARY_SOURCE_COMPARISON.md (Dittner 2019 + Behrens 2024, verified).
  - Motivating evidence (committed, results repo): relaxation_attempts/ash_guarded/
      production_batch/ENDPOINTS_PROVENANCE.md  (K2 reactant C1-C6 5.00 vs 3.12 bare),
      production_batch/ENDPOINT_VALIDITY_ASSESSMENT.md (endpoints valid; field pre-organises AWAY from TS),
      barrier_k2_scan/K2_SCAN_FINDING.md (1D scan +36.5, dissociative ridge C1-C6 6.08; suggestive only),
      barrier_k2_neb/ (NEB band corroborates dissociative character).

## Why the switch (the numbers backing it)
Certified Tier-1 optimises the FROZEN single-point proxy Dv = V_TS - V_R. Under relaxed validation the
certified-optimal K2 design (one +1, one -1) RELAXES the substrate's forming bond C1-C6 from 3.12 A (bare
substrate) to 5.00 A - pre-organising the substrate AWAY from the near-attack geometry (reactive window
2.6-3.5 A; TS ~2.5 A). Both barrier probes point the same way: the 1D O3-C4 scan crests at +36.5 kcal/mol
with C1-C6 = 6.08 A (dissociative, not the concerted TS; vs +17.47 bare baseline), and the NEB band shows a
dissociative/slack path. This is the concrete instance of the standing caveat "certified optimality is
optimality of the PROXY, not of catalysis" (DESIGN_DECISIONS.md). The fix is to make the charge
representation a credible enzyme-site model AND teach the proxy to reward reactant pre-organisation -
WITHOUT abandoning the certificate.

## GOCAT-verified positioning (do NOT over-claim)
GOCAT ALREADY has: bounded +-1 charges, minimum inter-charge distance (r_min=1A), net-charge/neutrality
constraint, vdW embedding surface, dipolar/H-bond/LJ representations (Dittner 2019 Eqs 2.15-2.19; Behrens
2024). So the realism INGREDIENTS are shared, NOT our novelty. Our novelty is encoding them as
LINEAR/MILP constraints under a LINEAR whole-path Dv objective solved to a PROVABLE gap=0.000 - which
GOCAT's nonlinear 11*dE NEB-barrier objective + ogolem GA cannot certify. The report's realism upgrades are
therefore adopted, but framed as "certifiable linear encoding of GOCAT's own charge model", not as new
physics.

## The four steps (implement in order; each references the report)

### Step 0 — Define the reactant-distorting mode  [prerequisite]
Report basis: section B(a). Compute the direction that carries C1-C6 3.12 -> 5.0 A from the BARE reactant.
Candidate definitions (decide before coding): (i) raw C1-C6 internal coordinate [lean: cleanest linear
statement]; (ii) low-frequency normal mode projecting onto C1-C6 opening; (iii) empirical difference vector
bare-reactant -> distorted-K2-reactant [cross-check]. Output: preorg/distorting_mode.dat + per-candidate-site
linear coefficients a_s = [V_s(C1) - V_s(C6)] (the differential potential each candidate charge site induces
ACROSS the C1<->C6 pair; linear in q by superposition).

### Step 1 — Realism constraints as LINEAR MILP rows  [Tier 1, certificate-preserving]
Report basis: Stage 1 / section A. Add to the sweep, each verified linear: bounded +-1 (present);
net-neutrality Sum q = const (linear; GOCAT-shared); big-M excluded-volume / min-distance (GOCAT r_min, as
MILP no-overlap rows); CPCM(eps=4)-screened per-site potential kernels (linear superposition); optional
net-neutral dipolar-surrogate library as one-hot placements. CONFIRM HiGHS gap=0.000 unchanged. State in
README: ingredients shared with GOCAT, the certifiable encoding is ours.

### Step 2 — Preorganisation objective term  [Tier 1, THE fix]
Report basis: Stage 2 / section C. Augment the linear whole-path Dv min-max with a linear preorganisation
term using Step 0: hard cap |V(C1)-V(C6)| <= tau AND/OR penalty mu * (distorting-mode projection). Re-solve;
CONFIRM gap=0.000 still holds. This is the genuine objective-level novelty (a preorganisation-penalised
CERTIFIABLE proxy - no prior certified catalytic-field objective couples path Dv to a preorganisation term).

### Step 3 — Ranking-preservation check  [make-or-break]
Report basis: Stage 3 + DESIGN_DECISIONS load-bearing test. Re-run certified sweep with realism+preorg; take
new top 3-5; run Tier-2 relaxed validation (08_relaxed_validation, the GOCAT-like nonlinear step). Metrics:
does relaxed C1-C6 stay 2.6-3.5 A? does Spearman rho(Tier-1 proxy, Tier-2 effect) improve? is the K2-type
inversion gone? The certificate only MATTERS if ranking survives.

### Step 4 — Enzyme yardstick  [interpretation, not fitting]
Report basis: Stage 4 / section D. Overlay realism-augmented certified-optimal design on BsCM
(Arg90/Arg7/Glu78); report convergence OR divergence honestly (07_validate machinery). Benchmarks: enzyme
DTSS ~7.3 kcal/mol (Claeyssens 2011), optimal catalytic field V_TS-V_R (Szefczyk/Sokalski 2004).

## Parallel, non-blocking
K2 Cartesian 2D scan / NEB barrier = the "bare K2, no realism" anchor (the BEFORE to 05b's AFTER). Finish
when convenient; does not gate the formulation work.

## Directory / branch
Branch: tier1-realism (off main). New dir: 05b_optimise_realism/ (parallel to 05_optimise/, which stays
untouched and reproducible - the DIFF between them is the contribution). preorg/ holds the distorting-mode
vector + per-site coefficients. Merge to main only when Step 3 ranking-preservation passes.
