# Step 0 — reactant-distorting mode definition (record)

PLAN Step 0 / report section B(a). Defines the linear per-site coefficients the preorganisation term
(Step 2) uses to penalise fields that drive C1-C6 apart (the observed 3.12 -> 5.0 A distortion).

## Chosen coefficient (PRIMARY): a_pot,s = V_s(C1) - V_s(C6)
The electrostatic potential a unit +1 charge at candidate site s induces across the C1<->C6 pair
(k=1, Angstrom). Linear in the source charge by Coulomb superposition, so a design's total distorting
drive = sum_s q_s * a_pot,s is LINEAR in the MILP placement variables -> encodable as a hard cap or
penalty WITHOUT breaking the certificate. Same differential-potential construction as Dv = V_TS - V_R,
but across the two reacting carbons rather than two path states -> consistent with the pipeline's
existing Sokalski-field language.

## Cross-check coefficient (SECONDARY, committed alongside): a_force,s = [F_s(C6) - F_s(C1)] . u_hat
Mechanical force-imbalance along the C1->C6 axis. Measures geometric push-apart directly (a gradient
quantity), one step removed from the potential the pipeline uses.

## Why a_pot is primary (decision record)
Considered three mode definitions (PLAN Step 0): (i) raw C1-C6 internal coordinate via potential
difference [CHOSEN]; (ii) low-freq normal mode; (iii) empirical bare->distorted difference vector.
Chose (i)/a_pot because:
  1. VALIDATED on ground truth: the K2 max-lowering design is KNOWN to distort C1-C6 apart (3.12->5.0).
     Scoring K2's two charges gives sum q*a_pot = +0.105 (>0 = distorting) -> a_pot correctly flags the
     one design we can check. (a_force also flags it, +0.014 -> corroborates the mode; see below.)
  2. Consistent with the differential-potential objects the MILP already uses (linear in the same V's).
  3. Cleanest linear statement of "do not separate charge across the reacting carbons".

## Honest caveat (do not hide)
a_pot and a_force AGREE on the aggregate verdict for K2 (both >0) but disagree in SIGN at the individual
site level (~50% agreement) - they measure different mechanisms (charge separation vs mechanical push)
and diverge per-site. We adopt a_pot, commit a_force as an independent cross-check, and let Step 3
(ranking preservation under REAL relaxed validation) be the final arbiter of whether the a_pot-based
penalty actually removes the distortion. If Step 3 shows a_pot fails to fix C1-C6 but a_force would,
revisit this choice - the decision is falsifiable at Step 3.

## Files
compute_distorting_coeffs.py : generator (reads reactant.xyz + dv_grid.tsv).
distorting_coeffs.tsv        : per-site idx, shell, a_pot, a_pot_norm (max|.|=1), a_force (xcheck).
Spread: 331 sites; |a_pot_norm|>0.3 on 39 sites, >0.5 on 9 - the term discriminates (not ~0 everywhere).
Strongest coefficients sit on the innermost shell (2.0 A), as expected (closest sites distort most).

## Consumed by
Step 2 (preorg objective term): hard cap |sum_s q_s a_pot,s| <= tau  AND/OR penalty mu*|sum_s q_s a_pot,s|,
added to the linear whole-path Dv min-max in 05b's augmented tier1 sweep. Both keep the MILP certifiable.
