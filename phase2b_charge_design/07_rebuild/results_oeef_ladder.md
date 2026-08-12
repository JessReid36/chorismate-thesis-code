# OEEF field-vs-integrity tradeoff (sparse, all-bond-protected designs)

The sparse OEEF optimiser (all 24 substrate bonds field-capped at FBOND) trades catalytic -z field
against per-bond field (integrity). Field scales ~linearly with FBOND:

  FBOND    avg -z field   charges   ~expected lowering (from finite-field scan)
  0.010    0.0124         4         ~-5.5 kcal/mol
  0.013    0.0156         4         ~-7   kcal/mol   <- currently screening (safe: < pair's 0.0137)
  0.016    0.0189         4         ~-9   kcal/mol   <- above pair ceiling; hold unknown
  0.020    0.0221         6         ~-11  kcal/mol   <- full catalytic field; hold unknown

Fragmentation threshold (bond-field the substrate tolerates) is bounded so far:
  - held:      pair oez0p50 at max bond-field 0.0137  (VALID both endpoints)
  - fragmented: dense 2927-charge design, bond-field uncontrolled (broke 7 bonds)
Threshold is somewhere > 0.0137. Screening the ladder (0.013 -> 0.016 -> 0.020) walks the field up until
something fragments, mapping the exact integrity ceiling AND the maximal catalytic field within it.
This is the quantitative electrostatic-ceiling characterisation (with the field correctly aligned to
Delta-mu, unlike the Sum q.dV proxy).
