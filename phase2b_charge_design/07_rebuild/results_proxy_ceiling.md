# Result: the frozen-geometry Sum q.dV proxy does NOT track the true barrier (ceiling)

## The decisive data point
Design cradleK4q0p30 (cradle + 4 catalytic charges at +-0.3; certified MILP objective Sum q.dV = -0.011,
i.e. the proxy calls it "stabilising"):
  - NEB-TS (IDPP, ActiveRegion, 10 images, ~99 iters): TRUE barrier = +37.6 kcal/mol
  - Relaxed scan (physical region, pre-derailment): +29 kcal/mol climbing
  - Bare reference: +17.47 kcal/mol
The field the proxy called stabilising DOUBLES the true barrier (+37.6 vs 17.5).

## Interpretation (matches Beker & Sokalski 2020 + GOCAT's own experience)
Sum q.dV is the first-order, additive, geometry-FROZEN "extreme approximation" limit of DTSS. It maps
WHERE charges go, not the relaxed-path barrier. A "good" (negative) objective can correspond to a much
HIGHER true barrier -- demonstrated here. Reweighting the proxy (incl. the new varied-magnitude QP)
cannot fix this: the objective doesn't see the relaxed path. GOCAT hit the same wall (their fixed-path
version gave "only small catalytic effects"; they needed on-the-fly path optimisation).

## Consequence
The certified-optimal-of-the-proxy programme has a ceiling for this substrate. Next: a barrier-CORRELATED
objective (field-response / OEEF) or a neutral steric cradle. FIRST measure whether the barrier is
field-sensitive at all (cheap diagnostic) before building a response optimiser.
