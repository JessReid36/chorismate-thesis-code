Stage 4 - polarization correction and the frozen-density validation.

The optimization (Stage 3) uses the frozen-density barrier change dDE = sum_i q_i*Dv_i, which
holds the substrate density fixed. Stage 4 tests this against the self-consistent result by
embedding a selected discrete design as external point charges (ORCA %pointcharges, inside the
eps=4 CPCM cavity with the substrate) and running four single points, {reactant, TS} x {bare,
with-charges}, all at B3LYP-D3BJ/def2-SVP. The polarized barrier change is the difference of
differences, [E_TS(pc) - E_R(pc)] - [E_TS(bare) - E_R(bare)]; the charge-charge self-energy is
identical at reactant and TS (fixed charge positions) and cancels, isolating the substrate-charge
interaction plus the density-relaxation (polarization) response.

For the placeholder geometry, two K=2 designs were tested. The net-neutral design (a +1 at the
most-negative-Dv point - the developing-negative-charge region, the enzyme's Arg90 analogue - plus
a counter -1) gives a frozen-density barrier lowering of -9.27 kcal/mol and a self-consistent value
of -8.35, a polarization correction of +0.92 (10.0%). The net-charge-free design (two -1 charges on
the two most-positive-Dv points, total system charge -4) gives -11.28 frozen and -10.10 polarized,
a correction of +1.18 (10.4%). Both SCF sets converged.

Two results follow. First, the frozen-density model is validated for both charge arrangements: the
polarization correction is ~10% in each case, so the linear q.Dv optimization captures about 90% of
the true self-consistent barrier change, and the correction is robust to whether the design is
net-neutral or net-charged. This justifies performing the full optimization on the frozen density.
Second, and contrary to an initial expectation that the net-free advantage would be a frozen-density
artifact, the net-free design outperforms the net-neutral one even after self-consistent relaxation
(-10.10 vs -8.35 kcal/mol). The two corrections being essentially equal shows this is a real effect,
not an artifact: the largest-magnitude values of Dv lie in charge-depletion regions (favouring
anions), and the net-neutral constraint sacrifices stabilisation by forcing a cation onto a
weaker-potential point and spending a charge on balance. The electrostatically optimal designed field
for this reaction is therefore net-charged, and outperforms the enzyme's net-neutral arrangement.
This is reported as a strong preliminary finding pending confirmation on the real isolated-substrate
geometry (the margin currently rests on placeholder Dv values) and subject to the realizability
considerations of a net-charged environment.
