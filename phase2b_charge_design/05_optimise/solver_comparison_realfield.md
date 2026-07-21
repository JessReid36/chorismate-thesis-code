# CBC vs HiGHS on the REAL Dv field (production, 331 pts, menu +/-1, K=1-10)

Both solvers on the actual dv_grid.tsv. Max-lowering closes at gap 0.000 for both (simpler problem).
The DISTRIBUTED min-max tail is where they diverge.

## HiGHS (highspy) - canonical
K  max-lowering(ddE/barrier/gap)   distributed(max|contrib|/total/barrier/gap)
 1   -5.77  +11.70  0.000            5.27  -5.27  +12.20  0.000
 2  -11.04   +6.43  0.000            2.65  -5.29  +12.18  0.000
 3  -16.07   +1.40  0.000            1.83  -5.30  +12.17  0.000
 4  -20.88   -3.41  0.000            1.34  -5.28  +12.19  0.000
 5  -25.54   -8.07  0.000            1.08  -5.27  +12.20  0.000
 6  -29.94  -12.47  0.000            0.90  -5.30  +12.17  0.000
 7  -34.15  -16.68  0.000            0.77  -5.27  +12.20  0.000
 8  -38.05  -20.58  0.000            0.69  -5.28  +12.19  0.000
 9  -41.89  -24.42  0.000            0.60  -5.27  +12.20  0.000
10  -45.00  -27.53  0.000            0.55  -5.28  +12.19  0.000
-> all rows certified gap 0.000, full sweep in ~8 s.

## CBC (python-mip, tmax=600 s) - distributed tail
K   distributed max|contrib| / gap
 6   0.90 / 0.004     (near-certified)
 7   0.80 / 0.046     (NOT certified AND SUBOPTIMAL: HiGHS optimum is 0.77)
 8   0.69 / 0.002     (near-certified)
 9   0.60 / 0.000
10   0.55 / 0.000
(K=1-5 CBC = 0.000, matching HiGHS.)

## The finding
CBC leaves K=6-8 distributed uncertified even at 600 s, and at K=7 its incumbent (0.80) is not
just unproven but WORSE than the true optimum (HiGHS 0.77, certified). So an open gap can conceal
a genuinely suboptimal design, not merely an unproven-correct one. Since Tier-1's contribution is
CERTIFIED optimality, this is decisive: HiGHS is used for all reported certificates; CBC is retained
only as an independent correctness cross-check (agrees where it closes). See SOLVER_CHOICE.md.
