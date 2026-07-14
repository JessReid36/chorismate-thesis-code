Stage 3 - charge-selection optimizer.

The frozen-density barrier change from external point charges is linear in the charges,
dDE = sum_i q_i * Dv(r_i), where Dv is the difference potential (Stage 1) evaluated at the
candidate grid points (Stage 2). Dv was computed on the 326-point grid by orca_vpot from the
reactant and TS densities (3.1); of 326 candidate positions, 194 are stabilising (Dv < 0).

Two selection problems are solved on this objective. (a) A charge-BUDGETED continuous solve:
minimise dDE subject to a total-charge budget sum|q_i| <= Q, net-neutrality, and |q_i| <= 1.
The barrier is unbounded without a budget (the naive box-constrained optimum pins every charge
to its bound), so the L1 budget is what makes the continuous problem well-posed; it also induces
sparsity, so the solver selects where a limited budget is most effective. Sweeping Q gives the
barrier-lowering-vs-charge curve (the continuous ceiling). (b) A discrete design: charges in
{-1,0,+1}, at most K active, optionally net-neutral, with a physical minimum-separation between
active charges (2.5 A, the floor below which real charged moieties overlap; stable to 3.5 A,
salt-bridge scale). The discrete objective is linear, so this is a mixed-integer linear program
solved to provable optimality with a certificate (CBC) - a global-optimality guarantee that the
heuristic global search of prior charge-optimisation methods (e.g. GOCAT's evolutionary algorithm)
does not provide.

For the placeholder geometry the net-neutral discrete designs match the continuous ceiling to
within a few percent, and the optimizer independently and repeatably places a +1 cation at the
grid point of most-negative Dv (the developing-negative-charge region on the breaking ether
oxygen) - the from-scratch analogue of the catalytic Arg90. This placement is stable across the
continuous and discrete solves and across separation values. A net-charge-free variant scores
marginally higher at low K by placing anions on the largest-positive-Dv points on the innermost
shell, but this advantage is concentrated in the regime where the frozen-density model is least
reliable and collapses to the net-neutral solution once realistic separation is enforced at
larger K; it is carried as a sensitivity case for the Stage 4 polarization correction to test.
