# Step 2 — preorganisation penalty sweep (findings)

PLAN Step 2 / report section C. Adds a LINEAR preorganisation penalty to the certified Tier-1 max-lowering
objective, using the Step-0 coefficients a_pot,i = V_i(C1) - V_i(C6). Distorting drive of a design
D = sum_i q_i * a_pot,i (linear); |D| linearised (1 aux var + 2 rows). Objective: minimise
[sum q*Dv*627.509] + mu*|D|. NET-FREE (optimal +-1 at each K; no net-charge constraint). Every mu solved
to HiGHS gap = 0.000 -> each point on the trade-off is individually CERTIFIED. Script: tier1_sweep_preorg.py.

## The trade-off (certified at every mu; barriers are the FROZEN Dv PROXY, not real catalysis)
As mu increases, the certified optimum moves OFF the C1-C6-distorting sites (D -> ~0 or negative) and the
proxy barrier rises. Distortion at mu=0 grows with K (more charges = more distorting drive):
  K1: D 0.058 (mu0, barrier +11.70) -> ~0.002 (barrier +13.58)
  K2: D 0.105 (mu0, +6.43)          -> 0.0006 (+11.15) ; original distorting K2 = the mu=0 design (sites 150/249)
  K3: D 0.158 (mu0, +1.40)          -> ~0      (+7.32)
  K4: D 0.213 (mu0, -3.41)          -> ~0      (+3.68)
Notably, HIGHER-K designs retain MORE barrier-lowering at low distortion (K4 barrier +3.68 at D~0 vs K2
+11.15) -> more charges let the design preorganise AND lower the barrier: the quantitative form of the
"higher K = fuller charge cavity = more enzyme-like" hypothesis. Emergent, not imposed.

## HONEST caveats (load-bearing - do not over-read the trade-off)
1. Barriers here are the FROZEN Dv PROXY. The mu=0 K2 "+6.43" is the SAME proxy that relaxed to a
   dissociative +36.5 kcal/mol (barrier_k2_scan). So NO barrier number in this sweep is trustworthy as
   catalysis. This is a SCREEN output, not a catalysis result.
2. Step 3 (relaxed validation) is the ARBITER: take a low-D design (e.g. K2 mu=100, sites 105/121) through
   Tier-2 and check whether relaxed C1-C6 actually stays in 2.6-3.5 A. If a low-D design STILL relaxes open,
   the a_pot mode was the wrong choice -> revisit with a_force (STEP0_MODE_DEFINITION falsifiable arbiter).
3. Odd motifs to scrutinise, NOT trust: K2 mu>=100 flips to TWO +1 (no dipole); K4 keeps a DISTORTING -1 at
   site 9 (a=+0.117) because its Dv payoff outweighs the penalty. May be real, may be proxy-gaming.

## Certificate status
Retained. At fixed mu the augmented objective is linear -> gap 0.000. The preorganisation-penalised
CERTIFIABLE proxy is the objective-level novelty (no prior certified catalytic-field objective couples the
path/Dv screen to a preorganisation term). GOCAT could not certify this (nonlinear 11*dE + GA).

## Next (Step 3)
Pick ~2-3 designs spanning the trade-off (mu=0 distorting, mu=100 preorganised, one higher-K) -> relaxed
validation in 08_relaxed_validation -> does the preorg penalty actually keep C1-C6 compact under relaxation,
and does Tier-1 ranking survive? That decides whether a_pot is the right mode and whether the penalty works.
