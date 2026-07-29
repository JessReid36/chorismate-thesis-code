# Shell control test: does moving charges OUTWARD reduce over-polarisation/fragmentation?

## The question (user)
Our optimiser always clusters charges on the CLOSEST shell (strongest |Dv|) -> harshest field -> fragments.
Would forcing them OUTWARD reduce the harshness? Naive translation is ill-defined (which vector? not
equidistant). Resolution: the optimiser naturally makes charges EQUIDISTANT on whatever shell it's allowed,
so RE-OPTIMISE per shell (no translation) = a clean uniform-distance sweep. Angular placement still optimises;
only radial distance is controlled.

## Designs (K=6, re-optimised on each shell separately, certified gap 0.000)
  K6_shell2 : 6 charges all on 2A shell (surr-substrate 2.21A)  proxy barrier -6.19  HARSH  -> expect fragment
  K6_shell3 : 6 charges all on 3A shell (surr-substrate 3.24A)  proxy barrier -0.04  mid
  K6_shell4 : 6 charges all on 4A shell (surr-substrate 4.37A)  proxy barrier +3.47  GENTLE -> expect intact
Tradeoff already visible in the proxy: catalysis FADES as charges move out (-6.19 -> +3.47). The relaxation
tests whether FRAGMENTATION also fades. Compare to the all-shell K6_apot baseline (4.24, valid, open).

## Run (own dir per job)
for t in shell2 shell3 shell4; do mkdir -p run_K6_$t; cp relax_shell.py K6_${t}_coords.xyz K6_${t}_charges.txt run_K6_${t}.pbs run_K6_$t/; done
Dry-run (compute node): for t in shell2 shell3 shell4; do (cd run_K6_$t && python3 relax_shell.py K6_$t); done
Submit: cd run_K6_shellN && qsub run_K6_shellN.pbs && cd ..   (each)
Analyse: python3 ../analyse_relaxed.py run_K6_$t/K6_${t}_relaxed.xyz

## Decision
- Outer shells fragment LESS (shell4 intact where shell2 fragments) -> DISTANCE is a real rupture factor;
  validates the proximity-weighted penalty (Stage-1 QP) as the principled fix. But catalysis is weak far out
  -> confirms the tradeoff, motivating the OPTIMISED middle ground (penalty finds it).
- All shells still fragment/open -> distance NOT sufficient; problem is field character, not just proximity
  -> go straight to certified-gentle QP/MISOCP.
Either outcome sharpens the optimiser design. Cheap (3 jobs), rules out the simple explanation.
