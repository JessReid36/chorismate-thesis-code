# Step 3 RESULT — K2 mu=100 (a_pot penalty): NEGATIVE. a_pot mode does NOT fix the distortion.

Relaxed the certified a_pot-preorganised K2 design (two +1 guanidiniums, sites 105/121, D_apot=+0.0006).
Job converged clean: ASH 191 min, 221 QM/MM steps, single connected molecule, E=-837.923 Eh.

## Measured relaxed geometry (analyse_relaxed.py)
  C1-C6 (forming) = 5.277 A   -> OPEN, essentially identical to the mu=0 distorting control (5.004)
  O3-C4 (breaking)= 1.462 A   -> intact ether (correct reactant)
  C4-C6           = 2.505 A
  connectivity    = SINGLE molecule (valid endpoint, not fragmented)

## Verdict: FALSIFIED a_pot as the preorganisation coordinate
The certified optimiser drove Sum q*a_pot -> ~0 (D=0.0006, gap 0.000), yet the substrate STILL relaxed
C1-C6 open to 5.28 A. Minimising the differential POTENTIAL across C1<->C6 (a_pot = V(C1)-V(C6)) does NOT
prevent the geometric opening. This confirms the Step-0 concern: a_pot measures charge-SEPARATION
tendency, but the pathology is GEOMETRIC separation. The ~50% a_pot/a_force sign disagreement (STEP0) was
the early warning; this is the empirical falsification the STEP0 record set up as the arbiter.

## This is a real, retained finding (not a failure to hide)
The frozen-proxy sweep + certificate worked exactly as designed; the certified optimum on the a_pot
objective was found provably. What failed is the CHOICE of a_pot as the mode - a scientific hypothesis
that Step 3 was built to test and correctly rejected. Reproducibility-first: negative result committed.

## Next (decision tree)
1. Swap to a_force (mechanical force-imbalance, already committed as the cross-check). a_force recommends
   DIFFERENT, dipole-preserving K2 designs (mu=100 -> sites 185(-1)/249(+1); mu=500 -> 117(+1)/150(-1)),
   NOT the two-+1 that failed. Relaxing an a_force design distinguishes two explanations:
     (a) a_force compact  -> coordinate was the problem; a_force is the right preorg mode; FIX WORKS.
     (b) a_force also open -> NO K=2 monopole pair can preorganise this dianion (confirms the original
         "two isolated monopoles cannot cradle a floppy dianion" intuition) -> redirect to higher-K or a
         neutral cradle, not a better penalty. Either outcome is scientifically informative.
2. The two-+1 motif (a_pot mu=100) is also now understood: it minimised a_pot by removing the -1, but
   without a balancing pull it could not hold the carbons - consistent with a_force keeping a dipole.

## Files
k2_mu100_relaxed.xyz : the relaxed geometry (C1-C6 5.277)
k2_mu100_relax.log   : ASH log (converged, 191 min)
