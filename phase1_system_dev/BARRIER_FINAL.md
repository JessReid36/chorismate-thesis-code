# Phase-1 final barrier — freq-verified, consistent-region  [retained-final]

Chorismate->prephenate Claisen in BsCM (frame 820), QM/MM (ORCA 6.0.1 / AMBER22),
102-atom reduced active region, common frozen bulk (18c ts_start splice).

## Energies (FINAL SINGLE POINT ENERGY (QM/MM), from the Opt/OptTS outputs)
Reactant  18e_reactant_reduced/reactant_reduced.out : -1076.130495 Ha
TS        18c_reduced_region/ts_reduced.out         : -1076.104998 Ha
Product   18f_product_reduced/product_reduced.out   : -1076.151454 Ha

## Result
Barrier (TS - R) = +16.00 kcal/mol
Reaction energy (P - R) = -13.15 kcal/mol

## Verification (freq, NumFreq on the 102-atom active region, 612 displacements each)
Reactant : ORCA TERMINATED NORMALLY, 612/612, 0 imaginary modes -> true minimum
Product  : ORCA TERMINATED NORMALLY, 612/612, 0 imaginary modes -> true minimum
TS       : one imaginary mode, -313.30 cm-1 (transition vector); IRC connects R and P
All three in the SAME 102-atom reduced region on the SAME frozen bulk (consistent-region).

## Cross-validation (barrier, kcal/mol)
reduced-region OptTS (this, retained-final) : 16.00
NEB-CI                                       : 15.94
relaxed scan                                 : 16.37
Spread ~0.4 kcal/mol; deviation from exp. dH! (12.7) ~ +3.3 kcal/mol.

## Supersedes
The earlier +15.08 kcal/mol barrier is SUPERSEDED: it subtracted the reduced-region TS from a
~1959-atom-region reactant/product (mixed active region -> apples-to-oranges). R and P were
re-optimised in the SAME 102-atom region on the common 18c bulk to fix this (see TIER2_KNOWN_ISSUES
B3). +16.00 is the consistent-region value and is the retained-final Phase-1 barrier.
