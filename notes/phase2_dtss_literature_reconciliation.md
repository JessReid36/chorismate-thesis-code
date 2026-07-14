# DTSS literature reconciliation — Szefczyk 2004 (JACS, ja049376t) + Ranaghan 2004 (OBC, B313759g)
# vs. our Dv (difference-potential) agreement map at enzyme residue positions.
# Retrieved directly from the two PDFs (Szefczyk Table 1, MP2/6-31G(d) DTSS column).

## Signed per-residue DTSS (Szefczyk 2004, Table 1; kcal/mol; negative = TS-stabilizing/catalytic)
Arg90   -9.06   strongly catalytic (dominant; mutagenesis: dead enzyme when removed)
Arg7    -5.90   strongly catalytic (Arg7Ala: 10^6-fold kcat/Km reduction)
Glu78   -3.57   catalytic (conformation-dependent; see Szefczyk Tables 2-3)
XSOL124 -2.48   catalytic (crystallographic water)
Arg116  -2.45   catalytic
Arg63   -1.40   WEAKLY catalytic
Cys75   -0.81   weakly catalytic
Tyr108  +0.29   ~neutral / slightly anticatalytic
Phe57   +0.85   slightly anticatalytic
Lys60   +1.39   ANTICATALYTIC (a cation here RAISES the barrier)
sum    -23.26   total DTSS (MP2)

## Reconciliation with our placeholder-geometry Dv agreement map
CONFIRMED  Arg90 (=our AMBER Arg217): dominant catalytic cation; our method rediscovers a +1
           2.06 A away. Our Dv at Arg217 was negative (wants cation) - AGREE. Solid.
CONFIRMED  Lys60 (=our Lys60): DTSS = +1.39 (anticatalytic). Our Dv at Lys60 was POSITIVE
           (wants anion; enzyme put cation) -> DISAGREE. The literature independently confirms
           Lys60 is a cationic residue sitting where a cation raises the barrier - a BINDING
           residue at an electrostatic cost to the barrier. THIS is the confirmed poster child
           for the binding-vs-catalysis decomposition, not Arg63.
CORRECTED  Arg63 (=our AMBER Arg63): DTSS = -1.40 (weakly CATALYTIC). Our placeholder Dv gave
           +0.0014 (wants anion) - WRONG SIGN. This is a placeholder-geometry artifact. Drop the
           Arg63 sign-inversion claim. Re-running on the real Phase 1 TS should flip Arg63 to
           agree with Szefczyk; if it does, that is a VALIDATION of the method against MP2 DTSS.
CONFIRMED  Arg7, Arg116 (=our Arg134, and Arg243 by the +127 map): catalytic cations, consistent
           with our cation placements in that region.

## Precedent and novelty (Szefczyk 2004, p.7, Fig. 6-7)
Szefczyk et al. construct the identical object we call Dv: "the difference of the ... electrostatic
potentials ... of the substrate and the transition state," multiplied by -1, with "+"/"-" signs on
the vdW surface "where positive or negative charges accelerate the reaction," and explicitly propose
it "for the design of a completely new catalytic environment for a known reaction." This is a direct
PRECEDENT for our objective function (validates the physics), and sharpens our novelty statement:
 - Szefczyk 2004: qualitative VISUALIZATION of the differential field as a manual design aid.
 - GOCAT: black-box GENETIC-ALGORITHM search over charge placements (no optimality guarantee).
 - THIS WORK: quantitative CONVEX / MIQP inverse design on the same differential field, with a
   provable-optimality certificate and physical constraints (net charge, discreteness {-1,0,+1},
   minimum inter-charge separation, total-charge budget). We compute the optimal charge arrangement,
   not just a good one, and we do it in milliseconds from two QM calculations.

## Thesis-narrative consequence
Lead with (1) Arg90 rediscovery [bulletproof] and (2) the total-vs-differential ranking divergence:
DeYonker 2024 F-SAPT TOTAL interaction ranks Arg7 > Arg63 > Arg90 > Lys60 (binding), whereas Szefczyk
DIFFERENTIAL DTSS ranks Arg90 > Arg7 > Glu78 > Arg116 > Arg63, with Lys60 positive - binding-strong
is not catalysis-strong. Use Lys60 (+1.39 DTSS) as the confirmed "binding cation at an electrostatic
cost" example. Present the Arg63 placeholder mis-sign honestly and use its correction on real geometry
as a method-validation checkpoint. Do NOT claim "net-charged designs beat the enzyme" as a headline
until the real-geometry Dv is in and the net-free polarized advantage is re-confirmed there.
