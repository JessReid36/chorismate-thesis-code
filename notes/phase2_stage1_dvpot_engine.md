Stage 1 - difference-potential engine.

The catalytic objective is the difference electrostatic potential Dv(r) = V_TS(r) - V_R(r),
the change in the substrate's own electrostatic potential between the transition-state and
reactant geometries. For a fixed density this is exactly the first-order (frozen-density)
interaction available to an external point charge, so the barrier change from a set of
charges {q_i} at points {r_i} is Sum_i q_i Dv(r_i), linear and pairwise-additive to first
order (polarization is a separate second-order correction, treated later).

V_R and V_TS were computed with orca_vpot from single-point densities at B3LYP-D3BJ/def2-SVP
(def2/J, RIJCOSX) with CPCM (epsilon = 4); the dianion is unbound in vacuum (Stage 1.1a) and
the continuum is applied identically at both geometries so it cancels in the difference. This
ORCA 6.0.1 build stores the SCF density as member <base>.scfp inside the <base>.densities
container; orca_vpot reads it via the container-basename argument.

The convention and extraction were validated by a sign check at physically motivated test
points: a +1 charge placed ~2 A beyond the breaking C4-O3 ether oxygen (along C4->O3) gives
q*Dv = -0.0122 Eh < 0, i.e. stabilizing - the developing negative charge on the ether oxygen
at the TS is preferentially stabilized by a cation, reproducing the catalytic role of Arg90.
A contrast point in the forming C6-C1 region gives q*Dv = +0.0027 Eh (destabilizing), the
opposite sign, consistent with charge depletion there on going to the TS. Dv therefore
reproduces the correct spatial structure of differential transition-state stabilization for
the chorismate rearrangement.
