# PIVOT RESULT: the chorismate Claisen barrier is strongly field-catalysable (OEEF)

## Finite-field barrier scan (fixed committed geometries, uniform field, CPCM eps=4)
F=0 barrier = 17.47 kcal/mol (reproduces the canonical bare barrier -- consistency check).
Barrier vs uniform field at |F|=0.02 a.u. (script verdict):
  -z  : 17.47 -> 6.85   (-10.63 kcal/mol)   <- max; the Delta-mu axis
  rxn : 17.47 -> 8.17   ( -9.31)            reaction-coordinate direction
  +x  : 17.47 -> 9.51   ( -7.97)
  -y  : 17.47 -> 11.82  ( -5.65)
Every axis is monotonic and sign-antisymmetric (one direction lowers, opposite raises) = textbook OEEF.

## Why this overturns the "uncatalysable" reading
Bare-dianion fragmentation and the proxy's +37 kcal/mol barrier were NOT the reaction's ceiling -- they
were the Sum q.dV proxy MIS-ALIGNING the field. Sum q.dV optimises differential POTENTIAL at grid
points, not alignment of the net FIELD with Delta-mu. A 0.02 a.u. field along Delta-mu lowers the
barrier ~10.6 kcal/mol. The electrostatic channel is wide open.

## Consequence: build the OEEF/response optimiser
Find the charge arrangement that maximises the field along the catalytic axis (-z / Delta-mu) at the
reacting region, subject to net-neutral, box, and a fragmentation trust-region ||B_frag q|| <= F_max.
Convex (linear obj + SOC), certifiable on Clarabel. Open design question: can buildable +-Q_max charges
at LJ standoff produce ~0.02 a.u. along Delta-mu at the reacting region WITHOUT fragmenting?

## Caveat
Fixed-geometry (electronic) sensitivity -- relaxation will shift the numbers, but F=0 reproduces the
canonical barrier and the direction/magnitude are robust. Confirm with a relaxed barrier under the best
achievable field once a design is built.
