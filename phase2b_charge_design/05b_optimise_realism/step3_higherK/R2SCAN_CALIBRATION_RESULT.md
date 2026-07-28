# r2SCAN-3c calibration RESULT: FAILS the geometry gate → DFT-only (pre-registered)

## Full table (DFT B3LYP-D3BJ/def2-SVP vs r2SCAN-3c, method-only difference, same designs)
| design     | DFT C1-C6  O3-C4  conn | r2SCAN C1-C6 O3-C4 conn | dC1-C6 | note |
| K6_apot    | 4.243  1.406  1mol(VALID) | 5.861  4.415  FRAG | 1.618 | **DECISIVE: DFT valid, r2SCAN ruptures** |
| K6_aforce  | 5.612  6.004  FRAG        | 6.947  6.799  FRAG | 1.334 | both fragment (agree on instability) |
| K10_apot   | 4.201  1.421  FRAG        | 4.505  7.483  FRAG | 0.304 | both fragment |
| K10_aforce | (DFT frag, 10.49/7.15)    | no relaxed xyz     |  --   | r2SCAN job incomplete; moot |

Summary: mean |dC1-C6| 1.086 A, max 1.618 A (>> 0.2 tolerance). Window-verdict 3/3 (misleading - "agreement"
is mostly both-say-open on already-broken designs). Connectivity 2/3, and the ONE disagreement is the load-
bearing case (K6_apot: DFT intact, r2SCAN fragmented). GATE: FAIL.

## Verdict
**r2SCAN-3c FAILS the geometry calibration gate. DFT (B3LYP-D3BJ/def2-SVP) only for geometry.**
The decisive case is K6_apot - the single design where DFT gives a VALID intact reactant (4.243, 24/24). There
r2SCAN over-delocalises and RUPTURES the substrate (5.861, O3-C4 4.415, 8/24). Same design, method-only diff
-> pure method effect (not charges/placement - those held fixed). On designs DFT itself fragments, r2SCAN also
fragments (agreement on instability, uninformative about reliability).

## Why (mechanism) + why no cheap fix
r2SCAN is a PURE meta-GGA (no exact exchange) -> self-interaction/delocalisation error, catastrophic for the
trifecta here: anion (-2) + strong external field + breaking bond. Density spuriously spills toward the +ve
embedding charges; the bond can't hold; dissociates. B3LYP holds it via 20% exact exchange (cancels much SIE).
- Dispersion is NOT the issue (r2SCAN-3c already has D4+gCP; dispersion is attractive, can't cause rupture).
- Bigger/diffuse basis wouldn't fix it (SIE is a functional problem, not basis).
- Only EXACT-EXCHANGE hybrids (r2SCANh/r2SCAN0) fix SIE - but cost ~ B3LYP, forfeiting the speedup.
=> No cheap cure. Stay at DFT. Pre-registered in DESIGN_DECISIONS as a mandatory gate ("do not skip").

## Consequence
No cheap-method speedup on this system. All geometry work stays B3LYP-D3BJ/def2-SVP. Timeline unaffected
(that's the level we've run throughout). Calibration did its job: caught an unreliable shortcut BEFORE any
result was built on it - a methods-chapter point in favour.
