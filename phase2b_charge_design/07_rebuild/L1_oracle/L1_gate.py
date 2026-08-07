#!/usr/bin/env python3
"""Layer 1 acceptance gate: zero-charge evaluate() must reproduce Layer 0."""
import os, sys
INP = sys.argv[1] if len(sys.argv) > 1 else "../inputs"

sub = open(os.path.join(INP, "reactant.xyz")).read().splitlines()[2:26]
open("zero_coords.xyz", "w").write("24\nzero-charge design\n" + "\n".join(l.rstrip() for l in sub) + "\n")
open("zero_charges.txt", "w").write(",".join(["0.0000"] * 24))
print("built zero-charge design (24 substrate atoms, 24 zero charges)")

import oracle
res = oracle.evaluate("zero", "reactant", workdir=".", inputs_dir=INP)
g = res["geometry"]
def d(i, j): return sum((g[i][1][k] - g[j][1][k]) ** 2 for k in range(3)) ** 0.5
o34 = d(7, 8); c16 = d(0, 12)
ok_int = res["integrity"] == "VALID"; ok_geo = o34 < 1.8
lines = ["L1 GATE (zero-charge evaluate vs L0):",
    "  integrity = %s   %s" % (res["integrity"], "PASS" if ok_int else "FAIL"),
    "  O3-C4 = %.3f (L0: 1.448)   %s" % (o34, "PASS" if ok_geo else "FAIL"),
    "  C1-C6 = %.3f (L0: 3.538)" % c16,
    "  energy = %.6f Eh" % res["energy"],
    "L1 %s" % ("PASS -- oracle trusted" if (ok_int and ok_geo) else "FAIL -- STOP")]
open("L1_RESULT.txt", "w").write("\n".join(lines) + "\n")
print("\n".join(lines))
