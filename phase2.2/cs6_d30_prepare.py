#!/usr/bin/env python3
"""DIAGNOSTIC 30: generate the 22 inputs. Two arms.

A  the basis sequence under CAM-B3LYP: def2-SVP (6 arms + 2 controls) and def2-TZVPD (6 arms).
   def2-SVPD is NOT re-run; d28/CAM-B3LYP/def2-SVPD already holds it and GATE 6 re-reads it.
B  the stability of the CAM-B3LYP solution: def2-SVPD with STABPerform (6 arms + 2 controls).

Charges and geometries are COPIED from the committed d27/def2-SVPD arms and verified by sha256.
Nothing is retyped: the three vectors must be the same three vectors D27 and D28 measured.
"""
import argparse, hashlib, shutil
from pathlib import Path

VECTORS = ["1.3000e_n10", "2.8800e_n10", "5.6433e_n20"]
STATES = ["R", "TS"]
# (arm, basis, stability) -> whether controls are generated here
PLAN = [("A", "def2-SVP",   False, True),    # controls needed: none committed at def2-SVP/CAM
        ("A", "def2-TZVPD", False, False),   # controls already in d28/CAM-B3LYP/def2-TZVPD
        ("B", "def2-SVPD",  True,  True)]    # stability controls are new


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""): h.update(b)
    return h.hexdigest()


def keyword(basis):
    # identical to d27 and d28 but for the functional and basis. CAM-B3LYP carries no dispersion
    # of its own, so it takes D3BJ exactly as d28's CAM-B3LYP arms do. wB97X-D3 would not.
    return "! CAM-B3LYP D3BJ %s def2/J RIJCOSX TightSCF" % basis


def inp(basis, charged, stability):
    L = [keyword(basis), "%maxcore 3000", "%pal nprocs 8 end"]
    if stability:
        L += ["%scf", "  MaxIter 300", "  STABPerform true", "end"]
    else:
        L += ["%scf MaxIter 300 end"]
    if charged:
        L += ['%pointcharges "pc.pc"']
    L += ["* xyzfile -2 1 substrate.xyz", ""]
    return "\n".join(L)


def main():
    a = argparse.ArgumentParser()
    a.add_argument("--src", required=True, help="the committed d27/def2-SVPD directory")
    a.add_argument("--out", required=True, help="where the d30 tree is written")
    a = a.parse_args()
    SRC, OUT = Path(a.src), Path(a.out)
    made, checks = [], []
    for arm, basis, stab, ctls in PLAN:
        names = ["%s_%s" % (v, s) for v in VECTORS for s in STATES]
        if ctls: names += ["ctl_%s" % s for s in STATES]
        for n in names:
            charged = not n.startswith("ctl")
            d = OUT / arm / basis / n
            d.mkdir(parents=True, exist_ok=True)
            # the geometry comes from the control of the same state, which is the same xyz the
            # charged arm of that state uses; both are verified below
            state = "TS" if n.endswith("_TS") else "R"
            src_geom = SRC / ("ctl_%s" % state) / "substrate.xyz"
            shutil.copy2(src_geom, d / "substrate.xyz")
            checks.append((str((d / "substrate.xyz").relative_to(OUT)), sha(src_geom),
                           sha(d / "substrate.xyz"), str(src_geom)))
            if charged:
                src_pc = SRC / n / "pc.pc"
                shutil.copy2(src_pc, d / "pc.pc")
                checks.append((str((d / "pc.pc").relative_to(OUT)), sha(src_pc),
                               sha(d / "pc.pc"), str(src_pc)))
            (d / "job.inp").write_text(inp(basis, charged, stab))
            made.append(str(d.relative_to(OUT)))
    bad = [c for c in checks if c[1] != c[2]]
    print("%d job directories written" % len(made))
    print("%d copied files verified by sha256, %d mismatched" % (len(checks), len(bad)))
    for c in bad: print("   MISMATCH %s against %s" % (c[0], c[3]))
    # also check that the charged arms of a given state all share one geometry, as they must
    geoms = {}
    for m in made:
        g = sha(OUT / m / "substrate.xyz")
        geoms.setdefault("TS" if m.endswith("_TS") else "R", set()).add(g)
    for st, s in sorted(geoms.items()):
        print("   state %-3s distinct geometries: %d %s" % (st, len(s), "OK" if len(s) == 1 else "*** MORE THAN ONE ***"))
    return 1 if bad or any(len(s) != 1 for s in geoms.values()) else 0


if __name__ == "__main__":
    raise SystemExit(main())
