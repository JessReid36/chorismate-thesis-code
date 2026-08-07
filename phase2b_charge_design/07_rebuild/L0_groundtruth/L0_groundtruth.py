#!/usr/bin/env python3
"""Layer 0 -- GROUND TRUTH. No external charges. No oracle machinery. No inherited code.

Establishes the anchor the whole rebuild rests on, from the committed geometries only:
  (a) the bare reactant relaxes and STAYS INTACT (O3-C4 ~1.45, no bond changes), and
  (b) the bare FORWARD BARRIER reproduces the canonical +17.47 kcal/mol.

Theory: B3LYP D3BJ / def2-SVP / def2/J / RIJCOSX / CPCM(eps=4), charge -2, singlet.
Barrier = E(TS) - E(R) as single points on the TRUSTED committed R/TS geometries (they are the
verified Phase-1/2b stationary points; we only confirm the theory reproduces the number). A bare
reactant relaxation confirms geometric stability. If either check fails, STOP.
"""
import os, subprocess

HERE   = os.path.dirname(os.path.abspath(__file__))
INP    = os.path.join(HERE, "..", "inputs")
ORCA   = "/home/apps2/ORCA/6.0.1/orca"
SIMPLE = "! B3LYP D3BJ def2-SVP def2/J RIJCOSX CPCM TightSCF"
BLOCKS = "%cpcm epsilon 4.0 end\n%scf MaxIter 300 end"
REF_BARRIER = 17.47
HARTREE_KCAL = 627.5094740631

def read_xyz_body(path):
    L = open(path).read().splitlines(); n = int(L[0].split()[0]); return L[2:2+n]

def write_inp(tag, body, opt=False):
    with open(os.path.join(HERE, tag + ".inp"), "w") as f:
        f.write(SIMPLE + (" Opt" if opt else "") + "\n" + BLOCKS + "\n* xyz -2 1\n")
        for ln in body: f.write(ln.rstrip() + "\n")
        f.write("*\n")

def run(tag):
    with open(os.path.join(HERE, tag + ".out"), "w") as fo:
        subprocess.run([ORCA, tag + ".inp"], cwd=HERE, stdout=fo, stderr=subprocess.STDOUT, check=False)

def final_energy(out):
    E = None
    for ln in open(out):
        if "FINAL SINGLE POINT ENERGY" in ln: E = float(ln.split()[-1])
    return E

def dist(body, i, j):
    a=[float(x) for x in body[i].split()[1:4]]; b=[float(x) for x in body[j].split()[1:4]]
    return sum((a[k]-b[k])**2 for k in range(3))**0.5

def main():
    react = read_xyz_body(os.path.join(INP, "reactant.xyz"))
    ts    = read_xyz_body(os.path.join(INP, "ts.xyz"))
    assert len(react) == 24 and len(ts) == 24, "expected 24-atom substrate"
    log = []

    write_inp("L0_react_sp", react, opt=False); run("L0_react_sp")
    write_inp("L0_ts_sp",    ts,    opt=False); run("L0_ts_sp")
    Er  = final_energy(os.path.join(HERE, "L0_react_sp.out"))
    Ets = final_energy(os.path.join(HERE, "L0_ts_sp.out"))
    assert Er and Ets, "an SCF did not complete"
    barrier = (Ets - Er) * HARTREE_KCAL
    ok_b = abs(barrier - REF_BARRIER) < 1.0
    log.append("BARRIER (single-point, committed geoms): %.2f kcal/mol  (ref %.2f)  %s"
               % (barrier, REF_BARRIER, "PASS" if ok_b else "FAIL"))

    o34_before = dist(react, 7, 8)
    write_inp("L0_react_opt", react, opt=True); run("L0_react_opt")
    optxyz = os.path.join(HERE, "L0_react_opt.xyz")
    if os.path.exists(optxyz):
        after = read_xyz_body(optxyz)
        o34 = dist(after, 7, 8); c16 = dist(after, 0, 12)
        ok_g = o34 < 1.8
        log.append("BARE REACTANT RELAX: O3-C4 %.3f -> %.3f, C1-C6 %.3f  %s"
                   % (o34_before, o34, c16, "PASS (intact)" if ok_g else "FAIL"))
    else:
        log.append("BARE REACTANT RELAX: no optimised xyz -- check L0_react_opt.out"); ok_g = False

    log.append("L0 PASS -- foundation trusted" if (ok_b and ok_g) else "L0 FAIL -- STOP, fix foundation")
    open(os.path.join(HERE, "L0_RESULT.txt"), "w").write("\n".join(log) + "\n")
    print("\n".join(log))

if __name__ == "__main__":
    main()
