#!/usr/bin/env python3
"""CONTROL TEST: is the continuum (CPCM eps=4) responsible for the substrate-collapse
under charge, or is the collapse intrinsic? Relax the SAME reactant geometry 4 ways:
  1. bare + CPCM eps=4   2. bare + vacuum   3. Arg90 charges + CPCM   4. Arg90 charges + vacuum
Readout: does C1-C6 (atoms 0-12) stay near-attack (~3.1) or open into the non-reactive well (~3.9)?
Committed refs: bare/CPCM opens 3.12->3.538 (intact); Arg90/CPCM free-relax opens 3.12->3.88 (collapse).
"""
import os, sys, subprocess, re
import numpy as np
ORCADIR = "/home/apps2/ORCA/6.0.1"
NCORES = 8
BASE = os.path.expanduser("~/07_rebuild")

# reacting-atom indices (0-based)
C1, C6, O3, C4 = 0, 12, 7, 8

def read_substrate_24(path):
    """First 24 atoms (element x y z) from an xyz (substrate)."""
    L = open(path).read().splitlines()
    n = int(L[0].split()[0])
    return [l for l in L[2:2+24]]

def read_charges_tail(coords_path, charges_path):
    """External charge lines (atoms 24..end) as (q, x, y, z) for the .pc file."""
    L = open(coords_path).read().splitlines()
    n = int(L[0].split()[0])
    ext_xyz = L[2+24:2+n]  # lines after the 24 substrate atoms
    q = [float(x) for x in open(charges_path).read().strip().split(",")]
    ext_q = q[24:]
    pc = []
    for i, line in enumerate(ext_xyz):
        p = line.split()
        pc.append((ext_q[i], float(p[1]), float(p[2]), float(p[3])))
    return pc

def dist(lines, i, j):
    a = np.array([float(x) for x in lines[i].split()[1:4]])
    b = np.array([float(x) for x in lines[j].split()[1:4]])
    return np.linalg.norm(a-b)

def write_pc_file(path, pc):
    with open(path, "w") as f:
        f.write("%d\n" % len(pc))
        for (q,x,y,z) in pc:
            f.write("%.6f %.8f %.8f %.8f\n" % (q,x,y,z))

def run_relax(tag, substrate_lines, use_cpcm, pc=None):
    """Constrained-free geometry opt of the 24-atom substrate. CPCM optional. External .pc optional.
    Returns the final geometry lines (element x y z)."""
    inp = "ctrl_%s.inp" % tag
    simple = "! B3LYP D3BJ def2-SVP def2/J RIJCOSX TightSCF Opt"
    blocks = "%scf MaxIter 300 end\n"
    if use_cpcm:
        blocks += "%cpcm epsilon 4.0 end\n"
    
    if pc is not None:
        write_pc_file("ctrl_%s.pc" % tag, pc)
        blocks = blocks + "%pointcharges \"ctrl_" + tag + ".pc\"\n"
    with open(inp, "w") as f:
        f.write(simple+"\n"+blocks+"* xyz -2 1\n")
        for l in substrate_lines:
            f.write(l+"\n")
        f.write("*\n")
    subprocess.run(["%s/orca"%ORCADIR, inp], stdout=open("ctrl_%s.out"%tag,"w"),
                   stderr=subprocess.STDOUT, env=dict(os.environ, PATH=ORCADIR+":"+os.environ["PATH"]))
    # read optimized geometry from the .xyz ORCA writes
    optxyz = "ctrl_%s.xyz" % tag
    if os.path.exists(optxyz):
        L = open(optxyz).read().splitlines()
        n = int(L[0].split()[0])
        return [l for l in L[2:2+n]], True
    # fallback: check if opt converged at all
    out = open("ctrl_%s.out"%tag).read()
    conv = "HURRAY" in out or "THE OPTIMIZATION HAS CONVERGED" in out
    return None, conv

def main():
    react = read_substrate_24(BASE+"/L0_groundtruth/L0_react_opt.xyz")
    print("=== starting reactant geometry ===")
    print("  C1-C6 = %.3f, O3-C4 = %.3f (0-based atoms %d-%d, %d-%d)" % (
        dist(react,C1,C6), dist(react,O3,C4), C1,C6,O3,C4))
    # load Arg90 charges (positions + magnitudes)
    pc = read_charges_tail(BASE+"/L4_barrier/design_cradleArg90_coords.xyz",
                           BASE+"/L4_barrier/design_cradleArg90_charges.txt")
    print("  Arg90 external charges: %d (sum q = %+.2f)" % (len(pc), sum(p[0] for p in pc)))
    print()
    conditions = [
        ("1_bare_cpcm",   True,  None, "bare + CPCM eps=4"),
        ("2_bare_vac",    False, None, "bare + VACUUM"),
        ("3_arg90_cpcm",  True,  pc,   "Arg90 charges + CPCM eps=4"),
        ("4_arg90_vac",   False, pc,   "Arg90 charges + VACUUM"),
    ]
    results = {}
    for tag, cpcm, charges, label in conditions:
        print(">> running %s (%s)..." % (tag, label), flush=True)
        geom, ok = run_relax(tag, react, cpcm, charges)
        if geom is None:
            print("   !! no optimized geometry (converged=%s) - check ctrl_%s.out" % (ok, tag))
            results[tag] = None
            continue
        c16 = dist(geom, C1, C6); o34 = dist(geom, O3, C4)
        results[tag] = (c16, o34)
        print("   final C1-C6 = %.3f, O3-C4 = %.3f" % (c16, o34))
    print()
    print("=" * 68)
    print("RESULTS: does the substrate stay near-attack or collapse?")
    print("  (near-attack ~3.1 A; collapse into non-reactive well ~3.9 A)")
    print("=" * 68)
    print("%-28s %10s %10s" % ("condition", "C1-C6", "O3-C4"))
    print("%-28s %10.3f %10.3f   (start)" % ("reactant (start)", dist(react,C1,C6), dist(react,O3,C4)))
    for tag, cpcm, charges, label in conditions:
        if results[tag]:
            c16, o34 = results[tag]
            flag = "COLLAPSE" if c16 > 3.6 else "near-attack"
            print("%-28s %10.3f %10.3f   %s" % (label, c16, o34, flag))
    print()
    print("=== INTERPRETATION ===")
    r = results
    if r.get("3_arg90_cpcm") and r.get("4_arg90_vac"):
        c_cpcm = r["3_arg90_cpcm"][0]; c_vac = r["4_arg90_vac"][0]
        print("Charged case: CPCM C1-C6 = %.3f vs VACUUM C1-C6 = %.3f" % (c_cpcm, c_vac))
        if c_cpcm > 3.6 and c_vac > 3.6:
            print(">>> Both collapse -> the continuum is NOT the cause; charge-driven collapse is INTRINSIC.")
        elif c_cpcm > 3.6 and c_vac < 3.6:
            print(">>> CPCM collapses but VACUUM does not -> the CONTINUUM (or charge-continuum mismatch)")
            print(">>> IS implicated in the collapse. Worth reconsidering the electrostatics setup.")
        elif c_cpcm < 3.6 and c_vac < 3.6:
            print(">>> Neither collapses here -> this relaxation didn't reproduce the earlier collapse;")
            print(">>> may need the exact earlier conditions (free vs constrained) to compare.")
        else:
            print(">>> VACUUM collapses but CPCM does not -> continuum is STABILISING, opposite of the worry.")
    if r.get("1_bare_cpcm") and r.get("2_bare_vac"):
        print("Bare case: CPCM C1-C6 = %.3f vs VACUUM C1-C6 = %.3f (are they similar?)" % (
            r["1_bare_cpcm"][0], r["2_bare_vac"][0]))

if __name__ == "__main__":
    main()
