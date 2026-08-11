#!/usr/bin/env python3
"""Diagnostic: is chorismate's Claisen barrier field-sensitive AT ALL? (OEEF ceiling test)"""
import os, sys, subprocess
import numpy as np

ORCADIR = "/home/apps2/ORCA/6.0.1"
SIMPLE  = "! B3LYP D3BJ def2-SVP def2/J RIJCOSX CPCM TightSCF"
BLOCKS_CPCM = "%cpcm epsilon 4.0 end"
NCORES  = 8
HARTREE_KCAL = 627.5094740631
INP = sys.argv[1] if len(sys.argv) > 1 else "."
F_ACHIEVABLE = 0.02
FIELDS = [0.0, 0.005, 0.01, 0.02]
AXES = {"x": (1,0,0), "y": (0,1,0), "z": (0,0,1)}


def reaction_axis():
    L = open(os.path.join(INP, "reactant.xyz")).read().splitlines()[2:26]
    R = np.array([[float(v) for v in l.split()[1:4]] for l in L])
    m_break = 0.5*(R[7]+R[8]); m_form = 0.5*(R[0]+R[12])
    u = m_form - m_break; return u/np.linalg.norm(u)


def write_inp(path, geomfile, F):
    body = open(os.path.join(INP, geomfile)).read().splitlines()[2:26]
    with open(path, "w") as f:
        f.write(SIMPLE + "\n" + BLOCKS_CPCM + "\n")
        f.write("%%scf\n  EField %.6f, %.6f, %.6f\nend\n" % (F[0], F[1], F[2]))
        f.write("%%pal nprocs %d end\n" % NCORES)
        f.write("* xyz -2 1\n")
        for l in body:
            f.write(l.rstrip() + "\n")
        f.write("*\n")


def run_orca(inp):
    base = inp[:-4]
    with open(base + ".out", "w") as o:
        subprocess.run(["%s/orca" % ORCADIR, inp], stdout=o, stderr=subprocess.STDOUT)
    E = None
    for ln in open(base + ".out"):
        if "FINAL SINGLE POINT ENERGY" in ln:
            E = float(ln.split()[-1])
    return E


def run_orca_for(geomfile, F, tag):
    inp = "ff_%s.inp" % tag
    write_inp(inp, geomfile, F)
    return run_orca(inp)


def barrier_at(F, tag):
    er = run_orca_for("reactant.xyz", F, "r_%s" % tag)
    et = run_orca_for("ts.xyz",       F, "t_%s" % tag)
    if er is None or et is None:
        return None
    return (et - er) * HARTREE_KCAL


def main():
    rax = reaction_axis()
    axes = dict(AXES); axes["rxn"] = tuple(rax)
    b0 = barrier_at((0,0,0), "F0")
    print("bare (F=0) fixed-geometry barrier = %.2f kcal/mol" % b0)
    print("\naxis   F(a.u.)   barrier(kcal)   dBarrier_vs_F0")
    results = {}
    for name, u in axes.items():
        for Fm in FIELDS:
            if Fm == 0.0:
                results[(name, 0.0)] = b0; continue
            for sgn in ([+1, -1] if Fm > 0 else [+1]):
                F = tuple(sgn*Fm*np.array(u))
                b = barrier_at(F, "%s_%+.3f" % (name, sgn*Fm))
                results[(name, sgn*Fm)] = b
                print("%-5s %+7.3f   %8.2f      %+8.2f" % (name, sgn*Fm, b, b-b0))
    best_drop = 0.0; best = None
    for (name, Fm), b in results.items():
        if b is None or abs(Fm) > F_ACHIEVABLE: continue
        if (b - b0) < best_drop:
            best_drop = b - b0; best = (name, Fm)
    print("\n=== VERDICT ===")
    print("max barrier LOWERING within |F|<=%.3f a.u.: %.2f kcal/mol (axis %s)"
          % (F_ACHIEVABLE, best_drop, best))
    if best_drop > -3.0:
        print("  -> barrier is FIELD-INSENSITIVE (<3 kcal/mol lowering). No external-charge scheme")
        print("     (proxy OR response) catalyses meaningfully -> PIVOT TO STERIC.")
    else:
        print("  -> barrier IS field-sensitive; a response-based (OEEF) optimiser is justified.")
    open("field_sensitivity_result.txt", "w").write(
        "b0_kcal=%.2f\nbest_drop_kcal=%.2f\nbest_axis=%s\nF_achievable=%.3f\n"
        % (b0, best_drop, str(best), F_ACHIEVABLE))


if __name__ == "__main__":
    main()
