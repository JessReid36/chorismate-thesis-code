#!/usr/bin/env python3
"""cs6_d51_expect_pyscf.py - what Claude expects of Stage 6 DIAGNOSTIC 51, computed in a second program.

NOT RUN ON hpc1 AND NOT RUN BY THE READER. It needs PySCF (written with 2.15.1), numpy and scipy, which the reader
does not. It was run once, on 23 November 2026, before D51's pre-registration was written, and what it wrote is the
committed cs6_d51_expected.py. It is kept so that the numbers in that file can be formed again.

What it computes, for one electron at a bare point charge q with the functions of D50 (def2-TZVPD hydrogen and an s
and a p function at each of 0.035, 0.012 [X2] and also 0.004 [X3]), the ghost atom 0.01 A from the charge as in D50:

  FULL   Hartree-Fock and CAM-B3LYP with exact integrals (no RI, no COSX) and an integration grid of 400 radial
         points to beyond 100000 bohr: what the functional gives in that set, free of the grid.
  LIMIT  the same in a set of 27 even-tempered s functions at the charge: what the functional gives in a complete
         basis. Hartree-Fock there is 13.605693 q^2 to 1e-6 eV.
  MODEL  CAM-B3LYP with a radial grid built as ORCA's manual describes its own (Gauss-Chebyshev points under the M3
         mapping of Treutler and Ahlrichs, r = xi/ln2 * ln(2/(1-x)), n_r = 15*IntAcc - 40 + b*ROW): n_r points, with
         xi fitted to ONE number, the electron count D50's committed output prints for set X2 at 0.25 e. Two models,
         n_r = 26 and 31 at IntAcc 4.388, since b is not known. Each then says what a larger IntAcc should give.

usage: python3 cs6_d51_expect_pyscf.py <results repository> <file to write>
"""
import csv
import os
import sys

import numpy as np
from pyscf import dft, gto, qmmm, scf
from scipy.integrate import lebedev_rule
from scipy.optimize import brentq

EV = 27.211386245988
RYD = EV / 2.0
LN2 = np.log(2.0)
SIZES = (10, 15, 20, 25, 30, 40, 50, 60, 75, 100)
ADD = {"X2": (0.035, 0.012), "X3": (0.035, 0.012, 0.004)}
RUNGS = (("g2", 4.388), ("g3", 4.959), ("i08", 8.0), ("i12", 12.0), ("i20", 20.0), ("i40", 40.0))
MODELS = (26, 31)
XC = "camb3lyp"


def mol_set(name):
    b = gto.basis.load("def2-tzvpd", "H") + [[0, [e, 1.0]] for e in ADD[name]] + [[1, [e, 1.0]] for e in ADD[name]]
    return gto.M(atom="GHOST-H 0.01 0 0", basis={"GHOST-H": b}, charge=-1, spin=1, unit="Angstrom", verbose=0)


def mol_large(lo=0.0004, hi=3000.0, ratio=1.8):
    ex, e = [], hi
    while e > lo:
        ex.append(e)
        e /= ratio
    return gto.M(atom="GHOST-H 0 0 0", basis={"GHOST-H": [[0, [e, 1.0]] for e in ex]}, charge=-1, spin=1, unit="Angstrom", verbose=0), len(ex)


def sphere(centre, r, dr, order=29):
    xyz, wl = lebedev_rule(order)
    co = (r[:, None, None] * xyz.T[None, :, :]).reshape(-1, 3) + np.asarray(centre)
    return co, ((r * r * dr)[:, None] * wl[None, :]).reshape(-1)


def grid_full(centre, n=400, scale=4.0):
    x, w = np.polynomial.legendre.leggauss(n)
    return sphere(centre, scale * (1 + x) / (1 - x), 2 * scale / (1 - x) ** 2 * w, 41)


def grid_m3(centre, n, xi):
    i = np.arange(1, n + 1)
    x = np.cos(i * np.pi / (n + 1))
    w = np.pi / (n + 1) * np.sin(i * np.pi / (n + 1)) ** 2 / np.sqrt(1 - x * x)
    return sphere(centre, xi / LN2 * np.log(2.0 / (1 - x)), xi / LN2 / (1 - x) * w)


def run(mol, q, grid=None):
    """(B in eV, the electron count the grid gives) for Hartree-Fock (grid None) or the functional on that grid"""
    mf = scf.UHF(mol) if grid is None else dft.UKS(mol)
    if grid is not None:
        mf.xc = XC
        mf.grids.coords, mf.grids.weights = grid
    mf = qmmm.mm_charge(mf, np.array([[0.0, 0.0, 0.0]]), np.array([q]), unit="Angstrom")
    mf.init_guess, mf.conv_tol, mf.max_cycle = "hcore", 1e-12, 300
    e = mf.kernel()
    if not mf.converged:
        raise SystemExit("not converged: q %.2f" % q)
    n = None
    if grid is not None:
        dm = mf.make_rdm1()
        ao = dft.numint.eval_ao(mol, grid[0])
        n = float((np.einsum("gi,ij,gj->g", ao, dm[0] + dm[1], ao) * grid[1]).sum())
    return -e * EV, n


def orca_count(s6, name, k):
    v = None
    for line in open(os.path.join(s6, "d50", "ONE_CAM_%sA_p%03d" % (name, k), "job.out"), errors="replace"):
        if line.startswith("N(Total)"):
            v = float(line.split()[2])
    return v


def main():
    if len(sys.argv) != 3:
        sys.stderr.write("usage: python3 cs6_d51_expect_pyscf.py <results repository> <file to write>\n")
        return 1
    s6 = os.path.join(sys.argv[1], "phase2.2/cs_stage6")
    if os.path.exists(sys.argv[2]):
        sys.stderr.write("STOP  %s exists already.\n" % sys.argv[2])
        return 1
    orca = dict(((r["set"], int(round(float(r["size_e"]) * 100))), r) for r in csv.DictReader(open(os.path.join(s6, "d50_calibration.tsv")), delimiter="\t"))
    mols = dict((name, mol_set(name)) for name in sorted(ADD))
    large, nl = mol_large()
    FULL, LIMIT, MODEL, FIT = {}, {}, {}, {}
    for k in SIZES:
        q = k / 100.0
        LIMIT[k] = (run(large, q)[0], run(large, q, grid_full(large.atom_coord(0)))[0])
        for name in sorted(ADD):
            m = mols[name]
            FULL[(name, k)] = (run(m, q)[0], run(m, q, grid_full(m.atom_coord(0)))[0])
    for n0 in MODELS:
        m = mols["X2"]
        xi = brentq(lambda x: run(m, 0.25, grid_m3(m.atom_coord(0), n0, x))[1] - orca_count(s6, "X2", 25), 0.5, 2.5, xtol=1e-6)
        dn, db = [], []
        for name in sorted(ADD):
            for rung, eps in RUNGS:
                n = n0 + int(round(15.0 * (eps - RUNGS[0][1])))
                for k in SIZES:
                    b, cnt = run(mols[name], k / 100.0, grid_m3(mols[name].atom_coord(0), n, xi))
                    MODEL[(n0, name, rung, k)] = (cnt, b)
                    if rung == "g2":
                        dn.append(cnt - orca_count(s6, name, k))
                        db.append(b - float(orca[(name, k)]["B_CAM_A_eV"]))
        FIT[n0] = (xi, float(np.sqrt(np.mean(np.square(dn)))), max(abs(v) for v in dn), float(np.sqrt(np.mean(np.square(db)))), max(abs(v) for v in db))
    L = ['"""cs6_d51_expected.py - what Claude expects of Stage 6 DIAGNOSTIC 51: numbers only.',
         "Written by cs6_d51_expect_pyscf.py (PySCF, a second program) on 23 November 2026, before the pre-registration.",
         "COMPUTED BY CLAUDE. Not ORCA and not a measurement. B is minus the energy of the one electron, eV.",
         "FULL[(set, size)] = (B at Hartree-Fock, B with CAM-B3LYP on a grid that counts the whole electron), exact integrals.",
         "LIMIT[size] = the same two in a set of %d even-tempered s functions at the charge." % nl,
         "MODEL[(n0, set, rung, size)] = (electron count, B with CAM-B3LYP) on a model of ORCA's radial grid with n0 radial",
         "points at IntAcc 4.388 and 15 more for each unit of IntAcc.",
         "FIT[n0] = (xi in bohr, then over D50's 20 runs with treatment A: rms and largest difference of the count,",
         'rms and largest difference of B in eV, model minus D50). Sizes are in hundredths of e."""',
         "LARGE_SET = %d" % nl, "MODELS = %r" % (MODELS,), "RUNGS = %r" % (RUNGS,), "FULL = {"]
    for key in sorted(FULL):
        L.append("    %r: (%.6f, %.6f)," % (key, FULL[key][0], FULL[key][1]))
    L += ["}", "LIMIT = {"]
    for key in sorted(LIMIT):
        L.append("    %r: (%.6f, %.6f)," % (key, LIMIT[key][0], LIMIT[key][1]))
    L += ["}", "MODEL = {"]
    order = dict((r, i) for i, (r, e) in enumerate(RUNGS))
    for key in sorted(MODEL, key=lambda t: (t[0], t[1], order[t[2]], t[3])):
        L.append("    %r: (%.4f, %.4f)," % (key, MODEL[key][0], MODEL[key][1]))
    L += ["}", "FIT = {"]
    for key in sorted(FIT):
        L.append("    %r: (%.4f, %.4f, %.4f, %.4f, %.4f)," % ((key,) + FIT[key]))
    L += ["}", ""]
    with open(sys.argv[2], "w") as f:
        f.write("\n".join(L))
    print("written: %d FULL, %d LIMIT, %d MODEL; fit %r" % (len(FULL), len(LIMIT), len(MODEL), FIT))
    return 0


if __name__ == "__main__":
    sys.exit(main())
