#!/usr/bin/env python3
"""cs6_d52_expect_pyscf.py - what Claude expects of Stage 6 DIAGNOSTIC 52, computed in a second program.

NOT RUN ON hpc1 AND NOT RUN BY THE READER. It needs PySCF (written with 2.15.1), numpy and scipy, which the reader
does not. It was run once, on 23 November 2026, before D52's pre-registration was written, and what it wrote is the
committed cs6_d52_expected.py. It is kept so that the numbers in that file can be formed again.

One electron in the field of charges of +0.25 e at positive sites of the design B_n11_q0.25: no substrate and no
negative charge. The positions are read from committed files: the charges from the design's point-charge file in
D50, the ghost atoms (0.01 A beyond each charge) from its input.

  LIMIT    the lowest level of the one-electron Hamiltonian in a large even-tempered basis of s, p and d functions
           at every charge of the system. This is the reference the runs are held against: for one centre it is
           13.605693 q^2. LIMIT2 is the same in a still larger basis, to show that LIMIT is converged.
  HF       the lowest level in the functions an ORCA run has: set X2 or X3 at ALL SIX ghost atoms. Hartree-Fock is
           exact for one electron in the functions it is given, so ORCA's Hartree-Fock run should give this number.
  CAM      CAM-B3LYP in those functions, with exact integrals (no RI, no COSX), on a grid built here: at each ghost
           atom 120 Gauss-Legendre radial points (r = 4 bohr (1+x)/(1-x)) on a Lebedev sphere of 302 points, the
           space shared out among the six centres by Becke's cells. With it, the electrons the grid counts.
           CAMCHECK repeats one case on a finer grid (160 radial points, 590 on the sphere).

usage: python3 cs6_d52_expect_pyscf.py <results repository> <file to write> [<log file>]
"""
import hashlib
import os
import sys
import time

import numpy as np
from pyscf import dft, gto, qmmm, scf
from scipy.integrate import lebedev_rule

EV = 27.211386245988
ANG = 0.52917721090380
S6 = "phase2.2/cs_stage6"
SRC = "d50/CAMTZX_B_n11_q0.25_X2A_di_R"
SITES = (1124, 1073, 1107, 1359, 1005, 1052)
SYSTEMS = (("s1", (1124,)), ("pc", (1073, 1052)), ("pf", (1124, 1052)), ("t3", (1124, 1073, 1107)),
           ("f4", (1124, 1073, 1107, 1052)), ("x6", SITES))
ADD = {"X2": (0.035, 0.012), "X3": (0.035, 0.012, 0.004)}
Q = 0.25
XC = "camb3lyp"
GRID = (120, 29)          # radial points; scipy's order of the Lebedev rule (29: 302 points)
GRID2 = (160, 41)         # the finer grid of CAMCHECK (41: 590 points)
CHECK = ("X3", "f4")


def sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


def positions(res):
    pcl = [l for l in open(os.path.join(res, S6, SRC, "pc.pc")).read().split("\n")[1:] if l.split()]
    pos = [tuple(float(v) for v in l.split()[1:4]) for l in pcl if float(l.split()[0]) > 0]
    gh = [tuple(float(v) for v in l.split()[1:4]) for l in open(os.path.join(res, S6, SRC, "job.inp")).read().split("\n") if l.startswith("H: ")]
    assert len(pos) == len(gh) == len(SITES)
    for c, g in zip(pos, gh):
        assert abs(np.linalg.norm(np.array(c) - np.array(g)) - 0.01) < 2e-6
    return dict(zip(SITES, pos)), dict(zip(SITES, gh))


def even(hi, lo, ratio):
    out, e = [], hi
    while e > lo:
        out.append(e)
        e /= ratio
    return out


def large(which):
    s, p, d = (even(60.0, 0.0015, 2.0), even(6.0, 0.004, 2.0), even(1.5, 0.01, 2.2)) if which == 1 else \
              (even(200.0, 0.0008, 1.6), even(20.0, 0.002, 1.7), even(5.0, 0.005, 1.8))
    return [[0, [e, 1.0]] for e in s] + [[1, [e, 1.0]] for e in p] + [[2, [e, 1.0]] for e in d]


def functions(name):
    return gto.basis.load("def2-tzvpd", "H") + [[0, [e, 1.0]] for e in ADD[name]] + [[1, [e, 1.0]] for e in ADD[name]]


def lowest(centres, basis, charged):
    """minus the lowest level, eV: functions 'basis' at 'centres', charges of +Q at 'charged' (all in angstrom)"""
    mol = gto.Mole()
    mol.atom = [["GHOST-H", c] for c in centres]
    mol.basis = {"GHOST-H": basis}
    mol.unit = "Angstrom"
    mol.charge, mol.spin, mol.verbose = -1, 1, 0
    mol.build()
    S, H = mol.intor("int1e_ovlp"), mol.intor("int1e_kin")
    for c in charged:
        with mol.with_rinv_origin(np.array(c) / ANG):
            H = H - Q * mol.intor("int1e_rinv")
    w, U = np.linalg.eigh(S)
    keep = w > 1e-7
    X = U[:, keep] / np.sqrt(w[keep])
    return -np.linalg.eigvalsh(X.T @ H @ X)[0] * EV


def becke(centres, nrad, order, scale=4.0):
    """a grid for several centres of one kind (bohr): Gauss-Legendre in x with r = scale (1+x)/(1-x), a Lebedev sphere,
    and Becke's cells (three iterations of his polynomial, equal sizes)"""
    x, w = np.polynomial.legendre.leggauss(nrad)
    r, dr = scale * (1 + x) / (1 - x), 2 * scale / (1 - x) ** 2 * w
    xyz, wl = lebedev_rule(order)
    C = np.asarray(centres)
    co, we = [], []
    for i in range(len(C)):
        p = (r[:, None, None] * xyz.T[None, :, :]).reshape(-1, 3) + C[i]
        wt = ((r * r * dr)[:, None] * wl[None, :]).reshape(-1)
        d = np.linalg.norm(p[:, None, :] - C[None, :, :], axis=2)
        P = np.ones_like(d)
        for a in range(len(C)):
            for b in range(len(C)):
                if a != b:
                    mu = (d[:, a] - d[:, b]) / np.linalg.norm(C[a] - C[b])
                    for _ in range(3):
                        mu = 1.5 * mu - 0.5 * mu ** 3
                    P[:, a] *= 0.5 * (1 - mu)
        co.append(p)
        we.append(wt * P[:, i] / P.sum(axis=1))
    return np.vstack(co), np.concatenate(we)


def functional(ghosts, name, charged, grid):
    """(minus the energy of the one electron in eV, the electrons the grid counts) for CAM-B3LYP in set 'name' at the six ghosts"""
    mol = gto.M(atom=[["GHOST-H", ghosts[s]] for s in SITES], basis={"GHOST-H": functions(name)}, charge=-1, spin=1, unit="Angstrom", verbose=0)
    mf = dft.UKS(mol)
    mf.xc = XC
    mf.grids.coords, mf.grids.weights = becke(mol.atom_coords(), grid[0], grid[1])
    mf = qmmm.mm_charge(mf, np.array(charged), np.array([Q] * len(charged)), unit="Angstrom")
    mf.init_guess, mf.conv_tol, mf.max_cycle = "hcore", 1e-10, 300
    e = mf.kernel()
    assert mf.converged
    dm = mf.make_rdm1()
    n = 0.0
    for i in range(0, len(mf.grids.weights), 100000):
        ao = dft.numint.eval_ao(mol, mf.grids.coords[i:i + 100000])
        n += float((np.einsum("gi,ij,gj->g", ao, dm[0] + dm[1], ao) * mf.grids.weights[i:i + 100000]).sum())
    return -e * EV, n


def hf_scf(ghosts, name, charged):
    """the same Hartree-Fock number by the SCF route, as a check of the route the functional takes"""
    mol = gto.M(atom=[["GHOST-H", ghosts[s]] for s in SITES], basis={"GHOST-H": functions(name)}, charge=-1, spin=1, unit="Angstrom", verbose=0)
    mf = qmmm.mm_charge(scf.UHF(mol), np.array(charged), np.array([Q] * len(charged)), unit="Angstrom")
    mf.init_guess, mf.conv_tol = "hcore", 1e-12
    e = mf.kernel()
    assert mf.converged
    return -e * EV


def main():
    res, out = sys.argv[1], sys.argv[2]
    log = open(sys.argv[3], "w") if len(sys.argv) > 3 else sys.stdout

    def say(t):
        log.write(t + "\n")
        log.flush()
    site, ghost = positions(res)
    import pyscf
    import scipy
    L = ['"""cs6_d52_expected.py - what Claude expects of DIAGNOSTIC 52: written by cs6_d52_expect_pyscf.py on 23 November 2026,',
         "before D52's pre-registration. COMPUTED IN A SECOND PROGRAM (PySCF %s, numpy %s, scipy %s); not ORCA, not a" % (pyscf.__version__, np.__version__, scipy.__version__),
         "measurement and not part of the record. One electron, charges of +0.25 e at positive sites of B_n11_q0.25.",
         "LIMIT: minus the lowest level in a large basis at the charges, eV. LIMIT2: the same in a still larger basis.",
         "HF: minus the lowest level in set X2 or X3 at all six ghost atoms. CAM: (minus the energy with CAM-B3LYP in those",
         "functions on a fine grid with exact integrals, the electrons that grid counts). CAMCHECK: one case on a finer grid.",
         'The positions were read from %s, sha256 of job.inp %s' % (SRC, sha(os.path.join(res, S6, SRC, "job.inp"))),
         'and of pc.pc %s."""' % sha(os.path.join(res, S6, SRC, "pc.pc")),
         "SITE = {%s}" % ", ".join("%d: (%s)" % (s, ", ".join(repr(v) for v in site[s])) for s in SITES),
         "GHOST = {%s}" % ", ".join("%d: (%s)" % (s, ", ".join(repr(v) for v in ghost[s])) for s in SITES)]
    lim, lim2, hf, cam = {}, {}, {}, {}
    for k, sites in SYSTEMS:
        ch = [site[s] for s in sites]
        lim[k], lim2[k] = lowest(ch, large(1), ch), lowest(ch, large(2), ch)
        say("%s limit %.6f larger basis %.6f" % (k, lim[k], lim2[k]))
    for name in sorted(ADD):
        for k, sites in SYSTEMS:
            ch = [site[s] for s in sites]
            hf[(name, k)] = lowest([ghost[s] for s in SITES], functions(name), ch)
            t0 = time.time()
            cam[(name, k)] = functional(ghost, name, ch, GRID)
            say("%s %s HF %.6f (by SCF %.6f) CAM %.6f count %.6f  %.0f s" % (name, k, hf[(name, k)], hf_scf(ghost, name, ch), cam[(name, k)][0], cam[(name, k)][1], time.time() - t0))
    chk = functional(ghost, CHECK[0], [site[s] for s in dict(SYSTEMS)[CHECK[1]]], GRID2)
    say("check %s %s finer grid: CAM %.6f count %.6f" % (CHECK[0], CHECK[1], chk[0], chk[1]))
    L.append("LIMIT = {%s}" % ", ".join('"%s": %.6f' % (k, lim[k]) for k, _ in SYSTEMS))
    L.append("LIMIT2 = {%s}" % ", ".join('"%s": %.6f' % (k, lim2[k]) for k, _ in SYSTEMS))
    L.append("HF = {")
    for name in sorted(ADD):
        L.append("    " + " ".join('("%s", "%s"): %.6f,' % (name, k, hf[(name, k)]) for k, _ in SYSTEMS))
    L.append("}")
    L.append("CAM = {")
    for name in sorted(ADD):
        L.append("    " + " ".join('("%s", "%s"): (%.6f, %.6f),' % (name, k, cam[(name, k)][0], cam[(name, k)][1]) for k, _ in SYSTEMS))
    L.append("}")
    L.append('CAMCHECK = {("%s", "%s"): (%.6f, %.6f)}' % (CHECK[0], CHECK[1], chk[0], chk[1]))
    with open(out, "w") as f:
        f.write("\n".join(L) + "\n")
    say("written")
    return 0


if __name__ == "__main__":
    sys.exit(main())
