#!/usr/bin/env python3
"""
wp4_esp_compare.py - does the AM1-BCC methylguanidinium (WP4 validation group, decision D4) present
the same electrostatics as the enzyme's own arginine?

Individual atomic charges are not unique: different charge sets can give nearly the same potential
outside a molecule, and the potential is what a neighbouring atom feels. So the two charge sets are
compared by the electrostatic potential (ESP) they produce, on the same geometry:
    A  the AM1-BCC charges from wp4a (charges_mgn_am1bcc.dat, atom order of mgn_input.mol2)
    B  Amber's own arginine charges for Arg90 (residue 217), read from the committed
       05_qmmm/13_bridge/complex_solvated.ORCAFF.prms (atoms CD..HH22). Amber's CD..HH22 sum to
       +0.8755 because the rest of the residue carries the remainder; the methyl cap HD1 is given
       that remainder so that B is also +1 (stated as a convention, not a force-field value).

Points: four shells at 1.4, 1.6, 1.8 and 2.0 times the Bondi van der Waals radii (C 1.70, N 1.55,
H 1.20 A) - the shell scales of Merz-Kollman ESP fitting, with Bondi's radii rather than theirs -
200 random directions per atom per shell (fixed seed),
keeping only points outside every atom's sphere at that scale. Also: the potential at the five
hydrogen-bond acceptor positions, 1.9 A beyond each N-H hydrogen along the N-H bond.
Potential in kcal/mol per unit charge, K = 332.0637133, no dielectric.

USAGE   python3 wp4_esp_compare.py <results_repo> <mgn_input.mol2> <charges_mgn_am1bcc.dat> [report]
Needs numpy.
"""
import sys
from pathlib import Path

import numpy as np

K = 332.0637133
ORDER = ["CD", "HD1", "HD2", "HD3", "NE", "HE", "CZ", "NH1", "HH11", "HH12", "NH2", "HH21", "HH22"]
VDW = {"C": 1.70, "N": 1.55, "H": 1.20}


def main():
    if len(sys.argv) < 4:
        sys.exit(__doc__)
    R, mol2, chg = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
    rep_path = Path(sys.argv[4]) if len(sys.argv) > 4 else None
    blk = mol2.read_text().split("@<TRIPOS>ATOM\n")[1].split("@<TRIPOS>BOND\n")[0].splitlines()
    names = [l.split()[1] for l in blk]
    if names != ORDER:
        sys.exit("FAIL atom order in {} is {}, expected {}".format(mol2, names, ORDER))
    X = np.array([[float(v) for v in l.split()[2:5]] for l in blk])
    qa = np.array([float(x) for x in chg.read_text().split()])
    if len(qa) != 13 or abs(qa.sum() - 1.0) > 1e-6:
        sys.exit("FAIL AM1-BCC charge file: {} values summing to {:+.6f}".format(len(qa), qa.sum()))

    # Amber's Arg90 charges from the committed force field: PDB serials 3498-3509 = prms rows CD..HH22
    L = (R / "05_qmmm/13_bridge/complex_solvated.ORCAFF.prms").read_text().splitlines()
    rows = {int(l.split()[0]): l.split() for l in L[4:4 + int(L[3].split()[0])]}
    pdbnames = {}
    for l in (R / "05_qmmm/18e_reactant_reduced/reactant_reduced.pdb").read_text().splitlines():
        if l[:6] in ("ATOM  ", "HETATM") and l[17:20] == "ARG" and int(l[22:26]) == 217:
            pdbnames[l[12:16].strip()] = int(l[6:11])
    amber = {n: float(rows[pdbnames[n]][2]) for n in ORDER if n != "HD1"}
    amber["HD1"] = 1.0 - sum(amber.values())
    qb = np.array([amber[n] for n in ORDER])

    rad = np.array([VDW[n[0]] for n in ORDER])
    rng = np.random.default_rng(20261006)
    P = []
    for s in (1.4, 1.6, 1.8, 2.0):
        for i in range(len(X)):
            v = rng.normal(size=(200, 3))
            v /= np.linalg.norm(v, axis=1)[:, None]
            pts = X[i] + s * rad[i] * v
            d = np.linalg.norm(pts[:, None, :] - X[None], axis=2)
            P.append(pts[(d >= s * rad[None, :] - 1e-9).all(1)])
    P = np.vstack(P)

    def esp(q, pts):
        return K * (q[None, :] / np.linalg.norm(pts[:, None, :] - X[None], axis=2)).sum(1)

    Va, Vb = esp(qa, P), esp(qb, P)
    rms_d = float(np.sqrt(((Va - Vb) ** 2).mean()))
    rms_b = float(np.sqrt((Vb ** 2).mean()))
    out = ["WP4 ESP comparison (wp4_esp_compare.py): AM1-BCC methylguanidinium against Amber's Arg90 charges",
           "", "charges (e):   atom   AM1-BCC    Amber Arg90"]
    for n, a, b in zip(ORDER, qa, qb):
        out.append("               {:<5} {:+.4f}    {:+.4f}{}".format(n, a, b, "  (cap: remainder to +1)" if n == "HD1" else ""))
    out += ["", "ESP on {} points (shells at 1.4-2.0 x Bondi radii, seed 20261006):".format(len(P)),
            "  RMS difference {:.2f} kcal/mol/e against an RMS potential of {:.2f}: {:.1f}%; mean difference {:+.2f}".format(
                rms_d, rms_b, 100 * rms_d / rms_b, float((Va - Vb).mean())),
            "", "potential at the hydrogen-bond acceptor positions (1.9 A beyond each N-H hydrogen):"]
    ratios = []
    for h, n in (("HE", "NE"), ("HH11", "NH1"), ("HH12", "NH1"), ("HH21", "NH2"), ("HH22", "NH2")):
        u = X[ORDER.index(h)] - X[ORDER.index(n)]
        p = (X[ORDER.index(h)] + 1.9 * u / np.linalg.norm(u))[None]
        a, b = esp(qa, p)[0], esp(qb, p)[0]
        ratios.append(a / b)
        out.append("  beyond {:<4}: AM1-BCC {:6.1f}  Amber {:6.1f} kcal/mol/e  ratio {:.3f}".format(h, a, b, a / b))
    out += ["", "RESULT: the AM1-BCC group's potential differs from Amber's arginine by {:.1f}% RMS overall and by "
            "{:.0f}-{:.0f}% at the acceptor positions".format(100 * rms_d / rms_b, 100 * min(abs(1 - r) for r in ratios),
                                                            100 * max(abs(1 - r) for r in ratios))]
    text = "\n".join(out) + "\n"
    if rep_path:
        rep_path.write_text(text)
    print(text, end="")
    return 0


if __name__ == "__main__":
    sys.exit(main())
