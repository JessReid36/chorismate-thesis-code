#!/usr/bin/env python3
"""
wp4_analyse.py - WP4 analysis: applies the criteria fixed before submission (WP4_CRITERIA.txt and
WP4_CRITERIA_CARTESIAN.txt, committed in 86e8f43) to the nine QM/MM optimisations in wp4_s13_rerun/.

Reads, per folder: combined.xyz (start), job.xyz (ORCA's final geometry), job.out (convergence),
site.ORCAFF.prms (the group's charges, arm c), and job_trj.xyz only if a run did not converge.
Atom order: substrate 0-23 (C1 H1 H2 C2 C3 O1 O2 O3 C4 H3 C5 H4 C6 C7 H5 C8 H6 C9 H7 O4 H8 C10 O5 O6),
then the site (arms a, b: atom 24) or methylguanidinium (arm c: 24-36, CD HD1 HD2 HD3 NE HE CZ NH1 HH11
HH12 NH2 HH21 HH22).

Definitions used (fixed here, before reading any result):
  hydrogen bond   N-H...O with H...O <= 2.5 A and N...O <= 3.5 A (a purely geometric definition, fixed
                  here; the D-WP4 windows constrain N...O more tightly where a verdict depends on it)
  bridge kept     arm c has a hydrogen bond to O3 AND one to a carboxylate oxygen (O1 O2 O5 O6)
  collapse        arm a: the site's closest substrate oxygen below 1.2 A
  integrity       arms b, c: every one of the substrate's 24 bonds within 0.15 A of its start length
  charge centre   sum(q r)/sum(q) over the group, with its force-field charges
Writes WP4_REPORT.txt and wp4_results.tsv into the folder; exits 0 (the criteria are reported, not
enforced by exit code), or 1 if a file is missing or malformed.

USAGE   python3 wp4_analyse.py [wp4_s13_rerun folder]
Python 3.6, standard library only (runs on the hpc1 login node).
"""
import math
import re
import sys
from pathlib import Path

W = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("/home/18660916/system_development/phase2.2/wp4_s13_rerun")
NAMES = "C1 H1 H2 C2 C3 O1 O2 O3 C4 H3 C5 H4 C6 C7 H5 C8 H6 C9 H7 O4 H8 C10 O5 O6".split()
GROUP = "CD HD1 HD2 HD3 NE HE CZ NH1 HH11 HH12 NH2 HH21 HH22".split()
BONDS = [(0, 1), (0, 2), (0, 3), (3, 4), (4, 5), (4, 6), (3, 7), (7, 8), (8, 9), (8, 10), (10, 11), (10, 12),
         (12, 13), (13, 14), (13, 15), (15, 16), (15, 17), (17, 18), (17, 19), (19, 20), (17, 8), (12, 21),
         (21, 22), (21, 23)]
OXY = [5, 6, 7, 19, 22, 23]                       # O1 O2 O3 O4 O5 O6
CARBOX = [5, 6, 22, 23]
DONORS = [("NE", "HE"), ("NH1", "HH11"), ("NH1", "HH12"), ("NH2", "HH21"), ("NH2", "HH22")]
ENZ = {"NE..O2": 2.71, "NH2..O3": 2.83, "CZ..O3": 3.55}      # the enzyme's QM/MM contacts at frame 41786
FOLDERS = ["arma_anch", "arma_free", "arma_freecart", "armb_anch", "armb_free", "armb_freecart",
           "armc_anch", "armc_free", "armc_freecart"]


def read_xyz(p):
    L = p.read_text().splitlines()
    n = int(L[0].split()[0])
    return [l.split()[0] for l in L[2:2 + n]], [[float(v) for v in l.split()[1:4]] for l in L[2:2 + n]]


def read_traj(p):
    L, out, i = p.read_text().splitlines(), [], 0
    while i < len(L):
        if not L[i].strip():
            i += 1
            continue
        n = int(L[i].split()[0])
        out.append([[float(v) for v in l.split()[1:4]] for l in L[i + 2:i + 2 + n]])
        i += n + 2
    return out


def d(a, b):
    return math.sqrt(sum((a[k] - b[k]) ** 2 for k in range(3)))


def centroid(X):
    return [sum(x[k] for x in X) / len(X) for k in range(3)]


def superpose(A, B):
    """rotation angle (deg) and RMSD after optimal superposition of B onto A (Horn's quaternion method)"""
    n = len(A)
    ca, cb = centroid(A), centroid(B)
    a = [[p[k] - ca[k] for k in range(3)] for p in A]
    b = [[p[k] - cb[k] for k in range(3)] for p in B]
    S = [[sum(a[i][r] * b[i][c] for i in range(n)) for c in range(3)] for r in range(3)]
    (xx, xy, xz), (yx, yy, yz), (zx, zy, zz) = S
    N = [[xx + yy + zz, yz - zy, zx - xz, xy - yx], [yz - zy, xx - yy - zz, xy + yx, zx + xz],
         [zx - xz, xy + yx, -xx + yy - zz, yz + zy], [xy - yx, zx + xz, yz + zy, -xx - yy + zz]]
    lam = max(abs(v) for row in N for v in row) * 4 + 1.0
    v = [1.0, 0.0, 0.0, 0.0]
    for _ in range(2000):
        w = [sum((N[i][j] + (lam if i == j else 0.0)) * v[j] for j in range(4)) for i in range(4)]
        s = math.sqrt(sum(x * x for x in w))
        v = [x / s for x in w]
    emax = sum(v[i] * sum(N[i][j] * v[j] for j in range(4)) for i in range(4))
    ga, gb = sum(x * x for p in a for x in p), sum(x * x for p in b for x in p)
    rmsd = math.sqrt(max(0.0, (ga + gb - 2.0 * emax) / n))
    angle = math.degrees(2.0 * math.acos(min(1.0, abs(v[0]))))
    return angle, rmsd


def converged(p):
    t = p.read_text(errors="replace")
    return "THE OPTIMIZATION HAS CONVERGED" in t, "ORCA TERMINATED NORMALLY" in t, t.count("GEOMETRY OPTIMIZATION CYCLE")


def group_charges(prms):
    L = prms.read_text().splitlines()
    i = L.index("$atoms")
    n = int(L[i + 1].split()[0])
    rows = [l.split() for l in L[i + 2:i + 2 + n]]
    return [float(r[2]) for r in rows[24:37]]


def main():
    rep, tsv = [], []
    say = rep.append
    say("WP4 report (wp4_analyse.py) - criteria from WP4_CRITERIA.txt and WP4_CRITERIA_CARTESIAN.txt (commit 86e8f43)")
    say("frame 41786; enzyme's own QM/MM contacts at the start: NE..O2 2.71, NH2..O3 2.83, CZ..O3 3.55 A\n")
    for f in FOLDERS:
        F = W / f
        for need in ("combined.xyz", "job.xyz", "job.out"):
            if not (F / need).exists():
                sys.exit("FAIL {}: missing {}".format(f, need))
        el0, X0 = read_xyz(F / "combined.xyz")
        el1, X1 = read_xyz(F / "job.xyz")
        if el0 != el1 or len(X0) not in (25, 37):
            sys.exit("FAIL {}: start and final geometries differ in atoms".format(f))
        for i, j in BONDS:
            if d(X0[i], X0[j]) > 1.6:
                sys.exit("FAIL {}: bond list does not match the start geometry ({}-{})".format(f, NAMES[i], NAMES[j]))
        conv, term, ncyc = converged(F / "job.out")
        S0, S1 = X0[:24], X1[:24]
        shift = d(centroid(S0), centroid(S1))
        ang, rms = superpose(S0, S1)
        dbond = max(abs(d(S1[i], S1[j]) - d(S0[i], S0[j])) for i, j in BONDS)
        wb = max(BONDS, key=lambda b: abs(d(S1[b[0]], S1[b[1]]) - d(S0[b[0]], S0[b[1]])))
        c16 = (d(S0[0], S0[12]), d(S1[0], S1[12]))
        arm, var = f[3], f.split("_", 1)[1]
        say("== {}  ({} cycles, {}{})".format(f, ncyc, "converged" if conv else "NOT converged",
                                                 "" if term else ", NOT terminated normally"))
        say("   substrate: centroid moved {:.3f} A, rotated {:.1f} deg, internal RMSD {:.3f} A; C1..C6 {:.2f} -> {:.2f} A".format(
            shift, ang, rms, c16[0], c16[1]))
        say("   largest bond change {} {:+.3f} A".format(NAMES[wb[0]] + "-" + NAMES[wb[1]],
                                                       d(S1[wb[0]], S1[wb[1]]) - d(S0[wb[0]], S0[wb[1]])))
        if not conv:
            # s13's rule for a run that stopped: the contact distance moves by at most 0.05 A over the
            # last third of the cycles (contact: site to nearest O for arms a/b, nearest N..O for arm c)
            T = read_traj(F / "job_trj.xyz") if (F / "job_trj.xyz").exists() else []
            def contact(X):
                if arm in "ab":
                    return min(d(X[24], X[i]) for i in OXY)
                return min(d(X[24 + GROUP.index(n)], X[i]) for n in ("NE", "NH1", "NH2") for i in OXY)
            if len(T) >= 3:
                tail = [contact(X) for X in T[-max(1, len(T) // 3):]]
                say("   stopped early: contact moved {:.3f} A over the last {} of {} cycles -> stability rule {}".format(
                    max(tail) - min(tail), len(tail), len(T), "MET" if max(tail) - min(tail) <= 0.05 else "NOT MET"))
            else:
                say("   stopped early and no usable trajectory: result not assessable")
        row = {"run": f, "cycles": ncyc, "converged": conv, "centroid_shift": shift, "rotation": ang,
               "internal_rmsd": rms, "C1C6_start": c16[0], "C1C6_final": c16[1], "max_bond_change": dbond}
        if arm in "ab":
            site = X1[24]
            do = sorted((d(site, S1[i]), NAMES[i]) for i in OXY)
            da = min((d(site, S1[i]), NAMES[i]) for i in range(24))
            say("   site to oxygens: " + ", ".join("{} {:.2f}".format(n, x) for x, n in do)
                + "; closest atom {} {:.2f} A".format(da[1], da[0]))
            row.update(nearest_O=do[0][1], d_nearest_O=do[0][0], d_O3=[x for x, n in do if n == "O3"][0])
            if arm == "a":
                co = {n: d(S1[i], S1[j]) for n, i, j in (("C3-O1", 4, 5), ("C3-O2", 4, 6), ("C10-O5", 21, 22),
                                                           ("C10-O6", 21, 23), ("C4-O3", 8, 7), ("C9-O4", 17, 19))}
                say("   C-O bonds: " + ", ".join("{} {:.3f}".format(k, v) for k, v in co.items()))
                ok = do[0][0] < 1.2
                say("   CRITERION collapse (closest O < 1.2 A): {} - {} at {:.3f} A".format(
                    "as expected" if ok else "NOT as expected", do[0][1], do[0][0]))
                row["verdict"] = "collapse" if ok else "no collapse"
            else:
                integ = dbond <= 0.15
                if var == "anch":
                    say("   CRITERION (reported, no pass/fail): stop between the LJ stop (~2.6) and the start; "
                        "nearest O {} {:.3f} A".format(do[0][1], do[0][0]))
                    row["verdict"] = "reported"
                else:
                    ok = do[0][1] in ("O2", "O3") and 2.28 <= do[0][0] <= 2.88
                    say("   CRITERION stop 2.28-2.88 A to O2 or O3 (rigid 2.584 to O2): {} - {} {:.3f} A{}".format(
                        "PASS" if ok else "OUTSIDE", do[0][1], do[0][0],
                        "  [internal coordinates: frame held by ORCA's coordinate choice]" if var == "free" else ""))
                    row["verdict"] = "pass" if ok else "outside"
                say("   integrity (bonds within 0.15 A): {}".format("PASS" if integ else "FAIL"))
                row["integrity"] = integ
        else:
            G0, G1 = X0[24:37], X1[24:37]
            gi = {n: i for i, n in enumerate(GROUP)}
            hb = []
            for nN, nH in DONORS:
                for o in OXY:
                    dh, dn = d(G1[gi[nH]], S1[o]), d(G1[gi[nN]], S1[o])
                    if dh <= 2.5 and dn <= 3.5:
                        hb.append((nN, nH, NAMES[o], dn, dh))
            for h in hb:
                say("   H-bond {}-{}...{}: N..O {:.2f}, H..O {:.2f} A".format(*h))
            if not hb:
                say("   no N-H...O hydrogen bond (H..O <= 2.5, N..O <= 3.5 A)")
            to_O3 = any(h[2] == "O3" for h in hb)
            to_carb = any(h[2] in ("O1", "O2", "O5", "O6") for h in hb)
            nO = min((d(G1[gi[n]], S1[o]), n, NAMES[o]) for n in ("NE", "NH1", "NH2") for o in OXY)
            cz = min((d(G1[gi["CZ"]], S1[o]), NAMES[o]) for o in OXY)
            q = group_charges(F / "site.ORCAFF.prms")
            cc0 = [sum(qq * g[k] for qq, g in zip(q, G0)) / sum(q) for k in range(3)]
            cc1 = [sum(qq * g[k] for qq, g in zip(q, G1)) / sum(q) for k in range(3)]
            ccmin = min((d(cc1, S1[i]), NAMES[i]) for i in range(24))
            dev = {"NE..O2": d(G1[gi["NE"]], S1[6]), "NH2..O3": d(G1[gi["NH2"]], S1[7]), "CZ..O3": d(G1[gi["CZ"]], S1[7])}
            say("   nearest N..O {} {} {:.2f} A; CZ..nearest O {} {:.2f} A".format(nO[1], nO[2], nO[0], cz[1], cz[0]))
            say("   enzyme contacts now: " + ", ".join("{} {:.2f} ({:+.2f} from the enzyme's {:.2f})".format(
                k, v, v - ENZ[k], ENZ[k]) for k, v in dev.items()))
            say("   charge centre moved {:.2f} A; closest approach of the charge centre to any substrate atom: {} {:.2f} A".format(
                d(cc0, cc1), ccmin[1], ccmin[0]))
            say("   group: CZ moved {:.2f} A (CD and HD1 held)".format(d(G0[gi["CZ"]], G1[gi["CZ"]])))
            integ = dbond <= 0.15
            mode = "+".join(sorted({h[2] for h in hb})) or "none"
            win = 2.62 <= nO[0] <= 3.11 and 3.18 <= cz[0] <= 3.86
            if var == "anch":
                ok = to_O3 and to_carb and win
                say("   CRITERION (validation): bridge kept {}; N..O in 2.62-3.11 {}; CZ..O in 3.18-3.86 {} -> {}".format(
                    "yes" if (to_O3 and to_carb) else "NO", "yes" if 2.62 <= nO[0] <= 3.11 else "NO",
                    "yes" if 3.18 <= cz[0] <= 3.86 else "NO", "PASS" if ok else "FAIL"))
                row["verdict"] = "pass" if ok else "fail"
            else:
                say("   informative: binding mode {}; contacts {} the D-WP4 windows".format(mode, "within" if win else "outside"))
                row["verdict"] = "mode " + mode
            say("   integrity (bonds within 0.15 A): {}".format("PASS" if integ else "FAIL"))
            row.update(integrity=integ, mode=mode, d_nearest_NO=nO[0], d_CZ_O=cz[0], cc_moved=d(cc0, cc1),
                       cc_closest=ccmin[0], cc_closest_atom=ccmin[1], **{k.replace("..", "_"): v for k, v in dev.items()})
        say("")
        tsv.append(row)
    keys = []
    for r in tsv:
        for k in r:
            if k not in keys:
                keys.append(k)
    lines = ["\t".join(keys)]
    for r in tsv:
        lines.append("\t".join(("{:.4f}".format(r[k]) if isinstance(r.get(k), float) else str(r.get(k, ""))) for k in keys))
    (W / "wp4_results.tsv").write_text("\n".join(lines) + "\n")
    (W / "WP4_REPORT.txt").write_text("\n".join(rep) + "\n")
    print("\n".join(rep))
    return 0


if __name__ == "__main__":
    sys.exit(main())
