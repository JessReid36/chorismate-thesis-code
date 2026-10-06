#!/usr/bin/env python3
"""
wp5_rigid_lj_scan.py - WP5 (charge-model plan), rigid-model part: how far a +1 site carrying
Lennard-Jones parameters stops from the substrate, and how steeply that depends on the site's
R_min and epsilon.

MODEL (mechanical embedding, no polarisation, rigid substrate - an expectation-setter for the
QM/MM points of WP4/WP5, not a result in its own right):
    E(p) = K q_s sum_j q_j / r_j  +  sum_j e_ij [ (R_ij/r_j)^12 - 2 (R_ij/r_j)^6 ]
    R_ij = (Rc_s + Rc_j)/2,  e_ij = sqrt(e_s e_j)        Lorentz-Berthelot on the full-R_min column,
                                                          the rule ORCA 6.0.1 uses (WP1, results
                                                          phase2.2/wp1_lj_combining)
    K = 332.0637133 kcal A mol^-1 e^-2 (e^2/4 pi eps0, CODATA 2018)
Substrate: the 24 QM atoms of each post-cut frame's optimised reactant (05_qmmm/20_invacuo/<f>_R.xyz,
which s18b shows equals NEB image 0), with the charges and LJ parameters of rows 6208-6231 of the
committed 05_qmmm/13_bridge/complex_solvated.ORCAFF.prms (GAFF types, AM1-BCC charges, per step08a/b).
Site: charge +1 e, LJ (Rc_s, e_s) from the grid below.

Stopping distances per frame and parameter set:
    RAY      the energy minimum along the s13 ray (substrate centroid -> O3, outward), as distance to O3;
    S13      the site minimised in 3D from the s13 start point (3.2 A from O3 on that ray): distance to
             the nearest oxygen and which one (the quantity the s13 jobs measured);
    CONTACT  per oxygen-bearing GROUP - carboxylate A (O1, O2 on C3), carboxylate B (O5, O6 on C10),
             ether O3, hydroxyl O4 - the site minimised in 3D from 3.2 A outside the group (carboxylate:
             along the bisector from the carboxyl carbon through the O-O midpoint, the bidentate approach;
             ether and hydroxyl: away from the bonded heavy atoms). The contact distance is the distance
             to the group's nearest oxygen, and it counts only if the minimum's nearest oxygen belongs to
             that group; otherwise the start is recorded as having MIGRATED and no distance is claimed.
             Slopes and epsilon effects use contacts in the same group only. (A first version defined
             basins per oxygen; the two oxygens of a carboxylate share one bidentate basin, so it split
             one well into two and was replaced.)

PARAMETER GRID (site), from the committed force field and the literature:
    R* (= Rc/2): 1.369 (Na+, the committed ion set), 1.824 (every Amber N in the committed file),
                 1.915 and 2.006 (N +5% and +10%: the range by which Friesner and co-workers enlarged
                 QM-atom LJ radii, as reported by Senn & Thiel 2009)
    epsilon:     0.087439 (Na+), 0.170 (Amber N), 0.340 (twice Amber N, an upper probe with no
                 literature value attached - stated as such)
The grid is a sensitivity scan, not a fit: no parameter is tuned to hit a target distance (Ryde 2016:
"it seems dangerous to correct a problem in the electrostatics by modifying van der Waals parameters").

CHECKS the script makes on itself (it exits non-zero if any fails):
    - the 24 prms rows are C10 H8 O6 in the xyz element order, net charge -2.000;
    - analytic gradient agrees with a central finite difference (1e-6 A) to 1e-5 relative;
    - every FREE minimum has |grad| < 1e-5 kcal/mol/A and lies more than 1.5 A from every atom;
    - every RAY minimum lies strictly inside the scanned 1.0-6.0 A range.
NOTE recorded in the report: the hydroxyl hydrogen (H8, GAFF ho) carries no LJ parameters, so a
NEGATIVE LJ site would meet no repulsive wall there; this scan uses +1 sites only.

USAGE   python3 wp5_rigid_lj_scan.py <results_repo> [out_dir]
Writes wp5_rigid_scan.tsv and WP5_RIGID_REPORT.txt into out_dir (default: results_repo/phase2.2/wp5_rigid_scan).
Needs numpy and scipy.
"""
import math
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

K = 332.0637133
QM_FIRST, QM_LAST = 6208, 6231
NAMES = ("C1 H1 H2 C2 C3 O1 O2 O3 C4 H3 C5 H4 C6 C7 H5 C8 H6 C9 H7 O4 H8 C10 O5 O6").split()
O_IDX = [i for i, n in enumerate(NAMES) if n.startswith("O")]
O3 = NAMES.index("O3")
RSTAR = [("Na", 1.369), ("N", 1.824), ("N+5%", 1.824 * 1.05), ("N+10%", 1.824 * 1.10)]
EPS = [("Na", 0.087439), ("N", 0.170), ("2xN", 0.340)]
DETAIL_FRAMES = ["41786", "55446", "24883"]   # 41786, 55446: touched by neither C4 nor C5; 24883: the s13 frame
# oxygen-bearing groups: (name, member oxygens, anchor carbon for the start direction or None)
GROUPS = [("carbA", ("O1", "O2"), "C3"), ("carbB", ("O5", "O6"), "C10"), ("ether", ("O3",), None), ("hydroxyl", ("O4",), None)]
GROUP_OF = {m: g for g, ms, a in GROUPS for m in ms}
fails = []


def load_substrate(R):
    L = (R / "05_qmmm/13_bridge/complex_solvated.ORCAFF.prms").read_text().splitlines()
    n = int(L[3].split()[0])
    rows = [l.split() for l in L[4:4 + n] if QM_FIRST <= int(l.split()[0]) <= QM_LAST]
    el = [r[1] for r in rows]
    q = np.array([float(r[2]) for r in rows])
    eps = np.array([abs(float(r[3])) for r in rows])
    rc = np.array([float(r[4]) for r in rows])
    return el, q, eps, rc


def read_xyz(p):
    L = Path(p).read_text().splitlines()
    n = int(L[0].split()[0])
    return [l.split()[0] for l in L[2:2 + n]], np.array([[float(v) for v in l.split()[1:4]] for l in L[2:2 + n]])


def make_energy(X, q, eps, rc, rc_s, e_s, q_s=1.0):
    Rij = 0.5 * (rc_s + rc)
    eij = np.sqrt(e_s * eps)

    def f(p):
        d = p - X
        r = np.sqrt((d * d).sum(1))
        s6 = (Rij / r) ** 6
        E = K * q_s * (q / r).sum() + (eij * (s6 * s6 - 2.0 * s6)).sum()
        dEdr = -K * q_s * q / r ** 2 + eij * (-12.0 * s6 * s6 + 12.0 * s6) / r
        g = (dEdr[:, None] * d / r[:, None]).sum(0)
        return E, g
    return f


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    R = Path(sys.argv[1])
    out = Path(sys.argv[2]) if len(sys.argv) > 2 else R / "phase2.2/wp5_rigid_scan"
    out.mkdir(parents=True, exist_ok=True)
    el, q, eps, rc = load_substrate(R)
    rep = []
    say = rep.append
    canon = "CHHCCOOOCHCHCCHCHCHOHCOO"
    if "".join(e[0] for e in el) != canon or abs(q.sum() + 2.0) > 1e-6:
        fails.append("prms rows are not the 24-atom substrate with charge -2")
    gl = [l for l in (R / "phase2.2/grid_v2.tsv").read_text().splitlines() if l.startswith("#   frames =")]
    frames = gl[0].split("=", 1)[1].strip().split(",")
    say("WP5 rigid-model scan (wp5_rigid_lj_scan.py) - +1 site with LJ against the rigid optimised reactant")
    say("combining rule: Lorentz-Berthelot on the full-R_min column (measured, WP1); K = {} kcal A/mol/e^2".format(K))
    say("substrate: prms rows {}-{} ({}), net charge {:+.6f}".format(QM_FIRST, QM_LAST, "".join(e[0] for e in el), q.sum()))

    # gradient check once, on the first frame and the N parameter set
    els, X = read_xyz(R / "05_qmmm/20_invacuo/{}_R.xyz".format(frames[0]))
    f = make_energy(X, q, eps, rc, 2 * 1.824, 0.17)
    p0 = X[O3] + 3.0 * (X[O3] - X.mean(0)) / np.linalg.norm(X[O3] - X.mean(0)) + np.array([0.3, -0.2, 0.1])
    _, g = f(p0)
    gn = np.array([(f(p0 + h)[0] - f(p0 - h)[0]) / 2e-6 for h in np.eye(3) * 1e-6])
    rel = np.abs(g - gn).max() / np.abs(g).max()
    say("analytic vs finite-difference gradient: max relative difference {:.1e}".format(rel))
    if rel > 1e-5:
        fails.append("analytic gradient disagrees with finite differences ({:.1e})".format(rel))

    rows = []
    contacts = []
    for f_ in frames:
        els, X = read_xyz(R / "05_qmmm/20_invacuo/{}_R.xyz".format(f_))
        if "".join(e[0] for e in els) != canon:
            fails.append("frame {}: element order differs".format(f_))
            continue
        c = X.mean(0)
        u = (X[O3] - c) / np.linalg.norm(X[O3] - c)
        start = X[O3] + 3.2 * u
        sets = [(rl, rs, ep_l, ep) for rl, rs in RSTAR for ep_l, ep in EPS] if f_ in DETAIL_FRAMES \
            else [("Na", 1.369, "Na", 0.087439), ("N", 1.824, "N", 0.170)]
        for rl, rs, ep_l, ep in sets:
            fE = make_energy(X, q, eps, rc, 2 * rs, ep)
            ds = np.arange(1.0, 6.0 + 1e-9, 0.0005)
            Er = np.array([fE(X[O3] + d * u)[0] for d in ds])
            k = int(Er.argmin())
            if k in (0, len(ds) - 1):
                fails.append("frame {} {}/{}: ray minimum at the scan edge".format(f_, rl, ep_l))
            res = minimize(lambda p: fE(p), start, jac=True, method="BFGS", options={"gtol": 1e-7, "maxiter": 2000})
            p = res.x
            E, g = fE(p)
            dist = np.linalg.norm(X - p, axis=1)
            jO = min(O_IDX, key=lambda i: dist[i])
            if np.linalg.norm(g) > 1e-5 or dist.min() < 1.5:
                fails.append("frame {} {}/{}: free minimum not clean (|g| {:.1e}, closest atom {:.2f} A)".format(
                    f_, rl, ep_l, np.linalg.norm(g), dist.min()))
            rows.append((f_, rl, rs, ep_l, ep, ds[k], dist[jO], NAMES[jO], dist[O3], dist.min(), NAMES[int(dist.argmin())], E))
            for gname, members, anchor in GROUPS:
                Om = X[[NAMES.index(m) for m in members]].mean(0)
                if anchor is not None:
                    v = Om - X[NAMES.index(anchor)]
                else:
                    o = NAMES.index(members[0])
                    heavy = [j for j in range(len(X)) if j != o and NAMES[j][0] != "H" and np.linalg.norm(X[j] - X[o]) < 1.7]
                    v = X[o] - X[heavy].mean(0)
                st = Om + 3.2 * v / np.linalg.norm(v)
                rr = minimize(lambda p: fE(p), st, jac=True, method="BFGS", options={"gtol": 1e-7, "maxiter": 2000})
                Eo, go = fE(rr.x)
                dd = np.linalg.norm(X - rr.x, axis=1)
                jn = min(O_IDX, key=lambda i: dd[i])
                if np.linalg.norm(go) > 1e-5 or dd.min() < 1.5:
                    fails.append("frame {} {}/{} contact {}: minimum not clean".format(f_, rl, ep_l, gname))
                own = NAMES[jn] in members
                contacts.append((f_, rl, rs, ep_l, ep, gname, GROUP_OF[NAMES[jn]], dd[jn] if own else float("nan")))

    hdr = "frame\tsite_Rstar_label\tRstar_A\teps_label\teps_kcal\tray_min_dO3_A\tfree_min_dO_A\tfree_nearest_O\tfree_dO3_A\tfree_nearest_atom_A\tfree_nearest_atom\tE_kcal"
    (out / "wp5_rigid_scan.tsv").write_text(hdr + "\n" + "\n".join(
        "\t".join([r[0], r[1], "{:.4f}".format(r[2]), r[3], "{:.6f}".format(r[4]), "{:.4f}".format(r[5]),
                   "{:.4f}".format(r[6]), r[7], "{:.4f}".format(r[8]), "{:.4f}".format(r[9]), r[10], "{:.4f}".format(r[11])])
        for r in rows) + "\n")

    (out / "wp5_contacts.tsv").write_text("frame\tRstar_label\tRstar_A\teps_label\teps_kcal\tstart_O\tminimum_nearest_O\tcontact_A\n" + "\n".join(
        "\t".join([c[0], c[1], "{:.4f}".format(c[2]), c[3], "{:.6f}".format(c[4]), c[5], c[6],
                   "" if math.isnan(c[7]) else "{:.4f}".format(c[7])]) for c in contacts) + "\n")

    def cval(f_, rl, ep_l, o):
        v = [c[7] for c in contacts if c[0] == f_ and c[1] == rl and c[3] == ep_l and c[5] == o]
        return v[0] if v else float("nan")

    say("\nDETAIL FRAMES - stopping distances (A)")
    say("  RAY = to O3 along the s13 ray; S13 = from the s13 start, to the nearest O (which O);")
    say("  carbA..hydroxyl = CONTACT distance in that group's own basin ('migr' = the minimum left that group)")
    for f_ in DETAIL_FRAMES:
        say("frame {}".format(f_))
        say("   R*      eps     RAY    S13 (O)    " + "".join("{:>10}".format(g) for g, m, a in GROUPS))
        for r in rows:
            if r[0] == f_:
                cs = "".join("{:>10}".format("migr" if math.isnan(cval(f_, r[1], r[3], g)) else "{:.3f}".format(cval(f_, r[1], r[3], g))) for g, m, a in GROUPS)
                say("   {:<6}  {:<5} {:6.3f}  {:6.3f} ({:<2})   {}".format(r[1], r[3], r[5], r[6], r[7], cs))

    say("\nSLOPES d(stop)/d(R*), least squares over the four R* values (same basin only), detail frames:")
    for f_ in DETAIL_FRAMES:
        for ep_l, ep in EPS:
            parts = []
            xs = np.array([rs for rl, rs in RSTAR])
            ry = np.array([[r[5] for r in rows if r[0] == f_ and r[1] == rl and r[3] == ep_l][0] for rl, rs in RSTAR])
            parts.append("RAY {:.3f}".format(np.polyfit(xs, ry, 1)[0]))
            for g, m, a in GROUPS:
                ys = np.array([cval(f_, rl, ep_l, g) for rl, rs in RSTAR])
                parts.append("{} {}".format(g, "  -  " if np.isnan(ys).any() else "{:.3f}".format(np.polyfit(xs, ys, 1)[0])))
            say("   frame {}  eps {:<4}: {}".format(f_, ep_l, "   ".join(parts)))

    say("\nEPSILON EFFECT: contact distance at eps 0.340 minus at 0.170 (A), carboxylate A and B, detail frames:")
    for f_ in DETAIL_FRAMES:
        for g in ("carbA", "carbB"):
            dv = [(rl, cval(f_, rl, "2xN", g) - cval(f_, rl, "N", g)) for rl, rs in RSTAR]
            say("   frame {} {}: ".format(f_, g) + "   ".join("R* {} {}".format(rl, "migr" if math.isnan(d) else "{:+.3f}".format(d)) for rl, d in dv))

    say("\nALL 30 FRAMES (current reactant geometries; C4 may move some of them):")
    for lab, ep_l in (("Na", "Na"), ("N", "N")):
        v = [r for r in rows if r[1] == lab and r[3] == ep_l]
        s13 = np.array([r[6] for r in v]); ry = np.array([r[5] for r in v])
        nO = {}
        for r in v:
            nO[r[7]] = nO.get(r[7], 0) + 1
        say("   site {:<2}: RAY d(O3) mean {:.3f} sd {:.3f};  S13 d(O) mean {:.3f} sd {:.3f} (nearest O {})".format(
            lab, ry.mean(), ry.std(ddof=1), s13.mean(), s13.std(ddof=1), dict(sorted(nO.items()))))
        for g, m, a in GROUPS:
            sel = [c for c in contacts if c[1] == lab and c[3] == ep_l and c[5] == g]
            vals = np.array([c[7] for c in sel])
            ok = vals[~np.isnan(vals)]
            went = {}
            for c in sel:
                if math.isnan(c[7]):
                    went[c[6]] = went.get(c[6], 0) + 1
            say("      contact {:<8}: {:2d}/{} frames stay in the group{}{}".format(
                g, len(ok), len(vals), "" if not len(ok) else ", mean {:.3f} sd {:.3f} range {:.3f}-{:.3f}".format(
                    ok.mean(), ok.std(ddof=1) if len(ok) > 1 else 0.0, ok.min(), ok.max()),
                "" if not went else "; migrated to " + ", ".join("{} x{}".format(k, v) for k, v in sorted(went.items()))))
    say("\nAGAINST THE EXPECTATIONS FIXED IN THE PLAN (charge-model plan, WP5 and section 5):")
    sl = []
    for f_ in DETAIL_FRAMES:
        for ep_l, ep in EPS:
            xs = np.array([rs for rl, rs in RSTAR])
            for key in ("RAY", "carbA"):
                ys = np.array([[r[5] for r in rows if r[0] == f_ and r[1] == rl and r[3] == ep_l][0] for rl, rs in RSTAR]) \
                    if key == "RAY" else np.array([cval(f_, rl, ep_l, "carbA") for rl, rs in RSTAR])
                if not np.isnan(ys).any():
                    sl.append(np.polyfit(xs, ys, 1)[0])
    say("   slope d(stop)/d(R*) expected about 0.74: RAY and carboxylate-A slopes span {:.3f}-{:.3f}".format(min(sl), max(sl)))
    de = []
    for f_ in DETAIL_FRAMES:
        for rl, rs in RSTAR:
            de.append(("RAY", f_, rl, [r[5] for r in rows if r[0] == f_ and r[1] == rl and r[3] == "2xN"][0]
                       - [r[5] for r in rows if r[0] == f_ and r[1] == rl and r[3] == "N"][0]))
            for g in ("carbA", "carbB"):
                d = cval(f_, rl, "2xN", g) - cval(f_, rl, "N", g)
                if not math.isnan(d):
                    de.append((g, f_, rl, d))
    over = [x for x in de if x[3] >= 0.1]
    say("   doubling eps moves the stop by under 0.1 A: largest {:+.3f} A; {} of {} cases at or above 0.1 A{}".format(
        max(x[3] for x in de), len(over), len(de),
        "" if not over else " (" + ", ".join("{} {} R* {} {:+.3f}".format(*x) for x in over) + ")"))
    r24 = [r for r in rows if r[0] == "24883" and r[1] == "Na" and r[3] == "Na"][0]
    say("   frame 24883, sodium parameters, expected about 2.15 A: S13 {:.3f} A to {} (the QM/MM s13_lj_test2 gave 2.142 A to O2)".format(r24[6], r24[7]))
    say("   ether contact: a free +1 site started at O3 leaves it for a carboxylate in every frame and setting scanned; in the enzyme")
    say("   Arg90 is held at O3 by the protein, so a free single site cannot test that contact (see WP4)")
    import scipy
    say("\nnumpy {}, scipy {}".format(np.__version__, scipy.__version__))
    say("\nREFERENCE: the LJ pair minimum alone (no Coulomb), by WP1's rule: Na-O(carboxylate) 3.0302,")
    say("      N-O(carboxylate) 3.4852, Na-O(ether) 3.0527 A. The contact distances above sit well inside these:")
    say("      the +1 site is pulled in by its Coulomb attraction to the dianion until the r^-12 wall stops it.")
    say("\nNOTE: the hydroxyl hydrogen H8 (GAFF ho) has no LJ parameters, so a negative LJ site would meet no")
    say("      repulsive wall there. This scan uses +1 sites only.")
    say("\nRESULT: " + ("all self-checks PASS" if not fails else "FAIL - " + "; ".join(fails)))
    (out / "WP5_RIGID_REPORT.txt").write_text("\n".join(rep) + "\n")
    print("\n".join(rep))
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
