#!/usr/bin/env python3
"""
cs2_contact_scan.py - Stage 2 of the charged-sphere development: contact distance as a function of the sphere's
charge q and size, in the rigid model of WP5 (mechanical embedding, rigid substrate, no polarisation). The energy
model, substrate charges and Lennard-Jones parameters, combining rule and minimiser are WP5's, imported unchanged
from wp5_rigid_lj_scan.py; Stage 1 showed this Lennard-Jones form is ORCA's to 2e-10 Eh.

Scope (STAGE2_PROTOCOL.txt):
  frames   the 30 post-cut frames of grid_v2 (read from its header, as WP5), each at its in vacuo R and TS geometry
           (05_qmmm/20_invacuo/<f>_R.xyz, <f>_TS.xyz)
  sizes    N (R* 1.824 A, eps 0.170), N+5% (1.9152), N+10% (2.0064) [Murphy et al. 2000 range], OPLS-AA sp3 carbon
           (sigma 3.500 A, eps 0.066; R* = 2^(1/6) sigma / 2 = 1.9644), Behrens' stated sigma 3.55 A (eps 0.066,
           R* 1.9924)
  q > 0    +0.25, +0.5, +0.75, +1 against the four oxygen groups, started and classified exactly as WP5's CONTACT
           (carboxylate A O1/O2, carboxylate B O5/O6, ether O3, hydroxyl O4; 'own' only if the minimum's nearest
           oxygen belongs to the group, otherwise 'migrated')
  q < 0    -0.25, -0.5, -0.75, -1 against each hydrogen (H1..H8), started 3.2 A from the hydrogen along heavy atom ->
           hydrogen; outcome 'bound' (|g| < 1e-5 kcal/mol/A, nearest atom within 6 A), 'escaped' (nearest atom
           beyond 6 A: no contact minimum) or 'unconverged' (reported, not used)
  H8 flag  any bound minimum with a negative sphere closer than 2.0 A to H8 (the zero-LJ hydroxyl hydrogen)
Self-checks (exit non-zero if any fails): substrate rows and charge as WP5; analytic vs finite-difference gradient
for a positive and a negative sphere; every positive-q minimum clean (|g| < 1e-5, nearest atom > 1.5 A).
USAGE  python3 cs2_contact_scan.py <results repo> [out dir]   (default out: <repo>/phase2.2/cs_stage2)
Needs numpy and scipy.
"""
import importlib.util, math, sys
from pathlib import Path
import numpy as np
from scipy.optimize import minimize

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("wp5", HERE / "wp5_rigid_lj_scan.py")
wp5 = importlib.util.module_from_spec(spec); spec.loader.exec_module(wp5)
NAMES, GROUPS, O_IDX = wp5.NAMES, wp5.GROUPS, wp5.O_IDX
SIZES = [("N", 1.824, 0.170), ("N+5%", 1.824 * 1.05, 0.170), ("N+10%", 1.824 * 1.10, 0.170),
         ("OPLS_sp3C", 2 ** (1 / 6) * 3.50 / 2, 0.066), ("Behrens_3.55", 2 ** (1 / 6) * 3.55 / 2, 0.066)]
QPOS = [0.25, 0.5, 0.75, 1.0]
QNEG = [-0.25, -0.5, -0.75, -1.0]
H_IDX = [i for i, n in enumerate(NAMES) if n.startswith("H")]
H8 = NAMES.index("H8")
OPT = {"gtol": 1e-7, "maxiter": 5000}
fails = []

def h_start(X, h):
    heavy = min((j for j in range(len(X)) if NAMES[j][0] != "H"), key=lambda j: np.linalg.norm(X[j] - X[h]))
    v = X[h] - X[heavy]
    return X[h] + 3.2 * v / np.linalg.norm(v)

def o_start(X, members, anchor):
    Om = X[[NAMES.index(m) for m in members]].mean(0)
    if anchor is not None:
        v = Om - X[NAMES.index(anchor)]
    else:
        o = NAMES.index(members[0])
        heavy = [j for j in range(len(X)) if j != o and NAMES[j][0] != "H" and np.linalg.norm(X[j] - X[o]) < 1.7]
        v = X[o] - X[heavy].mean(0)
    return Om + 3.2 * v / np.linalg.norm(v)

def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    R = Path(sys.argv[1])
    out = Path(sys.argv[2]) if len(sys.argv) > 2 else R / "phase2.2/cs_stage2"
    out.mkdir(parents=True, exist_ok=True)
    el, q, eps, rc = wp5.load_substrate(R)
    if "".join(e[0] for e in el) != "CHHCCOOOCHCHCCHCHCHOHCOO" or abs(q.sum() + 2.0) > 1e-6:
        fails.append("substrate rows are not the 24-atom dianion")
    gl = [l for l in (R / "phase2.2/grid_v2.tsv").read_text().splitlines() if l.startswith("#   frames =")]
    frames = gl[0].split("=", 1)[1].strip().split(",")
    if len(frames) != 30:
        fails.append("expected 30 frames in grid_v2, found {}".format(len(frames)))
    # gradient checks, first frame R, N sphere, q = +1 and -1
    _, X = wp5.read_xyz(R / "05_qmmm/20_invacuo/{}_R.xyz".format(frames[0]))
    for qs in (1.0, -1.0):
        f = wp5.make_energy(X, q, eps, rc, 2 * 1.824, 0.17, q_s=qs)
        p0 = X.mean(0) + np.array([4.1, -3.3, 2.2])
        _, g = f(p0)
        gn = np.array([(f(p0 + h)[0] - f(p0 - h)[0]) / 2e-6 for h in np.eye(3) * 1e-6])
        rel = np.abs(g - gn).max() / np.abs(g).max()
        if rel > 1e-5:
            fails.append("gradient check q={:+.0f}: {:.1e}".format(qs, rel))
    rows = []
    for fr in frames:
        for st in ("R", "TS"):
            els, X = wp5.read_xyz(R / "05_qmmm/20_invacuo/{}_{}.xyz".format(fr, st))
            if "".join(e[0] for e in els) != "CHHCCOOOCHCHCCHCHCHOHCOO":
                fails.append("frame {} {}: element order".format(fr, st)); continue
            for lab, rs, ep in SIZES:
                for qs in QPOS:
                    fE = wp5.make_energy(X, q, eps, rc, 2 * rs, ep, q_s=qs)
                    for gname, members, anchor in GROUPS:
                        r = minimize(lambda p: fE(p), o_start(X, members, anchor), jac=True, method="BFGS", options=OPT)
                        _, g = fE(r.x)
                        dd = np.linalg.norm(X - r.x, axis=1)
                        jn = min(O_IDX, key=lambda i: dd[i])
                        if np.linalg.norm(g) > 1e-5 or dd.min() < 1.5:
                            fails.append("{} {} {} q={} {}: minimum not clean".format(fr, st, lab, qs, gname))
                        own = NAMES[jn] in members
                        ja = int(dd.argmin())
                        rows.append((fr, st, lab, qs, gname, "own" if own else "migrated", NAMES[jn], dd[jn],
                                     NAMES[ja], dd[ja], dd[H8]))
                for qs in QNEG:
                    fE = wp5.make_energy(X, q, eps, rc, 2 * rs, ep, q_s=qs)
                    for h in H_IDX:
                        r = minimize(lambda p: fE(p), h_start(X, h), jac=True, method="BFGS", options=OPT)
                        _, g = fE(r.x)
                        dd = np.linalg.norm(X - r.x, axis=1)
                        ja = int(dd.argmin())
                        if dd.min() > 6.0:
                            kind = "escaped"
                        elif np.linalg.norm(g) > 1e-5:
                            kind = "unconverged"
                        else:
                            kind = "bound_own" if ja == h else "bound_migrated"
                        rows.append((fr, st, lab, qs, NAMES[h], kind, NAMES[ja], dd[ja], NAMES[ja], dd[ja], dd[H8]))
    hdr = "frame\tstate\tsize\tq\tstart\toutcome\tcontact_atom\tcontact_A\tnearest_atom\tnearest_A\tdist_H8_A"
    (out / "cs2_contacts.tsv").write_text(hdr + "\n" + "\n".join("\t".join([r[0], r[1], r[2], "{:+.2f}".format(r[3]), r[4], r[5],
        r[6], "{:.4f}".format(r[7]), r[8], "{:.4f}".format(r[9]), "{:.4f}".format(r[10])]) for r in rows) + "\n")
    rep = []
    say = rep.append
    say("Stage 2 report (cs2_contact_scan.py) - rigid-model contact distance by charge and size; protocol STAGE2_PROTOCOL.txt")
    say("frames {} x states R, TS = {} geometries; rows {}".format(len(frames), 2 * len(frames), len(rows)))
    say("self-checks: {}".format("all pass" if not fails else "{} FAILED".format(len(fails))))
    for lab, rs, ep in SIZES:
        say("\nSIZE {} (R* {:.4f} A, eps {:.3f} kcal/mol)".format(lab, rs, ep))
        say("  positive spheres - contact distance to the group's nearest oxygen, over all geometries (own basin only):")
        say("    q      " + "".join("{:>27}".format(g) for g, m, a in GROUPS))
        for qs in QPOS:
            cells = []
            for g, m, a in GROUPS:
                v = np.array([r[7] for r in rows if r[2] == lab and r[3] == qs and r[4] == g and r[5] == "own"])
                nmig = sum(1 for r in rows if r[2] == lab and r[3] == qs and r[4] == g and r[5] == "migrated")
                cells.append("{:2d}/{:2d} {:.3f} [{:.3f}-{:.3f}]".format(len(v), len(v) + nmig, v.mean(), v.min(), v.max()) if len(v) else "  0 own")
            say("    {:+.2f}  ".format(qs) + "".join("{:>27}".format(c) for c in cells))
        say("  negative spheres - outcome per hydrogen start (bound to that H / migrated / escaped / unconverged), and the")
        say("  closest approach to any atom among bound minima:")
        for qs in QNEG:
            sel = [r for r in rows if r[2] == lab and r[3] == qs]
            cnt = {k: sum(1 for r in sel if r[5] == k) for k in ("bound_own", "bound_migrated", "escaped", "unconverged")}
            bound = [r for r in sel if r[5].startswith("bound")]
            near = min(bound, key=lambda r: r[9]) if bound else None
            say("    {:+.2f}: own {} / migrated {} / escaped {} / unconverged {}{}".format(qs, cnt["bound_own"], cnt["bound_migrated"],
                cnt["escaped"], cnt["unconverged"], "" if near is None else "; closest {} at {:.3f} A".format(near[8], near[9])))
    # monotonicity in q (E2.1), N sphere, carboxylates
    mono = tot = 0
    for fr in frames:
        for st in ("R", "TS"):
            for g in ("carbA", "carbB"):
                v = [next((r[7] for r in rows if r[0] == fr and r[1] == st and r[2] == "N" and r[3] == qs and r[4] == g and r[5] == "own"), None) for qs in QPOS]
                if None in v:
                    continue
                tot += 1
                mono += all(v[k] > v[k + 1] for k in range(3))
    say("\nE2.1 positive N sphere, carboxylates: contact distance falls monotonically as q rises in {}/{} frame-state-groups "
        "(expectation: at least 90%) -> {}".format(mono, tot, "as expected" if tot and mono >= 0.9 * tot else "NOT as expected"))
    h8 = [r for r in rows if r[3] < 0 and r[5].startswith("bound") and r[10] < 2.0]
    say("E2.3 H8 flag: {} bound negative-sphere minima closer than 2.0 A to H8{}".format(
        len(h8), "" if not h8 else " (closest {:.3f} A, size {}, q {:+.2f})".format(min(r[10] for r in h8),
                                                                                  min(h8, key=lambda r: r[10])[2], min(h8, key=lambda r: r[10])[3])))
    if fails:
        say("\nSELF-CHECK FAILURES (first 10):")
        for x in fails[:10]:
            say("  " + x)
    (out / "STAGE2_REPORT.txt").write_text("\n".join(rep) + "\n")
    print("\n".join(rep))
    sys.exit(1 if fails else 0)

if __name__ == "__main__":
    main()
