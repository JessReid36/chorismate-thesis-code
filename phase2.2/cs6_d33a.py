#!/usr/bin/env python3
"""DIAGNOSTIC 33 stage A: read the A matrix at CAM-B3LYP/def2-SVPD, and apply the pre-registered
rule that picks stage B's designs.

The potentials are parsed by the COMMITTED assembler's own function, phase2.2/cs5_amatrix.vpot,
imported rather than reimplemented, so a_i^new is read exactly as the committed a_i was. The sites,
ranges, L and r_min come from the committed optimiser, phase2.2/cs_opt_milp.
"""
import argparse, importlib.util, re, sys, hashlib
from pathlib import Path
import numpy as np

H2KCAL = 627.5094740631
DELTA = 0.05                      # e, cs5_amatrix's finite-difference charge
FD_SITES = (1678, 1954, 1612)
G1_PINNED = {1678: 0.0728, 1954: 0.0611, 1612: 0.0237}   # D33.5 G1, as written
INV_LIMIT = 1e-4                  # Eh/e, the committed invariance limit
CONV = 2.88                       # e, largest Sigma|q| D30 confirmed converged
LADDER = ((3, 0.75, "bare"), (5, 0.50, "bare"), (11, 0.25, "bare"),
          (2, 1.00, "wall"), (3, 1.00, "wall"), (20, 0.25, "bare"))

def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m); return m

def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""): h.update(b)
    return h.hexdigest()

def energy(d):
    t = (Path(d) / "job.out").read_text(errors="replace")
    if "ORCA TERMINATED NORMALLY" not in t: return "did not terminate normally"
    e = re.findall(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", t)
    if len(e) != 1: return "%d FINAL SINGLE POINT ENERGY lines, expected 1" % len(e)
    return float(e[0])

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--d33a", required=True)
    ap.add_argument("--res", required=True)
    ap.add_argument("--code", required=True)
    ap.add_argument("--prereg", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--result")
    a = ap.parse_args()
    W, RES, COD, OUT = Path(a.d33a), Path(a.res), Path(a.code), Path(a.out)
    S6 = RES / "phase2.2" / "cs_stage6"
    OUT.mkdir(parents=True, exist_ok=True)
    gates, rep = [], []

    cs5 = load(COD / "phase2.2" / "cs5_amatrix.py", "cs5")
    opt = load(COD / "phase2.2" / "cs_opt_milp.py", "opt")
    if abs(cs5.DELTA - DELTA) > 1e-12:
        gates.append("the committed finite-difference charge is %r, not %r" % (cs5.DELTA, DELTA))
    if abs(cs5.H - H2KCAL) > 1e-9:
        gates.append("the committed conversion is %r, not %r" % (cs5.H, H2KCAL))

    # ---- every run must be present and clean
    runs = [l.split("\t")[0] for l in (W / "runs.tsv").read_text().splitlines()
            if l.strip() and not l.startswith("#")]
    if len(runs) != 16: gates.append("runs.tsv lists %d runs, expected 16" % len(runs))
    E, V = {}, {}
    for r in runs:
        e = energy(W / r)
        if isinstance(e, str): gates.append("%s: %s" % (r, e)); continue
        E[r] = e
    pts = {}
    for l in (W / "runs.tsv").read_text().splitlines():
        if not l.strip() or l.startswith("#"): continue
        c = l.split("\t")
        if len(c) >= 7 and c[6].strip() and c[6].strip() != "-": pts[c[0]] = c[6].strip()
    for r, p in pts.items():
        v = cs5.vpot(W / r, W / p)
        if isinstance(v, str): gates.append("%s: %s" % (r, v)); continue
        V[r] = v
    rep.append("read %d of %d energies and %d of %d potential sets" % (len(E), len(runs), len(V), len(pts)))
    if gates:
        print("STOP, before anything was computed:")
        for g in gates: print("   " + g)
        return 1

    # ---- the new A matrix
    a_new = (V["sp_42436_TS"] - V["sp_42436_R"]) * H2KCAL
    ap_new = (V["sp_42436_P"] - V["sp_42436_R"]) * H2KCAL
    bar = (E["sp_42436_TS"] - E["sp_42436_R"]) * H2KCAL

    # ---- G2, the invariance check
    dv = float(np.abs(V["inv_42436_R"] - V["sp_42436_R"]).max())
    de = abs(E["inv_42436_R"] - E["sp_42436_R"])
    g2 = dv <= INV_LIMIT
    rep.append("G2  frame 42436 in its own coordinates: max |dV| %.2e Eh/e (limit %.0e), |dE| %.2e Eh -> %s"
               % (dv, INV_LIMIT, de, "PASS" if g2 else "FAIL"))
    if not g2: gates.append("G2 invariance")

    # ---- G1, the finite-difference checks, under ERRATUM 1's resolution
    rep.append("G1  finite difference, |a_fd - a_new| against the stricter of G1's pinned limit and")
    rep.append("    the committed relative formula 0.01 + 0.01|a_new|  (ERRATUM 1)")
    for i in FD_SITES:
        g = lambda s, t: E["fd_42436_%s_%d_%s" % (s, i, t)]
        afd = (g("TS", "p") - g("TS", "m") - g("R", "p") + g("R", "m")) / (2 * DELTA) * H2KCAL
        an = float(a_new[i])
        rel = 0.01 + 0.01 * abs(an)
        lim = min(G1_PINNED[i], rel)
        which = "pinned" if lim == G1_PINNED[i] else "relative"
        ok = abs(afd - an) <= lim
        rep.append("    site %-5d a_new %+9.4f  a_fd %+9.4f  diff %+8.4f   pinned %.4f  relative %.4f"
                   "  bound by %-8s -> %s" % (i, an, afd, afd - an, G1_PINNED[i], rel, which,
                                              "PASS" if ok else "FAIL"))
        if not ok: gates.append("G1 site %d" % i)

    # ---- R33.1, against the committed A
    xyz, layer, ranges = opt.read_sites(str(S6 / "sites_f42436.tsv"))
    a_old = np.asarray(opt.read_matrix(str(S6 / "A_f42436.tsv"), len(xyz), "A")).ravel()
    L = np.asarray(opt.read_matrix(str(S6 / "L_f42436.tsv"), len(xyz), "L")).ravel()
    if len(a_new) != len(a_old):
        gates.append("a_new has %d entries, the committed A has %d" % (len(a_new), len(a_old)))
    if gates:
        print("STOP:")
        for g in gates: print("   " + g)
        for l in rep: print(l)
        return 1
    rk_old = np.argsort(np.argsort(-np.abs(a_old)))
    rk_new = np.argsort(np.argsort(-np.abs(a_new)))
    pear = float(np.corrcoef(a_old, a_new)[0, 1])
    spear = float(np.corrcoef(rk_old.astype(float), rk_new.astype(float))[0, 1])
    top_old = set(np.argsort(-np.abs(a_old))[:20].tolist())
    top_new = set(np.argsort(-np.abs(a_new))[:20].tolist())
    rep.append("")
    rep.append("R33.1  the new A matrix against the committed one, over all %d sites" % len(a_new))
    rep.append("    in vacuo barrier   %.4f kcal/mol at CAM-B3LYP/def2-SVPD, against 18.4610 committed" % bar)
    rep.append("    a range            %+.4f to %+.4f, against %+.4f to %+.4f committed"
               % (a_new.min(), a_new.max(), a_old.min(), a_old.max()))
    rep.append("    correlation        Pearson %.6f, rank %.6f" % (pear, spear))
    rep.append("    largest change     %+.4f kcal/mol per e at site %d"
               % ((a_new - a_old)[np.argmax(np.abs(a_new - a_old))], int(np.argmax(np.abs(a_new - a_old)))))
    rep.append("    median |change|    %.4f;  median |a_new/a_old| %.4f"
               % (float(np.median(np.abs(a_new - a_old))),
                  float(np.median(np.abs(a_new) / np.maximum(np.abs(a_old), 1e-12)))))
    rep.append("    top-20 overlap     %d of 20" % len(top_old & top_new))

    # ---- the pre-registered selection rule
    XYZ = np.asarray(xyz); RMIN = opt.CONFIG["r_min"]; BIG = 1e9
    def best_at(i, q):
        c = [a_new[i] * v + L[i] for lo, hi in ranges[i] for v in (-q, q) if lo <= v <= hi]
        if not c: return None, None
        vals = {a_new[i] * v + L[i]: v for lo, hi in ranges[i] for v in (-q, q) if lo <= v <= hi}
        m = min(vals); return m, vals[m]
    def greedy(vals, k):
        sel = []
        for i in np.argsort(vals):
            if vals[i] >= BIG or len(sel) == k: break
            if all(np.linalg.norm(XYZ[i] - XYZ[j]) >= RMIN for j in sel): sel.append(int(i))
        return sel
    designs = []
    rep.append("")
    rep.append("THE SIX DESIGNS, by D33.2's rule with a_new, r_min = %.2f A" % RMIN)
    for n, q, arm in LADDER:
        vals = np.array([(best_at(i, q)[0] if best_at(i, q)[0] is not None else BIG)
                         for i in range(len(XYZ))])
        sel = greedy(vals, n)
        if len(sel) != n:
            rep.append("    %-2d x %.2f  %-4s  only %d admissible sites found - G6" % (n, q, arm, len(sel)))
            gates.append("G6 %d x %.2f admits only %d sites" % (n, q, len(sel)))
            continue
        qs = [best_at(i, q)[1] for i in sel]
        sq = float(sum(abs(x) for x in qs))
        sc = float(vals[sel].sum())
        band = "confirmed" if sq <= CONV else "outside"
        designs.append(dict(n=n, q=q, arm=arm, sites=sel, charges=qs, sumq=sq, score=sc, band=band))
        rep.append("    %-2d x %.2f  %-4s  Sigma|q| %5.2f e  %-9s  score %+8.3f  sites %s"
                   % (n, q, arm, sq, band, sc, " ".join(str(s) for s in sel[:6]) + (" ..." if n > 6 else "")))
    for l in rep: print(l)
    if gates:
        print()
        print("STOP:")
        for g in gates: print("   " + g)
        return 1

    # ---- outputs
    (OUT / "A_f42436_d33a.tsv").write_text(
        "# a_if = V_TS(r_i) - V_R(r_i), kcal/mol per e, in vacuo CAM-B3LYP D3BJ/def2-SVPD at the\n"
        "# aligned geometries of frame 42436. DIAGNOSTIC 33 stage A, cs6_d33a.py.\n"
        "# row i = site i of sites_f42436.tsv (sha256 %s)\n" % sha(S6 / "sites_f42436.tsv")
        + "\n".join("%.6f" % v for v in a_new) + "\n")
    (OUT / "AP_f42436_d33a.tsv").write_text(
        "# ap_if = V_P(r_i) - V_R(r_i), kcal/mol per e, same conditions. Not used by stage B.\n"
        + "\n".join("%.6f" % v for v in ap_new) + "\n")
    with open(OUT / "designs_d33b.tsv", "w") as fh:
        fh.write("# DIAGNOSTIC 33 stage B designs, by D33.2's rule with a_new. band is D33.3(2).\n")
        fh.write("name\tn\tq\tarm\tsumq\tband\tscore\tsites\tcharges\n")
        for d in designs:
            fh.write("%s\t%d\t%.2f\t%s\t%.4f\t%s\t%.6f\t%s\t%s\n" % (
                "%s_n%02d_q%.2f" % (d["arm"], d["n"], d["q"]), d["n"], d["q"], d["arm"], d["sumq"],
                d["band"], d["score"], ",".join(str(s) for s in d["sites"]),
                ",".join("%+.4f" % c for c in d["charges"])))
    if a.result:
        p = Path(a.result)
        if p.exists(): print("the result document already exists; refusing to overwrite"); return 1
        p.write_text("\n".join([
            "STAGE 6 DIAGNOSTIC 33 STAGE A RESULT",
            "The A matrix at CAM-B3LYP/def2-SVPD, and the designs stage B will run",
            "9 November 2026",
            "",
            "Pre-registration STAGE6_DIAGNOSTIC33_20261109.txt sha256 %s" % sha(a.prereg),
            "Erratum 1 applies to G1's tolerance: STAGE6_DIAGNOSTIC33_ERRATUM1_20261109.txt",
            "Potentials parsed by the committed phase2.2/cs5_amatrix.py vpot, sha256 %s"
            % sha(COD / "phase2.2" / "cs5_amatrix.py"),
            "Sites, ranges, L and r_min from the committed phase2.2/cs_opt_milp.py, sha256 %s"
            % sha(COD / "phase2.2" / "cs_opt_milp.py"),
            "",
        ] + rep) + "\n")
        print()
        print("result document written, %d lines" % len(p.read_text().splitlines()))
    return 0

if __name__ == "__main__":
    sys.exit(main())
