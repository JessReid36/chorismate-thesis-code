#!/usr/bin/env python3
"""DIAGNOSTIC 33 stage B: the barrier effect of the six designs, at CAM-B3LYP/def2-SVPD.

Reliability is the committed classifier phase2.2/cs6_audit.py at its committed thresholds, read
against d31's charge-free control, with the two rules D33.3 fixed in advance: the dipole leg is not
applied to a walled arm, and stability is judged on the LAST analysis per AMENDMENT 5. The gap leg
is treated as unreadable in a level-shifted run, which is the convention DIAGNOSTIC 31 and
DIAGNOSTIC 32 used, so these verdicts are comparable with theirs.
"""
import argparse, importlib.util, re, sys, hashlib
from pathlib import Path

H2KCAL = 627.5094740631
CONV = 2.88
SHIFT = [r"Leaving SCF to start the TRAH-SCF procedure", r"\(TRAH MAcro\)", r"[Ll]evel ?[Ss]hift set to"]
STAB = r"The eigenvalues of the stability matrix:\s*((?:\s*E\(\s*\d+\)\s*=\s*-?\d+\.\d+ Eh\n)+)"
CTL_SHA = {"R": "2f5148d27e932aad2280f0ea7960bdf32fafdf5d769d9efe8911fef0ccadf3f2",
           "TS": "5a650d9ef9d41ff957af9018aa65b9d597dc584148aab78144b976b4a6ba1f80"}

def load(p, n):
    s = importlib.util.spec_from_file_location(n, p); m = importlib.util.module_from_spec(s)
    s.loader.exec_module(m); return m

def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""): h.update(b)
    return h.hexdigest()

def stability(t):
    bl = []
    for mm in re.finditer(STAB, t):
        eig = [float(x) for x in re.findall(r"E\(\s*\d+\)\s*=\s*(-?\d+\.\d+) Eh", mm.group(1))]
        tail = t[mm.end():mm.end() + 400]
        bl.append(dict(v="UNSTABLE" if "wavefunction is unstable" in tail
                       else ("STABLE" if "wavefunction is stable" in tail else None),
                       nneg=sum(1 for e in eig if e < 0), lo=min(eig) if eig else None))
    return dict(ran="WAVEFUNCTION STABILITY ANALYSIS" in t, n=len(bl),
                first=bl[0]["v"] if bl else None, last=bl[-1]["v"] if bl else None,
                nneg=bl[0]["nneg"] if bl else None, lo=bl[0]["lo"] if bl else None)

def run(d, audit):
    p = Path(d) / "job.out"
    if not p.is_file(): return None
    t = p.read_text(errors="replace")
    c = audit.read(t)
    if c is None: return None
    c["shifted"] = any(re.search(x, t) for x in SHIFT)
    c["stab"] = stability(t)
    return c

def main():
    ap = argparse.ArgumentParser()
    for k in ("d33b", "d31", "res", "code", "prereg", "designs", "out"): ap.add_argument("--" + k, required=True)
    ap.add_argument("--result")
    a = ap.parse_args()
    W, D31, RES, COD, OUT = Path(a.d33b), Path(a.d31), Path(a.res), Path(a.code), Path(a.out)
    OUT.mkdir(parents=True, exist_ok=True)
    audit = load(COD / "phase2.2" / "cs6_audit.py", "audit")
    gates, rep = [], []

    for st, h in CTL_SHA.items():
        p = D31 / ("ctl_%s" % st) / "job.out"
        if not p.is_file() or sha(p) != h:
            gates.append("d31/ctl_%s/job.out is missing or not the committed one" % st)
    if gates:
        print("STOP:"); [print("   " + g) for g in gates]; return 1
    ctl = {st: run(D31 / ("ctl_%s" % st), audit) for st in ("R", "TS")}
    if any(v is None for v in ctl.values()):
        print("STOP: the control did not read as a completed single point on this dianion"); return 1

    rows = [[c.strip() for c in l.split("\t")] for l in Path(a.designs).read_text().splitlines()
            if l.strip() and not l.startswith("#") and not l.startswith("name")]
    if len(rows) != 6: gates.append("designs file lists %d designs, expected 6" % len(rows))

    rep.append("PAIRED AGAINST d31's CHARGE-FREE CONTROL at the same level of theory")
    rep.append("   D_LW %.2f e   R_DIP %.1f   GAP_ABS %.1f eV   R_GAP %.1f   (committed thresholds)"
               % (audit.D_LW, audit.R_DIP, audit.GAP_ABS, audit.R_GAP))
    rep.append("   dipole leg withheld for walled arms (D33.3); stability on the LAST analysis (AMENDMENT 5)")
    rep.append("   gap treated as unreadable in a shifted run, as in DIAGNOSTIC 31 and 32")
    rep.append("")

    # ---- the wall-only controls, G7
    w0 = {}
    for d in sorted(W.iterdir()):
        if not d.is_dir() or not d.name.startswith("wall0_"): continue
        n = d.name.split("_")[1]; st = d.name.rsplit("_", 1)[1]
        r = run(d, audit)
        if r is None: gates.append("%s did not read" % d.name); continue
        w0.setdefault(n, {})[st] = r
    rep.append("WALL-ONLY CONTROLS (argon centres, zero charge) - D33.5 G7")
    for n in sorted(w0):
        if set(w0[n]) != {"R", "TS"}: gates.append("wall0 %s is missing a state" % n); continue
        rise = max(w0[n][st]["maxlw"] - ctl[st]["maxlw"] for st in ("R", "TS"))
        moved = rise >= audit.D_LW
        rep.append("   %s   rise %+.4f e   %s   stability %s/%s"
                   % (n, rise, "MOVES" if moved else "does not move",
                      w0[n]["R"]["stab"]["last"], w0[n]["TS"]["stab"]["last"]))
        if moved: gates.append("G7 %s moves density, so its walled designs are not quotable" % n)
    rep.append("")

    # ---- the designs
    out = []
    rep.append("THE SIX DESIGNS")
    rep.append("   design            n   Sigma|q|  band       verdict       ddE         score     ratio")
    for r in rows:
        nm, n, q, arm, sq, band, score = r[0], int(r[1]), float(r[2]), r[3], float(r[4]), r[5], float(r[6])
        chg = [float(x) for x in r[8].split(",")]
        sq_chk = sum(abs(c) for c in chg)
        if abs(sq_chk - sq) > 5e-4: gates.append("%s: Sigma|q| in the file is %.4f, its charges give %.4f" % (nm, sq, sq_chk))
        got = {st: run(W / ("%s_%s" % (nm, st)), audit) for st in ("R", "TS")}
        if any(v is None for v in got.values()):
            rep.append("   %-16s  not read" % nm); gates.append("%s did not read" % nm); continue
        flags, held = [], False
        for st in ("R", "TS"):
            g = got[st]; c = ctl[st]
            if not g["norm"]: flags.append("%s abnormal termination" % st)
            if g["notconv"]: flags.append("%s did not converge" % st)
            if g["nb"] != c["nb"]: gates.append("%s %s: %s basis functions against the control's %s" % (nm, st, g["nb"], c["nb"]))
            rise = g["maxlw"] - c["maxlw"]
            if rise >= audit.D_LW: flags.append("%s relocation: rise %+.4f" % (st, rise))
            if arm != "wall" and c["dip"] and g["dip"] and g["dip"] / c["dip"] >= audit.R_DIP:
                flags.append("%s dipole %.2fx control" % (st, g["dip"] / c["dip"]))
            if g["shifted"] or g["gap"] is None or c["gap"] is None: held = True
            elif g["gap"] < audit.GAP_ABS and g["gap"] < audit.R_GAP * c["gap"]:
                flags.append("%s gap %.4f eV" % (st, g["gap"]))
            s = g["stab"]
            if not s["ran"]: gates.append("%s %s ran no stability analysis" % (nm, st))
            if s["last"] == "UNSTABLE": flags.append("%s unstable at the last analysis, %s negative" % (st, s["nneg"]))
            elif s["first"] == "UNSTABLE":
                rep.append("      note: %s %s first solution unstable (%s negative, lowest %+.6f Eh), corrected"
                           % (nm, st, s["nneg"], s["lo"]))
        verdict = "UNRELIABLE" if flags else ("INCONCLUSIVE" if held else "RELIABLE")
        dde = ((got["TS"]["E"] - ctl["TS"]["E"]) - (got["R"]["E"] - ctl["R"]["E"])) * H2KCAL
        quot = (verdict == "RELIABLE") and (sq <= CONV)
        ratio = dde / score if score else float("nan")
        rep.append("   %-16s %-3d %6.2f e  %-9s  %-12s %+9.4f  %+9.4f  %5.2f%s"
                   % (nm, n, sq, band, verdict, dde, score, ratio, "" if quot else "   NOT QUOTABLE"))
        for f in flags: rep.append("        flag: %s" % f)
        out.append(dict(name=nm, n=n, q=q, arm=arm, sumq=sq, band=band, verdict=verdict,
                        dde=dde, score=score, ratio=ratio, quotable=quot, flags=flags))
    rep.append("")
    quot = [d for d in out if d["quotable"]]
    rep.append("R33.3  QUOTABLE: %s" % (", ".join("%s %+.4f kcal/mol" % (d["name"], d["dde"]) for d in quot)
                                        if quot else "none of the six"))
    if quot:
        best = min(quot, key=lambda d: d["dde"])
        rep.append("       the largest quotable lowering is %s at %+.4f kcal/mol" % (best["name"], best["dde"]))
    rep.append("")
    rep.append("R33.4  THE LINEAR ESTIMATE AGAINST THE MEASUREMENT")
    rep.append("       the record's two pairs, both B3LYP/def2-SVP at the unconstrained optimum:")
    rep.append("         S6-2  score -11.2184  ddE -12.3835  ratio 1.10")
    rep.append("         S6-20 score -19.6537  ddE -24.8261  ratio 1.26")
    rr = [d["ratio"] for d in out if d["verdict"] != "UNRELIABLE"]
    if rr: rep.append("       here, over the %d designs not flagged: %.2f to %.2f" % (len(rr), min(rr), max(rr)))
    for l in rep: print(l)
    if gates:
        print(); print("STOP:"); [print("   " + g) for g in gates]; return 1
    if a.result:
        p = Path(a.result)
        if p.exists(): print("the result already exists; refusing to overwrite"); return 1
        p.write_text("\n".join([
            "STAGE 6 DIAGNOSTIC 33 STAGE B RESULT",
            "The barrier effect of six reliability-validated designs, at CAM-B3LYP/def2-SVPD",
            "",
            "Pre-registration STAGE6_DIAGNOSTIC33_20261109.txt sha256 %s" % sha(a.prereg),
            "Classifier phase2.2/cs6_audit.py sha256 %s" % sha(COD / "phase2.2" / "cs6_audit.py"),
            "Designs designs_d33b.tsv sha256 %s" % sha(a.designs),
            "Control d31/ctl_R and ctl_TS, by the sha256 the pre-registration names",
            "",
        ] + rep) + "\n")
        print(); print("result document written, %d lines" % len(p.read_text().splitlines()))
    return 0

if __name__ == "__main__":
    sys.exit(main())
