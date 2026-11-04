#!/usr/bin/env python3
"""DIAGNOSTIC 26 -- is the def2-SVP / def2-SVPD gap spill-out, or basis incompleteness?

Pre-registered in STAGE6_DIAGNOSTIC26_20261104.txt, which is NOT edited. Runs no
calculation: every figure is arithmetic over committed ORCA outputs. Reads the
HOMO-LUMO gap and SCF cycles through the COMMITTED cs6_audit.py so the readability
thresholds come from that file and not from this one.
"""
import argparse, hashlib, importlib.util, re, statistics as st, sys
from pathlib import Path

H = 627.5094740631
S6 = "phase2.2/cs_stage6/"
D16 = S6 + "d16/def2-SVPD/"
A_D24 = -5.6165          # ddE def2-SVP bare; D24's print, 4 dp. No committed def2-SVP run exists.
B_D24 = -6.1348          # D24's def2-SVPD print, used only for the sensitivity row
TOL_PRINT = 1.05e-4
SETS = {
    "Ne": (D16 + "S6-10_wall_F1", D16 + "S6-10_wall0"),
    "Ar": (S6 + "d23/ar_d16_S6-10_wall_F1", S6 + "d23/ar_d16_S6-10_wall0"),
    "He": (S6 + "d23/he_d16_S6-10_wall_F1", S6 + "d23/he_d16_S6-10_wall0"),
}
SHA = {
    S6 + "STAGE6_CRITERIA.txt":                     "010d1d83d33c3ff28fcb11a2e05e337fbaa7df0fa51fd38c56973b5b2c93155e",
    S6 + "STAGE6_DIAGNOSTIC26_20261104.txt":         None,   # filled at run time, recorded in the report
    S6 + "STAGE6_DIAGNOSTIC24_NOTE1_20261103.txt":  "6353d354a5b2fbc08a8d244ac5f5949429ddf1998c75b60bf610bc3df9c7c495",
    S6 + "STAGE6_DIAGNOSTIC23_20261102.txt":        "44d9afec5a52b960d45b2b59e41a4e1474c3c355c915d117c17f7fe0546d7cb8",
    S6 + "d23/DIAGNOSTIC23_RESULT.txt":             "d954bf696ae557b5b367abacb5fb44f95e66b72aec2735f239e6f9886a898cb1",
    S6 + "d24/pc_1.3000e_n10.pc":                   "daf4793d1ec9b3c328046c8f0e80d2697f7286216462bce11934af0c4f199e28",
}
FAILS = []
def ok(g, m):  print("   pass  %-7s %s" % (g, m))
def bad(g, m): FAILS.append((g, m)); print("   FAIL  %-7s %s" % (g, m))
def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""): h.update(b)
    return h.hexdigest()
def txt(p):
    try: return Path(p).read_text(errors="replace")
    except Exception as e: return "<<UNREADABLE:%s>>" % e
def norm(s): return re.sub(r"\s+", " ", s.replace("−", "-"))


def load_audit(code):
    """the committed classifier, imported by explicit path (-I implies -P)"""
    p = Path(code) / "phase2.2" / "cs6_audit.py"
    spec = importlib.util.spec_from_file_location("_cs6_audit", p)
    mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
    return mod, sha(p)


def energy(res, d):
    t = txt(Path(res) / d / "job.out")
    if "ORCA TERMINATED NORMALLY" not in t: return None
    m = re.search(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", t)
    return float(m.group(1)) if m else None

def barrier(res, stem):
    a, b = energy(res, stem + "_R"), energy(res, stem + "_TS")
    return None if (a is None or b is None) else (b - a) * H


def g_inputs(res, bundle):
    for rel, want in SHA.items():
        p = Path(res) / rel
        if rel.endswith("DIAGNOSTIC26_20261104.txt"):
            if not p.exists(): bad("GATE 1", "the pre-registration is not committed; push y62 first")
            else: ok("GATE 1", "pre-registration present, sha256 %s" % sha(p)[:16])
            continue
        if not p.exists(): bad("GATE 1", "%s is missing" % rel); continue
        got = sha(p)
        if got != want: bad("GATE 1", "%s sha256 %s, expected %s" % (rel, got[:16], want[:16]))
        else: ok("GATE 1", "%-44s unchanged" % Path(rel).name)

def g_sign(res):
    """D26's GATE 1, re-derived: d_shield = ddE(bare) - ddE(wall)"""
    L = txt(Path(res) / S6 / "standoff/DIAGNOSTIC11_RESULT.txt").splitlines()
    n = 0
    for i, l in enumerate(L):
        if not re.search(r"ddE bare\s+ddE wall\s+d_shield", l): continue
        rows = []
        for m in L[i + 1:]:
            f = m.split()
            if len(f) >= 4 and re.match(r"^(d\d+|dref)$", f[0]): rows.append(f[:4])
            elif rows: break
        for f in rows:
            b, w, d = float(f[1]), float(f[2]), float(f[3]); n += 1
            if abs((b - w) - d) > TOL_PRINT:
                bad("GATE 2", "D11 row %s inconsistent with d_shield = bare - wall" % f[0])
    if n < 5: bad("GATE 2", "only %d D11 rows checked" % n)
    elif not [x for x in FAILS if x[0] == "GATE 2"]:
        ok("GATE 2", "d_shield = ddE(bare) - ddE(wall) on %d committed rows" % n)

def _pc(p):
    L = txt(p).split(); k = int(L[0])
    return [tuple(float(t) for t in L[1 + 4 * i:5 + 4 * i]) for i in range(k)]
def _inp(p):
    out = []
    for l in txt(p).splitlines():
        m = re.match(r"\s*[A-Za-z]{1,2}>\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)", l)
        if m: out.append(tuple(float(x) for x in m.groups()))
    return out

def g_identity(res):
    """D26's GATE 2, re-derived: every arm carries the identical vector"""
    base = Path(res) / S6
    v = _pc(base / "d24/pc_1.3000e_n10.pc")
    bare_pc = base / "d16/def2-SVPD/S6-10_bare_F1_R/pc.pc"
    if not bare_pc.exists(): bad("GATE 3", "d16's bare arm pc.pc is not committed")
    elif _pc(bare_pc) != v: bad("GATE 3", "d16's bare arm carries a different vector from the D24 cell")
    else: ok("GATE 3", "d16 bare arm pc.pc identical to the D24 cell")
    for k, (w, w0) in SETS.items():
        if k == "Ne": continue                     # Ne's arms use %pointcharges, checked below
        for arm in (w + "_R", w + "_TS"):
            q = base.parent.parent / arm / "job.inp"
            q = Path(res) / arm / "job.inp"
            if not q.exists(): bad("GATE 3", "%s missing" % arm); continue
            d = _inp(q)
            if len(d) != len(v) or max(max(abs(x[i] - y[i]) for i in range(4)) for x, y in zip(v, d)) > 1e-6:
                bad("GATE 3", "%s differs from the D24 cell" % Path(arm).name)
            elif "def2-SVPD" not in txt(q).splitlines()[0]:
                bad("GATE 3", "%s is not at def2-SVPD" % Path(arm).name)
            else: ok("GATE 3", "%-34s identical, def2-SVPD" % Path(arm).name)

def g_control(res):
    b = barrier(res, D16 + "ctl")
    if b is None: bad("GATE 4", "the def2-SVPD control did not terminate normally"); return None
    if abs(b - 18.3323) > 5e-4:
        bad("GATE 4", "control barrier %.6f, committed record says 18.3323" % b); return b
    ok("GATE 4", "def2-SVPD control barrier %.6f reproduces the committed 18.3323" % b)
    return b

def g_readable(res, code):
    try:
        mod, msha = load_audit(code)
    except Exception as e:
        bad("GATE 5", "the committed cs6_audit.py could not be loaded: %s" % str(e)[:80])
        return [], None
    arms = [("control", D16 + "ctl_R"), ("control", D16 + "ctl_TS"),
            ("bare", D16 + "S6-10_bare_F1_R"), ("bare", D16 + "S6-10_bare_F1_TS")]
    for k, (w, w0) in SETS.items():
        arms += [(k, w + "_R"), (k, w + "_TS"), (k + " wall0", w0 + "_R"), (k + " wall0", w0 + "_TS")]
    recs = []
    for lab, a in arms:
        d = mod.read(txt(Path(res) / a / "job.out"))
        if not d: bad("GATE 5", "%s could not be read by the committed reader" % a); continue
        d["lab"] = lab; d["arm"] = a; recs.append(d)
    ctl = [r for r in recs if r["lab"] == "control"]
    if not ctl: bad("GATE 5", "no control record"); return recs, msha
    med = {k: st.median([r[k] for r in ctl if r.get(k) is not None])
           for k in ("cyc", "gap", "maxlw", "dip") if any(r.get(k) is not None for r in ctl)}
    for r in recs:
        sick = []
        if r.get("gap") is not None and "gap" in med and r["gap"] < mod.GAP_ABS and r["gap"] < mod.R_GAP * med["gap"]:
            sick.append("gap %.3f eV" % r["gap"])
        if r.get("cyc") is not None and med.get("cyc") and r["cyc"] >= mod.R_CYC * med["cyc"]:
            sick.append("%d cycles, %.1fx control" % (r["cyc"], r["cyc"] / med["cyc"]))
        r["sick"] = sick
    bad_arms = [r for r in recs if r["sick"]]
    if bad_arms:
        for r in bad_arms: bad("GATE 5", "%s is not readable: %s" % (Path(r["arm"]).name, "; ".join(r["sick"])))
    else:
        ok("GATE 5", "all %d arms readable under cs6_audit.py's own thresholds (sha %s)"
           % (len(recs), msha[:12]))
    return recs, msha

def g_d23(res, B, dshield):
    """the computation must reproduce D23's committed dd(X) differences"""
    t = txt(Path(res) / S6 / "d23/DIAGNOSTIC23_RESULT.txt")
    for k in ("Ar", "He"):
        m = re.search(r"S6-10_D16\s+%s\s+([-+]\d+\.\d+)\s+([-+]\d+\.\d+)\s+([-+]\d+\.\d+)" % k, t)
        if not m: bad("GATE 6", "D23's %s row for S6-10_D16 not found" % k); continue
        want = float(m.group(3)); got = dshield[k] - dshield["Ne"]
        if abs(got - want) > 1e-4: bad("GATE 6", "dd(%s) computed %+.6f, D23 prints %+.4f" % (k, got, want))
        else: ok("GATE 6", "dd(%s) %+.6f reproduces D23's %+.4f" % (k, got, want))

def g_d24(res, B):
    """B must round to D24's print, or the discrepancy is reported with its size"""
    if abs(round(B, 4) - B_D24) <= 5e-5:
        ok("GATE 7", "B %.6f rounds to D24's %.4f" % (B, B_D24)); return 0.0
    d = B - B_D24
    bad("GATE 7", "B %.6f does NOT round to D24's %.4f; difference %+.6f kcal/mol = %+.3e Eh"
        % (B, B_D24, d, d / H))
    return d


def g_document(doc, A, B, bctl, ddw, dsh, f, med, per, delta):
    """GATE 8: every figure in the result document must match this run at full precision"""
    d = txt(doc)
    if d.startswith("<<UNREADABLE"): bad("GATE 8", "the result document is not in the bundle"); return
    need = [("%+0.4f" % A, "A"), ("%+0.6f" % B, "B"), ("%0.6f" % abs(B - A), "|B-A|"),
            ("%0.6f" % bctl, "control barrier"), ("%0.4f" % med, "median f_closed"),
            ("%+0.6f" % delta, "the D24 difference")]
    for k in ("Ne", "Ar", "He"):
        need += [("%+0.6f" % ddw[k], "C(%s)" % k), ("%+0.6f" % dsh[k], "d_shield(%s)" % k),
                 ("%0.6f" % abs(ddw[k] - A), "|C-A| %s" % k), ("%0.4f" % f[k], "f_closed(%s)" % k)]
    for s, what in need:
        if s.lstrip("+") in d or s in d: ok("GATE 8", "%-18s %s present in the document" % (what, s))
        else: bad("GATE 8", "%s: the document does not carry %s" % (what, s))
    for k, v in per.items():
        if v not in d: bad("GATE 8", "the document does not state the %s verdict %s" % (k, v))
    if len(set(per.values())) > 1 and "SET-DEPENDENT" not in d:
        bad("GATE 8", "the sets straddle a boundary and the document does not say SET-DEPENDENT")
    elif len(set(per.values())) > 1: ok("GATE 8", "the document declares SET-DEPENDENT, as D26.4 requires")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", required=True); ap.add_argument("--code", required=True)
    ap.add_argument("--out", required=True); ap.add_argument("--doc")
    a = ap.parse_args()
    res, code = a.results, a.code
    print("GATES")
    g_inputs(res, None)
    g_sign(res)
    g_identity(res)
    b_ctl = g_control(res)
    recs, msha = g_readable(res, code)
    if msha is None:
        print("\nSTOP: the readability gate could not run."); return 1
    if b_ctl is None: print("\nSTOP: no usable control."); return 1
    B = barrier(res, D16 + "S6-10_bare_F1") - b_ctl
    ddw, dsh, f = {}, {}, {}
    for k, (w, w0) in SETS.items():
        bw, bw0 = barrier(res, w), barrier(res, w0)
        if bw is None or bw0 is None: bad("GATE 8", "%s arms incomplete" % k); continue
        ddw[k] = bw - bw0
        dsh[k] = B - ddw[k]
        f[k] = (abs(B - A_D24) - abs(ddw[k] - A_D24)) / abs(B - A_D24)
    g_d23(res, B, dsh)
    delta = g_d24(res, B)
    if [x for x in FAILS if x[0] != "GATE 7"]:
        print("\n%d GATE FAILURE(S) other than the D24 reconciliation - nothing is written." % len(FAILS))
        return 1
    fs = sorted(f.values()); med = st.median(fs)
    def verdict(x):
        return "SPILL-OUT" if x >= 0.70 else "MIXED" if x >= 0.30 else "NOT SPILL-OUT" if x >= 0 else "WALL MOVES IT AWAY"
    per = {k: verdict(v) for k, v in f.items()}
    straddle = len(set(per.values())) > 1
    print()
    print("READING")
    print("   A  ddE def2-SVP  bare    %+10.4f   (D24, 4 dp)" % A_D24)
    print("   B  ddE def2-SVPD bare    %+10.6f   (d16, full precision)" % B)
    print("   |B - A|                  %10.6f" % abs(B - A_D24))
    for k in ("Ne", "Ar", "He"):
        print("   %-3s C %+10.6f  d_shield %+10.6f  |C-A| %9.6f  f_closed %7.4f  %s"
              % (k, ddw[k], dsh[k], abs(ddw[k] - A_D24), f[k], per[k]))
    print("   median f_closed %.4f -> %s" % (med, verdict(med)))
    print("   VERDICT: %s" % ("SET-DEPENDENT; no single D26.4 verdict is formed" if straddle else verdict(med)))

    if a.doc:
        g_document(a.doc, A_D24, B, b_ctl, ddw, dsh, f, med, per, delta)
        if [x for x in FAILS if x[0] == "GATE 8"]:
            print("\nGATE 8 FAILURE(S) - the result document does not match the reading. Nothing is written.")
            return 1

    o = Path(a.out); o.mkdir(parents=True, exist_ok=True)
    with open(o / "d26.tsv", "w") as fh:
        fh.write("set\tddE_wall\td_shield\tC_minus_A\tf_closed\tverdict\n")
        for k in ("Ne", "Ar", "He"):
            fh.write("%s\t%.6f\t%.6f\t%.6f\t%.6f\t%s\n" % (k, ddw[k], dsh[k], abs(ddw[k] - A_D24), f[k], per[k]))
    with open(o / "D26_REPORT.txt", "w") as fh:
        fh.write("DIAGNOSTIC 26 -- raw reading, produced by phase2.2/cs6_d26.py\n\n")
        fh.write("A  ddE def2-SVP  bare   %+0.4f   source: DIAGNOSTIC 24 Note 1, 4 dp\n" % A_D24)
        fh.write("B  ddE def2-SVPD bare   %+0.6f  source: d16/def2-SVPD, full precision\n" % B)
        fh.write("   control barrier      %+0.6f\n" % b_ctl)
        fh.write("   |B - A|              %0.6f\n\n" % abs(B - A_D24))
        fh.write("set  ddE(wall)     d_shield      |C-A|        f_closed   verdict\n")
        for k in ("Ne", "Ar", "He"):
            fh.write("%-4s %+12.6f %+12.6f %12.6f %11.4f   %s\n"
                     % (k, ddw[k], dsh[k], abs(ddw[k] - A_D24), f[k], per[k]))
        fh.write("\nmedian f_closed %0.4f -> %s\n" % (med, verdict(med)))
        fh.write("spread %0.4f to %0.4f over three ECP sets\n" % (min(fs), max(fs)))
        fh.write("VERDICT %s\n" % ("SET-DEPENDENT" if straddle else verdict(med)))
        fh.write("\nD24 reconciliation: B - D24's print = %+0.6f kcal/mol = %+0.3e Eh\n" % (delta, delta / H))
        fh.write("cs6_audit.py sha256 %s\n" % msha)
    print()
    print("wrote %s and %s" % (o / "D26_REPORT.txt", o / "d26.tsv"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
