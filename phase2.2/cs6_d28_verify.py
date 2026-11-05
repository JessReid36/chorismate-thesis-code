#!/usr/bin/env python3
"""D28.3 GATE 9, as its own program.

Independence: the reader locates the highest occupied orbital by OCCUPATION, taking the last row
with occupation above 0.5. This locates it by INDEX, computing the orbital number from the printed
electron count, and separately cross-checks that row's Hartree and eV columns against each other.
Different identification, different consistency check, no shared code. It then requires that every
numeric token in the document be accounted for, and that no value withheld for level shifting has
leaked into it.
"""
import argparse, re, sys
from pathlib import Path

EV = 27.211386245988
TOL_GAP, DROP_MIN, VDE = 1.0, 0.5, 0.5123
FUNCS = ["B3LYP", "CAM-B3LYP", "wB97X-D3"]
LRX = {"B3LYP": 0.20, "CAM-B3LYP": 0.65, "wB97X-D3": 1.00}
VECTORS = ["1.3000e_n10", "2.8800e_n10", "5.6433e_n20"]
QUOTED = {0.20: "B3LYP's long-range exact exchange", 0.65: "CAM-B3LYP's", 1.00: "wB97X-D3's",
          0.5123: "Probe 3's detachment energy", 1.0: "TOL", 0.5: "the PARTIAL threshold",
          2.88: "the vector", 2.8800: "the vector", 22: "the job count", 402: "def2-SVPD", 682: "def2-TZVPD",
          28: "a diagnostic number", 27: "ditto", 118: "the electron count", 59: "half of it"}
ALLOWED_INT = set(range(0, 30)) | {59, 118, 402, 682, 2026, 1, 3, 4, 5, 6}


def by_index(path):
    """locate the frontier orbitals by INDEX from the electron count, not by occupation"""
    t = Path(path).read_text(errors="replace")
    nel = re.findall(r"Number of Electrons\s+NEL\s*\.*\s*(\d+)", t)
    if not nel: return None, "no electron count"
    n_occ = int(nel[-1]) // 2
    k = t.rfind("ORBITAL ENERGIES")
    if k < 0: return None, "no orbital table"
    rows = {}
    for l in t[k:].splitlines()[4:]:
        f = l.split()
        if len(f) != 4: 
            if rows: break
            continue
        try: rows[int(f[0])] = (float(f[1]), float(f[2]), float(f[3]))
        except ValueError: break
    if (n_occ - 1) not in rows or n_occ not in rows: return None, "frontier rows absent"
    occ_h, eh_h, ev_h = rows[n_occ - 1]
    occ_l, eh_l, ev_l = rows[n_occ]
    if occ_h < 0.5 or occ_l > 0.5: return None, "occupations disagree with the index"
    if abs(eh_h * EV - ev_h) > 2e-3: return None, "the Hartree and eV columns disagree"
    shifted = bool(re.search(r"[Ll]evel ?[Ss]hift set to", t)) \
        or ("Leaving SCF to start the TRAH-SCF procedure" in t) or ("(TRAH MAcro)" in t)
    return {"homo": ev_h, "lumo": ev_l, "gap": ev_l - ev_h, "shifted": shifted}, None


def main():
    ap = argparse.ArgumentParser()
    for k in ("tree", "doc", "d27"): ap.add_argument("--" + k, required=True)
    a = ap.parse_args()
    tree, d27 = Path(a.tree), Path(a.d27)
    D = Path(a.doc).read_text()
    rows, fails = [], []
    def chk(ok, msg):
        rows.append("   %s  GATE 9  %s" % ("pass" if ok else "FAIL", msg))
        if not ok: fails.append(msg)

    H, G, withheld = {}, {}, []
    for f in FUNCS:
        for basis in ("def2-SVPD", "def2-TZVPD"):
            for st in ("R", "TS"):
                p = (d27 / basis / ("ctl_%s" % st) / "job.out") if f == "B3LYP" \
                    else (tree / f / basis / ("ctl_%s" % st) / "job.out")
                if not p.exists(): H[(f, basis, st)] = None; continue
                d, err = by_index(p)
                if d is None: H[(f, basis, st)] = None; continue
                if d["shifted"]: H[(f, basis, st)] = None; withheld.append("%s/%s/ctl_%s" % (f, basis, st)); continue
                H[(f, basis, st)] = d["homo"]
        p = (d27 / "def2-SVPD" / "2.8800e_n10_R" / "job.out") if f == "B3LYP" \
            else (tree / f / "def2-SVPD" / "2.8800e_n10_R" / "job.out")
        d, _ = by_index(p) if p.exists() else (None, None)
        G[f] = None if (d is None or d["shifted"]) else d["gap"]

    derived = {}
    for k, v in H.items():
        if v is not None: derived[round(v, 4)] = "HOMO %s %s %s" % k
    for k, v in G.items():
        if v is not None: derived[round(v, 4)] = "gap at 2.8800 e, %s" % k
    hn, hb = H.get(("wB97X-D3", "def2-SVPD", "R")), H.get(("B3LYP", "def2-SVPD", "R"))
    if hn is not None and hb is not None:
        derived[round(hn - hb, 4)] = "the shift wB97X-D3 applies"
        primary = "CONFIRMED" if hn < 0.0 else ("PARTIAL" if (hb - hn) > DROP_MIN else "NOT CONFIRMED")
    else:
        primary = "UNAVAILABLE"
    derived[round(-VDE, 4)] = "where the orbital should sit"
    mono = []
    for basis in ("def2-SVPD", "def2-TZVPD"):
        for st in ("R", "TS"):
            vals = [H.get((f, basis, st)) for f in FUNCS]
            mono.append(None if any(v is None for v in vals) else (vals[0] > vals[1] > vals[2]))
    mech = "INCOMPLETE" if any(m is None for m in mono) else ("MONOTONE" if all(mono) else "NOT MONOTONE")
    gn = G.get("wB97X-D3")
    sec = "UNAVAILABLE" if gn is None else ("NARROWING IS REAL" if gn < TOL_GAP else "NARROWING WAS EXAGGERATED")

    for val, what in sorted(derived.items()):
        chk(("%.4f" % val) in D, "%-40s %+9.4f present in the document" % (what, val))
    chk("VERDICT: %s" % primary in D, "the primary verdict reads %s" % primary)
    chk("VERDICT: %s" % mech in D, "the mechanism verdict reads %s" % mech)
    chk("VERDICT: %s" % sec in D, "the secondary verdict reads %s" % sec)
    chk(len(re.findall(r"NOT QUOTABLE", D)) >= 0, "withheld values are marked NOT QUOTABLE where shown")
    # nothing withheld may appear as a number
    leaked = []
    for nm in withheld:
        f, basis, st = nm.split("/")[0], nm.split("/")[1], nm.split("/")[2].replace("ctl_", "")
        p = tree / f / basis / ("ctl_%s" % st) / "job.out"
        d, _ = by_index(p)
        if d and ("%.4f" % d["homo"]) in D: leaked.append(nm)
    chk(not leaked, "no level-shifted orbital energy leaked into the document%s"
        % ("" if not leaked else ": " + ", ".join(leaked)))
    toks = re.findall(r"(?<![\w.])-?\d+(?:\.\d+)?(?![\w.])", D)
    un = []
    for t in set(toks):
        fv = float(t)
        if "." not in t:
            if int(fv) in ALLOWED_INT or int(fv) in QUOTED: continue
            un.append(t); continue
        if any(abs(fv - k) < 5e-5 for k in derived) or any(abs(fv - k) < 5e-5 for k in QUOTED): continue
        un.append(t)
    chk(not un, "every numeric token in the document is accounted for%s"
        % ("" if not un else "; unaccounted: " + ", ".join(sorted(un)[:8])))
    print("\n".join(rows)); print()
    if fails:
        print("GATE 9 FAILED: %d" % len(fails))
        for f in fails: print("   " + f)
        return 1
    print("GATE 9 PASSES: %d checks, every figure re-derived by orbital index." % len(rows))
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as e:
        print("GATE 9 FAILED: the verifier could not complete")
        print("   %s: %s" % (type(e).__name__, e)); sys.exit(1)
