#!/usr/bin/env python3
"""Stage 6 DIAGNOSTIC 28's reader.

Reads the 22 committed outputs and applies the decision rules D28.4 fixed before the run.

THE RULE THAT MATTERS MOST: an orbital energy from a run that engaged ORCA's TRAH converger is
reported as UNAVAILABLE, never quoted. In such a run the printed ORBITAL ENERGIES table carries
level-shifted energies, and the negative-gap warnings ORCA emits precede SCF convergence rather
than describing the converged state. DIAGNOSTIC 27 did not have this rule and a claim had to be
retracted because of it. D28.3 GATE 8 is that rule.
"""
import argparse, hashlib, json, re, sys
from pathlib import Path

EV = 27.211386245988
TOL_GAP = 1.0           # eV, D28.4 secondary reading
DROP_MIN = 0.5          # eV, D28.4 PARTIAL threshold
BOUND = 0.0             # eV, D28.4 CONFIRMED threshold
VDE = 0.5123            # eV, Probe 3's measured detachment energy
FUNCS = ["B3LYP", "CAM-B3LYP", "wB97X-D3"]      # in order of long-range exact exchange
LRX = {"B3LYP": 0.20, "CAM-B3LYP": 0.65, "wB97X-D3": 1.00}
VECTORS = ["1.3000e_n10", "2.8800e_n10", "5.6433e_n20"]
NBASIS = {"def2-SVPD": 402, "def2-TZVPD": 682}


class Gates:
    def __init__(self): self.rows, self.failed = [], []
    def __call__(self, ok, gate, msg):
        self.rows.append("   %s  %-7s %s" % ("pass" if ok else "FAIL", gate, msg))
        if not ok: self.failed.append("%s  %s" % (gate, msg))
        return ok
    def report(self): return "\n".join(self.rows)


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""): h.update(b)
    return h.hexdigest()


def read_out(path):
    """Everything D28 needs from one output, with the TRAH question answered explicitly."""
    t = Path(path).read_text(errors="replace")
    d = {"normal": "****ORCA TERMINATED NORMALLY****" in t}
    m = re.search(r"SCF CONVERGED AFTER\s+(\d+)\s+CYCLES", t)
    d["cycles"] = int(m.group(1)) if m else None
    nb = sorted(set(int(x) for x in re.findall(r"Number of basis functions\s*\.*\s*(\d+)", t)))
    d["nbasis"] = nb[0] if len(nb) == 1 else nb
    nel = sorted(set(int(x) for x in re.findall(r"Number of Electrons\s+NEL\s*\.*\s*(\d+)", t)))
    d["nel"] = nel[0] if len(nel) == 1 else nel
    fsp = re.findall(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", t)
    d["energy"] = float(fsp[-1]) if fsp else None
    # --- the functional and dispersion ORCA reports, for GATE 5
    d["exchange"] = (re.findall(r"Exchange Functional\s+Exchange\s*\.*\s*(\S+)", t) or [None])[-1]
    d["correlation"] = (re.findall(r"Correlation Functional\s+Correlation\s*\.*\s*(\S+)", t) or [None])[-1]
    d["disp_blocks"] = len(re.findall(r"DFT DISPERSION CORRECTION", t))
    dv = re.findall(r"Dispersion correction\s+(-?\d+\.\d+)", t)
    d["disp"] = float(dv[-1]) if dv else None
    # --- TRAH: engaged at all, and does any negative-gap warning precede convergence?
    L = t.split("\n")
    w = [i for i, l in enumerate(L) if "negative HOMO - LUMO gap" in l]
    c = [i for i, l in enumerate(L) if "SCF CONVERGED" in l]
    # a bare search for "TRAH" matches the author credits ("Benjamin Helmich-Paris : ... TRAH-")
    # and the settings listing ("Trust-Rad. Augm. Hess. CNVTRAH .... auto"), both present in every
    # ORCA output. These are the markers of TRAH actually ENGAGING:
    d["trah"] = ("Leaving SCF to start the TRAH-SCF procedure" in t) or ("(TRAH MAcro)" in t)
    d["neg_warn"] = len(w)
    d["warn_before_conv"] = bool(w and c and max(w) < max(c))
    d["level_shifted"] = bool(re.search(r"[Ll]evel ?[Ss]hift set to", t))
    # the level shift is what corrupts the printed table, and TRAH is what applies it
    d["shifted"] = d["level_shifted"] or d["trah"]
    # --- orbital energies, located by OCCUPATION
    k = t.rfind("ORBITAL ENERGIES")
    occ, vir = [], []
    if k >= 0:
        for l in t[k:].splitlines()[4:]:
            f = l.split()
            if len(f) != 4:
                if occ and vir: break
                continue
            try: o, eh, ev = float(f[1]), float(f[2]), float(f[3])
            except ValueError: break
            (occ if o > 0.5 else vir).append((eh, ev))
    nocc_expected = (d["nel"] // 2) if isinstance(d["nel"], int) else None
    d["nocc"] = len(occ)
    d["nocc_ok"] = (nocc_expected is not None and len(occ) == nocc_expected)
    if occ and vir and d["nocc_ok"]:
        d["homo_raw"], d["lumo_raw"] = occ[-1][1], vir[0][1]
        d["gap_raw"] = vir[0][1] - occ[-1][1]
    else:
        d["homo_raw"] = d["lumo_raw"] = d["gap_raw"] = None
    # --- and THEN the rule: a level-shifted run yields no quotable orbital energy
    usable = (d["homo_raw"] is not None) and (not d["shifted"])
    d["orbitals_usable"] = usable
    d["homo"] = d["homo_raw"] if usable else None
    d["lumo"] = d["lumo_raw"] if usable else None
    d["gap"] = d["gap_raw"] if usable else None
    return d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tree", required=True); ap.add_argument("--out", required=True)
    ap.add_argument("--d27", required=True)
    a = ap.parse_args()
    tree, out, d27 = Path(a.tree), Path(a.out), Path(a.d27)
    out.mkdir(parents=True, exist_ok=True)
    g = Gates()
    NEW = ["wB97X-D3", "CAM-B3LYP"]

    # ---------------- read everything ---------------------------------------
    R, probes = {}, {}
    for fn in NEW:
        p = tree / "probe" / fn / "job.out"
        probes[fn] = read_out(p) if p.exists() else None
        for basis in ("def2-SVPD", "def2-TZVPD"):
            for st in ("R", "TS"):
                q = tree / fn / basis / ("ctl_%s" % st) / "job.out"
                if q.exists(): R[(fn, basis, "ctl", st)] = read_out(q)
        for v in VECTORS:
            for st in ("R", "TS"):
                q = tree / fn / "def2-SVPD" / ("%s_%s" % (v, st)) / "job.out"
                if q.exists(): R[(fn, "def2-SVPD", v, st)] = read_out(q)
    # B3LYP comes from the committed d27 tree, re-read rather than quoted
    B = {}
    for basis in ("def2-SVPD", "def2-TZVPD"):
        for st in ("R", "TS"):
            q = d27 / basis / ("ctl_%s" % st) / "job.out"
            if q.exists(): B[(basis, "ctl", st)] = read_out(q)
    for v in VECTORS:
        for st in ("R", "TS"):
            q = d27 / "def2-SVPD" / ("%s_%s" % (v, st)) / "job.out"
            if q.exists(): B[("def2-SVPD", v, st)] = read_out(q)

    # ---------------- GATE 1  the probes ------------------------------------
    for fn in NEW:
        d = probes[fn]
        g(d is not None and d["normal"], "GATE 1",
          "%-10s keyword probe terminated normally%s" % (fn, "" if d and d["normal"] else
          " - ORCA does not accept this keyword, so %s is not read" % fn))

    # ---------------- GATE 6  termination and convergence -------------------
    g(len(R) == 20, "GATE 6", "all 20 substrate runs produced an output (%d read)" % len(R))
    nn = [k for k, v in R.items() if not v["normal"]]
    nc = [k for k, v in R.items() if v["cycles"] is None]
    g(not nn, "GATE 6", "every run reports ORCA TERMINATED NORMALLY%s"
      % ("" if not nn else "; did not: " + ", ".join("/".join(map(str, k)) for k in nn[:3])))
    g(not nc, "GATE 6", "every run reports SCF convergence%s"
      % ("" if not nc else "; did not: " + ", ".join("/".join(map(str, k)) for k in nc[:3])))

    # ---------------- GATE 5  one dispersion, and not B3LYP's functional ----
    dd = [k for k, v in R.items() if v["disp_blocks"] != 1]
    g(not dd, "GATE 5", "every run reports exactly one dispersion correction%s"
      % ("" if not dd else "; wrong count: " + ", ".join("/".join(map(str, k)) for k in dd[:3])))
    b3 = [k for k, v in R.items() if v["exchange"] == "B88" and v["correlation"] == "LYP"]
    g(not b3, "GATE 5", "no run reports B88/LYP, which is what ORCA prints for B3LYP%s"
      % ("" if not b3 else "; does: " + ", ".join("/".join(map(str, k)) for k in b3[:3])))
    for fn in NEW:
        ex = sorted(set(v["exchange"] for k, v in R.items() if k[0] == fn))
        co = sorted(set(v["correlation"] for k, v in R.items() if k[0] == fn))
        g(len(ex) == 1 and len(co) == 1, "GATE 5",
          "%-10s reports one functional throughout: exchange %s, correlation %s" % (fn, ex, co))

    # ---------------- GATE 7  basis and electron count ---------------------
    bn = [k for k, v in R.items() if v["nbasis"] != NBASIS.get(k[1])]
    g(not bn, "GATE 7", "every run reports its basis's own function count, 402 / 682%s"
      % ("" if not bn else "; wrong: " + ", ".join("/".join(map(str, k)) for k in bn[:3])))
    en = [k for k, v in R.items() if v["nel"] != 118]
    g(not en, "GATE 7", "every run reports 118 electrons%s"
      % ("" if not en else "; wrong: " + ", ".join("/".join(map(str, k)) for k in en[:3])))
    no = [k for k, v in R.items() if not v["nocc_ok"]]
    g(not no, "GATE 7", "every run's occupied-orbital count is NEL/2 = 59%s"
      % ("" if not no else "; wrong: " + ", ".join("/".join(map(str, k)) for k in no[:3])))

    # ---------------- GATE 8  TRAH is declared, and its orbitals withheld ---
    tr = sorted(k for k, v in R.items() if v["shifted"])
    g(True, "GATE 8", "runs whose orbital energies are level-shifted and therefore NOT quoted: %d%s"
      % (len(tr), "" if not tr else " - " + ", ".join("/".join(map(str, k)) for k in tr[:4])))
    wb = sorted(k for k, v in R.items() if v["warn_before_conv"])
    g(True, "GATE 8", "runs emitting a negative-gap warning BEFORE convergence: %d%s"
      % (len(wb), "" if not wb else " - " + ", ".join("/".join(map(str, k)) for k in wb[:4])))
    leak = [k for k, v in R.items() if v["shifted"] and v["homo"] is not None]
    g(not leak, "GATE 8", "no level-shifted run contributes a quoted orbital energy")

    # ---------------- the readings, by D28.4's rules ------------------------
    H = {}
    for basis in ("def2-SVPD", "def2-TZVPD"):
        for st in ("R", "TS"):
            H[("B3LYP", basis, st)] = B.get((basis, "ctl", st), {}).get("homo")
            for fn in NEW:
                H[(fn, basis, st)] = R.get((fn, basis, "ctl", st), {}).get("homo")
    key = ("wB97X-D3", "def2-SVPD", "R")
    h_new, h_b3 = H.get(key), H.get(("B3LYP", "def2-SVPD", "R"))
    if h_new is None or h_b3 is None:
        primary = "UNAVAILABLE"
        drop = None
    else:
        drop = h_b3 - h_new
        if h_new < BOUND: primary = "CONFIRMED"
        elif drop > DROP_MIN: primary = "PARTIAL"
        else: primary = "NOT CONFIRMED"
    # mechanism: monotone in long-range exact exchange, both states and bases
    mech, mech_rows = True, []
    for basis in ("def2-SVPD", "def2-TZVPD"):
        for st in ("R", "TS"):
            vals = [H.get((f, basis, st)) for f in FUNCS]
            if any(v is None for v in vals): mech = None; mech_rows.append((basis, st, vals, None)); continue
            okm = vals[0] > vals[1] > vals[2]
            mech_rows.append((basis, st, vals, okm))
            if not okm: mech = False
    mech_verdict = ("MONOTONE" if mech is True else ("NOT MONOTONE" if mech is False else "INCOMPLETE"))
    # secondary: the gap at 2.8800 e
    sec = {}
    for fn in NEW + ["B3LYP"]:
        src = R.get((fn, "def2-SVPD", "2.8800e_n10", "R")) if fn in NEW else B.get(("def2-SVPD", "2.8800e_n10", "R"))
        sec[fn] = src.get("gap") if src else None
    gnew = sec.get("wB97X-D3")
    secondary = ("UNAVAILABLE" if gnew is None else
                 ("NARROWING IS REAL" if gnew < TOL_GAP else "NARROWING WAS EXAGGERATED"))

    payload = {"primary": primary, "drop": drop, "H": dict(("|".join(k), v) for k, v in H.items()),
               "mechanism": mech_verdict,
               "mech_rows": [{"basis": b, "state": s, "values": v, "monotone": o} for b, s, v, o in mech_rows],
               "secondary": secondary, "gap_2880": sec, "VDE": VDE, "LRX": LRX,
               "shifted": ["|".join(map(str, k)) for k in tr],
               "warn_before_conv": ["|".join(map(str, k)) for k in wb],
               "cycles": dict(("|".join(map(str, k)), v["cycles"]) for k, v in R.items()),
               "maxlw_dip": dict(("|".join(map(str, k)), None) for k in R),
               "gates_failed": g.failed}
    (out / "d28_values.json").write_text(json.dumps(payload, indent=1, sort_keys=True))
    with open(out / "d28.tsv", "w") as f:
        f.write("functional\tlong_range_exact_exchange\tbasis\tstate\tHOMO_eV\tquotable\n")
        for basis in ("def2-SVPD", "def2-TZVPD"):
            for st in ("R", "TS"):
                for fn in FUNCS:
                    v = H.get((fn, basis, st))
                    f.write("%s\t%.2f\t%s\t%s\t%s\t%s\n" % (fn, LRX[fn], basis, st,
                            ("%.4f" % v) if v is not None else "NA", "yes" if v is not None else "no"))
    with open(out / "D28_REPORT.txt", "w") as f:
        f.write("DIAGNOSTIC 28: the functional test, read from the 22 committed outputs\n\n")
        f.write("GATES\n%s\n\n" % g.report())
        f.write("PRIMARY READING - highest occupied orbital of the BARE substrate, def2-SVPD, R\n")
        for fn in FUNCS:
            v = H.get((fn, "def2-SVPD", "R"))
            f.write("   %-10s long-range exact exchange %.2f   HOMO %s\n"
                    % (fn, LRX[fn], ("%+8.4f eV" % v) if v is not None else "NOT QUOTABLE (level-shifted)"))
        f.write("   Probe 3's detachment energy implies it should sit near %+.4f eV\n" % -VDE)
        if drop is not None: f.write("   wB97X-D3 moves it by %+.4f eV\n" % -drop)
        f.write("   VERDICT: %s\n\n" % primary)
        f.write("MECHANISM CHECK - monotone in long-range exact exchange?\n")
        for b, s, v, o in mech_rows:
            f.write("   %-11s %-3s  %s   %s\n" % (b, s,
                    "  ".join(("%+8.4f" % x) if x is not None else "      NA" for x in v),
                    {True: "monotone", False: "NOT monotone", None: "incomplete"}[o]))
        f.write("   VERDICT: %s\n\n" % mech_verdict)
        f.write("SECONDARY READING - gap at sum|q| = 2.8800 e, def2-SVPD, R\n")
        for fn in FUNCS:
            v = sec.get(fn)
            f.write("   %-10s %s\n" % (fn, ("%7.4f eV" % v) if v is not None else "NOT QUOTABLE"))
        f.write("   TOL %.1f eV   VERDICT: %s\n" % (TOL_GAP, secondary))
        if g.failed:
            f.write("\nFAILED GATES, committed as failures:\n")
            for x in g.failed: f.write("   %s\n" % x)
    print(g.report()); print()
    for fn in FUNCS:
        v = H.get((fn, "def2-SVPD", "R"))
        print("   %-10s LRX %.2f  HOMO %s" % (fn, LRX[fn], ("%+8.4f eV" % v) if v is not None else "NOT QUOTABLE"))
    print("   PRIMARY: %s     MECHANISM: %s     SECONDARY: %s" % (primary, mech_verdict, secondary))
    return 2 if g.failed else 0


if __name__ == "__main__":
    sys.exit(main())
