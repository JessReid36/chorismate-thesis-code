#!/usr/bin/env python3
"""
cs6_audit_diff.py - what the two cs6_audit defects cost, measured, not argued.

  python3 cs6_audit_diff.py <results-repo> <code-phase2.2>

Reads every committed ORCA output twice and prints the verdict difference:

  BEFORE  max|Loewdin| over every row of the Loewdin block, including ORCA's inline ECP
          centres (symbol Q, which carry their ASSIGNED charge); HOMO-LUMO from the
          ORBITAL ENERGIES table, which prints LEVEL-SHIFTED energies when TRAH ran.
  AFTER   Q rows excluded; TRAH's own converged gap preferred where it exists.

Both readings carry its own copy of the logic, so the comparison does not depend on which
cs6_audit.py is installed and stays reproducible after the fix lands.

Walled runs are scored with AMENDMENT 1 A1.2's dipole correction applied, because that is what
DIAGNOSTICS 14-17 do. Without it the printed dipole of a walled run carries the ECP centres'
own classical dipole and the 2x rule fires on nearly all of them, which masks the Loewdin
defect and makes the diff report a far smaller number than the diagnostics would see.

EXIT 1 if any verdict moves TOWARD COLLAPSED. That is the only direction that changes a
published readability verdict, and it is the line between a note and an erratum.
"""
import re, sys, tarfile
import numpy as np
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import cs6_audit as AUD
import cs6_transfer as TR
import cs6_negatives as NG

EV = AUD.EV


def rows_of(s):
    k = s.rfind("LOEWDIN ATOMIC CHARGES")
    if k < 0:
        return []
    out = []
    for l in s[k:].splitlines()[2:]:
        m = re.match(r"\s*\d+\s+([A-Za-z]{1,2})\s*:\s*(-?\d+\.\d+)", l)
        if not m:
            break
        out.append((m.group(1), float(m.group(2))))
    return out


def gaps(s):
    k = s.rfind("ORBITAL ENERGIES")
    occ, vir = [], []
    for m in re.finditer(r"^\s+(\d+)\s+([0-9.]+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s*$", s[k:], re.M):
        (occ if float(m.group(2)) > 0.5 else vir).append(float(m.group(3)))
    table = (vir[0] - occ[-1]) * EV if (occ and vir) else None
    tr = re.findall(r"HOMO - LUMO gap\s*:\s*(-?\d+\.\d+)", s)
    return table, (float(tr[-1]) * EV if tr else None)


def geom_and_mem(s):
    """The substrate's 24 atoms and the inline ECP centres with their assigned charges."""
    k = s.find("CARTESIAN COORDINATES (ANGSTROEM)")
    if k < 0:
        return None, None
    g = []
    for l in s[k:].splitlines()[2:]:
        m = re.match(r"\s*([A-Za-z]{1,3}>?)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s*$", l)
        if not m:
            break
        g.append((m.group(1), np.array([float(m.group(i)) for i in (2, 3, 4)])))
    if len(g) < 24:
        return None, None
    ecp = [r for r in g[24:] if r[0].endswith(">")]
    q = [v for sym, v in rows_of(s) if sym == AUD.ECP_SYMBOL] if hasattr(AUD, "ECP_SYMBOL") \
        else [v for sym, v in rows_of(s) if sym == "Q"]
    if not ecp or len(q) != len(ecp):
        return None, None
    M = np.array([NG.MASS[a[0]] for a in g[:24]])
    X = np.array([a[1] for a in g[:24]])
    return (M[:, None] * X).sum(0) / M.sum(), [("Q", q[i], ecp[i][1]) for i in range(len(ecp))]


def verdict(nb, maxlw, dip, gap, cyc, C):
    c = C.get(nb)
    if not c or "gap" not in c:
        return "NO CONTROL"
    moved, sick = [], []
    if maxlw is not None and "maxlw" in c and maxlw - c["maxlw"] >= AUD.D_LW:
        moved.append(1)
    if dip is not None and c.get("dip") and dip / c["dip"] >= AUD.R_DIP:
        moved.append(1)
    if gap is not None and gap < AUD.GAP_ABS and gap < AUD.R_GAP * c["gap"]:
        sick.append(1)
    if cyc is not None and c.get("cyc") and cyc >= AUD.R_CYC * c["cyc"]:
        sick.append(1)
    if moved and sick:
        return "COLLAPSED"
    return "STRAINED" if sick else ("POLARISED" if moved else "CLEAN")


def outputs(root):
    root = Path(root)
    for p in sorted(root.rglob("job.out")):
        if ".git" in p.parts:
            continue
        yield str(p.relative_to(root)), p.read_text(errors="replace"), p.parent
    for a in sorted(x for x in root.rglob("*.tar.xz") if ".git" not in x.parts):
        try:
            t = tarfile.open(str(a))
        except Exception:
            continue
        for m in t.getmembers():
            if not m.isfile() or not m.name.endswith(".out") or m.name.endswith(".pbs.out"):
                continue
            try:
                yield "%s::%s" % (a.relative_to(root), m.name), \
                      t.extractfile(m).read().decode(errors="replace"), None
            except Exception:
                continue


def main(res, cod):
    T6 = Path(res) / "phase2.2" / "cs_stage6"
    C = TR.controls(T6)
    print("control table, from the committed DIAGNOSTIC6_RESULT.txt:")
    for nb in sorted(C):
        c = C[nb]
        print("   nb=%-4d %-12s n=%-4d cyc %-3d gap %.4f maxlw %.4f dip %.2f"
              % (nb, c["basis"], c["n"], c["cyc"], c["gap"], c["maxlw"], c["dip"]))
    print()
    tot = walled = changed = toward = 0
    trans = {}
    examples = []
    for ident, s, d in outputs(res):
        a = AUD.read(s)
        if not a or not a["norm"] or a["notconv"] or a["nb"] not in C:
            continue
        tot += 1
        rows = rows_of(s)
        lw_before = max((abs(v) for sym, v in rows if sym not in AUD.METALS), default=None)
        lw_after = max((abs(v) for sym, v in rows if sym not in AUD.METALS and sym != "Q"), default=None)
        g_table, g_trah = gaps(s)
        dip = a["dip"]
        O, mem = (geom_and_mem(s) if d is not None else (None, None))
        if mem:
            walled += 1
            dip = NG.corrected_dip(d, mem, O)       # A1.2, as D14-D17 apply it
        vb = verdict(a["nb"], lw_before, dip, g_table, a["cyc"], C)
        va = verdict(a["nb"], lw_after, dip, (g_trah if g_trah is not None else g_table), a["cyc"], C)
        if vb != va:
            changed += 1
            trans[(vb, va)] = trans.get((vb, va), 0) + 1
            if va == "COLLAPSED":
                toward += 1
                examples.append("%s  %s -> %s" % (ident, vb, va))
            elif len(examples) < 5:
                examples.append("%s  %s -> %s" % (ident, vb, va))
    print("outputs classified                        : %d" % tot)
    print("   of which walled (A1.2 correction used) : %d" % walled)
    print("verdicts changed by the two fixes         : %d" % changed)
    for (x, y), n in sorted(trans.items()):
        print("   %-10s -> %-10s %d" % (x, y, n))
    print("MOVED TOWARD COLLAPSED                    : %d" % toward)
    for e in examples:
        print("   %s" % e)
    print()
    if toward:
        print("STOP: a published readability verdict would change. This is an erratum, not a note.")
        return 1
    print("No verdict moves toward COLLAPSED. POLARISED and CLEAN are both readable (A1.3),")
    print("so no reading is withheld or granted that was not before, and no Dset changes.")
    return 0


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print(__doc__); raise SystemExit(2)
    raise SystemExit(main(sys.argv[1], sys.argv[2]))
