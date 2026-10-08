#!/usr/bin/env python3
"""
cs4b_analyse.py - Stage 4b: applies STAGE4B_CRITERIA.txt, as amended by STAGE4B_AMENDMENT1.txt, to
phase2.2/cs_stage4b/ (and, for the bare sphere at 2.674 A, Stage 4's committed single points in phase2.2/cs_stage4/).
Python 3.6, standard library only. Writes STAGE4B_REPORT.txt.
Amendment 1 (Arm B only; criteria and thresholds unchanged): ORCA would not attach the %basis pseudopotential to a
zero-charge Ne> centre and refuses a Ne> centre without one, so q = 0 comes from B2_inline_<basis>_q0 (or, if that
did not run, the mean of B2_inline_<basis>_q{p,m}00001), and the check-1 controls are B2_pc_svp_q{p,m}1 (the sphere as
an external point-charge file, no pseudopotential).
USAGE  python3 cs4b_analyse.py [phase2.2 folder]
"""
import os, re, sys
P2 = sys.argv[1] if len(sys.argv) > 1 else "/home/18660916/system_development/phase2.2"
B4, S4 = os.path.join(P2, "cs_stage4b"), os.path.join(P2, "cs_stage4")
H2KCAL = 627.5094740631
QS = ((-1.0, "qm1"), (-0.5, "qm05"), (0.0, "q0"), (0.5, "qp05"), (0.75, "qp075"), (1.0, "qp1"))
DIST = (("d2674", 2.674), ("d290", 2.90), ("d320", 3.20), ("d350", 3.50))

def text(p):
    return open(p, errors="replace").read() if os.path.exists(p) else ""
def parse(p, qmmm):
    t = text(p)
    if "ORCA TERMINATED NORMALLY" not in t:
        return None
    pat = r"FINAL SINGLE POINT ENERGY \(QM/MM\)\s+(-?\d+\.\d+)" if qmmm else r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)"
    e = re.findall(pat, t)
    i = t.rfind("ORBITAL ENERGIES"); orbs = []
    for l in t[i:i + 40000].splitlines()[4:]:
        m = re.match(r"\s*(\d+)\s+(\d\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)", l)
        if m: orbs.append((float(m.group(2)), float(m.group(4))))
        elif orbs: break
    occ = [o[1] for o in orbs if o[0] > 0.5]; vir = [o[1] for o in orbs if o[0] < 0.5]
    nel = re.findall(r"Number of Electrons\s+NEL\s+\.+\s+(\d+)", t)
    nbf = re.findall(r"Number of basis functions\s+\.+\s+(\d+)", t)
    d3 = re.findall(r"Dispersion correction\s+(-?\d+\.\d+)", t)
    lj = re.findall(r"LJ interaction\s+(-?\d+\.\d+)", t)
    if not e:
        return None
    return {"E": float(e[-1]), "homo": occ[-1] if occ else None, "lumo": vir[0] if vir else None,
            "nel": int(nel[-1]) if nel else None, "nbf": int(nbf[0]) if nbf else None,
            "d3": float(d3[-1]) if d3 else None, "lj": float(lj[-1]) if lj else None, "text": t}

def fits(E):
    """E: dict q -> energy (Eh). Returns a, b (symmetric, exact through q^4), and the negative-side quadratic
    reference (from q = -1, -0.5, 0) deviations at +0.5, +0.75, +1 in kcal/mol."""
    D1 = (E[1.0] - E[-1.0]) * H2KCAL; D5 = (E[0.5] - E[-0.5]) * H2KCAL
    S1 = (E[1.0] + E[-1.0] - 2 * E[0.0]) * H2KCAL; S5 = (E[0.5] + E[-0.5] - 2 * E[0.0]) * H2KCAL
    a = (8 * D5 - D1) / 6; b = (16 * S5 - S1) / 6
    x1 = (E[-1.0] - E[0.0]) * H2KCAL; x5 = (E[-0.5] - E[0.0]) * H2KCAL
    beta = 2 * (x1 - 2 * x5); alpha = beta - x1
    dev = {q: (E[q] - E[0.0]) * H2KCAL - (alpha * q + beta * q * q) for q in (0.5, 0.75, 1.0) if q in E}
    return a, b, dev

rep = []; say = rep.append
say("Stage 4b report (cs4b_analyse.py) - criteria from STAGE4B_CRITERIA.txt")
# ---------------- Arm A
say("\nARM A - bare sphere (QM/MM), def2-SVP vs def2-SVPD, sphere moved out along O2 -> sphere")
A = {}
for bt in ("svp", "svpd"):
    for dt, d in DIST:
        R = {}
        for q, qt in QS:
            if dt == "d2674" and qt != "qp075":
                R[q] = parse(os.path.join(S4, "sp_%s_%s" % (bt, qt.replace("qp075", "")), "job.out"), True)
            else:
                R[q] = parse(os.path.join(B4, "A_%s_%s_%s" % (bt, dt, qt), "job.out"), True)
        A[bt, dt] = R
safe = {}
for dt, d in DIST:
    ok_all = True
    res = {}
    for bt in ("svp", "svpd"):
        R = A[bt, dt]
        if any(R[q] is None for q, _ in QS):
            ok_all = False; continue
        E = {q: R[q]["E"] for q, _ in QS}
        res[bt] = fits(E) + ({q: R[q]["lumo"] - R[q]["homo"] for q, _ in QS},)
    if not ok_all or len(res) < 2:
        say("  %.3f A: a run is missing or did not terminate normally" % d); continue
    (a0, b0, dev0, gap0), (a1, b1, dev1, gap1) = res["svp"], res["svpd"]
    brel = abs(b1 - b0) / abs(b0)
    say("  %.3f A from O2:  a SVP %.3f / SVPD %.3f;  b SVP %.3f / SVPD %.3f (%.1f%%) -> b criterion %s" % (
        d, a0, a1, b0, b1, 100 * brel, "met" if brel <= 0.10 else "NOT met"))
    for q in (0.5, 0.75, 1.0):
        s = gap1[q] >= 2.0 and abs(dev1[q]) <= 0.5
        safe[dt, q] = s and brel <= 0.10
        say("     q=%+.2f  SVPD gap %.2f eV, deviation from trend %+.3f kcal/mol (SVP %.2f eV, %+.3f) -> %s" % (
            q, gap1[q], dev1[q], gap0[q], dev0[q], "safe" if safe[dt, q] else "NOT safe"))
say("  safe distance by charge (smallest tested distance from O2 at which that charge, and every larger tested distance, is safe):")
for q in (0.5, 0.75, 1.0):
    ds = [d for dt, d in DIST if all(safe.get((dt2, q), False) for dt2, d2 in DIST if d2 >= d)]
    say("     q=%+.2f: %s" % (q, "%.3f A" % min(ds) if ds else "none of the tested distances"))
# ---------------- Arm B
say("\nARM B - sphere = coreless ECP centre Ne> with charge q and the Ne-type pseudopotential (Marefat Khah 2020, SI S2)")
B = {(bt, q): parse(os.path.join(B4, "B_ecp_%s_%s" % (bt, qt), "job.out"), False) for bt in ("svp", "svpd") for q, qt in QS if qt != "q0"}
extra = []                                     # replacement runs used, for check 0
for bt in ("svp", "svpd"):                     # Amendment 1: q = 0 with the pseudopotential given inline
    z = parse(os.path.join(B4, "B2_inline_%s_q0" % bt, "job.out"), False)
    zp = parse(os.path.join(B4, "B2_inline_%s_qp00001" % bt, "job.out"), False)
    zm = parse(os.path.join(B4, "B2_inline_%s_qm00001" % bt, "job.out"), False)
    if z is not None:
        B[bt, 0.0] = z; extra.append(z); src = "B2_inline_%s_q0" % bt
        if zp is not None and zm is not None:
            extra += [zp, zm]
            src += "; the +-0.0001 pair's mean differs by %+.2e Eh" % ((zp["E"] + zm["E"]) / 2 - z["E"])
    elif zp is not None and zm is not None:
        z = dict(zp); z["E"] = (zp["E"] + zm["E"]) / 2; B[bt, 0.0] = z; extra += [zp, zm]
        src = "mean of B2_inline_%s_q+-0.0001 (the q = 0 run did not complete)" % bt
    else:
        B[bt, 0.0] = None; src = "missing"
    say("  q = 0 (%s): %s" % (bt, src))
c1 = parse(os.path.join(B4, "B2_pc_svp_qp1", "job.out"), False)
c2 = parse(os.path.join(B4, "B2_pc_svp_qm1", "job.out"), False)
ci = parse(os.path.join(B4, "B_inline_svp_qp1", "job.out"), False)
runs = [v for v in B.values()] + [c1, c2, ci] + extra
if any(v is None for v in runs):
    say("  a run is missing or did not terminate normally -> Arm B not assessable")
else:
    nel = {v["nel"] for v in runs}; nbf = {(k[0], v["nbf"]) for k, v in B.items()}
    ck0 = nel == {118} and nbf == {("svp", 264), ("svpd", 402)} and c1["nbf"] == c2["nbf"] == ci["nbf"] == 264 \
        and {v["nbf"] for v in extra} <= {264, 402}
    say("  check 0 (118 electrons in every run; the centre adds no basis functions: 264 / 402): %s  [NEL %s; basis %s]" % (
        "PASS" if ck0 else "FAIL", sorted(nel), sorted(nbf)))
    ok1 = True
    for c, qt in ((c1, "qp1"), (c2, "qm1")):
        s4 = parse(os.path.join(S4, "sp_svp_%s" % qt, "job.out"), False)   # the QM part of Stage 4's QM/MM energy
        dd3 = (c["d3"] or 0.0) - (s4["d3"] or 0.0)
        diff = c["E"] - s4["E"]
        ok1 = ok1 and abs(diff - dd3) <= 1e-6
        say("  check 1 control (point-charge file, no ECP, %s) vs the QM part of Stage 4's QM/MM energy: %+.2e Eh; dispersion-correction "
            "difference %+.2e Eh; residual %+.2e Eh" % (qt, diff, dd3, diff - dd3))
    say("  check 1 -> %s" % ("PASS" if ok1 else "FAIL"))
    dinl = ci["E"] - B["svp", 1.0]["E"]
    say("  check 2 inline ECP vs %%basis ECP (q=+1, SVP): %+.2e Eh -> %s" % (dinl, "PASS" if abs(dinl) <= 1e-7 else "FAIL"))
    applied = abs(B["svp", 1.0]["E"] - c1["E"]) > 1e-4
    say("  check 3 the ECP acts (ECP vs no-ECP at q=+1 differ by %.4f kcal/mol) -> %s" % (
        (B["svp", 1.0]["E"] - c1["E"]) * H2KCAL, "PASS" if applied else "FAIL (ECP not applied)"))
    say("  ECP lines in the inline run's output (first 12 matching, for inspection):")
    blk = [l for l in ci["text"].splitlines() if re.search(r"ECP|N_core|lmax|NewECP", l)][:12]
    for l in blk: say("     | " + l.rstrip()[:110])
    E0 = {q: B["svp", q]["E"] for q, _ in QS}; E1 = {q: B["svpd", q]["E"] for q, _ in QS}
    a0, b0, dev0 = fits(E0); a1, b1, dev1 = fits(E1)
    gap1 = B["svpd", 1.0]["lumo"] - B["svpd", 1.0]["homo"]
    brel = abs(b1 - b0) / abs(b0)
    B1 = gap1 >= 2.0 and abs(dev1[1.0]) <= 0.5 and brel <= 0.10
    say("  B1 (SVPD, ECP sphere at contact): gap at q=+1 %.2f eV; deviation from trend at q=+1 %+.3f kcal/mol "
        "(+0.75: %+.3f); b SVP %.3f / SVPD %.3f (%.1f%%) -> %s" % (gap1, dev1[1.0], dev1[0.75], b0, b1, 100 * brel, "PASS" if B1 else "FAIL"))
    S4E = {q: parse(os.path.join(S4, "sp_svp_%s" % qt, "job.out"), True)["E"] for q, qt in QS if qt != "qp075"}
    S4E[0.75] = A["svp", "d2674"][0.75]["E"] if A["svp", "d2674"][0.75] else S4E[1.0]
    ab, bb, _ = fits(S4E)
    B2 = abs(a0 - ab) <= 0.01 * abs(ab) and abs(b0 - bb) <= 0.10 * abs(bb)
    say("  B2 (decision input; SVP, ECP sphere vs bare sphere): a %.3f vs %.3f (%+.2f%%); b %.3f vs %.3f (%+.1f%%) -> %s" % (
        a0, ab, 100 * (a0 - ab) / abs(ab), b0, bb, 100 * (b0 - bb) / abs(bb),
        "within 1% / 10%: adopting the ECP leaves the A matrix and Stage 2-3 conclusions unchanged" if B2 else
        "outside 1% / 10%: adopting the ECP changes design energies; Stage 2-3 conclusions to be revisited"))
    s40 = parse(os.path.join(S4, "sp_svp_q0", "job.out"), False)    # QM part of Stage 4's q = 0 energy
    s40m = parse(os.path.join(S4, "sp_svp_q0", "job.out"), True)
    pauli = (B["svp", 0.0]["E"] - s40["E"]) * H2KCAL
    pd3 = ((B["svp", 0.0]["d3"] or 0.0) - (s40["d3"] or 0.0)) * H2KCAL
    say("  informative: energy the pseudopotential alone adds at contact (q=0, SVP) %+.3f kcal/mol, of which %+.3f is a "
        "change in the dispersion correction (if ORCA's D3 counts the centre); the sphere's LJ energy at the same "
        "geometry %+.3f kcal/mol" % (pauli, pd3, s40m["lj"] * H2KCAL))
open(os.path.join(B4, "STAGE4B_REPORT.txt"), "w").write("\n".join(rep) + "\n")
print("\n".join(rep))
