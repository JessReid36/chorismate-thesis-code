#!/usr/bin/env python3
"""DIAGNOSTIC 31 reader: the (n, q) envelope, by the rules pre-registered in
STAGE6_DIAGNOSTIC31_20261106.txt.

Every threshold comes from the committed classifier phase2.2/cs6_audit.py, loaded by explicit
path. Nothing is transcribed here. A condition is judged in BOTH states and is never called
RELIABLE on a gap that could not be read.
"""
import argparse, hashlib, importlib.util, json, re
from pathlib import Path

NS = (1, 2, 3, 5, 10, 20)
LADDER = (0.25, 0.50, 0.75, 1.00)
H2KCAL = 627.5094740631
# the narrow level-shift markers. A bare "TRAH" matches the author credits in every ORCA output.
SHIFT = [r"Leaving SCF to start the TRAH-SCF procedure", r"\(TRAH MAcro\)", r"[Ll]evel ?[Ss]hift set to"]


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""): h.update(b)
    return h.hexdigest()


def load_audit(code):
    p = Path(code) / "phase2.2" / "cs6_audit.py"
    if not p.exists(): raise SystemExit("the committed classifier is not at %s" % p)
    spec = importlib.util.spec_from_file_location("_a", p)
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
    return m, sha(p)


def stability(t):
    blocks = []
    for mm in re.finditer(r"The eigenvalues of the stability matrix:\s*((?:\s*E\(\s*\d+\)\s*=\s*-?\d+\.\d+ Eh\n)+)", t):
        eig = [float(x) for x in re.findall(r"E\(\s*\d+\)\s*=\s*(-?\d+\.\d+) Eh", mm.group(1))]
        tail = t[mm.end():mm.end() + 400]
        v = "UNSTABLE" if "wavefunction is unstable" in tail else ("STABLE" if "wavefunction is stable" in tail else None)
        blocks.append({"eigenvalues": eig[:6], "n_negative": sum(1 for e in eig if e < 0), "verdict": v})
    return {"ran": "WAVEFUNCTION STABILITY ANALYSIS" in t, "n_blocks": len(blocks),
            "first": blocks[0]["verdict"] if blocks else None,
            "last": blocks[-1]["verdict"] if blocks else None,
            "n_negative": blocks[0]["n_negative"] if blocks else None,
            "eigenvalues": blocks[0]["eigenvalues"] if blocks else []}


def read_run(d, audit):
    jo = Path(d) / "job.out"
    if not jo.exists(): return {"present": False}
    t = jo.read_text(errors="replace")
    e = re.findall(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", t)
    c = audit.read(t)
    return {"present": True, "normal": "ORCA TERMINATED NORMALLY" in t,
            "converged": bool(e) and not (re.search(r"MAXIMUM NUMBER OF .*ITERATIONS", t, re.I)
                                          or "SCF NOT CONVERGED AFTER" in t),
            "E": float(e[-1]) if e else None,
            "gap": c.get("gap"), "maxlw": c.get("maxlw"), "dip": c.get("dip"),
            "shifted": any(re.search(p, t) for p in SHIFT),
            "stability": stability(t), "sha": sha(jo)}


def main():
    a = argparse.ArgumentParser()
    a.add_argument("--d31", required=True); a.add_argument("--code", required=True)
    a.add_argument("--prereg", required=True); a.add_argument("--out", required=True)
    a.add_argument("--result")
    a = a.parse_args()
    audit, audit_sha = load_audit(a.code)
    D, OUT = Path(a.d31), Path(a.out)
    gates, notes = [], []
    D_LW, R_DIP, GAP_ABS, R_GAP = audit.D_LW, audit.R_DIP, audit.GAP_ABS, audit.R_GAP

    ctl = {}
    for st in ("R", "TS"):
        r = read_run(D / ("ctl_%s" % st), audit)
        if not r["present"]: gates.append("the ctl_%s control is missing" % st)
        ctl[st] = r
    if gates: return bail(gates, OUT)

    CONDS = ([("A/n%02d" % n, "A", n, 1.00) for n in NS]
             + [("B/n%02d_q%.2f" % (n, q), "B", n, q) for n in NS for q in LADDER])
    rows = {}
    for name, arm, n, q in CONDS:
        states = {}
        for st in ("R", "TS"):
            r = read_run(D / ("%s_%s" % (name, st)), audit)
            if not r["present"]: gates.append("%s_%s is missing" % (name, st)); continue
            c = ctl[st]
            flags = []
            if not r["converged"]: flags.append("did not converge")
            if not r["normal"]: flags.append("abnormal termination")
            rise = (r["maxlw"] - c["maxlw"]) if (r["maxlw"] is not None and c["maxlw"] is not None) else None
            ratio = (r["dip"] / c["dip"]) if (r["dip"] is not None and c["dip"]) else None
            if rise is not None and rise >= D_LW: flags.append("relocation: Loewdin rise %+.4f" % rise)
            if ratio is not None and ratio >= R_DIP: flags.append("relocation: dipole ratio %.2f" % ratio)
            # the gap is withheld where a level shift makes the printed orbital energies the floor
            gap_read = (not r["shifted"]) and r["gap"] is not None and c["gap"] is not None
            if gap_read and r["gap"] < GAP_ABS and r["gap"] < R_GAP * c["gap"]:
                flags.append("gap %.4f eV, %.2f of control" % (r["gap"], r["gap"] / c["gap"]))
            if r["stability"]["first"] == "UNSTABLE":
                flags.append("unstable, %s negative eigenvalue(s)" % r["stability"]["n_negative"])
            if not r["stability"]["ran"]: gates.append("%s_%s ran no stability analysis" % (name, st))
            states[st] = {"flags": flags, "rise": rise, "ratio": ratio, "gap": r["gap"],
                          "gap_read": gap_read, "shifted": r["shifted"], "E": r["E"],
                          "stability": r["stability"], "converged": r["converged"]}
        if len(states) != 2: continue
        allflags = sorted(set(f for s in states.values() for f in s["flags"]))
        anyheld = any(not s["gap_read"] for s in states.values())
        verdict = "UNRELIABLE" if allflags else ("INCONCLUSIVE" if anyheld else "RELIABLE")
        bar = None
        if all(s["E"] is not None for s in states.values()):
            bar = (states["TS"]["E"] - states["R"]["E"]) * H2KCAL
        rows[name] = {"arm": arm, "n": n, "q": q, "sumabsq": n * q, "maxabsq": q,
                      "sumq2": n * q * q, "verdict": verdict, "flags": allflags,
                      "gap_withheld": anyheld, "barrier": bar, "states": states}

    # ---- D31.4  the two envelopes -----------------------------------------
    armA = {}
    seq = [(n, rows["A/n%02d" % n]["verdict"]) for n in NS if ("A/n%02d" % n) in rows]
    env, broke, mono = None, False, True
    for n, v in seq:
        if v == "RELIABLE" and not broke: env = n
        elif v != "RELIABLE": broke = True
        elif broke: mono = False
    armA = {"per_n": dict(seq), "envelope": env if mono else None, "monotone": mono}
    armB = {}
    for n in NS:
        s = [(q, rows["B/n%02d_q%.2f" % (n, q)]["verdict"]) for q in LADDER
             if ("B/n%02d_q%.2f" % (n, q)) in rows]
        e, br, mo = None, False, True
        for q, v in s:
            if v == "RELIABLE" and not br: e = q
            elif v != "RELIABLE": br = True
            elif br: mo = False
        armB[n] = {"per_q": {("%.2f" % q): v for q, v in s}, "envelope": e if mo else None, "monotone": mo}

    # ---- D31.5  the iso-total test ----------------------------------------
    fam = {}
    for name, r in rows.items(): fam.setdefault(round(r["sumabsq"], 4), []).append(name)
    iso = {}
    for s, mem in sorted(fam.items()):
        if len(mem) < 2: continue
        vs = set(rows[m]["verdict"] for m in mem)
        iso["%.2f" % s] = {"members": sorted(mem), "verdicts": {m: rows[m]["verdict"] for m in sorted(mem)},
                           "agree": len(vs) == 1}
    n_dis = sum(1 for v in iso.values() if not v["agree"])
    iso_verdict = ("ALL %d FAMILIES AGREE" % len(iso) if n_dis == 0 else
                   "%d of %d FAMILIES DISAGREE" % (n_dis, len(iso)))

    # ---- D31.6  the abandonment conditions, evaluated by their own rules --
    B = [r for r in rows.values() if r["arm"] == "B"]
    nconv = sum(1 for r in rows.values() for s in r["states"].values() if not s["converged"])
    fb = {}
    fb["F1"] = {"fired": bool(B) and not any(r["verdict"] == "UNRELIABLE" for r in B),
                "meaning": "no Arm B condition is UNRELIABLE: the restricted pool never reached the "
                           "boundary, and the magnitude ladder should be re-run on the strongest "
                           "sites with the floors treated as advisory"}
    fb["F2"] = {"fired": bool(B) and all(r["verdict"] == "UNRELIABLE" for r in B),
                "meaning": "every Arm B condition is UNRELIABLE, including n=1 at q=0.25: diagnose "
                           "the setup before concluding anything"}
    fb["F3"] = {"fired": nconv > 0.25 * 2 * len(rows),
                "meaning": "more than a quarter of the runs did not converge: the finding is about "
                           "the protocol at these field strengths, not about the designs"}
    fb["F4"] = {"fired": n_dis > 0,
                "meaning": "the iso-total families disagree, which D31.5 expects and which "
                           "triggers nothing"}

    payload = {"audit_sha": audit_sha, "prereg_sha": sha(a.prereg),
               "D_LW": D_LW, "R_DIP": R_DIP, "GAP_ABS": GAP_ABS, "R_GAP": R_GAP,
               "n_conditions": len(rows), "n_nonconverged_runs": nconv,
               "control": {k: {kk: vv for kk, vv in v.items() if kk != "sha"} for k, v in ctl.items()},
               "conditions": rows, "arm_A": armA, "arm_B": armB,
               "iso_total": iso, "iso_verdict": iso_verdict, "fallback": fb,
               "notes": notes, "gates_failed": gates}
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "d31_values.json").write_text(json.dumps(payload, indent=1, sort_keys=True, default=str))

    print("ARM A  unit charges, n ascending        D_LW %.2f  R_DIP %.1f  GAP_ABS %.1f  R_GAP %.1f"
          % (D_LW, R_DIP, GAP_ABS, R_GAP))
    for n in NS:
        k = "A/n%02d" % n
        if k not in rows: print("   n=%-3d MISSING" % n); continue
        r = rows[k]
        print("   n=%-3d S|q| %6.2f  %-12s %s" % (n, r["sumabsq"], r["verdict"],
                                                  "; ".join(r["flags"][:2]) or ""))
    print("   ENVELOPE  %s%s" % (armA["envelope"] if armA["envelope"] is not None else "none",
                                 "" if armA["monotone"] else "   (NOT MONOTONE: reported per n)"))
    print()
    print("ARM B  magnitude ladder per n")
    print("   n     q=0.25       q=0.50       q=0.75       q=1.00       envelope")
    for n in NS:
        b = armB.get(n, {})
        cells = "".join("%-13s" % b.get("per_q", {}).get("%.2f" % q, "-")[:12] for q in LADDER)
        print("   %-5d %s%s%s" % (n, cells, b.get("envelope") if b.get("envelope") is not None else "none",
                                  "" if b.get("monotone", True) else "  (not monotone)"))
    print()
    print("ISO-TOTAL  %s" % iso_verdict)
    for s, v in iso.items():
        if not v["agree"]:
            print("   S|q| = %-6s %s" % (s, "  ".join("%s=%s" % (m.split('/')[-1], w) for m, w in v["verdicts"].items())))
    print()
    for k in ("F1", "F2", "F3", "F4"):
        print("   %s %s" % (k, "FIRED - " + fb[k]["meaning"] if fb[k]["fired"] else "not fired"))
    if gates:
        print(); print("   GATES FAILED INSIDE THE READER:")
        for g in gates: print("     " + g)
    if a.result and not gates:
        n = write_result(a.result, payload)
        print("\nresult document written, %d lines" % n)
    elif a.result:
        print("\nresult document NOT written: the reader's own gates failed")
    return 1 if gates else 0


def bail(gates, OUT):
    print("   GATES FAILED INSIDE THE READER:")
    for g in gates: print("     " + g)
    return 1


def write_result(path, p):
    import datetime
    A, B, iso, fb = p["arm_A"], p["arm_B"], p["iso_total"], p["fallback"]
    L = []; w = L.append
    w("STAGE 6  DIAGNOSTIC 31  RESULT")
    w("How many charges, and how large, before the model stops being reliable")
    w(datetime.date.today().strftime("%-d %B %Y"))
    w("")
    w("Pre-registered in STAGE6_DIAGNOSTIC31_20261106.txt")
    w("  sha256 %s" % p["prereg_sha"])
    w("Classifier: phase2.2/cs6_audit.py sha256 %s" % p["audit_sha"])
    w("Thresholds read from it, not transcribed: D_LW %.2f e, R_DIP %.1f, GAP_ABS %.1f eV, R_GAP %.1f"
      % (p["D_LW"], p["R_DIP"], p["GAP_ABS"], p["R_GAP"]))
    w("")
    w("")
    w("R31.1  ARM A: THE UNIT-CHARGE FAMILY")
    w("")
    w("  n     Sigma|q|   verdict        what fired")
    for n in NS:
        k = "A/n%02d" % n
        if k not in p["conditions"]: continue
        r = p["conditions"][k]
        w("  %-5d %8.2f   %-14s %s" % (n, r["sumabsq"], r["verdict"], "; ".join(r["flags"]) or "-"))
    w("")
    w("  ENVELOPE IN n: %s" % (A["envelope"] if A["envelope"] is not None else
                               "none formed; the arm is reported per n"))
    if not A["monotone"]:
        w("  NOT MONOTONE: a larger n is RELIABLE where a smaller one is not, so by D31.4 no")
        w("  envelope is formed. That is an admissible outcome and not a failure.")
    w("")
    w("")
    w("R31.2  ARM B: THE MAGNITUDE LADDER AT EACH n")
    w("")
    w("  n      q=0.25         q=0.50         q=0.75         q=1.00         envelope in q")
    for n in NS:
        b = B.get(n) or B.get(str(n)) or {}
        cells = "".join("%-15s" % (b.get("per_q", {}).get("%.2f" % q, "-"))[:15] for q in LADDER)
        w("  %-6d %s%s" % (n, cells, b.get("envelope") if b.get("envelope") is not None else "none"))
    w("")
    w("")
    w("R31.3  DOES THE TOTAL GOVERN IT? THE ISO-TOTAL FAMILIES")
    w("")
    w("  %s" % p["iso_verdict"])
    w("")
    for s, v in iso.items():
        w("  Sigma|q| = %-6s %s" % (s, "   ".join("%s %s" % (m, w2) for m, w2 in v["verdicts"].items())))
        w("  %-17s %s" % ("", "agree" if v["agree"] else "DISAGREE"))
    w("")
    w("  By D31.5, disagreement confirms DIAGNOSTIC 29 on a purpose-built family and means the")
    w("  envelope must be stated in (n, q) and never as a total. Agreement would be reported as")
    w("  agreement on these families and nothing more; it would not overturn DIAGNOSTIC 29, whose")
    w("  counterexample is one charge vector in two geometries.")
    w("")
    w("")
    w("R31.4  THE PRE-REGISTERED ABANDONMENT CONDITIONS")
    w("")
    for k in ("F1", "F2", "F3", "F4"):
        w("  %s  %s" % (k, "FIRED" if fb[k]["fired"] else "not fired"))
        if fb[k]["fired"]: w("      %s" % fb[k]["meaning"])
    w("")
    w("  %d of %d runs did not converge." % (p["n_nonconverged_runs"], 2 * p["n_conditions"]))
    w("")
    w("")
    w("R31.5  WHAT THIS LICENSES")
    w("")
    w("  The envelopes above are the largest tested designs that passed, on these two site pools,")
    w("  at def2-SVPD under CAM-B3LYP. They are not a rule a new design must obey: DIAGNOSTIC 29")
    w("  showed no function of the charges orders relocation, so nothing here may be written into")
    w("  the optimiser as a cap. The withdrawn Sigma q^2 <= 0.58 e^2 stays withdrawn.")
    w("")
    w("  D31.8's six limitations stand as written. In particular DIAGNOSTIC 30 found def2-SVPD NOT")
    w("  converged at Sigma|q| = 5.6433 e under CAM-B3LYP, so every condition here above about 3 e")
    w("  inherits that doubt and no def2-TZVPD arm was run.")
    w("")
    Path(path).write_text("\n".join(L) + "\n")
    return len(L)


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except SystemExit:
        raise
    except Exception as e:
        print("READER FAILED: %s: %s" % (type(e).__name__, e))
        raise SystemExit(1)
