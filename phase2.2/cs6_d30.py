#!/usr/bin/env python3
"""DIAGNOSTIC 30 reader. Two arms, both pre-registered in STAGE6_DIAGNOSTIC30_20261105.txt.

A  the basis sequence under CAM-B3LYP: ddE, gap_SD, gap_DT, by D27's rule with D27's tolerance.
B  the stability of the CAM-B3LYP solution, by D19 Reading 2's rule.

TOL is parsed out of the pre-registration rather than written here, so the figure the reader uses
and the figure the pre-registration committed cannot drift apart.
"""
import argparse, hashlib, json, re
from pathlib import Path

H2KCAL = 627.5094740631          # the project constant, as cs1_analyse.py and cs13_qscale.py use
VECTORS = ["1.3000e_n10", "2.8800e_n10", "5.6433e_n20"]
STATES = ["R", "TS"]
# the narrow TRAH markers. A bare "TRAH" matches the author credits and the settings listing in
# every ORCA output, so it must not be used.
TRAH = [r"Leaving SCF to start the TRAH-SCF procedure", r"\(TRAH MAcro\)", r"[Ll]evel ?[Ss]hift set to"]


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""): h.update(b)
    return h.hexdigest()


def read_tol(prereg):
    p = Path(prereg)
    if not p.exists():
        raise SystemExit("the pre-registration is not at %s. TOL is read from it and never assumed." % prereg)
    m = re.search(r"TOL = ([\d.]+) kcal/mol", p.read_text())
    if not m: raise SystemExit("could not find TOL in the pre-registration; refusing to assume one")
    return float(m.group(1))


def energy(d):
    """the LAST FINAL SINGLE POINT ENERGY, plus what else the gates need."""
    jo = Path(d) / "job.out"
    if not jo.exists(): return None
    t = jo.read_text(errors="replace")
    fsp = re.findall(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", t)
    if not fsp: return None
    nb = re.search(r"Number of basis functions\s*\.*\s*(\d+)", t)
    return {"E": float(fsp[-1]), "n_fsp": len(fsp),
            "normal": "ORCA TERMINATED NORMALLY" in t,
            "maxiter": bool(re.search(r"MAXIMUM NUMBER OF .*ITERATIONS", t, re.I))
                       or "SCF NOT CONVERGED AFTER" in t,
            "nbasis": int(nb.group(1)) if nb else None,
            "trah": any(re.search(p, t) for p in TRAH),
            "sha": sha(jo)}


def stability(d):
    """every WAVEFUNCTION STABILITY ANALYSIS block, in order, with its eigenvalues and verdict.
    The format is ORCA's own, read from a committed d20 output before this was written:
        The eigenvalues of the stability matrix:
           E( 0) =  -0.00753321 Eh
        The stability analysis indicates that the wavefunction is unstable
        The stability analysis shows that the wavefunction is stable
    """
    t = (Path(d) / "job.out").read_text(errors="replace")
    blocks = []
    for m in re.finditer(r"The eigenvalues of the stability matrix:\s*((?:\s*E\(\s*\d+\)\s*=\s*-?\d+\.\d+ Eh\n)+)", t):
        eig = [float(x) for x in re.findall(r"E\(\s*\d+\)\s*=\s*(-?\d+\.\d+) Eh", m.group(1))]
        tail = t[m.end():m.end() + 400]
        if "wavefunction is unstable" in tail: v = "UNSTABLE"
        elif "wavefunction is stable" in tail: v = "STABLE"
        else: v = None
        blocks.append({"eigenvalues": eig, "n_negative": sum(1 for e in eig if e < 0), "verdict": v})
    return {"n_blocks": len(blocks), "blocks": blocks,
            "ran": "WAVEFUNCTION STABILITY ANALYSIS" in t,
            "restarted": "Restart requested" in t,
            # the first block judges the solution the plain run would have kept; the last judges
            # the one this run finished on
            "first": blocks[0]["verdict"] if blocks else None,
            "last": blocks[-1]["verdict"] if blocks else None}


def write_result(path, pay, prereg_sha, audit_note):
    """the result document: the pre-registered structure, the measured numbers, and the verdicts
    the pre-registered rules give. No interpretation beyond what those rules state."""
    A, B, b3, TOL = pay["arm_A"], pay["arm_B"], pay["b3lyp_comparison"], pay["TOL"]
    L = []
    w = L.append
    w("STAGE 6  DIAGNOSTIC 30  RESULT")
    w("Is def2-SVPD converged under CAM-B3LYP, and is the CAM-B3LYP solution stable?")
    w(audit_note)
    w("")
    w("Pre-registered in STAGE6_DIAGNOSTIC30_20261105.txt")
    w("  sha256 %s" % prereg_sha)
    w("Reader: phase2.2/cs6_d30.py, committed with this document.")
    w("TOL = %.1f kcal/mol, read from the pre-registration and not written into the reader." % TOL)
    w("Energies are the last FINAL SINGLE POINT ENERGY of each run; the conversion is the project")
    w("constant %.10f kcal/mol per Hartree." % pay["H2KCAL"])
    w("")
    w("")
    w("R30.1  ARM A: THE BASIS SEQUENCE UNDER CAM-B3LYP")
    w("")
    w("  ddE(basis, v) = barrier(v) - barrier(control of that basis), in kcal/mol")
    w("")
    w("  vector            def2-SVP    def2-SVPD   def2-TZVPD     gap_SD     gap_DT   verdict")
    for v, r in A["vectors"].items():
        if r.get("status") == "INCOMPLETE":
            w("  %-14s   INCOMPLETE" % v); continue
        w("  %-14s %+10.4f  %+10.4f   %+10.4f  %9.4f  %9.4f   %s"
          % (v, r["ddE"]["def2-SVP"], r["ddE"]["def2-SVPD"], r["ddE"]["def2-TZVPD"],
             r["gap_SD"], r["gap_DT"], r["verdict"]))
        if "arbitration" in r:
            w("  %-14s   arbitration, by D27's rule: %s" % ("", r["arbitration"]))
    w("")
    w("  VERDICT   %s" % A["verdict"])
    w("")
    w("  B3LYP's OWN SEQUENCE, RE-READ FROM d27 AND REPORTED BESIDE, NOT ENTERING THE VERDICT:")
    w("")
    w("  vector            def2-SVP    def2-SVPD   def2-TZVPD     gap_SD     gap_DT")
    for v in A["vectors"]:
        if not b3.get(v): w("  %-14s   not available" % v); continue
        r = b3[v]
        w("  %-14s %+10.4f  %+10.4f   %+10.4f  %9.4f  %9.4f"
          % (v, r["ddE"]["def2-SVP"], r["ddE"]["def2-SVPD"], r["ddE"]["def2-TZVPD"],
             r["gap_SD"], r["gap_DT"]))
    w("")
    w("  These reproduce DIAGNOSTIC 27's committed figures from the raw outputs rather than")
    w("  quoting them, which is also what validates this reader's arithmetic.")
    w("")
    w("")
    w("R30.2  ARM B: THE STABILITY OF THE CAM-B3LYP SOLUTION")
    w("")
    w("  The verdict of the FIRST stability analysis is the one that judges the solution a plain")
    w("  run would have kept. Where a run restarted, the second analysis judges what it finished on.")
    w("")
    w("  condition           first     blocks  negative eigenvalues   dE vs the committed plain run")
    for k, c in sorted(B["conditions"].items()):
        s = c["stability"]
        w("  %-18s  %-9s %-6d %-22s %+9.4f kcal/mol"
          % (k, s["first"] or "?", s["n_blocks"],
             (str(s["blocks"][0]["n_negative"]) if s["blocks"] else "?"),
             c["dE_kcal"] if c["dE_kcal"] is not None else float("nan")))
    w("")
    w("  vector              barrier, committed   barrier, stability-aware   change")
    for k, bb in sorted(B["barriers"].items()):
        w("  %-18s %18.4f %26.4f %+9.4f   %s"
          % (k, bb["barrier_committed"], bb["barrier_stab"], bb["change"],
             "MATERIAL" if bb["material"] else "IMMATERIAL"))
    w("")
    w("  VERDICT   %s" % B["verdict"])
    if B["control_itself_unstable"]:
        w("")
        w("  THE CHARGE-FREE CONTROL IS ITSELF UNSTABLE. By D30.4 this finding is about the method")
        w("  and the bare dianion, not about the designed field, and is not attributed to the charges.")
    w("")
    w("  The lowering is measured against the committed plain def2-SVPD run of the same name,")
    w("  because ORCA does not print the pre-restart SCF energy in a labelled block and because")
    w("  the committed value is the one every committed quantity uses.")
    w("")
    w("")
    w("R30.3  WHAT FOLLOWS, BY THE PRE-REGISTERED RULES AND NOTHING ELSE")
    w("")
    notc = [v for v, r in A["vectors"].items() if r.get("verdict") == "NOT CONVERGED"]
    if notc:
        w("  ARM A. def2-SVPD is not sufficient under CAM-B3LYP at %s." % ", ".join(notc))
        w("  By D30.3, every CAM-B3LYP quantity at or above that Sum|q| is provisional until re-run.")
    elif A["verdict"] == "CONVERGED AT def2-SVPD":
        w("  ARM A. def2-SVPD is the production basis under CAM-B3LYP. Nothing in the record needs")
        w("  re-running on account of the basis.")
    else:
        w("  ARM A. %s. The consequence is taken per vector, as D30.3 states." % A["verdict"])
    w("")
    if B["verdict"] == "ALL STABLE":
        w("  ARM B. No condition reports an instability. The committed CAM-B3LYP energies stand.")
    elif B["verdict"] == "UNSTABLE, IMMATERIAL":
        w("  ARM B. Instabilities are reported but no barrier moves by more than TOL. They are")
        w("  recorded and the committed energies stand.")
    else:
        w("  ARM B. At least one barrier moves by more than TOL: %s." % ", ".join(B["material_vectors"]))
        w("  By D30.4, every CAM-B3LYP barrier at or above that vector is provisional until the")
        w("  reconverged solution is used.")
    w("")
    w("  D30.6's six limitations stand as written. In particular this says nothing about +-1 e")
    w("  designs and nothing about density relocation, and it does not by itself make CAM-B3LYP")
    w("  the production functional.")
    w("")
    if pay["notes"]:
        w("")
        w("NOTES")
        for n in pay["notes"]: w("  - " + n)
        w("")
    Path(path).write_text("\n".join(L) + "\n")
    return len(L)


def main():
    a = argparse.ArgumentParser()
    a.add_argument("--d30", required=True)
    a.add_argument("--d28", required=True, help="d28/CAM-B3LYP, for the def2-SVPD arms and the def2-TZVPD controls")
    a.add_argument("--d27", required=True, help="d27, for B3LYP's own sequence, reported beside")
    a.add_argument("--prereg", required=True)
    a.add_argument("--out", required=True)
    a.add_argument("--result", help="also write the result document to this path")
    a = a.parse_args()
    D30, D28, D27, OUT = Path(a.d30), Path(a.d28), Path(a.d27), Path(a.out)
    TOL = read_tol(a.prereg)
    gates, notes = [], []

    # ---- where each basis's runs live -------------------------------------
    SRC = {"def2-SVP":   {"arm": D30 / "A/def2-SVP",  "ctl": D30 / "A/def2-SVP"},
           "def2-SVPD":  {"arm": D28 / "def2-SVPD",   "ctl": D28 / "def2-SVPD"},
           "def2-TZVPD": {"arm": D30 / "A/def2-TZVPD", "ctl": D28 / "def2-TZVPD"}}
    NB = {"def2-SVP": 264, "def2-SVPD": 402, "def2-TZVPD": 682}

    # ---- ARM A ------------------------------------------------------------
    E, barrier, ddE = {}, {}, {}
    for b, s in SRC.items():
        for st in STATES:
            r = energy(s["ctl"] / ("ctl_%s" % st))
            if r is None: gates.append("ARM A: %s ctl_%s is missing or has no energy" % (b, st)); continue
            E[(b, "ctl", st)] = r
        for v in VECTORS:
            for st in STATES:
                r = energy(s["arm"] / ("%s_%s" % (v, st)))
                if r is None: gates.append("ARM A: %s %s_%s is missing or has no energy" % (b, v, st)); continue
                E[(b, v, st)] = r
    for (b, k, st), r in E.items():
        if r["nbasis"] != NB[b]:
            gates.append("ARM A: %s %s %s reports %s basis functions, expected %d"
                         % (b, k, st, r["nbasis"], NB[b]))
        if not r["normal"]: gates.append("ARM A: %s %s %s did not terminate normally" % (b, k, st))
        if r["maxiter"]: gates.append("ARM A: %s %s %s hit the iteration limit" % (b, k, st))
    for b in SRC:
        for k in ["ctl"] + VECTORS:
            if (b, k, "R") in E and (b, k, "TS") in E:
                barrier[(b, k)] = (E[(b, k, "TS")]["E"] - E[(b, k, "R")]["E"]) * H2KCAL
    for b in SRC:
        for v in VECTORS:
            if (b, v) in barrier and (b, "ctl") in barrier:
                ddE[(b, v)] = barrier[(b, v)] - barrier[(b, "ctl")]

    armA = {"vectors": {}, "verdict": None, "per_vector": True}
    verdicts = {}
    for v in VECTORS:
        got = all((b, v) in ddE for b in SRC)
        if not got:
            armA["vectors"][v] = {"status": "INCOMPLETE"}; verdicts[v] = None; continue
        sd = abs(ddE[("def2-SVPD", v)] - ddE[("def2-SVP", v)])
        dt = abs(ddE[("def2-TZVPD", v)] - ddE[("def2-SVPD", v)])
        shrinking = dt < sd
        if dt <= TOL and shrinking: vv = "CONVERGED"
        elif dt <= TOL: vv = "WITHIN TOL BUT NOT SHRINKING"
        else: vv = "NOT CONVERGED"
        row = {"ddE": {b: ddE[(b, v)] for b in SRC}, "gap_SD": sd, "gap_DT": dt,
               "shrinking": shrinking, "verdict": vv}
        # D27's third-vector arbitration, applied to whichever vector fails
        if vv == "NOT CONVERGED":
            a1 = abs(ddE[("def2-TZVPD", v)] - ddE[("def2-SVPD", v)])
            a2 = abs(ddE[("def2-TZVPD", v)] - ddE[("def2-SVP", v)])
            row["arbitration"] = ("def2-SVP IS THE ARTEFACT" if a1 <= TOL else
                                  "def2-SVPD IS THE ARTEFACT" if a2 <= TOL else "BOTH ARE ARTEFACTS")
        armA["vectors"][v] = row; verdicts[v] = vv
    vs = set(x for x in verdicts.values() if x)
    if len(vs) == 1 and None not in verdicts.values():
        armA["verdict"] = "CONVERGED AT def2-SVPD" if vs == {"CONVERGED"} else vs.pop()
        armA["per_vector"] = False
    else:
        armA["verdict"] = "PER VECTOR; NO SINGLE VERDICT IS FORMED"

    # B3LYP's own sequence, re-read from d27 and reported beside, not entering the verdict
    b3 = {}
    for v in VECTORS:
        try:
            bb = {}
            for b, p in (("def2-SVP", D27 / "def2-SVP"), ("def2-SVPD", D27 / "def2-SVPD"),
                         ("def2-TZVPD", D27 / "def2-TZVPD")):
                er, et = energy(p / ("%s_R" % v)), energy(p / ("%s_TS" % v))
                cr, ct = energy(p / "ctl_R"), energy(p / "ctl_TS")
                if None in (er, et, cr, ct): raise KeyError(b)
                bb[b] = ((et["E"] - er["E"]) - (ct["E"] - cr["E"])) * H2KCAL
            b3[v] = {"ddE": bb, "gap_SD": abs(bb["def2-SVPD"] - bb["def2-SVP"]),
                     "gap_DT": abs(bb["def2-TZVPD"] - bb["def2-SVPD"])}
        except KeyError:
            b3[v] = None
            notes.append("B3LYP comparison unavailable for %s; reported as absent, not assumed" % v)

    # ---- ARM B ------------------------------------------------------------
    armB = {"conditions": {}, "verdict": None}
    for k in ["ctl"] + VECTORS:
        for st in STATES:
            d = D30 / "B/def2-SVPD" / ("%s_%s" % (k, st) if k == "ctl" else "%s_%s" % (k, st))
            e = energy(d)
            if e is None: gates.append("ARM B: %s_%s is missing or has no energy" % (k, st)); continue
            if not e["normal"]: gates.append("ARM B: %s_%s did not terminate normally" % (k, st))
            s = stability(d)
            if not s["ran"]: gates.append("ARM B: %s_%s ran no stability analysis" % (k, st))
            plain = energy(D28 / "def2-SVPD" / ("%s_%s" % (k, st)))
            # the lowering is measured against the COMMITTED plain run, because the pre-restart SCF
            # energy is not printed in a labelled block, and because the committed value is the one
            # every committed quantity uses
            dE = (e["E"] - plain["E"]) * H2KCAL if plain else None
            armB["conditions"]["%s_%s" % (k, st)] = {
                "E_stab": e["E"], "E_committed": plain["E"] if plain else None,
                "dE_kcal": dE, "stability": s, "trah": e["trah"]}
    bars = {}
    for k in ["ctl"] + VECTORS:
        a_, b_ = armB["conditions"].get("%s_R" % k), armB["conditions"].get("%s_TS" % k)
        if a_ and b_ and a_["E_committed"] is not None and b_["E_committed"] is not None:
            new = (b_["E_stab"] - a_["E_stab"]) * H2KCAL
            old = (b_["E_committed"] - a_["E_committed"]) * H2KCAL
            bars[k] = {"barrier_stab": new, "barrier_committed": old, "change": new - old,
                       "material": abs(new - old) > TOL}
    armB["barriers"] = bars
    unstable = sorted(k for k, c in armB["conditions"].items()
                      if c["stability"]["first"] == "UNSTABLE")
    mat = sorted(k for k, b in bars.items() if b["material"])
    armB["unstable"] = unstable
    armB["material_vectors"] = mat
    if not unstable: armB["verdict"] = "ALL STABLE"
    elif not mat: armB["verdict"] = "UNSTABLE, IMMATERIAL"
    else: armB["verdict"] = "UNSTABLE, MATERIAL"
    armB["control_itself_unstable"] = any(k.startswith("ctl") for k in unstable)

    payload = {"TOL": TOL, "H2KCAL": H2KCAL, "prereg_sha": sha(a.prereg),
               "arm_A": armA, "arm_B": armB, "b3lyp_comparison": b3,
               "barriers": {"|".join(k): v for k, v in barrier.items()},
               "ddE": {"|".join(k): v for k, v in ddE.items()},
               "energies": {"|".join(k): v for k, v in E.items()},
               "notes": notes, "gates_failed": gates}
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "d30_values.json").write_text(json.dumps(payload, indent=1, sort_keys=True, default=str))
    if a.result and not gates:
        import datetime
        n = write_result(a.result, payload, payload["prereg_sha"],
                         datetime.date.today().strftime("%-d %B %Y"))
        print("result document written, %d lines" % n)
    elif a.result:
        print("result document NOT written: the reader's own gates failed")

    print("ARM A  the basis sequence under CAM-B3LYP      TOL = %.1f kcal/mol" % TOL)
    for v in VECTORS:
        r = armA["vectors"][v]
        if r.get("status") == "INCOMPLETE": print("   %-14s INCOMPLETE" % v); continue
        print("   %-14s ddE  SVP %+9.4f  SVPD %+9.4f  TZVPD %+9.4f" %
              (v, r["ddE"]["def2-SVP"], r["ddE"]["def2-SVPD"], r["ddE"]["def2-TZVPD"]))
        print("   %-14s gap_SD %8.4f   gap_DT %8.4f   %s%s" %
              ("", r["gap_SD"], r["gap_DT"], r["verdict"],
               "   -> " + r["arbitration"] if "arbitration" in r else ""))
        if b3.get(v):
            print("   %-14s B3LYP, beside: gap_SD %8.4f   gap_DT %8.4f"
                  % ("", b3[v]["gap_SD"], b3[v]["gap_DT"]))
    print("   VERDICT      %s" % armA["verdict"])
    print()
    print("ARM B  the stability of the CAM-B3LYP solution")
    for k, c in sorted(armB["conditions"].items()):
        s = c["stability"]
        eig = s["blocks"][0]["eigenvalues"][:3] if s["blocks"] else []
        print("   %-18s %-8s blocks %d  negatives %s  dE vs committed %+8.4f kcal/mol  %s"
              % (k, s["first"] or "?", s["n_blocks"],
                 s["blocks"][0]["n_negative"] if s["blocks"] else "?",
                 c["dE_kcal"] if c["dE_kcal"] is not None else float("nan"),
                 " ".join("%+.5f" % e for e in eig)))
    for k, b in sorted(bars.items()):
        print("   %-18s barrier %8.4f -> %8.4f   change %+8.4f   %s"
              % (k, b["barrier_committed"], b["barrier_stab"], b["change"],
                 "MATERIAL" if b["material"] else "IMMATERIAL"))
    print("   VERDICT      %s" % armB["verdict"])
    if armB["control_itself_unstable"]:
        print("   NOTE         the charge-free control is itself unstable: this is about the method")
        print("                and the bare dianion, not about the designed field")
    if notes:
        print()
        for n in notes: print("   note: " + n)
    if gates:
        print(); print("   GATES FAILED INSIDE THE READER:")
        for g in gates: print("     " + g)
    return 1 if gates else 0


if __name__ == "__main__":
    # a traceback from a reader once swallowed an entire gate report in this project; any
    # exception becomes a clean failure line instead
    try:
        raise SystemExit(main())
    except SystemExit:
        raise
    except Exception as e:
        print("READER FAILED: %s: %s" % (type(e).__name__, e))
        raise SystemExit(1)
