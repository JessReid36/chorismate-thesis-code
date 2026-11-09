#!/usr/bin/env python3
"""DIAGNOSTIC 32 reader: the walled (n, q) grid, PAIRED against DIAGNOSTIC 31.

Every threshold comes from the committed classifier, loaded by explicit path. Relocation is judged
on the Loewdin leg ALONE, as D32.3 declared in advance: ORCA counts an ECP centre as a nucleus, so
a walled arm's dipole is bookkeeping. The dipole ratio is recorded for every condition and gates
nothing.

The pairing recomputes each d31 twin's verdict UNDER THE SAME LOEWDIN-ONLY RULE rather than taking
its published verdict, so the two sides are judged alike. The published verdict is recorded beside
it, and a gate requires that no d31 condition was flagged on the dipole alone - if one were, the
two rules would differ and the pairing would be unfair.
"""
import argparse, hashlib, importlib.util, json, re
from pathlib import Path

NS = (1, 2, 3, 5, 10, 20)
LADDER = (0.25, 0.50, 0.75, 1.00)
WALL0_NS = (1, 5, 20)
H2KCAL = 627.5094740631
SHIFT = [r"Leaving SCF to start the TRAH-SCF procedure", r"\(TRAH MAcro\)", r"[Ll]evel ?[Ss]hift set to"]
# the twelve conditions that answer the +-1 question, and which job they are in
JOB_A = ["A/n%02d" % n for n in NS] + ["B/n%02d_q1.00" % n for n in NS]


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
        blocks.append({"n_negative": sum(1 for e in eig if e < 0), "verdict": v})
    return {"ran": "WAVEFUNCTION STABILITY ANALYSIS" in t, "n_blocks": len(blocks),
            "first": blocks[0]["verdict"] if blocks else None,
            "n_negative": blocks[0]["n_negative"] if blocks else None}


def read_run(d, audit):
    jo = Path(d) / "job.out"
    if not jo.exists(): return None
    t = jo.read_text(errors="replace")
    e = re.findall(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", t)
    c = audit.read(t)
    nb = re.search(r"Number of basis functions\s*\.*\s*(\d+)", t)
    nel = re.search(r"Number of Electrons\s*NEL\s*\.*\s*(\d+)", t)
    return {"normal": "ORCA TERMINATED NORMALLY" in t,
            "converged": bool(e) and not (re.search(r"MAXIMUM NUMBER OF .*ITERATIONS", t, re.I)
                                          or "SCF NOT CONVERGED AFTER" in t),
            "E": float(e[-1]) if e else None, "gap": c.get("gap"),
            "maxlw": c.get("maxlw"), "dip": c.get("dip"),
            "shifted": any(re.search(p, t) for p in SHIFT), "stability": stability(t),
            "nbasis": int(nb.group(1)) if nb else None, "nel": int(nel.group(1)) if nel else None}


def judge(rows, ctl, audit, loewdin_only):
    """the verdict for one condition, from its two states."""
    flags, held = [], False
    for st, r in rows.items():
        c = ctl[st]
        if not r["converged"]: flags.append("%s did not converge" % st)
        if not r["normal"]: flags.append("%s abnormal termination" % st)
        rise = r["maxlw"] - c["maxlw"] if (r["maxlw"] is not None and c["maxlw"] is not None) else None
        ratio = r["dip"] / c["dip"] if (r["dip"] is not None and c["dip"]) else None
        if rise is not None and rise >= audit.D_LW: flags.append("%s relocation: rise %+.4f" % (st, rise))
        if not loewdin_only and ratio is not None and ratio >= audit.R_DIP:
            flags.append("%s relocation: dipole %.2f" % (st, ratio))
        gap_read = (not r["shifted"]) and r["gap"] is not None and c["gap"] is not None
        if not gap_read: held = True
        elif r["gap"] < audit.GAP_ABS and r["gap"] < audit.R_GAP * c["gap"]:
            flags.append("%s gap %.4f eV" % (st, r["gap"]))
        if r["stability"]["first"] == "UNSTABLE":
            flags.append("%s unstable, %s negative" % (st, r["stability"]["n_negative"]))
    v = "UNRELIABLE" if flags else ("INCONCLUSIVE" if held else "RELIABLE")
    return v, sorted(set(flags)), held


def main():
    a = argparse.ArgumentParser()
    a.add_argument("--d32", required=True); a.add_argument("--d31", required=True)
    a.add_argument("--code", required=True); a.add_argument("--prereg", required=True)
    a.add_argument("--out", required=True); a.add_argument("--result")
    a.add_argument("--allow-partial", action="store_true",
                   help="read whatever is present and refuse to write a result document")
    a = a.parse_args()
    audit, audit_sha = load_audit(a.code)
    W, B, OUT = Path(a.d32), Path(a.d31), Path(a.out)
    gates, notes = [], []

    # the baseline is d31's own charge-free control, re-read not re-run, so every rise here sits on
    # the same baseline as every rise there
    ctl = {}
    for st in ("R", "TS"):
        r = read_run(B / ("ctl_%s" % st), audit)
        if r is None: gates.append("d31's ctl_%s is missing; it is the baseline for every rise" % st)
        ctl[st] = r
    if gates:
        for g in gates: print("   " + g)
        return 1

    # ---- the fairness check the pairing depends on ------------------------
    dip_only = []
    for c in ["A/n%02d" % n for n in NS] + ["B/n%02d_q%.2f" % (n, q) for n in NS for q in LADDER]:
        lw = dp = False
        ok = True
        for st in ("R", "TS"):
            r = read_run(B / ("%s_%s" % (c, st)), audit)
            if r is None: ok = False; break
            if r["maxlw"] - ctl[st]["maxlw"] >= audit.D_LW: lw = True
            if ctl[st]["dip"] and r["dip"] / ctl[st]["dip"] >= audit.R_DIP: dp = True
        if ok and dp and not lw: dip_only.append(c)
    if dip_only:
        gates.append("%d d31 condition(s) were flagged on the dipole leg alone (%s); the "
                     "Loewdin-only rule used here would not reproduce their published verdict, so "
                     "the pairing is not like for like" % (len(dip_only), dip_only[:3]))

    CONDS = (["A/n%02d" % n for n in NS]
             + ["B/n%02d_q%.2f" % (n, q) for n in NS for q in LADDER])
    rows, missing = {}, []
    for c in CONDS:
        got = {st: read_run(W / ("%s_%s" % (c, st)), audit) for st in ("R", "TS")}
        if any(v is None for v in got.values()): missing.append(c); continue
        v32, f32, held32 = judge(got, ctl, audit, True)
        src = {st: read_run(B / ("%s_%s" % (c, st)), audit) for st in ("R", "TS")}
        if any(v is None for v in src.values()):
            gates.append("%s has no d31 twin to pair against" % c); continue
        v31, f31, _ = judge(src, ctl, audit, True)
        if v32 == "INCONCLUSIVE" or v31 == "INCONCLUSIVE": pr = "INCONCLUSIVE"
        elif v31 == "UNRELIABLE" and v32 == "RELIABLE": pr = "RESCUED"
        elif v31 == "UNRELIABLE" and v32 == "UNRELIABLE": pr = "STILL UNRELIABLE"
        elif v31 == "RELIABLE" and v32 == "UNRELIABLE": pr = "BROKEN BY THE WALL"
        else: pr = "UNCHANGED"
        rise32 = max(got[st]["maxlw"] - ctl[st]["maxlw"] for st in ("R", "TS"))
        rise31 = max(src[st]["maxlw"] - ctl[st]["maxlw"] for st in ("R", "TS"))
        rows[c] = {"d32": v32, "d31": v31, "pairing": pr, "flags_d32": f32, "flags_d31": f31,
                   "rise_d32": rise32, "rise_d31": rise31,
                   "factor": (rise31 / rise32) if abs(rise32) > 1e-9 else None,
                   "dip_d32": {st: (got[st]["dip"] / ctl[st]["dip"] if ctl[st]["dip"] else None)
                               for st in ("R", "TS")},
                   "nbasis": [got[st]["nbasis"] for st in ("R", "TS")],
                   "nel": [got[st]["nel"] for st in ("R", "TS")],
                   "barrier": ((got["TS"]["E"] - got["R"]["E"]) * H2KCAL
                               if all(got[st]["E"] is not None for st in ("R", "TS")) else None)}
    # the wall-only controls
    w0 = {}
    for n in WALL0_NS:
        got = {st: read_run(W / ("wall0/n%02d_%s" % (n, st)), audit) for st in ("R", "TS")}
        if any(v is None for v in got.values()): missing.append("wall0/n%02d" % n); continue
        rise = max(got[st]["maxlw"] - ctl[st]["maxlw"] for st in ("R", "TS"))
        w0[n] = {"rise": rise, "moved": rise >= audit.D_LW,
                 "stability": {st: got[st]["stability"]["first"] for st in ("R", "TS")},
                 "nbasis": [got[st]["nbasis"] for st in ("R", "TS")],
                 "nel": [got[st]["nel"] for st in ("R", "TS")]}
    if missing:
        notes.append("not present, so not read: %s" % ", ".join(sorted(missing)))

    # GATE 6: the wall added neither electrons nor basis functions
    badnb = [c for c, r in list(rows.items()) + [("wall0/n%02d" % k, v) for k, v in w0.items()]
             if r["nbasis"] != [402, 402] or r["nel"] != [118, 118]]
    if badnb: gates.append("the wall changed NEL or the basis size in: %s" % badnb[:3])

    # ---- the primary verdict, over the twelve +-1 conditions --------------
    twelve = [c for c in JOB_A if c in rows]
    resc = [c for c in twelve if rows[c]["pairing"] == "RESCUED"]
    unch = [c for c in twelve if rows[c]["pairing"] == "UNCHANGED"]
    still = [c for c in twelve if rows[c]["pairing"] == "STILL UNRELIABLE"]
    inc = [c for c in twelve if rows[c]["pairing"] == "INCONCLUSIVE"]
    if len(twelve) < len(JOB_A):
        primary = "INCOMPLETE: %d of %d +-1 conditions read" % (len(twelve), len(JOB_A))
    elif not resc and not unch: primary = "THE WALL DOES NOT RESCUE THE UNIT CHARGE"
    elif len(resc) + len(unch) == len(twelve): primary = "THE WALL RESCUES THE UNIT CHARGE"
    else:
        def nset(cs): return sorted(set(int(re.search(r"n(\d\d)", c).group(1)) for c in cs))
        primary = ("THE WALL RESCUES IT UP TO SOME n: rescued or unchanged at n = %s, "
                   "still unreliable at n = %s" % (nset(resc + unch), nset(still)))

    broken = [c for c, r in rows.items() if r["pairing"] == "BROKEN BY THE WALL"]
    nonconv = sum(1 for c in rows for st in ("R", "TS")
                  if "did not converge" in " ".join(rows[c]["flags_d32"]))
    fb = {"G1": {"fired": any(v["moved"] for v in w0.values()),
                 "meaning": "a wall-only control moves density: the wall is not a neutral addition, "
                            "and no rescue is claimed at that n"},
          "G2": {"fired": bool(twelve) and not resc,
                 "meaning": "nothing among the twelve is RESCUED: the wall does not solve this, "
                            "DIAGNOSTIC 31's envelope stands, and the optimiser works in partial charges"},
          "G3": {"fired": bool(broken),
                 "meaning": "a condition reliable without the wall is unreliable with it: the wall is "
                            "not a pure improvement and no blanket recommendation is made"},
          "G4": {"fired": nonconv > 0.25 * 2 * max(len(rows), 1),
                 "meaning": "more than a quarter of the runs did not converge: the finding is about "
                            "the protocol with ECP centres at these field strengths"}}

    payload = {"audit_sha": audit_sha, "prereg_sha": sha(a.prereg), "partial": bool(missing),
               "D_LW": audit.D_LW, "R_DIP": audit.R_DIP, "GAP_ABS": audit.GAP_ABS, "R_GAP": audit.R_GAP,
               "dipole_only_in_d31": dip_only, "conditions": rows, "wall0": w0,
               "primary": primary, "rescued": resc, "unchanged": unch, "still": still,
               "inconclusive": inc, "broken": broken, "fallback": fb,
               "notes": notes, "gates_failed": gates}
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "d32_values.json").write_text(json.dumps(payload, indent=1, sort_keys=True, default=str))

    print("PAIRED AGAINST DIAGNOSTIC 31   D_LW %.2f   relocation on the Loewdin leg only" % audit.D_LW)
    print("   d31 conditions flagged on the dipole leg alone: %d  %s"
          % (len(dip_only), "(so the pairing is like for like)" if not dip_only else dip_only))
    print()
    print("THE +-1 QUESTION")
    print("   condition          d31            d32            pairing              rise: bare -> walled")
    for c in JOB_A:
        if c not in rows: print("   %-18s not read yet" % c); continue
        r = rows[c]
        print("   %-18s %-14s %-14s %-20s %+7.4f -> %+7.4f%s"
              % (c, r["d31"], r["d32"], r["pairing"], r["rise_d31"], r["rise_d32"],
                 "   (x%.1f)" % r["factor"] if r["factor"] and r["factor"] > 1 else ""))
    print()
    print("   PRIMARY   %s" % primary)
    if w0:
        print()
        print("WALL-ONLY CONTROLS (argon centres, zero charge)")
        for n, v in sorted(w0.items()):
            print("   n=%-3d rise %+7.4f   %s   stability %s/%s" % (n, v["rise"],
                  "MOVES" if v["moved"] else "does not move", v["stability"]["R"], v["stability"]["TS"]))
    rest = [c for c in CONDS if c not in JOB_A and c in rows]
    if rest:
        print()
        print("THE REST OF THE LADDER")
        print("   n      q=0.25               q=0.50               q=0.75")
        for n in NS:
            cells = ""
            for q in (0.25, 0.50, 0.75):
                k = "B/n%02d_q%.2f" % (n, q)
                cells += "%-21s" % (rows[k]["pairing"] if k in rows else "-")
            print("   %-6d %s" % (n, cells))
    print()
    for k in ("G1", "G2", "G3", "G4"):
        print("   %s %s" % (k, "FIRED - " + fb[k]["meaning"] if fb[k]["fired"] else "not fired"))
    for nt in notes: print("\n   note: " + nt)
    if gates:
        print(); print("   GATES FAILED INSIDE THE READER:")
        for g in gates: print("     " + g)
        return 1
    if a.result:
        if missing or a.allow_partial:
            print("\nNO RESULT DOCUMENT WRITTEN: %s" %
                  ("not every condition has run yet" if missing else "--allow-partial was given"))
            print("Nothing above is committed. Rerun when both jobs have finished.")
        else:
            n = write_result(a.result, payload)
            print("\nresult document written, %d lines" % n)
    return 0


def write_result(path, p):
    import datetime
    rows, w0, fb = p["conditions"], p["wall0"], p["fallback"]
    L = []; w = L.append
    w("STAGE 6  DIAGNOSTIC 32  RESULT")
    w("Does giving the spheres a size rescue the unit charge?")
    w(datetime.date.today().strftime("%-d %B %Y"))
    w("")
    w("Pre-registered in STAGE6_DIAGNOSTIC32_20261108.txt")
    w("  sha256 %s" % p["prereg_sha"])
    w("Classifier: phase2.2/cs6_audit.py sha256 %s" % p["audit_sha"])
    w("Thresholds read from it: D_LW %.2f e, GAP_ABS %.1f eV, R_GAP %.1f." % (p["D_LW"], p["GAP_ABS"], p["R_GAP"]))
    w("RELOCATION IS JUDGED ON THE LOEWDIN LEG ALONE, as D32.3 declared before the run: ORCA counts")
    w("an ECP centre as a nucleus, so a walled arm's dipole is bookkeeping. The dipole ratio is")
    w("recorded for every condition below and gates nothing.")
    w("")
    w("THE PAIRING IS LIKE FOR LIKE. Each d31 twin's verdict is recomputed here under the same")
    w("Loewdin-only rule rather than taken from its published verdict. %d d31 conditions were"
      % len(p["dipole_only_in_d31"]))
    w("flagged on the dipole leg alone, so the two rules agree and the comparison is unambiguous.")
    w("")
    w("")
    w("R32.1  THE +-1 QUESTION")
    w("")
    w("  condition          d31            d32            pairing")
    for c in JOB_A:
        if c not in rows: continue
        r = rows[c]
        w("  %-18s %-14s %-14s %s" % (c, r["d31"], r["d32"], r["pairing"]))
    w("")
    w("  PRIMARY VERDICT   %s" % p["primary"])
    w("")
    w("  the Loewdin rise, without the wall and with it:")
    w("")
    w("  condition            bare      walled    factor")
    for c in JOB_A:
        if c not in rows: continue
        r = rows[c]
        w("  %-18s %+9.4f %+9.4f %s" % (c, r["rise_d31"], r["rise_d32"],
                                        ("%.1fx" % r["factor"]) if r["factor"] else "-"))
    w("")
    w("  DIAGNOSTIC 23 measured a factor of nine at the strongest field then available. The")
    w("  factors above are reported for comparison and gate nothing.")
    w("")
    w("")
    w("R32.2  THE WALL-ONLY CONTROLS")
    w("")
    w("  argon centres at the same sites, carrying zero charge.")
    w("")
    w("  n      Loewdin rise   verdict              stability R / TS")
    for n, v in sorted(w0.items()):
        w("  %-6d %+12.4f   %-20s %s / %s" % (n, v["rise"], "MOVES" if v["moved"] else "does not move",
                                              v["stability"]["R"], v["stability"]["TS"]))
    w("")
    if any(v["moved"] for v in w0.values()):
        w("  AT LEAST ONE MOVES. By D32.5 the wall is not a neutral addition at that n, the paired")
        w("  comparison there is confounded, and no rescue is claimed for it.")
    else:
        w("  None moves, so at the n values tested the wall does not by itself perturb the density")
        w("  and the paired comparison is not confounded.")
    w("")
    w("")
    w("R32.3  THE REST OF THE LADDER")
    w("")
    w("  n      q=0.25               q=0.50               q=0.75")
    for n in NS:
        cells = ""
        for q in (0.25, 0.50, 0.75):
            k = "B/n%02d_q%.2f" % (n, q)
            cells += "%-21s" % (rows[k]["pairing"] if k in rows else "-")
        w("  %-6d %s" % (n, cells))
    w("")
    if p["broken"]:
        w("  BROKEN BY THE WALL: %s" % ", ".join(p["broken"]))
        w("  By D32.6 G3 the wall is not a pure improvement, and no blanket recommendation to wall")
        w("  everything is made.")
    w("")
    w("")
    w("R32.4  THE PRE-REGISTERED STOPPING CONDITIONS")
    w("")
    for k in ("G1", "G2", "G3", "G4"):
        w("  %s  %s" % (k, "FIRED" if fb[k]["fired"] else "not fired"))
        if fb[k]["fired"]: w("      %s" % fb[k]["meaning"])
    w("")
    w("")
    w("R32.5  WHAT THIS LICENSES")
    w("")
    w("  D32.8's six limitations stand as written. It is one ECP species, chosen in DIAGNOSTIC 23 by")
    w("  matching the Amber nitrogen Lennard-Jones wall; it does not validate the wall against")
    w("  anything physical; no def2-TZVPD arm was run, so the strong-field conditions inherit")
    w("  DIAGNOSTIC 30's basis doubt; one leg of the relocation criterion was unavailable")
    w("  throughout; and it is one frame, one grid, in vacuo, at fixed geometries.")
    w("")
    w("  IT SETS NO BUDGET AND NO CAP. DIAGNOSTIC 29 still rules out a limit written as a function")
    w("  of the charges, and the withdrawn Sigma q^2 <= 0.58 e^2 stays withdrawn.")
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
