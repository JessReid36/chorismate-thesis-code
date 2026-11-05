#!/usr/bin/env python3
"""DIAGNOSTIC 29. Re-derive the ordering of density relocation, under the committed classifier's
RELOCATION criterion alone, after Sigma q^2 <= 0.58 e^2 was withdrawn.

Pre-registered in STAGE6_DIAGNOSTIC29_20261105.txt. This reader implements that document and
nothing else. In particular it does NOT read the HOMO-LUMO gap or the SCF cycle count: both enter
the classifier's full verdict, and DIAGNOSTIC 28 showed both to be functional-dependent.
"""
import argparse, hashlib, importlib.util, json, re
from pathlib import Path

NBASIS = {402: "def2-SVPD", 682: "def2-TZVPD"}
# ORCA's own printed exchange functional, so the functional is verified and not taken from the path
XPRINT = {"B88": "B3LYP", "WB88": "CAM-B3LYP", "WB97X-D3": "wB97X-D3"}


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""): h.update(b)
    return h.hexdigest()


def load_audit(code_root):
    """the COMMITTED classifier: its D_LW and R_DIP, and its own Loewdin and dipole readers."""
    p = Path(code_root) / "phase2.2" / "cs6_audit.py"
    spec = importlib.util.spec_from_file_location("_cs6_audit", p)
    mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
    return mod, sha(p)


def read_pc(p):
    """point charges: Sigma q^2, Sigma |q|, max |q|, n, and the site geometry's own fingerprint."""
    rows = []
    for ln in Path(p).read_text().splitlines()[1:]:
        f = ln.split()
        if len(f) >= 4:
            try: rows.append((float(f[0]), f[1], f[2], f[3]))
            except ValueError: pass
    q = [r[0] for r in rows]
    sites = "|".join(sorted("%s,%s,%s" % r[1:] for r in rows))
    return {"n": len(q), "sumq2": sum(x * x for x in q), "sumabsq": sum(abs(x) for x in q),
            "maxabsq": max((abs(x) for x in q), default=0.0), "sumq": sum(q),
            "arrangement": hashlib.sha256(sites.encode()).hexdigest()[:16]}


def classify_dir(d, audit):
    """one job directory, or None if it is not an admissible condition or control.
    Only maxlw and dip are taken from the classifier. The gap and the cycle count are not read."""
    jo = d / "job.out"
    if not jo.exists(): return None
    t = jo.read_text(errors="replace")
    if "FINAL SINGLE POINT ENERGY" not in t: return None
    nb = re.findall(r"Number of basis functions\s*\.*\s*(\d+)", t)
    if not nb or int(nb[0]) not in NBASIS: return None
    if "LOEWDIN ATOMIC CHARGES" not in t or "Magnitude (Debye)" not in t: return None
    xf = re.findall(r"Exchange Functional\s+Exchange\s*\.*\s*(\S+)", t)
    fn = XPRINT.get(xf[0]) if xf else None
    if fn is None: return None
    st = "TS" if d.name.endswith("_TS") else ("R" if d.name.endswith("_R") else None)
    if st is None: return None
    c = audit.read(t)
    if c.get("maxlw") is None or c.get("dip") is None: return None
    return {"path": str(d), "fn": fn, "basis": NBASIS[int(nb[0])], "state": st,
            "maxlw": float(c["maxlw"]), "dip": float(c["dip"])}


def main():
    a = argparse.ArgumentParser()
    a.add_argument("--tree", required=True, help="phase2.2/cs_stage6 in the results repository")
    a.add_argument("--code", required=True, help="the code repository, for the committed classifier")
    a.add_argument("--out", required=True)
    a = a.parse_args()
    audit, audit_sha = load_audit(a.code)
    T, OUT = Path(a.tree), Path(a.out)
    gates, notes = [], []

    # ---- 1. gather. A condition carries a pc.pc. A directory WITHOUT one is not thereby
    # charge-free: the ECP-walled arms carry their charges as Ne> lines in job.inp, so a
    # charge-free control is identified by name (ctl_*, or d20's lam000 zero-charge arm) AND by
    # the absence of charges in both places.
    conds, ctls = [], []
    # The whole tree, not a glob over directory names. An earlier version of this reader walked
    # T.glob("d*/**/"), which silently excluded every top-level directory not beginning with "d" -
    # transfer, additivity, negatives, negatives_armA among them - and so dropped 72 of 182
    # conditions, including the pair that settles this diagnostic. Erratum 1 records that.
    seen = set()
    cand = sorted({p.parent for p in T.rglob("pc.pc")} | {p.parent for p in T.rglob("job.out")})
    for d in cand:
        if d in seen: continue
        seen.add(d)
        r = classify_dir(d, audit)
        if r is None: continue
        pc = d / "pc.pc"
        inp = (d / "job.inp").read_text(errors="replace") if (d / "job.inp").exists() else ""
        has_ecp_charges = bool(re.search(r"^\s*\w+>\s", inp, re.M))
        if pc.exists():
            r.update(read_pc(pc)); r["kind"] = "condition"; conds.append(r)
        elif not has_ecp_charges and (d.name.startswith("ctl") or "lam000" in d.name):
            r["kind"] = "control"; ctls.append(r)
    notes.append("conditions with a pc.pc in a diffuse basis, both blocks present: %d" % len(conds))
    notes.append("charge-free controls found: %d" % len(ctls))

    # ---- 2. one canonical control per (functional, basis, state), and PROVE the duplicates agree
    by_key = {}
    for c in ctls: by_key.setdefault((c["fn"], c["basis"], c["state"]), []).append(c)
    CTL, spread = {}, []
    for k, v in sorted(by_key.items()):
        mw = [x["maxlw"] for x in v]; dp = [x["dip"] for x in v]
        dm, dd = max(mw) - min(mw), max(dp) - min(dp)
        spread.append({"key": "|".join(k), "n": len(v), "maxlw_spread": dm, "dip_spread": dd,
                       "paths": [x["path"].split("cs_stage6/")[-1] for x in v]})
        if dm > 1e-3 or dd > 1e-2:
            gates.append("GATE C: the %s controls do not agree (maxlw spread %.4f, dip spread %.4f); "
                         "they are not the same control and cannot be pooled" % ("|".join(k), dm, dd))
        CTL[k] = {"maxlw": sum(mw) / len(mw), "dip": sum(dp) / len(dp), "n": len(v)}

    # ---- 3. the criterion, and nothing else
    for c in conds:
        k = (c["fn"], c["basis"], c["state"])
        ctl = CTL.get(k)
        if ctl is None:
            c["moved"] = None; c["why"] = "no same-functional, same-basis, same-state control"
            continue
        c["rise"] = c["maxlw"] - ctl["maxlw"]
        c["ratio"] = c["dip"] / ctl["dip"] if ctl["dip"] else None
        by_lw = c["rise"] >= audit.D_LW
        by_dip = (c["ratio"] is not None and c["ratio"] >= audit.R_DIP)
        c["moved"] = bool(by_lw or by_dip)
        c["why"] = ("loewdin" if by_lw else "") + ("+dipole" if by_dip else "")
    used = [c for c in conds if c["moved"] is not None]
    orphan = [c for c in conds if c["moved"] is None]
    if orphan: notes.append("conditions with no matching control, excluded: %d" % len(orphan))

    primary = [c for c in used if c["fn"] == "B3LYP"]
    held = [c for c in used if c["fn"] != "B3LYP"]
    notes.append("primary set (B3LYP): %d    held back for the cross-functional check: %d"
                 % (len(primary), len(held)))
    notes.append("distinct arrangements in the primary set: %d"
                 % len(set(c["arrangement"] for c in primary)))

    # ---- 4. the decision rule, exactly as pre-registered
    SCALARS = [("sumq2", "Sigma q^2"), ("sumabsq", "Sigma |q|"), ("maxabsq", "max |q|")]
    sep = {}
    for key, label in SCALARS:
        still = [c[key] for c in primary if not c["moved"]]
        move = [c[key] for c in primary if c["moved"]]
        if not still or not move:
            sep[key] = {"label": label, "verdict": "NOT TESTABLE",
                        "n_still": len(still), "n_moved": len(move)}
            continue
        C, B = max(still), min(move)
        sep[key] = {"label": label, "C": C, "B": B, "margin": B - C,
                    "rel_margin": (B - C) / C if C else None,
                    "verdict": "SEPARATES" if B > C else "OVERLAPS",
                    "overlap": (C - B) if B <= C else None,
                    "n_still": len(still), "n_moved": len(move)}
    separating = [k for k, v in sep.items() if v["verdict"] == "SEPARATES"]
    if not separating:
        verdict = "NO SCALAR IN THE SET ORDERS RELOCATION"
        preferred = None
    else:
        preferred = max(separating, key=lambda k: (sep[k]["rel_margin"] if sep[k]["rel_margin"]
                                                   is not None else -1e9))
        verdict = "SEPARATES: " + ", ".join(sep[k]["label"] for k in separating)

    # ---- 5. monotonicity, able to fail on its own
    mono = {"applicable": preferred is not None, "groups": [], "verdict": None}
    if preferred:
        bad = 0
        for gk in sorted(set((c["basis"], c["state"], c["fn"]) for c in primary)):
            g = sorted([c for c in primary if (c["basis"], c["state"], c["fn"]) == gk],
                       key=lambda c: c[preferred])
            rises = [c["rise"] for c in g]
            ok = all(rises[i + 1] >= rises[i] - 1e-9 for i in range(len(rises) - 1))
            mono["groups"].append({"group": "|".join(gk), "n": len(g), "non_decreasing": ok})
            if not ok: bad += 1
        mono["verdict"] = "MONOTONE" if bad == 0 else "NOT MONOTONE in %d of %d groups" % (
            bad, len(mono["groups"]))

    # ---- 6. the cross-functional check, also separate
    cross = {"applicable": preferred is not None, "conditions": [], "verdict": None}
    if preferred:
        C, B = sep[preferred]["C"], sep[preferred]["B"]
        bad = []
        for c in sorted(held, key=lambda c: (c["fn"], c["basis"], c["state"], c[preferred])):
            expect = None
            if c[preferred] <= C: expect = False      # below the boundary: should not have moved
            elif c[preferred] >= B: expect = True     # at or above it: should have moved
            consistent = (expect is None) or (c["moved"] == expect)
            cross["conditions"].append({
                "path": c["path"].split("cs_stage6/")[-1], "fn": c["fn"], "state": c["state"],
                preferred: c[preferred], "moved": c["moved"], "rise": c["rise"],
                "ratio": c["ratio"], "expected": expect,
                "status": "CONSISTENT" if consistent else "INCONSISTENT"})
            if not consistent: bad.append(c["path"].split("cs_stage6/")[-1])
        cross["verdict"] = "CONSISTENT" if not bad else "INCONSISTENT in %d condition(s)" % len(bad)
        cross["inconsistent"] = bad

    payload = {"audit_sha": audit_sha, "D_LW": audit.D_LW, "R_DIP": audit.R_DIP,
               "n_conditions": len(conds), "n_used": len(used), "n_primary": len(primary),
               "n_held": len(held), "control_spread": spread, "separation": sep,
               "separating": separating, "preferred": preferred, "verdict": verdict,
               "monotonicity": mono, "cross_functional": cross, "notes": notes,
               "gates_failed": gates,
               "conditions": [{k: v for k, v in c.items() if k != "path"} |
                              {"path": c["path"].split("cs_stage6/")[-1]} for c in used]}
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "d29_values.json").write_text(json.dumps(payload, indent=1, sort_keys=True))
    for n in notes: print("   " + n)
    print()
    for key, label in SCALARS:
        v = sep[key]
        if v["verdict"] in ("NOT TESTABLE",):
            print("   %-12s %s (did not move: %d, moved: %d)" % (label, v["verdict"], v["n_still"], v["n_moved"]))
        else:
            print("   %-12s C = %.4f   B = %.4f   margin %+.4f   %s" % (label, v["C"], v["B"], v["margin"], v["verdict"]))
    print()
    print("   VERDICT      %s" % verdict)
    if preferred: print("   PREFERRED    %s (relative margin %.4f)" % (sep[preferred]["label"], sep[preferred]["rel_margin"]))
    print("   MONOTONE     %s" % (mono["verdict"] or "not applicable"))
    print("   CROSS-FUNC   %s" % (cross["verdict"] or "not applicable"))
    if gates:
        print(); print("   GATES FAILED INSIDE THE READER:")
        for g in gates: print("     " + g)
    return 1 if gates else 0


if __name__ == "__main__":
    raise SystemExit(main())
