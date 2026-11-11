#!/usr/bin/env python3
"""DIAGNOSTIC 37's reader: six candidate functionals on the record's own two tests.

Bare single points, so the energy is the last FINAL SINGLE POINT ENERGY, the convention the
committed cs6_audit.py uses. Orbital energies come from D28's committed reader, cs6_d28.py
read_out, imported by path, under D28's rule: an output in which TRAH engaged or a level
shift was set yields no quotable orbital energy.

The three functionals already in the record are recomputed here from their committed raw
outputs (d27 for B3LYP, d28 for CAM-B3LYP and wB97X-D3) and are never typed in.

usage: cs6_d37.py --res <results-repo> --code <code-repo> --prereg <file> --result <file> --tsv <file>
"""
import argparse, glob, hashlib, importlib.util, os, re, subprocess, sys

H2K = 627.5094740631
NEL, NOCC = 118, 59
VDE = 0.5123                      # eV, Probe 3, as D28 uses it
BOUND = 0.0                       # eV, D28.4's criterion
TOL_EH, TOL_KCAL = 1e-6, 0.0010   # G37.4
BASES = [("SV", "def2-SVPD"), ("TZ", "def2-TZVPD")]
CTL_E = {("SV", "R"): -836.31486034, ("SV", "TS"): -836.27211231,
         ("TZ", "R"): -837.20938385, ("TZ", "TS"): -837.16842778}
CTL_BAR = {"SV": 26.8248, "TZ": 25.7003}
CAND = [("wB97MV", "wB97M-V"), ("wB97XV", "wB97X-V"), ("PW6B95", "PW6B95-D3BJ"),
        ("M062X", "M06-2X-D3(0)"), ("PWPB95", "PWPB95-D3BJ"), ("revDSD", "revDSD-PBEP86-D4")]
DOUBLE = ("PWPB95", "revDSD")
# the three already in the record: name -> where its bare pair sits, by basis
OLD = [("B3LYP-D3BJ", "d27/%s/ctl_%s"), ("CAM-B3LYP-D3BJ", "d28/CAM-B3LYP/%s/ctl_%s"),
       ("wB97X-D3", "d28/wB97X-D3/%s/ctl_%s")]


def stop(m):
    sys.stderr.write("STOP  %s\n" % m)
    raise SystemExit(1)


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 16), b""):
            h.update(b)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    for a in ("res", "code", "prereg", "result", "tsv"):
        ap.add_argument("--" + a, required=True)
    a = ap.parse_args()
    S6 = os.path.join(a.res, "phase2.2/cs_stage6")
    D = os.path.join(S6, "d37")
    p28 = os.path.join(a.code, "phase2.2", "cs6_d28.py")
    sp = importlib.util.spec_from_file_location("d37_cs6_d28", p28)
    d28 = importlib.util.module_from_spec(sp)
    sp.loader.exec_module(d28)

    def read(d):
        p = os.path.join(d, "job.out")
        if not os.path.isfile(p):
            return None
        r = d28.read_out(p)
        t = open(p, errors="replace").read()
        r["notconv"] = bool(re.search(r"SCF NOT CONVERGED AFTER|SERIOUS PROBLEM", t))
        r["text"] = t
        return r

    rep, gates = [], []
    say = rep.append
    say("STAGE 6, DIAGNOSTIC 37 - RESULT")
    say("Six candidate functionals on the record's own two tests.")
    say("")
    say("Read by phase2.2/cs6_d37.py            sha256 %s" % sha(__file__))
    say("   orbital energies by cs6_d28.py read_out sha256 %s" % sha(p28))
    say("   pre-registration %-38s sha256 %s" % (os.path.basename(a.prereg), sha(a.prereg)))
    say("")
    say("Bare single points on frame 42436: no charges, no MM region, no wall. Barriers in")
    say("kcal/mol, orbital energies in eV.")
    say("")

    # ---- G37.1 to G37.3 and G37.5, on the files as they are now
    geo = {leg: os.path.join(S6, "d33a", "sp_42436_%s" % leg, "substrate.xyz") for leg in ("R", "TS")}
    names = ["CAM"] + [c for c, _ in CAND]
    for tag, basis in BASES:
        ctl = {}
        for leg in ("R", "TS"):
            d = os.path.join(D, "CAM_%s_%s" % (tag, leg))
            src = os.path.join(S6, "d28/CAM-B3LYP/%s/ctl_%s/job.inp" % (basis, leg))
            if not os.path.isfile(os.path.join(d, "job.inp")):
                gates.append("G37.2 CAM_%s_%s has no input" % (tag, leg)); continue
            if sha(src) != sha(os.path.join(d, "job.inp")):
                gates.append("G37.2 CAM_%s_%s's input is no longer D28's" % (tag, leg))
            ctl[leg] = open(os.path.join(d, "job.inp")).read().split("\n")
        for n in names:
            for leg in ("R", "TS"):
                d = os.path.join(D, "%s_%s_%s" % (n, tag, leg))
                g = os.path.join(d, "substrate.xyz")
                if not os.path.isfile(g) or sha(g) != sha(geo[leg]):
                    gates.append("G37.1 %s_%s_%s's geometry is not frame 42436's" % (n, tag, leg))
                if n == "CAM" or leg not in ctl:
                    continue
                p = os.path.join(d, "job.inp")
                t = open(p).read().split("\n") if os.path.isfile(p) else []
                if len(t) != len(ctl[leg]) or t[1:] != ctl[leg][1:] or t[:1] == ctl[leg][:1]:
                    gates.append("G37.3 %s_%s_%s does not differ from the control in its first line alone" % (n, tag, leg))
    try:
        q = subprocess.run(["git", "-C", a.res, "diff", "--quiet", "HEAD", "--", "phase2.2/cs_stage6/d28",
                            "phase2.2/cs_stage6/d33a", "phase2.2/cs_stage6/d36"])
        if q.returncode != 0:
            gates.append("G37.5 a committed file under d28, d33a or d36 differs from HEAD")
    except OSError:
        gates.append("G37.5 could not be checked: git did not run")

    # ---- G37.4, the control
    C = {}
    for tag, basis in BASES:
        for leg in ("R", "TS"):
            r = read(os.path.join(D, "CAM_%s_%s" % (tag, leg)))
            if r is None or not r["normal"] or r["energy"] is None:
                gates.append("G37.4 CAM_%s_%s has no normally terminated output" % (tag, leg))
            else:
                C[(tag, leg)] = r
    if gates:
        say("GATES")
        for g in gates:
            say("   FAILED  %s" % g)
        print("\n".join(rep))
        stop("%d gate(s) failed before any reading. Nothing is written." % len(gates))

    say("=" * 78)
    say("R37.1 and G37.4  THE CONTROL. D28's own CAM-B3LYP inputs, run again.")
    say("=" * 78)
    say("")
    for tag, basis in BASES:
        for leg in ("R", "TS"):
            e = C[(tag, leg)]["energy"]
            dd = abs(e - CTL_E[(tag, leg)])
            say("   %-10s %-2s  %.8f Eh   committed %.8f   differ %.2e" % (basis, leg, e, CTL_E[(tag, leg)], dd))
            if dd > TOL_EH:
                gates.append("G37.4 %s %s differs from D28's energy by %.2e Eh, over the %.0e limit" % (basis, leg, dd, TOL_EH))
        b = (C[(tag, "TS")]["energy"] - C[(tag, "R")]["energy"]) * H2K
        say("   %-10s barrier %.4f   committed %.4f   differ %.2e" % (basis, b, CTL_BAR[tag], abs(b - CTL_BAR[tag])))
        if abs(b - CTL_BAR[tag]) > TOL_KCAL:
            gates.append("G37.4 the %s control barrier is %.4f, not %.4f within %.4f" % (basis, b, CTL_BAR[tag], TOL_KCAL))
    say("")
    if gates:
        say("   THE CONTROL FAILED. Nothing below is readable.")
        for g in gates:
            say("   FAILED  %s" % g)
        print("\n".join(rep))
        stop("the control failed.")
    say("   The control holds, so the pipeline and the geometry are what D37 says they are.")
    say("")

    # ---- the nine functionals
    rows = []      # name, kind, {tag: {leg: run}}, status
    for name, pat in OLD:
        R = {}
        for tag, basis in BASES:
            R[tag] = {leg: read(os.path.join(S6, pat % (basis, leg))) for leg in ("R", "TS")}
        bad = [1 for tag, _ in BASES for leg in ("R", "TS") if R[tag][leg] is None or R[tag][leg]["energy"] is None]
        if bad:
            stop("the committed bare pair of %s could not be read; the comparison rows would be typed in, and are not." % name)
        rows.append((name, "committed", R, "committed"))
    for c, name in CAND:
        R, why = {}, []
        for tag, basis in BASES:
            R[tag] = {}
            for leg in ("R", "TS"):
                r = read(os.path.join(D, "%s_%s_%s" % (c, tag, leg)))
                R[tag][leg] = r
                k = "%s_%s_%s" % (c, tag, leg)
                if r is None:
                    why.append("%s produced no output" % k); continue
                if not r["normal"]:
                    err = [l.strip() for l in r["text"].splitlines() if re.search(r"error|unrecognized|unknown|abort", l, re.I)]
                    why.append("%s did not terminate normally%s" % (k, (": " + err[0][:90]) if err else ""))
                    continue
                if r["notconv"]:
                    why.append("%s did not converge its SCF" % k)
                if r["nel"] != NEL:
                    why.append("%s reports %s electrons" % (k, r["nel"]))
                if not r["nocc_ok"] or r["nocc"] != NOCC:
                    why.append("%s has %s occupied orbitals" % (k, r["nocc"]))
                if r["energy"] is None:
                    why.append("%s has no final energy" % k)
                if c in DOUBLE and not why:
                    m = re.findall(r"RI-MP2 CORRELATION ENERGY\s*:?\s*(-?\d+\.\d+)", r["text"])
                    if not m:
                        tail = r["text"][max(r["text"].rfind("SCF CONVERGED"), 0):]
                        seen = [l.strip() for l in tail.splitlines() if "MP2" in l][:6]
                        why.append("G37.7 %s: no 'RI-MP2 CORRELATION ENERGY' line was found. Lines naming MP2: %s"
                                   % (k, " // ".join(seen) if seen else "none"))
        rows.append((name, "candidate", R, "; ".join(why) if why else "measured"))

    def bar(R, tag):
        return (R[tag]["TS"]["energy"] - R[tag]["R"]["energy"]) * H2K

    def homo(R, tag, leg):
        r = R[tag][leg]
        return None if r is None else r["homo"]

    def fh(v):
        return "UNAVAILABLE" if v is None else "%+.4f" % v

    ok = [r for r in rows if r[3] in ("committed", "measured")]
    failed = [r for r in rows if r[3] not in ("committed", "measured")]

    say("=" * 78)
    say("G37.6 and G37.7  WHICH CANDIDATES WERE MEASURED")
    say("=" * 78)
    say("")
    for name, kind, R, st in rows:
        if kind == "candidate":
            say("   %-18s %s" % (name, "measured" if st == "measured" else "NOT MEASURED"))
            if st != "measured":
                for w in st.split("; "):
                    say("        %s" % w)
    say("")

    say("=" * 78)
    say("R37.2  TEST 1. The highest occupied orbital, against D28's criterion and target.")
    say("=" * 78)
    say("")
    say("   BOUND if below %.1f eV at def2-SVPD, reactant (D28.4). Target %.4f eV (Probe 3)." % (BOUND, -VDE))
    say("   The double hybrids' orbital energies are those of their hybrid SCF step.")
    say("")
    say("   %-18s %-10s %12s %10s %12s %12s %12s %12s" % ("functional", "", "SVPD, R", "", "from target", "SVPD, TS", "TZVPD, R", "TZVPD, TS"))
    for name, kind, R, st in ok:
        h = homo(R, "SV", "R")
        say("   %-18s %-10s %12s %10s %12s %12s %12s %12s"
            % (name, kind, fh(h), "" if h is None else ("BOUND" if h < BOUND else "NOT BOUND"),
               "" if h is None else "%+.4f" % (h + VDE), fh(homo(R, "SV", "TS")), fh(homo(R, "TZ", "R")), fh(homo(R, "TZ", "TS"))))
    say("")

    say("=" * 78)
    say("R37.3 and R37.4  The barrier at both basis sets, and the shift between them.")
    say("=" * 78)
    say("")
    say("   %-18s %-10s %12s %12s %12s" % ("functional", "", "def2-SVPD", "def2-TZVPD", "shift"))
    for name, kind, R, st in ok:
        say("   %-18s %-10s %12.4f %12.4f %+12.4f" % (name, kind, bar(R, "SV"), bar(R, "TZ"), bar(R, "TZ") - bar(R, "SV")))
    say("")

    # ---- R37.5, only if D36 is committed
    say("=" * 78)
    say("R37.5  TEST 2. Each def2-TZVPD barrier against D36's DLPNO-CCSD(T) barrier.")
    say("=" * 78)
    say("")
    bench, note = None, None
    d36res = sorted(glob.glob(os.path.join(S6, "STAGE6_DIAGNOSTIC36_RESULT_*.txt")))
    bR, bT = (read(os.path.join(S6, "d36", "B_%s" % leg)) for leg in ("R", "TS"))
    if not d36res:
        note = "D36 has no committed result document."
    elif not (bR and bT and bR["normal"] and bT["normal"] and bR["energy"] is not None and bT["energy"] is not None):
        note = "D36's coupled-cluster outputs are not both present and normally terminated."
    else:
        b = (bT["energy"] - bR["energy"]) * H2K
        if ("%.4f" % b) not in open(d36res[-1], errors="replace").read():
            note = ("the barrier recomputed from d36/B_{R,TS}, %.4f, is not printed in %s; the two must agree before it is used."
                    % (b, os.path.basename(d36res[-1])))
        else:
            bench = b
    if bench is None:
        say("   NOT EVALUATED. %s" % note)
        say("   By R37.5 it is added by a dated addendum once D36 is committed.")
    else:
        say("   benchmark %.4f, recomputed from d36/B_{R,TS}/job.out and found in %s" % (bench, os.path.basename(d36res[-1])))
        say("   (sha256 %s). Whether D36's gates passed, and whether its D36.5 caveat" % sha(d36res[-1]))
        say("   applies, is stated in that document and is not re-derived here.")
        say("")
        say("   %-18s %-10s %12s %12s" % ("functional", "", "def2-TZVPD", "minus bench"))
        for name, kind, R, st in sorted(ok, key=lambda r: abs(bar(r[2], "TZ") - bench)):
            say("   %-18s %-10s %12.4f %+12.4f" % (name, kind, bar(R, "TZ"), bar(R, "TZ") - bench))
    say("")

    say("=" * 78)
    say("R37.6  The two tests side by side. No score combines them and none is defined.")
    say("=" * 78)
    say("")
    say("   %-18s %-10s %12s %14s %14s" % ("functional", "", "orbital", "from target", "from benchmark"))
    for name, kind, R, st in ok:
        h = homo(R, "SV", "R")
        say("   %-18s %-10s %12s %14s %14s"
            % (name, kind, "UNAVAILABLE" if h is None else ("BOUND" if h < BOUND else "NOT BOUND"),
               "" if h is None else "%+.4f eV" % (h + VDE),
               "not evaluated" if bench is None else "%+.4f" % (bar(R, "TZ") - bench)))
    say("")
    say("   This does not adopt a level of theory (D37.1). It bears on no design (D37.8).")
    if failed:
        say("   %d of 6 candidates were not measured and are absent from the tables above." % len(failed))

    L = ["functional\tkind\tbasis\tE_R_Eh\tE_TS_Eh\tbarrier_kcal\thomo_R_eV\thomo_TS_eV\tstatus"]
    for name, kind, R, st in rows:
        for tag, basis in BASES:
            if st in ("committed", "measured"):
                L.append("%s\t%s\t%s\t%.8f\t%.8f\t%.4f\t%s\t%s\t%s" % (
                    name, kind, basis, R[tag]["R"]["energy"], R[tag]["TS"]["energy"], bar(R, tag),
                    fh(homo(R, tag, "R")), fh(homo(R, tag, "TS")), st))
            else:
                L.append("%s\t%s\t%s\t\t\t\t\t\tNOT MEASURED" % (name, kind, basis))
    open(a.tsv, "w").write("\n".join(L) + "\n")
    open(a.result, "w").write("\n".join(rep) + "\n")
    print("\n".join(rep))
    return 0


if __name__ == "__main__":
    sys.exit(main())
