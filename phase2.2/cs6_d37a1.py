#!/usr/bin/env python3
"""DIAGNOSTIC 37, AMENDMENT 1's reader: revDSD-PBEP86-D4, run two ways.

Energies are the last FINAL SINGLE POINT ENERGY; orbital energies come from D28's committed
reader, cs6_d28.py read_out, imported by path, under D28's rule (no quotable orbital energy
from an output in which TRAH engaged or a level shift was set). That is how D37 reads them.

Which arm is entered for revDSD-PBEP86-D4 is decided by the rule of R-A1.4 and by nothing
in the numbers: arm RI if it passed its gates and its outputs say RI was used; otherwise
arm K if it passed its gates; otherwise NOT MEASURED.

D37's own rows are taken from its committed result document and checked against d37.tsv.
They are not recomputed and not changed.

Test 2 for the entry (R-A1.6) takes D36's barrier from D36's committed result document and
requires it to agree with the barrier recomputed from d36a2/B_{R,TS}/job.out, which is where
D36's result says its coupled-cluster arm is read from.

usage: cs6_d37a1.py --res <results-repo> --code <code-repo> --prereg <amendment> --result <file> --tsv <file>
"""
import argparse, hashlib, importlib.util, os, re, subprocess, sys

H2K = 627.5094740631
NEL, NOCC = 118, 59
VDE, BOUND = 0.5123, 0.0
TOL_EH, TOL_KCAL = 1e-6, 0.0010
BASES = [("SV", "def2-SVPD"), ("TZ", "def2-TZVPD")]
CTL_E = {"R": -836.31486034, "TS": -836.27211231}
CTL_BAR = 26.8248
ARMS = [("RI", "the keyword, with RI asked for in an %mp2 block"), ("K", "the keyword alone")]
MP2 = "%mp2 RI true end"
NAME = "revDSD-PBEP86-D4"
D37_RESULT = "STAGE6_DIAGNOSTIC37_RESULT_20261111.txt"
D36_RESULT = "STAGE6_DIAGNOSTIC36_RESULT_20261114.txt"
GUARDED = ("d28", "d33a", "d36", "d36a2", "d37")


def stop(m):
    sys.stderr.write("STOP  %s\n" % m)
    raise SystemExit(1)


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 16), b""):
            h.update(b)
    return h.hexdigest()


def orca_words(t):
    """ORCA's own words for why a run stopped: the INPUT ERROR box if there is one."""
    L = t.splitlines()
    for i, l in enumerate(L):
        if "INPUT ERROR" in l:
            box = [x.strip() for x in L[i:i + 4] if x.strip() and not set(x.strip()) <= set("!")]
            return " / ".join(box)[:200]
    err = [l.strip() for l in L if re.search(r"error|aborting|unrecognized", l, re.I)]
    return (err[-1][:200] if err else "no error line found; the output ends: " + (L[-1].strip()[:120] if L else "(empty)"))


def main():
    ap = argparse.ArgumentParser()
    for a in ("res", "code", "prereg", "result", "tsv"):
        ap.add_argument("--" + a, required=True)
    a = ap.parse_args()
    S6 = os.path.join(a.res, "phase2.2/cs_stage6")
    D = os.path.join(S6, "d37a1")
    p28 = os.path.join(a.code, "phase2.2", "cs6_d28.py")
    sp = importlib.util.spec_from_file_location("d37a1_cs6_d28", p28)
    d28 = importlib.util.module_from_spec(sp)
    sp.loader.exec_module(d28)

    def read(d):
        p = os.path.join(d, "job.out")
        if not os.path.isfile(p):
            return None
        r = d28.read_out(p)
        t = open(p, errors="replace").read()
        r["notconv"] = bool(re.search(r"SCF NOT CONVERGED AFTER|SERIOUS PROBLEM", t))
        r["second_order"] = bool(re.search(r"^\s*ORCA\s+MP2\s*$|MP2 TOTAL ENERGY|MP2 CORRELATION ENERGY", t, re.M))
        r["ri"] = bool(re.search(r"RI-MP2 CORRELATION ENERGY|RI-MP2 ENERGY", t))
        r["text"] = t
        return r

    rep, gates = [], []
    say = rep.append
    say("STAGE 6, DIAGNOSTIC 37, AMENDMENT 1 - RESULT")
    say("revDSD-PBEP86-D4, which D37 could not run, run two ways.")
    say("")
    say("Read by phase2.2/cs6_d37a1.py          sha256 %s" % sha(__file__))
    say("   orbital energies by cs6_d28.py read_out sha256 %s" % sha(p28))
    say("   amendment %-45s sha256 %s" % (os.path.basename(a.prereg), sha(a.prereg)))
    say("")
    say("Bare single points on frame 42436: no charges, no MM region, no wall. Barriers in")
    say("kcal/mol, orbital energies in eV.")
    say("")

    # ---- GA1.1 to GA1.3 and GA1.5
    geo = {leg: os.path.join(S6, "d33a", "sp_42436_%s" % leg, "substrate.xyz") for leg in ("R", "TS")}
    runs = [("CAM_SV_%s" % leg, leg) for leg in ("R", "TS")]
    runs += [("revDSD_%s_%s_%s" % (arm, tag, leg), leg) for arm, _ in ARMS for tag, _ in BASES for leg in ("R", "TS")]
    for name, leg in runs:
        g = os.path.join(D, name, "substrate.xyz")
        if not os.path.isfile(g) or sha(g) != sha(geo[leg]):
            gates.append("GA1.1 %s's geometry is not frame 42436's" % name)
    for leg in ("R", "TS"):
        src = os.path.join(S6, "d28/CAM-B3LYP/def2-SVPD/ctl_%s/job.inp" % leg)
        p = os.path.join(D, "CAM_SV_%s" % leg, "job.inp")
        if not os.path.isfile(p) or sha(src) != sha(p):
            gates.append("GA1.2 CAM_SV_%s's input is not D28's" % leg)
    for tag, basis in BASES:
        for leg in ("R", "TS"):
            c = open(os.path.join(S6, "d28/CAM-B3LYP/%s/ctl_%s/job.inp" % (basis, leg))).read().split("\n")
            pk = os.path.join(D, "revDSD_K_%s_%s" % (tag, leg), "job.inp")
            pr = os.path.join(D, "revDSD_RI_%s_%s" % (tag, leg), "job.inp")
            k = open(pk).read().split("\n") if os.path.isfile(pk) else []
            r = open(pr).read().split("\n") if os.path.isfile(pr) else []
            if len(k) != len(c) or k[1:] != c[1:] or k[:1] == c[:1]:
                gates.append("GA1.3 revDSD_K_%s_%s does not differ from D28's control in its first line alone" % (tag, leg))
            if [l for l in r if l != MP2] != k or r.count(MP2) != 1:
                gates.append("GA1.3 revDSD_RI_%s_%s does not differ from arm K in the one added line alone" % (tag, leg))
    try:
        q = subprocess.run(["git", "-C", a.res, "diff", "--quiet", "HEAD", "--"] +
                           ["phase2.2/cs_stage6/" + x for x in GUARDED])
        if q.returncode != 0:
            gates.append("GA1.5 a committed file under %s differs from HEAD" % ", ".join(GUARDED))
    except OSError:
        gates.append("GA1.5 could not be checked: git did not run")

    C = {}
    for leg in ("R", "TS"):
        r = read(os.path.join(D, "CAM_SV_%s" % leg))
        if r is None or not r["normal"] or r["energy"] is None:
            gates.append("GA1.4 CAM_SV_%s has no normally terminated output" % leg)
        else:
            C[leg] = r
    if gates:
        say("GATES")
        for g in gates:
            say("   FAILED  %s" % g)
        print("\n".join(rep))
        stop("%d gate(s) failed before any reading. Nothing is written." % len(gates))

    say("=" * 78)
    say("R-A1.1 and GA1.4  THE CONTROL. D28's own CAM-B3LYP input at def2-SVPD, run again.")
    say("=" * 78)
    say("")
    for leg in ("R", "TS"):
        e = C[leg]["energy"]
        dd = abs(e - CTL_E[leg])
        say("   def2-SVPD  %-2s  %.8f Eh   committed %.8f   differ %.2e" % (leg, e, CTL_E[leg], dd))
        if dd > TOL_EH:
            gates.append("GA1.4 %s differs from D28's energy by %.2e Eh, over the %.0e limit" % (leg, dd, TOL_EH))
    b = (C["TS"]["energy"] - C["R"]["energy"]) * H2K
    say("   def2-SVPD  barrier %.4f   committed %.4f   differ %.2e" % (b, CTL_BAR, abs(b - CTL_BAR)))
    if abs(b - CTL_BAR) > TOL_KCAL:
        gates.append("GA1.4 the control barrier is %.4f, not %.4f within %.4f" % (b, CTL_BAR, TOL_KCAL))
    say("")
    if gates:
        say("   THE CONTROL FAILED. Nothing below is readable.")
        for g in gates:
            say("   FAILED  %s" % g)
        print("\n".join(rep))
        stop("the control failed.")
    say("   The control holds.")
    say("")

    # ---- GA1.6, per arm
    A = {}
    for arm, what in ARMS:
        R, why = {}, []
        for tag, basis in BASES:
            R[tag] = {}
            for leg in ("R", "TS"):
                k = "revDSD_%s_%s_%s" % (arm, tag, leg)
                r = read(os.path.join(D, k))
                R[tag][leg] = r
                if r is None:
                    why.append("%s produced no output" % k); continue
                if not r["normal"]:
                    why.append("%s did not terminate normally: %s" % (k, orca_words(r["text"]))); continue
                if r["notconv"]:
                    why.append("%s did not converge its SCF" % k)
                if r["nel"] != NEL:
                    why.append("%s reports %s electrons" % (k, r["nel"]))
                if not r["nocc_ok"] or r["nocc"] != NOCC:
                    why.append("%s has %s occupied orbitals" % (k, r["nocc"]))
                if r["energy"] is None:
                    why.append("%s has no final energy" % k)
                if not r["second_order"]:
                    why.append("%s contains no second-order section" % k)
        ok = not why
        ri = None
        if ok:
            flags = [R[tag][leg]["ri"] for tag, _ in BASES for leg in ("R", "TS")]
            ri = "yes" if all(flags) else ("no" if not any(flags) else "in some runs only")
        A[arm] = dict(R=R, why=why, ok=ok, ri=ri, what=what)

    def bar(R, tag):
        return (R[tag]["TS"]["energy"] - R[tag]["R"]["energy"]) * H2K

    def homo(R, tag, leg):
        return R[tag][leg]["homo"]

    def fh(v):
        return "UNAVAILABLE" if v is None else "%+.4f" % v

    say("=" * 78)
    say("GA1.6 and R-A1.2  THE TWO ARMS")
    say("=" * 78)
    say("")
    for arm, what in ARMS:
        x = A[arm]
        say("   ARM %-2s  %s" % (arm, what))
        if not x["ok"]:
            say("      NOT MEASURED")
            for w in x["why"]:
                say("         %s" % w)
            say("")
            continue
        h = homo(x["R"], "SV", "R")
        say("      its outputs say RI was used in the second-order part: %s" % x["ri"])
        say("      barrier   def2-SVPD %.4f   def2-TZVPD %.4f   shift %+.4f"
            % (bar(x["R"], "SV"), bar(x["R"], "TZ"), bar(x["R"], "TZ") - bar(x["R"], "SV")))
        say("      highest occupied orbital, def2-SVPD, reactant: %s eV%s" % (
            fh(h), "" if h is None else "   %s   %+.4f eV from the target" % ("BOUND" if h < BOUND else "NOT BOUND", h + VDE)))
        say("      the same orbital: def2-SVPD TS %s   def2-TZVPD R %s   def2-TZVPD TS %s"
            % (fh(homo(x["R"], "SV", "TS")), fh(homo(x["R"], "TZ", "R")), fh(homo(x["R"], "TZ", "TS"))))
        say("      (a double hybrid's orbital energies are those of its hybrid SCF step)")
        say("")

    say("=" * 78)
    say("R-A1.3  ARM RI MINUS ARM K")
    say("=" * 78)
    say("")
    if A["RI"]["ok"] and A["K"]["ok"]:
        for tag, basis in BASES:
            say("   %-10s barrier %+.4f kcal/mol" % (basis, bar(A["RI"]["R"], tag) - bar(A["K"]["R"], tag)))
        say("   Arm K's outputs say RI was used: %s. If they say yes, the two arms ran the same" % A["K"]["ri"])
        say("   calculation and this difference is not the effect of the RI approximation.")
    else:
        say("   NOT EVALUATED: both arms are needed.")
    say("")

    # ---- R-A1.4, the rule
    say("=" * 78)
    say("R-A1.4  THE ENTRY FOR revDSD-PBEP86-D4, by the rule and not by the numbers")
    say("=" * 78)
    say("")
    if A["RI"]["ok"] and A["RI"]["ri"] == "yes":
        entry, label = "RI", "arm RI: it passed its gates and its outputs say RI was used"
    elif A["K"]["ok"]:
        entry = "K"
        label = "arm K, because arm RI %s. Arm K's outputs say RI was used: %s" % (
            "was not measured" if not A["RI"]["ok"] else "ran but its outputs do not say RI was used (%s)" % A["RI"]["ri"], A["K"]["ri"])
    else:
        entry, label = None, "NOT MEASURED: neither arm qualifies"
    say("   %s." % label)
    say("")

    # ---- R-A1.6, test 2: D36's number, from its result document, required to agree with d36a2
    bench, note = None, None
    p36 = os.path.join(S6, D36_RESULT)
    bR, bT = (read(os.path.join(S6, "d36a2", "B_%s" % leg)) for leg in ("R", "TS"))
    m36 = re.search(r"THE BARRIER IS (\d+\.\d{4}) kcal/mol\.", open(p36, errors="replace").read()) if os.path.isfile(p36) else None
    if m36 is None:
        note = "%s is missing or does not print its barrier." % D36_RESULT
    elif not (bR and bT and bR["normal"] and bT["normal"] and bR["energy"] is not None and bT["energy"] is not None):
        note = "d36a2/B_R and B_TS are not both present and normally terminated."
    else:
        bb = (bT["energy"] - bR["energy"]) * H2K
        if "%.4f" % bb != m36.group(1):
            note = "the barrier recomputed from d36a2/B_{R,TS}, %.4f, is not the %s that %s prints." % (bb, m36.group(1), D36_RESULT)
        else:
            bench = bb

    # ---- R-A1.5, D37's table with the entry added
    say("=" * 78)
    say("R-A1.5  D37's TABLE OF THE TWO TESTS, WITH THE ENTRY ADDED")
    say("=" * 78)
    say("")
    p37, t37 = os.path.join(S6, D37_RESULT), os.path.join(S6, "d37.tsv")
    if not (os.path.isfile(p37) and os.path.isfile(t37)):
        stop("D37's result document or d37.tsv is missing; its rows cannot be read.")
    lines = open(p37).read().splitlines()
    try:
        i0 = next(i for i, l in enumerate(lines) if l.startswith("R37.6"))
        i1 = next(i for i in range(i0, len(lines)) if lines[i].startswith("   functional"))
    except StopIteration:
        stop("D37's result document has no R37.6 table.")
    table = []
    for l in lines[i1:]:
        if not l.strip():
            break
        table.append(l)
    tsv = {}
    for l in open(t37).read().splitlines()[1:]:
        f = l.split("\t")
        tsv[(f[0], f[2])] = f
    for l in table[1:]:
        f = l.split()
        key = (f[0], "def2-SVPD")
        if key not in tsv or tsv[key][8] not in ("committed", "measured"):
            stop("D37's table names %s, which d37.tsv does not carry as measured." % f[0])
        h = float(tsv[key][6])
        if ("BOUND" if h < BOUND else "NOT BOUND") not in l or ("%+.4f eV" % (h + VDE)) not in l:
            stop("D37's table row for %s does not agree with d37.tsv." % f[0])
    say("   The first %d rows are D37's own, copied from %s" % (len(table) - 1, D37_RESULT))
    say("   and checked against d37.tsv. They are not recomputed.")
    say("")
    for l in table:
        say(l)
    if entry:
        x = A[entry]
        h = homo(x["R"], "SV", "R")
        say("   %-18s %-10s %12s %14s %14s"
            % (NAME, "arm " + entry, "UNAVAILABLE" if h is None else ("BOUND" if h < BOUND else "NOT BOUND"),
               "" if h is None else "%+.4f eV" % (h + VDE),
               "not evaluated" if bench is None else "%+.4f" % (bar(x["R"], "TZ") - bench)))
        say("")
        say("   For the barriers beside D37's: %s %.4f at def2-SVPD and" % (NAME, bar(x["R"], "SV")))
        say("   %.4f at def2-TZVPD." % bar(x["R"], "TZ"))
    else:
        say("")
        say("   %s is NOT MEASURED and no row is added." % NAME)
    say("")

    say("=" * 78)
    say("R-A1.6  TEST 2 FOR THE ENTRY, against D36's DLPNO-CCSD(T) barrier")
    say("=" * 78)
    say("")
    if bench is None:
        say("   NOT EVALUATED. %s" % note)
    elif entry is None:
        say("   NOT EVALUATED: there is no entry.")
    else:
        say("   benchmark %.4f kcal/mol, as %s prints it," % (bench, D36_RESULT))
        say("   and recomputed here from d36a2/B_{R,TS}/job.out to the same four decimals.")
        say("   %s, arm %s, def2-TZVPD: %.4f, which is %+.4f from the benchmark." % (NAME, entry, bar(A[entry]["R"], "TZ"), bar(A[entry]["R"], "TZ") - bench))
        say("   D37's own rows carry the benchmark column their document printed; their test 2 is")
        say("   D37's R37.5 and is not evaluated here.")
        say("   DA1.3: for this one functional, test 2 was fixed before its own number existed and")
        say("   after the benchmark's did.")
    say("")
    say("   This does not adopt a level of theory and bears on no design (DA1.4).")

    L = ["arm\tbasis\tE_R_Eh\tE_TS_Eh\tbarrier_kcal\thomo_R_eV\thomo_TS_eV\tRI_said\tstatus\tentered"]
    for arm, _ in ARMS:
        x = A[arm]
        for tag, basis in BASES:
            if x["ok"]:
                L.append("%s\t%s\t%.8f\t%.8f\t%.4f\t%s\t%s\t%s\tmeasured\t%s" % (
                    arm, basis, x["R"][tag]["R"]["energy"], x["R"][tag]["TS"]["energy"], bar(x["R"], tag),
                    fh(homo(x["R"], tag, "R")), fh(homo(x["R"], tag, "TS")), x["ri"], "yes" if arm == entry else "no"))
            else:
                L.append("%s\t%s\t\t\t\t\t\t\tNOT MEASURED\tno" % (arm, basis))
    open(a.tsv, "w").write("\n".join(L) + "\n")
    open(a.result, "w").write("\n".join(rep) + "\n")
    print("\n".join(rep))
    return 0


if __name__ == "__main__":
    sys.exit(main())
