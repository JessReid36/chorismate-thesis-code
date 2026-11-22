#!/usr/bin/env python3
"""cs6_d49.py - the reader of Stage 6 DIAGNOSTIC 49: the calibration of DIAGNOSTIC 48, run again so that ORCA ends
normally. Committed before any job of D49 is submitted.

It reads the forty runs of D49 with the committed reader of outputs (cs6_d40.py), as D48's reader read D48's
calibration, and forms:

  B(q)   what the def2-TZVPD functions of hydrogen bind at a bare charge q: minus the energy of one electron, in eV,
         at Hartree-Fock (HF) and with the functional of D48's runs (CAM). At each size it is taken from variant N (the
         D48 input with the dipole moment switched off) where that run is measured, and otherwise from variant Z (a
         hydrogen atom whose nuclear charge is q, the charge at the centre of the functions), and it says which.
  S'     for each of D48's 44 runs with a charge: S with 13.605693122994 q_i^2 replaced by B(CAM) at q_i, as D48's
         pre-registration defines it. S is formed again here from the committed files, by D48's formula, and must be
         the S that d48.tsv holds. What each run did is taken from d48.tsv and is not read again.

S' PREDICTS NOTHING: D48's result was known when D49 was written. It is set beside S.

The inputs are not taken on trust: every file of d49 is formed again by cs6_d49_prepare.py from D48's committed files
and compared, every output must print the committed input beside it, and the files of D48 this leans on must be the
committed ones. A gate that fails stops the reader before anything is written.

usage: python3 cs6_d49.py --res <results repository> --code <code repository> --prereg <the pre-registration>
                          --result <the result to write> --tsv <the table to write>
"""
import argparse
import hashlib
import importlib.util
import os
import re
import subprocess
import sys
import textwrap

S = "phase2.2/cs_stage6"
PREREG = "STAGE6_DIAGNOSTIC49_20261121.txt"
D48PRE = "STAGE6_DIAGNOSTIC48_20261121.txt"
D48RES = "STAGE6_DIAGNOSTIC48_RESULT_20261121.txt"
EV = 27.211386245988
RYD = EV / 2.0
BOHR = 0.5291772083
D_LW = 0.25
TOL_S = 6e-5          # eV: S formed here against the S of d48.tsv, which is written to four decimals
NBF = 9               # the def2-TZVPD functions of one hydrogen atom, as D48's outputs report them
# quoted sentences: each must be in the document it is quoted from, and the first two in the pre-registration given
Q_A = "S' PREDICTS NOTHING. It is set beside S"                                                         # D48, section 2
Q_B = "NOT MEASURED at a size: S' is not formed for the arrangements that hold that size"                 # D48, section 6
Q_C = "Calibration measured: 0 of 20."                                                                    # D48's result
GUARDED = ["d31", "d33a", "d37", "d40", "d43", "d44", "d47", "d48"]


def stop(m):
    sys.stderr.write("STOP  %s\n" % m)
    raise SystemExit(1)


def load(path, name):
    sp = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(sp)
    sp.loader.exec_module(m)
    return m


def sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


def echoed(text):
    m = re.search(r"INPUT FILE\n=+\nNAME = [^\n]*\n((?:\|\s*\d+>[^\n]*\n)+)", text)
    if not m:
        return None
    ls = [re.sub(r"^\|\s*\d+> ?", "", l) for l in m.group(1).split("\n") if l]
    while ls and (not ls[-1].strip() or "****END OF INPUT****" in ls[-1]):
        ls.pop()
    return ls


def table(p):
    ls = open(p).read().splitlines()
    h = ls[0].split("\t")
    return [dict(zip(h, l.split("\t"))) for l in ls[1:]]


def ranks(v):
    """ranks from 1, ties given the mean of their places"""
    o = sorted(range(len(v)), key=lambda i: v[i])
    r = [0.0] * len(v)
    i = 0
    while i < len(o):
        j = i
        while j + 1 < len(o) and v[o[j + 1]] == v[o[i]]:
            j += 1
        for k in range(i, j + 1):
            r[o[k]] = (i + j) / 2.0 + 1.0
        i = j + 1
    return r


def spearman(a, b):
    if len(a) < 3:
        return None
    ra, rb = ranks(a), ranks(b)
    ma, mb = sum(ra) / len(ra), sum(rb) / len(rb)
    sa = sum((x - ma) ** 2 for x in ra)
    sb = sum((y - mb) ** 2 for y in rb)
    if sa == 0 or sb == 0:
        return None
    return sum((x - ma) * (y - mb) for x, y in zip(ra, rb)) / (sa * sb) ** 0.5


def main():
    ap = argparse.ArgumentParser()
    for a in ("res", "code", "prereg", "result", "tsv"):
        ap.add_argument("--" + a, required=True)
    a = ap.parse_args()
    S6 = os.path.join(a.res, S)
    D = os.path.join(S6, "d49")
    paths = dict((k, os.path.join(a.code, "phase2.2", f)) for k, f in (
        ("d40", "cs6_d40.py"), ("d33b", "cs6_d33b.py"), ("p48", "cs6_d48_prepare.py"), ("p49", "cs6_d49_prepare.py")))
    for k, p in paths.items():
        if not os.path.isfile(p):
            stop("the committed script %s is missing." % p)
    for f in (a.prereg, os.path.join(D, "RUNS.tsv"), os.path.join(S6, "d48.tsv"), os.path.join(S6, D48PRE), os.path.join(S6, D48RES)):
        if not os.path.isfile(f):
            stop("%s is missing." % f)
    for f in (a.result, a.tsv):
        if os.path.exists(f):
            stop("%s exists already; this reader overwrites nothing." % f)
    d40, d33b = load(paths["d40"], "d49_cs6_d40"), load(paths["d33b"], "d49_cs6_d33b")
    p48, p49 = load(paths["p48"], "d49_cs6_d48_prepare"), load(paths["p49"], "d49_cs6_d49_prepare")
    gates = []

    def git(*args):
        try:
            q = subprocess.run(["git", "-C", a.res] + list(args), stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
            return q.returncode, q.stdout.strip()
        except OSError:
            return 127, ""

    # ------------------------------------------------------------ G49.3: the committed files this leans on
    guard = ["%s/%s" % (S, g) for g in GUARDED] + ["%s/%s" % (S, f) for f in ("d48.tsv", D48PRE, D48RES)]
    rc, st = git("status", "--porcelain", "--untracked-files=all", "--", *guard)
    if rc != 0:
        gates.append("G49.3 could not be checked: git did not run")
    elif st:
        gates.append("G49.3 a file under %s, or a table or document of D48, is not the committed one" % ", ".join(GUARDED))
    rc, pcm = git("log", "--diff-filter=A", "--format=%H", "-1", "--", "%s/%s" % (S, PREREG))
    if rc != 0 or not pcm:
        gates.append("G49.3 the commit that added %s could not be found" % PREREG)
        pcm = ""
    else:
        rc1, b1 = git("rev-parse", "%s:%s/%s" % (pcm, S, PREREG))
        rc2, b2 = git("hash-object", a.prereg)
        if rc1 != 0 or rc2 != 0 or b1 != b2:
            gates.append("G49.3 the pre-registration given is not the committed %s" % PREREG)

    # ------------------------------------------------------------ G49.2: the sources
    try:
        F, runs = p49.build(S6)
    except SystemExit:
        stop("cs6_d49_prepare.py could not form the files of d49 again from the committed files of D48.")
    nsrc = 0
    psrc = os.path.join(D, "SOURCES.sha256")
    if not os.path.isfile(psrc) or open(psrc).read() != F["SOURCES.sha256"]:
        gates.append("G49.2 d49/SOURCES.sha256 is not the list formed again from the committed files of D48")
    for lst, base in ((psrc, S6), (os.path.join(S6, "d48", "SOURCES.sha256"), S6)):
        if not os.path.isfile(lst):
            gates.append("G49.2 %s is missing" % os.path.relpath(lst, S6))
            continue
        for l in open(lst).read().splitlines():
            if l.strip():
                h, f = l.split(None, 1)
                nsrc += 1
                p = os.path.join(base, f)
                if not os.path.isfile(p) or sha(p) != h:
                    gates.append("G49.2 %s is not the file listed in %s" % (f, os.path.relpath(lst, S6)))

    # ------------------------------------------------------------ G49.1: the inputs
    nsum = 0
    pic = os.path.join(D, "INPUT_CHECKSUMS.sha256")
    if not os.path.isfile(pic):
        gates.append("G49.1 d49/INPUT_CHECKSUMS.sha256 is missing")
    else:
        if open(pic).read() != F["INPUT_CHECKSUMS.sha256"]:
            gates.append("G49.1 d49/INPUT_CHECKSUMS.sha256 is not the list formed again here")
        for l in open(pic).read().splitlines():
            if l.strip():
                h, f = l.split(None, 1)
                nsum += 1
                p = os.path.join(D, f)
                if not os.path.isfile(p) or sha(p) != h:
                    gates.append("G49.1 %s is not the file listed in INPUT_CHECKSUMS.sha256" % f)
    for k in sorted(F):
        p = os.path.join(D, k)
        if not os.path.isfile(p) or open(p).read() != F[k]:
            gates.append("G49.1 d49/%s is not the file formed again from the committed files of D48" % k)
    for r in runs:
        for fn in ("substrate.xyz",) + (("pc.pc",) if r["variant"] == "Z" else ()):
            if os.path.exists(os.path.join(D, r["run"], fn)):
                gates.append("G49.1 %s holds a file %s that its input does not read" % (r["run"], fn))

    # ------------------------------------------------------------ G49.4: the quoted sentences
    def has(path, sent):
        return " ".join(sent.split()) in " ".join(open(path, errors="replace").read().split())
    for doc, q, what in ((D48PRE, Q_A, "D48's sentence that S' predicts nothing"), (D48PRE, Q_B, "D48's sentence on a calibration not measured"),
                         (D48RES, Q_C, "D48's count of the calibration")):
        if not has(os.path.join(S6, doc), q):
            gates.append("G49.4 %s is not in %s" % (what, doc))
    for q in (Q_A, Q_B):
        if not has(a.prereg, q):
            gates.append("G49.4 the words \"%s\" are not in the pre-registration given" % (q[:60] + " ..."))

    # ------------------------------------------------------------ the runs of D49, and G49.5, G49.6
    RUN = {}
    for r in runs:
        d = os.path.join(D, r["run"])
        x = d40.read(d, "pc" if r["variant"] == "N" else "plain", 1, d33b)
        x["status"] = "measured"
        if not x.get("exists"):
            x["status"] = "NOT MEASURED"
        elif x["E"] is None:
            x["status"] = "DID NOT CONVERGE" if "SCF NOT CONVERGED" in x["text"] else "NOT MEASURED"
        x.update(r)
        RUN[r["run"]] = x
        if x.get("exists"):
            e = echoed(x["text"])
            want = open(os.path.join(d, "job.inp")).read().split("\n") if os.path.isfile(os.path.join(d, "job.inp")) else None
            while want and not want[-1].strip():
                want.pop()
            if e is None:
                if x.get("normal"):
                    gates.append("G49.5 %s/job.out ends normally and does not print the input it ran" % r["run"])
            elif want is not None and e != want:
                gates.append("G49.5 the input printed at the top of %s/job.out is not the committed %s/job.inp" % (r["run"], r["run"]))
            if x["status"] == "measured":
                m = re.findall(r"Number of basis functions\s*\.+\s*(\d+)", x["text"])
                if not m or any(int(v) != NBF for v in m):
                    gates.append("G49.6 %s reports %s basis functions, not %d" % (r["run"], ", ".join(sorted(set(m))) or "no count of", NBF))

    # ------------------------------------------------------------ G49.7: S formed again is the S of d48.tsv
    T48 = table(os.path.join(S6, "d48.tsv"))
    need = ("run", "condition", "part", "state", "S_eV", "sum_24_minus_control_e", "what_the_run_did", "how_it_took", "S_predicts", "prediction")
    sel = None
    if T48 and any(k not in T48[0] for k in need):
        gates.append("G49.7 d48.tsv does not hold the columns this reader needs")
    else:
        try:
            sel = p48.select(S6)
        except SystemExit:
            gates.append("G49.7 cs6_d48_prepare.py could not form D48's conditions again from the committed files")
    SP = {}
    if sel:
        M = sel["M"]
        cmap = dict((x["name"], x) for x in sel["conds"])
        ct = dict((st, p48.energy(os.path.join(S6, "%s_%s" % (p48.CTL, st), "job.out"))) for st in p48.STATES)
        bm = dict((st, p48.energy(os.path.join(S6, "%s_%s" % (p48.BARE, st), "job.out"))) for st in p48.STATES)

        def s_of(x, st, bind=None):
            mo, di = M[("mo", st)], M[("di", st)]
            vde = (bm[st] - ct[st]) * EV + sum(q * (mo[p] - di[p]) * EV for (q, xyz), p in zip(x["charges"], x["rows"]))
            best = None
            for i, ((q, xyz), p) in enumerate(zip(x["charges"], x["rows"])):
                if q <= 0:
                    continue
                other = sum(qj * BOHR / p48.dist(xyz, xj) for j, (qj, xj) in enumerate(x["charges"]) if j != i)
                last = RYD * q * q if bind is None else bind(q)
                if last is None:
                    return None, None
                s = vde - EV * (mo[p] + other) - last
                if best is None or s < best[0]:
                    best = (s, q)
            return best
        if len(T48) != 2 * len(cmap):
            gates.append("G49.7 d48.tsv holds %d runs and D48's conditions give %d" % (len(T48), 2 * len(cmap)))
        for row in T48:
            x = cmap.get(row["condition"])
            if x is None:
                gates.append("G49.7 d48.tsv holds a condition, %s, that D48's conditions do not" % row["condition"])
                continue
            s, q = s_of(x, row["state"])
            if abs(s - float(row["S_eV"])) > TOL_S:
                gates.append("G49.7 S of %s formed here is %+.4f eV and d48.tsv holds %s" % (row["run"], s, row["S_eV"]))
            SP[row["run"]] = dict(x=x, S=s, q=q)

    if gates:
        sys.stdout.write("\n".join("FAILED  " + g for g in gates) + "\n")
        stop("%d gate(s) failed. Nothing was written." % len(gates))

    # ------------------------------------------------------------ the calibration
    one = p49.one_runs(S6)
    sizes = []
    for run, lv, size, q in one:
        if size not in sizes:
            sizes.append(size)
    CAL = {}
    for run, lv, size, q in one:
        n, z = RUN["%s_N" % run], RUN["%s_Z" % run]
        t48 = open(os.path.join(S6, "d48", run, "job.out"), errors="replace").read()
        e48 = re.findall(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", t48)
        ok48 = bool(re.search(r"SCF CONVERGED AFTER\s+\d+\s+CYCLES", t48)) and "SCF NOT CONVERGED" not in t48 and len(e48) == 1
        c = dict(run=run, lv=lv, size=size, q=q, eN=n["E"] if n["status"] == "measured" else None, eZ=z["E"] if z["status"] == "measured" else None,
                 e48=float(e48[0]) if ok48 else None, exact=RYD * q * q)
        c["src"] = "N" if c["eN"] is not None else ("Z" if c["eZ"] is not None else "")
        c["E"] = c["eN"] if c["eN"] is not None else c["eZ"]
        c["B"] = None if c["E"] is None else -c["E"] * EV
        CAL[(lv, size)] = c

    def bind_cam(q):
        c = CAL.get(("CAM", "%.2f" % q))
        return None if c is None else c["B"]
    for row in T48:
        y = SP[row["run"]]
        sp, q = s_of(y["x"], row["state"], bind_cam)
        y.update(Sp=sp, left=float(row["sum_24_minus_control_e"]) if row["sum_24_minus_control_e"] else None, did=row["what_the_run_did"],
                 how=row["how_it_took"], cond=row["condition"], st=row["state"], part=row["part"], run=row["run"], pred=row["S_predicts"], verdict=row["prediction"])

    # ------------------------------------------------------------ the result
    rep = []
    say = rep.append

    def fs(v, fmt="%.6f", none="-"):
        return none if v is None else fmt % v

    def wrapsay(text, lead="   ", width=96, hang=""):
        for i, l in enumerate(textwrap.wrap(text, width - len(lead))):
            say(lead + (hang if i else "") + l)
    rc, head = git("rev-parse", "HEAD")
    say("STAGE 6 - DIAGNOSTIC 49 - RESULT")
    say("The calibration of DIAGNOSTIC 48, run again so that ORCA ends normally")
    say("=" * 96)
    say("")
    say("Read by cs6_d49.py under %s, which was committed before any job was submitted." % PREREG)
    say("   results repository at %s; the pre-registration was added by %s" % (head[:12] if rc == 0 else "?", pcm[:12]))
    say("   reader               cs6_d49.py           sha256 %s" % sha(os.path.abspath(__file__)))
    say("   inputs formed again by cs6_d49_prepare.py  sha256 %s" % sha(paths["p49"]))
    say("   S formed again by    cs6_d48_prepare.py   sha256 %s" % sha(paths["p48"]))
    say("   outputs read by      cs6_d40.py           sha256 %s" % sha(paths["d40"]))
    say("")
    say("One electron, the def2-TZVPD functions of hydrogen, a bare charge q. No substrate. D49 adopts nothing")
    say("and edits nothing of D48: D48's own reading of its calibration stays NOT MEASURED. S' PREDICTS")
    say("NOTHING: D48's result was known when D49 was written.")
    say("")
    say("THE GATES, as they held:")
    say("   G49.1  all %d files of d49/INPUT_CHECKSUMS.sha256 are the ones listed, and every file of d49 is the" % nsum)
    say("          one formed again here from the committed files of D48.")
    say("   G49.2  the %d files of d49/SOURCES.sha256 and of d48/SOURCES.sha256 are the ones listed." % nsrc)
    say("   G49.3  no file under %s, and no table or document of D48," % ", ".join(GUARDED))
    say("          differs from the committed one; the pre-registration given is the committed one.")
    say("   G49.4  the sentences quoted here are in the documents they are quoted from.")
    say("   G49.5  every output that exists prints, as the input it ran, the committed input beside it.")
    say("   G49.6  every measured run reports %d basis functions." % NBF)
    say("   G49.7  S of each of D48's %d runs with a charge, formed again here, is the S that d48.tsv holds." % len(T48))
    say("")
    # ---- R49.1
    say("=" * 96)
    say("R49.1  THE RUNS")
    say("=" * 96)
    say("")
    say("   run                level  variant  q, e   status             energy, Eh")
    tail = []
    for r in runs:
        x = RUN[r["run"]]
        say("   %-17s  %-5s  %-7s  %s   %-17s  %s" % (r["run"], r["level"], r["variant"], r["size"], x["status"], fs(x["E"] if x["status"] == "measured" else None, "%.10f")))
        if x["status"] != "measured":
            tail.append("%s: %s: %s" % (r["run"], x["status"], "; ".join(x.get("why") or ["no readable output"])))
    say("")
    for w in tail:
        wrapsay("NOTE  " + w, "   ", 96, "      ")
    nN = sum(1 for r in runs if r["variant"] == "N" and RUN[r["run"]]["status"] == "measured")
    nZ = sum(1 for r in runs if r["variant"] == "Z" and RUN[r["run"]]["status"] == "measured")
    say("   Measured: %d of %d of variant N, %d of %d of variant Z." % (nN, len(one), nZ, len(one)))
    say("")
    # ---- R49.2
    say("=" * 96)
    say("R49.2  THE CALIBRATION: WHAT THE FUNCTIONS BIND AT A CHARGE q")
    say("=" * 96)
    say("")
    say("B is minus the energy of the one electron, eV. Exact, in a complete basis: 13.6057 q^2. HF is exact for")
    say("one electron but for the functions; CAM is the functional of D48's runs without its dispersion term.")
    say("B is taken from variant N where it is measured, otherwise from variant Z; 'from' says which.")
    say("")
    say("   q, e   13.6057 q^2   B, HF      from   B, CAM     from   exact minus HF   CAM minus HF   exact minus CAM (what S' adds to S)")
    for size in sizes:
        h, c = CAL[("HF", size)], CAL[("CAM", size)]
        say("   %s   %-11s   %-9s  %-4s   %-9s  %-4s   %-14s   %-12s   %s" % (
            size, "%.4f" % h["exact"], fs(h["B"], "%+.4f"), h["src"] or "-", fs(c["B"], "%+.4f"), c["src"] or "-",
            fs(None if h["B"] is None else h["exact"] - h["B"], "%+.4f"), fs(None if (h["B"] is None or c["B"] is None) else c["B"] - h["B"], "%+.4f"),
            fs(None if c["B"] is None else c["exact"] - c["B"], "%+.4f")))
    say("")
    ncal = sum(1 for c in CAL.values() if c["B"] is not None)
    say("   Measured: %d of %d (from variant N %d, from variant Z %d)." % (ncal, len(CAL), sum(1 for c in CAL.values() if c["src"] == "N"), sum(1 for c in CAL.values() if c["src"] == "Z")))
    neg = [c for c in CAL.values() if c["B"] is not None and c["B"] <= 0]
    for c in neg:
        say("   B at or below zero: %s at %s e. These functions do not hold an electron at a charge of that size." % (c["lv"], c["size"]))
    say("")
    say("THE TWO VARIANTS AND D48'S OWN OUTPUTS, energies in Eh. Variant Z has the charge at the centre of the")
    say("functions, variant N 0.01 A from it. D48's outputs did not end normally: the last column is the final")
    say("energy they hold, set here for comparison only, and it is not a measured value of D48.")
    say("")
    say("   level  q, e   variant N        variant Z        D48's output     Z minus N, eV   N minus D48's, Eh")
    dzn, dn48 = [], []
    for lv in p49.LEVELS:
        for size in sizes:
            c = CAL[(lv, size)]
            a1 = None if (c["eN"] is None or c["eZ"] is None) else (c["eZ"] - c["eN"]) * EV
            a2 = None if (c["eN"] is None or c["e48"] is None) else c["eN"] - c["e48"]
            if a1 is not None:
                dzn.append((abs(a1), a1, lv, size))
            if a2 is not None:
                dn48.append((abs(a2), a2, lv, size))
            say("   %-5s  %s   %-15s  %-15s  %-15s  %-14s  %s" % (lv, size, fs(c["eN"], "%.10f"), fs(c["eZ"], "%.10f"), fs(c["e48"], "%.10f"), fs(a1, "%+.5f"), fs(a2, "%+.2e")))
    say("")
    if dzn:
        m = max(dzn)
        say("   Z minus N: formed at %d of %d; the largest is %+.5f eV (%s, %s e)." % (len(dzn), len(CAL), m[1], m[2], m[3]))
        if m[0] > 0.02:
            say("   MORE THAN 0.02 eV at %d size(s): the two variants are not one calibration there. Where B is taken" % sum(1 for v in dzn if v[0] > 0.02))
            say("   from variant Z it is the charge at the centre of the functions, and S' from it is that.")
    else:
        say("   Z minus N: not formed; no size has both variants measured.")
    if dn48:
        m = max(dn48)
        say("   N minus the final energy in D48's output: formed at %d of %d; the largest is %+.2e Eh (%s, %s e)." % (len(dn48), len(CAL), m[1], m[2], m[3]))
        if m[0] > 1e-6:
            say("   MORE THAN 0.000001 hartree at %d size(s): variant N does not give back the final energy of D48's" % sum(1 for v in dn48 if v[0] > 1e-6))
            say("   output there, though the inputs differ by the dipole switch alone.")
    else:
        say("   N minus the final energy in D48's output: not formed.")
    say("")
    # ---- R49.3
    say("=" * 96)
    say("R49.3  S' BESIDE S, FOR EACH OF D48'S RUNS WITH A CHARGE")
    say("=" * 96)
    say("")
    say("S and what each run did are D48's (d48.tsv). S' is S with 13.6057 q^2 replaced by B(CAM) at the size of")
    say("the charge S names. 'left' is the electron charge that left the 24 atoms against the control, e.")
    say("")
    say("   run                                    q, e   S, eV    S', eV   left     the run did")
    rows = [SP[row["run"]] for row in T48]
    for y in rows:
        say("   %-37s  %.2f   %+.2f    %-6s   %-7s  %s%s" % (y["run"], y["q"], y["S"], fs(y["Sp"], "%+.2f"), fs(y["left"], "%+.3f"), y["did"], ", %s" % y["how"] if y["how"] == "A LEG ONLY" else ""))
    say("")
    nsp = sum(1 for y in rows if y["Sp"] is not None)
    say("   S' is formed for %d of %d runs." % (nsp, len(rows)))
    say("")
    # ---- R49.4
    say("=" * 96)
    say("R49.4  WHETHER S' SEPARATES WHAT D48'S RUNS DID, BESIDE S")
    say("=" * 96)
    say("")
    say("As D48 defines it: it separates if the largest value at which a run TOOK ELECTRONS is below the smallest")
    say("at which a run DID NOT TAKE. Over the runs for which S' is formed, with the runs A LEG ONLY set aside.")
    say("")
    use = [y for y in rows if y["Sp"] is not None and y["how"] != "A LEG ONLY"]
    nleg = sum(1 for y in rows if y["how"] == "A LEG ONLY")
    sepw = {}
    for key, lab in (("S", "S"), ("Sp", "S'")):
        yes = [(y[key], y["run"]) for y in use if y["did"] == "TOOK ELECTRONS"]
        no = [(y[key], y["run"]) for y in use if y["did"] == "DID NOT TAKE"]
        if not yes or not no:
            sepw[key] = "NOT FORMED"
            say("   %s: the separation is NOT FORMED: %s." % (lab, "S' is formed for no run" if not use else "the runs did not do both things"))
            say("")
            continue
        hi, lo = max(yes), min(no)
        say("   %s: the largest at which a run took electrons is %+.2f eV (%s);" % (lab, hi[0], hi[1]))
        say("   %s  the smallest at which a run did not take is %+.2f eV (%s)." % (" " * len(lab), lo[0], lo[1]))
        if hi[0] < lo[0]:
            sepw[key] = "SEPARATES"
            say("   %s SEPARATES what the runs did: any threshold between %+.2f and %+.2f eV gives every one of them." % (lab, hi[0], lo[0]))
        else:
            sepw[key] = "DOES NOT SEPARATE"
            inside = sorted([(v, n, "TOOK") for v, n in yes if v >= lo[0]] + [(v, n, "did not") for v, n in no if v <= hi[0]])
            say("   %s DOES NOT SEPARATE what the runs did. Between %+.2f and %+.2f eV lie:" % (lab, lo[0], hi[0]))
            for v, n, w in inside:
                say("      %+.2f  %-8s %s" % (v, w, n))
        say("")
    say("   Runs in this reading: %d of %d (S' not formed for %d; A LEG ONLY, set aside: %d)." % (len(use), len(rows), len(rows) - nsp, nleg))
    say("")
    # ---- R49.5
    say("=" * 96)
    say("R49.5  HOW WELL S AND S' ORDER THE CHARGE THAT LEFT THE ATOMS")
    say("=" * 96)
    say("")
    say("Spearman's rank correlation between the estimate and the electron charge that left the 24 atoms, over")
    say("the runs for which S' is formed and the charge is read. -1 would be a perfect order: the lower the")
    say("estimate, the more charge left. On D48's committed table; no threshold enters it.")
    say("")
    RHO = {}
    for lab, pick in (("one positive charge (part A)", lambda y: y["part"] == "A"), ("every run with a charge", lambda y: True)):
        u = [y for y in rows if pick(y) and y["Sp"] is not None and y["left"] is not None]
        rs, rp = spearman([y["S"] for y in u], [y["left"] for y in u]), spearman([y["Sp"] for y in u], [y["left"] for y in u])
        RHO[lab] = (len(u), rs, rp)
        say("   %-29s %2d runs:   S %s    S' %s" % (lab + ",", len(u), fs(rs, "%+.4f"), fs(rp, "%+.4f")))
    say("")
    # ---- R49.6
    say("=" * 96)
    say("R49.6  WHAT THE TEST OF B_n11_q0.25 WITH FUNCTIONS WAS WORTH")
    say("=" * 96)
    say("")
    c25 = CAL.get(("CAM", "0.25"))
    b11 = [y for y in rows if y["cond"] == "B_n11_q0.25"]
    if c25 is None or c25["B"] is None:
        worth = "NOT KNOWN"
        say("   B(CAM) at 0.25 e is not measured: NOT KNOWN.")
    else:
        worth = "NOT A TEST OF CAPTURE" if c25["B"] <= 0 else ("A WEAK TEST" if c25["B"] < 0.5 * c25["exact"] else "A TEST AS IT STANDS")
        say("   B(CAM) at 0.25 e is %+.4f eV; in a complete basis it would be %.4f eV." % (c25["B"], c25["exact"]))
        say("   By the pre-registration: %s." % worth)
    for y in b11:
        say("   %s: S %+.2f eV, S' %s eV; the run %s (%s e left the atoms)." % (y["run"], y["S"], fs(y["Sp"], "%+.2f"), y["did"], fs(y["left"], "%+.3f")))
    say("")
    # ---- cannot say
    say("=" * 96)
    say("WHAT THIS CANNOT SAY (section 6 of the pre-registration, in short; its own words govern)")
    say("=" * 96)
    say("")
    for w in ("S' is formed after D48's result was known. It predicts nothing and D49 tests nothing.",
              "B is what these nine functions bind. It says what D48's runs with functions could show; it does not say what a complete basis would.",
              "Variant Z has the charge at the centre of the functions and not 0.01 A from it; where B is taken from Z the table says so.",
              "The final energies in D48's outputs are set beside for comparison. D48's reading of its calibration stays NOT MEASURED.",
              "One electron and a bare charge: no substrate, no other charge, nothing relaxed.",
              "D49 adopts nothing. A threshold on S or on S' is for a pre-registration of its own, on conditions not yet run."):
        wrapsay(w, "   ", 96, "")
    say("")
    line = ("Calibration measured: %d of %d (variant N %d of %d runs, variant Z %d of %d). S' formed for %d of %d runs. S' %s what D48's runs did; S %s. "
            "Rank correlation with the charge that left, single charges: S %s, S' %s. B_n11_q0.25 with functions: %s." % (
                ncal, len(CAL), nN, len(one), nZ, len(one), nsp, len(rows),
                {"SEPARATES": "separates", "DOES NOT SEPARATE": "does not separate", "NOT FORMED": "is not held against"}[sepw["Sp"]],
                {"SEPARATES": "separates them", "DOES NOT SEPARATE": "does not", "NOT FORMED": "is not held against them"}[sepw["S"]],
                fs(RHO["one positive charge (part A)"][1], "%+.3f"), fs(RHO["one positive charge (part A)"][2], "%+.3f"), worth))
    say(line)
    # ---- tables
    T = ["run\tlevel\tvariant\tsize_e\tstatus\tE_Eh\tminus_E_eV\twhy"]
    for r in runs:
        x = RUN[r["run"]]
        ok = x["status"] == "measured"
        T.append("\t".join([r["run"], r["level"], r["variant"], r["size"], x["status"], fs(x["E"] if ok else None, "%.10f", ""), fs(-x["E"] * EV if ok else None, "%+.6f", ""),
                            "" if ok else "; ".join(x.get("why") or ["no readable output"])]))
    C = ["level\tsize_e\texact_eV\tB_eV\tfrom\tE_N_Eh\tE_Z_Eh\tE_D48_output_Eh\tZ_minus_N_eV\tN_minus_D48_output_Eh\texact_minus_B_eV"]
    for lv in p49.LEVELS:
        for size in sizes:
            c = CAL[(lv, size)]
            C.append("\t".join([lv, size, "%.6f" % c["exact"], fs(c["B"], "%+.6f", ""), c["src"], fs(c["eN"], "%.10f", ""), fs(c["eZ"], "%.10f", ""), fs(c["e48"], "%.10f", ""),
                                fs(None if (c["eN"] is None or c["eZ"] is None) else (c["eZ"] - c["eN"]) * EV, "%+.6f", ""),
                                fs(None if (c["eN"] is None or c["e48"] is None) else c["eN"] - c["e48"], "%+.3e", ""), fs(None if c["B"] is None else c["exact"] - c["B"], "%+.6f", "")]))
    P = ["run\tcondition\tpart\tstate\tsize_of_the_charge_S_names_e\tS_eV\tS_prime_eV\tsum_24_minus_control_e\twhat_the_run_did\thow_it_took"]
    for y in rows:
        P.append("\t".join([y["run"], y["cond"], y["part"], y["st"], "%.2f" % y["q"], "%+.4f" % y["S"], fs(y["Sp"], "%+.4f", ""), fs(y["left"], "%+.4f", ""), y["did"], y["how"]]))
    rep[:] = [l.rstrip() for l in rep]
    base = os.path.splitext(a.tsv)[0]
    open(a.tsv, "w").write("\n".join(T) + "\n")
    open(base + "_calibration.tsv", "w").write("\n".join(C) + "\n")
    open(base + "_sprime.tsv", "w").write("\n".join(P) + "\n")
    open(a.result, "w").write("\n".join(rep) + "\n")
    sys.stdout.write("\n".join(rep) + "\n")


if __name__ == "__main__":
    main()
