#!/usr/bin/env python3
"""cs6_d45.py - DIAGNOSTIC 45's reader: the approximations of the coupled-cluster benchmark, measured.

Fourteen new single points on the bare substrate of frame 42436 (d45/), beside the four committed
ones the benchmark is made of: d38/N and d41/T at def2-SVPD, d36a2/B_R and B_TS at def2-TZVPD.

Every coupled-cluster output, new or committed, is read by D38's committed reader, cs6_d38.py
read_out, imported by path, so a committed leg is read here exactly as D38 and D41 read it. The
inputs are rebuilt by cs6_d45_prepare.py, imported by path, and compared with the files that ran.
What ORCA prints for a pure RI-MP2 run, for DLPNO-CCSD(T1) and for a TCutPNO other than 1e-7 has no
committed precedent: the patterns below were written from the committed outputs of the
non-iterative runs and of the double hybrid d37/PWPB95_SV_R, and a run whose lines they do not
find is NOT READ, with the line named.

Nothing is typed in except committed reference values, each named with its source and checked
against the committed file it comes from before it is used, and the constants of the two
extrapolation formulas, which the pre-registration names with their sources.

A failure of G45.1 to G45.4 stops the reader and nothing is written. G45.5 is per run.

usage: cs6_d45.py --res <results-repo> --code <code-repo> --prereg <file> --result <file> --tsv <file>
"""
import argparse, hashlib, importlib.util, math, os, re, subprocess, sys

H2K = 627.5094740631
S = "phase2.2/cs_stage6"
NEL, NCORE, NOCC = 118, 16, 59              # electrons; frozen core orbitals (32 el, d38/N); occupied orbitals
TOL_ID = 2e-9                               # identities inside one output, Eh (D41, G41.5)
TOL_REP = 1e-6                              # the same input reproduces a committed energy, Eh (D44, section 2)
ALPHA34, BETA34 = 7.88, 2.97                # ORCA 6.0 manual, automatic extrapolation, def2-TZVPP/def2-QZVPP
CPS_F = 1.5                                 # Altun, Neese, Bistoni, J. Chem. Theory Comput. 2020, 16, 6142
B_SV, B_TZ, SHIFT = "22.8968", "23.8603", "+0.9635"       # D41's result, D36's result, D41's result
CAM_DIST = "+1.8400"                        # AMENDMENT 7, A7.1: CAM-B3LYP-D3BJ at def2-TZVPD against D36
GUARDED = ("d33a", "d36a2", "d38", "d41")
D36R, D41P, D41R, A7 = ("STAGE6_DIAGNOSTIC36_RESULT_20261114.txt", "STAGE6_DIAGNOSTIC41_20261115.txt",
                        "STAGE6_DIAGNOSTIC41_RESULT_20261115.txt", "STAGE6_AMENDMENT7_20261117.txt")
QUOTES = ((D36R, "THE BARRIER IS %s kcal/mol." % B_TZ),
          (D36R, "TightPNO was requested; the residual PNO error is not measured here."),
          (D36R, "23.8603 is a reference for comparison at this basis and must not be quoted as the true in vacuo barrier."),
          (D41R, "THE BARRIER AT def2-SVPD IS %s kcal/mol." % B_SV),
          (D41R, "THE CORRELATED BARRIER SHIFTS BY %s kcal/mol from def2-SVPD to def2-TZVPD." % SHIFT),
          (D41P, "DLPNO is an approximation; TightPNO is requested and its residual error is not measured."),
          (D41P, "Two basis sets do not make an extrapolation, and none is attempted."),
          (A7, "WHICH BENCHMARK IS NEARER THE TRUTH."),
          (A7, "Nothing in the record measures it. Agreement of three functionals is not that measurement."))
CTL = {("def2-SVPD", "R"): "d38/N", ("def2-SVPD", "TS"): "d41/T",
       ("def2-TZVPD", "R"): "d36a2/B_R", ("def2-TZVPD", "TS"): "d36a2/B_TS"}
KEEP = ("TCutPairs", "TCutMKN", "TCutDO")   # TightPNO's other three thresholds, which no run here may change
COMP = (("Hartree-Fock", "e0"), ("CCSD correlation", "corr"), ("triples", "tcorr"), ("CCSD(T)", "ccsdt"))


def stop(m):
    sys.stderr.write("STOP  %s\n" % m)
    raise SystemExit(1)


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 16), b""):
            h.update(b)
    return h.hexdigest()


def flat(s):
    return " ".join(s.split())


def load(path, name):
    sp = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(sp)
    sp.loader.exec_module(m)
    return m


def echoed(text):
    """The input ORCA prints at the top of its output, as lines; None if there is none. As cs6_d44.py reads it."""
    m = re.search(r"INPUT FILE\n=+\nNAME = [^\n]*\n((?:\|\s*\d+>[^\n]*\n)+)", text)
    if not m:
        return None
    ls = [re.sub(r"^\|\s*\d+> ?", "", l) for l in m.group(1).split("\n") if l]
    while ls and (not ls[-1].strip() or "****END OF INPUT****" in ls[-1]):
        ls.pop()
    return ls


def settings(text):
    """ORCA's DLPNO SETTINGS block, as the lines it prints, in its order; [] if there is none."""
    m = re.search(r"DLPNO SETTINGS[^\n]*\n-+\n\n(.*?)\n\n", text, re.S)
    return [flat(l) for l in m.group(1).split("\n")] if m else []


def setting(block, name):
    for l in block:
        k = re.match(r"^%s\s*:\s*(\S+)$" % re.escape(name), l)
        if k:
            return k.group(1)
    return None


def triples_times(text):
    a = re.findall(r"Total time for T0\s+\.\.\.\s+(\d+\.\d+) sec", text)
    b = re.findall(r"Total time for Iterative \(T\)\s+\.\.\.\s+(\d+\.\d+) sec", text)
    return (float(a[-1]) if a else None, float(b[-1]) if b else None)


def nbf(text):
    m = re.search(r"Basis Dimension\s+Dim\s*\.+\s*(\d+)", text)
    return int(m.group(1)) if m else None


def one(vals, label, bad):
    """the one value a label prints; a label printed more than once must print one value"""
    if not vals:
        return None
    if len(set(vals)) > 1:
        bad.append("%s is printed with %d different values" % (label, len(set(vals))))
    return float(vals[-1])


def mp2_of(text):
    bad = []
    scf = one(re.findall(r"^Total Energy\s+:\s+(-?\d+\.\d+) Eh", text, re.M), "the SCF Total Energy", bad)
    cor = one(re.findall(r"RI-MP2 CORRELATION ENERGY:\s+(-?\d+\.\d+) Eh", text), "RI-MP2 CORRELATION ENERGY", bad)
    tot = one(re.findall(r"MP2 TOTAL ENERGY:\s+(-?\d+\.\d+) Eh", text), "MP2 TOTAL ENERGY", bad)
    rng = re.search(r"Orbital Range Operator\s+0:\s+(\d+)-\s*(\d+) to\s+(\d+)-\s*(\d+)", text)
    return dict(scf=scf, cor=cor, tot=tot, rng=tuple(int(x) for x in rng.groups()) if rng else None, bad=bad)


def orca_words(t):
    L = t.splitlines()
    err = [l.strip() for l in L if re.search(r"error|aborting|not converged|No space left|cannot", l, re.I)]
    if err:
        return err[-1][:160]
    return "no error line found; the output ends: " + (L[-1].strip()[:120] if L else "(empty)")


def job_facts(D):
    """R45.8 and the per-run lines, from every log of the job script that is there, the earliest first."""
    logs = [f for f in os.listdir(D) if re.match(r"^run\.pbs\.out(\.part\d+)?$", f)]
    logs.sort(key=lambda f: (0, int(f.rsplit("part", 1)[1])) if "part" in f else (1, 0))
    places, per = [], {}
    for f in logs:
        t = open(os.path.join(D, f), errors="replace").read()
        h, j = re.search(r"host=(\S+)", t), re.search(r"jobid=(\S+)", t)
        m = re.search(r"R45\.8 scratch=(\S+) avail=(\d+)G margin=(\S+)", t)
        fin = re.search(r"^finished: (\d+) of (\d+)", t, re.M)
        places.append(dict(log=f, host=h.group(1) if h else None, jobid=j.group(1) if j else None,
                           scratch=m.group(1) if m else None, avail=int(m.group(2)) if m else None,
                           margin=m.group(3) if m else None, finished=fin.group(0) if fin else None,
                           refused=re.findall(r"^   candidate (\S+)\s+(.*)$", t, re.M)))
        for r, rc, n, pk in re.findall(r"^(\S+) orca_exit=(\d+) normal=(\d+) peak_mb=(\d+)", t, re.M):
            per[r] = dict(rc=int(rc), peak=int(pk), log=f)
    return places, per


def f4(x, plus=True):
    if x is None:
        return "NOT FORMED"
    return ("%+.4f" if plus else "%.4f") % x


def main():
    ap = argparse.ArgumentParser()
    for a in ("res", "code", "prereg", "result", "tsv"):
        ap.add_argument("--" + a, required=True)
    a = ap.parse_args()
    S6 = os.path.join(a.res, S)
    D = os.path.join(S6, "d45")
    p38 = os.path.join(a.code, "phase2.2", "cs6_d38.py")
    pprep = os.path.join(a.code, "phase2.2", "cs6_d45_prepare.py")
    for p in (p38, pprep, a.prereg):
        if not os.path.isfile(p):
            stop("%s is missing." % p)
    d38 = load(p38, "d45_cs6_d38")
    prep = load(pprep, "d45_cs6_d45_prepare")

    rep, gates = [], []
    say = rep.append
    say("STAGE 6, DIAGNOSTIC 45 - RESULT")
    say("The approximations of the coupled-cluster benchmark, measured: the triples, the PNO")
    say("threshold and the basis set. Bare substrate of frame 42436, no charges.")
    say("")
    say("Read by phase2.2/cs6_d45.py               sha256 %s" % sha(__file__))
    say("   coupled-cluster outputs by cs6_d38.py   sha256 %s" % sha(p38))
    say("   inputs rebuilt by cs6_d45_prepare.py    sha256 %s" % sha(pprep))
    say("   pre-registration %-22s sha256 %s" % (os.path.basename(a.prereg), sha(a.prereg)))
    say("")

    # ------------------------------------------------------------------ G45.1 the inputs
    try:
        prep.select(S6)
        runs = prep.runs_of()
    except prep.Refused as e:
        stop("G45.1/G45.2 the committed files are not what the inputs were built from: %s" % e)
    want = {}
    for run, kind, basis, st in runs:
        want["%s/job.inp" % run] = prep.build_input(kind, basis).encode()
        want["%s/substrate.xyz" % run] = open(os.path.join(S6, prep.GEOM[st]), "rb").read()
    want["RUNS.tsv"] = prep.runs_text(runs).encode()
    want["run.pbs"] = prep.pbs_text(runs, prep.HPC).encode()
    pc = os.path.join(D, "INPUT_CHECKSUMS.sha256")
    listed = {}
    if not os.path.isfile(pc):
        gates.append("G45.1 d45/INPUT_CHECKSUMS.sha256 is missing")
    else:
        for l in open(pc).read().splitlines():
            f = l.split(None, 1)
            if len(f) == 2:
                listed[f[1].strip()] = f[0]
        if sorted(listed) != sorted(list(want) + ["SOURCES.sha256"]):
            gates.append("G45.1 d45/INPUT_CHECKSUMS.sha256 does not list exactly the %d files of the diagnostic" % (len(want) + 1))
        for n, h in sorted(listed.items()):
            p = os.path.join(D, n)
            if not os.path.isfile(p) or sha(p) != h:
                gates.append("G45.1 d45/%s is missing or is not the file whose sha256 is listed" % n)
            elif n in want and open(p, "rb").read() != want[n]:
                gates.append("G45.1 d45/%s is not what cs6_d45_prepare.py writes from the committed files" % n)

    # ------------------------------------------------------------------ G45.2 the sources
    ps = os.path.join(D, "SOURCES.sha256")
    if not os.path.isfile(ps):
        gates.append("G45.2 d45/SOURCES.sha256 is missing")
    else:
        src = {}
        for l in open(ps).read().splitlines():
            f = l.split(None, 1)
            if len(f) == 2:
                src[f[1].strip()] = f[0]
        if sorted(src) != sorted(prep.sources()):
            gates.append("G45.2 d45/SOURCES.sha256 does not list exactly the files this diagnostic reads")
        for n, h in sorted(src.items()):
            p = os.path.join(S6, n)
            if not os.path.isfile(p) or sha(p) != h:
                gates.append("G45.2 %s is missing or is not the file whose sha256 d45/SOURCES.sha256 records" % n)
    try:
        q = subprocess.run(["git", "-C", a.res, "diff", "--quiet", "HEAD", "--"] + ["%s/%s" % (S, g) for g in GUARDED])
        if q.returncode != 0:
            gates.append("G45.2 a committed file under %s differs from HEAD" % ", ".join(GUARDED))
    except OSError:
        gates.append("G45.2 could not be checked: git did not run")

    # ------------------------------------------------------------------ G45.3 the controls
    C, CT = {}, {}
    for key, rel in CTL.items():
        r = d38.read_out(os.path.join(S6, rel))
        C[key] = r
        if r is None:
            gates.append("G45.3 %s/job.out is missing" % rel)
            continue
        t = open(r["path"], errors="replace").read()
        CT[key] = dict(block=settings(t), t0=triples_times(t)[0], it=triples_times(t)[1], nbf=nbf(t), echo=echoed(t))
        if not r["norm"] or r["notconv"] or r["nel"] != NEL:
            gates.append("G45.3 %s did not end normally with a converged SCF and %d electrons" % (rel, NEL))
        if any(r.get(k) is None for _, k in COMP):
            gates.append("G45.3 %s does not print every component of its energy" % rel)
            continue
        for b in r["inconsistent"]:
            gates.append("G45.3 %s: %s" % (rel, b))
        if abs(r["e0"] + r["corr"] - r["ccsd"]) > TOL_ID or abs(r["ccsd"] + r["tcorr"] - r["ccsdt"]) > TOL_ID:
            gates.append("G45.3 %s: its components do not add up to its E(CCSD(T)) within %.0e Eh" % (rel, TOL_ID))
        if setting(CT[key]["block"], "TCutPNO") != "1.000e-07":
            gates.append("G45.3 %s does not print TCutPNO: 1.000e-07" % rel)
        if CT[key]["it"] != 0.0 or not CT[key]["t0"]:
            gates.append("G45.3 %s does not print a time for T0 and 0.000 sec for the iterative (T)" % rel)
        if CT[key]["nbf"] is None:
            gates.append("G45.3 %s does not print its number of basis functions" % rel)
    bar = {}
    if not any(g.startswith("G45.3") for g in gates):
        for basis, lab in (("def2-SVPD", B_SV), ("def2-TZVPD", B_TZ)):
            bar[basis] = (C[(basis, "TS")]["ccsdt"] - C[(basis, "R")]["ccsdt"]) * H2K
            if "%.4f" % bar[basis] != lab:
                gates.append("G45.3 the committed %s barrier recomputes to %.4f, not %s" % (basis, bar[basis], lab))
            if CT[(basis, "R")]["nbf"] != CT[(basis, "TS")]["nbf"]:
                gates.append("G45.3 the two committed %s legs do not report one number of basis functions" % basis)
            if CT[(basis, "R")]["block"] != CT[(basis, "TS")]["block"]:
                gates.append("G45.3 the two committed %s legs do not print one DLPNO SETTINGS block" % basis)
        if not gates and "%+.4f" % (bar["def2-TZVPD"] - bar["def2-SVPD"]) != SHIFT:
            gates.append("G45.3 the committed basis shift recomputes to %+.4f, not %s" % (bar["def2-TZVPD"] - bar["def2-SVPD"], SHIFT))

    # ------------------------------------------------------------------ G45.4 the quoted sentences
    pre = flat(open(a.prereg, errors="replace").read())
    for fn, q in QUOTES:
        p = os.path.join(S6, fn)
        if not os.path.isfile(p) or flat(q) not in flat(open(p, errors="replace").read()):
            gates.append("G45.4 %s does not hold the sentence: %s" % (fn, q))
        if flat(q) not in pre:
            gates.append("G45.4 the pre-registration does not hold the sentence it is read under: %s" % q)

    # the nine functionals at def2-TZVPD, from the committed tables
    FN = []
    p37, p37a = os.path.join(S6, "d37.tsv"), os.path.join(S6, "d37a1.tsv")
    if os.path.isfile(p37) and os.path.isfile(p37a):
        for l in open(p37).read().splitlines()[1:]:
            f = l.split("\t")
            if len(f) >= 6 and f[2] == "def2-TZVPD" and f[5]:
                FN.append((f[0], f[5]))
        for l in open(p37a).read().splitlines()[1:]:
            f = l.split("\t")
            if len(f) >= 10 and f[0] == "RI" and f[1] == "def2-TZVPD" and f[4] and f[9] == "yes":
                FN.append(("revDSD-PBEP86-D4", f[4]))
    if len(FN) != 9 or len(set(n for n, _ in FN)) != 9:
        gates.append("G45.3 d37.tsv and d37a1.tsv do not give nine functionals at def2-TZVPD")
    elif "def2-TZVPD" in bar:
        cam = [b for n, b in FN if n == "CAM-B3LYP-D3BJ"]
        if len(cam) != 1 or "%+.4f" % (float(cam[0]) - bar["def2-TZVPD"]) != CAM_DIST:
            gates.append("G45.3 CAM-B3LYP-D3BJ's committed def2-TZVPD barrier minus D36's is not the %s AMENDMENT 7 prints" % CAM_DIST)

    if gates:
        say("GATES")
        for g in gates:
            say("   FAILED  %s" % g)
        say("")
        say("A FAILURE OF G45.1 TO G45.4 STOPS THE DIAGNOSTIC. No reading is made.")
        print("\n".join(rep))
        stop("%d gate(s) failed. Nothing is written." % len(gates))

    NBF = {b: CT[(b, "R")]["nbf"] for b in ("def2-SVPD", "def2-TZVPD")}
    geo = {st: d38.xyz(os.path.join(S6, prep.GEOM[st])) for st in prep.STATES}
    places, per = job_facts(D)

    # ------------------------------------------------------------------ G45.5 each run
    R = {}
    for run, kind, basis, st in runs:
        z = dict(run=run, kind=kind, basis=basis, st=st, why=[], status="NOT RUN", nbf=None, hf=None, cor=None,
                 trip=None, tot=None, tcut=None, t0=None, it=None, runtime=None, diff=[], dhf=None, dccsd=None)
        R[run] = z
        r = d38.read_out(os.path.join(D, run))
        if r is None:
            z["why"].append("no output")
            continue
        t = open(r["path"], errors="replace").read()
        w = z["why"]
        z["runtime"], z["nbf"] = r["runtime"], nbf(t)
        if not r["norm"]:
            w.append("did not terminate normally. ORCA's words: %s" % orca_words(t))
        if r["notconv"]:
            w.append("did not converge its SCF")
        if r["nel"] != NEL:
            w.append("reports %s electrons, not %d" % (r["nel"], NEL))
        if r["atoms"] != geo[st]:
            w.append("ran on coordinates that are not frame 42436's %s geometry" % st)
        e = echoed(t)
        built = prep.build_input(kind, basis).split("\n")
        while built and not built[-1].strip():
            built.pop()
        if e is None or [l.rstrip() for l in e] != built:
            w.append("the input ORCA prints is not the committed input")
        if basis in NBF and z["nbf"] != NBF[basis]:
            w.append("reports %s basis functions; the committed %s runs report %d" % (z["nbf"], basis, NBF[basis]))
        ctl = C.get((basis, st))
        if kind == "MP2":
            m = mp2_of(t)
            z["hf"], z["cor"], z["tot"] = m["scf"], m["cor"], m["tot"]
            w.extend(m["bad"])
            for lab, v in (("an SCF Total Energy", m["scf"]), ("RI-MP2 CORRELATION ENERGY", m["cor"]), ("MP2 TOTAL ENERGY", m["tot"])):
                if v is None:
                    w.append("the line for %s was not found" % lab)
            if None not in (m["scf"], m["cor"], m["tot"]):
                if abs(m["scf"] + m["cor"] - m["tot"]) > TOL_ID:
                    w.append("its SCF and correlation energies do not add up to its MP2 total within %.0e Eh" % TOL_ID)
                if r["E"] is None or abs(r["E"] - m["tot"]) > TOL_ID:
                    w.append("its last FINAL SINGLE POINT ENERGY is not its MP2 total energy")
            if m["rng"] is None:
                w.append("the line for the orbital range of the correlation treatment was not found")
            elif m["rng"][0] != NCORE or m["rng"][1] != NOCC - 1 or m["rng"][2] != NOCC:
                w.append("correlates orbitals %d to %d, not %d to %d: the core is not frozen as in the benchmark"
                         % (m["rng"][0], m["rng"][1], NCORE, NOCC - 1))
        else:
            z["hf"], z["cor"], z["trip"], z["tot"] = r["e0"], r["corr"], r["tcorr"], r["ccsdt"]
            w.extend(r["inconsistent"])
            if r["ccsdt"] is None or None in (r["e0"], r["corr"], r["ccsd"], r["tcorr"]):
                w.append("does not print every component of a CCSD(T) energy under the labels the committed outputs use")
            else:
                if r["E"] is None or abs(r["E"] - r["ccsdt"]) > TOL_ID:
                    w.append("its last FINAL SINGLE POINT ENERGY is not its E(CCSD(T))")
                if abs(r["e0"] + r["corr"] - r["ccsd"]) > TOL_ID or abs(r["ccsd"] + r["tcorr"] - r["ccsdt"]) > TOL_ID:
                    w.append("its components do not add up to its E(CCSD(T)) within %.0e Eh" % TOL_ID)
                z["dccsd"] = r["ccsd"] - ctl["ccsd"]
            blk, cb = settings(t), CT[(basis, st)]["block"]
            z["tcut"] = setting(blk, "TCutPNO")
            z["t0"], z["it"] = triples_times(t)
            if not blk:
                w.append("prints no DLPNO SETTINGS block")
            else:
                if z["tcut"] != prep.TCUT[kind]:
                    w.append("prints TCutPNO: %s, not the %s asked for" % (z["tcut"], prep.TCUT[kind]))
                for k in KEEP:
                    if setting(blk, k) != setting(cb, k):
                        w.append("prints %s: %s; the committed run prints %s" % (k, setting(blk, k), setting(cb, k)))
                z["diff"] = [l for l in blk if l not in cb]
            if kind == "T1":
                if z["it"] is None:
                    w.append("the line for the time of the iterative (T) was not found")
                elif z["it"] <= 0.0:
                    w.append("prints 0.000 sec for the iterative (T): by ORCA's own timing the iterative triples did not run")
            else:
                if z["it"] is None or z["t0"] is None:
                    w.append("the lines for the times of the triples were not found")
                elif z["it"] != 0.0 or z["t0"] <= 0.0:
                    w.append("its triples are not the non-iterative ones of the committed runs, by ORCA's own timing")
        if ctl is not None and z["hf"] is not None:
            z["dhf"] = z["hf"] - ctl["e0"]
            if abs(z["dhf"]) > TOL_REP:
                w.append("its Hartree-Fock energy differs from the committed %s's by %+.2e Eh" % (CTL[(basis, st)], z["dhf"]))
        z["status"] = "NOT READ" if w else "measured"

    def ok(*names):
        return all(R[n]["status"] == "measured" for n in names)

    def dE(key, top, bot):
        return (R[top][key] - R[bot][key]) * H2K

    nmeas = sum(1 for z in R.values() if z["status"] == "measured")
    say("=" * 78)
    say("GATES G45.1 TO G45.4: ALL PASS")
    say("=" * 78)
    say("")
    say("   G45.1  all %d files listed in d45/INPUT_CHECKSUMS.sha256 have the sha256 listed, and" % len(listed))
    say("          every input, RUNS.tsv and run.pbs are what cs6_d45_prepare.py writes from the")
    say("          committed files.")
    say("   G45.2  all %d files of d45/SOURCES.sha256 have the sha256 listed, and no committed" % len(prep.sources()))
    say("          file under %s differs from HEAD." % ", ".join(GUARDED))
    say("   G45.3  the four committed legs end normally, add up, print TCutPNO: 1.000e-07 and")
    say("          0.000 sec for the iterative (T), and give %s and %s kcal/mol," % (B_SV, B_TZ))
    say("          %s apart." % SHIFT)
    say("          The nine functionals' def2-TZVPD barriers are in d37.tsv and d37a1.tsv.")
    say("   G45.4  the %d sentences quoted from committed documents are in them word for word." % len(QUOTES))
    say("   G45.5  is per run, below.")
    say("")

    # ------------------------------------------------------------------ R45.1
    say("=" * 78)
    say("R45.1  WHAT RAN")
    say("=" * 78)
    say("")
    say("   run       basis         functions  status      TCutPNO     T0 s   iter (T) s  run time")
    for run, kind, basis, st in runs:
        z = R[run]
        rt = ""
        if z["runtime"]:
            k = re.match(r"(\d+) days (\d+) hours (\d+) minutes (\d+) seconds", z["runtime"])
            rt = "%dh%02dm" % (int(k.group(1)) * 24 + int(k.group(2)), int(k.group(3))) if k else z["runtime"][:12]
        say("   %-9s %-13s %9s  %-10s  %-9s %8s %10s   %s" % (
            run, basis, z["nbf"] if z["nbf"] is not None else "", z["status"], z["tcut"] or "",
            "%.0f" % z["t0"] if z["t0"] is not None else "", "%.0f" % z["it"] if z["it"] is not None else "", rt))
    say("")
    say("   %d of %d runs are measured." % (nmeas, len(runs)))
    for run, kind, basis, st in runs:
        for wl in R[run]["why"]:
            say("   %s %s: %s" % (run, R[run]["status"], wl))
    say("   The four committed legs:")
    for k in CTL:
        say("      %-11s %4d functions, TCutPNO %s, T0 %.0f s, iterative (T) %.0f s" % (
            CTL[k], CT[k]["nbf"], setting(CT[k]["block"], "TCutPNO"), CT[k]["t0"], CT[k]["it"]))
    odd = [(run, R[run]["diff"]) for run, kind, _, _ in runs if R[run]["diff"]]
    if odd:
        say("   Lines of ORCA's DLPNO SETTINGS block that are not in the committed run's block:")
        for run, ls in odd:
            say("      %-8s %s" % (run, "; ".join(ls)))
    say("")

    # ------------------------------------------------------------------ R45.2
    say("=" * 78)
    say("R45.2  THE SAME REFERENCE: each new Hartree-Fock energy beside the committed one")
    say("=" * 78)
    say("")
    say("   A new run is used only if its Hartree-Fock energy is within 1e-6 Eh of the")
    say("   committed run's at the same basis set and state (the figure of D44, section 2).")
    say("")
    say("   run       committed leg   new minus committed, Eh    word")
    for run, kind, basis, st in runs:
        z = R[run]
        if (basis, st) not in CTL:
            say("   %-9s %-15s %-26s %s" % (run, "none", "", "no committed counterpart at this basis set"))
        elif z["dhf"] is None:
            say("   %-9s %-15s %-26s %s" % (run, CTL[(basis, st)], "", "not formed"))
        else:
            say("   %-9s %-15s %+.3e %15s %s" % (run, CTL[(basis, st)], z["dhf"], "", "REPRODUCES" if abs(z["dhf"]) <= TOL_REP else "DIFFERS"))
    say("")
    say("   E(CCSD) of each new coupled-cluster run minus the committed leg's, Eh. For T1 the")
    say("   thresholds are the committed run's, so this should be zero; for P6 and P8 it is")
    say("   the effect of the threshold on CCSD.")
    for run, kind, basis, st in runs:
        if kind != "MP2" and R[run]["dccsd"] is not None:
            say("      %-8s %+.9f" % (run, R[run]["dccsd"]))
    say("")

    Q = []                                    # d45_barriers.tsv

    def keep(name, v, label, src):
        Q.append((name, v, label, src))

    b7, b36 = bar["def2-SVPD"], bar["def2-TZVPD"]
    keep("barrier CCSD(T0) TightPNO def2-SVPD", b7, "measured, committed (D41)", "d38/N, d41/T")
    keep("barrier CCSD(T0) TightPNO def2-TZVPD", b36, "measured, committed (D36)", "d36a2/B_R, d36a2/B_TS")
    cparts = {k: (C[("def2-SVPD", "TS")][k] - C[("def2-SVPD", "R")][k]) * H2K for _, k in COMP}

    # ------------------------------------------------------------------ R45.3
    say("=" * 78)
    say("R45.3  THE TRIPLES: iterative (T1) beside the committed non-iterative (T0), def2-SVPD")
    say("=" * 78)
    say("")
    dT = None
    if ok("T1_R", "T1_TS"):
        say("   state   (T0), committed, Eh     (T1), here, Eh      T1 minus T0, Eh")
        for st in prep.STATES:
            c, z = C[("def2-SVPD", st)], R["T1_%s" % st]
            say("   %-5s %20.9f %20.9f %+18.9f" % (st, c["tcorr"], z["trip"], z["trip"] - c["tcorr"]))
        bT1 = dE("tot", "T1_TS", "T1_R")
        tT1 = dE("trip", "T1_TS", "T1_R")
        dT = bT1 - b7
        say("")
        say("   the barrier with (T0), committed            %10.4f kcal/mol" % b7)
        say("   the barrier with (T1)                       %10.4f kcal/mol" % bT1)
        say("   T1 MINUS T0 ON THE BARRIER                  %+10.4f kcal/mol   measured" % dT)
        say("   of which the triples' own part              %+10.4f   (the rest is E(CCSD), R45.2)" % (tT1 - cparts["tcorr"]))
        say("   the triples' part of the barrier: %+.4f with (T0), %+.4f with (T1)" % (cparts["tcorr"], tT1))
        keep("barrier CCSD(T1) TightPNO def2-SVPD", bT1, "measured", "d45/T1_R, d45/T1_TS")
        keep("triples, T1 minus T0, def2-SVPD", dT, "measured", "d45/T1_*, d38/N, d41/T")
    else:
        say("   NOT FORMED: T1_R and T1_TS are not both measured.")
    say("   Reported, not interpreted. No threshold is set on it.")
    say("")

    # ------------------------------------------------------------------ R45.4
    say("=" * 78)
    say("R45.4  THE PNO THRESHOLD: TCutPNO 1e-6, 1e-7 (committed) and 1e-8, def2-SVPD, (T0)")
    say("=" * 78)
    say("")
    b6 = dE("tot", "P6_TS", "P6_R") if ok("P6_R", "P6_TS") else None
    b8 = dE("tot", "P8_TS", "P8_R") if ok("P8_R", "P8_TS") else None
    say("   TCutPNO     E(CCSD(T)) R, Eh    E(CCSD(T)) TS, Eh    HF    CCSD corr   triples    barrier")
    for lab, tag in (("1e-6", "P6"), ("1e-7", None), ("1e-8", "P8")):
        if tag is None:
            c = C[("def2-SVPD", "R")], C[("def2-SVPD", "TS")]
            say("   %-6s %19.9f %20.9f %9.4f %9.4f %9.4f %10.4f   committed" % (
                lab, c[0]["ccsdt"], c[1]["ccsdt"], cparts["e0"], cparts["corr"], cparts["tcorr"], b7))
        elif ok(tag + "_R", tag + "_TS"):
            say("   %-6s %19.9f %20.9f %9.4f %9.4f %9.4f %10.4f" % (
                lab, R[tag + "_R"]["tot"], R[tag + "_TS"]["tot"], dE("hf", tag + "_TS", tag + "_R"),
                dE("cor", tag + "_TS", tag + "_R"), dE("trip", tag + "_TS", tag + "_R"), dE("tot", tag + "_TS", tag + "_R")))
        else:
            say("   %-6s NOT FORMED: %s_R and %s_TS are not both measured" % (lab, tag, tag))
    d67 = b7 - b6 if b6 is not None else None
    d78 = b8 - b7 if b8 is not None else None
    cps = b6 + CPS_F * (b7 - b6) if b6 is not None else None
    say("")
    say("   1e-7 minus 1e-6 on the barrier, d67         %s kcal/mol   measured" % f4(d67))
    say("   1e-8 minus 1e-7 on the barrier, d78         %s kcal/mol   measured" % f4(d78))
    say("   two-point estimate from 1e-6 and 1e-7       %s kcal/mol   ESTIMATE" % f4(cps, False))
    say("      (b6 + %.1f x (b7 - b6), the formula of Altun, Neese and Bistoni 2020)" % CPS_F)
    say("   the barrier at 1e-8 minus that estimate     %s kcal/mol" % f4(b8 - cps if None not in (b8, cps) else None))
    if None in (d67, d78):
        word = "NOT ESTABLISHED"
    else:
        word = "CONVERGING" if abs(d78) < abs(d67) else "NOT CONVERGING"
    say("")
    say("   THE PNO SERIES IS %s." % word)
    if word == "NOT ESTABLISHED":
        say("   The rule needs both differences, and at least one is not formed.")
    else:
        say("   The rule fixed in advance: CONVERGING if |d78| is smaller than |d67|, and NOT")
        say("   CONVERGING otherwise. |d78| = %.4f and |d67| = %.4f." % (abs(d78), abs(d67)))
        say("   The two steps are %s." % ("in the same direction" if d67 * d78 > 0 else "NOT in the same direction"))
    for st in prep.STATES:
        if ok("P6_" + st, "P8_" + st):
            c = C[("def2-SVPD", st)]
            e = [R["P6_" + st]["tot"] - R["P6_" + st]["hf"], c["ccsdt"] - c["e0"], R["P8_" + st]["tot"] - R["P8_" + st]["hf"]]
            say("   %s: the correlation energy recovered is %.6f, %.6f, %.6f Eh at 1e-6, 1e-7, 1e-8:" % (st, e[0], e[1], e[2]))
            say("      it %s at each tighter threshold." % ("is larger in size" if e[0] > e[1] > e[2] else "IS NOT larger in size"))
    if b6 is not None:
        keep("barrier CCSD(T0) TCutPNO 1e-6 def2-SVPD", b6, "measured", "d45/P6_R, d45/P6_TS")
        keep("PNO, 1e-7 minus 1e-6, def2-SVPD", d67, "measured", "d45/P6_*, d38/N, d41/T")
        keep("barrier, two-point PNO estimate, def2-SVPD", cps, "ESTIMATE", "b6 + 1.5 x (b7 - b6)")
    if b8 is not None:
        keep("barrier CCSD(T0) TCutPNO 1e-8 def2-SVPD", b8, "measured", "d45/P8_R, d45/P8_TS")
        keep("PNO, 1e-8 minus 1e-7, def2-SVPD", d78, "measured", "d45/P8_*, d38/N, d41/T")
    say("")

    # ------------------------------------------------------------------ R45.5
    say("=" * 78)
    say("R45.5  THE BASIS SET, AT MP2")
    say("=" * 78)
    say("")
    M = {}
    say("   basis          functions     HF      MP2 corr    MP2 barrier, kcal/mol")
    for key, basis in prep.MP2:
        if ok(key + "_R", key + "_TS"):
            M[key] = dict(hf=dE("hf", key + "_TS", key + "_R"), cor=dE("cor", key + "_TS", key + "_R"),
                          tot=dE("tot", key + "_TS", key + "_R"))
            say("   %-13s %9s %10.4f %10.4f %12.4f" % (basis, R[key + "_R"]["nbf"], M[key]["hf"], M[key]["cor"], M[key]["tot"]))
            keep("barrier RI-MP2 %s" % basis, M[key]["tot"], "measured", "d45/%s_R, d45/%s_TS" % (key, key))
        else:
            say("   %-13s NOT FORMED: %s_R and %s_TS are not both measured" % (basis, key, key))
    say("")
    steps = (("M_SV", "M_TZ"), ("M_TZ", "M_TP"), ("M_TP", "M_QZ"), ("M_TZ", "M_QZ"))
    names = dict(prep.MP2)
    for lo, hi in steps:
        v = M[hi]["tot"] - M[lo]["tot"] if lo in M and hi in M else None
        say("   MP2, %-11s minus %-11s %s kcal/mol%s" % (names[hi], names[lo], f4(v), "   measured" if v is not None else ""))
    say("   CCSD(T), def2-TZVPD minus def2-SVPD     %s kcal/mol   measured, committed (D41)" % SHIFT)
    cbs = None
    if "M_TP" in M and "M_QZ" in M:
        ex, ey = math.exp(-ALPHA34 * math.sqrt(3.0)), math.exp(-ALPHA34 * math.sqrt(4.0))
        hf_c = (M["M_TP"]["hf"] * ey - M["M_QZ"]["hf"] * ex) / (ey - ex)
        co_c = (3.0 ** BETA34 * M["M_TP"]["cor"] - 4.0 ** BETA34 * M["M_QZ"]["cor"]) / (3.0 ** BETA34 - 4.0 ** BETA34)
        cbs = hf_c + co_c
        say("")
        say("   EXTRAPOLATED from def2-TZVPPD and def2-QZVPPD, with the exponents ORCA's manual")
        say("   gives for def2-TZVPP and def2-QZVPP, alpha = %.2f and beta = %.2f:" % (ALPHA34, BETA34))
        say("      HF %10.4f   MP2 corr %10.4f   MP2 barrier %10.4f kcal/mol   ESTIMATE" % (hf_c, co_c, cbs))
        say("      that estimate minus MP2/def2-QZVPPD   %+.4f kcal/mol" % (cbs - M["M_QZ"]["tot"]))
        keep("barrier RI-MP2, (3/4) extrapolation", cbs, "ESTIMATE", "def2-TZVPPD, def2-QZVPPD; alpha 7.88, beta 2.97")
    say("")
    say("   THE PREMISE of adding an MP2 basis step to a coupled-cluster barrier is that")
    say("   CCSD(T) minus MP2 does not depend on the basis set. Measured over the one step")
    say("   the record has:")
    dl = {}
    for key, basis, b in (("M_SV", "def2-SVPD", b7), ("M_TZ", "def2-TZVPD", b36)):
        dl[key] = b - M[key]["tot"] if key in M else None
        say("      CCSD(T) minus MP2 at %-11s %s kcal/mol" % (basis, f4(dl[key])))
    prem = dl["M_TZ"] - dl["M_SV"] if None not in (dl["M_SV"], dl["M_TZ"]) else None
    say("      its change from def2-SVPD to def2-TZVPD   %s kcal/mol" % f4(prem))
    if prem is not None:
        keep("CCSD(T) minus MP2, change from def2-SVPD to def2-TZVPD", prem, "measured", "d45/M_SV_*, d45/M_TZ_*, committed legs")
    say("")

    # ------------------------------------------------------------------ R45.6
    say("=" * 78)
    say("R45.6  D36'S BARRIER WITH THE MEASURED CORRECTIONS BESIDE IT")
    say("=" * 78)
    say("")
    c_pno_m = d78
    c_pno_e = (cps - b7) if cps is not None else None
    c_bas_m = M["M_QZ"]["tot"] - M["M_TZ"]["tot"] if "M_QZ" in M and "M_TZ" in M else None
    c_bas_e = cbs - M["M_TZ"]["tot"] if cbs is not None and "M_TZ" in M else None
    say("   D36's barrier, DLPNO-CCSD(T0), TightPNO, def2-TZVPD        %10.4f   measured, committed" % b36)
    say("   triples, T1 minus T0, at def2-SVPD (R45.3)                 %10s   %s" % (f4(dT), "measured" if dT is not None else ""))
    say("   PNO threshold, 1e-8 minus 1e-7, at def2-SVPD (R45.4)       %10s   %s" % (f4(c_pno_m), "measured" if c_pno_m is not None else ""))
    say("   PNO threshold, two-point estimate minus 1e-7, def2-SVPD    %10s   %s" % (f4(c_pno_e), "ESTIMATE" if c_pno_e is not None else ""))
    say("   basis set, MP2/def2-QZVPPD minus MP2/def2-TZVPD (R45.5)    %10s   %s" % (f4(c_bas_m), "measured" if c_bas_m is not None else ""))
    say("   basis set, MP2 extrapolated minus MP2/def2-TZVPD           %10s   %s" % (f4(c_bas_e), "ESTIMATE" if c_bas_e is not None else ""))
    cA = b36 + dT + c_pno_m + c_bas_m if None not in (dT, c_pno_m, c_bas_m) else None
    cB = b36 + dT + c_pno_e + c_bas_e if None not in (dT, c_pno_e, c_bas_e) else None
    say("")
    say("   COMPOSITE A, D36's barrier plus the three measured steps    %10s kcal/mol" % f4(cA, False))
    say("   COMPOSITE B, D36's barrier plus the triples and the two")
    say("                extrapolated corrections                       %10s kcal/mol" % f4(cB, False))
    say("")
    say("   BOTH ARE ESTIMATES. Each adds corrections that were measured one at a time, two")
    say("   of them at the smaller basis set and one at MP2, and whether they add is not")
    say("   measured. Neither replaces D36's number, and no conclusion rests on either.")
    if cA is not None:
        keep("composite A", cA, "ESTIMATE", "23.8603 + triples + (1e-8 minus 1e-7) + MP2 (QZVPPD minus TZVPD)")
    if cB is not None:
        keep("composite B", cB, "ESTIMATE", "23.8603 + triples + two-point PNO + MP2 (extrapolated minus TZVPD)")
    say("")

    # ------------------------------------------------------------------ R45.7
    say("=" * 78)
    say("R45.7  THE NINE FUNCTIONALS' COMMITTED def2-TZVPD BARRIERS BESIDE EACH")
    say("=" * 78)
    say("")
    say("   functional          barrier    minus D36's    minus composite A    minus composite B")
    for n, b in FN:
        v = float(b)
        say("   %-18s %8s %+13.4f %18s %20s" % (n, b, v - b36, f4(v - cA if cA is not None else None), f4(v - cB if cB is not None else None)))
        keep("%s def2-TZVPD minus D36's barrier" % n, v - b36, "derived from committed files", "d37.tsv, d37a1.tsv, d36a2")
        if cA is not None:
            keep("%s def2-TZVPD minus composite A" % n, v - cA, "ESTIMATE", "")
        if cB is not None:
            keep("%s def2-TZVPD minus composite B" % n, v - cB, "ESTIMATE", "")
    say("")
    say("   By distance, nearest first:")
    for lab, ref in (("from D36's barrier", b36), ("from composite A", cA), ("from composite B", cB)):
        if ref is None:
            say("      %-20s NOT FORMED" % lab)
            continue
        o = sorted((n for n, _ in FN), key=lambda n: abs(float(dict(FN)[n]) - ref))
        say("      %-20s %s" % (lab, ", ".join(o[:5]) + ","))
        say("      %-20s %s" % ("", ", ".join(o[5:])))
    say("   The first column of differences is the record's. The other two are differences")
    say("   from estimates and are reported, not interpreted: the list of AMENDMENT 7, A7.3 is")
    say("   not redrawn on them, and no functional is added to it or dropped from it.")
    say("")

    # ------------------------------------------------------------------ R45.8
    say("=" * 78)
    say("R45.8  WHERE IT RAN, AND WHAT IT USED")
    say("=" * 78)
    say("")
    if not places:
        say("   no log of the job script was found.")
    for pl in places:
        say("   %-22s node %s, job %s" % (pl["log"], pl["host"] or "not reported", pl["jobid"] or "not reported"))
        say("      scratch %s" % (pl["scratch"] or "not reported"))
        say("      %s free when it was chosen, margin %s" % (
            "%d GB" % pl["avail"] if pl["avail"] is not None else "not reported", pl["margin"] or "not reported"))
        say("      its last line: %s" % (pl["finished"] or "NOT REACHED: the job did not run to its end"))
        for c, why in pl["refused"]:
            say("      candidate refused: %s %s" % (c, why))
    say("   The maximum of the scratch in use, sampled every 60 s while ORCA ran; a lower")
    say("   bound on the true peak:")
    for run, kind, basis, st in runs:
        if run in per:
            say("      %-8s %7d MB   ORCA's exit status %d   (%s)" % (run, per[run]["peak"], per[run]["rc"], per[run]["log"]))
        else:
            say("      %-8s not reported" % run)
    say("")

    # ------------------------------------------------------------------ what follows, and what it cannot say
    say("=" * 78)
    say("WHAT FOLLOWS (section 6 of the pre-registration)")
    say("=" * 78)
    say("")
    if word == "CONVERGING":
        say("   F1. THE PNO SERIES IS CONVERGING. DIAGNOSTIC 46 may be submitted once this result")
        say("   is committed, as AMENDMENT 8 says.")
    else:
        say("   F2. THE PNO SERIES IS %s. DIAGNOSTIC 46 is not submitted. What is done" % word)
        say("   about it is a dated amendment, written first, and a decision for Jessica.")
    if nmeas < len(runs):
        say("   F3. %d run(s) are not measured. The cause is diagnosed first. A job that was cut" % (len(runs) - nmeas))
        say("   short is submitted again as it is; anything that changes an input is a dated")
        say("   amendment.")
    say("   F5. D36's 23.8603 and D41's 22.8968 stand as committed. Nothing here edits them.")
    say("")
    say("=" * 78)
    say("WHAT THIS CANNOT SAY (section 7 of the pre-registration, unchanged)")
    say("=" * 78)
    say("")
    for l in ("D45.1  One frame, frozen geometries, in vacuo, no charges. It is not the experimental",
              "       number and not a prediction of a physical barrier.",
              "D45.2  It names no design, says nothing about ddE and adopts no level of theory.",
              "D45.3  (T1) is still a perturbative, local treatment of the triples. Nothing here",
              "       measures the error of CCSD(T) itself.",
              "D45.4  The triples and the PNO threshold are checked at def2-SVPD only, and TCutPairs,",
              "       TCutMKN and TCutDO are not varied.",
              "D45.5  MP2 is used for its dependence on the basis set and for nothing else. The two",
              "       exponents were made for def2-TZVPP and def2-QZVPP, without diffuse functions.",
              "D45.6  The auxiliary basis is AutoAux and the core is frozen in every run; neither",
              "       error is measured.",
              "D45.7  A composite is an estimate and never carries a conclusion.",
              "D45.8  Nothing here is shaped toward an outcome. Every reading is printed whichever",
              "       way it came out."):
        say("   %s" % l)
    say("")
    say("D45: %d of %d runs measured; the PNO series is %s; triples %s, composite A %s, composite B %s kcal/mol."
        % (nmeas, len(runs), word, f4(dT), f4(cA, False), f4(cB, False)))

    L = ["run\tkind\tbasis\tstate\tstatus\tn_basis\tE_HF_Eh\tE_corr_Eh\tE_triples_Eh\tE_total_Eh\tTCutPNO\tT0_s\titerative_T_s\t"
         "run_time\tpeak_scratch_MB\twhy_not_measured"]
    for key, rel in CTL.items():
        c = C[key]
        L.append("\t".join([rel, "committed", key[0], key[1], "committed", str(CT[key]["nbf"]), "%.9f" % c["e0"], "%.9f" % c["corr"],
                            "%.9f" % c["tcorr"], "%.9f" % c["ccsdt"], setting(CT[key]["block"], "TCutPNO"), "%.3f" % CT[key]["t0"],
                            "%.3f" % CT[key]["it"], c["runtime"] or "", "", ""]))
    for run, kind, basis, st in runs:
        z = R[run]
        g = z["status"] == "measured"
        L.append("\t".join([run, kind, basis, st, z["status"], str(z["nbf"]) if z["nbf"] is not None else "",
                            "%.9f" % z["hf"] if g else "", "%.9f" % z["cor"] if g else "",
                            "%.9f" % z["trip"] if g and z["trip"] is not None else "", "%.9f" % z["tot"] if g else "",
                            z["tcut"] or "", "%.3f" % z["t0"] if z["t0"] is not None else "", "%.3f" % z["it"] if z["it"] is not None else "",
                            z["runtime"] or "", str(per[run]["peak"]) if run in per else "", " | ".join(z["why"])]))
    B = ["quantity\tkcal_mol\tlabel\tformed_from"] + ["%s\t%.4f\t%s\t%s" % q for q in Q]
    rep[:] = [l.rstrip() for l in rep]
    open(a.tsv, "w").write("\n".join(L) + "\n")
    open(os.path.splitext(a.tsv)[0] + "_barriers.tsv", "w").write("\n".join(B) + "\n")
    open(a.result, "w").write("\n".join(rep) + "\n")
    print("\n".join(rep))
    return 0


if __name__ == "__main__":
    sys.exit(main())
