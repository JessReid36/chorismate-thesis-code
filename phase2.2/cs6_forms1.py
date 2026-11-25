#!/usr/bin/env python3
"""cs6_forms1.py - the reader of STAGE 6, FORMS TEST 1, committed before any of its outputs exists.

It answers one question for each of fifteen runs: did ORCA 6.0.1 on hpc1 accept this form of input,
and does the output show the form did what it was written to do (the core potential read on the
right centre, the functions of a ghost present, the right count of electrons and of point charges,
a gradient file, a potential file). It names no design and forms no change of a barrier. The one
energy it compares is that of the control with one electron removed, against D56's committed one.

usage (as cs6_forms1_harvest.sh calls it):
   python3 cs6_forms1.py --res <results repo> --code <code repo> --prereg <pre-registration>
                         --result <result document to write> --tsv <table to write>
If a gate fails it prints why and writes nothing.
"""
import argparse, csv, hashlib, importlib.util, os, re, subprocess, sys

PIN_PREREG = "71cb659f1e31c3242525e79f76939a75041206787adf757bc85fff308188ba4f"
S6REL = "phase2.2/cs_stage6"
DIR = "forms1"
NBF = {"SV": 402, "TZ": 682}          # D56's G56.6: the basis functions of the substrate at each level
TOL_E = 1.0e-6                        # Eh: the control with one electron removed reproduces D56's
W = 96

# run -> what is expected: electrons, point charges read, the centre with a core potential, functions beyond the substrate's
EXPECT = {
    "SV_wallpos_R": dict(nel=118, npc=1, ecp="He>", extra=False),
    "SV_shi_R": dict(nel=118, npc=None, ecp="Na>", extra=False),
    "TZ_shi_R": dict(nel=118, npc=None, ecp="Na>", extra=False),
    "TZ_wallHe_R": dict(nel=118, npc=None, ecp="He>", extra=False),
    "TZ_wallNe_R": dict(nel=118, npc=None, ecp="Ne>", extra=False),
    "SV_Na_R": dict(nel=128, npc=None, ecp=None, extra=True),
    "TZ_Na_R": dict(nel=128, npc=None, ecp=None, extra=True),
    "SV_K_R": dict(nel=136, npc=None, ecp=None, extra=True),
    "TZ_K_R": dict(nel=136, npc=None, ecp=None, extra=True),
    "SV_ghostNa_R": dict(nel=118, npc=None, ecp=None, extra=True),
    "TZ_ghostNa_R": dict(nel=118, npc=None, ecp=None, extra=True),
    "SV_ghostK_R": dict(nel=118, npc=None, ecp=None, extra=True),
    "TZ_ghostK_R": dict(nel=118, npc=None, ecp=None, extra=True),
    "SV_engrad_R": dict(nel=118, npc=2, ecp=None, extra=False),
    "SV_mono_R": dict(nel=117, npc=None, ecp=None, extra=False),
}


def sha(p):
    with open(p, "rb") as f:
        return hashlib.sha256(f.read()).hexdigest()


def rd(p):
    with open(p, errors="replace") as f:
        return f.read()


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def stop(why):
    print("THE READER STOPPED AT A GATE: " + why)
    print("Nothing was written.")
    sys.exit(1)


def echoed_input(text):
    """the input lines as ORCA echoes them, up to END OF INPUT"""
    out = []
    for l in text.split("\n"):
        if "****END OF INPUT****" in l:
            break
        m = re.match(r"^\|\s*\d+>( (.*))?$", l)
        if m:
            out.append((m.group(2) or "").rstrip())
    return out


def num(pattern, text):
    m = re.search(pattern, text, re.M)
    return int(m.group(1)) if m else None


def engrad(path):
    """(number of atoms, number of gradient values) of an ORCA .engrad file, or None"""
    if not os.path.isfile(path):
        return None
    L = [l.strip() for l in rd(path).split("\n")]
    try:
        i = L.index("# Number of atoms")
        n = int(L[i + 2])
        j = L.index("# The current gradient in Eh/bohr")
        g = []
        k = j + 2
        while k < len(L) and not L[k].startswith("#") and L[k]:
            g.append(float(L[k]))
            k += 1
        return n, len(g)
    except (ValueError, IndexError):
        return None


def pcgrad(path):
    if not os.path.isfile(path):
        return "NOT WRITTEN"
    L = [l.split() for l in rd(path).split("\n") if l.strip()]
    try:
        if int(L[0][0]) == 2 and len(L) >= 3 and all(len(r) == 3 for r in L[1:3]):
            [float(x) for r in L[1:3] for x in r]
            return "WRITTEN, 2 ROWS OF 3"
    except (ValueError, IndexError):
        pass
    return "WRITTEN, NOT IN THE FORM EXPECTED"


def vpot(d, pts):
    p = os.path.join(d, "vpot.out")
    if not os.path.isfile(p):
        return "NOT WRITTEN"
    L = [l.split() for l in rd(p).split("\n") if l.strip()]
    P = [l.split() for l in pts.split("\n")[1:] if l.strip()]
    try:
        if int(L[0][0]) != 5 or len(L) != 6 or any(len(r) != 4 for r in L[1:]):
            return "WRITTEN, NOT IN THE FORM EXPECTED"
        for r, q in zip(L[1:], P):
            if max(abs(float(a) - float(b)) for a, b in zip(r[:3], q)) > 1e-6:
                return "WRITTEN, POINTS DO NOT MATCH"
            float(r[3])
    except (ValueError, IndexError):
        return "WRITTEN, NOT IN THE FORM EXPECTED"
    return "WRITTEN, 5 POINTS"


def main():
    ap = argparse.ArgumentParser()
    for a in ("res", "code", "prereg", "result", "tsv"):
        ap.add_argument("--" + a, required=True)
    a = ap.parse_args()
    S6 = os.path.join(a.res, S6REL)
    D = os.path.join(S6, DIR)

    # GF.1 the pre-registration is the committed one
    if not os.path.isfile(a.prereg) or sha(a.prereg) != PIN_PREREG:
        stop("the pre-registration given is not the committed STAGE6_FORMS_TEST1_20261125.txt.")
    # GF.2 the programs are the committed ones
    for f in ("cs6_forms1_prepare.py", "cs6_forms1.py"):
        q = subprocess.run(["git", "-C", a.code, "diff", "--quiet", "HEAD", "--", "phase2.2/" + f])
        r = subprocess.run(["git", "-C", a.code, "cat-file", "-e", "HEAD:phase2.2/" + f], stderr=subprocess.DEVNULL)
        if q.returncode or r.returncode:
            stop("%s is not committed, or differs from its committed version." % f)
    prep = load(os.path.join(a.code, "phase2.2/cs6_forms1_prepare.py"), "prep")
    # GF.3 the committed inputs are the ones the committed program forms from committed files
    F = prep.build(S6)
    tracked = subprocess.run(["git", "-C", a.res, "ls-files", S6REL + "/" + DIR], stdout=subprocess.PIPE, check=True).stdout.decode().split("\n")
    tracked = set(os.path.relpath(t, S6REL + "/" + DIR) for t in tracked if t)
    inputs = set(F)
    if not inputs <= tracked:
        stop("files of forms1/ that the program forms are not committed: %s" % ", ".join(sorted(inputs - tracked)[:5]))
    for p in sorted(F):
        if rd(os.path.join(D, p)) != F[p]:
            stop("forms1/%s is not the file the committed program forms from the committed files." % p)
    # GF.4 every output that exists echoes, as the input it ran, the committed input beside it
    runs = [r[0] for r in prep.RUNS]
    for r in runs:
        o = os.path.join(D, r, "job.out")
        if os.path.isfile(o):
            want = rd(os.path.join(D, r, "job.inp")).split("\n")
            got = echoed_input(rd(o))
            strip = lambda L: [x.rstrip() for x in L][:max([i + 1 for i, x in enumerate(L) if x.strip()] or [0])]
            if got and strip(got) != strip(want):
                stop("%s/job.out does not echo the committed job.inp beside it." % r)
    # D56's committed energy of the control with one electron removed
    E56 = None
    for row in csv.DictReader(open(os.path.join(S6, "d56.tsv")), delimiter="\t"):
        if row["run"] == "SV_C_bare_R":
            E56 = float(row["E_Eh"])
    if E56 is None:
        stop("d56.tsv holds no energy for SV_C_bare_R.")
    # the job logs
    exitc, vexit = {}, None
    for j in prep.JOBS:
        lp = os.path.join(D, "run_%s.pbs.out" % j)
        if os.path.isfile(lp):
            for l in rd(lp).split("\n"):
                m = re.match(r"^(\S+) orca_exit=(\d+) normal=(\d+)", l)
                if m:
                    exitc[m.group(1)] = int(m.group(2))
                m = re.match(r"^SV_mono_R vpot_exit=(\d+)", l)
                if m:
                    vexit = int(m.group(1))

    rows = []
    for r, job, cap, lv, form in prep.RUNS:
        e = EXPECT[r]
        d = os.path.join(D, r)
        o = os.path.join(d, "job.out")
        row = dict(run=r, job=job, level=lv, cap_s=cap, orca_exit=exitc.get(r, ""), status="", nbf="", nel="", npc="", ecp_centre="",
                   checks="", extra_file="", energy_check="", word="", form=form)
        if not os.path.isfile(o):
            row.update(status="NOT RUN", word="NOT SETTLED")
            rows.append(row)
            continue
        t = rd(o)
        began = ("ORCA LEAN-SCF" in t) or re.search(r"^Iteration\s+Energy \(Eh\)", t, re.M) is not None
        normal = "ORCA TERMINATED NORMALLY" in t
        if not began:
            row["status"] = "REFUSED"
        elif normal:
            row["status"] = "ACCEPTED AND ENDED NORMALLY"
        elif row["orca_exit"] == 124:
            row["status"] = "ACCEPTED, STOPPED AT ITS CAP"
        else:
            row["status"] = "ACCEPTED, DID NOT END NORMALLY"
        nbf, nel = num(r"^Number of basis functions\s+\.\.\.\s+(\d+)", t), num(r"Number of Electrons\s+NEL\s+\.\.\.\.\s+(\d+)", t)
        npc = num(r"Reading point charge file\s+\.\.\. ok \((\d+) point charges?\)", t)
        ecp = re.findall(r"^Atom\s+\d+\s*([A-Z][a-z]?>)\s+ECP group", t, re.M)
        row.update(nbf="" if nbf is None else nbf, nel="" if nel is None else nel, npc="" if npc is None else npc, ecp_centre=" ".join(ecp))
        bad = []
        if nel != e["nel"]:
            bad.append("electrons %s, expected %d" % (nel, e["nel"]))
        if e["npc"] is not None and npc != e["npc"]:
            bad.append("point charges read %s, expected %d" % (npc, e["npc"]))
        if e["npc"] is None and npc is not None:
            bad.append("point charges read where none were given")
        if e["ecp"] is not None and ecp != [e["ecp"]]:
            bad.append("core potential on %s, expected one on %s" % (ecp or "no centre", e["ecp"]))
        if e["ecp"] is None and ecp:
            bad.append("a core potential where none was given")
        if nbf is None:
            bad.append("no count of basis functions")
        elif e["extra"] and not nbf > NBF[lv]:
            bad.append("basis functions %d, expected more than the substrate's %d" % (nbf, NBF[lv]))
        elif not e["extra"] and nbf != NBF[lv]:
            bad.append("basis functions %d, expected the substrate's %d" % (nbf, NBF[lv]))
        need_end = r in ("SV_engrad_R", "SV_mono_R")
        if r == "SV_engrad_R" and normal:
            g = engrad(os.path.join(d, "job.engrad"))
            row["extra_file"] = "job.engrad: " + ("NOT READ" if g is None else "%d atoms, %d gradient values" % g) + "; job.pcgrad: " + pcgrad(os.path.join(d, "job.pcgrad"))
            if g != (24, 72):
                bad.append("job.engrad does not hold 24 atoms and 72 gradient values")
        if r == "SV_mono_R" and normal:
            E = re.findall(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", t)
            if E:
                dE = float(E[-1]) - E56
                row["energy_check"] = "%s D56's SV_C_bare_R (%+.1e Eh)" % ("REPRODUCES" if abs(dE) <= TOL_E else "DIFFERS FROM", dE)
                if abs(dE) > TOL_E:
                    bad.append("the energy differs from D56's committed control by %.1e Eh" % dE)
            else:
                bad.append("no final energy")
            v = vpot(d, F["points5.xyz"])
            row["extra_file"] = "vpot.out: %s (orca_vpot exit %s)" % (v, "not run" if vexit is None else vexit)
            if v != "WRITTEN, 5 POINTS":
                bad.append("orca_vpot gave no potential at the five points")
        row["checks"] = ("not read: no SCF began" if not began else ("HOLD" if not bad else "; ".join(bad)))
        if row["status"] == "REFUSED":
            row["word"] = "REFUSED"
        elif bad:
            row["word"] = "CHECK FAILED"
        elif normal:
            row["word"] = "WORKS"
        else:
            row["word"] = "NOT SETTLED" if need_end else "ACCEPTED, NOT FINISHED"
        rows.append(row)

    # the document
    import textwrap
    P = lambda s, ind=0: "\n".join(textwrap.wrap(" ".join(s.split()), W, initial_indent=" " * ind, subsequent_indent=" " * ind, break_long_words=False))
    head = subprocess.run(["git", "-C", a.res, "rev-parse", "HEAD"], stdout=subprocess.PIPE).stdout.decode().strip()
    L = ["STAGE 6 - FORMS TEST 1 - RESULT", "Does ORCA 6.0.1 on hpc1 accept each new form of input (AMENDMENT 13, A13.4 item 9)", "=" * W, ""]
    L.append(P("Read by cs6_forms1.py under STAGE6_FORMS_TEST1_20261125.txt, which was committed with the inputs before any job was submitted."))
    L.append("   results repository at %s" % head[:12])
    L.append("   reader               cs6_forms1.py  sha256 %s" % sha(os.path.join(a.code, "phase2.2/cs6_forms1.py")))
    L.append("   inputs formed by     cs6_forms1_prepare.py  sha256 %s" % sha(os.path.join(a.code, "phase2.2/cs6_forms1_prepare.py")))
    L += ["", P("THE GATES, as they held: GF.1 the pre-registration is the committed one; GF.2 the two programs are the committed ones; GF.3 every input of forms1/ is the file the committed program forms from the committed files; GF.4 every output that exists echoes the committed input beside it."), ""]
    L += ["=" * W, "RF.1  EACH RUN", "=" * W, ""]
    for row in rows:
        L.append("%-14s %s  %s" % (row["run"], row["level"], row["word"]))
        L.append(P("form: %s" % row["form"], 3))
        L.append(P("status: %s%s; electrons %s; basis functions %s; point charges read %s; core potential on %s" % (
            row["status"], "" if row["orca_exit"] == "" else " (orca exit %s)" % row["orca_exit"], row["nel"] or "-", row["nbf"] or "-",
            row["npc"] if row["npc"] != "" else "-", row["ecp_centre"] or "-"), 3))
        L.append(P("checks: %s" % (row["checks"] or "-"), 3))
        if row["extra_file"]:
            L.append(P(row["extra_file"], 3))
        if row["energy_check"]:
            L.append(P(row["energy_check"], 3))
        L.append("")
    cnt = {}
    for row in rows:
        cnt[row["word"]] = cnt.get(row["word"], 0) + 1
    L += ["=" * W, "RF.2  WHAT EACH WORD SENDS THE NEXT STEP TO (pre-registration, section 4; its words govern)", "=" * W, ""]
    L.append(P("WORKS: the form is used as it was tested. ACCEPTED, NOT FINISHED: ORCA took the input and the checks that come before the SCF hold; the later package sets its cap with that run's time in view and says so. REFUSED, CHECK FAILED or NOT SETTLED: the form is not used in any later package until a dated note says how it is mended and a run shows it."))
    L += ["", P("This test names no design, forms no change of a barrier and adopts nothing.")]
    L += ["", "Runs: %d. %s." % (len(rows), "; ".join("%s: %d" % (k, cnt[k]) for k in ("WORKS", "ACCEPTED, NOT FINISHED", "CHECK FAILED", "REFUSED", "NOT SETTLED") if k in cnt))]
    txt = "\n".join(L) + "\n"
    with open(a.result, "w") as f:
        f.write(txt)
    cols = ["run", "job", "level", "cap_s", "orca_exit", "status", "nel", "nbf", "npc", "ecp_centre", "checks", "extra_file", "energy_check", "word", "form"]
    with open(a.tsv, "w") as f:
        f.write("\t".join(cols) + "\n")
        for row in rows:
            f.write("\t".join(str(row[c]) for c in cols) + "\n")
    print(txt, end="")


if __name__ == "__main__":
    main()
