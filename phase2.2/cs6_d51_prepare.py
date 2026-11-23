#!/usr/bin/env python3
"""cs6_d51_prepare.py - Stage 6 DIAGNOSTIC 51: the calibration of D50 done again, with the integration grid measured.

D50 held the functional's value for one electron at a bare charge against the exact one and its rule gave NONE. Read
afterwards, its committed outputs show that ORCA's default integration grid does not count the whole of so diffuse
an electron on a lone ghost atom (0.92 of it at 0.25 e), and a calculation in a second program shows that the
functional itself lies 0.05 to 0.06 eV from exact at 0.20 to 0.60 e. D51 measures the two apart, in ORCA.

One electron, a bare charge q, the functions of D50 (sets X2 and X3), ten sizes:

  HF     D50's committed Hartree-Fock inputs, byte for byte.
  g2     D50's committed inputs with the functional and auxiliary treatment A, byte for byte (ORCA's default grid).
  g3     the same with DefGrid3.
  i08, i12, i20, i40   the same as g2 with IntAcc 8, 12, 20, 40: more radial points in the grid of the functional.
  c40    i40 with the three grids of the exchange (COSX) raised to IntAcc 40 as well.
  n2     the functional with no RI and no COSX (NORI), default grid.
  n40    the same with IntAcc 40.

Every input is D50's committed input with at most one line changed and one block added; this script forms them
through the committed cs6_d50_prepare.py, which it loads by path and whose sha256 it checks, and refuses unless
what it forms without a change is the committed input of D50. It reads committed files only and edits none.

usage: python3 cs6_d51_prepare.py <results repository> <code repository> <output directory> [<directory of d51 on hpc1>]
"""
import hashlib
import importlib.util
import os
import sys

S6 = "phase2.2/cs_stage6"
P50 = "cs6_d50_prepare.py"
P50_SHA = "bdbfce5202644a9e2d31c557377d5139007eccc099a38b00225150af2a4b3f82"
EXPF = "cs6_d51_expected.py"
EXP_SHA = "69340f4b7d373156a34e29bdaccb49b1b0227252a16fac08db8bba5d309c91a4"
HPC = "/home/18660916/system_development/phase2.2/cs_stage6/d51"
AUX = "A"                                                  # the auxiliary treatment of the ladder: D50's first, ORCA's AutoAux
RI_WORDS, NORI_WORD = "def2/J RIJCOSX", "NORI"             # what the line of the functional loses and gains in rungs n2 and n40
# rung: (kind, what is added to the "!" line, IntAcc of the functional's grid or None, IntAcc of the three COSX grids or None)
RUNGS = (("g2", "ri", "", None, None), ("g3", "ri", "DefGrid3", None, None), ("i08", "ri", "", 8.0, None), ("i12", "ri", "", 12.0, None),
         ("i20", "ri", "", 20.0, None), ("i40", "ri", "", 40.0, None), ("c40", "ri", "", 40.0, 40.0), ("n2", "nori", "", None, None),
         ("n40", "nori", "", 40.0, None))
LADDER = ("g2", "g3", "i08", "i12", "i20", "i40")           # the grid of the functional, from ORCA's default upward
INTACC = {"g2": 4.388, "g3": 4.959, "i08": 8.0, "i12": 12.0, "i20": 20.0, "i40": 40.0}   # g2 and g3: what ORCA prints for them
TOL_ADEQ = 0.05                                            # eV, D50's: within this of exact (HF), or of the larger set (the functional)
TOL_GRID = 0.005                                           # eV, a tenth of it: the two highest rungs measured must agree to this
TOL_EH = 1e-6                                              # hartree: HF and g2 against D50's committed energies
FROM = 25                                                  # hundredths of e: the rule asks for every size from this one up
JOBS = [("K1", "cm_d51_k1")]
EXTRA = ["STAGE6_DIAGNOSTIC50_20261122.txt", "STAGE6_DIAGNOSTIC50_RESULT_20261123.txt", "d50_calibration.tsv", "d50.tsv", "d48/run_E1.pbs"]


class Refused(Exception):
    pass


def stop(m):
    raise Refused(m)


def sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


def sha_text(t):
    return hashlib.sha256(t.encode()).hexdigest()


def load_by_path(folder, fn, want, name):
    p = os.path.join(folder, fn)
    if not os.path.isfile(p):
        stop("%s is not in %s." % (fn, folder))
    if sha(p) != want:
        stop("%s is not the committed file this script was written for." % fn)
    spec = importlib.util.spec_from_file_location(name, p)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def method_block(xc, cosx):
    L = ["%method", "  IntAcc %.1f" % xc]
    if cosx is not None:
        L.append("  IntAccX %.1f, %.1f, %.1f" % (cosx, cosx, cosx))
    return L + ["end"]


def hf_name(name, k):
    return "ONE_HF_%s_p%03d" % (name, k)


def cam_name(name, rung, k):
    return "ONE_CAM_%s_%s_p%03d" % (name, rung, k)


def cam_input(p50, p48, name, rung):
    """D50's committed input of the functional with treatment A, changed as the rung says"""
    kind, word, xc, cosx = dict((r[0], r[1:]) for r in RUNGS)[rung]
    if kind == "ri":
        ls = p50.one_input(p48, "CAM", name, AUX).split("\n")
    else:
        plain = p50.one_plain(p48, "CAM").split("\n")
        if plain[0].count(RI_WORDS) != 1:
            stop("the line of the functional in D49's variant N does not hold '%s' once." % RI_WORDS)
        plain[0] = plain[0].replace(RI_WORDS, NORI_WORD)
        ls = p50.with_block("\n".join(plain), name, None, 1).split("\n")
    if not ls[0].startswith("! "):
        stop("an input does not begin with its line of keywords.")
    if word:
        ls[0] = ls[0] + " " + word
    if xc is not None:
        if ls.count(p48.PC_LINE) != 1:
            stop("an input does not hold the line of the point charges once.")
        k = ls.index(p48.PC_LINE)
        ls = ls[:k] + method_block(xc, cosx) + ls[k:]
    return "\n".join(ls)


def check_d50(s6, p50, p48):
    """what is formed here without a change is what D50 committed and ran"""
    for name in sorted(p50.ADD):
        for k in p50.SIZES:
            for run, text in ((p50.one_name("HF", name, None, k), p50.one_input(p48, "HF", name, None)),
                              (p50.one_name("CAM", name, AUX, k), cam_input(p50, p48, name, "g2"))):
                d = os.path.join(s6, "d50", run)
                for fn, want in (("job.inp", text), ("pc.pc", p48.one_pc(k))):
                    p = os.path.join(d, fn)
                    if not os.path.isfile(p) or open(p).read() != want:
                        stop("d50/%s/%s is not the file formed here from the committed scripts." % (run, fn))
                if not os.path.isfile(os.path.join(d, "job.out")):
                    stop("d50/%s/job.out, D50's committed output, is missing." % run)


def runs_of(p50):
    out = []
    for name in sorted(p50.ADD):
        for k in p50.SIZES:
            out.append(dict(run=hf_name(name, k), job="K1", level="ONE_HF", set=name, rung="", size=k, d50=p50.one_name("HF", name, None, k)))
        for rung in [r[0] for r in RUNGS]:
            for k in p50.SIZES:
                out.append(dict(run=cam_name(name, rung, k), job="K1", level="ONE_CAM", set=name, rung=rung, size=k,
                                d50=p50.one_name("CAM", name, AUX, k) if rung == "g2" else ""))
    if len(set(r["run"] for r in out)) != len(out):
        stop("two runs carry one name.")
    return out


def input_text(p50, p48, r):
    if r["level"] == "ONE_HF":
        return p50.one_input(p48, "HF", r["set"], None)
    return cam_input(p50, p48, r["set"], r["rung"])


def runs_text(runs):
    L = ["run\tjob\tlevel\tset\trung\tsize_e\tthe_same_input_in_d50"]
    for r in runs:
        L.append("\t".join([r["run"], r["job"], r["level"], r["set"], r["rung"], "%.2f" % (r["size"] / 100.0), r["d50"]]))
    return "\n".join(L) + "\n"


def rungs_text():
    L = ["rung\tintegrals\tadded_to_the_line_of_keywords\tIntAcc_of_the_grid_of_the_functional\tIntAcc_of_the_three_COSX_grids"]
    for rung, kind, word, xc, cosx in RUNGS:
        L.append("\t".join([rung, "RI-J with AutoAux on the ghost atom, COSX" if kind == "ri" else "no RI, no COSX (NORI)", word or "-",
                            ("%.3f, ORCA's own for this keyword" % INTACC[rung]) if xc is None and rung in INTACC else ("default" if xc is None else "%.1f" % xc),
                            "default" if cosx is None else "%.1f" % cosx]))
    return "\n".join(L) + "\n"


def expected_text(p50, p48, ex):
    L = ["set\tsize_e\texact_eV\tB_HF_full_eV\tB_CAM_full_eV\tB_HF_limit_eV\tB_CAM_limit_eV\t" +
         "\t".join("count_%s_model%d\tB_CAM_%s_model%d_eV" % (rung, n0, rung, n0) for rung in LADDER for n0 in ex.MODELS)]
    for name in sorted(p50.ADD):
        for k in p50.SIZES:
            q = k / 100.0
            row = [name, "%.2f" % q, "%.6f" % (p48.RYD * q * q), "%.6f" % ex.FULL[(name, k)][0], "%.6f" % ex.FULL[(name, k)][1],
                   "%.6f" % ex.LIMIT[k][0], "%.6f" % ex.LIMIT[k][1]]
            for rung in LADDER:
                for n0 in ex.MODELS:
                    row += ["%.4f" % ex.MODEL[(n0, name, rung, k)][0], "%.4f" % ex.MODEL[(n0, name, rung, k)][1]]
            L.append("\t".join(row))
    return "\n".join(L) + "\n"


def sources(p50):
    out = list(EXTRA)
    for name in sorted(p50.ADD):
        for k in p50.SIZES:
            for run in (p50.one_name("HF", name, None, k), p50.one_name("CAM", name, AUX, k)):
                out += ["d50/%s/job.inp" % run, "d50/%s/pc.pc" % run, "d50/%s/job.out" % run]
    return out


def pbs_text(s6, runs, hpc):
    """the job script: D48's committed run_E1.pbs (one electron, each run stopped after 900 s), with its names and its
    runs changed and nothing else, as D50's run_K1.pbs was made"""
    jk, pbs = JOBS[0]
    old = "E1"
    mine = [r["run"] for r in runs]
    ls = open(os.path.join(s6, "d48", "run_%s.pbs" % old)).read().split("\n")
    i0 = next(i for i, l in enumerate(ls) if l.startswith("run "))
    i1 = next(i for i, l in enumerate(ls) if l.startswith('echo "end='))
    head, mid, tail = ls[:i0], ls[i0:i1], ls[i1:]
    names = [l[4:] for l in mid]
    loop = "n=0; for r in %s; do grep -q 'ORCA TERMINATED NORMALLY' \"$W/$r/job.out\" 2>/dev/null && n=$((n+1)); done"
    if any(not l.startswith("run ") for l in mid) or len(tail) != 4 or tail[3] != "" or tail[1] != loop % " ".join(names) \
            or tail[2] != 'echo "finished: $n of %d"' % len(names):
        stop("d48/run_%s.pbs is not the job script this script was written for." % old)
    d48 = "/home/18660916/system_development/phase2.2/cs_stage6/d48"
    rep = [("#PBS -N cm_d48_%s" % old.lower(), "#PBS -N %s" % pbs),
           ("#PBS -o %s/run_%s.pbs.out" % (d48, old), "#PBS -o %s/run_%s.pbs.out" % (hpc, jk)),
           ("# Stage 6 DIAGNOSTIC 48: is the estimate S of DIAGNOSTIC 47 a rule?",
            "# Stage 6 DIAGNOSTIC 51: the calibration of D50 done again, with the integration grid measured."),
           ("# This job is one of 7. It makes %d runs in the order below; a run that fails does not stop the ones after it." % len(names),
            "# This job is the only one. It makes %d runs in the order below; a run that fails does not stop the ones after it." % len(mine)),
           ("W=%s" % d48, "W=%s" % hpc),
           ('if ! mkdir .running_%s 2>/dev/null; then echo "FAIL: another job holds .running_%s in $W"; exit 1; fi' % (old, old),
            'if ! mkdir .running_%s 2>/dev/null; then echo "FAIL: another job holds .running_%s in $W"; exit 1; fi' % (jk, jk)),
           ("trap 'rmdir \"$W/.running_%s\" 2>/dev/null' EXIT" % old, "trap 'rmdir \"$W/.running_%s\" 2>/dev/null' EXIT" % jk)]
    for a, b in rep:
        if head.count(a) != 1:
            stop("d48/run_%s.pbs does not hold the line: %s" % (old, a))
        head[head.index(a)] = b
    if any(("d48" in l or old in l or "of 7" in l) for l in head):
        stop("a line of d48/run_%s.pbs still names D48 after the names were changed." % old)
    return "\n".join(head + ["run %s" % n for n in mine] + [tail[0], loop % " ".join(mine), 'echo "finished: $n of %d"' % len(mine), ""])


def build(s6, code, hpc=HPC):
    """every file of d51, as {path under d51: text}; the checksum list last"""
    p50 = load_by_path(os.path.join(code, "phase2.2"), P50, P50_SHA, "d51_cs6_d50_prepare")      # committed: from the code repository given
    ex = load_by_path(os.path.dirname(os.path.abspath(__file__)), EXPF, EXP_SHA, "d51_cs6_d51_expected")   # committed beside this script
    try:
        p48 = p50.load48(code)
        check_d50(s6, p50, p48)
    except p50.Refused as e:
        stop("cs6_d50_prepare.py: %s" % e)
    runs = runs_of(p50)
    F = {}
    for r in runs:
        F["%s/job.inp" % r["run"]] = input_text(p50, p48, r)
        F["%s/pc.pc" % r["run"]] = p48.one_pc(r["size"])
    F["RUNS.tsv"] = runs_text(runs)
    F["RUNGS.tsv"] = rungs_text()
    F["EXPECTED.tsv"] = expected_text(p50, p48, ex)
    F["SOURCES.sha256"] = "".join("%s  %s\n" % (sha(os.path.join(s6, f)), f) for f in sources(p50))
    F["run_K1.pbs"] = pbs_text(s6, runs, hpc)
    F["INPUT_CHECKSUMS.sha256"] = "".join("%s  %s\n" % (sha_text(F[k]), k) for k in sorted(F))
    return F, runs, p50, p48, ex


def main():
    if len(sys.argv) not in (4, 5):
        sys.stderr.write("usage: python3 cs6_d51_prepare.py <results repository> <code repository> <output directory> [<directory of d51 on hpc1>]\n")
        return 1
    s6 = os.path.join(sys.argv[1], S6)
    code, out = sys.argv[2], sys.argv[3]
    hpc = sys.argv[4] if len(sys.argv) == 5 else HPC
    if os.path.exists(out):
        sys.stderr.write("STOP  %s exists already; this script overwrites nothing.\n" % out)
        return 1
    try:
        F, runs, p50, p48, ex = build(s6, code, hpc)
    except Refused as e:
        sys.stderr.write("STOP  %s\n" % e)
        return 1
    for k in sorted(F):
        p = os.path.join(out, k)
        os.makedirs(os.path.dirname(p), exist_ok=True)
        with open(p, "w") as f:
            f.write(F[k])
        if k.endswith(".pbs"):
            os.chmod(p, 0o755)
    print("%d runs of one electron (%d Hartree-Fock, %d with the functional in %d rungs), %d files checksummed, %d job" % (
        len(runs), sum(1 for r in runs if r["level"] == "ONE_HF"), sum(1 for r in runs if r["level"] == "ONE_CAM"), len(RUNGS), len(F) - 1, len(JOBS)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
