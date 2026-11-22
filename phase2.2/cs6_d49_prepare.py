#!/usr/bin/env python3
"""cs6_d49_prepare.py - Stage 6 DIAGNOSTIC 49: the calibration of DIAGNOSTIC 48, run again so that ORCA ends normally.

D48's part C is one electron, the def2-TZVPD functions of a hydrogen atom (a ghost atom, no nucleus) and a bare point
charge q 0.01 A from their centre, at Hartree-Fock (HF) and with the functional of D48's runs (CAM), at ten sizes. All
twenty runs converged and wrote a final energy, and ORCA then stopped in its property step, where the dipole moment
needs a nucleus and there is none. D48's committed reader counts a run that does not end normally NOT MEASURED.

D49 runs each of the twenty in two variants, written from the committed inputs of D48:
  N  the committed D48 input with three lines added that switch the dipole moment off. Same electron, same functions,
     same charge, same place.
  Z  no point charge and no ghost atom: a hydrogen atom whose nuclear charge is set to q ("Z = q" on its coordinate
     line), neutral by ORCA's count, one electron. The charge then sits AT the centre of the functions and not 0.01 A
     from it, and ORCA has a nucleus for its property step.
Neither variant has been run in ORCA before in this record.

Writes, under the output directory: one folder for each run (job.inp, and pc.pc for variant N), RUNS.tsv,
SOURCES.sha256, INPUT_CHECKSUMS.sha256 and the job script run_K1.pbs. It reads committed files only and edits none.

usage: python3 cs6_d49_prepare.py <results repository> <output directory> [<directory of d49 on hpc1>]
"""
import hashlib
import os
import re
import sys

S6 = "phase2.2/cs_stage6"
D48 = "d48"
PC_LINE = '%pointcharges "pc.pc"'
GHOST_LINE = "H: 0.010000 0.000000 0.000000"
ELPROP = ["%elprop", "  Dipole false", "end"]
VARIANTS = ("N", "Z")
LEVELS = ("CAM", "HF")
NRUN48 = 20
JOB, PBSNAME = "K1", "cm_d49_k1"
TIMEOUT = 900
HPC = "/home/18660916/system_development/phase2.2/cs_stage6/d49"
OTHER = ["STAGE6_DIAGNOSTIC48_20261121.txt", "STAGE6_DIAGNOSTIC48_RESULT_20261121.txt", "d48.tsv", "d48/RUNS.tsv", "d48/CONDITIONS.tsv",
         "d48/SOURCES.sha256", "d48/INPUT_CHECKSUMS.sha256", "d48/run_E1.pbs"]


def stop(m):
    sys.stderr.write("STOP  %s\n" % m)
    raise SystemExit(1)


def sha_text(t):
    return hashlib.sha256(t.encode()).hexdigest()


def sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


def one_runs(s6):
    """D48's twenty calibration runs, in D48's order: (run, level, size as text, q)"""
    ls = open(os.path.join(s6, D48, "RUNS.tsv")).read().splitlines()
    h = ls[0].split("\t")
    out = []
    for l in ls[1:]:
        r = dict(zip(h, l.split("\t")))
        if r["kind"] != "one":
            continue
        lv = r["level"].replace("ONE_", "")
        if lv not in LEVELS or not re.fullmatch(r"\d\.\d\d", r["size_e"]) or r["run"] != "ONE_%s_p%03d" % (lv, round(100 * float(r["size_e"]))):
            stop("d48/RUNS.tsv holds a calibration run this script does not know: %s" % l)
        out.append((r["run"], lv, r["size_e"], float(r["size_e"])))
    if len(out) != NRUN48 or len(set(x[0] for x in out)) != NRUN48:
        stop("d48/RUNS.tsv holds %d calibration runs, not %d." % (len(out), NRUN48))
    return out


def inputs(s6, run, q):
    """the two inputs written from D48's committed input of one calibration run, and the charge file of variant N"""
    t = open(os.path.join(s6, D48, run, "job.inp")).read()
    pc = open(os.path.join(s6, D48, run, "pc.pc")).read()
    ls = t.split("\n")
    if ls.count(PC_LINE) != 1 or ls.count(GHOST_LINE) != 1 or ls.count("* xyz -1 2") != 1 or not ls[0].startswith("! ") or "%elprop" in t:
        stop("d48/%s/job.inp is not the input this script was written for." % run)
    w = pc.split()
    if len(w) != 5 or w[0] != "1" or abs(float(w[1]) - q) > 1e-9 or any(abs(float(v)) > 1e-12 for v in w[2:]):
        stop("d48/%s/pc.pc is not one charge of %.2f e at the origin." % (run, q))
    k = ls.index(PC_LINE)
    n = "\n".join(ls[:k] + ELPROP + ls[k:])
    body = [l for l in ls[1:k]]                                   # %maxcore and %scf of D48's input, unchanged
    z = "\n".join([ls[0]] + body + ["* xyz 0 2", "H 0.000000 0.000000 0.000000 Z = %.6f" % q, "*", ""])
    return {"N": n, "Z": z}, pc


def runs_of(s6):
    out = []
    one = one_runs(s6)
    for v in VARIANTS:
        for run, lv, size, q in one:
            out.append(dict(run="%s_%s" % (run, v), job=JOB, level=lv, variant=v, size=size, q=q, d48=run))
    return out


def runs_text(runs):
    L = ["run\tjob\tlevel\tvariant\tsize_e\trun_of_D48"]
    for r in runs:
        L.append("\t".join([r["run"], r["job"], r["level"], r["variant"], r["size"], r["d48"]]))
    return "\n".join(L) + "\n"


def source_files(s6):
    out = list(OTHER)
    for run, lv, size, q in one_runs(s6):
        out += ["%s/%s/job.inp" % (D48, run), "%s/%s/pc.pc" % (D48, run), "%s/%s/job.out" % (D48, run)]
    return out


def pbs_text(s6, runs, hpc):
    """the job script: D48's committed run_E1.pbs with its names and its runs changed, and nothing else"""
    t = open(os.path.join(s6, D48, "run_E1.pbs")).read()
    ls = t.split("\n")
    i0 = next(i for i, l in enumerate(ls) if l.startswith("run ONE_"))
    i1 = next(i for i, l in enumerate(ls) if l.startswith('echo "end='))
    head, tail = ls[:i0], ls[i1:]
    if [l for l in ls[i0:i1] if not l.startswith("run ONE_")] or len(ls[i0:i1]) != NRUN48 or len(tail) != 4 or tail[3] != "" \
            or not tail[1].startswith("n=0; for r in ONE_") or tail[2] != 'echo "finished: $n of %d"' % NRUN48:
        stop("d48/run_E1.pbs is not the job script this script was written for.")
    rep = [("#PBS -N cm_d48_e1", "#PBS -N %s" % PBSNAME),
           ("#PBS -o /home/18660916/system_development/phase2.2/cs_stage6/d48/run_E1.pbs.out", "#PBS -o %s/run_%s.pbs.out" % (hpc, JOB)),
           ("# Stage 6 DIAGNOSTIC 48: is the estimate S of DIAGNOSTIC 47 a rule?", "# Stage 6 DIAGNOSTIC 49: the calibration of DIAGNOSTIC 48, run again so that ORCA ends normally."),
           ("# This job is one of 7. It makes 20 runs in the order below; a run that fails does not stop the ones after it.",
            "# This is the only job. It makes %d runs in the order below; a run that fails does not stop the ones after it." % len(runs)),
           ("W=/home/18660916/system_development/phase2.2/cs_stage6/d48", "W=%s" % hpc),
           ('if ! mkdir .running_E1 2>/dev/null; then echo "FAIL: another job holds .running_E1 in $W"; exit 1; fi',
            'if ! mkdir .running_%s 2>/dev/null; then echo "FAIL: another job holds .running_%s in $W"; exit 1; fi' % (JOB, JOB)),
           ("trap 'rmdir \"$W/.running_E1\" 2>/dev/null' EXIT", "trap 'rmdir \"$W/.running_%s\" 2>/dev/null' EXIT" % JOB)]
    for a, b in rep:
        if head.count(a) != 1:
            stop("d48/run_E1.pbs does not hold the line: %s" % a)
        head[head.index(a)] = b
    if any(("d48" in l or "E1" in l or "of 7" in l) for l in head):
        stop("a line of d48/run_E1.pbs still names D48 after the names were changed.")
    names = [r["run"] for r in runs]
    mid = ["run %s" % n for n in names]
    end = [tail[0], "n=0; for r in %s; do grep -q 'ORCA TERMINATED NORMALLY' \"$W/$r/job.out\" 2>/dev/null && n=$((n+1)); done" % " ".join(names),
           'echo "finished: $n of %d"' % len(names), ""]
    if tail[1] != "n=0; for r in %s; do grep -q 'ORCA TERMINATED NORMALLY' \"$W/$r/job.out\" 2>/dev/null && n=$((n+1)); done" % " ".join(l[4:] for l in ls[i0:i1]):
        stop("the last loop of d48/run_E1.pbs is not the one this script was written for.")
    return "\n".join(head + mid + end)


def build(s6, hpc=HPC):
    """every file of d49, as {path under d49: text}; the checksum list last"""
    runs = runs_of(s6)
    F = {}
    for r in runs:
        ins, pc = inputs(s6, r["d48"], r["q"])
        F["%s/job.inp" % r["run"]] = ins[r["variant"]]
        if r["variant"] == "N":
            F["%s/pc.pc" % r["run"]] = pc
    F["RUNS.tsv"] = runs_text(runs)
    F["SOURCES.sha256"] = "".join("%s  %s\n" % (sha(os.path.join(s6, f)), f) for f in source_files(s6))
    F["run_%s.pbs" % JOB] = pbs_text(s6, runs, hpc)
    F["INPUT_CHECKSUMS.sha256"] = "".join("%s  %s\n" % (sha_text(F[k]), k) for k in sorted(F))
    return F, runs


def main():
    if len(sys.argv) not in (3, 4):
        stop("usage: python3 cs6_d49_prepare.py <results repository> <output directory> [<directory of d49 on hpc1>]")
    s6 = os.path.join(sys.argv[1], S6)
    out = sys.argv[2]
    hpc = sys.argv[3] if len(sys.argv) == 4 else HPC
    for f in source_files(s6):
        if not os.path.isfile(os.path.join(s6, f)):
            stop("the committed file %s is missing." % f)
    if os.path.exists(out):
        stop("%s exists already; this script overwrites nothing." % out)
    F, runs = build(s6, hpc)
    for k in sorted(F):
        p = os.path.join(out, k)
        os.makedirs(os.path.dirname(p), exist_ok=True)
        with open(p, "w") as f:
            f.write(F[k])
        if k.endswith(".pbs"):
            os.chmod(p, 0o755)
    print("%d runs (%d systems in %d variants), %d files checksummed, 1 job" % (len(runs), NRUN48, len(VARIANTS), len(F) - 1))


if __name__ == "__main__":
    main()
