#!/usr/bin/env python3
"""Stage 6, DIAGNOSTIC 54: the design at the second level, from two starts. THE INPUTS.

Five runs of the dianion of frame 42436 with the eleven charges of B_n11_q0.25 and the functions of D50 at the six
positive sites, each started from orbitals that a run of D53 left on hpc1:

  O  the run of D53 at wB97X-V again, started from its own converged orbitals (R and TS);
  S  the run at R started from the solution D53 found at CAM-B3LYP-D3BJ, in which half an electron sits on the
     functions: at wB97X-V (the question), at CAM-B3LYP-D3BJ (the control of the way of starting) and at Hartree-Fock.

Every input is the committed input of its run of D53 with two changes and no other: the first line asks for the
orbitals to be read (and, for S, for separate orbitals of the two spins, since the start has them), and one line names
the file, start.gbw. The file is not in the repository; its sha256 is, in D53's committed HPC_MANIFEST.sha256, and
stage_orbitals.sh, written here, copies it into place on hpc1 and checks it.

usage: python3 cs6_d54_prepare.py <results repository> <code repository> <output directory> [<directory of d54 on hpc1>]
"""
import hashlib
import os
import sys

S6 = "phase2.2/cs_stage6"
HPC = "/home/18660916/system_development/phase2.2/cs_stage6/d54"
D53 = "d53"
MANIFEST = "d53/HPC_MANIFEST.sha256"
PBS53 = "d53/run_K3.pbs"
PBS53_HPC = "/home/18660916/system_development/phase2.2/cs_stage6/d53"
NB = 784
CAP53 = 10200       # the cap of one run in D53's committed job script, s
CAP = 18000         # the cap of one run here, s: D53's longest run, the one in which the electrons were taken, took 7509 s,
                    # and a run stopped by its cap is not made again
LIMIT_H = 24
LEVEL_NAME = {"CAM": "CAM-B3LYP-D3BJ", "WXV": "wB97X-V", "HF": "Hartree-Fock"}
# the first line of a run that starts from the capture solution: separate orbitals for the two spins, spelled as D52's committed inputs spell them
# ("! UKS wB97X-V ..." for a functional; "! UHF ..." for Hartree-Fock, where UHF takes the place of HF)
FIRST_S = {"CAM": ("! ", "! UKS "), "WXV": ("! ", "! UKS "), "HF": ("! HF ", "! UHF ")}
START = "start.gbw"
# run, job, level, kind, state, the run of D53 whose input it is, the run of D53 whose orbitals it starts from
RUNS = (("DI7_WXV_O_d11_R", "K1", "WXV", "O", "R", "DI6_WXV_F_d11_R", "DI6_WXV_F_d11_R"),
        ("DI7_WXV_O_d11_TS", "K1", "WXV", "O", "TS", "DI6_WXV_F_d11_TS", "DI6_WXV_F_d11_TS"),
        ("DI7_WXV_S_d11_R", "K1", "WXV", "S", "R", "DI6_WXV_F_d11_R", "DI6_CAM_F_d11_R"),
        ("DI7_CAM_S_d11_R", "K2", "CAM", "S", "R", "DI6_CAM_F_d11_R", "DI6_CAM_F_d11_R"),
        ("DI7_HF_S_d11_R", "K2", "HF", "S", "R", "DI6_HF_F_d11_R", "DI6_CAM_F_d11_R"))
JOBS = (("K1", "cm_d54_k1", "the second level: its own start at R and at TS, then the capture solution as start at R"),
        ("K2", "cm_d54_k2", "the capture solution as start at R: the working functional (the control), then Hartree-Fock"))
# the runs of D53 each run of D54 is read against: its control with functions, the plain control, the partner without functions
REFS = {"WXV": {"R": ("DI6_WXV_GF_R", "DI6_WXV_GN_R", "DI6_WXV_N_d11_R"), "TS": ("DI6_WXV_GF_TS", "DI6_WXV_GN_TS", "DI6_WXV_N_d11_TS")},
        "CAM": {"R": ("DI6_CAM_GF_R", "DI6_CAM_GN_R", "DI6_CAM_N_d11_R")},
        "HF": {"R": ("DI6_HF_GF_R", "DI6_HF_GN_R", None)}}


class Refused(Exception):
    pass


def stop(m):
    raise Refused(m)


def sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


def sha_text(t):
    return hashlib.sha256(t.encode()).hexdigest()


def runs_of():
    return [dict(run=r[0], job=r[1], level=r[2], kind=r[3], state=r[4], source=r[5], start=r[6], nb=NB) for r in RUNS]


def manifest(s6):
    """{file under d53: sha256}, from the manifest hpc1 wrote over D53's files and D53's harvest committed"""
    m = {}
    for l in open(os.path.join(s6, MANIFEST)).read().split("\n"):
        if l.strip():
            h, f = l.split(None, 1)
            m[f.strip()] = h
    return m


def input_of(s6, r):
    """the committed input of the run of D53, with the first line changed and one line added"""
    ls = open(os.path.join(s6, D53, r["source"], "job.inp")).read().split("\n")
    if not ls[0].startswith("! ") or "MORead" in ls[0] or any("moinp" in l.lower() for l in ls) or " UKS" in ls[0] or " UHF" in ls[0]:
        stop("%s/%s/job.inp is not the input this script was written for." % (D53, r["source"]))
    if r["kind"] == "O":
        first = ls[0] + " MORead"
    else:
        a, b = FIRST_S[r["level"]]
        if not ls[0].startswith(a):
            stop("the first line of %s/%s/job.inp is not the one this script was written for." % (D53, r["source"]))
        first = b + ls[0][len(a):] + " MORead"
    return "\n".join([first, '%%moinp "%s"' % START] + ls[1:])


def start_table(s6):
    m = manifest(s6)
    rows = []
    for r in runs_of():
        f = "%s/job.gbw" % r["start"]
        if f not in m:
            stop("%s does not list %s." % (MANIFEST, f))
        rows.append((r["run"], "%s/%s" % (D53, f), m[f]))
    return rows


def stage_text(s6):
    rows = start_table(s6)
    ls = ["#!/bin/bash",
          "# Stage 6 DIAGNOSTIC 54: puts the orbitals each run starts from into its folder as %s." % START,
          "# Run on hpc1 in the folder d54, which lies beside d53, before any job is submitted (block P4 runs it).",
          "# Each file is a copy of an orbital file that a run of D53 left on hpc1, d53/<run>/job.gbw, and is",
          "# checked against the sha256 that D53's committed HPC_MANIFEST.sha256 holds for that file, before it is",
          "# copied and again after.",
          "# It copies nothing over a file that is there already and changes nothing in d53.",
          'cd "$(dirname "$0")" || exit 1']
    for src in sorted(set(s for _, s, _ in rows)):
        ls.append('[ -f "../%s" ] || { echo "STOP: %s is not on hpc1. Nothing was copied."; exit 1; }' % (src, src))
    for src, h in sorted(set((s, h) for _, s, h in rows)):
        ls.append('[ "$(sha256sum "../%s" | cut -c1-64)" = %s ] || { echo "STOP: %s on hpc1 is not the file the committed manifest of D53 lists. Nothing was copied."; exit 1; }' % (src, h, src))
    for run, src, h in rows:
        ls.append('[ -e "%s/%s" ] || cp "../%s" "%s/%s" || { echo "STOP: the copy into %s failed."; exit 1; }' % (run, START, src, run, START, run))
    ls.append('sha256sum -c --quiet START_ORBITALS.sha256 || { echo "STOP: a start file is not the file the committed manifest of D53 lists. Submit nothing."; exit 1; }')
    ls.append('echo "start orbitals: %d files in place, each the file the committed manifest of D53 lists"' % len(rows))
    return "\n".join(ls) + "\n"


def pbs_text(s6, job, hpc):
    """the job script: D53's committed run_K3.pbs with its names, its runs and the cap of one run changed"""
    jk, pbs, what = job
    mine = [r["run"] for r in runs_of() if r["job"] == jk]
    ls = open(os.path.join(s6, PBS53)).read().split("\n")
    i0 = next(i for i, l in enumerate(ls) if l.startswith("run "))
    i1 = next(i for i, l in enumerate(ls) if l.startswith('echo "end='))
    head, mid, tail = ls[:i0], ls[i0:i1], ls[i1:]
    names = [l[4:] for l in mid]
    loop = "n=0; for r in %s; do grep -q 'ORCA TERMINATED NORMALLY' \"$W/$r/job.out\" 2>/dev/null && n=$((n+1)); done"
    if any(not l.startswith("run ") for l in mid) or len(tail) != 4 or tail[3] != "" or tail[1] != loop % " ".join(names) \
            or tail[2] != 'echo "finished: $n of %d"' % len(names):
        stop("%s is not the job script this script was written for." % PBS53)
    rep = [("#PBS -N cm_d53_k3", "#PBS -N %s" % pbs),
           ("#PBS -o %s/run_K3.pbs.out" % PBS53_HPC, "#PBS -o %s/run_%s.pbs.out" % (hpc, jk)),
           ("# Stage 6 DIAGNOSTIC 53: charges of 0.25 e together, on the dianion, at two levels.",
            "# Stage 6 DIAGNOSTIC 54: the design at the second level, from two starts."),
           ("# Job K3 of 9: the second level, with the functions, at R. It makes %d runs in the order below; a run that fails does not stop the ones after it." % len(names),
            "# Job %s of %d: %s. It makes %d runs in the order below; a run that fails does not stop the ones after it." % (jk, len(JOBS), what, len(mine))),
           ("W=%s" % PBS53_HPC, "W=%s" % hpc),
           ('if ! mkdir .running_K3 2>/dev/null; then echo "FAIL: another job holds .running_K3 in $W"; exit 1; fi',
            'if ! mkdir .running_%s 2>/dev/null; then echo "FAIL: another job holds .running_%s in $W"; exit 1; fi' % (jk, jk)),
           ("trap 'rmdir \"$W/.running_K3\" 2>/dev/null' EXIT", "trap 'rmdir \"$W/.running_%s\" 2>/dev/null' EXIT" % jk)]
    for a, b in rep:
        if head.count(a) != 1:
            stop("%s does not hold the line: %s" % (PBS53, a))
        head[head.index(a)] = b
    cap53, cap = 'timeout -k 10 %d "$ORCA/orca" job.inp' % CAP53, 'timeout -k 10 %d "$ORCA/orca" job.inp' % CAP
    if sum(cap53 in l for l in head) != 1:
        stop("%s does not hold the cap of %d s." % (PBS53, CAP53))
    head = [l.replace(cap53, cap) for l in head]
    if any(("d53" in l or "D53" in l or "DIAGNOSTIC 53" in l or "K3" in l) for l in head):
        stop("a line of %s still names D53 or its job after the names were changed." % PBS53)
    if len(mine) * (CAP + 10) / 3600.0 >= LIMIT_H:
        stop("job %s could reach its limit of %d hours." % (jk, LIMIT_H))
    return "\n".join(head + ["run %s" % n for n in mine] + [tail[0], loop % " ".join(mine), 'echo "finished: $n of %d"' % len(mine), ""])


def sources_list():
    f = [MANIFEST, PBS53]
    for r in runs_of():
        for x in ("%s/%s/job.inp" % (D53, r["source"]), "%s/%s/pc.pc" % (D53, r["source"])):
            if x not in f:
                f.append(x)
    return f


def runs_text(s6):
    st = dict((run, (src, h)) for run, src, h in start_table(s6))
    ls = ["\t".join(["run", "job", "level", "kind", "state", "its_input_is_that_of", "its_start_orbitals", "sha256_of_the_start", "basis_functions"])]
    for r in runs_of():
        ls.append("\t".join([r["run"], r["job"], r["level"], r["kind"], r["state"], "%s/%s" % (D53, r["source"]), st[r["run"]][0], st[r["run"]][1], str(r["nb"])]))
    return "\n".join(ls) + "\n"


def build(s6, code=None, hpc=HPC):
    """every file of d54, as {path under d54: text}; the checksum list last"""
    for f in sources_list():
        if not os.path.isfile(os.path.join(s6, f)):
            stop("the committed file %s is missing." % f)
    pcs = set(open(os.path.join(s6, D53, r["source"], "pc.pc")).read() for r in runs_of())
    if len(set(p.split("\n", 1)[0].strip() for p in pcs)) != 1 or [p.split("\n", 1)[0].strip() for p in pcs][0] != "11":
        stop("a pc.pc of the source runs does not hold eleven charges.")
    if len(set(tuple(l for l in p.split("\n")[1:12]) for p in pcs)) != 1:
        stop("the pc.pc files of the source runs do not hold the same eleven charges.")
    F = {}
    for r in runs_of():
        F["%s/job.inp" % r["run"]] = input_of(s6, r)
        F["%s/pc.pc" % r["run"]] = open(os.path.join(s6, D53, r["source"], "pc.pc")).read()
    F["RUNS.tsv"] = runs_text(s6)
    F["START_ORBITALS.sha256"] = "".join("%s  %s/%s\n" % (h, run, START) for run, src, h in start_table(s6))
    F["stage_orbitals.sh"] = stage_text(s6)
    F["SOURCES.sha256"] = "".join("%s  %s\n" % (sha(os.path.join(s6, f)), f) for f in sources_list())
    for job in JOBS:
        F["run_%s.pbs" % job[0]] = pbs_text(s6, job, hpc)
    F["INPUT_CHECKSUMS.sha256"] = "".join("%s  %s\n" % (sha_text(F[k]), k) for k in sorted(F))
    return F, runs_of()


def main():
    if len(sys.argv) not in (4, 5):
        sys.stderr.write("usage: python3 cs6_d54_prepare.py <results repository> <code repository> <output directory> [<directory of d54 on hpc1>]\n")
        return 1
    s6 = os.path.join(sys.argv[1], S6)
    out = sys.argv[3]
    hpc = sys.argv[4] if len(sys.argv) == 5 else HPC
    if os.path.exists(out):
        sys.stderr.write("STOP  %s exists already; this script overwrites nothing.\n" % out)
        return 1
    try:
        F, runs = build(s6, sys.argv[2], hpc)
    except Refused as e:
        sys.stderr.write("STOP  %s\n" % e)
        return 1
    for k in sorted(F):
        p = os.path.join(out, k)
        os.makedirs(os.path.dirname(p), exist_ok=True)
        with open(p, "w") as f:
            f.write(F[k])
        if k.endswith(".pbs") or k.endswith(".sh"):
            os.chmod(p, 0o755)
    print("%d runs of the dianion (%s), %d files checksummed, %d jobs" % (
        len(runs), ", ".join("%s %d" % (lv, sum(1 for r in runs if r["level"] == lv)) for lv in ("WXV", "CAM", "HF")), len(F) - 1, len(JOBS)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
