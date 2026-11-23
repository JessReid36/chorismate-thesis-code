#!/usr/bin/env python3
"""Stage 6, DIAGNOSTIC 53: charges of 0.25 e together, on the dianion, at two levels. THE INPUTS.

Forms every file of phase2.2/cs_stage6/d53 from committed files and from nothing else:

  F   the charges of an arrangement, and set X2 with treatment A at all six positive sites of B_n11_q0.25:
      D50's committed input of the design with those functions, its pc.pc cut down to the arrangement;
  GF  no charge, the same functions: that input without its line of point charges (D50's committed control);
  N   the same charges and no functions (the partner, for the gain): D44's committed input of the design;
  GN  no charge and no functions: D44's committed control.

Three levels: CAM, the committed first line as it stands; WXV, the same line with D37's committed words for the
functional; HF, the same line with HF. Every run of CAM and WXV gives hydrogen IntAcc 12 (the block D52 found to be
what a ghost atom then gets). The nine job scripts are D52's committed run_K1.pbs with its names, its runs and its cap
changed. Nothing is run here and nothing is read from an output but the number of basis functions two committed
outputs print, which the reader later asks of every run.

usage: python3 cs6_d53_prepare.py <results repository> <code repository> <output directory> [<directory of d53 on hpc1>]
"""
import hashlib
import os
import re
import sys

S6 = "phase2.2/cs_stage6"
HPC = "/home/18660916/system_development/phase2.2/cs_stage6/d53"
STATES = ("R", "TS")
SRC = {"F": "d50/CAMTZX_B_n11_q0.25_X2A_di_%s", "GF": "d50/CAMTZX_G_n11_q0.25_X2A_di_%s",
       "N": "d44/CAMTZ_B_n11_q0.25_di_%s", "GN": "d44/CAMTZ_C_ctl_%s"}
WXV_INP = "d37/wB97XV_TZ_R/job.inp"
PBS52 = "d52/run_K1.pbs"
D52 = "/home/18660916/system_development/phase2.2/cs_stage6/d52"
CAM_LINE = "! CAM-B3LYP D3BJ def2-TZVPD def2/J RIJCOSX TightSCF"
CAM_WORDS = "CAM-B3LYP D3BJ"
LEVELS = (("CAM", CAM_WORDS), ("WXV", "wB97X-V"), ("HF", "HF"))
LEVEL_NAME = {"CAM": "CAM-B3LYP-D3BJ", "WXV": "wB97X-V", "HF": "Hartree-Fock"}
GRID = ["%method", "  SpecialGridAtoms 1", "  SpecialGridIntAcc 12.0", "end"]      # D52's rung z1, word for word
GRID_LEVELS = ("CAM", "WXV")
PC_SITES = (557, 1124, 534, 1073, 1107, 1359, 1421, 417, 1460, 1005, 1052)       # the order of the lines of the committed pc.pc
POS = (1124, 1073, 1107, 1359, 1005, 1052)                                       # the order of the ghost atoms of the committed input
ARR = (("s1", (1124,)), ("pc", (1073, 1052)), ("pf", (1124, 1052)), ("t3", (1124, 1073, 1107)),
       ("f4", (1124, 1073, 1107, 1052)), ("x6", POS), ("d11", PC_SITES))
Q = 0.25
NB = {"F": 784, "GF": 784, "N": 682, "GN": 682}          # checked below against what the committed outputs print
NB_OUT = {"F": "d50/CAMTZX_B_n11_q0.25_X2A_di_R/job.out", "N": "d44/CAMTZ_C_ctl_R/job.out"}
KINDS = {"CAM": ("F", "N", "GF", "GN"), "WXV": ("F", "N", "GF", "GN"), "HF": ("F", "GF", "GN")}
# job, its PBS name, level, the kinds it holds (controls first), states, cap in seconds of one run, what it is
JOBS = (("K1", "cm_d53_k1", "CAM", ("GF", "F"), ("R",), 10200, "the working functional, with the functions, at R"),
        ("K2", "cm_d53_k2", "CAM", ("GF", "F"), ("TS",), 10200, "the working functional, with the functions, at TS"),
        ("K3", "cm_d53_k3", "WXV", ("GF", "F"), ("R",), 10200, "the second level, with the functions, at R"),
        ("K4", "cm_d53_k4", "WXV", ("GF", "F"), ("TS",), 10200, "the second level, with the functions, at TS"),
        ("K5", "cm_d53_k5", "CAM", ("GN", "N"), ("R",), 9000, "the working functional, without the functions, at R"),
        ("K6", "cm_d53_k6", "CAM", ("GN", "N"), ("TS",), 9000, "the working functional, without the functions, at TS"),
        ("K7", "cm_d53_k7", "WXV", ("GN", "N"), ("R",), 9000, "the second level, without the functions, at R"),
        ("K8", "cm_d53_k8", "WXV", ("GN", "N"), ("TS",), 9000, "the second level, without the functions, at TS"),
        ("K9", "cm_d53_k9", "HF", ("GN", "GF", "F"), ("R", "TS"), 4500, "Hartree-Fock, the bracket"))
LIMIT_H = 24.0


class Refused(Exception):
    pass


def stop(m):
    raise Refused(m)


def sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


def sha_text(t):
    return hashlib.sha256(t.encode()).hexdigest()


def lines_of(path):
    ls = open(path).read().split("\n")
    while ls and not ls[-1].strip():
        ls.pop()
    return ls


def name_of(level, kind, arr, state):
    return "DI6_%s_%s_%s" % (level, kind, state) if kind in ("GF", "GN") else "DI6_%s_%s_%s_%s" % (level, kind, arr, state)


def first_line(s6, level):
    words = dict(LEVELS)[level]
    if CAM_LINE.count(CAM_WORDS) != 1:
        stop("the committed first line does not hold its functional once.")
    line = CAM_LINE.replace(CAM_WORDS, words)
    if level == "WXV" and lines_of(os.path.join(s6, WXV_INP))[0] != line:
        stop("the first line formed for wB97X-V is not the first line of the committed %s." % WXV_INP)
    return line


def source(s6, kind, state):
    ls = lines_of(os.path.join(s6, SRC[kind] % state, "job.inp"))
    if ls[0] != CAM_LINE:
        stop("%s/job.inp does not begin with the committed first line." % (SRC[kind] % state))
    return ls


def input_of(s6, level, kind, state):
    ls = source(s6, kind, state)
    ls[0] = first_line(s6, level)
    if level in GRID_LEVELS:
        i = next((k for k, l in enumerate(ls) if l.startswith("%pointcharges")), None)
        if i is None:
            i = next(k for k, l in enumerate(ls) if l.startswith("*"))
        ls[i:i] = GRID
    return "\n".join(ls) + "\n"


def pc_lines(s6, state):
    ls = lines_of(os.path.join(s6, SRC["F"] % state, "pc.pc"))
    if ls[0].strip() != str(len(PC_SITES)) or len(ls) != len(PC_SITES) + 1:
        stop("the committed pc.pc of the design does not hold %d charges." % len(PC_SITES))
    return ls[1:]


def pc_of(s6, sites, state):
    body = pc_lines(s6, state)
    keep = [body[PC_SITES.index(s)] for s in PC_SITES if s in sites]
    if len(keep) != len(sites):
        stop("an arrangement names a site that is not a site of the design.")
    return "%d\n" % len(keep) + "\n".join(keep) + "\n"


def charges_of(s6, sites, state):
    """[(charge, (x, y, z))] of an arrangement, in the order of the committed pc.pc"""
    out = []
    for s, l in zip(PC_SITES, pc_lines(s6, state)):
        f = l.split()
        if s in sites:
            out.append((float(f[0]), tuple(float(v) for v in f[1:4])))
    return out


def geometry(s6, state):
    """the 24 atoms and the six ghost centres of the committed input with functions: [(element, x, y, z)], [(x, y, z)]"""
    atoms, ghosts = [], []
    ls = source(s6, "F", state)
    i0 = next(k for k, l in enumerate(ls) if l.startswith("*"))
    for l in ls[i0 + 1:]:
        m = re.match(r"\s*([A-Z][a-z]?)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s*$", l)
        g = re.match(r"\s*H\s*:\s*(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s*$", l)
        if m:
            atoms.append((m.group(1),) + tuple(float(m.group(k)) for k in (2, 3, 4)))
        elif g:
            ghosts.append(tuple(float(g.group(k)) for k in (1, 2, 3)))
    return atoms, ghosts


def check_sources(s6):
    """the committed files agree with one another where this script assumes they do"""
    for st in STATES:
        f, gf = source(s6, "F", st), source(s6, "GF", st)
        if [l for l in f if not l.startswith("%pointcharges")] != gf or sum(l.startswith("%pointcharges") for l in f) != 1:
            stop("D50's control with functions is not D50's design input without its line of point charges (%s)." % st)
        if lines_of(os.path.join(s6, SRC["N"] % st, "pc.pc")) != lines_of(os.path.join(s6, SRC["F"] % st, "pc.pc")):
            stop("the two committed pc.pc of the design differ (%s)." % st)
        n, gn = source(s6, "N", st), source(s6, "GN", st)
        if [l for l in n if not l.startswith("%pointcharges")] != gn:
            stop("D44's control is not D44's design input without its line of point charges (%s)." % st)
        sx = [open(os.path.join(s6, SRC[k] % st, "substrate.xyz")).read() for k in ("N", "GN")]
        if sx[0] != sx[1]:
            stop("the two committed substrate.xyz differ (%s)." % st)
        atoms, ghosts = geometry(s6, st)
        xyz = [l.split() for l in sx[0].split("\n")[2:26]]
        if len(atoms) != 24 or len(ghosts) != len(POS) or len(xyz) != 24:
            stop("the committed input with functions does not hold 24 atoms and six ghost centres (%s)." % st)
        for a, b in zip(atoms, xyz):
            if a[0] != b[0] or max(abs(a[k + 1] - float(b[k + 1])) for k in range(3)) > 1e-9:
                stop("the atoms of the input with functions are not those of the committed substrate.xyz (%s)." % st)
        body = pc_lines(s6, st)
        for s, l in zip(PC_SITES, body):
            if abs(abs(float(l.split()[0])) - Q) > 1e-12 or (float(l.split()[0]) > 0) != (s in POS):
                stop("the committed pc.pc does not hold +0.25 e at the positive sites and -0.25 e at the others (%s)." % st)
        for s, g in zip(POS, ghosts):
            c = [float(v) for v in body[PC_SITES.index(s)].split()[1:4]]
            d = sum((c[k] - g[k]) ** 2 for k in range(3)) ** 0.5
            if not 0.005 < d < 0.02:
                stop("a ghost centre does not lie 0.01 A from the charge at site %d (%s): %.4f A." % (s, st, d))
    for k, rel in NB_OUT.items():
        t = open(os.path.join(s6, rel), errors="replace").read()
        nb = sorted(set(int(v) for v in re.findall(r"Number of basis functions\s*\.+\s*(\d+)", t)))
        if nb != [NB[k]]:
            stop("%s does not print %d basis functions." % (rel, NB[k]))


def runs_of():
    runs = []
    for jk, pbs, level, kinds, states, cap, what in JOBS:
        for st in states:
            for kind in kinds:
                for arr, sites in (ARR if kind in ("F", "N") else (("-", ()),)):
                    r = dict(run=name_of(level, kind, arr, st), job=jk, level=level, kind=kind, arr=arr, state=st, nb=NB[kind], ncharge=len(sites),
                             control={"F": name_of(level, "GF", "-", st), "N": name_of(level, "GN", "-", st), "GF": name_of(level, "GN", "-", st), "GN": "-"}[kind],
                             partner=name_of(level, "N", arr, st) if (kind == "F" and "N" in KINDS[level]) else "-")
                    runs.append(r)
    want = sum(len(KINDS[lv]) for lv, _ in LEVELS)
    names = [r["run"] for r in runs]
    if len(set(names)) != len(names) or len(set((r["level"], r["kind"]) for r in runs)) != want:
        stop("the list of runs is not the one this script was written for.")
    return runs


def runs_text(runs):
    cols = ("run", "job", "level", "kind", "arr", "state", "ncharge", "nb", "control", "partner")
    head = ("run", "job", "level", "kind", "arrangement", "state", "point_charges", "basis_functions", "its_control", "its_partner_without_functions")
    return "\t".join(head) + "\n" + "".join("\t".join(str(r[c]) for c in cols) + "\n" for r in runs)


def arrangements_text(s6):
    ls = ["arrangement\tsites\tcharges_e\tsum_e\n"]
    for arr, sites in ARR:
        order = [s for s in PC_SITES if s in sites]
        qs = [Q if s in POS else -Q for s in order]
        ls.append("%s\t%s\t%s\t%+.2f\n" % (arr, ",".join(str(s) for s in order), ",".join("%+.2f" % q for q in qs), sum(qs)))
    return "".join(ls)


def levels_text(s6):
    ls = ["level\tfirst_line\tblock_added\n"]
    for lv, _ in LEVELS:
        ls.append("%s\t%s\t%s\n" % (lv, first_line(s6, lv), " | ".join(GRID) if lv in GRID_LEVELS else "-"))
    return "".join(ls)


def pbs_text(s6, job, runs, hpc):
    """the job script: D52's committed run_K1.pbs with its names, its runs and its cap changed"""
    jk, pbs, level, kinds, states, cap, what = job
    mine = [r["run"] for r in runs if r["job"] == jk]
    ls = open(os.path.join(s6, PBS52)).read().split("\n")
    i0 = next(i for i, l in enumerate(ls) if l.startswith("run "))
    i1 = next(i for i, l in enumerate(ls) if l.startswith('echo "end='))
    head, mid, tail = ls[:i0], ls[i0:i1], ls[i1:]
    names = [l[4:] for l in mid]
    loop = "n=0; for r in %s; do grep -q 'ORCA TERMINATED NORMALLY' \"$W/$r/job.out\" 2>/dev/null && n=$((n+1)); done"
    if any(not l.startswith("run ") for l in mid) or len(tail) != 4 or tail[3] != "" or tail[1] != loop % " ".join(names) \
            or tail[2] != 'echo "finished: $n of %d"' % len(names):
        stop("%s is not the job script this script was written for." % PBS52)
    rep = [("#PBS -N cm_d52_k1", "#PBS -N %s" % pbs),
           ("#PBS -o %s/run_K1.pbs.out" % D52, "#PBS -o %s/run_%s.pbs.out" % (hpc, jk)),
           ("# Stage 6 DIAGNOSTIC 52: charges together, one electron first.",
            "# Stage 6 DIAGNOSTIC 53: charges of 0.25 e together, on the dianion, at two levels."),
           ("# This job is the only one. It makes %d runs in the order below; a run that fails does not stop the ones after it." % len(names),
            "# Job %s of %d: %s. It makes %d runs in the order below; a run that fails does not stop the ones after it." % (jk, len(JOBS), what, len(mine))),
           ("W=%s" % D52, "W=%s" % hpc),
           ('if ! mkdir .running_K1 2>/dev/null; then echo "FAIL: another job holds .running_K1 in $W"; exit 1; fi',
            'if ! mkdir .running_%s 2>/dev/null; then echo "FAIL: another job holds .running_%s in $W"; exit 1; fi' % (jk, jk)),
           ("trap 'rmdir \"$W/.running_K1\" 2>/dev/null' EXIT", "trap 'rmdir \"$W/.running_%s\" 2>/dev/null' EXIT" % jk),
           ('  else ( cd "$d" && timeout -k 10 3600 "$ORCA/orca" job.inp > job.out 2>&1 </dev/null ); echo "$1 orca_exit=$? normal=$(grep -c \'ORCA TERMINATED NORMALLY\' "$d/job.out") $(date)"; fi',
            '  else ( cd "$d" && timeout -k 10 %d "$ORCA/orca" job.inp > job.out 2>&1 </dev/null ); echo "$1 orca_exit=$? normal=$(grep -c \'ORCA TERMINATED NORMALLY\' "$d/job.out") $(date)"; fi' % cap)]
    for a, b in rep:
        if head.count(a) != 1:
            stop("%s does not hold the line: %s" % (PBS52, a))
        head[head.index(a)] = b
    if any(("d52" in l or "D52" in l or "DIAGNOSTIC 52" in l or "K1" in l.replace(jk, "")) for l in head):
        stop("a line of %s still names D52 or its job after the names were changed." % PBS52)
    if len(mine) * (cap + 10) / 3600.0 >= LIMIT_H:
        stop("job %s could reach its limit of %d hours: %d runs of %d s." % (jk, LIMIT_H, len(mine), cap))
    return "\n".join(head + ["run %s" % n for n in mine] + [tail[0], loop % " ".join(mine), 'echo "finished: $n of %d"' % len(mine), ""])


def sources_list():
    fs = [PBS52, WXV_INP, "d52_rules.tsv", "d50_dianion.tsv"] + list(NB_OUT.values())
    for st in STATES:
        fs += [SRC["F"] % st + "/job.inp", SRC["F"] % st + "/pc.pc", SRC["GF"] % st + "/job.inp",
               SRC["N"] % st + "/job.inp", SRC["N"] % st + "/pc.pc", SRC["N"] % st + "/substrate.xyz",
               SRC["GN"] % st + "/job.inp", SRC["GN"] % st + "/substrate.xyz"]
    return sorted(set(fs))


def build(s6, code=None, hpc=HPC):
    """every file of d53, as {path under d53: text}; the checksum list last"""
    for f in sources_list():
        if not os.path.isfile(os.path.join(s6, f)):
            stop("the committed file %s is missing." % f)
    check_sources(s6)
    runs = runs_of()
    arr = dict(ARR)
    F = {}
    for r in runs:
        F["%s/job.inp" % r["run"]] = input_of(s6, r["level"], r["kind"], r["state"])
        if r["kind"] in ("F", "N"):
            F["%s/pc.pc" % r["run"]] = pc_of(s6, arr[r["arr"]], r["state"])
        if r["kind"] in ("N", "GN"):
            F["%s/substrate.xyz" % r["run"]] = open(os.path.join(s6, SRC["GN"] % r["state"], "substrate.xyz")).read()
    F["RUNS.tsv"] = runs_text(runs)
    F["ARRANGEMENTS.tsv"] = arrangements_text(s6)
    F["LEVELS.tsv"] = levels_text(s6)
    F["SOURCES.sha256"] = "".join("%s  %s\n" % (sha(os.path.join(s6, f)), f) for f in sources_list())
    for job in JOBS:
        F["run_%s.pbs" % job[0]] = pbs_text(s6, job, runs, hpc)
    F["INPUT_CHECKSUMS.sha256"] = "".join("%s  %s\n" % (sha_text(F[k]), k) for k in sorted(F))
    return F, runs


def main():
    if len(sys.argv) not in (4, 5):
        sys.stderr.write("usage: python3 cs6_d53_prepare.py <results repository> <code repository> <output directory> [<directory of d53 on hpc1>]\n")
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
        if k.endswith(".pbs"):
            os.chmod(p, 0o755)
    print("%d runs of the dianion (%s), %d files checksummed, %d jobs" % (
        len(runs), ", ".join("%s %d" % (lv, sum(1 for r in runs if r["level"] == lv)) for lv, _ in LEVELS), len(F) - 1, len(JOBS)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
