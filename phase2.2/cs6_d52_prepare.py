#!/usr/bin/env python3
"""cs6_d52_prepare.py - Stage 6 DIAGNOSTIC 52: charges together, one electron first.

D50 found that the design B_n11_q0.25 gave up half an electron at R to the functions at its positive charges, and
that one charge of +0.25 e alone did not. D51 measured, for ONE charge, which functions and which grid let ORCA
describe an electron held there, and how far the functional lies from exact. D52 measures the same for one electron
held by SEVERAL charges of +0.25 e at once, before any run of the dianion with charges together (AMENDMENT 10).

One electron and nothing else: no substrate, no negative charge. Six systems of charges of +0.25 e at positive sites
of B_n11_q0.25, with the functions at ALL SIX positive sites in every run, as the runs of the dianion will have them:

  s1   site 1124                         pc   1073 and 1052, the closest two        pf   1124 and 1052, the farthest two
  t3   1124, 1073 and 1107               f4   1124, 1073, 1107 and 1052             x6   all six

  HF     Hartree-Fock (exact for one electron in the functions given; no grid), sets X2 and X3.
  CAM    the working functional without its dispersion term, auxiliary treatment A, on IntAcc 12, 20 and 40, sets X2
         and X3. At (X2, s1) also on ORCA's default grid (g2) and with IntAcc 12 given to hydrogen alone (z1).
  WXV, WMV, RDS   three candidates for a second level: wB97X-V and wB97M-V (the functionals of D37's committed
         inputs) and revDSD-PBEP86-D4/2021 (the committed line of D44's level RDTZ), on IntAcc 12 and 20 in set X2
         and on IntAcc 20 in set X3.

Every input is D49's committed one-electron input (formed through the committed cs6_d50_prepare.py, loaded by path
and checked by sha256) with its one ghost atom replaced by the six ghost atoms of the committed input of the design
in D50, and for a functional the line of keywords of its level and at most two blocks added. The point charges are
lines of the committed point-charge file of the design. This script reads committed files only and edits none.

usage: python3 cs6_d52_prepare.py <results repository> <code repository> <output directory> [<directory of d52 on hpc1>]
"""
import hashlib
import importlib.util
import os
import sys

S6 = "phase2.2/cs_stage6"
P50 = "cs6_d50_prepare.py"
P50_SHA = "bdbfce5202644a9e2d31c557377d5139007eccc099a38b00225150af2a4b3f82"
EXPF = "cs6_d52_expected.py"
EXP_SHA = "4ab741567ba717331d5c4c7eb2ec41d0394d19d8c0be5796f87bda638e2aa3dc"
HPC = "/home/18660916/system_development/phase2.2/cs_stage6/d52"
AUX = "A"                                                  # the auxiliary treatment D51's rule gives
SRC = "d50/CAMTZX_B_n11_q0.25_X2A_di_R"                    # committed: the design's input (its six ghost atoms) and its point charges
RD = "d44/RDTZ_C_ctl_R"                                    # committed: the line of keywords of revDSD-PBEP86-D4/def2-TZVPD
D37 = {"WXV": "d37/wB97XV_TZ_R", "WMV": "d37/wB97MV_TZ_R"}  # committed: where each of these functionals ran on the dianion
D37_LINE = "! %s def2-TZVPD def2/J RIJCOSX TightSCF"
RD_LINE = "! revDSD-PBEP86-D4/2021 def2-TZVPD def2/J def2-TZVPD/C RIJCOSX TightSCF"
RD_MP2 = "%mp2 RI true end"
PBS51 = "d51/run_K1.pbs"                                   # committed: the job script this one is made from
SITES = (1124, 1073, 1107, 1359, 1005, 1052)               # the positive charges of B_n11_q0.25, in the order of its point-charge file
SYSTEMS = (("s1", (1124,)), ("pc", (1073, 1052)), ("pf", (1124, 1052)), ("t3", (1124, 1073, 1107)),
           ("f4", (1124, 1073, 1107, 1052)), ("x6", SITES))
SETS = ("X2", "X3")
Q = 0.25                                                   # e, every charge
GHOST_AT = 0.01                                            # A: a ghost atom sits this far beyond its charge (D47's rule)
# functional: (key, the word that stands for CAM-B3LYP in the line of keywords, or None for the committed line of RDTZ)
FUNCS = (("CAM", "CAM-B3LYP"), ("WXV", "wB97X-V"), ("WMV", "wB97M-V"), ("RDS", None))
CANDIDATES = ("WXV", "WMV", "RDS")
# rung: the block added before the line of the point charges
RUNGS = (("g2", None), ("z1", ("SpecialGridAtoms 1", "SpecialGridIntAcc 12.0")), ("i12", ("IntAcc 12.0",)), ("i20", ("IntAcc 20.0",)),
         ("i40", ("IntAcc 40.0",)))
INTACC = {"i12": 12.0, "i20": 20.0, "i40": 40.0}
LADDER = ("i12", "i20", "i40")                             # the working functional's ladder
CAND_RUNGS = {"X2": ("i12", "i20"), "X3": ("i20",)}        # the candidates' rungs
TOP = "i20"                                                # the rung at which functionals are set against the limit and against one another
FORM = ("X2", "s1")                                        # where the default grid and the grid for hydrogen alone are run
TOL_ADEQ = 0.05                                            # eV, D50's and D51's: the functions hold within this of the limit
TOL_GRID = 0.005                                           # eV, D51's: two rungs agree to this
TOL_TWO = 0.001                                            # eV: ORCA's Hartree-Fock against PySCF's in the same functions
TOL_EH = 1e-6                                              # hartree: two runs give the same energy
D51_ERR = 0.0631                                           # eV: D51's result, the functional's largest distance from exact at one charge
FEW = 3                                                    # RULE C: the limit is this many times the functional's largest distance
JOBS = [("K1", "cm_d52_k1")]
CAP = 3600                                                 # s: a run is stopped after this (D51's job stopped a run after 900)
EXTRA = ["STAGE6_AMENDMENT10_20261123.txt", "STAGE6_DIAGNOSTIC50_20261122.txt", "STAGE6_DIAGNOSTIC50_RESULT_20261123.txt",
         "STAGE6_DIAGNOSTIC51_20261123.txt", "STAGE6_DIAGNOSTIC51_RESULT_20261123.txt", "STAGE6_DIAGNOSTIC44_NOTE3_20261123.txt",
         "d43/SITES.tsv", "d51.tsv", "d51_rule.tsv", "%s/job.inp" % SRC, "%s/pc.pc" % SRC, "%s/job.inp" % RD,
         "%s/job.inp" % D37["WXV"], "%s/job.inp" % D37["WMV"], PBS51]


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


def geometry(s6):
    """the six positive charges of the design (lines of its committed point-charge file) and its six ghost atoms
    (lines of its committed input), each checked against the committed table of sites"""
    pcl = open(os.path.join(s6, SRC, "pc.pc")).read().split("\n")
    if pcl[-1] != "" or int(pcl[0].split()[0]) != len(pcl) - 2:
        stop("%s/pc.pc is not a point-charge file of the committed form." % SRC)
    pos = [l for l in pcl[1:-1] if float(l.split()[0]) > 0]
    gh = [l for l in open(os.path.join(s6, SRC, "job.inp")).read().split("\n") if l.startswith("H: ")]
    if len(pos) != len(SITES) or len(gh) != len(SITES):
        stop("%s does not hold %d positive charges and %d ghost atoms." % (SRC, len(SITES), len(SITES)))
    rows = [l.split("\t") for l in open(os.path.join(s6, "d43", "SITES.tsv")).read().split("\n")[1:] if l]
    xyz = dict((int(r[0]), tuple(float(v) for v in r[2:5])) for r in rows)
    out = {}
    for site, pl, gl in zip(SITES, pos, gh):
        q, c, g = float(pl.split()[0]), tuple(float(v) for v in pl.split()[1:4]), tuple(float(v) for v in gl.split()[1:4])
        if abs(q - Q) > 1e-9 or site not in xyz or max(abs(a - b) for a, b in zip(c, xyz[site])) > 1e-6:
            stop("positive charge %d of %s/pc.pc is not +%.2f e at site %d of d43/SITES.tsv." % (SITES.index(site) + 1, SRC, Q, site))
        d = sum((a - b) ** 2 for a, b in zip(c, g)) ** 0.5
        if abs(d - GHOST_AT) > 2e-6:
            stop("ghost atom %d of %s/job.inp is %.6f A from its charge, not %.2f." % (SITES.index(site) + 1, SRC, d, GHOST_AT))
        out[site] = dict(pc_line=pl, ghost_line=gl, xyz=c, ghost=g)
    return out


def method_block(lines):
    return ["%method"] + ["  " + l for l in lines] + ["end"]


def hf_name(name, sysk):
    return "ONE6_HF_%s_%s" % (name, sysk)


def fn_name(func, name, rung, sysk):
    return "ONE6_%s_%s_%s_%s" % (func, name, rung, sysk)


def first_line(s6, p48, func):
    cam = dict(p48.ONE)["CAM"]
    word = dict(FUNCS)[func]
    if func == "CAM":
        return cam
    if word is not None:
        if cam.count("CAM-B3LYP") != 1:
            stop("the committed line of the functional does not hold CAM-B3LYP once.")
        line = cam.replace("CAM-B3LYP", word)
        if open(os.path.join(s6, D37[func], "job.inp")).read().split("\n")[0] != D37_LINE % word or line != "! UKS " + (D37_LINE % word)[2:]:
            stop("%s/job.inp does not begin with the committed line of %s this script was written for." % (D37[func], word))
        return line
    ls = open(os.path.join(s6, RD, "job.inp")).read().split("\n")
    if ls[0] != RD_LINE or ls.count(RD_MP2) != 1:
        stop("%s/job.inp is not the committed input of level RDTZ this script was written for." % RD)
    return "! UKS " + RD_LINE[2:]


def input_of(s6, p50, p48, geo, level, name, rung=None):
    """level: HF or a key of FUNCS. D49's one-electron input with the six ghost atoms of the design, the functions of
    the set under each, and for a functional the block of its rung"""
    ls = p50.one_plain(p48, "HF" if level == "HF" else "CAM").split("\n")
    one = "H: %.6f 0.000000 0.000000" % p48.GHOST_OFF
    if ls.count(one) != 1 or not ls[0].startswith("! "):
        stop("D49's one-electron input does not have the form this script was written for.")
    k = ls.index(one)
    ls = ls[:k] + [geo[s]["ghost_line"] for s in SITES] + ls[k + 1:]
    add = []
    if level != "HF":
        ls[0] = first_line(s6, p48, level)
        if dict(FUNCS)[level] is None:
            add.append(RD_MP2)
        block = dict(RUNGS)[rung]
        if block is not None:
            add += method_block(block)
    if add:
        if ls.count(p48.PC_LINE) != 1:
            stop("an input does not hold the line of the point charges once.")
        k = ls.index(p48.PC_LINE)
        ls = ls[:k] + add + ls[k:]
    return p50.with_block("\n".join(ls), name, None if level == "HF" else AUX, len(SITES))


def pc_of(geo, sites):
    return "%d\n" % len(sites) + "".join(geo[s]["pc_line"] + "\n" for s in sites)


def runs_of():
    out = []

    def add(level, name, rung, sysk):
        out.append(dict(run=hf_name(name, sysk) if level == "HF" else fn_name(level, name, rung, sysk), job="K1", level=level, set=name,
                        rung=rung or "", system=sysk))
    for name in SETS:
        for sysk, _ in SYSTEMS:
            add("HF", name, None, sysk)
    for rung in ("i12", "i20"):
        for name in SETS:
            for sysk, _ in SYSTEMS:
                add("CAM", name, rung, sysk)
    for rung in ("g2", "z1"):
        add("CAM", FORM[0], rung, FORM[1])
    for func in CANDIDATES:
        for name in SETS:
            for rung in CAND_RUNGS[name]:
                for sysk, _ in SYSTEMS:
                    add(func, name, rung, sysk)
    for name in SETS:
        for sysk, _ in SYSTEMS:
            add("CAM", name, "i40", sysk)                    # the slowest runs last
    if len(set(r["run"] for r in out)) != len(out):
        stop("two runs carry one name.")
    return out


def runs_text(runs):
    L = ["run\tjob\tlevel\tset\trung\tsystem\tsites_charged"]
    sy = dict(SYSTEMS)
    for r in runs:
        L.append("\t".join([r["run"], r["job"], r["level"], r["set"], r["rung"], r["system"], ",".join(str(s) for s in sy[r["system"]])]))
    return "\n".join(L) + "\n"


def systems_text(geo, ex):
    L = ["system\tcharges\tsites_charged\tlimit_eV\tlimit_in_a_larger_basis_eV\tHF_in_X2_expected_eV\tHF_in_X3_expected_eV\t"
         "CAM_in_X2_expected_eV\tCAM_in_X3_expected_eV"]
    for sysk, sites in SYSTEMS:
        L.append("\t".join([sysk, str(len(sites)), ",".join(str(s) for s in sites), "%.6f" % ex.LIMIT[sysk], "%.6f" % ex.LIMIT2[sysk],
                            "%.6f" % ex.HF[("X2", sysk)], "%.6f" % ex.HF[("X3", sysk)], "%.6f" % ex.CAM[("X2", sysk)][0], "%.6f" % ex.CAM[("X3", sysk)][0]]))
    L.append("")
    L.append("site\tcharge_x_A\tcharge_y_A\tcharge_z_A\tghost_x_A\tghost_y_A\tghost_z_A")
    for s in SITES:
        L.append("\t".join([str(s)] + ["%.6f" % v for v in geo[s]["xyz"]] + ["%.6f" % v for v in geo[s]["ghost"]]))
    return "\n".join(L) + "\n"


def levels_text(s6, p48):
    L = ["level\tline_of_keywords\tblock_added"]
    L.append("HF\t%s\t-" % dict(p48.ONE)["HF"])
    for func, _ in FUNCS:
        L.append("%s\t%s\t%s" % (func, first_line(s6, p48, func), RD_MP2 if dict(FUNCS)[func] is None else "-"))
    L.append("")
    L.append("rung\tblock_added_before_the_point_charges")
    for rung, block in RUNGS:
        L.append("%s\t%s" % (rung, "none (ORCA's default grid)" if block is None else " | ".join(method_block(block))))
    return "\n".join(L) + "\n"


def pbs_text(s6, runs, hpc):
    """the job script: D51's committed run_K1.pbs with its names and its runs changed, a run stopped after CAP seconds
    in place of 900, and no run made a second time in its folder (D44 NOTE 3)"""
    jk, pbs = JOBS[0]
    mine = [r["run"] for r in runs]
    ls = open(os.path.join(s6, PBS51)).read().split("\n")
    i0 = next(i for i, l in enumerate(ls) if l.startswith("run "))
    i1 = next(i for i, l in enumerate(ls) if l.startswith('echo "end='))
    head, mid, tail = ls[:i0], ls[i0:i1], ls[i1:]
    names = [l[4:] for l in mid]
    loop = "n=0; for r in %s; do grep -q 'ORCA TERMINATED NORMALLY' \"$W/$r/job.out\" 2>/dev/null && n=$((n+1)); done"
    if any(not l.startswith("run ") for l in mid) or len(tail) != 4 or tail[3] != "" or tail[1] != loop % " ".join(names) \
            or tail[2] != 'echo "finished: $n of %d"' % len(names):
        stop("%s is not the job script this script was written for." % PBS51)
    d51 = "/home/18660916/system_development/phase2.2/cs_stage6/d51"
    rep = [("#PBS -N cm_d51_k1", "#PBS -N %s" % pbs),
           ("#PBS -o %s/run_K1.pbs.out" % d51, "#PBS -o %s/run_%s.pbs.out" % (hpc, jk)),
           ("# Stage 6 DIAGNOSTIC 51: the calibration of D50 done again, with the integration grid measured.",
            "# Stage 6 DIAGNOSTIC 52: charges together, one electron first."),
           ("# This job is the only one. It makes %d runs in the order below; a run that fails does not stop the ones after it." % len(names),
            "# This job is the only one. It makes %d runs in the order below; a run that fails does not stop the ones after it." % len(mine)),
           ("# No -q: the allocation is automatic and walltime 24:00:00 routes to day. Restartable:",
            "# No -q: the allocation is automatic and walltime 24:00:00 routes to day. No run is made twice:"),
           ("#   submitted again, it skips every run that ended normally.",
            "#   submitted again, it makes only the runs whose folder holds no output yet (D44 NOTE 3)."),
           ("W=%s" % d51, "W=%s" % hpc),
           ('  if [ -s "$d/job.out" ] && grep -q "ORCA TERMINATED NORMALLY" "$d/job.out"; then echo "$1 done before"',
            '  if [ -e "$d/job.out" ]; then echo "$1 has an output already: not run again"'),
           ('  else ( cd "$d" && timeout -k 10 900 "$ORCA/orca" job.inp > job.out 2>&1 </dev/null ); echo "$1 orca_exit=$? normal=$(grep -c \'ORCA TERMINATED NORMALLY\' "$d/job.out") $(date)"; fi',
            '  else ( cd "$d" && timeout -k 10 %d "$ORCA/orca" job.inp > job.out 2>&1 </dev/null ); echo "$1 orca_exit=$? normal=$(grep -c \'ORCA TERMINATED NORMALLY\' "$d/job.out") $(date)"; fi' % CAP)]
    for a, b in rep:
        if head.count(a) != 1:
            stop("%s does not hold the line: %s" % (PBS51, a))
        head[head.index(a)] = b
    if any(("d51" in l or "D51" in l or "DIAGNOSTIC 51" in l) for l in head):
        stop("a line of %s still names D51 after the names were changed." % PBS51)
    return "\n".join(head + ["run %s" % n for n in mine] + [tail[0], loop % " ".join(mine), 'echo "finished: $n of %d"' % len(mine), ""])


def build(s6, code, hpc=HPC):
    """every file of d52, as {path under d52: text}; the checksum list last"""
    p50 = load_by_path(os.path.join(code, "phase2.2"), P50, P50_SHA, "d52_cs6_d50_prepare")      # committed: from the code repository given
    ex = load_by_path(os.path.dirname(os.path.abspath(__file__)), EXPF, EXP_SHA, "d52_cs6_d52_expected")   # committed beside this script
    try:
        p48 = p50.load48(code)
    except p50.Refused as e:
        stop("cs6_d50_prepare.py: %s" % e)
    for name in SETS:
        if name not in p50.ADD:
            stop("set %s is not a set of cs6_d50_prepare.py." % name)
    for f in EXTRA:
        if not os.path.isfile(os.path.join(s6, f)):
            stop("the committed file %s is missing." % f)
    geo = geometry(s6)
    for s in SITES:
        if max(abs(a - b) for a, b in zip(geo[s]["xyz"], ex.SITE[s])) > 1e-9 or max(abs(a - b) for a, b in zip(geo[s]["ghost"], ex.GHOST[s])) > 1e-9:
            stop("what Claude expects was computed for other positions than the committed ones (site %d)." % s)
    runs = runs_of()
    sy = dict(SYSTEMS)
    F = {}
    for r in runs:
        try:
            F["%s/job.inp" % r["run"]] = input_of(s6, p50, p48, geo, r["level"], r["set"], r["rung"] or None)
        except p50.Refused as e:
            stop("cs6_d50_prepare.py: %s" % e)
        F["%s/pc.pc" % r["run"]] = pc_of(geo, sy[r["system"]])
    F["RUNS.tsv"] = runs_text(runs)
    F["SYSTEMS.tsv"] = systems_text(geo, ex)
    F["LEVELS.tsv"] = levels_text(s6, p48)
    F["SOURCES.sha256"] = "".join("%s  %s\n" % (sha(os.path.join(s6, f)), f) for f in EXTRA)
    F["run_K1.pbs"] = pbs_text(s6, runs, hpc)
    F["INPUT_CHECKSUMS.sha256"] = "".join("%s  %s\n" % (sha_text(F[k]), k) for k in sorted(F))
    return F, runs, p50, p48, ex, geo


def main():
    if len(sys.argv) not in (4, 5):
        sys.stderr.write("usage: python3 cs6_d52_prepare.py <results repository> <code repository> <output directory> [<directory of d52 on hpc1>]\n")
        return 1
    s6 = os.path.join(sys.argv[1], S6)
    code, out = sys.argv[2], sys.argv[3]
    hpc = sys.argv[4] if len(sys.argv) == 5 else HPC
    if os.path.exists(out):
        sys.stderr.write("STOP  %s exists already; this script overwrites nothing.\n" % out)
        return 1
    try:
        F, runs, p50, p48, ex, geo = build(s6, code, hpc)
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
    print("%d runs of one electron (%d Hartree-Fock, %d with the working functional, %d with a candidate), %d files checksummed, %d job" % (
        len(runs), sum(1 for r in runs if r["level"] == "HF"), sum(1 for r in runs if r["level"] == "CAM"),
        sum(1 for r in runs if r["level"] in CANDIDATES), len(F) - 1, len(JOBS)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
