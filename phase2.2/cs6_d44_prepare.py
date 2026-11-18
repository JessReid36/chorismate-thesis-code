#!/usr/bin/env python3
"""cs6_d44_prepare.py - DIAGNOSTIC 44: the cross-check of AMENDMENT 7 on the six survivors of D43.

196 single points on frame 42436: the six designs that DIAGNOSTIC 43 lists as survivors of the
search, at the five levels AMENDMENT 7 names in A7.3, with bare point charges.

  the five levels   CAM-B3LYP-D3BJ/def2-TZVPD                       (CAMTZ)   28 runs
                    PWPB95-D3BJ/def2-SVPD and /def2-TZVPD           (PWSV, PWTZ)   42 runs each
                    revDSD-PBEP86-D4/def2-SVPD and /def2-TZVPD      (RDSV, RDTZ)   42 runs each

  at CAMTZ, a hybrid, every run is made as the working level makes it: the dianion (di) and the
  monoanion (mo), each with the stability analysis, at the reactant and at the transition state.

  at a double hybrid the dianion is run twice at each state. di is the committed charge-free
  input of the level with the point charges added and nothing else, so its energy is on the
  restricted solution that the record's comparisons of functionals used. ds is the same dianion
  with the stability analysis: the search for a lower self-consistent solution. The monoanion
  (mo) carries the stability analysis, as AMENDMENT 6 requires.

  each level has its controls with no charges, in the same kinds: C_ctl (di), C_ctlS (ds, at a
  double hybrid only) and C_bare (mo), at both states.

Nothing about a level or a design is typed in to be trusted.

  THE METHOD LINE of a level is line 1 of the committed charge-free input AMENDMENT 7 names for
  it, read from that file. The strings typed in below are a check on what was read, and nothing
  is written unless the file holds exactly that line at both states.

  WHAT IS ADDED to it. A run WITHOUT the analysis is the committed charge-free input of its
  level, every line, with the %pointcharges line of DIAGNOSTIC 31's committed charged input put
  before the geometry line; a control of that kind is the committed input byte for byte. A run
  WITH the analysis is DIAGNOSTIC 31's committed charged input, the input of the working level,
  with the method line of its level in place of its own, and the %mp2 block where the committed
  input of the level holds one. The monoanion changes the charge-and-multiplicity line from
  "-2 1" to "-1 2" and nothing else. A control has no %pointcharges line.

  THE SIX DESIGNS are the rows of the committed d43.tsv that say survivor, in the order of
  d43/CONDITIONS.tsv. Each point-charge file is the committed file of D43, byte for byte.

  THE GEOMETRIES are d33a's, byte for byte, which are also the geometries of the committed inputs
  AMENDMENT 7 names.

select() and everything it uses need the standard library only, so the reader imports them and
rebuilds every input on its own.

Usage:  cs6_d44_prepare.py <results-repo> <out-dir> <hpc-dir>
"""
import hashlib, os, shutil, sys

S6 = "phase2.2/cs_stage6"
GEOM = {"R": "d33a/sp_42436_R/substrate.xyz", "TS": "d33a/sp_42436_TS/substrate.xyz"}
TEMPLATE = "d31/B/n01_q0.25_R/job.inp"       # the committed charged dianion input of the working level
TEMPLATE_LINES = ["! CAM-B3LYP D3BJ def2-SVPD def2/J RIJCOSX TightSCF", "%maxcore 3000", "%pal nprocs 8 end", "%scf",
                  "  MaxIter 300", "  STABPerform true", "end", '%pointcharges "pc.pc"', "* xyzfile -2 1 substrate.xyz", ""]
PC_LINE = '%pointcharges "pc.pc"'
SCF_PLAIN = "%scf MaxIter 300 end"           # the %scf line of the committed charge-free inputs of the five levels
MP2_LINE = "%mp2 RI true end"
DI, MONO = "-2 1", "-1 2"
STATES = ("R", "TS")
# key, the level as AMENDMENT 7 writes it, the committed inputs that carry its method line (_R and _TS),
# whether it is a double hybrid, the name of its job, the method line as a check on what is read
LEVELS = [
    ("CAMTZ", "CAM-B3LYP-D3BJ/def2-TZVPD", "d37/CAM_TZ", False, "cm_d44_ct",
     "! CAM-B3LYP D3BJ def2-TZVPD def2/J RIJCOSX TightSCF"),
    ("PWSV", "PWPB95-D3BJ/def2-SVPD", "d37/PWPB95_SV", True, "cm_d44_ps",
     "! RI-PWPB95 D3BJ def2-SVPD def2/J def2-SVPD/C RIJCOSX TightSCF"),
    ("PWTZ", "PWPB95-D3BJ/def2-TZVPD", "d37/PWPB95_TZ", True, "cm_d44_pt",
     "! RI-PWPB95 D3BJ def2-TZVPD def2/J def2-TZVPD/C RIJCOSX TightSCF"),
    ("RDSV", "revDSD-PBEP86-D4/def2-SVPD", "d37a1/revDSD_RI_SV", True, "cm_d44_rs",
     "! revDSD-PBEP86-D4/2021 def2-SVPD def2/J def2-SVPD/C RIJCOSX TightSCF"),
    ("RDTZ", "revDSD-PBEP86-D4/def2-TZVPD", "d37a1/revDSD_RI_TZ", True, "cm_d44_rt",
     "! revDSD-PBEP86-D4/2021 def2-TZVPD def2/J def2-TZVPD/C RIJCOSX TightSCF"),
]
MP2_LEVELS = ("RDSV", "RDTZ")                # the levels whose committed input carries the %mp2 block
D43_COND = "d43/CONDITIONS.tsv"
D43_TSV = "d43.tsv"
D43_RESULT = "STAGE6_DIAGNOSTIC43_RESULT_20261117.txt"
# the kinds of run: tag, charge and multiplicity, whether the input carries the stability analysis
KINDS_HYBRID = (("di", DI, True), ("mo", MONO, True))
KINDS_DH = (("di", DI, False), ("ds", DI, True), ("mo", MONO, True))
CTL_NAME = {"di": "C_ctl", "ds": "C_ctlS", "mo": "C_bare"}
NDESIGN, NLEVEL, NRUN = 6, 5, 196
NPER = {False: 28, True: 42}                 # runs of a level: a hybrid, a double hybrid
NSUM = 568                                   # files checksummed: 168 x 3 + 28 x 2, three tables, five job scripts
WALL, MEM = "24:00:00", "48gb"
HPC = "/home/18660916/system_development/phase2.2/cs_stage6/d44"


class Refused(Exception):
    pass


def stop(m):
    raise Refused(m)


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 16), b""):
            h.update(b)
    return h.hexdigest()


def kinds(level):
    return KINDS_DH if level["dh"] else KINDS_HYBRID


def rname(level, cond, kind, state):
    """a run directory: <level>_<design>_<di|ds|mo>_<R|TS>"""
    return "%s_%s_%s_%s" % (level, cond, kind, state)


def cname(level, kind, state):
    """a control's run directory: <level>_C_ctl_<R|TS>, <level>_C_ctlS_<R|TS> or <level>_C_bare_<R|TS>"""
    return "%s_%s_%s" % (level, CTL_NAME[kind], state)


def table(p):
    ls = open(p).read().splitlines()
    hdr = ls[0].split("\t")
    return [dict(zip(hdr, l.split("\t"))) for l in ls[1:] if l.strip()]


def sources(designs):
    """every committed file this script or the reader reads, relative to phase2.2/cs_stage6"""
    out = [GEOM["R"], GEOM["TS"], TEMPLATE, D43_COND, D43_TSV, D43_RESULT,
           "STAGE6_CRITERIA.txt", "STAGE6_AMENDMENT6_20261117.txt", "STAGE6_AMENDMENT7_20261117.txt",
           "vde_probe/PROBE2_RESULT.txt",
           "d31/ctl_R/job.out", "d31/ctl_TS/job.out", "d40/B_bare_R/job.out", "d40/B_bare_TS/job.out"]
    for key, name, src, dh, job, line in LEVELS:
        for st in STATES:
            out += ["%s_%s/job.inp" % (src, st), "%s_%s/job.out" % (src, st)]
    for d in designs:
        out.append("d43/%s_di_R/pc.pc" % d)
        for sp in ("di", "mo"):
            for st in STATES:
                out.append("d43/%s_%s_%s/job.out" % (d, sp, st))
    return out


def committed_input(level):
    """the lines of the committed charge-free input of a level, as this script requires them to be"""
    return [level["line"], TEMPLATE_LINES[1], TEMPLATE_LINES[2], SCF_PLAIN] + ([MP2_LINE] if level["mp2"] else []) \
        + ["* xyzfile %s substrate.xyz" % DI, ""]


def select(s6):
    """The levels and the designs, from the committed files alone. Raises Refused where a committed
    file is not what this script was written for."""
    for f in (GEOM["R"], GEOM["TS"], TEMPLATE, D43_COND, D43_TSV):
        if not os.path.isfile(os.path.join(s6, f)):
            stop("%s is missing." % f)
    tl = open(os.path.join(s6, TEMPLATE)).read().split("\n")
    if tl != TEMPLATE_LINES or tl[7] != PC_LINE:
        stop("%s is not the committed charged input of the working level, line for line." % TEMPLATE)

    levels = []
    for key, name, src, dh, job, line in LEVELS:
        lv = dict(key=key, name=name, src=src, dh=dh, job=job, line=line, mp2=key in MP2_LEVELS)
        for st in STATES:
            p = os.path.join(s6, "%s_%s" % (src, st), "job.inp")
            if not os.path.isfile(p):
                stop("%s_%s/job.inp, which AMENDMENT 7 names for %s, is missing." % (src, st, name))
            if open(p).read().split("\n") != committed_input(lv):
                stop("%s_%s/job.inp is not the committed charge-free input of %s, line for line." % (src, st, name))
            g = os.path.join(s6, "%s_%s" % (src, st), "substrate.xyz")
            if not os.path.isfile(g) or sha(g) != sha(os.path.join(s6, GEOM[st])):
                stop("%s_%s/substrate.xyz is not frame 42436's committed %s geometry." % (src, st, st))
        levels.append(lv)
    if [l["key"] for l in levels] != ["CAMTZ", "PWSV", "PWTZ", "RDSV", "RDTZ"] or len(set(l["job"] for l in levels)) != NLEVEL \
            or [l["dh"] for l in levels] != [False, True, True, True, True]:
        stop("the five levels are not the five of AMENDMENT 7, in its order.")

    # the designs: D43's survivors, in the order of its CONDITIONS.tsv
    cond = table(os.path.join(s6, D43_COND))
    by = {}
    for r in table(os.path.join(s6, D43_TSV)):
        by.setdefault(r.get("condition"), []).append(r)
    designs = []
    for c in cond:
        rs = by.get(c.get("condition"), [])
        if sorted(r.get("state") for r in rs) != ["R", "TS"]:
            stop("d43.tsv does not hold one row at each state for %s." % c.get("condition"))
        surv = set(r.get("survivor") for r in rs)
        if surv == {"no"}:
            continue
        if surv != {"yes"} or c.get("part") != "B" or any(r.get("four_legs") != "RELIABLE" or r.get("fifth_leg") != "PASSED" for r in rs):
            stop("d43.tsv does not list %s as a design of part B that is RELIABLE and passes the fifth leg at both rows." % c.get("condition"))
        name = c["condition"]
        pcs = [os.path.join(s6, "d43", "%s_%s_%s" % (name, sp, st), "pc.pc") for sp in ("di", "mo") for st in STATES]
        if not all(os.path.isfile(p) for p in pcs) or len(set(sha(p) for p in pcs)) != 1:
            stop("the four committed point-charge files of %s are not one file." % name)
        n = int(c["n"])
        first = open(pcs[0]).read().split("\n")
        if int(first[0].split()[0]) != n or len([l for l in first[1:] if l.strip()]) != n:
            stop("d43/%s_di_R/pc.pc does not hold the %d charges CONDITIONS.tsv lists." % (name, n))
        qs = [float(l.split()[0]) for l in first[1:1 + n]]
        if ",".join("%+.4f" % v for v in qs) != c["charges"]:
            stop("d43/%s_di_R/pc.pc does not hold the charges CONDITIONS.tsv lists, in its order." % name)
        designs.append(dict(name=name, n=n, q=c["q_e"], sum=c["sum_q_e"], total=c["sum_abs_q_e"], sites=c["sites"],
                            charges=c["charges"], pc="d43/%s_di_R/pc.pc" % name, pcsha=sha(pcs[0])))
    if len(designs) != NDESIGN:
        stop("d43.tsv lists %d survivors, not %d." % (len(designs), NDESIGN))
    return dict(levels=levels, designs=designs, tl=tl)


def build_input(level, cm, charged, stab):
    """The text of one input.
    stab: DIAGNOSTIC 31's charged input with the level's method line, and its %mp2 block if it has one.
    not stab: the committed charge-free input of the level, with the %pointcharges line if charged."""
    t = TEMPLATE_LINES
    if stab:
        lines = [level["line"]] + t[1:7]
    else:
        lines = committed_input(level)[:4]
    lines = lines + ([MP2_LINE] if level["mp2"] else []) + ([PC_LINE] if charged else []) + ["* xyzfile %s substrate.xyz" % cm, ""]
    return "\n".join(lines)


def runs_of(sel):
    """every run, in the order it is made: (run, level key, design or None, kind, state, charge-and-multiplicity, stab)"""
    out = []
    for lv in sel["levels"]:
        for kind, cm, stab in kinds(lv):
            for st in STATES:
                out.append((cname(lv["key"], kind, st), lv["key"], None, kind, st, cm, stab))
        for d in sel["designs"]:
            for kind, cm, stab in kinds(lv):
                for st in STATES:
                    out.append((rname(lv["key"], d["name"], kind, st), lv["key"], d["name"], kind, st, cm, stab))
    if len(out) != NRUN or len(set(r[0] for r in out)) != NRUN:
        stop("%d runs were listed, not %d distinct ones." % (len(out), NRUN))
    for lv in sel["levels"]:
        if sum(1 for r in out if r[1] == lv["key"]) != NPER[lv["dh"]]:
            stop("%s does not have %d runs." % (lv["key"], NPER[lv["dh"]]))
    return out


def conditions_text(designs):
    """d44/CONDITIONS.tsv"""
    L = ["design\tn\tq_e\tsum_q_e\tsum_abs_q_e\tsites\tcharges\tpoint_charge_file\tsha256"]
    for d in designs:
        L.append("%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (d["name"], d["n"], d["q"], d["sum"], d["total"], d["sites"], d["charges"], d["pc"], d["pcsha"]))
    return "\n".join(L) + "\n"


def levels_text(levels):
    """d44/LEVELS.tsv"""
    L = ["key\tlevel\tdouble_hybrid\tmethod_line\tmp2_block\tmethod_line_from\tkinds_of_run\tjob\truns"]
    for l in levels:
        L.append("%s\t%s\t%s\t%s\t%s\t%s_R/job.inp and %s_TS/job.inp\t%s\t%s\t%d" % (
            l["key"], l["name"], "yes" if l["dh"] else "no", l["line"], MP2_LINE if l["mp2"] else "", l["src"], l["src"],
            ", ".join("%s %s" % (k, "with the stability analysis" if s else "without it") for k, cm, s in kinds(l)), l["job"], NPER[l["dh"]]))
    return "\n".join(L) + "\n"


def pbs_text(level, runs, hpc):
    """the job script of one level"""
    k = level["key"]
    mine = [r for r in runs if r[1] == k]
    L = ["#!/bin/bash", "#PBS -N %s" % level["job"],
         "#PBS -l select=1:ncpus=8:mem=%s" % MEM, "#PBS -l walltime=%s" % WALL,
         "#PBS -m ae", "#PBS -M 18660916@sun.ac.za", "#PBS -j oe",
         "#PBS -o %s/run_%s.pbs.out" % (hpc, k),
         "# Stage 6 DIAGNOSTIC 44: the cross-check of AMENDMENT 7 on the six survivors of D43.",
         "# This job is one of five, one for each level. Its level is %s." % level["name"]]
    if level["dh"]:
        L += ["# Order: the six controls with no charges: the dianion without the stability analysis (C_ctl),",
              "#   the dianion with it (C_ctlS) and the monoanion (C_bare), each at the reactant and at the",
              "#   transition state; then each of the six designs as di R, di TS, ds R, ds TS, mo R, mo TS."]
    else:
        L += ["# Order: the four controls with no charges: the dianion (C_ctl) and the monoanion (C_bare) at",
              "#   the reactant and at the transition state; then each of the six designs as di R, di TS,",
              "#   mo R, mo TS."]
    L += ["# A run that fails does not stop the ones after it.",
          "# No -q: the allocation is automatic and walltime %s routes to day. Restartable:" % WALL,
          "#   submitted again, it skips every run that ended normally.",
          'export PATH="/apps/openmpi/4.1.1/bin:$PATH"',
          'export LD_LIBRARY_PATH="/home/apps2/ORCA/6.0.1/lib:/apps/openmpi/4.1.1/lib:'
          '/apps/mambaforge/envs/medaka/lib:${LD_LIBRARY_PATH:-}"',
          "export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1",
          "ORCA=/home/apps2/ORCA/6.0.1",
          "W=%s" % hpc, 'cd "$W" || exit 1',
          'if ! mkdir .running_%s 2>/dev/null; then echo "FAIL: another job holds .running_%s in $W"; exit 1; fi' % (k, k),
          "trap 'rmdir \"$W/.running_%s\" 2>/dev/null' EXIT" % k,
          'echo "host=$(hostname) start=$(date)"',
          'echo "free memory at start:"; free -g | head -2',
          "run() {",
          '  d="$W/$1"',
          '  if [ -s "$d/job.out" ] && grep -q "ORCA TERMINATED NORMALLY" "$d/job.out"; then echo "$1 done before"',
          '  else ( cd "$d" && "$ORCA/orca" job.inp > job.out 2>&1 </dev/null ); '
          'echo "$1 orca_exit=$? normal=$(grep -c \'ORCA TERMINATED NORMALLY\' "$d/job.out") $(date)"; fi',
          "}"]
    L += ["run %s" % r[0] for r in mine]
    L += ['echo "end=$(date)"',
          'echo "finished: $(grep -l \'ORCA TERMINATED NORMALLY\' %s_*/job.out 2>/dev/null | wc -l) of %d"' % (k, len(mine))]
    return "\n".join(L) + "\n"


def main():
    if len(sys.argv) != 4:
        stop("usage: cs6_d44_prepare.py <results-repo> <out-dir> <hpc-dir>")
    res, out, hpc = sys.argv[1:4]
    s6 = os.path.join(res, S6)
    if os.path.exists(out) and os.listdir(out):
        stop("%s exists and is not empty; this script writes into a new directory only." % out)
    sel = select(s6)
    levels, designs = sel["levels"], sel["designs"]
    lmap = dict((l["key"], l) for l in levels)
    dmap = dict((d["name"], d) for d in designs)
    runs = runs_of(sel)
    for f in sources([d["name"] for d in designs]):
        if not os.path.isfile(os.path.join(s6, f)):
            stop("%s, which this diagnostic reads, is missing." % f)

    os.makedirs(out, exist_ok=True)
    sums = []
    for run, lk, dn, kind, st, cm, stab in runs:
        d = os.path.join(out, run)
        os.makedirs(d)
        shutil.copyfile(os.path.join(s6, GEOM[st]), os.path.join(d, "substrate.xyz"))
        with open(os.path.join(d, "job.inp"), "w") as f:
            f.write(build_input(lmap[lk], cm, dn is not None, stab))
        if dn is not None:
            shutil.copyfile(os.path.join(s6, dmap[dn]["pc"]), os.path.join(d, "pc.pc"))
        for fn in sorted(os.listdir(d)):
            sums.append((sha(os.path.join(d, fn)), "%s/%s" % (run, fn)))

    # ---- the checks, on what was actually written, against the committed files and not against build_input
    tl = sel["tl"]
    for run, lk, dn, kind, st, cm, stab in runs:
        d = os.path.join(out, run)
        lv = lmap[lk]
        if sha(os.path.join(d, "substrate.xyz")) != sha(os.path.join(s6, GEOM[st])):
            stop("%s/substrate.xyz is not frame 42436's committed %s geometry." % (run, st))
        want = ["job.inp", "pc.pc", "substrate.xyz"] if dn is not None else ["job.inp", "substrate.xyz"]
        if sorted(os.listdir(d)) != want:
            stop("%s does not hold exactly %s." % (run, ", ".join(want)))
        b = open(os.path.join(d, "job.inp")).read().split("\n")
        psrc = os.path.join(s6, "%s_%s" % (lv["src"], st), "job.inp")
        src = open(psrc).read().split("\n")
        if b[0] != src[0]:
            stop("%s: the method line is not line 1 of %s_%s/job.inp." % (run, lv["src"], st))
        if stab:
            # the working level's committed charged input, with: the method line; the %mp2 block where the
            # level's committed input holds one; no %pointcharges line in a control; the charge-and-multiplicity line
            exp = list(tl)
            exp[0] = src[0]
            if dn is None:
                exp.remove(PC_LINE)
            if MP2_LINE in src:
                exp.insert(7, MP2_LINE)
            exp[-2] = "* xyzfile %s substrate.xyz" % cm
            if b != exp:
                stop("%s is not the working level's input with the method line of %s and nothing else changed." % (run, lv["name"]))
        else:
            # the committed charge-free input of the level, with the %pointcharges line before the geometry line
            if cm != DI:
                stop("%s: only a dianion is run without the stability analysis." % run)
            if dn is None:
                if sha(os.path.join(d, "job.inp")) != sha(psrc):
                    stop("%s is not %s_%s/job.inp byte for byte." % (run, lv["src"], st))
            else:
                exp = list(src)
                exp.insert(len(exp) - 2, PC_LINE)
                if b != exp:
                    stop("%s is not %s_%s/job.inp with the %%pointcharges line added and nothing else." % (run, lv["src"], st))
        if b.count(MP2_LINE) != (1 if MP2_LINE in src else 0) or b.count(PC_LINE) != (1 if dn is not None else 0) \
                or ("STABPerform true" in "\n".join(b)) != stab:
            stop("%s does not carry the %%mp2 block, the point charges and the stability analysis exactly where it should." % run)
        if dn is not None and sha(os.path.join(d, "pc.pc")) != dmap[dn]["pcsha"]:
            stop("%s/pc.pc is not D43's committed point-charge file of %s." % (run, dn))

    with open(os.path.join(out, "CONDITIONS.tsv"), "w") as f:
        f.write(conditions_text(designs))
    with open(os.path.join(out, "LEVELS.tsv"), "w") as f:
        f.write(levels_text(levels))
    with open(os.path.join(out, "SOURCES.sha256"), "w") as f:
        for s in sources([d["name"] for d in designs]):
            f.write("%s  %s\n" % (sha(os.path.join(s6, s)), s))
    for fn in ("CONDITIONS.tsv", "LEVELS.tsv", "SOURCES.sha256"):
        sums.append((sha(os.path.join(out, fn)), fn))
    for l in levels:
        p = os.path.join(out, "run_%s.pbs" % l["key"])
        with open(p, "w") as f:
            f.write(pbs_text(l, runs, hpc))
        os.chmod(p, 0o755)
        sums.append((sha(p), "run_%s.pbs" % l["key"]))
    if len(sums) != NSUM:
        stop("%d files were written, not %d." % (len(sums), NSUM))
    with open(os.path.join(out, "INPUT_CHECKSUMS.sha256"), "w") as f:
        for h, n in sums:
            f.write("%s  %s\n" % (h, n))
    for l in levels:
        print("   %-5s  %-28s  %2d runs  %s%s" % (l["key"], l["name"], NPER[l["dh"]], l["line"], "   + " + MP2_LINE if l["mp2"] else ""))
    for d in designs:
        print("   %-13s n %2d  q %s  sum %s e  total %s e  sites %s" % (d["name"], d["n"], d["q"], d["sum"], d["total"], d["sites"]))
    print("   %d run directories, %d files checksummed, five jobs, %s, walltime %s each" % (NRUN, len(sums), MEM, WALL))


if __name__ == "__main__":
    try:
        main()
    except Refused as e:
        sys.stderr.write("REFUSED  %s\n" % e)
        raise SystemExit(2)
