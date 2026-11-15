#!/usr/bin/env python3
"""cs6_d42_prepare.py - DIAGNOSTIC 42: is the dianion bound with the charges of DIAGNOSTIC 31's RELIABLE conditions?

27 single points on frame 42436 at CAM-B3LYP-D3BJ/def2-SVPD with the stability analysis.

  the control   D40's committed monoanion input for the bare substrate at the reactant, d40/B_bare_R, byte for byte
  26 monoanions the 13 conditions of D31's arm B that D31's result calls RELIABLE, at the reactant and at the
                transition state: each is the COMMITTED dianion input of d31/B/<condition>_<state> with its
                charge-and-multiplicity line changed from "-2 1" to "-1 2" and nothing else touched, beside
                byte-for-byte copies of that run's own geometry and point-charge file

The 13 conditions are not typed in as a list to be trusted: they are read from the table R31.2 of D31's committed
result, and this script refuses to finish unless the list below is exactly the set that table calls RELIABLE.
It also writes CONDITIONS.tsv: each condition's n, q, the sum of its charges as its committed pc.pc holds them,
and the class the pre-registration reads it in (H: the sum is zero or more; N: the sum is below zero).

Usage:  cs6_d42_prepare.py <results-repo> <out-dir> <hpc-dir>
"""
import hashlib, os, re, shutil, sys

S6 = "phase2.2/cs_stage6"
GEOM = {"R": "d33a/sp_42436_R/substrate.xyz", "TS": "d33a/sp_42436_TS/substrate.xyz"}
DI, MONO = "-2 1", "-1 2"
D31R = "STAGE6_DIAGNOSTIC31_RESULT_20261107.txt"
METHOD = "! CAM-B3LYP D3BJ def2-SVPD def2/J RIJCOSX TightSCF"
STATES = ("R", "TS")
CONTROL = ("C_bare_R", "d40/B_bare_R", "R", ["substrate.xyz"])
# (n, q) of the conditions of D31's arm B that R31.2 calls RELIABLE; checked against that table below
RELIABLE = [(1, "0.25"), (2, "0.25"), (3, "0.25"), (5, "0.25"), (10, "0.25"), (20, "0.25"),
            (1, "0.50"), (2, "0.50"), (3, "0.50"), (5, "0.50"),
            (1, "0.75"), (2, "0.75"), (3, "0.75")]
WALL, MEM = "24:00:00", "48gb"


def stop(m):
    sys.stderr.write("REFUSED  %s\n" % m)
    raise SystemExit(2)


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 16), b""):
            h.update(b)
    return h.hexdigest()


def cname(n, q):
    return "n%02d_q%s" % (n, q)


def charges(p):
    ls = open(p).read().split("\n")
    n = int(ls[0].split()[0])
    q = [float(l.split()[0]) for l in ls[1:1 + n]]
    if len(q) != n:
        stop("%s does not hold the %d charges its first line announces." % (p, n))
    return q


def reliable_in_r312(text):
    """The (n, q) cells that the table of R31.2 calls RELIABLE."""
    i = text.find("R31.2")
    j = text.find("R31.3")
    if i < 0 or j < 0:
        stop("%s does not hold R31.2 and R31.3." % D31R)
    block = text[i:j]
    head = [l for l in block.splitlines() if re.match(r"^\s+n\s+q=", l)]
    if len(head) != 1:
        stop("R31.2 of %s does not hold exactly one heading row." % D31R)
    qs = re.findall(r"q=(\d\.\d\d)", head[0])
    out = []
    for l in block.splitlines():
        m = re.match(r"^\s+(\d+)\s+((?:(?:RELIABLE|UNRELIABLE|INCONCLUSIVE)\s+){%d})" % len(qs), l)
        if m:
            for q, v in zip(qs, m.group(2).split()):
                if v == "RELIABLE":
                    out.append((int(m.group(1)), q))
    return out


def main():
    if len(sys.argv) != 4:
        stop("usage: cs6_d42_prepare.py <results-repo> <out-dir> <hpc-dir>")
    res, out, hpc = sys.argv[1], sys.argv[2], sys.argv[3]
    s6 = os.path.join(res, S6)
    if os.path.exists(out) and os.listdir(out):
        stop("%s exists and is not empty; this script writes into a new directory only." % out)

    # ---- the 13 conditions are the ones D31's committed result calls RELIABLE, and no others
    p31 = os.path.join(s6, D31R)
    if not os.path.isfile(p31):
        stop("%s is missing." % D31R)
    found = reliable_in_r312(open(p31).read())
    if sorted(found) != sorted(RELIABLE) or len(set(RELIABLE)) != len(RELIABLE):
        stop("the conditions listed here are not exactly those R31.2 calls RELIABLE: listed %d, the table has %d."
             % (len(RELIABLE), len(found)))

    runs = [CONTROL]
    cond = []
    for n, q in RELIABLE:
        c = cname(n, q)
        qs = {st: charges(os.path.join(s6, "d31", "B", "%s_%s" % (c, st), "pc.pc")) for st in STATES}
        if qs["R"] != qs["TS"]:
            stop("%s does not carry the same charges at the two states." % c)
        if len(qs["R"]) != n or any(abs(abs(x) - float(q)) > 1e-9 for x in qs["R"]):
            stop("%s does not hold %d charges of magnitude %s e." % (c, n, q))
        s = sum(qs["R"])
        cond.append((c, n, q, sum(1 for x in qs["R"] if x > 0), sum(1 for x in qs["R"] if x < 0), s, "H" if s > -5e-7 else "N"))
    # the order of running: the control; the conditions whose charges sum to zero or more, which are the test; then the rest
    for cls in ("H", "N"):
        for c, n, q, npos, nneg, s, k in cond:
            if k == cls:
                runs += [("%s_%s" % (c, st), "d31/B/%s_%s" % (c, st), st, ["substrate.xyz", "pc.pc"]) for st in STATES]
    nrun = 1 + 2 * len(RELIABLE)
    if len(runs) != nrun or len(set(r[0] for r in runs)) != nrun:
        stop("%d runs were listed, expected %d." % (len(runs), nrun))

    # ---- every source, before anything is written
    for name, src, st, extra in runs:
        d = os.path.join(s6, src)
        for fn in ["job.inp", "job.out"] + extra:
            if not os.path.isfile(os.path.join(d, fn)):
                stop("%s/%s is missing; %s is built from it." % (src, fn, name))
        t = open(os.path.join(d, "job.inp")).read()
        want = MONO if name == CONTROL[0] else DI
        if sum(1 for l in t.split("\n") if l.startswith("* xyzfile ")) != 1 or ("* xyzfile %s " % want) not in t:
            stop("%s/job.inp does not have exactly one geometry line at charge and multiplicity %s." % (src, want))
        if "****ORCA TERMINATED NORMALLY****" not in open(os.path.join(d, "job.out"), errors="replace").read():
            stop("%s/job.out did not terminate normally." % src)
        if sha(os.path.join(d, "substrate.xyz")) != sha(os.path.join(s6, GEOM[st])):
            stop("%s/substrate.xyz is not frame 42436's committed %s geometry." % (src, st))

    os.makedirs(out, exist_ok=True)
    sums, first = [], {}
    for name, src, st, extra in runs:
        d = os.path.join(out, name)
        os.makedirs(d)
        sd = os.path.join(s6, src)
        for fn in extra:
            shutil.copyfile(os.path.join(sd, fn), os.path.join(d, fn))
            if sha(os.path.join(sd, fn)) != sha(os.path.join(d, fn)):
                stop("%s/%s did not copy byte for byte into %s" % (src, fn, name))
        dst = os.path.join(d, "job.inp")
        if name == CONTROL[0]:
            shutil.copyfile(os.path.join(sd, "job.inp"), dst)
        else:
            lines = open(os.path.join(sd, "job.inp")).read().split("\n")
            new = [l.replace("* xyzfile %s " % DI, "* xyzfile %s " % MONO) if l.startswith("* xyzfile ") else l for l in lines]
            with open(dst, "w") as f:
                f.write("\n".join(new))
        for fn in sorted(os.listdir(d)):
            sums.append((sha(os.path.join(d, fn)), "%s/%s" % (name, fn)))
        first[name] = open(dst).read().splitlines()[0]

    # ---- the comparability check, on what was actually written
    for name, src, st, extra in runs:
        a = open(os.path.join(s6, src, "job.inp")).read().split("\n")
        b = open(os.path.join(out, name, "job.inp")).read().split("\n")
        diff = [(x, y) for x, y in zip(a, b) if x != y]
        if len(a) != len(b):
            stop("%s: the input has a different number of lines from its committed source." % name)
        if name == CONTROL[0]:
            if diff or sha(os.path.join(s6, src, "job.inp")) != sha(os.path.join(out, name, "job.inp")):
                stop("%s is not D40's committed monoanion input byte for byte." % name)
        else:
            if len(diff) != 1:
                stop("%s differs from %s/job.inp in %d lines; exactly one is allowed." % (name, src, len(diff)))
            x, y = diff[0]
            if not (x.startswith("* xyzfile %s " % DI) and y == x.replace(DI, MONO, 1)):
                stop("%s: the one changed line is not the charge-and-multiplicity line: %r -> %r" % (name, x, y))
        if sorted(os.listdir(os.path.join(out, name))) != sorted(["job.inp"] + extra):
            stop("%s holds files other than its input and the %d files copied beside it." % (name, len(extra)))
        t = open(os.path.join(out, name, "job.inp")).read()
        if first[name] != METHOD or "STABPerform true" not in t or ("* xyzfile %s substrate.xyz" % MONO) not in t:
            stop("%s is not a CAM-B3LYP-D3BJ/def2-SVPD monoanion input with the stability analysis." % name)
        if ('%pointcharges "pc.pc"' in t) != ("pc.pc" in extra):
            stop("%s: the point-charge line and the point-charge file do not go together." % name)

    with open(os.path.join(out, "CONDITIONS.tsv"), "w") as f:
        f.write("condition\tn\tq_e\tn_positive\tn_negative\tsum_q_e\tclass\tsource_R\tsource_TS\n")
        for c, n, q, npos, nneg, s, k in cond:
            f.write("%s\t%d\t%s\t%d\t%d\t%+.4f\t%s\td31/B/%s_R\td31/B/%s_TS\n" % (c, n, q, npos, nneg, s, k, c, c))
    sums.append((sha(os.path.join(out, "CONDITIONS.tsv")), "CONDITIONS.tsv"))

    pbs = os.path.join(out, "run.pbs")
    with open(pbs, "w") as f:
        f.write("#!/bin/bash\n#PBS -N cm_d42\n")
        f.write("#PBS -l select=1:ncpus=8:mem=%s\n#PBS -l walltime=%s\n" % (MEM, WALL))
        f.write("#PBS -m ae\n#PBS -M 18660916@sun.ac.za\n#PBS -j oe\n")
        f.write("#PBS -o %s/run.pbs.out\n" % hpc)
        f.write("# Stage 6 DIAGNOSTIC 42: is the dianion bound with the charges of D31's RELIABLE conditions?\n")
        f.write("# C_bare_R is D40's committed monoanion input for the bare substrate, byte for byte: the CONTROL.\n")
        f.write("# Every other run is a monoanion, charge -1 multiplicity 2: a committed dianion input of D31's\n")
        f.write("#   arm B with that one line changed, beside its own geometry and point-charge file.\n")
        f.write("# Order: the control; the 8 conditions whose charges sum to zero or more, both states; then\n")
        f.write("#   the 5 whose charges sum to less than zero. A run that fails does not stop the ones after it.\n")
        f.write("# No -q: the allocation is automatic and walltime %s routes to day. Restartable.\n" % WALL)
        f.write('export PATH="/apps/openmpi/4.1.1/bin:$PATH"\n')
        f.write('export LD_LIBRARY_PATH="/home/apps2/ORCA/6.0.1/lib:/apps/openmpi/4.1.1/lib:'
                '/apps/mambaforge/envs/medaka/lib:${LD_LIBRARY_PATH:-}"\n')
        f.write("export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1\n")
        f.write("ORCA=/home/apps2/ORCA/6.0.1\n")
        f.write('W=%s\ncd "$W" || exit 1\n' % hpc)
        f.write('if ! mkdir .running 2>/dev/null; then echo "FAIL: another job holds .running in $W"; exit 1; fi\n')
        f.write("trap 'rmdir \"$W/.running\" 2>/dev/null' EXIT\n")
        f.write('echo "host=$(hostname) start=$(date)"\n')
        f.write('echo "free memory at start:"; free -g | head -2\n')
        f.write("run() {\n")
        f.write('  d="$W/$1"\n')
        f.write('  if [ -s "$d/job.out" ] && grep -q "ORCA TERMINATED NORMALLY" "$d/job.out"; then echo "$1 done before"\n')
        f.write('  else ( cd "$d" && "$ORCA/orca" job.inp > job.out 2>&1 </dev/null ); '
                'echo "$1 orca_exit=$? normal=$(grep -c \'ORCA TERMINATED NORMALLY\' "$d/job.out") $(date)"; fi\n')
        f.write("}\n")
        for r in runs:
            f.write("run %s\n" % r[0])
        f.write('echo "end=$(date)"\n')
        f.write('echo "finished: $(grep -l \'ORCA TERMINATED NORMALLY\' */job.out 2>/dev/null | wc -l) of %d"\n' % nrun)
    os.chmod(pbs, 0o755)
    sums.append((sha(pbs), "run.pbs"))

    with open(os.path.join(out, "INPUT_CHECKSUMS.sha256"), "w") as f:
        for h, n in sums:
            f.write("%s  %s\n" % (h, n))
    for c, n, q, npos, nneg, s, k in cond:
        print("   %-10s n %2d  q %s  %2d positive %2d negative  sum %+.4f e  class %s" % (c, n, q, npos, nneg, s, k))
    print("   %d run directories, %d files checksummed, %s, walltime %s" % (nrun, len(sums), MEM, WALL))


if __name__ == "__main__":
    main()
