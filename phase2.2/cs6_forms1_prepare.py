#!/usr/bin/env python3
"""cs6_forms1_prepare.py - form the inputs of STAGE 6, FORMS TEST 1 from committed files.

AMENDMENT 13, A13.4 item 9: "A short test of the new forms of input on hpc1, which reads no number
of a design, comes first." This program writes those inputs. It computes nothing and reads no
output. Every geometry, site, wall block and method line is taken from a committed file; the only
numbers written here are the block of the published capping potential (AMENDMENT 12, A12.3 item 5:
N_core 0, lmax f, s 1.732 / 14.676, p 1.115 / 5.1757, d 1.203 / -1.816, f zero) and the caps.

usage:
   python3 cs6_forms1_prepare.py --res <results repository> --out <folder to write forms1 into>
The reader cs6_forms1.py calls build() and stops unless the committed forms1/ is what build() gives.
"""
import argparse, hashlib, os, sys

S6REL = "phase2.2/cs_stage6"
DIR = "forms1"
POS, NEG = 1954, 1612          # the two finite-difference sites of A_REPORT.txt; NOT site 1690 of D3
HPC = "/home/18660916/system_development/phase2.2/cs_stage6/" + DIR
SHI = ["  N_core 0", "  lmax f", "  s 1", "    1  1.732000000  14.676000000  2", "  p 1", "    1  1.115000000  5.175700000  2",
       "  d 1", "    1  1.203000000  -1.816000000  2", "  f 1", "    1  1.000000000  0.000000000  2", "end"]

# run, job, cap in seconds, level, what form it tests
RUNS = [
    ("SV_wallpos_R", "F1", 1200, "SV", "helium wall on a positive charge, a bare negative charge beside it (%pointcharges)"),
    ("SV_shi_R", "F1", 1200, "SV", "the published capping block (N_core 0) on a coreless centre that carries +1 e, written Na> as D57 will write it"),
    ("SV_Na_R", "F1", 1200, "SV", "Na+ as an atom of the calculation"),
    ("SV_K_R", "F1", 1200, "SV", "K+ as an atom of the calculation"),
    ("SV_ghostNa_R", "F1", 1200, "SV", "the functions of Na with no nucleus and no electrons"),
    ("SV_ghostK_R", "F1", 1200, "SV", "the functions of K with no nucleus and no electrons"),
    ("SV_engrad_R", "F1", 2400, "SV", "energy and gradient with %pointcharges (expects job.engrad and job.pcgrad)"),
    ("SV_mono_R", "F1", 1800, "SV", "the control with one electron removed, as D56 ran it, then orca_vpot at five sites"),
    ("TZ_wallHe_R", "F2", 3000, "TZ", "helium wall carrying +1 e at def2-TZVPD"),
    ("TZ_wallNe_R", "F2", 3000, "TZ", "neon wall carrying +1 e at def2-TZVPD"),
    ("TZ_shi_R", "F2", 3000, "TZ", "the published capping block at def2-TZVPD"),
    ("TZ_Na_R", "F3", 3000, "TZ", "Na+ as an atom at def2-TZVPD"),
    ("TZ_K_R", "F3", 3000, "TZ", "K+ as an atom at def2-TZVPD"),
    ("TZ_ghostNa_R", "F3", 3000, "TZ", "the functions of Na alone at def2-TZVPD"),
    ("TZ_ghostK_R", "F3", 3000, "TZ", "the functions of K alone at def2-TZVPD"),
]
JOBS = ("F1", "F2", "F3")
BASIS = {"SV": "def2-SVPD", "TZ": "def2-TZVPD"}


def rd(p):
    with open(p) as f:
        return f.read()


def head(level, extra="", stab=True):
    scf = "%scf\n  MaxIter 300\n  STABPerform true\nend\n" if stab else "%scf MaxIter 300 end\n"
    return "! CAM-B3LYP D3BJ %s def2/J RIJCOSX TightSCF%s\n%%maxcore 3000\n%%pal nprocs 8 end\n%s" % (BASIS[level], extra, scf)


def wall_block(path, tag):
    """the lines of a committed wall block, from the line after '<tag>' to its 'end'"""
    L = rd(path).split("\n")
    i = [k for k, l in enumerate(L) if l.split()[:1] == [tag]][0]
    j = i + 1
    while L[j].strip() != "end":
        j += 1
    return L[i + 1:j + 1]


def build(S6):
    """every file of forms1/, as {relative path: text}"""
    sub = rd(os.path.join(S6, "d56/SV_C_ctl_R/substrate.xyz"))
    atoms = [l for l in sub.split("\n")[2:] if l.strip()]
    assert int(sub.split("\n")[0]) == 24 and len(atoms) == 24
    site = {}
    for l in rd(os.path.join(S6, "sites_f42436.tsv")).split("\n"):
        c = l.split("\t")
        if c[0] in (str(POS), str(NEG)):
            site[int(c[0])] = tuple(float(x) for x in c[2:5])
    assert set(site) == {POS, NEG}
    P, N = site[POS], site[NEG]
    he = wall_block(os.path.join(S6, "d23/he_d16_S6-10_wall0_R/job.inp"), "He>")
    ne = wall_block(os.path.join(S6, "additivity/def2-SVP/n10_wall0_R/job.inp"), "Ne>")
    pts = rd(os.path.join(S6, "d33a/points.xyz")).split("\n")
    mono = rd(os.path.join(S6, "d56/SV_C_bare_R/job.inp"))
    assert rd(os.path.join(S6, "d56/SV_C_bare_R/substrate.xyz")) == sub
    xyz = lambda p: "%16.10f  %18.10f  %18.10f" % p
    pcl = lambda q, p: " %9.6f  %s\n" % (q, xyz(p))

    def centre(tag, q, block):
        return "* xyz -2 1\n" + "".join("%s\n" % a for a in atoms) + "  %s  %.6f  %s  NewECP\n" % (tag, q, xyz(P)) + "".join("%s\n" % b for b in block) + "*\n"

    def system(sym):
        return "25\nsubstrate of frame 42436 at R, then %s at site %d\n" % (sym, POS) + "".join("%s\n" % a for a in atoms) + "  %-3s %s\n" % (sym, xyz(P))

    F = {}
    F["SV_wallpos_R/job.inp"] = head("SV") + '%pointcharges "pc.pc"\n' + centre("He>", 0.25, he)
    F["SV_wallpos_R/pc.pc"] = "1\n" + pcl(-0.25, N)
    for lv in ("SV", "TZ"):
        F["%s_shi_R/job.inp" % lv] = head(lv) + centre("Na>", 1.0, SHI)
        for sym, name, chg in (("Na", "Na", -1), ("K", "K", -1), ("Na:", "ghostNa", -2), ("K:", "ghostK", -2)):
            F["%s_%s_R/job.inp" % (lv, name)] = head(lv) + "* xyzfile %d 1 system.xyz\n" % chg
            F["%s_%s_R/system.xyz" % (lv, name)] = system(sym)
    F["TZ_wallHe_R/job.inp"] = head("TZ") + centre("He>", 1.0, he)
    F["TZ_wallNe_R/job.inp"] = head("TZ") + centre("Ne>", 1.0, ne)
    F["SV_engrad_R/job.inp"] = head("SV", " EnGrad", stab=False) + '%pointcharges "pc.pc"\n* xyzfile -2 1 substrate.xyz\n'
    F["SV_engrad_R/pc.pc"] = "2\n" + pcl(0.25, P) + pcl(-0.25, N)
    F["SV_engrad_R/substrate.xyz"] = sub
    F["SV_mono_R/job.inp"] = mono
    F["SV_mono_R/substrate.xyz"] = sub
    F["points5.xyz"] = "5\n" + "".join("%s\n" % l for l in pts[1:6])
    F["RUNS.tsv"] = "run\tjob\tcap_s\tlevel\tform\n" + "".join("%s\t%s\t%d\t%s\t%s\n" % r for r in RUNS)
    for j in JOBS:
        F["run_%s.pbs" % j] = pbs(j)
    assert sorted(set(p.split("/")[0] for p in F if "/" in p)) == sorted(r[0] for r in RUNS)
    F["INPUT_CHECKSUMS.sha256"] = "".join("%s  %s\n" % (hashlib.sha256(F[p].encode()).hexdigest(), p) for p in sorted(F))
    return F


def pbs(j):
    mine = [r for r in RUNS if r[1] == j]
    t = """#!/bin/bash
#PBS -N cm_forms1_%s
#PBS -l select=1:ncpus=8:mem=48gb
#PBS -l walltime=06:00:00
#PBS -m ae
#PBS -M 18660916@sun.ac.za
#PBS -j oe
#PBS -o %s/run_%s.pbs.out
# Stage 6 FORMS TEST 1: does ORCA accept each new form of input (AMENDMENT 13, A13.4 item 9).
# Job %s of 3. It makes %d runs in the order below (name, cap in s); a run that fails does not stop the ones after it.
# No -q: the allocation is automatic. The job script is D56's committed run_C0.pbs with its names, runs and caps
# changed, and for F1 the potential step of D33's committed run.pbs. No run is made twice:
#   submitted again, it makes only the runs whose folder holds no output yet (D44 NOTE 3).
export PATH="/apps/openmpi/4.1.1/bin:$PATH"
export LD_LIBRARY_PATH="/home/apps2/ORCA/6.0.1/lib:/apps/openmpi/4.1.1/lib:/apps/mambaforge/envs/medaka/lib:${LD_LIBRARY_PATH:-}"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
ORCA=/home/apps2/ORCA/6.0.1
W=%s
cd "$W" || exit 1
if ! mkdir .running_%s 2>/dev/null; then echo "FAIL: another job holds .running_%s in $W"; exit 1; fi
trap 'rmdir "$W/.running_%s" 2>/dev/null' EXIT
echo "host=$(hostname) start=$(date)"
echo "free memory at start:"; free -g | head -2
run() {
  d="$W/$1"
  if [ -e "$d/job.out" ]; then echo "$1 has an output already: not run again"
  else ( cd "$d" && timeout -k 10 "$2" "$ORCA/orca" job.inp > job.out 2>&1 </dev/null ); echo "$1 orca_exit=$? normal=$(grep -c 'ORCA TERMINATED NORMALLY' "$d/job.out") $(date)"; fi
}
""" % (j.lower(), HPC, j, j, len(mine), HPC, j, j, j)
    if j == "F1":
        t += """vpot() {
  d="$W/$1"; pts="$W/$2"
  if [ -e "$d/vpot.log" ]; then echo "$1 vpot has a log already: not run again"; return; fi
  if ! grep -q "ORCA TERMINATED NORMALLY" "$d/job.out" 2>/dev/null; then echo "$1 vpot skipped, run not normal"; return; fi
  ( cd "$d" && "$ORCA/orca_vpot" job.gbw job.scfp "$pts" vpot.out job > vpot.log 2>&1 )
  echo "$1 vpot_exit=$? lines=$(wc -l < "$d/vpot.out" 2>/dev/null || echo 0) expected=6"
}
"""
    for r in mine:
        t += "run %s %d\n" % (r[0], r[2])
    if j == "F1":
        t += "vpot SV_mono_R points5.xyz\n"
    t += 'echo "end=$(date)"\n'
    t += "n=0; for r in %s; do grep -q 'ORCA TERMINATED NORMALLY' \"$W/$r/job.out\" 2>/dev/null && n=$((n+1)); done\n" % " ".join(r[0] for r in mine)
    t += 'echo "finished: $n of %d"\n' % len(mine)
    return t


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--res", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    F = build(os.path.join(a.res, S6REL))
    dst = os.path.join(a.out, DIR)
    if os.path.exists(dst):
        sys.exit("REFUSED: %s exists already." % dst)
    for p in sorted(F):
        os.makedirs(os.path.dirname(os.path.join(dst, p)), exist_ok=True)
        with open(os.path.join(dst, p), "w") as f:
            f.write(F[p])
    print("wrote %d files into %s: %d runs in %d jobs" % (len(F), dst, len(RUNS), len(JOBS)))


if __name__ == "__main__":
    main()
