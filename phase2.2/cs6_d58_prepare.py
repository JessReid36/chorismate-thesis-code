#!/usr/bin/env python3
"""cs6_d58_prepare.py - form the inputs of STAGE 6, DIAGNOSTIC 58 from committed files.

AMENDMENT 13, A13.4 items 2 to 4: the force each survivor of D56 adds to the substrate, measured
directly; the first-order map of that force (F); and the first-order map of the hold on the last
electron (E). This program writes the inputs. It computes nothing and reads no output.

   G runs (38): energy and gradient at def2-SVPD with D56's own point charges, for the 18 survivors
                named by the committed d56_survivors.tsv, at R and at TS; and the substrate without
                charges at R and at TS. The form is FORMS TEST 1's SV_engrad_R.
   M runs (144): the substrate without charges at R with one atom moved by STEP_BOHR along one
                axis, either way; each is followed by orca_vpot at the 2626 sites (d33a/points.xyz).
                The form is D33's committed sp_42436_R with its potential step.
   E runs (4):  D56's committed SV_C_ctl_R, SV_C_ctl_TS, SV_C_bare_R and SV_C_bare_TS inputs,
                unchanged, each followed by orca_vpot at the 2626 sites.

usage:
   python3 cs6_d58_prepare.py --res <results repository> --out <folder to write d58 into>
The reader cs6_d58.py calls build() and stops unless the committed d58/ is what build() gives.
"""
import argparse, csv, hashlib, os, sys

S6REL = "phase2.2/cs_stage6"
DIR = "d58"
HPC = "/home/18660916/system_development/phase2.2/cs_stage6/" + DIR
BOHR = 0.5291772083                    # A per bohr, as cs6_analyse.py has it
STEP_BOHR = 0.01                       # the step of the map F, bohr (fixed in the pre-registration)
CAP = {"G": 3600, "M": 900, "E": 1800}
AXES = "xyz"


def rd(p):
    with open(p) as f:
        return f.read()


def survivors(S6):
    rows = list(csv.DictReader(open(os.path.join(S6, "d56_survivors.tsv")), delimiter="\t"))
    s = [r["design"] for r in rows if r["the_design"] == "SURVIVOR"]
    assert len(s) == 18, len(s)
    return s


def runs(S6):
    """(run, job, cap) in the order the jobs make them"""
    G = ["G_C_R", "G_C_TS"] + ["G_%s_%s" % (d, st) for d in survivors(S6) for st in ("R", "TS")]
    jobs = [("G1", G[0:10]), ("G2", G[10:20]), ("G3", G[20:30]), ("G4", G[30:38])]
    M = [["M_a%02d_%s_%s" % (a, ax, sg) for a in range(lo, lo + 8) for ax in AXES for sg in ("p", "m")] for lo in (0, 8, 16)]
    jobs += [("M1", M[0]), ("M2", M[1]), ("M3", M[2])]
    jobs += [("E1", ["E_ctl_R", "E_ctl_TS", "E_bare_R", "E_bare_TS"])]
    out = [(r, j, CAP[r[0]]) for j, rs in jobs for r in rs]
    assert len(out) == 186 and len(set(r for r, _, _ in out)) == 186
    return out, [j for j, _ in jobs]


def head(engrad):
    if engrad:
        return "! CAM-B3LYP D3BJ def2-SVPD def2/J RIJCOSX TightSCF EnGrad\n%maxcore 3000\n%pal nprocs 8 end\n%scf MaxIter 300 end\n"
    return "! CAM-B3LYP D3BJ def2-SVPD def2/J RIJCOSX TightSCF\n%maxcore 3000\n%pal nprocs 8 end\n%scf MaxIter 300 end\n"


def displaced(sub, atom, axis, sign):
    L = sub.split("\n")
    n = int(L[0])
    assert n == 24
    out = [L[0], "d58: the substrate of d56/SV_C_ctl_R with atom %d moved %s%.2f bohr along %s" % (atom, "+" if sign > 0 else "-", STEP_BOHR, axis)]
    for k in range(n):
        c = L[2 + k].split()
        x = [float(v) for v in c[1:4]]
        if k == atom:
            x[AXES.index(axis)] += sign * STEP_BOHR * BOHR
        out.append("  %-2s %24.14f %21.14f %21.14f" % (c[0], x[0], x[1], x[2]))
    return "\n".join(out) + "\n"


def build(S6):
    """every file of d58/, as {relative path: text}"""
    R, J = runs(S6)
    d56 = os.path.join(S6, "d56")
    subR, subTS = rd(os.path.join(d56, "SV_C_ctl_R/substrate.xyz")), rd(os.path.join(d56, "SV_C_ctl_TS/substrate.xyz"))
    F = {}
    for r, j, c in R:
        if r.startswith("G_"):
            st = r.rsplit("_", 1)[1]
            sub = subR if st == "R" else subTS
            if r.startswith("G_C_"):
                F[r + "/job.inp"] = head(True) + "* xyzfile -2 1 substrate.xyz\n"
            else:
                d = r[2:].rsplit("_", 1)[0]
                src = os.path.join(d56, "SV_%s_di_%s" % (d, st))
                assert rd(os.path.join(src, "substrate.xyz")) == sub
                F[r + "/job.inp"] = head(True) + '%pointcharges "pc.pc"\n* xyzfile -2 1 substrate.xyz\n'
                F[r + "/pc.pc"] = rd(os.path.join(src, "pc.pc"))
            F[r + "/substrate.xyz"] = sub
        elif r.startswith("M_"):
            a, ax, sg = r.split("_")[1:]
            F[r + "/job.inp"] = head(False) + "* xyzfile -2 1 substrate.xyz\n"
            F[r + "/substrate.xyz"] = displaced(subR, int(a[1:]), ax, 1 if sg == "p" else -1)
        else:
            src = os.path.join(d56, {"E_ctl_R": "SV_C_ctl_R", "E_ctl_TS": "SV_C_ctl_TS", "E_bare_R": "SV_C_bare_R", "E_bare_TS": "SV_C_bare_TS"}[r])
            F[r + "/job.inp"] = rd(os.path.join(src, "job.inp"))
            F[r + "/substrate.xyz"] = rd(os.path.join(src, "substrate.xyz"))
    F["points.xyz"] = rd(os.path.join(S6, "d33a/points.xyz"))
    F["RUNS.tsv"] = "run\tjob\tcap_s\n" + "".join("%s\t%s\t%d\n" % x for x in R)
    for j in J:
        F["run_%s.pbs" % j] = pbs(j, [(r, c) for r, jj, c in R if jj == j])
    F["INPUT_CHECKSUMS.sha256"] = "".join("%s  %s\n" % (hashlib.sha256(F[p].encode()).hexdigest(), p) for p in sorted(F))
    return F


def pbs(j, mine):
    pot = j[0] in "ME"
    t = """#!/bin/bash
#PBS -N cm_d58_%s
#PBS -l select=1:ncpus=8:mem=48gb
#PBS -l walltime=24:00:00
#PBS -m ae
#PBS -M 18660916@sun.ac.za
#PBS -j oe
#PBS -o %s/run_%s.pbs.out
# Stage 6 DIAGNOSTIC 58: the force each survivor of D56 adds, and the first-order maps of the force and of the
# hold on the last electron (AMENDMENT 13, A13.4 items 2 to 4).
# Job %s of 8. It makes %d runs in the order below (name, cap in s)%s; a run that fails does not stop the ones after it.
# No -q: the allocation is automatic and walltime 24:00:00 routes to day. The job script is FORMS TEST 1's committed
# run_F1.pbs with its names, runs and caps changed. No run is made twice:
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
""" % (j.lower(), HPC, j, j, len(mine), ", each followed by orca_vpot at the 2626 sites" if pot else "", HPC, j, j, j)
    if pot:
        t += """vpot() {
  d="$W/$1"; pts="$W/points.xyz"
  if [ -e "$d/vpot.log" ]; then echo "$1 vpot has a log already: not run again"; return; fi
  if ! grep -q "ORCA TERMINATED NORMALLY" "$d/job.out" 2>/dev/null; then echo "$1 vpot skipped, run not normal"; return; fi
  ( cd "$d" && "$ORCA/orca_vpot" job.gbw job.scfp "$pts" vpot.out job > vpot.log 2>&1 )
  echo "$1 vpot_exit=$? lines=$(wc -l < "$d/vpot.out" 2>/dev/null || echo 0) expected=2627"
}
"""
    for r, c in mine:
        t += "run %s %d\n" % (r, c)
        if pot:
            t += "vpot %s\n" % r
    t += 'echo "end=$(date)"\n'
    t += "n=0; for r in %s; do grep -q 'ORCA TERMINATED NORMALLY' \"$W/$r/job.out\" 2>/dev/null && n=$((n+1)); done\n" % " ".join(r for r, _ in mine)
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
        os.makedirs(os.path.dirname(os.path.join(dst, p)) or dst, exist_ok=True)
        with open(os.path.join(dst, p), "w") as f:
            f.write(F[p])
    print("wrote %d files into %s: 186 runs in 8 jobs" % (len(F), dst))


if __name__ == "__main__":
    main()
