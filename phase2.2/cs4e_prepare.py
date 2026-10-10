#!/usr/bin/env python3
"""
cs4e_prepare.py - Stage 4e inputs (method dependence of the sign beside the carboxylates; STAGE4E_CRITERIA.txt).
Runs on the PC from committed files only; Python 3.6, standard library. Reads, under RESULTS/phase2.2/:
  cs_stage4d/bare_{R,TS}/substrate.xyz      frame 41786's reactant and TS in their own coordinates (Stage 4d)
  cs_stage4d/pc/ALONE_<res>.pc.xz           the six residues' point charges (Stage 4d, unchanged)
  cs_stage4d/ALONE_R90_R/job.inp            Stage 4d's input, so that L1's header is checked to be identical
  grid_v2.tsv, aligned/frame_41786/transform.txt, A_v2.tsv   probe sites, taken back to frame 41786's coordinates
Writes RESULTS/phase2.2/cs_stage4e (must not exist): <run>/{job.inp,substrate.xyz,env.pc} for 170 runs, runs.tsv,
probes.tsv, batch_<level>.pbs (5). Deterministic: a second run into another folder gives identical files.
USAGE  python3 cs4e_prepare.py RESULTS_REPO [OUT]
"""
import lzma, math, os, sys
RES = sys.argv[1]
P2 = os.path.join(RES, "phase2.2")
D4 = os.path.join(P2, "cs_stage4d")
OUT = sys.argv[2] if len(sys.argv) > 2 else os.path.join(P2, "cs_stage4e")
if os.path.exists(OUT): sys.exit("STOP: %s exists" % OUT)
N = "C1 H1 H2 C2 C3 O1 O2 O3 C4 H3 C5 H4 C6 C7 H5 C8 H6 C9 H7 O4 H8 C10 O5 O6".split()
COMMON = "%maxcore 3000\n%pal nprocs 8 end\n%scf MaxIter 300 end\n%output\n  Print[P_Hirshfeld] 1\nend\n"
LEVELS = [("L1", "! B3LYP D3BJ def2-SVP def2/J RIJCOSX TightSCF\n"),
          ("L2", "! B3LYP 6-31G(d) TightSCF\n"),
          ("L3", "! HF 6-31G(d) TightSCF\n"),
          ("L4", "! MP2 6-31G(d) TightSCF\n"),
          ("L5", "! WB97X-D3 def2-SVP def2/J RIJCOSX TightSCF\n")]
RESIDUES = ["R90", "R7", "E78", "R116", "R63p", "K60p"]
SITES = [37, 52, 87, 209, 211]
def stop(s): sys.exit("STOP: " + s)
# L1 must be Stage 4d's input exactly
ref = open(os.path.join(D4, "ALONE_R90_R", "job.inp")).read()
mine = LEVELS[0][1] + COMMON + '%pointcharges "env.pc"\n* xyzfile -2 1 substrate.xyz\n'
if ref != mine: stop("L1's input differs from Stage 4d's")
SUBTXT = {g: open(os.path.join(D4, "bare_" + g, "substrate.xyz")).read() for g in ("R", "TS")}
SUB = {g: [tuple(float(v) for v in l.split()[1:4]) for l in SUBTXT[g].splitlines()[2:26]] for g in SUBTXT}
PC = {}
for r in RESIDUES:
    PC["ALONE_" + r] = lzma.open(os.path.join(D4, "pc", "ALONE_%s.pc.xz" % r), "rt").read()
# probe sites: aligned = original @ R + t, so original = (aligned - t) @ R.T
T = [[float(v) for v in l.split()] for l in open(os.path.join(P2, "aligned/frame_41786/transform.txt")) if not l.startswith("#")]
Rm, t = T[:3], T[3]
G = [l.rstrip("\n").split("\t") for l in open(os.path.join(P2, "grid_v2.tsv")) if not l.startswith("#") and not l.startswith("idx")]
fr = [l for l in open(os.path.join(P2, "A_v2.tsv")) if l.startswith("# frames")][0].split()[3:]
A = [l.rstrip("\n").split("\t") for l in open(os.path.join(P2, "A_v2.tsv")) if not l.startswith("#")]
if len(G) != 233 or len(A) != 233: stop("grid_v2 / A_v2 size")
prow = []
for s in SITES:
    if int(G[s][0]) != s: stop("grid row %d" % s)
    ga = [float(v) - t[k] for k, v in enumerate(G[s][1:4])]
    o = tuple(sum(ga[j] * Rm[i][j] for j in range(3)) for i in range(3))
    d = sorted((math.sqrt(sum((o[k] - x[k]) ** 2 for k in range(3))), N[j]) for j, x in enumerate(SUB["R"]))
    prow.append("%d\t%s\t%.6f %.6f %.6f\t%s %.3f, %s %.3f\t%s" % (s, " ".join(G[s][1:4]), o[0], o[1], o[2],
                d[0][1], d[0][0], d[1][1], d[1][0], A[s][fr.index("41786")]))
    for sign, tag in ((1.0, "p"), (-1.0, "m")):
        PC["P%d%s" % (s, tag)] = "1\n%10.6f %14.8f %14.8f %14.8f\n" % (sign, o[0], o[1], o[2])
ENVS = ["bare"] + ["ALONE_" + r for r in RESIDUES] + ["P%d%s" % (s, tag) for s in SITES for tag in "pm"]
os.makedirs(OUT)
rows = []
for lev, head in LEVELS:
    for e in ENVS:
        for g in ("R", "TS"):
            run = "%s_%s_%s" % (lev, e, g); d = os.path.join(OUT, run); os.makedirs(d)
            open(os.path.join(d, "substrate.xyz"), "w").write(SUBTXT[g])
            inp = head + COMMON
            if e != "bare":
                open(os.path.join(d, "env.pc"), "w").write(PC[e]); inp += '%pointcharges "env.pc"\n'
            open(os.path.join(d, "job.inp"), "w").write(inp + "* xyzfile -2 1 substrate.xyz\n")
            npc = 0 if e == "bare" else int(PC[e].split("\n")[0])
            rows.append("%s\t%s\t%s\t%s\t%d" % (run, lev, e, g, npc))
open(os.path.join(OUT, "runs.tsv"), "w").write("# run\tlevel\tenvironment\tgeometry\tn_point_charges\n" + "\n".join(rows) + "\n")
open(os.path.join(OUT, "probes.tsv"), "w").write(
    "# site\taligned xyz (grid_v2)\tframe 41786 xyz\tnearest substrate atoms at R (A)\ta_41786 (A_v2, kcal/mol per e)\n"
    + "\n".join(prow) + "\n")
PBS = """#!/bin/bash
#PBS -N cm34e_%s
#PBS -l select=1:ncpus=8:mem=24gb
#PBS -l walltime=168:00:00
#PBS -m ae
#PBS -M 18660916@sun.ac.za
#PBS -j oe
#PBS -o %s/batch_%s.pbs.out
export PATH="/apps/openmpi/4.1.1/bin:$PATH"
export LD_LIBRARY_PATH="/home/apps2/ORCA/6.0.1/lib:/apps/openmpi/4.1.1/lib:/apps/mambaforge/envs/medaka/lib:${LD_LIBRARY_PATH:-}"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
echo "host=$(hostname) start=$(date)"
"""
SDH = "/home/18660916/system_development/phase2.2/cs_stage4e"
for lev, _ in LEVELS:
    s = PBS % (lev, SDH, lev)
    for r in rows:
        f = r.split("\t")
        if f[1] != lev: continue
        s += "cd %s/%s\n/home/apps2/ORCA/6.0.1/orca job.inp > job.out 2>&1 </dev/null\n" % (SDH, f[0])
        s += 'echo "%s orca_exit=$? $(grep -c \'ORCA TERMINATED NORMALLY\' job.out) $(date)"\n' % f[0]
        s += "rm -f job.gbw job*.tmp 2>/dev/null\n"
    s += 'echo "end=$(date)"\n'
    open(os.path.join(OUT, "batch_%s.pbs" % lev), "w").write(s)
print("cs_stage4e: %d runs at %d levels, probes.tsv, runs.tsv, %d batch jobs" % (len(rows), len(LEVELS), len(LEVELS)))
