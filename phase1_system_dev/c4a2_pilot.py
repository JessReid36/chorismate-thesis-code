#!/usr/bin/env python3
"""
c4a2_pilot.py - inputs of the C4 Amendment 2 pilot (exploratory; pilot_a2/PILOT_A2.txt): the TS-anchored forth and back on
one finished frame. RUN ON hpc1 (login node; builds inputs only, runs no QM). Python 3.6, standard library; deterministic.

  T  FRAME STAGING OUT   T1: the environment relaxed around the climbing-image substrate. Start: the final band's climbing
                         image (neb_MEP.allxyz, full system, rows copied verbatim). Active atoms: C4's 1959 minus the
                         substrate (QM atoms 6207-6230): 1935 atoms; the substrate does not move (ORCA's L-Opt moves only
                         active atoms and ignores constraints and fragments). Input otherwise S2_FRAME/reopt.inp.
  R  FRAME STAGING OUT   R2: back to the reactant, after T1 has finished (gradients met, energy plateau, c4_extract ok).
                         Start: T1's last geometry (reopt.xyz) with the substrate rows replaced by image k_low's
                         (neb_MEP.allxyz). Active atoms: all 1959. Input file: S2_FRAME/reopt.inp, byte for byte.
STAGING/expected.sha256 lists "<sha256>  <absolute hpc1 path>" for every source file (built on the PC from the results and
code repos); each source is checked before use. OUT is where the inputs are written (PD for the real inputs, a scratch
folder for the second, identical run); the job scripts always point at PD. T1's outputs are always read from PD.
Writes OUT/T/{active_atoms.txt, T1_FRAME/{start.xyz,reopt.inp,reopt.pbs,prms link}, GENERATOR_LOG_T1_FRAME.txt} or the
same under OUT/R for R2. Stops (exit 1) on any failed check; nothing is written unless every check passes. OUT "-" runs
every check and writes nothing.
"""
import hashlib, math, os, re, sys
ROOT = "/home/18660916/system_development"
SD = ROOT + "/05_qmmm/23_c4_reactant"
PD = SD + "/pilot_a2"
QM0, NQM, NATOMS = 6207, 24, 55680
H = 627.5094740631
PLATEAU_EH, PLATEAU_N = 2.5e-4, 50
CLASH_STOP = 1.0          # A: a new substrate atom closer than this to any other atom stops R
if len(sys.argv) != 5 or sys.argv[1] not in ("T", "R"): sys.exit(__doc__)
MODE, FRAME, STAGING, OUT = sys.argv[1:]
# canonical hpc1 paths (job scripts, checksum list) are read through loc(): the home folder is taken from STAGING, so the
# same code runs in a copy of the tree elsewhere (the sandbox); on hpc1 loc() changes nothing
if "/system_development/" not in STAGING: sys.exit("STAGING must lie under .../system_development/")
HOME_REAL = os.path.abspath(STAGING)[:os.path.abspath(STAGING).index("/system_development/")]
def loc(p): return HOME_REAL + p[len("/home/18660916"):] if p.startswith("/home/18660916/") else p
LOG = []
def say(s): LOG.append(s); print(s)
def stop(s): print("STOP: " + s); sys.exit(1)
def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""): h.update(b)
    return h.hexdigest()
EXP = {}
for l in open(os.path.join(STAGING, "expected.sha256")):
    if l.strip(): EXP[l[66:].rstrip("\n")] = l[:64]
def src(p):
    if p not in EXP: stop("no expected checksum for " + p)
    if not os.path.isfile(loc(p)): stop("missing " + p)
    if sha(loc(p)) != EXP[p]: stop("checksum differs from the repos' record: " + p)
    return loc(p)
# ---- sources common to T and R
BL = [l.rstrip("\n").split("\t") for l in open(src(SD + "/c4_bands.tsv")) if l.strip() and not l.startswith("#")]
BAND = {r[0]: dict(zip(BL[0], r)) for r in BL[1:]}
if FRAME not in BAND: stop("frame %s is not in c4_bands.tsv" % FRAME)
B = BAND[FRAME]; CI, KL = int(B["ci"]), int(B["k_low"])
if not 0 < KL < CI: stop("frame %s: k_low %d, climbing image %d (the pilot needs 0 < k_low < ci)" % (FRAME, KL, CI))
ACT = [int(x) for x in open(src(SD + "/active_atoms.txt")).read().split()]
QMS = list(range(QM0, QM0 + NQM))
if len(ACT) != 1959 or not set(QMS) <= set(ACT): stop("active_atoms.txt")
ENV = [i for i in ACT if not QM0 <= i < QM0 + NQM]
S2 = "%s/S2_%s" % (SD, FRAME)
TEMPL = open(src(S2 + "/reopt.inp")).read()
TL = TEMPL.splitlines()
want = ["! QMMM B3LYP D3BJ def2-SVP def2/J RIJCOSX L-Opt TightSCF", "%maxcore 3000", "%pal nprocs 8 end",
        "%scf MaxIter 200 end", "%geom", "  MaxIter 2000", "  TolMaxG 1.0e-4", "  TolRMSG 3.0e-5", "end", "%qmmm",
        "  QMAtoms {6207:6230} end"]
if TL[:11] != want or TL[12:] != ['  ORCAFFFilename "complex_solvated.ORCAFF.prms"', "end", "*xyzfile -2 1 start.xyz"] \
        or TL[11] != "  ActiveAtoms {" + " ".join(map(str, ACT)) + "} end":
    stop("S2_%s/reopt.inp is not the C4 step-2 template with the 1959 active atoms" % FRAME)
src(SD + "/c4_extract.py")
PRMS = src(SD + "/complex_solvated.ORCAFF.prms")
def prms_elements(p):
    L = open(p).read().splitlines(); i = [k for k, l in enumerate(L) if l.strip() == "$atoms"][0]
    n = int(L[i + 1].split()[0])
    return [l.split()[1] for l in L[i + 2:i + 2 + n]]
EL = prms_elements(PRMS)
if len(EL) != NATOMS: stop("force field has %d atoms" % len(EL))
def read_allxyz(p):
    L = open(p).read().splitlines(); S = []; i = 0
    while i < len(L):
        s = L[i].strip()
        if not s or s.startswith(">"): i += 1; continue
        n = int(s.split()[0])
        if i + 2 + n > len(L): stop("truncated structure in " + p)
        S.append((L[i + 1], i + 2, n)); i += n + 2
    return L, S
def rows_of(L, s):
    R = [l.split() for l in L[s[1]:s[1] + s[2]]]
    if any(len(r) < 4 for r in R): stop("a short atom line")
    return R
def read_xyz(p):
    L = open(p).read().splitlines(); n = int(L[0].split()[0]); R = [l.split() for l in L[2:2 + n]]
    if n != NATOMS or len(R) != NATOMS or any(len(r) < 4 for r in R): stop("%s holds %d atoms, expected %d" % (p, n, NATOMS))
    return R
def fl(r): return (float(r[1]), float(r[2]), float(r[3]))
def dmax(a, b): return max(abs(a[k] - b[k]) for k in range(3))
ALLXYZ = "%s/05_qmmm/%s/frame_%s/neb_MEP.allxyz" % (ROOT, B["band"], FRAME)
XL, XS = read_allxyz(src(ALLXYZ))
if len(XS) != 10 or any(s[2] != NATOMS for s in XS): stop("neb_MEP.allxyz is not 10 structures of %d atoms" % NATOMS)
X0, XK, XC = rows_of(XL, XS[0]), rows_of(XL, XS[KL]), rows_of(XL, XS[CI])
del XL, XS
for nm, X in (("0", X0), ("k_low", XK), ("ci", XC)):
    if any(X[i][0].lower() != EL[i].lower() for i in range(NATOMS)): stop("image %s: elements differ from the force field" % nm)
AS = set(ACT)
fix = max(max(dmax(fl(X[i]), fl(X0[i])) for i in range(NATOMS) if i not in AS) for X in (XK, XC))
if fix > 1e-6: stop("an atom outside the active region differs between band images (%.1e A)" % fix)
say("c4a2_pilot.py %s, frame %s: band %s, climbing image %d, k_low %d; E_ci %s, E_klow %s Eh (final band)" % (
    MODE, FRAME, B["band"], CI, KL, B["E_ci"], B["E_klow"]))
say("sources: each matches its sha256 from the repos (%d listed); S2_%s/reopt.inp is the step-2 template; band images 0,"
    " k_low and ci: elements = force field, atoms outside the active region identical" % (len(EXP), FRAME))
PBS = """#!/bin/bash
#PBS -N c4p_{run}
#PBS -l select=1:ncpus=8:mem=120gb
#PBS -l walltime=176:00:00
#PBS -m ae
#PBS -M 18660916@sun.ac.za
#PBS -j oe
#PBS -o {pd}/{leg}/{run}/reopt.pbs.out
cd {pd}/{leg}/{run}
export PATH="/apps/openmpi/4.1.1/bin:$PATH"
export LD_LIBRARY_PATH="/home/apps2/ORCA/6.0.1/lib:/apps/openmpi/4.1.1/lib:/apps/mambaforge/envs/medaka/lib:${{LD_LIBRARY_PATH:-}}"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
for f in reopt.inp start.xyz complex_solvated.ORCAFF.prms; do [[ -s "$f" ]] || {{ echo "FAIL missing $f"; exit 1; }}; done
[[ -e reopt.out ]] && {{ echo "FAIL reopt.out exists: this run was already started"; exit 1; }}
[[ $(df -Pk . | awk 'NR==2{{print $4}}') -ge 3000000 ]] || {{ echo "FAIL less than 3 GB free here (the trajectory needs up to ~1.3 GB)"; exit 1; }}
echo "host=$(hostname) start=$(date)"
/home/apps2/ORCA/6.0.1/orca reopt.inp > reopt.out 2>&1 </dev/null
echo "orca_exit=$? $(date)"
grep -q "ORCA TERMINATED NORMALLY" reopt.out && grep -q "The minimization has converged" reopt.out && echo C4_REOPT_CONVERGED || echo C4_REOPT_INCOMPLETE
env -u LD_LIBRARY_PATH python3 {sd}/c4_extract.py {pd}/{leg}/{run} && echo "extract=ok" || echo "extract=FAILED"
# c4_extract.py checks that every atom outside ../active_atoms.txt (for T1 this includes the substrate) is where start.xyz
# put it. A finished run's trajectory is not needed: its last geometry is reopt.xyz. A run killed by the walltime never
# reaches this line, so its trajectory stays for a continuation
rm -f reopt_trj.dcd reopt*.tmp 2>/dev/null
echo "end=$(date)"
"""
def write(leg, run, active, start_text, inp_text, extra):
    if OUT == "-":
        for s in extra: say(s)
        say("check only: every check passed, nothing written"); return
    od = os.path.join(OUT, leg); rd = os.path.join(od, run)
    if os.path.exists(rd): stop(rd + " exists")
    os.makedirs(rd)
    files = {os.path.join(od, "active_atoms.txt"): "\n".join(str(i) for i in active) + "\n",
             os.path.join(rd, "start.xyz"): start_text, os.path.join(rd, "reopt.inp"): inp_text,
             os.path.join(rd, "reopt.pbs"): PBS.format(run=run, leg=leg, pd=PD, sd=SD)}
    ap = os.path.join(od, "active_atoms.txt")
    if os.path.exists(ap) and open(ap).read() != files[ap]: stop(ap + " exists and differs")
    for p, t in sorted(files.items()):
        with open(p, "w") as fh: fh.write(t)
    os.symlink("../../../complex_solvated.ORCAFF.prms", os.path.join(rd, "complex_solvated.ORCAFF.prms"))
    if os.path.abspath(OUT) == os.path.abspath(loc(PD)) and sha(os.path.join(rd, "complex_solvated.ORCAFF.prms")) != EXP[SD + "/complex_solvated.ORCAFF.prms"]:
        stop("force-field link")
    for s in extra: say(s)
    say("written (sha256, path relative to %s):" % leg)
    for p in sorted(files): say("  %s  %s" % (sha(p), os.path.relpath(p, od)))
    say("  link %s/complex_solvated.ORCAFF.prms -> ../../../complex_solvated.ORCAFF.prms" % run)
    open(os.path.join(od, "GENERATOR_LOG_%s.txt" % run), "w").write("\n".join(LOG) + "\n")
def mind(sub_rows, rows, sub_idx):
    """smallest distance from the substrate atoms to any other atom; (distance, substrate atom, other atom)"""
    S = [(i, fl(sub_rows[i - QM0])) for i in sub_idx]; best = (1e9, -1, -1)
    for j in range(NATOMS):
        if QM0 <= j < QM0 + NQM: continue
        b = fl(rows[j])
        for i, a in S:
            d = (a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2
            if d < best[0]: best = (d, i, j)
    return (best[0] ** 0.5, best[1], best[2])
if MODE == "T":
    run = "T1_" + FRAME
    inp = "\n".join(TL[:11] + ["  ActiveAtoms {" + " ".join(map(str, ENV)) + "} end"] + TL[12:]) + "\n"
    start = ("%d\nframe %s, image %d (climbing image) of the final band (05_qmmm/%s/frame_%s/neb_MEP.allxyz): start of T1\n"
             % (NATOMS, FRAME, CI, B["band"], FRAME) + "".join("%s %s %s %s\n" % tuple(r[:4]) for r in XC))
    write("T", run, ENV, start, inp, [
        "T1 input: S2_%s/reopt.inp with ActiveAtoms = the 1959 C4 active atoms minus 6207-6230 (%d atoms); every other"
        " line identical" % (FRAME, len(ENV)),
        "T1 start: band image %d (climbing image), 55680 atoms, coordinates copied verbatim" % CI])
else:
    run = "R2_" + FRAME; T1 = os.path.join(loc(PD), "T", "T1_" + FRAME)
    po = os.path.join(T1, "reopt.pbs.out")
    if not os.path.isfile(po): stop("T1 has not run (no %s)" % po)
    pt = open(po).read()
    if "C4_REOPT_CONVERGED" not in pt or "extract=ok" not in pt: stop("T1 did not finish with converged gradients and extract=ok; see " + po)
    cp = os.path.join(T1, "reopt-minimize-ener.csv")
    E = [float(l.split(";")[5]) for l in open(cp) if l.strip() and not l.lstrip().startswith("#")]
    slide = E[-1 - PLATEAU_N] - E[-1] if len(E) > PLATEAU_N else 0.0
    if slide > PLATEAU_EH: stop("T1's energy fell by %.2e Eh over its last %d steps (more than %.1e): continue it before R2" % (slide, PLATEAU_N, PLATEAU_EH))
    TS = read_xyz(os.path.join(T1, "start.xyz")); TF = read_xyz(os.path.join(T1, "reopt.xyz"))
    if any(TS[i][:4] != XC[i][:4] for i in range(NATOMS)): stop("T1's start.xyz is not the band's climbing image")
    ES = set(ENV)
    moved = max(dmax(fl(TF[i]), fl(TS[i])) for i in range(NATOMS) if i not in ES)
    if moved > 2e-6: stop("in T1 an atom outside its active atoms (the substrate included) moved by %.1e A" % moved)
    if any(TF[i][0].lower() != EL[i].lower() for i in range(NATOMS)): stop("T1's reopt.xyz: elements")
    new = [XK[i] if QM0 <= i < QM0 + NQM else TF[i] for i in range(NATOMS)]
    c_new = mind([XK[i] for i in QMS], TF, QMS)
    c_band = mind([XK[i] for i in QMS], XK, QMS)
    c_t1 = mind([TF[i] for i in QMS], TF, QMS)
    extra = ["T1: %d optimiser steps, E_start %.8f, E_min %.8f Eh; energy fall over its last %d steps %.2e Eh;"
             " atoms outside its active atoms (the substrate included) unmoved (max %.1e A)" % (len(E), E[0], E[-1], PLATEAU_N, slide, moved),
             "R2 start: T1's last geometry with atoms 6207-6230 (the substrate) taken from band image %d (k_low)" % KL,
             "closest contact, substrate to any other atom (A): R2 start %.3f (atoms %d-%d); band image k_low %.3f;"
             " T1's last geometry (climbing-image substrate) %.3f" % (c_new + (c_band[0], c_t1[0]))]
    if c_new[0] < CLASH_STOP:
        for s in extra: say(s)
        stop("a substrate atom of the R2 start is %.3f A from atom %d (below %.1f A): nothing written" % (c_new[0], c_new[2], CLASH_STOP))
    start = ("%d\nframe %s: T1's last geometry (pilot_a2/T/T1_%s/reopt.xyz), substrate 6207-6230 from band image %d"
             " (k_low): start of R2\n" % (NATOMS, FRAME, FRAME, KL) + "".join("%s %s %s %s\n" % tuple(r[:4]) for r in new))
    write("R", run, ACT, start, TEMPL, extra + ["R2 input: S2_%s/reopt.inp byte for byte (sha256 %s)" % (FRAME, EXP[S2 + "/reopt.inp"])])
