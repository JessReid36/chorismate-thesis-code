#!/usr/bin/env python3
"""
spa_prepare.py - inputs of the stationary-point barriers under SP Amendment 1, the alternating scheme (results
05_qmmm/24_sp_barriers/SP_AMENDMENT1.txt; everything not changed there is SP_PLAN.txt). RUN ON hpc1 (login node; builds
inputs only, runs no QM). Python 3.6, standard library; deterministic (no dates, no random). It is sp_prepare.py with
the four optimisations (ts, r_irc, p_irc, r0) replaced by alternations; regions and start files are made exactly as
sp_prepare.py makes them.

USAGE  python3 spa_prepare.py EXPECTED OUT RUN [RUN ...]
  EXPECTED  lines "<sha256>  <absolute hpc1 path>" for every source file (built on the PC from the repos' records);
            every source is checked against it before use
  OUT       where the run folders are written: SD for the real inputs, any other folder for the second, identical
            generation (the job scripts always point at SD); "-" runs every check and writes nothing
  RUN       FRAME:RADIUS:BULK, e.g. 42436:8:md. RADIUS 6 or 8 (A). BULK md = every atom outside the region stays at the
            frame's MD structure (frame_F_CHA2.pdb); ci = it stays at the final band's climbing image
Writes, per RUN, OUT/f<FRAME>_r<RADIUS>[_ci]_a1/ (must not exist):
  region.txt                    region atoms (0-based ORCA indices), one per line; the environment is the region
                                without the 24 substrate atoms 6207-6230 (the core)
  starts.sha256                 sha256 of the two start files the generator writes (ts, r0)
  ts/start.xyz                  TS start: region atoms from the band's climbing image, the rest from the bulk
  r0/start.xyz                  R0 start: region atoms from the band's image 0 (the Phase 1 reactant), the rest from the bulk
  ts/ r_irc/ p_irc/ r0/         env.inp (L-Opt of the environment, core fixed) and core.inp (the core alone, environment
                                fixed: OptTS with an exact substrate Hessian for ts, Opt for the minima); the job runs
                                them in turn, cycle K in STEP/cKK_env and STEP/cKK_core (spa_stage.py alt_prep)
  ts_freq/ r_irc_freq/ p_irc_freq/ r0_freq/ freq.inp   frequencies, substrate atoms only displaced (partial Hessian)
  irc/irc.inp                   IRC from the TS, substrate atoms moving, Hessian from ts_freq
  job_ts.pbs job_ts_freq_irc.pbs job_r_irc.pbs job_r_irc_freq.pbs job_p_irc.pbs job_p_irc_freq.pbs job_r0.pbs
  job_r0_freq.pbs, submit.sh (qsub with dependencies), GENERATOR_LOG.txt
  every step folder links complex_solvated.ORCAFF.prms -> ../../complex_solvated.ORCAFF.prms (SD's copy)
Start files of the later steps are written in the jobs by spa_stage.py, after each step's checks.
Stops (exit 1) on any failed check; nothing is written unless every check of every RUN passes.
"""
import hashlib, os, re, sys

ROOT = "/home/18660916/system_development"
SD = ROOT + "/05_qmmm/24_sp_barriers"
QM0, NQM, NATOMS = 6207, 24, 55680
QMS = list(range(QM0, QM0 + NQM))
BRK, FRM = (6214, 6215), (6207, 6219)      # breaking O-C and forming C-C (0-based; step16_restraint_scan_driver.sh)
PRMTOP = ROOT + "/03_amber/tleap_build/complex_solvated.prmtop"
PRMS = ROOT + "/05_qmmm/19_ensemble/frame_41786/complex_solvated.ORCAFF.prms"
BANDS = ROOT + "/05_qmmm/23_c4_reactant/c4_bands.tsv"
R12 = ROOT + "/05_qmmm/23_c4_reactant/active_atoms.txt"
CLASH = 1.0          # A: a non-bonded contact between a region atom and an outside atom closer than this stops the RUN
TOL_BAND = 6e-4      # A: band rows (allxyz) against the 3-decimal MD PDB, for atoms the band held fixed

if len(sys.argv) < 4: sys.exit(__doc__)
EXPECTED, OUT, RUNS = sys.argv[1], sys.argv[2], sys.argv[3:]
LOG = []
def say(s): LOG.append(s); print(s)
def stop(s): print("STOP: " + s); sys.exit(1)
def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""): h.update(b)
    return h.hexdigest()
def shas(t): return hashlib.sha256(t.encode()).hexdigest()
EXP = {}
for l in open(EXPECTED):
    if l.strip():
        if len(l) < 67 or l[64:66] != "  ": stop("malformed line in " + EXPECTED + ": " + l.rstrip())
        EXP[l[66:].rstrip("\n")] = l[:64]
def src(p):
    if p not in EXP: stop("no expected checksum for " + p)
    if not os.path.isfile(p): stop("missing " + p)
    if sha(p) != EXP[p]: stop("checksum differs from the repos' record: " + p)
    return p

# ---- runs requested
SPEC = []
for r in RUNS:
    m = re.match(r"^(\d{5}):(6|8):(md|ci)$", r)
    if not m: stop("RUN must be FRAME:RADIUS:BULK with RADIUS 6 or 8 and BULK md or ci: " + r)
    f, rad, bulk = m.group(1), int(m.group(2)), m.group(3)
    name = "f%s_r%d%s_a1" % (f, rad, "_ci" if bulk == "ci" else "")
    if any(s[3] == name for s in SPEC): stop("RUN given twice: " + r)
    SPEC.append((f, rad, bulk, name))

# ---- common sources
def read_prmtop(p):
    F = {}; cur = None; w = None
    for l in open(p):
        if l.startswith("%FLAG"): cur = l.split()[1]; F[cur] = []; continue
        if l.startswith("%FORMAT"):
            m = re.match(r"\((\d+)([aIE])(\d+)", l[7:].strip())
            if not m: stop("prmtop format line: " + l.strip())
            w = int(m.group(3)); continue
        if l.startswith("%") or cur is None: continue
        s = l.rstrip("\n")
        F[cur].extend(s[i:i + w] for i in range(0, len(s), w))
    return F
P = read_prmtop(src(PRMTOP))
PTR = [int(x) for x in P["POINTERS"] if x.strip()]
if PTR[0] != NATOMS: stop("prmtop holds %d atoms" % PTR[0])
NRES = PTR[11]
ANAME = [x.strip() for x in P["ATOM_NAME"]][:NATOMS]
RLAB = [x.strip() for x in P["RESIDUE_LABEL"]][:NRES]
RPTR = [int(x) - 1 for x in P["RESIDUE_POINTER"] if x.strip()]
if len(ANAME) != NATOMS or len(RLAB) != NRES or len(RPTR) != NRES or RPTR[0] != 0 or RPTR != sorted(RPTR):
    stop("prmtop residue tables")
RPTR.append(NATOMS)
RESOF = [0] * NATOMS
for r in range(NRES):
    for i in range(RPTR[r], RPTR[r + 1]): RESOF[i] = r
BONDS = []
for key in ("BONDS_INC_HYDROGEN", "BONDS_WITHOUT_HYDROGEN"):
    v = [int(x) for x in P[key] if x.strip()]
    if len(v) % 3: stop("prmtop " + key)
    BONDS.extend((v[k] // 3, v[k + 1] // 3) for k in range(0, len(v), 3))
del P
if len(set(RESOF[i] for i in QMS)) != 1 or RLAB[RESOF[QM0]] != "CHA" or RPTR[RESOF[QM0]] != QM0 \
        or RPTR[RESOF[QM0] + 1] != QM0 + NQM:
    stop("atoms 6207-6230 are not exactly one CHA residue in the prmtop")
def prms_elements(p):
    L = open(p).read().splitlines(); i = [k for k, l in enumerate(L) if l.strip() == "$atoms"][0]
    n = int(L[i + 1].split()[0])
    return [l.split()[1] for l in L[i + 2:i + 2 + n]]
EL = prms_elements(src(PRMS))
if len(EL) != NATOMS: stop("force field has %d atoms" % len(EL))
SDPRMS = SD + "/complex_solvated.ORCAFF.prms"
if OUT != "-" and os.path.abspath(OUT) == SD and (not os.path.isfile(SDPRMS) or sha(SDPRMS) != EXP[PRMS]):
    stop(SDPRMS + " is missing or is not the force field " + PRMS)
BL = [l.rstrip("\n").split("\t") for l in open(src(BANDS)) if l.strip() and not l.startswith("#")]
BAND = {r[0]: dict(zip(BL[0], r)) for r in BL[1:]}
if len(BAND) != 30: stop("c4_bands.tsv holds %d frames" % len(BAND))
ACT12 = [int(x) for x in open(src(R12)).read().split()]
A12 = set(ACT12)
if len(ACT12) != 1959 or len(A12) != 1959 or not set(QMS) <= A12: stop("the 1959-atom active list (R12)")
say("spa_prepare.py: %d run(s) in this call: %s" % (len(SPEC), " ".join(s[3] for s in SPEC)))
say("prmtop: %d atoms, %d residues, %d bonds; atoms 6207-6230 are one CHA residue; force field: %d atoms" % (
    NATOMS, NRES, len(BONDS), len(EL)))

# ---- geometry helpers
def read_pdb(p):
    R = []
    for l in open(p):
        if l[:6] in ("ATOM  ", "HETATM"):
            R.append((l[12:16].strip(), l[17:20].strip(), l[30:38].strip(), l[38:46].strip(), l[46:54].strip()))
    return R
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
    R = [l.split()[:4] for l in L[s[1]:s[1] + s[2]]]
    if any(len(r) < 4 for r in R): stop("a short atom line")
    return R
def read_xyz_frames(p):
    L = open(p).read().splitlines(); S = []; i = 0
    while i < len(L):
        if not L[i].strip(): i += 1; continue
        n = int(L[i].split()[0]); S.append([[float(v) for v in l.split()[1:4]] for l in L[i + 2:i + 2 + n]]); i += n + 2
    return S
def fl(r): return (float(r[1]), float(r[2]), float(r[3]))
def dist(a, b): return ((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2) ** 0.5
def dmax(a, b): return max(abs(a[k] - b[k]) for k in range(3))
def braces(t, key):
    m = re.search(key + r"\s*\{([^}]*)\}", t)
    if not m: return None
    out = []
    for tok in m.group(1).split():
        if ":" in tok:
            a, b = tok.split(":"); out.extend(range(int(a), int(b) + 1))
        else: out.append(int(tok))
    return out
def neighbours(X, A, B, cut):
    """pairs (i in A, j in B, d) with d < cut, by a cell grid"""
    g = {}
    for j in B:
        x = X[j]; g.setdefault((int(x[0] // cut), int(x[1] // cut), int(x[2] // cut)), []).append(j)
    out = []
    for i in A:
        x = X[i]; c = (int(x[0] // cut), int(x[1] // cut), int(x[2] // cut))
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    for j in g.get((c[0] + dx, c[1] + dy, c[2] + dz), ()):
                        d = dist(x, X[j])
                        if d < cut: out.append((i, j, d))
    return out
BONDED = set()
for i, j in BONDS: BONDED.add((i, j)); BONDED.add((j, i))

# ---- input templates (ORCA 6.0.1). core: sp_prepare.py's ts.inp and opt.inp (OptTS with Calc_Hess and Hybrid_Hess over
# the substrate, as in the frame-820 input step18c_ts_reduced.inp; Opt) with the substrate alone active; env: the C4
# pilot's T1 input (L-Opt, TolMaxG 1e-4, TolRMSG 3e-5) with the region minus the substrate active and MaxIter 5000 in
# place of 2000 (SP_AMENDMENT1.txt, A1.3); freq and irc: unchanged from sp_prepare.py
HEAD = "! QMMM B3LYP D3BJ def2-SVP def2/J RIJCOSX "
def qmmm(active):
    return ("%qmmm\n  QMAtoms {6207:6230} end\n  ActiveAtoms {" + " ".join(map(str, active)) + "} end\n"
            '  ORCAFFFilename "complex_solvated.ORCAFF.prms"\nend\n*xyzfile -2 1 start.xyz\n')
def inp_env(env):
    return (HEAD + "L-Opt TightSCF\n%maxcore 6000\n%pal nprocs 16 end\n%scf MaxIter 250 end\n"
            "%geom\n  MaxIter 5000\n  TolMaxG 1.0e-4\n  TolRMSG 3.0e-5\nend\n" + qmmm(env))
def inp_core_ts():
    return (HEAD + "OptTS TightOpt TightSCF\n%maxcore 12000\n%pal nprocs 16 end\n%scf MaxIter 250 end\n"
            "%geom\n  Calc_Hess true\n  Hybrid_Hess {" + " ".join(map(str, QMS)) + "} end\n  MaxIter 500\nend\n"
            + qmmm(QMS))
def inp_core_min():
    return (HEAD + "Opt TightOpt TightSCF\n%maxcore 6000\n%pal nprocs 16 end\n%scf MaxIter 250 end\n"
            "%geom\n  MaxIter 500\nend\n" + qmmm(QMS))
def inp_freq():
    return HEAD + "NumFreq TightSCF\n%maxcore 8000\n%pal nprocs 16 end\n%scf MaxIter 250 end\n" + qmmm(QMS)
def inp_irc():
    return (HEAD + "IRC TightSCF\n%maxcore 8000\n%pal nprocs 16 end\n%scf MaxIter 250 end\n"
            '%irc\n  MaxIter 100\n  Direction both\n  InitHess read\n  Hess_Filename "freq.hess"\n  hessMode 0\nend\n'
            + qmmm(QMS))

# ---- job scripts. Each prints SP_<STEP>_PASS or SP_<STEP>_FAIL per step and ends with "end_job=" when it was not
# killed; a failed step exits 1, so the jobs that depend on it (afterok) never start.
PBS_HEAD = """#!/bin/bash
#PBS -N {name}
#PBS -l select=1:ncpus=16:mem={mem}gb
#PBS -l walltime=176:00:00
#PBS -m a
#PBS -M 18660916@sun.ac.za
#PBS -j oe
#PBS -o {rd}/{job}.pbs.out
set -u
RD={rd}
SD={sd}
ORCA=/home/apps2/ORCA/6.0.1/orca
export PATH="/apps/openmpi/4.1.1/bin:$PATH"
export LD_LIBRARY_PATH="/home/apps2/ORCA/6.0.1/lib:/apps/openmpi/4.1.1/lib:/apps/mambaforge/envs/medaka/lib:${{LD_LIBRARY_PATH:-}}"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
PY() {{ env -u LD_LIBRARY_PATH python3 "$SD/spa_stage.py" "$@"; }}
used() {{ qstat -f "${{PBS_JOBID:-none}}" 2>/dev/null | grep -E "resources_used\\.(mem|vmem|walltime|cput|ncpus)" | sed 's/^ */pbs_/'; }}
quit() {{ used; echo "end_job=$(date)"; exit "$1"; }}
run_step() {{   # folder, ORCA basename, GB of free disk needed
  cd "$RD/$1" || {{ echo "SP_$1_FAIL no folder $RD/$1"; quit 1; }}
  for f in "$2.inp" start.xyz complex_solvated.ORCAFF.prms; do [[ -s "$f" ]] || {{ echo "SP_$1_FAIL missing $1/$f"; quit 1; }}; done
  [[ -e "$2.out.xz" ]] && {{ echo "SP_$1_FAIL $1/$2.out.xz exists: this step was already run"; quit 1; }}
  [[ $(df -Pk . | awk 'NR==2{{print $4}}') -ge $(( $3 * 1000000 )) ]] || {{ echo "SP_$1_FAIL less than $3 GB free in $RD"; quit 1; }}
  ( set -o noclobber; : > "$2.out" ) 2>/dev/null || {{ echo "SP_$1_FAIL $1/$2.out exists: this step was already started"; quit 1; }}
  echo "step=$1 start=$(date)"
  "$ORCA" "$2.inp" > "$2.out" 2>&1 </dev/null
  echo "step=$1 orca_exit=$? end=$(date)"
  if grep -q -i -e "No space left on device" -e "Disk quota exceeded" "$2.out" || [[ $(df -Pk . | awk 'NR==2{{print $4}}') -lt 1000000 ]]; then
    echo "SP_$1_FAIL machinery stop: the disk filled during ORCA (message in $1/$2.out, or under 1 GB free after it); no check was run"; quit 1
  fi
  if grep -q -e "INPUT ERROR" -e "UNRECOGNIZED OR DUPLICATED KEYWORD" "$2.out" && ! grep -q "FINAL SINGLE POINT ENERGY" "$2.out" && [[ ! -s "$2-minimize-ener.csv" ]]; then
    echo "SP_$1_FAIL machinery stop: ORCA reported an input error before its first energy (see $1/$2.out); no check was run"; quit 1
  fi
}}
finish_step() {{   # folder, ORCA basename, 1 if the step passed its checks: sizes recorded (bytes, name); after a pass, full-system
                  # trajectories and scratch files are deleted (in a cycle folder also the full-system PDB, the densities
                  # and start.xyz, which is the previous part's recorded output); any output over 50 MB is compressed with xz
  cd "$RD/$1" || return
  find . -maxdepth 1 \\( -type f -o -type l \\) -printf "%s %f\\n" | LC_ALL=C sort -k2 > FILES_BEFORE_CLEANUP.txt
  if [[ "$3" == 1 ]]; then
    rm -f "$2"_trj.xyz "$2"_IRC_F_trj.xyz "$2"_IRC_B_trj.xyz "$2"_IRC_Full_trj.xyz ./*.tmp ./*.tmp.*
    [[ "$1" == irc ]] && rm -f freq.hess
    [[ "$1" == */c[0-9][0-9]_* ]] && rm -f "$2"_trj.dcd "$2".pdb ./*.densities start.xyz
  fi
  if command -v xz >/dev/null; then
    for f in "$2.out" ./*activeRegion_trj.xyz; do
      [[ -f "$f" ]] && [[ $(stat -c %s "$f") -gt 50000000 ]] && xz -T1 "$f"
    done
  fi
  du -sk . | awk '{{print "folder_kB_after_cleanup=" $1}}'
}}
alternate() {{   # step, GB of free disk needed per part: cycles K = 1..8 of environment step then core step, until
                # spa_stage.py alt_decide finds the alternation converged (exit 0) or not converged after cycle 8 (exit 1)
  local S=$1 K P C rc
  for K in 1 2 3 4 5 6 7 8; do
    for P in env core; do
      C=$(printf "c%02d_%s" "$K" "$P")
      PY alt_prep "$RD" "$S" "$K" "$P" || {{ echo "SP_${{S}}_FAIL alt_prep $S/$C"; quit 1; }}
      run_step "$S/$C" "$P" "$2"
      PY alt_check "$RD" "$S" "$K" "$P" && OK=1 || OK=0
      finish_step "$S/$C" "$P" "$OK"
      [[ $OK == 1 ]] || {{ PY alt_fail "$RD" "$S" "$C did not pass its checks"; echo "SP_${{S}}_FAIL"; quit 1; }}
    done
    PY alt_decide "$RD" "$S" "$K"; rc=$?
    [[ $rc == 0 ]] && {{ echo "SP_${{S}}_PASS"; return 0; }}
    [[ $rc == 3 ]] || {{ echo "SP_${{S}}_FAIL"; quit 1; }}
  done
  echo "SP_${{S}}_FAIL no decision after cycle 8"; quit 1
}}
echo "job={job} host=$(hostname) start=$(date)"
"""
def step_block(step, base, gb, check):
    s = ("PY prep \"$RD\" {step} || {{ echo SP_{step}_FAIL; quit 1; }}\n"
         "run_step {step} {base} {gb}\n"
         "PY {check} && OK=1 || OK=0\n"
         "finish_step {step} {base} $OK\n"
         "[[ $OK == 1 ]] && echo SP_{step}_PASS || {{ echo SP_{step}_FAIL; quit 1; }}\n")
    return s.format(step=step, base=base, gb=gb, check=check)
def alt_block(step, gb):
    return ("PY prep \"$RD\" {step} || {{ echo SP_{step}_FAIL; quit 1; }}\n"
            "alternate {step} {gb}\n").format(step=step, gb=gb)
JOBS = [  # job, memory GB, [(step, ORCA basename, free GB, check command)]; basename None: an alternation
    ("job_ts", 250, [("ts", None, 20, None)]),
    ("job_ts_freq_irc", 250, [("ts_freq", "freq", 3, 'check_freq "$RD" ts_freq ts'),
                              ("irc", "irc", 5, 'check_irc "$RD"')]),
    ("job_r_irc", 180, [("r_irc", None, 20, None)]),
    ("job_r_irc_freq", 180, [("r_irc_freq", "freq", 3, 'check_freq "$RD" r_irc_freq min')]),
    ("job_p_irc", 180, [("p_irc", None, 20, None)]),
    ("job_p_irc_freq", 180, [("p_irc_freq", "freq", 3, 'check_freq "$RD" p_irc_freq min')]),
    ("job_r0", 180, [("r0", None, 20, None)]),
    ("job_r0_freq", 180, [("r0_freq", "freq", 3, 'check_freq "$RD" r0_freq min')])]
SHORT = {"job_ts": "ts", "job_ts_freq_irc": "fi", "job_r_irc": "r", "job_r_irc_freq": "rf", "job_p_irc": "p",
         "job_p_irc_freq": "pf", "job_r0": "r0", "job_r0_freq": "r0f"}
SUBMIT = """#!/bin/bash
# submit.sh - the eight jobs of {name} (SP Amendment 1), with dependencies (PBS afterok: a job starts only if the one it
# needs ended with exit 0, i.e. every check passed). Run once, on hpc1: bash {rd}/submit.sh
set -euo pipefail
cd {rd}
[[ -e SUBMITTED ]] && {{ echo "STOP: already submitted (SUBMITTED lists the job ids)"; exit 1; }}
[[ -e SUBMITTED.partial ]] && {{ echo "STOP: an earlier submission stopped part-way; these jobs are queued (tell Claude, do not re-run):"; cat SUBMITTED.partial; exit 1; }}
[[ -e DO_NOT_SUBMIT ]] && {{ echo "STOP: DO_NOT_SUBMIT: $(cat DO_NOT_SUBMIT)"; exit 1; }}
for j in job_ts job_ts_freq_irc job_r_irc job_r_irc_freq job_p_irc job_p_irc_freq job_r0 job_r0_freq; do
  [[ -s $j.pbs ]] || {{ echo "STOP: missing $j.pbs"; exit 1; }}
  [[ -e $j.pbs.out ]] && {{ echo "STOP: $j.pbs.out exists"; exit 1; }}
done
sub() {{ local id; id=$(qsub "$@") || {{ echo "STOP: qsub $* failed; submitted so far (also in SUBMITTED.partial):" >&2; cat SUBMITTED.partial >&2; exit 1; }}; echo "$id"; }}
: > SUBMITTED.partial
A1=$(sub job_ts.pbs);                                  echo "job_ts $A1" >> SUBMITTED.partial
A2=$(sub -W depend=afterok:$A1 job_ts_freq_irc.pbs);  echo "job_ts_freq_irc $A2" >> SUBMITTED.partial
B1=$(sub -W depend=afterok:$A2 job_r_irc.pbs);        echo "job_r_irc $B1" >> SUBMITTED.partial
B2=$(sub -W depend=afterok:$B1 job_r_irc_freq.pbs);   echo "job_r_irc_freq $B2" >> SUBMITTED.partial
C1=$(sub -W depend=afterok:$A2 job_p_irc.pbs);        echo "job_p_irc $C1" >> SUBMITTED.partial
C2=$(sub -W depend=afterok:$C1 job_p_irc_freq.pbs);   echo "job_p_irc_freq $C2" >> SUBMITTED.partial
D1=$(sub job_r0.pbs);                                  echo "job_r0 $D1" >> SUBMITTED.partial
D2=$(sub -W depend=afterok:$D1 job_r0_freq.pbs);      echo "job_r0_freq $D2" >> SUBMITTED.partial
mv SUBMITTED.partial SUBMITTED
echo "{name}: submitted"; cat SUBMITTED
"""

# ---- per run: checks, region, starts
def fmt_rows(rows): return "".join("%s %s %s %s\n" % tuple(r) for r in rows)
WRITE = []   # (run name, {relative path: text}, [(relative link, target)], log lines of this run)
HDR = None
for f, rad, bulk, name in SPEC:
    if HDR is None: HDR = len(LOG)
    L0 = len(LOG)
    if f not in BAND: stop("frame %s is not one of the 30 frames of c4_bands.tsv" % f)
    b = BAND[f]; band = b["band"]; ci = int(b["ci"])
    say("")
    say("== %s: frame %s, radius %d A, bulk %s; band %s, climbing image %d" % (name, f, rad, bulk, band, ci))
    FD = ROOT + "/05_qmmm/19_ensemble/frame_%s" % f
    BD = ROOT + "/05_qmmm/%s/frame_%s" % (band, f)
    M = read_pdb(src(FD + "/frame_%s_CHA2.pdb" % f))
    if len(M) != NATOMS: stop("frame_%s_CHA2.pdb holds %d atoms" % (f, len(M)))
    bad = [i for i in range(NATOMS) if M[i][0] != ANAME[i]]
    if bad: stop("frame_%s_CHA2.pdb: %d atom names differ from the prmtop (first: atom %d %s vs %s)" % (
        f, len(bad), bad[0], M[bad[0]][0], ANAME[bad[0]]))
    if any(M[i][1] != "CHA" for i in QMS): stop("frame_%s_CHA2.pdb: atoms 6207-6230 are not CHA" % f)
    MX = [(float(r[2]), float(r[3]), float(r[4])) for r in M]
    MROW = [[EL[i], M[i][2], M[i][3], M[i][4]] for i in range(NATOMS)]
    tn = open(src(BD + "/neb.inp")).read()
    if braces(tn, "ActiveAtoms") != ACT12 or braces(tn, "QMAtoms") != QMS:
        stop("frame %s: the band's neb.inp does not have the R12 active atoms and the QM atoms 6207-6230" % f)
    XL, XS = read_allxyz(src(BD + "/neb_MEP.allxyz"))
    if len(XS) != 10 or any(s[2] != NATOMS for s in XS): stop("frame %s: neb_MEP.allxyz is not 10 structures of %d atoms" % (f, NATOMS))
    X0, XC = rows_of(XL, XS[0]), rows_of(XL, XS[ci])
    del XL, XS
    for nm, X in (("0", X0), ("ci", XC)):
        if any(X[i][0].lower() != EL[i].lower() for i in range(NATOMS)): stop("frame %s image %s: elements differ from the force field" % (f, nm))
    fix = max(max(dmax(fl(X[i]), MX[i]) for i in range(NATOMS) if i not in A12) for X in (X0, XC))
    if fix > TOL_BAND: stop("frame %s: an atom outside R12 differs between the band and the MD structure by %.1e A" % (f, fix))
    T = read_xyz_frames(src(BD + "/neb_MEP.QMRegion_trj.xyz"))
    if len(T) != 10: stop("frame %s: neb_MEP.QMRegion_trj.xyz is not 10 structures" % f)
    dq = max(max(dmax(fl(X[QM0 + j]), T[k][j]) for j in range(NQM)) for X, k in ((X0, 0), (XC, ci)))
    if dq > 2e-6: stop("frame %s: band images 0/ci substrate differs from the committed written path (%.1e A)" % (f, dq))
    say("sources: every file matches its sha256; MD structure names = prmtop; band images 0 and %d: elements ="
        " force field, atoms outside R12 = MD structure (max %.1e A), substrate = committed written path (max %.1e A);"
        " band neb.inp active atoms = R12" % (ci, fix, dq))
    # region: every residue with an atom within rad of an atom of the substrate, on the MD structure
    S = [MX[i] for i in QMS]
    r2 = rad * rad; res = set()
    for i in range(NATOMS):
        x = MX[i]
        for s in S:
            if (x[0] - s[0]) ** 2 + (x[1] - s[1]) ** 2 + (x[2] - s[2]) ** 2 <= r2: res.add(RESOF[i]); break
    region = sorted(i for r in res for i in range(RPTR[r], RPTR[r + 1]))
    RS = set(region)
    if not set(QMS) <= RS: stop("region lacks a substrate atom")
    labs = {}
    for r in res: labs[RLAB[r]] = labs.get(RLAB[r], 0) + 1
    nwat = labs.get("WAT", 0)
    say("region (whole residues with any atom within %d A of any of atoms 6207-6230, MD structure): %d atoms, %d"
        " residues (%d waters, %d others incl. the substrate)" % (rad, len(region), len(res), nwat, len(res) - nwat))
    say("  residues: " + " ".join("%s:%d" % (k, labs[k]) for k in sorted(labs)))
    say("  non-water residues (prmtop residue number, 1-based): " + " ".join(
        "%s%d" % (RLAB[r], r + 1) for r in sorted(res) if RLAB[r] != "WAT"))
    out12 = len(RS - A12); in12 = len(A12 - RS)
    say("  against R12 (the band's 1959 active atoms): %d region atoms were held fixed in the band (%d of them"
        " water atoms); %d R12 atoms lie outside the region" % (out12, sum(1 for i in RS - A12 if RLAB[RESOF[i]] == "WAT"), in12))
    # bulk and starts
    if bulk == "md":
        BULKROW = MROW; bdesc = "the MD structure (frame_%s_CHA2.pdb, coordinates as written)" % f
    else:
        BULKROW = XC; bdesc = "the band's climbing image (neb_MEP.allxyz image %d)" % ci
    start_ts = [XC[i] if i in RS else BULKROW[i] for i in range(NATOMS)]
    start_r0 = [X0[i] if i in RS else BULKROW[i] for i in range(NATOMS)]
    XT = [fl(r) for r in start_ts]; XR = [fl(r) for r in start_r0]
    XCF = [fl(r) for r in XC]; X0F = [fl(r) for r in X0]
    reset = [i for i in range(NATOMS) if i not in RS and dmax(XT[i], XCF[i]) > 1e-3]
    if reset:
        dd = sorted(dist(XT[i], XCF[i]) for i in reset)
        say("bulk: %s. Atoms outside the region whose start position differs from the climbing image: %d (max %.3f A,"
            " median %.3f A)" % (bdesc, len(reset), dd[-1], dd[len(dd) // 2]))
    else:
        say("bulk: %s. Every atom outside the region is at its climbing-image position" % bdesc)
    cross = [(i, j) for i, j in BONDS if (i in RS) != (j in RS)]
    if cross:
        dev_t = max(abs(dist(XT[i], XT[j]) - dist(XCF[i], XCF[j])) for i, j in cross)
        dev_r = max(abs(dist(XR[i], XR[j]) - dist(X0F[i], X0F[j])) for i, j in cross)
        say("covalent bonds across the region boundary: %d; largest bond-length change made by the splice: TS start"
            " %.3f A (against the climbing image), R0 start %.3f A (against image 0)" % (len(cross), dev_t, dev_r))
    else:
        say("covalent bonds across the region boundary: none")
    outside = [i for i in range(NATOMS) if i not in RS]
    for nm, X in (("TS start", XT), ("R0 start", XR)):
        pr = [(d, i, j) for i, j, d in neighbours(X, region, outside, 3.0) if (i, j) not in BONDED]
        dmin = min(pr) if pr else (9.9, -1, -1)
        say("closest non-bonded contact, region to outside, %s: %.3f A (atoms %d-%d)" % (nm, dmin[0], dmin[1], dmin[2]))
        if dmin[0] < CLASH: stop("%s: %s has a contact of %.3f A (below %.1f A)" % (name, nm, dmin[0], CLASH))
    for nm, X in (("TS start", XT), ("R0 start", XR), ("MD structure", MX)):
        say("substrate %s: breaking O-C (6214-6215) %.3f A, forming C-C (6207-6219) %.3f A" % (
            nm, dist(X[BRK[0]], X[BRK[1]]), dist(X[FRM[0]], X[FRM[1]])))
    rd = SD + "/" + name
    files = {}
    files["region.txt"] = "\n".join(map(str, region)) + "\n"
    files["ts/start.xyz"] = "%d\n%s: TS start - region atoms from band %s image %d (climbing image), other atoms from %s\n" % (
        NATOMS, name, band, ci, "the MD structure" if bulk == "md" else "the same image") + fmt_rows(start_ts)
    files["r0/start.xyz"] = "%d\n%s: R0 start - region atoms from band %s image 0 (Phase 1 reactant), other atoms from %s\n" % (
        NATOMS, name, band, "the MD structure" if bulk == "md" else "band image %d (climbing image)" % ci) + fmt_rows(start_r0)
    files["starts.sha256"] = "%s  ts/start.xyz\n%s  r0/start.xyz\n" % (shas(files["ts/start.xyz"]), shas(files["r0/start.xyz"]))
    env = [i for i in region if i not in set(QMS)]
    if len(env) + NQM != len(region): stop("region minus the substrate")
    say("alternation (SP Amendment 1): core = the 24 substrate atoms 6207-6230; environment = the region without them,"
        " %d atoms (active in the env steps)" % len(env))
    for st in ("ts", "r_irc", "p_irc", "r0"): files[st + "/env.inp"] = inp_env(env)
    files["ts/core.inp"] = inp_core_ts()
    for st in ("r_irc", "p_irc", "r0"): files[st + "/core.inp"] = inp_core_min()
    files["irc/irc.inp"] = inp_irc()
    for st in ("ts_freq", "r_irc_freq", "p_irc_freq", "r0_freq"): files[st + "/freq.inp"] = inp_freq()
    for job, mem, steps in JOBS:
        body = PBS_HEAD.format(name="sa%s%s%s" % (f, rad, SHORT[job]) + ("c" if bulk == "ci" else ""), mem=mem, rd=rd,
                               sd=SD, job=job)
        for st, base, gb, chk in steps: body += alt_block(st, gb) if base is None else step_block(st, base, gb, chk)
        files[job + ".pbs"] = body + "quit 0\n"
    files["submit.sh"] = SUBMIT.format(name=name, rd=rd)
    links = [(st + "/complex_solvated.ORCAFF.prms", "../../complex_solvated.ORCAFF.prms") for st in
             ("ts", "ts_freq", "irc", "r_irc", "p_irc", "r0", "r_irc_freq", "p_irc_freq", "r0_freq")]
    WRITE.append((name, files, links, LOG[L0:]))
    del X0, XC, M, MX, MROW, XT, XR, XCF, X0F, start_ts, start_r0

# ---- write
if OUT == "-":
    say("check only: every check passed, nothing written"); sys.exit(0)
for name, files, links, rl in WRITE:
    if os.path.exists(os.path.join(OUT, name)): stop(os.path.join(OUT, name) + " exists")
for name, files, links, rl in WRITE:
    od = os.path.join(OUT, name)
    for st in sorted(set(os.path.dirname(p) for p in files) | set(os.path.dirname(l) for l, t in links)):
        os.makedirs(os.path.join(od, st), exist_ok=True)
    for p in sorted(files):
        with open(os.path.join(od, p), "w") as fh: fh.write(files[p])
    for l, t in links: os.symlink(t, os.path.join(od, l))
    for p in ("submit.sh",) + tuple(j + ".pbs" for j, m, s in JOBS): os.chmod(os.path.join(od, p), 0o755)
    lg = LOG[:HDR] + rl
    lg.append("")
    lg.append("written in %s (sha256, path relative to the run folder):" % name)
    lg.extend("  %s  %s" % (shas(files[p]), p) for p in sorted(files))
    lg.extend("  link %s -> %s" % (l, t) for l, t in links)
    with open(os.path.join(od, "GENERATOR_LOG.txt"), "w") as fh: fh.write("\n".join(lg) + "\n")
    print("written: %s (%d files, %d links)" % (od, len(files) + 1, len(links)))
