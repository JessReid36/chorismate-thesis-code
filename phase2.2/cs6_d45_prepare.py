#!/usr/bin/env python3
"""cs6_d45_prepare.py - DIAGNOSTIC 45: the approximations of the coupled-cluster benchmark, measured.

14 single points on the bare substrate of frame 42436, no charges, in one job.

  the triples       T1_R, T1_TS    the committed def2-SVPD benchmark input with DLPNO-CCSD(T1), the
                                   iterative triples, in place of DLPNO-CCSD(T)
  the PNO threshold P6_R, P6_TS    the committed input with TCutPNO 1e-6
                    P8_R, P8_TS    the committed input with TCutPNO 1e-8
                                   (the committed runs, d38/N and d41/T, are TightPNO: 1e-7)
  the basis set     M_SV, M_TZ, M_TP, M_QZ, each at R and TS: RI-MP2 at def2-SVPD, def2-TZVPD,
                                   def2-TZVPPD and def2-QZVPPD

Nothing about a run is typed in to be trusted. THE TEMPLATE is the committed input of the def2-SVPD
benchmark, d38/N/job.inp, which d41/T/job.inp repeats byte for byte; the strings typed in below
are a check on what is read, and nothing is written unless the file holds exactly those lines.
Every input here is that file with ONE thing changed, and the check at the end is made on the
files that were written, against the committed file and not against the function that wrote them.
THE GEOMETRIES are d33a's, byte for byte.

select() and everything it uses need the standard library only, so the reader imports them and
rebuilds every input on its own.

Usage:  cs6_d45_prepare.py <results-repo> <out-dir> <hpc-dir>
"""
import hashlib, os, shutil, sys

S6 = "phase2.2/cs_stage6"
GEOM = {"R": "d33a/sp_42436_R/substrate.xyz", "TS": "d33a/sp_42436_TS/substrate.xyz"}
TEMPLATE = "d38/N/job.inp"                    # the committed input of the def2-SVPD benchmark, reactant
TEMPLATE_TS = "d41/T/job.inp"                 # the same file, transition state
TEMPLATE_LINES = ["! DLPNO-CCSD(T) def2-SVPD AutoAux TightSCF TightPNO", "%maxcore 6000", "%pal nprocs 16 end",
                  "%scf MaxIter 300 end", "* xyzfile -2 1 substrate.xyz", ""]
STATES = ("R", "TS")
T1_LINE = "! DLPNO-CCSD(T1) def2-SVPD AutoAux TightSCF TightPNO"
PNO_BLOCK = {"P6": ["%mdci", "  TCutPNO 1e-6", "end"], "P8": ["%mdci", "  TCutPNO 1e-8", "end"]}
TCUT = {"P6": "1.000e-06", "P8": "1.000e-08", "T1": "1.000e-07"}      # as ORCA prints TCutPNO
MP2 = (("M_SV", "def2-SVPD"), ("M_TZ", "def2-TZVPD"), ("M_TP", "def2-TZVPPD"), ("M_QZ", "def2-QZVPPD"))
# the order the job makes them in: the cheap and the certain first, so that a job that is cut
# short has lost the least
ORDER = ["M_SV", "M_TZ", "M_TP", "P6", "T1", "M_QZ", "P8"]
NRUN, NSUM = 14, 31                           # 14 x 2 files, RUNS.tsv, SOURCES.sha256, run.pbs
CONTROLS = ("d38/N", "d41/T", "d36a2/B_R", "d36a2/B_TS")
JOB, CORES, MEM, WALL = "cm_d45", 16, "120gb", "48:00:00"
WANT_GB, MIN_GB = 250, 120
HPC = "/home/18660916/system_development/phase2.2/cs_stage6/d45"


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


def sources():
    """every committed file this script or the reader reads, relative to phase2.2/cs_stage6"""
    out = [GEOM["R"], GEOM["TS"], "STAGE6_DIAGNOSTIC36_RESULT_20261114.txt", "STAGE6_DIAGNOSTIC41_20261115.txt",
           "STAGE6_DIAGNOSTIC41_RESULT_20261115.txt", "STAGE6_AMENDMENT7_20261117.txt", "d37.tsv", "d37a1.tsv"]
    for c in CONTROLS:
        out += ["%s/job.inp" % c, "%s/job.out" % c, "%s/substrate.xyz" % c]
    return out


def select(s6):
    """What the committed files are required to be. Raises Refused where one is not."""
    for f in sources():
        if not os.path.isfile(os.path.join(s6, f)):
            stop("%s, which this diagnostic reads, is missing." % f)
    for f in (TEMPLATE, TEMPLATE_TS):
        if open(os.path.join(s6, f)).read().split("\n") != TEMPLATE_LINES:
            stop("%s is not the committed input of the def2-SVPD benchmark, line for line." % f)
    for c, st in (("d38/N", "R"), ("d41/T", "TS"), ("d36a2/B_R", "R"), ("d36a2/B_TS", "TS")):
        if sha(os.path.join(s6, c, "substrate.xyz")) != sha(os.path.join(s6, GEOM[st])):
            stop("%s/substrate.xyz is not frame 42436's committed %s geometry." % (c, st))
    return dict(tl=list(TEMPLATE_LINES))


def build_input(kind, basis=None):
    """The text of one input: the committed template with one thing changed."""
    t = list(TEMPLATE_LINES)
    if kind == "T1":
        t[0] = T1_LINE
    elif kind in PNO_BLOCK:
        t = t[:4] + PNO_BLOCK[kind] + t[4:]
    elif kind == "MP2":
        t[0] = "! RI-MP2 %s AutoAux TightSCF" % basis
    else:
        stop("no such kind of run: %s" % kind)
    return "\n".join(t)


def runs_of():
    """every run, in the order it is made: (run, kind, basis, state)"""
    bas = dict(MP2)
    out = []
    for k in ORDER:
        for st in STATES:
            if k in bas:
                out.append(("%s_%s" % (k, st), "MP2", bas[k], st))
            else:
                out.append(("%s_%s" % (k, st), k, "def2-SVPD", st))
    if len(out) != NRUN or len(set(r[0] for r in out)) != NRUN:
        stop("%d runs were listed, not %d distinct ones." % (len(out), NRUN))
    return out


def runs_text(runs):
    """d45/RUNS.tsv"""
    what = {"T1": "the committed input with DLPNO-CCSD(T1)", "P6": "the committed input with TCutPNO 1e-6",
            "P8": "the committed input with TCutPNO 1e-8", "MP2": "the committed input with RI-MP2 at this basis set, and no TightPNO"}
    L = ["run\tkind\tbasis\tstate\tmethod_line\twhat_differs_from_%s" % TEMPLATE]
    for run, kind, basis, st in runs:
        L.append("%s\t%s\t%s\t%s\t%s\t%s" % (run, kind, basis, st, build_input(kind, basis).split("\n")[0], what[kind]))
    return "\n".join(L) + "\n"


def pbs_text(runs, hpc):
    """the job script: DIAGNOSTIC 41's, for fourteen runs"""
    names = " ".join(r[0] for r in runs)
    return r"""#!/bin/bash
#PBS -N %(job)s
#PBS -l select=1:ncpus=%(cores)d:mem=%(mem)s
#PBS -l walltime=%(wall)s
#PBS -m ae
#PBS -M 18660916@sun.ac.za
#PBS -j oe
#PBS -o %(hpc)s/run.pbs.out
# Stage 6 DIAGNOSTIC 45: the approximations of the coupled-cluster benchmark, measured.
#   Fourteen single points on the bare substrate of frame 42436, one after another, in the
#   order of RUNS.tsv: the cheap and the certain first.
# No -q: the allocation is automatic. The longest runs here have never been made, so the
#   walltime is %(wall)s and not 24 hours. Restartable: submitted again, it skips every run
#   that ended normally. A run that fails does not stop the ones after it.
#
# SCRATCH. This is DIAGNOSTIC 41's script. Every run is made in node-local scratch and its
#   text output is copied home; nothing is computed in /home. The def2-SVPD benchmark needed
#   31 to 46 GB of scratch (D38 Erratum 1, D41). What a tighter PNO threshold needs is not
#   measured, so %(want)d GB is asked for and %(min)d GB is the least a candidate may offer.
#   The peak is sampled every 60 s WHILE ORCA runs. Scratch is removed after each run.
export PATH="/apps/openmpi/4.1.1/bin:$PATH"
export LD_LIBRARY_PATH="/home/apps2/ORCA/6.0.1/lib:/apps/openmpi/4.1.1/lib:/apps/mambaforge/envs/medaka/lib:${LD_LIBRARY_PATH:-}"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
ORCA=/home/apps2/ORCA/6.0.1
H=%(hpc)s
P="$H/PROGRESS.txt"
PK="$H/.peak_mb"
RUNS="%(names)s"

cd "$H" || exit 1
if mkdir .running 2>/dev/null; then :; else
  echo "FAIL: another job holds .running in $H"; exit 1
fi

MARK="d45.${PBS_JOBID:-nojobid}"
WANT_GB=%(want)d
MIN_GB=%(min)d
S=""
AVAIL=""
SMARGIN=""

pick() {
  for c in /scratch-small-local/18660916 /scratch/18660916 "${TMPDIR:-}"; do
    if [ -z "$c" ]; then continue; fi
    t="$c/$MARK"
    if mkdir -p "$t" 2>/dev/null; then :; else
      echo "   candidate $c  mkdir failed"; continue; fi
    if echo ok > "$t/.probe" 2>/dev/null; then :; else
      echo "   candidate $c  created but not writable"; rmdir "$t" 2>/dev/null; continue; fi
    rm -f "$t/.probe"
    a=$(df -P -BG "$t" 2>/dev/null | awk 'NR==2{sub(/G$/,"",$4); print $4}')
    case "${a:-}" in
      ''|*[!0-9]*) echo "   candidate $c  free space unreadable"; rmdir "$t" 2>/dev/null; continue ;;
    esac
    if [ "$a" -ge "$1" ]; then S="$t"; AVAIL="$a"; return 0; fi
    echo "   candidate $c  ${a}G free, under the ${1}G asked for"
    rmdir "$t" 2>/dev/null       # a refused candidate is left as it was found
  done
  return 1
}

cleanup() {
  kill "${WATCH:-}" 2>/dev/null
  if [ -n "$S" ] && [ -d "$S" ]; then
    for r in $RUNS; do
      if [ -f "$S/$r/job.out" ] && [ -d "$H/$r" ]; then
        cp -f "$S/$r/job.out" "$H/$r/job.out" 2>/dev/null
      fi
    done
    echo "scratch at exit: $(du -sh "$S" 2>/dev/null | cut -f1)"
    case "$S" in
      */d45.*) rm -rf "$S" ;;
      *) echo "REFUSING to remove an unexpected scratch path: $S" ;;
    esac
  fi
  rmdir "$H/.running" 2>/dev/null
}
trap cleanup EXIT

echo "host=$(hostname) start=$(date) jobid=${PBS_JOBID:-none}"
echo "cores=$(nproc) ; memory:"; free -g | head -2

echo "choosing scratch, in order of preference:"
if pick "$WANT_GB"; then
  SMARGIN="full"
elif pick "$MIN_GB"; then
  SMARGIN="REDUCED"
  echo "   NOTE: reduced margin. What the tighter PNO threshold needs is not measured."
else
  echo "FAIL: no candidate on $(hostname) is writable with at least ${MIN_GB}G free."
  echo "      Nothing is wrong with the calculation. Resubmit and PBS will usually"
  echo "      place the job on another node."
  exit 1
fi
echo "R45.8 scratch=$S avail=${AVAIL}G margin=$SMARGIN"
echo 0 > "$PK"

printf '%%s  d45 on %%s, scratch %%s, %%sG free, sampling every %%ds\n' \
  "$(date '+%%F %%T')" "$(hostname)" "$S" "$AVAIL" 60 >> "$P"

# The watcher carries the maximum it sees while ORCA is running and writes it to $PK, because
# a measurement taken after ORCA returns records the residue (D38's result, Erratum 1). It is
# set back to zero before each run, so each run has its own maximum.
( sleep 20
  while [ -d "$S" ]; do
    u=$(du -sm "$S" 2>/dev/null | cut -f1)
    case "${u:-}" in
      ''|*[!0-9]*) u=0 ;;
    esac
    M=$(cat "$PK" 2>/dev/null || echo 0)
    case "${M:-}" in
      ''|*[!0-9]*) M=0 ;;
    esac
    if [ "$u" -gt "$M" ]; then echo "$u" > "$PK"; fi
    printf '%%s  used=%%sM  running=%%s\n' "$(date '+%%F %%T')" "$u" "$(ls "$S" 2>/dev/null | tr '\n' ' ')" >> "$P"
    sleep 60
  done ) &
WATCH=$!

for r in $RUNS; do
  if [ -s "$H/$r/job.out" ] && grep -q "ORCA TERMINATED NORMALLY" "$H/$r/job.out"; then
    echo "$r done before"; continue
  fi
  echo 0 > "$PK"
  mkdir -p "$S/$r" || { echo "$r scratch directory could not be made"; continue; }
  cp -p "$H/$r/job.inp" "$H/$r/substrate.xyz" "$S/$r/" || { echo "$r inputs could not be copied to scratch"; continue; }
  ( cd "$S/$r" && "$ORCA/orca" job.inp > job.out 2>&1 </dev/null )
  rc=$?
  cp -f "$S/$r/job.out" "$H/$r/job.out" 2>/dev/null
  if [ -f "$S/$r/job.property.txt" ]; then cp -f "$S/$r/job.property.txt" "$H/$r/job.property.txt"; fi
  echo "$r orca_exit=$rc normal=$(grep -c 'ORCA TERMINATED NORMALLY' "$H/$r/job.out") peak_mb=$(cat "$PK" 2>/dev/null || echo 0) $(date)"
  printf '%%s  %%s orca_exit=%%s\n' "$(date '+%%F %%T')" "$r" "$rc" >> "$P"
  case "$S" in
    */d45.*) rm -rf "$S/$r" ;;
    *) echo "REFUSING to remove an unexpected scratch path: $S/$r" ;;
  esac
done
echo "end=$(date)"
echo "finished: $(grep -l 'ORCA TERMINATED NORMALLY' "$H"/*/job.out 2>/dev/null | wc -l) of %(n)d"
""" % dict(job=JOB, cores=CORES, mem=MEM, wall=WALL, hpc=hpc, names=names, want=WANT_GB, min=MIN_GB, n=len(runs))


def main():
    if len(sys.argv) != 4:
        stop("usage: cs6_d45_prepare.py <results-repo> <out-dir> <hpc-dir>")
    res, out, hpc = sys.argv[1:4]
    s6 = os.path.join(res, S6)
    if os.path.exists(out) and os.listdir(out):
        stop("%s exists and is not empty; this script writes into a new directory only." % out)
    select(s6)
    runs = runs_of()
    os.makedirs(out, exist_ok=True)
    sums = []
    for run, kind, basis, st in runs:
        d = os.path.join(out, run)
        os.makedirs(d)
        shutil.copyfile(os.path.join(s6, GEOM[st]), os.path.join(d, "substrate.xyz"))
        with open(os.path.join(d, "job.inp"), "w") as f:
            f.write(build_input(kind, basis))
        for fn in sorted(os.listdir(d)):
            sums.append((sha(os.path.join(d, fn)), "%s/%s" % (run, fn)))

    # ---- the checks, on what was actually written, against the committed input and not against build_input
    src = open(os.path.join(s6, TEMPLATE)).read().split("\n")
    for run, kind, basis, st in runs:
        d = os.path.join(out, run)
        if sorted(os.listdir(d)) != ["job.inp", "substrate.xyz"]:
            stop("%s does not hold exactly job.inp and substrate.xyz." % run)
        if sha(os.path.join(d, "substrate.xyz")) != sha(os.path.join(s6, GEOM[st])):
            stop("%s/substrate.xyz is not frame 42436's committed %s geometry." % (run, st))
        b = open(os.path.join(d, "job.inp")).read().split("\n")
        if kind == "T1":
            if b[1:] != src[1:] or b[0] != src[0].replace("DLPNO-CCSD(T)", "DLPNO-CCSD(T1)") or b[0] == src[0]:
                stop("%s is not the committed input with DLPNO-CCSD(T1) and nothing else changed." % run)
        elif kind in PNO_BLOCK:
            if b[:4] != src[:4] or b[4:7] != PNO_BLOCK[kind] or b[7:] != src[4:]:
                stop("%s is not the committed input with the TCutPNO block added and nothing else changed." % run)
        else:
            w = src[0].split()
            if b[1:] != src[1:] or b[0].split() != ["!", "RI-MP2", basis] + [x for x in w[3:] if x != "TightPNO"] \
                    or w[1] != "DLPNO-CCSD(T)" or w[2] != "def2-SVPD":
                stop("%s is not the committed input with RI-MP2 and %s in the method line and nothing else changed." % (run, basis))
        if "pointcharges" in "\n".join(b).lower():
            stop("%s carries point charges; D45 has none." % run)

    with open(os.path.join(out, "RUNS.tsv"), "w") as f:
        f.write(runs_text(runs))
    with open(os.path.join(out, "SOURCES.sha256"), "w") as f:
        for s in sources():
            f.write("%s  %s\n" % (sha(os.path.join(s6, s)), s))
    p = os.path.join(out, "run.pbs")
    with open(p, "w") as f:
        f.write(pbs_text(runs, hpc))
    os.chmod(p, 0o755)
    for fn in ("RUNS.tsv", "SOURCES.sha256", "run.pbs"):
        sums.append((sha(os.path.join(out, fn)), fn))
    if len(sums) != NSUM:
        stop("%d files were written, not %d." % (len(sums), NSUM))
    with open(os.path.join(out, "INPUT_CHECKSUMS.sha256"), "w") as f:
        for h, n in sums:
            f.write("%s  %s\n" % (h, n))
    for run, kind, basis, st in runs:
        print("   %-8s %-4s %-12s %s" % (run, kind, basis, build_input(kind, basis).split("\n")[0]))
    print("   %d run directories, %d files checksummed, one job, %d cores, %s, walltime %s" % (NRUN, len(sums), CORES, MEM, WALL))


if __name__ == "__main__":
    try:
        main()
    except Refused as e:
        sys.stderr.write("REFUSED  %s\n" % e)
        raise SystemExit(2)
