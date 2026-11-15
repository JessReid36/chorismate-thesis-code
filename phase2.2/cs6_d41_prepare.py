#!/usr/bin/env python3
"""cs6_d41_prepare.py - DIAGNOSTIC 41: the correlated barrier at the production basis.

One single point: the bare substrate at frame 42436's TRANSITION STATE, DLPNO-CCSD(T)/def2-SVPD.
Its input is DIAGNOSTIC 38's run N, byte for byte; only the geometry file beside it differs. The
reactant leg is not run: it is d38/N/job.out, committed, and its sha256 is written out here so
that the reader can require the same file.

Usage:  cs6_d41_prepare.py <results-repo> <out-dir> <hpc-dir>
"""
import hashlib, os, shutil, sys

S6 = "phase2.2/cs_stage6"
GEOM = {"R": "d33a/sp_42436_R/substrate.xyz", "TS": "d33a/sp_42436_TS/substrate.xyz"}
SRC = "d38/N"
LINE = "! DLPNO-CCSD(T) def2-SVPD AutoAux TightSCF TightPNO"
WANT_GB, MIN_GB = 100, 60
WALL, MEM, NCPU = "24:00:00", "120gb", 16

PBS = r'''#!/bin/bash
#PBS -N cm_d41
#PBS -l select=1:ncpus=@NCPU@:mem=@MEM@
#PBS -l walltime=@WALL@
#PBS -m ae
#PBS -M 18660916@sun.ac.za
#PBS -j oe
#PBS -o @HPC@/run.pbs.out
# DIAGNOSTIC 41: one single point, the bare substrate at frame 42436's transition
#   state, DLPNO-CCSD(T)/def2-SVPD. The input is DIAGNOSTIC 38's run N byte for
#   byte; the reactant leg is d38/N, committed, and is not run again.
# walltime @WALL@ routes to day automatically. No -q.
#
# SCRATCH. The job runs in node-local scratch and copies its output home. The
#   script is DIAGNOSTIC 39's, cut down to one closed-shell run. D38's result
#   Erratum 1 bounds a closed-shell def2-SVPD run of the reactant between 31 and
#   43 GB. The transition state is not measured; at def2-TZVPD it needed more
#   than the reactant (80 GB against 74, D36). So @WANT_GB@ GB is asked for and
#   @MIN_GB@ GB is the least a candidate may offer, which rules out $TMPDIR here
#   (38 to 43 GB) on purpose. The peak is sampled every 60 s WHILE ORCA runs.
export PATH="/apps/openmpi/4.1.1/bin:$PATH"
export LD_LIBRARY_PATH="/home/apps2/ORCA/6.0.1/lib:/apps/openmpi/4.1.1/lib:/apps/mambaforge/envs/medaka/lib:${LD_LIBRARY_PATH:-}"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
ORCA=/home/apps2/ORCA/6.0.1
H=@HPC@
P="$H/PROGRESS.txt"
PK="$H/.peak_mb"
RUNS="T"

cd "$H" || exit 1
if mkdir .running 2>/dev/null; then :; else
  echo "FAIL: another job holds .running in $H"; exit 1
fi

MARK="d41.${PBS_JOBID:-nojobid}"
WANT_GB=@WANT_GB@
MIN_GB=@MIN_GB@
S=""
AVAIL=""
SMARGIN=""
RESIDUE=0
echo 0 > "$PK"

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
      */d41.*) rm -rf "$S" ;;
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
  echo "   NOTE: reduced margin. The reactant needed 31 to 43 GB at this basis and"
  echo "   the transition state is not measured. R41.6 reports the sampled maximum."
else
  echo "FAIL: no candidate on $(hostname) is writable with at least ${MIN_GB}G free."
  echo "      Nothing is wrong with the calculation. Resubmit and PBS will usually"
  echo "      place the job on another node. \$TMPDIR is deliberately too small to"
  echo "      be accepted here."
  exit 1
fi
echo "R41.6 scratch=$S avail=${AVAIL}G margin=$SMARGIN"

printf '%s  d41 on %s, scratch %s, %sG free, sampling every %ds\n' \
  "$(date '+%F %T')" "$(hostname)" "$S" "$AVAIL" 60 > "$P"

# The watcher carries the maximum it sees while ORCA is running and writes it to
# $PK, because a measurement taken after ORCA returns records the residue (D38's
# result, Erratum 1).
( sleep 20
  M=0
  while [ -d "$S" ]; do
    u=$(du -sm "$S" 2>/dev/null | cut -f1)
    case "${u:-}" in
      ''|*[!0-9]*) u=0 ;;
    esac
    if [ "$u" -gt "$M" ]; then M="$u"; echo "$M" > "$PK"; fi
    printf '%s  used=%sM  max=%sM  storage=%s  running=%s\n' \
      "$(date '+%F %T')" "$u" "$M" \
      "$(grep -hm1 'PNO-Integral Storage' "$S"/*/job.out 2>/dev/null | sed 's/.*\.\.\.[[:space:]]*//' | head -1)" \
      "$(ls "$S" 2>/dev/null | tr '\n' ' ')" >> "$P"
    sleep 60
  done ) &
WATCH=$!

RC=0
for r in $RUNS; do
  if [ -s "$H/$r/job.out" ] && grep -q "ORCA TERMINATED NORMALLY" "$H/$r/job.out"; then
    echo "=== $r done before, skipping"; continue
  fi
  echo "=== $r  start $(date)"
  mkdir -p "$S/$r" || { RC=1; break; }
  cp -p "$H/$r/job.inp" "$H/$r/substrate.xyz" "$S/$r/" || { RC=1; break; }
  ( cd "$S/$r" && "$ORCA/orca" job.inp > job.out 2>&1 </dev/null )
  rc=$?
  cp -f "$S/$r/job.out" "$H/$r/job.out" 2>/dev/null
  for f in job.property.txt job.gbw; do
    if [ -f "$S/$r/$f" ]; then cp -f "$S/$r/$f" "$H/$r/$f"; fi
  done
  echo "$r orca_exit=$rc normal=$(grep -c 'ORCA TERMINATED NORMALLY' "$H/$r/job.out")"
  echo "$r storage: $(grep -m1 'PNO-Integral Storage' "$H/$r/job.out")"
  echo "$r E(0):    $(grep -m1 'E(0)' "$H/$r/job.out")"
  echo "$r E(CCSD): $(grep -m1 'E(CCSD)' "$H/$r/job.out")"
  echo "$r CCSD(T): $(grep -m1 'E(CCSD(T))' "$H/$r/job.out")"
  grep -m1 "TOTAL RUN TIME" "$H/$r/job.out"
  printf '%s  %s orca_exit=%s\n' "$(date '+%F %T')" "$r" "$rc" >> "$P"
  # The residue, kept under its own name. This is NOT the peak; $PK holds that.
  u=$(du -sm "$S" 2>/dev/null | cut -f1)
  case "${u:-}" in
    ''|*[!0-9]*) u=0 ;;
  esac
  if [ "$u" -gt "$RESIDUE" ]; then RESIDUE="$u"; fi
  echo "$r residue after this run: ${u} MB  (not the peak; see R41.6)"
  case "$S" in
    */d41.*) rm -rf "$S/$r" ;;
    *) echo "REFUSING to remove an unexpected scratch path: $S/$r" ;;
  esac
  if [ "$rc" -ne 0 ]; then RC=$rc; fi
done
PEAK=$(cat "$PK" 2>/dev/null || echo 0)
echo "R41.6 sampled maximum while running: ${PEAK} MB, every 60s, on $S at the $SMARGIN margin"
echo "R41.6 largest residue after the run:  ${RESIDUE} MB"
echo "R41.6 the sampled maximum is a lower bound on the true peak"
echo "all runs finished, worst exit=$RC $(date)"
exit $RC
'''


def stop(m):
    sys.stderr.write("REFUSED  %s\n" % m)
    raise SystemExit(2)


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 16), b""):
            h.update(b)
    return h.hexdigest()


def main():
    if len(sys.argv) != 4:
        stop("usage: cs6_d41_prepare.py <results-repo> <out-dir> <hpc-dir>")
    res, out, hpc = sys.argv[1], sys.argv[2], sys.argv[3]
    s6 = os.path.join(res, S6)
    if os.path.exists(out) and os.listdir(out):
        stop("%s exists and is not empty; this script writes into a new directory only." % out)
    src = os.path.join(s6, SRC)
    for fn in ("job.inp", "job.out", "substrate.xyz"):
        if not os.path.isfile(os.path.join(src, fn)):
            stop("%s/%s is missing; it is the reactant leg and the template." % (SRC, fn))
    t = open(os.path.join(src, "job.inp")).read()
    if t.splitlines()[0] != LINE:
        stop("%s/job.inp does not carry the method line this diagnostic names: %s" % (SRC, t.splitlines()[0]))
    if "* xyzfile -2 1 substrate.xyz" not in t or "pointcharges" in t.lower():
        stop("%s/job.inp is not a charge-free dianion input reading substrate.xyz." % SRC)
    if sha(os.path.join(src, "substrate.xyz")) != sha(os.path.join(s6, GEOM["R"])):
        stop("%s/substrate.xyz is not frame 42436's committed reactant." % SRC)
    o = open(os.path.join(src, "job.out"), errors="replace").read()
    if "****ORCA TERMINATED NORMALLY****" not in o or "E(CCSD(T))" not in o:
        stop("%s/job.out is not a normally terminated coupled-cluster output." % SRC)
    if sha(os.path.join(s6, GEOM["R"])) == sha(os.path.join(s6, GEOM["TS"])):
        stop("the two committed geometries are the same file.")

    d = os.path.join(out, "T")
    os.makedirs(d)
    shutil.copyfile(os.path.join(src, "job.inp"), os.path.join(d, "job.inp"))
    shutil.copyfile(os.path.join(s6, GEOM["TS"]), os.path.join(d, "substrate.xyz"))
    if sha(os.path.join(d, "job.inp")) != sha(os.path.join(src, "job.inp")):
        stop("the input did not copy byte for byte.")
    if sha(os.path.join(d, "substrate.xyz")) != sha(os.path.join(s6, GEOM["TS"])):
        stop("the transition-state geometry did not copy byte for byte.")

    pbs = os.path.join(out, "run.pbs")
    text = PBS
    for k, v in (("@HPC@", hpc), ("@WALL@", WALL), ("@MEM@", MEM), ("@NCPU@", str(NCPU)), ("@WANT_GB@", str(WANT_GB)), ("@MIN_GB@", str(MIN_GB))):
        if k not in text:
            stop("the job script template has no %s" % k)
        text = text.replace(k, v)
    if "@" in text.replace("18660916@sun.ac.za", ""):
        stop("the job script still holds an unfilled field.")
    with open(pbs, "w") as f:
        f.write(text)
    os.chmod(pbs, 0o755)
    sums = [(sha(os.path.join(d, fn)), "T/%s" % fn) for fn in sorted(os.listdir(d))] + [(sha(pbs), "run.pbs")]
    with open(os.path.join(out, "INPUT_CHECKSUMS.sha256"), "w") as f:
        for h, n in sums:
            f.write("%s  %s\n" % (h, n))
    # the reactant leg, by hash, for the reader: a path relative to cs_stage6
    with open(os.path.join(out, "REACTANT_LEG.sha256"), "w") as f:
        f.write("%s  %s/job.out\n" % (sha(os.path.join(src, "job.out")), SRC))
    print("   T   %s   (input of %s, geometry %s)" % (LINE, SRC, GEOM["TS"]))
    print("   reactant leg %s/job.out  sha256 %s" % (SRC, sha(os.path.join(src, "job.out"))))
    print("   1 run directory, %d files checksummed, %d cores, %s, walltime %s, scratch %d GB asked, %d GB least" % (len(sums), NCPU, MEM, WALL, WANT_GB, MIN_GB))


if __name__ == "__main__":
    main()
