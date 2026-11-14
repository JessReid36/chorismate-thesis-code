#!/usr/bin/env python3
"""DIAGNOSTIC 36, AMENDMENT 2.  Attempt 3 of arm B.

Attempt 2's input is copied and changed in exactly two ways:

  1. the %mdci block of the ORCA 6.0 manual, section 7.14.4, is inserted, so ORCA is
     asked to hold the much-used DLPNO integrals in shared memory;
  2. %maxcore 12000 becomes %maxcore 8000, leaving room inside the same 240 GB PBS
     request for that shared segment, because the manual does not say whether the
     segment is counted inside %maxcore or allocated on top of it.

Both are storage directives. Neither decides what is computed, and G36.8 is the
control that measures that rather than asserting it.

The job script is where the real change is. Attempt 2 ran with its working directory
on /home, which is NFS, shared and 97 per cent full, and wrote 116 GB of intermediates
there. Attempt 3 stages into /scratch on the compute node - 1.9 TB, 1.5 TB free,
local - runs there, and copies the output back. The manual says ORCA falls back to
on-disk storage when memory is short, so the fallback is planned for instead of hoped
against.

Refuses, and writes nothing, if:
  - the method line is not the pre-registered one
  - %maxcore 12000 is not found exactly once
  - the %mdci or StorageType keyword is already present
  - the geometry is not frame 42436's committed substrate, byte for byte
  - the new input does not round-trip back to attempt 2's, byte for byte

Usage: cs6_d36_amend2_prepare.py <results-repo> <out-dir> <hpc-dir>
"""
import hashlib, os, sys

S6 = "phase2.2/cs_stage6"
SRC = "d36a1"                       # attempt 2's committed inputs
METHOD = "! DLPNO-CCSD(T) def2-TZVPD AutoAux TightSCF TightPNO"
OLD_MAXCORE = "%maxcore 12000"
NEW_MAXCORE = "%maxcore 8000"
MDCI = ["%mdci", "     StorageType Shared", "end"]   # manual 7.14.4, verbatim layout
ANCHOR = "* xyzfile"
NEED_GB = 300
LEGS = ("R", "TS")


def stop(msg):
    print("REFUSED  %s" % msg)
    sys.exit(1)


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 16), b""):
            h.update(b)
    return h.hexdigest()


def transform(src_text):
    """attempt 2's input -> attempt 3's input, and nothing else."""
    lines = src_text.splitlines()
    if lines.count(OLD_MAXCORE) != 1:
        stop("%r appears %d times in the input, expected once."
             % (OLD_MAXCORE, lines.count(OLD_MAXCORE)))
    if any("mdci" in l.lower() or "storagetype" in l.lower() for l in lines):
        stop("the input already mentions mdci or StorageType; this script would "
             "duplicate it.")
    anchor = [i for i, l in enumerate(lines) if l.startswith(ANCHOR)]
    if len(anchor) != 1:
        stop("found %d %r lines, expected exactly one." % (len(anchor), ANCHOR))
    out = []
    for i, l in enumerate(lines):
        if i == anchor[0]:
            out.extend(MDCI)
        out.append(NEW_MAXCORE if l == OLD_MAXCORE else l)
    return "\n".join(out) + "\n"


def untransform(new_text):
    """attempt 3's input -> attempt 2's, for the round-trip check."""
    lines = new_text.splitlines()
    for k in MDCI:
        if lines.count(k) != 1:
            stop("round trip: %r appears %d times in the generated input."
                 % (k, lines.count(k)))
    out = [OLD_MAXCORE if l == NEW_MAXCORE else l for l in lines if l not in MDCI]
    return "\n".join(out) + "\n"


PBS = """#!/bin/bash
#PBS -N cm_d36a2_%(leg)s
#PBS -l select=1:ncpus=16:mem=240gb
#PBS -l walltime=168:00:00
#PBS -m ae
#PBS -M 18660916@sun.ac.za
#PBS -j oe
#PBS -o %(hpc)s/run_B_%(leg)s.pbs.out
# DIAGNOSTIC 36 AMENDMENT 2: arm B, %(leg)s leg, attempt 3.
# Attempt 2 failed at exit 1, not at the walltime. ORCA reported PNO-Integral
#   Storage ON DISK again, because %%maxcore is not the control for it; the control
#   is the %%mdci StorageType block, which was never set. The intermediates - 74 GB
#   for R and 42 GB for TS - went to /home, which is NFS and 97%% full, and both
#   jobs stopped 96 seconds apart on different nodes.
# This script does two things differently. The input now asks for shared-memory
#   storage. And the run happens in /scratch on the compute node, which is local
#   and has 1.5 TB free, so a fallback to disk is survivable instead of fatal.
# The method line, AutoAux, TightPNO, def2-TZVPD and the geometry are unchanged.
# walltime 168:00:00 routes to week automatically. No -q.
export PATH="/apps/openmpi/4.1.1/bin:$PATH"
export LD_LIBRARY_PATH="/home/apps2/ORCA/6.0.1/lib:/apps/openmpi/4.1.1/lib:/apps/mambaforge/envs/medaka/lib:${LD_LIBRARY_PATH:-}"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
ORCA=/home/apps2/ORCA/6.0.1
H=%(hpc)s/B_%(leg)s
P="$H/PROGRESS.txt"
NEED_GB=%(need)d

cd "$H" || exit 1
if [ -s job.out ] && grep -q "ORCA TERMINATED NORMALLY" job.out; then
  echo "B_%(leg)s done before"; exit 0; fi
if ! mkdir .running 2>/dev/null; then echo "FAIL: another job holds .running in $H"; exit 1; fi

S="/scratch/18660916/d36a2_B_%(leg)s.${PBS_JOBID:-nojobid}"

cleanup() {
  kill "${WATCH:-}" 2>/dev/null
  if [ -d "$S" ]; then
    echo "scratch at exit: $(du -sh "$S" 2>/dev/null | cut -f1)"
    if cp -f "$S/job.out" "$H/job.out" 2>/dev/null; then rm -f "$H/job.out.partial"; fi
    cp -f "$S/job.property.txt" "$H/" 2>/dev/null
    cp -f "$S/job.gbw" "$H/" 2>/dev/null
    case "$S" in
      /scratch/18660916/d36a2_B_*) rm -rf "$S" ;;
      *) echo "REFUSING to remove an unexpected scratch path: $S" ;;
    esac
  fi
  rmdir "$H/.running" 2>/dev/null
}
trap cleanup EXIT

echo "host=$(hostname) start=$(date) jobid=${PBS_JOBID:-none}"
echo "cores=$(nproc) ; memory:"; free -g | head -2

# GUARD (a): enough local scratch, checked now and not seventeen hours in.
mkdir -p "$S" || { echo "FAIL: cannot create $S on $(hostname)"; exit 1; }
AVAIL=$(df -P -BG "$S" 2>/dev/null | awk 'NR==2{sub(/G$/,"",$4); print $4}')
echo "scratch=$S avail=${AVAIL:-unknown}G need=${NEED_GB}G"
case "${AVAIL:-}" in
  ''|*[!0-9]*) echo "FAIL: could not read free space on $S"; exit 1 ;;
esac
if [ "$AVAIL" -lt "$NEED_GB" ]; then
  echo "FAIL: $(hostname) has ${AVAIL}G free on /scratch and this needs ${NEED_GB}G."
  echo "      Nothing is wrong with the calculation. Resubmit and PBS will usually"
  echo "      place the job on another node."
  exit 1
fi

cp -p "$H/job.inp" "$H/substrate.xyz" "$S/" || { echo "FAIL: staging into $S"; exit 1; }
cd "$S" || exit 1

printf '%%s  B_%(leg)s attempt 3 on %%s, scratch %%s, %%sG free\\n' \\
  "$(date '+%%F %%T')" "$(hostname)" "$S" "$AVAIL" > "$P"

# GUARDS (b) and (c): the storage mode and the tail of job.out, on /home, while the
# job is still alive. Attempt 2 printed the storage line too and nobody could read it
# until the job was dead, and both legs left their last ninety minutes unrecorded.
(
  sleep 120
  while [ -d "$S" ]; do
    cp -f "$S/job.out" "$H/job.out.partial" 2>/dev/null
    {
      printf '%%s  used=%%s  storage=%%s\\n' \\
        "$(date '+%%F %%T')" \\
        "$(du -sh "$S" 2>/dev/null | cut -f1)" \\
        "$(grep -m1 'PNO-Integral Storage' "$S/job.out" 2>/dev/null | sed 's/.*\\.\\.\\.[[:space:]]*//')"
      printf '    last: %%s\\n' \\
        "$(grep -v '^[[:space:]]*$' "$S/job.out" 2>/dev/null | tail -1 | cut -c1-100)"
    } >> "$P"
    sleep 600
  done
) &
WATCH=$!

"$ORCA/orca" job.inp > job.out 2>&1 </dev/null
RC=$?
echo "orca_exit=$RC normal=$(grep -c 'ORCA TERMINATED NORMALLY' job.out) $(date)"
echo "PNO integral storage: $(grep -m1 'PNO-Integral Storage' job.out)"
echo "ccsdt energy: $(grep -m1 'E(CCSD(T))' job.out)"
grep -m1 "TOTAL RUN TIME" job.out
printf '%%s  B_%(leg)s orca_exit=%%s\\n' "$(date '+%%F %%T')" "$RC" >> "$P"
exit $RC
"""


def main():
    if len(sys.argv) != 4:
        stop("usage: cs6_d36_amend2_prepare.py <results-repo> <out-dir> <hpc-dir>")
    res, out, hpc = sys.argv[1], sys.argv[2], sys.argv[3]
    s6 = os.path.join(res, S6)
    if os.path.exists(out):
        stop("%s already exists; this script does not overwrite." % out)

    print("DIAGNOSTIC 36 AMENDMENT 2, attempt 3 of arm B")
    print("  source inputs : %s/%s" % (S6, SRC))
    print("  the two input changes, and only these two:")
    print("      %s  ->  %s" % (OLD_MAXCORE, NEW_MAXCORE))
    print("      insert before %r:" % ANCHOR)
    for k in MDCI:
        print("          %s" % k)
    print()

    for leg in LEGS:
        src_inp = os.path.join(s6, SRC, "B_%s" % leg, "job.inp")
        src_xyz = os.path.join(s6, SRC, "B_%s" % leg, "substrate.xyz")
        ref_xyz = os.path.join(s6, "d33a", "sp_42436_%s" % leg, "substrate.xyz")
        for p in (src_inp, src_xyz, ref_xyz):
            if not os.path.isfile(p):
                stop("missing %s" % p)
        src = open(src_inp).read()
        if src.splitlines()[0] != METHOD:
            stop("B_%s: the first line of the source input is not the pre-registered "
                 "method line.\n         found    %r\n         expected %r"
                 % (leg, src.splitlines()[0], METHOD))
        if sha(src_xyz) != sha(ref_xyz):
            stop("B_%s: the source geometry is not d33a's frame 42436 substrate." % leg)

        new = transform(src)
        back = untransform(new)
        if back != src:
            a, b = src.splitlines(), back.splitlines()
            stop("B_%s: the generated input does not round-trip back to the source. "
                 "%d lines in, %d out." % (leg, len(a), len(b)))

        d = os.path.join(out, "B_%s" % leg)
        os.makedirs(d)
        open(os.path.join(d, "job.inp"), "w").write(new)
        with open(src_xyz, "rb") as f, open(os.path.join(d, "substrate.xyz"), "wb") as g:
            g.write(f.read())
        if sha(os.path.join(d, "substrate.xyz")) != sha(ref_xyz):
            stop("B_%s: the copied geometry does not match d33a's." % leg)
        open(os.path.join(out, "run_B_%s.pbs" % leg), "w").write(
            PBS % {"leg": leg, "hpc": hpc, "need": NEED_GB})
        print("  B_%s written. its input, in full:" % leg)
        for l in new.splitlines():
            print("      %s" % l)
        print()

    names = []
    for leg in LEGS:
        names += ["B_%s/job.inp" % leg, "B_%s/substrate.xyz" % leg]
    names += ["run_B_%s.pbs" % leg for leg in LEGS]
    with open(os.path.join(out, "INPUT_CHECKSUMS.sha256"), "w") as f:
        for n in names:
            f.write("%s  %s\n" % (sha(os.path.join(out, n)), n))
    print("  INPUT_CHECKSUMS.sha256 over %d files" % len(names))
    print()
    print("  the two inputs are identical to each other apart from the geometry:")
    a = open(os.path.join(out, "B_R", "job.inp")).read()
    b = open(os.path.join(out, "B_TS", "job.inp")).read()
    if a != b:
        stop("the two inputs differ from each other; they must not.")
    print("      yes, byte for byte")
    print()
    print("DONE. %s" % out)


if __name__ == "__main__":
    main()
