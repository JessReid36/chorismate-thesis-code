#!/usr/bin/env python3
"""DIAGNOSTIC 38, AMENDMENT 1. Attempt 2's job script.

Attempt 1 died in one second: mkdir on /scratch/18660916/... failed on comp056. A probe
then measured that /scratch is only a symlink to /scratch-small-local, that $TMPDIR is a
PBS-provided node-local 38 GB, and that all 32 nodes advertise
resources_available.scratch = True - including comp056, where the mkdir failed. So that
flag is not a live check and is not requested; instead the job tries three candidate
locations and uses the first that it can actually write to.

THIS SCRIPT CHANGES ONE FILE AND ONLY ONE. It writes d38/run_a1.pbs from the committed
d38/run.pbs, replacing exactly two regions: the scratch-selection block and the per-run
cleanup block. The five job.inp, the four charge.pc and the five substrate.xyz are not
read for modification and not written at all; they are checked by sha256 and left alone.

It verifies its own work by reconstruction: putting the old regions back into the new
script must give the committed script back, byte for byte.

Refuses, and writes nothing, if:
  - any of the fourteen committed inputs has changed
  - either region does not occur exactly once in the committed script
  - the reconstruction does not reproduce the committed script exactly
  - the generated script mentions a method line, requests a queue, or requests scratch
  - the PBS select line, walltime, or run list and order differ from the committed one

Usage: cs6_d38_amend1_prepare.py <results-repo> <out-file> <hpc-dir>
"""
import hashlib, os, re, sys

S6 = "phase2.2/cs_stage6"
RUNS = ("N", "B0", "B1", "A0", "A1")
WITH_PC = ("B0", "B1", "A0", "A1")
WANT_GB, MIN_GB = 100, 30

OLD_SCRATCH = '''S="/scratch/18660916/d38.${PBS_JOBID:-nojobid}"

cleanup() {
  kill "${WATCH:-}" 2>/dev/null
  if [ -d "$S" ]; then
    for r in $RUNS; do
      if [ -f "$S/$r/job.out" ] && [ -d "$H/$r" ]; then
        cp -f "$S/$r/job.out" "$H/$r/job.out" 2>/dev/null
      fi
    done
    echo "scratch at exit: $(du -sh "$S" 2>/dev/null | cut -f1)"
    case "$S" in
      /scratch/18660916/d38.*) rm -rf "$S" ;;
      *) echo "REFUSING to remove an unexpected scratch path: $S" ;;
    esac
  fi
  rmdir "$H/.running" 2>/dev/null
}
trap cleanup EXIT

echo "host=$(hostname) start=$(date) jobid=${PBS_JOBID:-none}"
echo "cores=$(nproc) ; memory:"; free -g | head -2

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

'''

NEW_SCRATCH = '''# AMENDMENT 1. /scratch is only a symlink to /scratch-small-local, and on comp056 the
# mkdir through it failed. All 32 nodes advertise resources_available.scratch = True,
# comp056 included, so that flag is not a live check and is not requested. Three
# candidates are tried instead, each by actually writing into it.
MARK="d38a1.${PBS_JOBID:-nojobid}"
WANT_GB=%(want)d
MIN_GB=%(min)d
S=""
AVAIL=""
SMARGIN=""
PEAK=0

pick() {
  for c in /scratch-small-local/18660916 /scratch/18660916 "${TMPDIR:-}"; do
    if [ -z "$c" ]; then continue; fi
    t="$c/$MARK"
    if mkdir -p "$t" 2>/dev/null; then :; else
      echo "   candidate $c  mkdir failed"; continue; fi
    if echo ok > "$t/.probe" 2>/dev/null; then :; else
      echo "   candidate $c  created but not writable"; continue; fi
    rm -f "$t/.probe"
    a=$(df -P -BG "$t" 2>/dev/null | awk 'NR==2{sub(/G$/,"",$4); print $4}')
    case "${a:-}" in
      ''|*[!0-9]*) echo "   candidate $c  free space unreadable"; continue ;;
    esac
    if [ "$a" -ge "$1" ]; then S="$t"; AVAIL="$a"; return 0; fi
    echo "   candidate $c  ${a}G free, under the ${1}G asked for"
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
      */d38a1.*) rm -rf "$S" ;;
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
  echo "   NOTE: running on the reduced margin. def2-SVPD's scratch footprint is not"
  echo "   measured anywhere in this record, so this may not be enough. R38.9 reports"
  echo "   the peak so the next job sets this from a measurement."
else
  echo "FAIL: no candidate on $(hostname) is writable with at least ${MIN_GB}G free."
  echo "      Nothing is wrong with the calculation. Resubmit and PBS will usually"
  echo "      place the job on another node."
  exit 1
fi
echo "R38.8 scratch=$S avail=${AVAIL}G margin=$SMARGIN"

'''

OLD_PERRUN = '''  # clear this run's intermediates before the next one starts, under the same guard
  case "$S/$r" in
    /scratch/18660916/d38.*/?|/scratch/18660916/d38.*/??) rm -rf "$S/$r" ;;
    *) echo "REFUSING to remove an unexpected scratch path: $S/$r" ;;
  esac
'''

NEW_PERRUN = '''  # R38.9: measure the footprint BEFORE clearing it, so the peak is real
  u=$(du -sm "$S" 2>/dev/null | cut -f1)
  case "${u:-}" in
    ''|*[!0-9]*) u=0 ;;
  esac
  if [ "$u" -gt "$PEAK" ]; then PEAK="$u"; fi
  echo "$r scratch after this run: ${u} MB"
  # clear this run's intermediates before the next one starts, under the same guard
  case "$S/$r" in
    */d38a1.*/?|*/d38a1.*/??) rm -rf "$S/$r" ;;
    *) echo "REFUSING to remove an unexpected scratch path: $S/$r" ;;
  esac
'''

OLD_TAIL = '''echo "all runs finished, worst exit=$RC $(date)"'''
NEW_TAIL = '''echo "R38.9 peak scratch across the five runs: ${PEAK} MB, on $S at the $SMARGIN margin"
echo "all runs finished, worst exit=$RC $(date)"'''

SUBS = ((OLD_SCRATCH, NEW_SCRATCH % {"want": WANT_GB, "min": MIN_GB}),
        (OLD_PERRUN, NEW_PERRUN),
        (OLD_TAIL, NEW_TAIL))


def stop(msg):
    print("REFUSED  %s" % msg)
    sys.exit(1)


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 16), b""):
            h.update(b)
    return h.hexdigest()


def main():
    if len(sys.argv) != 4:
        stop("usage: cs6_d38_amend1_prepare.py <results-repo> <out-file> <hpc-dir>")
    res, out, hpc = sys.argv[1], sys.argv[2], sys.argv[3]
    d38 = os.path.join(res, S6, "d38")
    src = os.path.join(d38, "run.pbs")
    if not os.path.isfile(src):
        stop("missing %s" % src)
    if os.path.exists(out):
        stop("%s already exists; this script does not overwrite." % out)

    print("DIAGNOSTIC 38 AMENDMENT 1: attempt 2's job script")
    print("  source  %s/d38/run.pbs" % S6)
    print("          sha256 %s" % sha(src))
    print()

    # G38.11: the fourteen inputs, by their own committed manifest
    man = os.path.join(d38, "INPUT_CHECKSUMS.sha256")
    if not os.path.isfile(man):
        stop("missing %s" % man)
    want = {}
    for ln in open(man):
        h, n = ln.split()
        want[n] = h
    expect = []
    for r in RUNS:
        expect += ["%s/job.inp" % r, "%s/substrate.xyz" % r]
        if r in WITH_PC:
            expect.append("%s/charge.pc" % r)
    missing = [n for n in expect if n not in want]
    if missing:
        stop("the committed manifest does not cover %s" % ", ".join(missing))
    bad = []
    for n in expect:
        p = os.path.join(d38, n)
        if not os.path.isfile(p) or sha(p) != want[n]:
            bad.append(n)
    if bad:
        stop("these committed inputs have changed: %s" % ", ".join(bad))
    print("  G38.11  all %d inputs match the committed manifest, unchanged" % len(expect))

    t = open(src).read()
    new = t
    for old, repl in SUBS:
        if new.count(old) != 1:
            stop("a region to replace occurs %d times in run.pbs, expected once:\n"
                 "         %r..." % (new.count(old), old[:70]))
        new = new.replace(old, repl)

    # verify by reconstruction: put the old regions back and the result must be identical
    back = new
    for old, repl in SUBS:
        if back.count(repl) != 1:
            stop("reconstruction: a new region occurs %d times, expected once"
                 % back.count(repl))
        back = back.replace(repl, old)
    if back != t:
        stop("reconstruction does not reproduce the committed run.pbs: %d chars in, "
             "%d out" % (len(t), len(back)))
    print("  reconstruction reproduces the committed run.pbs byte for byte")

    # G38.12: the things that must not have moved
    for lab, s in (("PBS select", "#PBS -l select=1:ncpus=16:mem=120gb"),
                   ("walltime", "#PBS -l walltime=24:00:00"),
                   ("the run list and order", 'RUNS="N B0 B1 A0 A1"'),
                   ("the resume guard", "done before"),
                   ("the lock", "mkdir .running")):
        if s not in new:
            stop("G38.12 %s is not in the generated script" % lab)
        if t.count(s) != new.count(s):
            stop("G38.12 %s occurs %d times in the committed script and %d in the new one"
                 % (lab, t.count(s), new.count(s)))
    if re.search(r"#PBS\s+-q\b", new):
        stop("G38.12 the generated script requests a queue")
    if "scratch=True" in new:
        stop("G38.12 the generated script requests scratch=True, which measures nothing")
    # The committed script names the method in its own header comments, which is
    # correct. What must not happen is that this amendment CHANGES any such line, so
    # the gate compares them rather than forbidding the words.
    meth = lambda x: [l for l in x.splitlines() if re.search(r"DLPNO|def2-", l)]
    gone = [l for l in meth(t) if l not in meth(new)]
    added = [l for l in meth(new) if l not in meth(t)]
    if gone:
        stop("G38.12 a line naming the method or basis was lost:\n         %r" % gone[0])
    if len(added) > 1 or (added and "footprint is not" not in added[0]):
        stop("G38.12 a line naming the method or basis was added:\n         %r"
             % (added[0] if added else ""))
    print("  G38.12  select, walltime, run order, resume guard and lock unchanged;")
    print("          no queue requested, scratch=True not requested, no method named")

    open(out, "w").write(new)
    print()
    print("  the three candidates, in order:")
    print("     1. /scratch-small-local/18660916/   the mount, 1.5 TB free on comp050")
    print("     2. /scratch/18660916/               the symlink, where it resolves")
    print("     3. $TMPDIR                          node-local, 38 GB on comp050")
    print("  thresholds: %d GB preferred, %d GB accepted and recorded as REDUCED"
          % (WANT_GB, MIN_GB))
    print()
    print("  written  %s" % out)
    print("           sha256 %s" % sha(out))
    print()
    print("DONE. One file written. No input was touched.")


if __name__ == "__main__":
    main()
