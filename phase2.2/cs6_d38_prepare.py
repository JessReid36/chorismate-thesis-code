#!/usr/bin/env python3
"""DIAGNOSTIC 38. Five DLPNO-CCSD(T)/def2-SVPD single points on frame 42436's reactant,
differing only in the embedding point charge.

The question is whether an external point charge reaches the correlated step or only the
Hartree-Fock reference. The ORCA manual says a point-charge file "adds the point charge
contribution to the one-electron matrix and the nuclear repulsion", which is
method-agnostic by construction, but nothing in the manual states that DLPNO-CCSD(T)
honours it and nothing states that it does not. The dangerous answer is not an error: it
is a normal-looking output whose correlation energy is the bare one.

Every charge and position is READ from its design's committed site.ORCAFF.prms and
combined.xyz. Nothing is invented, and the charge is cross-checked against what the
design's own name encodes.

Refuses, and writes nothing, if:
  - the geometry is not d33a's frame 42436 reactant, byte for byte
  - a design's ORCAFF file does not have exactly 24 QM atoms plus one sphere
  - the sphere's element in combined.xyz is not the one the ORCAFF row names
  - the charge disagrees with the design's name
  - the first 24 atoms of combined.xyz are not the committed substrate
  - the five generated inputs differ in anything but the %pointcharges line

Usage: cs6_d38_prepare.py <results-repo> <out-dir> <hpc-dir>
"""
import hashlib, os, re, sys

S6 = "phase2.2/cs_stage6"
METHOD = "! DLPNO-CCSD(T) def2-SVPD AutoAux TightSCF TightPNO"
NQM = 24
NEED_GB = 100
# run name -> (design or None, charge override or None)
# The B arm runs before the A arm: the record says B1's charge makes this dianion more
# bound and A1's sign makes it metastable, so the interpretable probe goes first.
RUNS = (("N",  None,                 None),
        ("B0", "ADD-s1690_qp10000",  0.0),
        ("B1", "ADD-s1690_qp10000",  None),
        ("A0", "ADD-s678_qm10000",   0.0),
        ("A1", "ADD-s678_qm10000",   None))


def stop(msg):
    print("REFUSED  %s" % msg)
    sys.exit(1)


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 16), b""):
            h.update(b)
    return h.hexdigest()


def name_charge(design):
    """qm06500 -> -0.6500, qp10000 -> +1.0000. The design's own name, decoded."""
    m = re.search(r"_q([mp])(\d+)$", design)
    if not m:
        stop("cannot decode a charge from the design name %r" % design)
    return (1.0 if m.group(1) == "p" else -1.0) * int(m.group(2)) / 10000.0


def read_design(s6, design):
    """The sphere's charge and position, read from the committed QM/MM inputs."""
    d = os.path.join(s6, "d34", "a2_%s_R" % design)
    ff = os.path.join(d, "site.ORCAFF.prms")
    cx = os.path.join(d, "combined.xyz")
    for p in (ff, cx):
        if not os.path.isfile(p):
            stop("missing %s" % p)
    L = open(ff).read().splitlines()
    if "$atoms" not in L:
        stop("%s has no $atoms section" % ff)
    i = L.index("$atoms")
    n = int(L[i + 1].split()[0])
    if n != NQM + 1:
        stop("%s declares %d atoms; this probe needs exactly %d QM plus one sphere"
             % (ff, n, NQM))
    row = L[i + 2 + NQM].split()
    # the AMBER ORCAFF row is: index element CHARGE eps rmin eps/2 rmin
    if len(row) < 3:
        stop("%s: the sphere's row has %d fields, expected at least 3" % (ff, len(row)))
    q = float(row[2])
    want = name_charge(design)
    if abs(q - want) > 1e-9:
        stop("%s: the force field says %+.6f e and the design name encodes %+.6f e"
             % (design, q, want))
    C = open(cx).read().splitlines()
    if int(C[0].split()[0]) != NQM + 1:
        stop("%s declares %s atoms, expected %d" % (cx, C[0].split()[0], NQM + 1))
    sph = C[2 + NQM].split()
    if sph[0] != row[1]:
        stop("%s: the sphere is element %s in combined.xyz and %s in the force field"
             % (design, sph[0], row[1]))
    return dict(q=q, el=sph[0], xyz=tuple(float(x) for x in sph[1:4]),
                ff=ff, cx=cx, qm=C[2:2 + NQM])


PBS = """#!/bin/bash
#PBS -N cm_d38
#PBS -l select=1:ncpus=16:mem=120gb
#PBS -l walltime=24:00:00
#PBS -m ae
#PBS -M 18660916@sun.ac.za
#PBS -j oe
#PBS -o %(hpc)s/run.pbs.out
# DIAGNOSTIC 38: five DLPNO-CCSD(T)/def2-SVPD single points on one geometry, differing
#   only in the embedding point charge. The question is whether the charge reaches the
#   correlated step or only the Hartree-Fock reference.
# The five run sequentially in one job and each run's scratch is cleared before the next
#   starts, so the footprint is one run at a time. This is a shared cluster.
# D36's guards carry over: free space checked before the run, the PNO storage mode
#   written to /home while the job is alive, and each job.out copied back as its run
#   finishes rather than everything at the end.
# walltime 24:00:00 routes to day automatically. No -q.
export PATH="/apps/openmpi/4.1.1/bin:$PATH"
export LD_LIBRARY_PATH="/home/apps2/ORCA/6.0.1/lib:/apps/openmpi/4.1.1/lib:/apps/mambaforge/envs/medaka/lib:${LD_LIBRARY_PATH:-}"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
ORCA=/home/apps2/ORCA/6.0.1
H=%(hpc)s
P="$H/PROGRESS.txt"
NEED_GB=%(need)d
RUNS="%(runs)s"

cd "$H" || exit 1
if ! mkdir .running 2>/dev/null; then echo "FAIL: another job holds .running in $H"; exit 1; fi

S="/scratch/18660916/d38.${PBS_JOBID:-nojobid}"

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

printf '%%s  d38 on %%s, scratch %%s, %%sG free\\n' \\
  "$(date '+%%F %%T')" "$(hostname)" "$S" "$AVAIL" > "$P"

( sleep 120
  while [ -d "$S" ]; do
    printf '%%s  used=%%s  storage=%%s  running=%%s\\n' \\
      "$(date '+%%F %%T')" \\
      "$(du -sh "$S" 2>/dev/null | cut -f1)" \\
      "$(grep -hm1 'PNO-Integral Storage' "$S"/*/job.out 2>/dev/null | sed 's/.*\\.\\.\\.[[:space:]]*//' | head -1)" \\
      "$(ls "$S" 2>/dev/null | tr '\\n' ' ')" >> "$P"
    sleep 600
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
  if [ -f "$H/$r/charge.pc" ]; then cp -p "$H/$r/charge.pc" "$S/$r/"; fi
  ( cd "$S/$r" && "$ORCA/orca" job.inp > job.out 2>&1 </dev/null )
  rc=$?
  cp -f "$S/$r/job.out" "$H/$r/job.out" 2>/dev/null
  for f in job.property.txt job.gbw; do
    if [ -f "$S/$r/$f" ]; then cp -f "$S/$r/$f" "$H/$r/$f"; fi
  done
  echo "$r orca_exit=$rc normal=$(grep -c 'ORCA TERMINATED NORMALLY' "$H/$r/job.out")"
  echo "$r storage: $(grep -m1 'PNO-Integral Storage' "$H/$r/job.out")"
  echo "$r pcharges: $(grep -im1 'point charge' "$H/$r/job.out")"
  echo "$r E(0):    $(grep -m1 'E(0)' "$H/$r/job.out")"
  echo "$r E(CORR): $(grep -m1 'E(CORR)(corrected)' "$H/$r/job.out")"
  echo "$r CCSD(T): $(grep -m1 'E(CCSD(T))' "$H/$r/job.out")"
  grep -m1 "TOTAL RUN TIME" "$H/$r/job.out"
  printf '%%s  %%s orca_exit=%%s\\n' "$(date '+%%F %%T')" "$r" "$rc" >> "$P"
  # clear this run's intermediates before the next one starts, under the same guard
  case "$S/$r" in
    /scratch/18660916/d38.*/?|/scratch/18660916/d38.*/??) rm -rf "$S/$r" ;;
    *) echo "REFUSING to remove an unexpected scratch path: $S/$r" ;;
  esac
  if [ "$rc" -ne 0 ]; then RC=$rc; fi
done
echo "all runs finished, worst exit=$RC $(date)"
exit $RC
"""


def main():
    if len(sys.argv) != 4:
        stop("usage: cs6_d38_prepare.py <results-repo> <out-dir> <hpc-dir>")
    res, out, hpc = sys.argv[1], sys.argv[2], sys.argv[3]
    s6 = os.path.join(res, S6)
    if os.path.exists(out):
        stop("%s already exists; this script does not overwrite." % out)
    geo = os.path.join(s6, "d33a", "sp_42436_R", "substrate.xyz")
    if not os.path.isfile(geo):
        stop("missing %s" % geo)
    gsha = sha(geo)
    gqm = open(geo).read().splitlines()[2:2 + NQM]

    print("DIAGNOSTIC 38: does an embedding point charge reach the correlated step?")
    print("  geometry : %s/d33a/sp_42436_R/substrate.xyz" % S6)
    print("             sha256 %s" % gsha)
    print("  method   : %s" % METHOD)
    print()

    D = {}
    for design in ("ADD-s1690_qp10000", "ADD-s678_qm10000"):
        d = read_design(s6, design)
        D[design] = d
        # the QM half of combined.xyz must be the committed substrate
        for a, b in zip(d["qm"], gqm):
            fa, fb = a.split(), b.split()
            if fa[0] != fb[0] or any(abs(float(x) - float(y)) > 1e-9
                                     for x, y in zip(fa[1:4], fb[1:4])):
                stop("%s: combined.xyz's QM atoms are not the committed substrate"
                     % design)
        print("  %s" % design)
        print("     name encodes          %+.6f e" % name_charge(design))
        print("     force field says      %+.6f e      %s" % (d["q"], d["ff"]))
        print("     sphere is %s at        %.8f %.8f %.8f" % ((d["el"],) + d["xyz"]))
        print("     combined.xyz QM half matches the committed substrate, all %d atoms"
              % NQM)
        print()

    texts = {}
    for run, design, override in RUNS:
        rd = os.path.join(out, run)
        os.makedirs(rd)
        with open(geo, "rb") as f, open(os.path.join(rd, "substrate.xyz"), "wb") as g:
            g.write(f.read())
        if sha(os.path.join(rd, "substrate.xyz")) != gsha:
            stop("%s: the copied geometry does not match d33a's" % run)
        lines = [METHOD, "%maxcore 6000", "%pal nprocs 16 end", "%scf MaxIter 300 end"]
        q = None
        if design is not None:
            q = D[design]["q"] if override is None else override
            lines.append('%pointcharges "charge.pc"')
            with open(os.path.join(rd, "charge.pc"), "w") as f:
                # first line the count; then magnitude in atomic units, position in
                # Angstrom, which is the format the manual specifies
                f.write("1\n%12.6f %16.8f %16.8f %16.8f\n" % ((q,) + D[design]["xyz"]))
        lines.append("* xyzfile -2 1 substrate.xyz")
        t = "\n".join(lines) + "\n"
        open(os.path.join(rd, "job.inp"), "w").write(t)
        texts[run] = t
        print("  %-2s  %s" % (run, "no embedding charge" if design is None
                              else "%+.6f e from %s" % (q, design)))
        for l in t.splitlines():
            print("        %s" % l)
        if design is not None:
            print("        charge.pc:")
            for l in open(os.path.join(rd, "charge.pc")).read().splitlines():
                print("          %s" % l)
        print()

    # G38.7 recomputed here: the five differ only in the %pointcharges line
    base = texts["N"].splitlines()
    for run, _d, _o in RUNS:
        got = [l for l in texts[run].splitlines() if not l.startswith("%pointcharges")]
        if got != base:
            stop("%s differs from N in more than the %%pointcharges line" % run)
    print("  the five inputs are identical apart from the %pointcharges line: yes")

    open(os.path.join(out, "run.pbs"), "w").write(
        PBS % {"hpc": hpc, "need": NEED_GB,
               "runs": " ".join(r for r, _d, _o in RUNS)})

    names = []
    for run, design, _o in RUNS:
        names += ["%s/job.inp" % run, "%s/substrate.xyz" % run]
        if design is not None:
            names.append("%s/charge.pc" % run)
    names.append("run.pbs")
    with open(os.path.join(out, "INPUT_CHECKSUMS.sha256"), "w") as f:
        for n in names:
            f.write("%s  %s\n" % (sha(os.path.join(out, n)), n))
    print("  INPUT_CHECKSUMS.sha256 over %d files" % len(names))
    print()
    print("DONE. %s" % out)


if __name__ == "__main__":
    main()
