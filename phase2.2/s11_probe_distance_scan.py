#!/usr/bin/env python3
"""
s11_probe_distance_scan.py - where does the bare point-charge model break, and what
does one +1 charge at Burschowsky's distance actually buy?

WHY THIS, AND NOT A SMEARING-WIDTH CALIBRATION
A smearing width cannot be chosen before the failure is measured. Two reasons:

  1. [DTHESIS] Ch. 3.6 names the missing ingredient as REPULSION, not the singularity
     alone: Coulomb implosion arises "due to the non-convexity around the singularity
     of the charge center WITHOUT ANY OTHER REPULSIONS SUCH AS DISPERSION INTERACTIONS
     OF REAL ATOMS". A Gaussian removes the 1/r divergence but still has no Pauli
     wall, so an anionic density can still pile onto it. Smearing may treat a symptom.
  2. If the model is well behaved at the distances the design actually uses, no fix is
     needed and the added parameter would be unjustified.

So: measure first. This scan answers three questions at once - the Burschowsky
calibration, whether smearing is needed at the working distance, and what minimum
approach distance the optimiser's constraint set should enforce.

WHAT IS COMPUTED
For each calibration frame, a +1 point charge is placed on the ray from the substrate
centroid through the ether oxygen O3, at a series of distances d from O3, and single
points are run on the committed reactant and transition-state QM geometries.

The quantity of interest is the DIFFERENTIAL stabilisation

    ddE(d) = [ E_TS(field) - E_TS(bare) ] - [ E_R(field) - E_R(bare) ]

NEGATIVE ddE means the transition state is stabilised MORE than the reactant, i.e. the
barrier falls. This is the differential quantity, not the raw stabilisation of either
state: a charge that stabilises both equally does nothing catalytically, which is also
why a continuum cannot be a source of catalysis.

THREE DIAGNOSTICS, each with behaviour stated in advance so the result is falsifiable

  A. DISTANCE SCALING. Far from the substrate the charge couples to the CHANGE in
     charge distribution between reactant and TS, which is dipole-like, so
     ddE ~ d^-2 and a log-log slope of about -2 is expected. A slope that departs
     sharply from -2 as d shrinks marks the onset of the density responding
     nonlinearly to the charge, i.e. the breakdown of the bare point-charge model.
  B. LOEWDIN CHARGE on the substrate atom nearest the probe. A monotone, modest drift
     is polarisation. A sharp acceleration is density collapsing onto the probe.
  C. HOMO ENERGY of the reactant. The bare dianion is already electronically UNBOUND
     at def2-SVP (HOMO +0.082 Eh, from refstate_diagnostic via
     SOLVATION_DECISION_NOTE.md). A +1 probe should pull it DOWN, i.e. bind it. If the
     HOMO drops far below where a real cation at that distance would put it, the
     probe is over-binding the density.

THE CALIBRATION
Burschowsky et al., PNAS 2014 (pnas.1408512111) replaced Arg90 in B. subtilis
chorismate mutase with CITRULLINE - isosteric, cationic guanidinium to neutral urea -
and found the arginine worth "up to 5.9 kcal/mol to transition state stabilization but
only 0.6 kcal/mol to the binding energy of the ground state", with the citrulline NH2
3.2 A from the ether oxygen of the transition-state analogue.

THE NUMBER TO CALIBRATE AGAINST IS THE DIFFERENTIAL, 5.9 - 0.6 = 5.3 kcal/mol, NOT
5.9. Expected sign here is NEGATIVE ddE of magnitude about 5.3 at d = 3.2 A.

TWO CAVEATS THAT MUST TRAVEL WITH ANY COMPARISON
  - Burschowsky's number is an EXPERIMENTAL FREE ENERGY difference from a protein
    mutation. Ours is an ELECTRONIC energy on a bare substrate with one point charge.
    Agreement to better than a kcal/mol would be luck, not validation. Order of
    magnitude and sign are what is being tested.
  - A guanidinium is not a point charge. Its formal +1 is delocalised over CZ and two
    NH2 groups, and it carries Pauli repulsion and hydrogen bonds that a bare charge
    does not. If a point +1 at 3.2 A gives much MORE than 5.3 kcal/mol, that is itself
    evidence of over-polarisation rather than of a better catalyst.

USAGE
    python3 s11_probe_distance_scan.py generate <ensemble_dir> <workdir> [frames...]
    qsub <workdir>/s11_probe_scan.pbs
    python3 s11_probe_distance_scan.py analyse <workdir>

Default calibration frames are post-cut, full_NAC, with both a barrier and an in vacuo
reference, and barriers nearest the post-cut ensemble mean of 13.07 kcal/mol.
"""
import sys
import math
import re
from pathlib import Path

HARTREE = 627.5095
CHARGE, MULT = -2, 1
O3_INDEX = 7                      # ether oxygen in the cha_gaff atom order
PROBE_Q = +1.0
# Distances from O3, in Angstrom. 3.2 is Burschowsky's citrulline NH2 separation.
DISTANCES = [1.5, 2.0, 2.5, 2.8, 3.2, 3.6, 4.0, 5.0, 6.0, 8.0, 10.0]
BURSCHOWSKY_D = 3.2
BURSCHOWSKY_TS = 5.9              # kcal/mol, TS stabilisation
BURSCHOWSKY_GS = 0.6              # kcal/mol, ground-state binding
BURSCHOWSKY_DIFFERENTIAL = BURSCHOWSKY_TS - BURSCHOWSKY_GS
DEFAULT_FRAMES = ["24883", "20000", "43087"]
LEVEL = "! B3LYP D3BJ def2-SVP def2/J RIJCOSX TightSCF"
# Every job here is at the production level and is compute-bound: measured on the
# equivalent b97-3c jobs, 202 s serial against 71-89 s on 8 ranks, with 81 per cent of
# the time in SCF iterations and 12 per cent in startup. That is the profile where
# parallelism pays. Cheap methods with second-scale runtimes are startup-dominated and
# must stay serial - see METHODS in s10_cheap_method_benchmark.py.
NPROCS = 8


def read_xyz(p):
    L = Path(p).read_text().splitlines()
    n = int(L[0].split()[0])
    return [(f[0], float(f[1]), float(f[2]), float(f[3]))
            for f in (l.split() for l in L[2:2 + n])]


def probe_at(atoms, d):
    cx = sum(a[1] for a in atoms) / len(atoms)
    cy = sum(a[2] for a in atoms) / len(atoms)
    cz = sum(a[3] for a in atoms) / len(atoms)
    ox, oy, oz = atoms[O3_INDEX][1:4]
    vx, vy, vz = ox - cx, oy - cy, oz - cz
    n = math.sqrt(vx * vx + vy * vy + vz * vz)
    return (ox + d * vx / n, oy + d * vy / n, oz + d * vz / n)


def _dist(a, b):
    """Euclidean distance. math.dist needs Python 3.8; the cluster login node is older."""
    return math.sqrt((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2)


def nearest_atom(atoms, pt):
    best = min(range(len(atoms)), key=lambda i: _dist(pt, atoms[i][1:4]))
    return best, _dist(pt, atoms[best][1:4])


def generate(ens, work, frames):
    ens, work = Path(ens), Path(work)
    work.mkdir(parents=True, exist_ok=True)
    n_inp = 0
    manifest = ["# frame\tstate\tdistance_A\tprobe_x\tprobe_y\tprobe_z\t"
                "nearest_atom\tnearest_dist_A"]
    for fr in frames:
        d = ens / f"frame_{fr}"
        for state, src in (("R", "reactant_qm.xyz"), ("TS", "transition_state_qm.xyz")):
            sp = d / src
            if not sp.exists():
                print(f"  SKIP {fr} {state}: {src} not found")
                continue
            atoms = read_xyz(sp)
            if len(atoms) != 24:
                print(f"  SKIP {fr} {state}: {len(atoms)} atoms, expected 24")
                continue
            gx = work / f"{fr}_{state}.xyz"
            gx.write_text(f"{len(atoms)}\n{fr} {state}\n" +
                          "".join(f"{a[0]:<3}{a[1]:>15.8f}{a[2]:>15.8f}{a[3]:>15.8f}\n"
                                  for a in atoms))
            # bare reference, no probe
            (work / f"{fr}_{state}_bare.inp").write_text(
                f"{LEVEL}\n%maxcore 3000\n%pal nprocs 8 end\n%scf MaxIter 300 end\n"
                f"%output Print[P_Loewdin] 1 end\n"
                f"* xyzfile {CHARGE} {MULT} {fr}_{state}.xyz\n")
            n_inp += 1
            for dd in DISTANCES:
                px, py, pz = probe_at(atoms, dd)
                ia, ndist = nearest_atom(atoms, (px, py, pz))
                manifest.append(f"{fr}\t{state}\t{dd}\t{px:.6f}\t{py:.6f}\t{pz:.6f}\t"
                                f"{atoms[ia][0]}{ia}\t{ndist:.4f}")
                tag = f"{fr}_{state}_d{dd:g}"
                (work / f"{tag}.pc").write_text(
                    f"1\n{PROBE_Q:.4f} {px:.8f} {py:.8f} {pz:.8f}\n")
                (work / f"{tag}.inp").write_text(
                    f"{LEVEL}\n%maxcore 3000\n%pal nprocs 8 end\n%scf MaxIter 300 end\n"
                    f"%output Print[P_Loewdin] 1 end\n"
                    f'%pointcharges "{tag}.pc"\n'
                    f"* xyzfile {CHARGE} {MULT} {fr}_{state}.xyz\n")
                n_inp += 1
    (work / "probe_manifest.tsv").write_text("\n".join(manifest) + "\n")

    pbs = work / "s11_probe_scan.pbs"
    pbs.write_text(f"""#!/usr/bin/env bash
#PBS -N cm_s11probe
#PBS -l select=1:ncpus=8:mem=24gb
#PBS -l walltime=24:00:00
#PBS -m ae
#PBS -M 18660916@sun.ac.za
#PBS -j oe
#PBS -o {work.resolve()}/s11_probe_scan.pbs.out
#
# s11_probe_scan.pbs - probe-distance scan.
# ORCA is NOT on PATH on this cluster and cannot run on a login node.
# Conventions copied from s8_invacuo_new.pbs.

set -uo pipefail
ORCA=/home/apps2/ORCA/6.0.1
WORK={work.resolve()}

export PATH="/apps/openmpi/4.1.1/bin:$PATH"
export LD_LIBRARY_PATH="$ORCA/lib:/apps/openmpi/4.1.1/lib:/apps/mambaforge/envs/medaka/lib:${{LD_LIBRARY_PATH:-}}"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1

cd "$WORK" || exit 1
if ! mkdir .running 2>/dev/null; then echo "FAIL: another job holds .running"; exit 1; fi
echo "host=$(hostname) start=$(date)"

for f in *.out; do
  [ -e "$f" ] || continue
  grep -q "TERMINATED NORMALLY" "$f" 2>/dev/null || rm -f "$f"
done
""" + r"""
n=0; ok=0; bad=0
SCRATCH="${TMPDIR:-/tmp}/orca_$$"
mkdir -p "$SCRATCH" || { echo "FAIL: cannot create $SCRATCH"; exit 1; }
trap 'rm -rf "$SCRATCH"; rmdir .running 2>/dev/null' EXIT
echo "scratch: $SCRATCH"

# MEASURED 2026-09-30: the same 24-atom job takes 28.7 s real on the shared
# filesystem and 6.6 s in node-local scratch, user time 4.9 s vs 3.6 s. The
# difference is ORCA temporary-file I/O, not compute. Each job therefore runs in
# node-local scratch and only the .out is copied back.
for inp in *.inp; do
  out="${inp%.inp}.out"
  if [ -f "$out" ] && grep -q "TERMINATED NORMALLY" "$out"; then continue; fi
  n=$((n+1))
  base="${inp%.inp}"
  rm -rf "$SCRATCH/$base"; mkdir -p "$SCRATCH/$base"
  cp "$inp" "$SCRATCH/$base/" || continue
  for dep in $(grep -oE '[A-Za-z0-9_.-]+\.(xyz|pc)' "$inp" | sort -u); do
    [ -f "$dep" ] && cp "$dep" "$SCRATCH/$base/"
  done
  ( cd "$SCRATCH/$base" && "$ORCA/orca" "$inp" > "$out" 2>&1 </dev/null )
  cp "$SCRATCH/$base/$out" "$out" 2>/dev/null
  rm -rf "$SCRATCH/$base"
  if grep -q "TERMINATED NORMALLY" "$out" 2>/dev/null; then ok=$((ok+1))
  else bad=$((bad+1)); echo "  WARN did not terminate: $inp"; fi
  if [ $((n % 50)) -eq 0 ]; then echo "  ... $n run, $ok ok, $bad failed  $(date +%H:%M)"; fi
done
echo
echo "run=$n ok=$ok failed=$bad"
echo "end=$(date)"
""")
    pbs.chmod(0o755)
    print(f"wrote {n_inp} inputs for {len(frames)} frames")
    print(f"  distances from O3: {DISTANCES}")
    print(f"  Burschowsky distance {BURSCHOWSKY_D} A is in the set")
    print(f"  probe {PROBE_Q:+.1f} e on the centroid->O3 ray")
    print(f"\n  ORCA is not on PATH here; submit with:  qsub {pbs.resolve()}")


def energy(p):
    t = Path(p).read_text(errors="replace")
    if "ORCA TERMINATED NORMALLY" not in t:
        return None, None, None
    e = re.findall(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", t)
    homo = None
    m = re.findall(r"^\s*\d+\s+2\.0000\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s*$", t, re.M)
    if m:
        homo = float(m[-1][0])          # last doubly occupied orbital, Eh
    low = {}
    blk = re.search(r"LOEWDIN ATOMIC CHARGES\n-+\n(.*?)\n\n", t, re.S)
    if blk:
        for line in blk.group(1).splitlines():
            f = line.split()
            if len(f) >= 4:
                try:
                    low[int(f[0])] = float(f[-1])
                except ValueError:
                    pass
    return (float(e[-1]) if e else None), homo, low


def analyse(work):
    work = Path(work)
    man = {}
    for l in (work / "probe_manifest.tsv").read_text().splitlines()[1:]:
        f = l.split("\t")
        man[(f[0], f[1], float(f[2]))] = (f[6], float(f[7]))
    data = {}
    for out in work.glob("*.out"):
        m = re.match(r"(\d+)_(R|TS)_(bare|d[\d.]+)\.out$", out.name)
        if not m:
            continue
        key = (m.group(1), m.group(2),
               None if m.group(3) == "bare" else float(m.group(3)[1:]))
        data[key] = energy(out)
    frames = sorted({k[0] for k in data})
    print(f"probe-distance scan, {len(frames)} frame(s)\n")
    print("ddE = [E_TS(field) - E_TS(bare)] - [E_R(field) - E_R(bare)]")
    print("NEGATIVE ddE means the transition state is stabilised more than the"
          " reactant.\n")
    for fr in frames:
        eb_r = data.get((fr, "R", None), (None,) * 3)
        eb_t = data.get((fr, "TS", None), (None,) * 3)
        if eb_r[0] is None or eb_t[0] is None:
            print(f"  frame {fr}: bare reference missing, skipping")
            continue
        print(f"  frame {fr}   bare HOMO(R) {eb_r[1]:+.6f} Eh "
              f"({'UNBOUND' if eb_r[1] and eb_r[1] > 0 else 'bound'})")
        print(f"    {'d/A':>6}{'ddE':>10}{'nearest':>10}{'gap/A':>8}"
              f"{'HOMO(R)':>11}{'dq(near)':>10}{'slope':>8}")
        rows = []
        for d in DISTANCES:
            fr_r = data.get((fr, "R", d), (None,) * 3)
            fr_t = data.get((fr, "TS", d), (None,) * 3)
            if fr_r[0] is None or fr_t[0] is None:
                continue
            dd = ((fr_t[0] - eb_t[0]) - (fr_r[0] - eb_r[0])) * HARTREE
            nm, gap = man.get((fr, "R", d), ("?", float("nan")))
            idx = int(re.sub(r"^\D+", "", nm)) if re.search(r"\d", nm) else None
            dq = (fr_r[2].get(idx, float("nan")) - eb_r[2].get(idx, float("nan"))
                  if idx is not None and fr_r[2] and eb_r[2] else float("nan"))
            rows.append((d, dd, nm, gap, fr_r[1], dq))
        for i, (d, dd, nm, gap, homo, dq) in enumerate(rows):
            slope = ""
            if i > 0:
                d0, dd0 = rows[i - 1][0], rows[i - 1][1]
                if dd0 * dd > 0 and d0 > 0:
                    slope = f"{(math.log(abs(dd)) - math.log(abs(dd0))) / (math.log(d) - math.log(d0)):+.2f}"
            star = "  <== Burschowsky" if abs(d - BURSCHOWSKY_D) < 1e-6 else ""
            print(f"    {d:>6.1f}{dd:>10.2f}{nm:>10}{gap:>8.2f}"
                  f"{homo if homo is None else f'{homo:>+11.5f}'}"
                  f"{dq:>10.3f}{slope:>8}{star}")
        at32 = [r for r in rows if abs(r[0] - BURSCHOWSKY_D) < 1e-6]
        if at32:
            v = at32[0][1]
            print(f"\n    at {BURSCHOWSKY_D} A: ddE = {v:+.2f} kcal/mol")
            print(f"    Burschowsky differential (5.9 TS - 0.6 GS) = "
                  f"-{BURSCHOWSKY_DIFFERENTIAL:.1f} kcal/mol expected sign negative")
            if v > 0:
                print("    SIGN WRONG - the probe raises the barrier. Check the ray "
                      "direction and the O3 index before interpreting anything else.")
            elif abs(v) > 2 * BURSCHOWSKY_DIFFERENTIAL:
                print("    MUCH LARGER than the measured differential. A bare point "
                      "charge has no Pauli repulsion and no delocalisation, so this is "
                      "evidence of OVER-POLARISATION rather than of a better catalyst.")
            else:
                print("    Same sign and order of magnitude as the measured value.")
        far = [r for r in rows if r[0] >= 5.0]
        if len(far) >= 2:
            sl = [(math.log(abs(far[i][1])) - math.log(abs(far[i-1][1]))) /
                  (math.log(far[i][0]) - math.log(far[i-1][0]))
                  for i in range(1, len(far)) if far[i][1] * far[i-1][1] > 0]
            if sl:
                print(f"    far-field log-log slope {sum(sl)/len(sl):+.2f} "
                      f"(dipole-like coupling predicts about -2)")
        print()
    print("READING THE SCAN")
    print("  The largest d at which the slope is near -2, the Loewdin drift is smooth")
    print("  and the HOMO is merely pulled down is the shortest distance a BARE point")
    print("  charge may be trusted. That distance is the floor for the optimiser's")
    print("  minimum approach constraint.")
    print("  If 3.2 A sits comfortably inside the well-behaved region, smearing is not")
    print("  needed for this design and adding it would be an unjustified parameter.")
    print("  If 3.2 A is already past the breakdown, a bare point charge cannot")
    print("  represent Arg90 and the charge model must change - but note that")
    print("  [DTHESIS] Ch.3.6 attributes the failure to missing REPULSION, so a")
    print("  Gaussian alone may not be sufficient.")


if __name__ == "__main__":
    if len(sys.argv) < 3 or sys.argv[1] not in ("generate", "analyse"):
        sys.exit(__doc__)
    if sys.argv[1] == "generate":
        fr = sys.argv[4:] if len(sys.argv) > 4 else DEFAULT_FRAMES
        generate(sys.argv[2], sys.argv[3], fr)
    else:
        analyse(sys.argv[2])
