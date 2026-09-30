#!/usr/bin/env python3
"""
s10_cheap_method_benchmark.py - can a cheap method screen charge designs for THIS reaction?

WHY THIS EXISTS
The Phase 2 design loop needs many barrier evaluations. B3LYP-D3BJ/def2-SVP is too
expensive to put inside an optimiser. Behrens & Hartke (Top. Catal. 2021) run GFN-xTB
inside their GOCAT loop and justify it as: "As long as the relative energies are accurate
enough between different solution candidates, our global optimization scheme will be able
to find optima that will be transferable to higher quality ab initio calculations."

That warrant is about RELATIVE energies, so it has to be tested on relative energies for
OUR system, not assumed. The literature is specifically unfavourable here. Claeyssens et
al. 2011 p.10: "both AM1 and HF overestimate the reaction barrier severely, while standard
non-hybrid density functional techniques (e.g. PBE and SCC-DFTB) underestimate the
reaction barrier". p.11: AM1 "overestimates substrate compression in chorismate mutase" -
a geometry error, not only an energy one. Woodcock's SCC-DFTB barrier is 6.0 kcal/mol
against Claeyssens' 11.3 at B3LYP. GFN-xTB is in the same tight-binding family.

So: an absolute offset is tolerable, a ranking failure is not.

WHAT IT TESTS

  Test A, the substrate PES. For each frame, a bare single point on the committed
  reactant and transition-state QM geometries. This reproduces exactly the quantity
  already in ensemble_barriers.tsv as barrier_vac, which is a bare 24-atom gas-phase
  barrier at the QM/MM geometries. Cheap method against DFT, frame by frame, across the
  ensemble. What matters is the CORRELATION, not the offset.

  Test B, the field response. The same two geometries with a single +1 point charge
  placed 3.2 A from the ether oxygen O3, outward along the centroid->O3 vector. The
  distance and the target are taken from Burschowsky et al., PNAS 2014, who replaced
  Arg90 with citrulline and found the cationic guanidinium worth 5.9 kcal/mol of TS
  stabilisation with its NH2 3.2 A from the ether oxygen of the TS analogue.
  A screening method must reproduce the RESPONSE to a field, not just the bare barrier.
  A method can rank bare barriers well and still get the response wrong, and the response
  is the entire design signal.

  VERDICT RULE, fixed before running so the result cannot be rationalised afterwards:
    - Test A Pearson r >= 0.95 across frames  -> usable as a screen
    - Test B response within 25 per cent of the DFT response, and same sign, every frame
    - An absolute offset of any size is acceptable provided both hold.
  A method failing A is unusable. A method passing A but failing B is usable only for
  pre-filtering geometry, never for scoring a design.

CANDIDATES
  XTB2        GFN2-xTB, what Behrens use
  HF-3c       already proven to run on this system, step 13 smoke test
  B97-3c      composite GGA, modern, no Claeyssens-era verdict applies
  r2SCAN-3c   composite meta-GGA, currently the best accuracy-per-cost composite
  reference   B3LYP D3BJ def2-SVP def2/J RIJCOSX, the production level

CAVEATS TO CARRY WHATEVER WINS
  - The substrate is a DIANION. Minimal-basis and tight-binding methods handle diffuse
    anionic density poorly. The in vacuo reference is already electronically unbound at
    def2-SVP (HOMO +0.082 Eh), so the cheap methods are being asked to describe an
    unbound species.
  - AM1's documented failure for this reaction is GEOMETRIC. This benchmark only
    EVALUATES on fixed DFT geometries. If a cheap method is ever allowed to optimise
    rather than evaluate, that failure mode returns and this benchmark says nothing
    about it.

USAGE
  Generate inputs and the PBS script:
      python3 s10_cheap_method_benchmark.py generate <ensemble_dir> <workdir>
      qsub <workdir>/s10_benchmark.pbs
  ORCA is not on PATH on this cluster and must be invoked by full path with specific
  library settings; the benchmark cannot run on a login node.
  After the jobs finish, analyse:
      python3 s10_cheap_method_benchmark.py analyse <workdir> <ensemble_barriers.tsv>

  <ensemble_dir> is 05_qmmm/19_ensemble_barriers
"""
import sys
import math
import re
from pathlib import Path

HARTREE = 627.5095
CHARGE, MULT = -2, 1
# Atom order in *_qm.xyz is the cha_gaff order. O3 is the ether oxygen, index 7 (0-based),
# and is the atom typed os; C4 at index 8 is the c3 it bonds to. The breaking bond is
# C4-O3 per the step 11 atom mapping.
O3_INDEX = 7
BURSCHOWSKY_R = 3.2      # Angstrom, citrulline NH2 to ether oxygen, PNAS 2014
PROBE_Q = +1.0

METHODS = {
    "xtb2":      "! XTB2",
    "hf3c":      "! HF-3c",
    "b97-3c":    "! B97-3c",
    "r2scan-3c": "! r2SCAN-3c",
    "ref":       "! B3LYP D3BJ def2-SVP def2/J RIJCOSX TightSCF",
}


def read_xyz(p):
    lines = Path(p).read_text().splitlines()
    n = int(lines[0].split()[0])
    at = []
    for l in lines[2:2 + n]:
        f = l.split()
        at.append((f[0], float(f[1]), float(f[2]), float(f[3])))
    return at


def probe_position(atoms):
    """+1 charge at BURSCHOWSKY_R from O3, outward along centroid->O3."""
    cx = sum(a[1] for a in atoms) / len(atoms)
    cy = sum(a[2] for a in atoms) / len(atoms)
    cz = sum(a[3] for a in atoms) / len(atoms)
    ox, oy, oz = atoms[O3_INDEX][1], atoms[O3_INDEX][2], atoms[O3_INDEX][3]
    vx, vy, vz = ox - cx, oy - cy, oz - cz
    n = math.sqrt(vx * vx + vy * vy + vz * vz)
    if n < 1e-6:
        raise SystemExit("O3 sits on the centroid; cannot define an outward direction")
    return (ox + BURSCHOWSKY_R * vx / n,
            oy + BURSCHOWSKY_R * vy / n,
            oz + BURSCHOWSKY_R * vz / n)


def generate(ens_dir, work):
    ens_dir, work = Path(ens_dir), Path(work)
    work.mkdir(parents=True, exist_ok=True)
    frames = sorted(d for d in ens_dir.glob("frame_*")
                    if (d / "reactant_qm.xyz").exists()
                    and (d / "transition_state_qm.xyz").exists())
    print(f"{len(frames)} frames with both a reactant and a transition-state geometry")
    jobs = []
    for d in frames:
        fr = d.name.replace("frame_", "")
        for state, src in (("R", "reactant_qm.xyz"), ("TS", "transition_state_qm.xyz")):
            atoms = read_xyz(d / src)
            if len(atoms) != 24:
                print(f"  SKIP {fr} {state}: {len(atoms)} atoms, expected 24")
                continue
            xyz = work / f"{fr}_{state}.xyz"
            xyz.write_text(f"{len(atoms)}\n{fr} {state}\n" +
                           "".join(f"{a[0]:<3}{a[1]:>15.8f}{a[2]:>15.8f}{a[3]:>15.8f}\n"
                                   for a in atoms))
            px, py, pz = probe_position(atoms)
            pc = work / f"{fr}_{state}.pc"
            pc.write_text(f"1\n{PROBE_Q:.4f} {px:.8f} {py:.8f} {pz:.8f}\n")
            for mname, kw in METHODS.items():
                for field in ("bare", "field"):
                    tag = f"{fr}_{state}_{mname}_{field}"
                    # NO %pal. These are 24-atom single points: ORCA reports 2-4 s of
                    # compute, but with %pal nprocs 8 each job took ~180 s wall, i.e.
                    # 98% MPI startup and teardown. Measured from file timestamps: the
                    # xtb2 jobs, the only ones without %pal, ran 8 s apart while every
                    # other method ran 2.5-4.5 min apart. Serial is ~40x faster here.
                    body = [kw, "%maxcore 3000", "%scf MaxIter 300 end"]
                    if field == "field":
                        body.append(f'%pointcharges "{fr}_{state}.pc"')
                    body.append(f"* xyzfile {CHARGE} {MULT} {fr}_{state}.xyz")
                    (work / f"{tag}.inp").write_text("\n".join(body) + "\n")
                    jobs.append(tag)
    pbs = work / "s10_benchmark.pbs"
    pbs.write_text(f"""#!/usr/bin/env bash
#PBS -N cm_s10bm
#PBS -l select=1:ncpus=1:mem=8gb
#PBS -l walltime=24:00:00
#PBS -m ae
#PBS -M 18660916@sun.ac.za
#PBS -j oe
#PBS -o {work.resolve()}/s10_benchmark.pbs.out
#
# s10_benchmark.pbs - cheap-method screening benchmark.
# Conventions copied from s8_invacuo_new.pbs: ORCA is NOT on PATH and must be
# invoked by full path with these library settings. This cannot run on a login node.

set -uo pipefail
ORCA=/home/apps2/ORCA/6.0.1
WORK={work.resolve()}

export PATH="/apps/openmpi/4.1.1/bin:$PATH"
export LD_LIBRARY_PATH="$ORCA/lib:/apps/openmpi/4.1.1/lib:/apps/mambaforge/envs/medaka/lib:${{LD_LIBRARY_PATH:-}}"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1

cd "$WORK" || exit 1
if ! mkdir .running 2>/dev/null; then echo "FAIL: another job holds .running"; exit 1; fi
trap 'rmdir .running 2>/dev/null' EXIT
echo "host=$(hostname) start=$(date)"

# Clear any output left by a failed attempt, so they are not mistaken for results.
for f in *.out; do
  [ -e "$f" ] || continue
  grep -q "TERMINATED NORMALLY" "$f" 2>/dev/null || rm -f "$f"
done

n=0; ok=0; bad=0
for inp in *.inp; do
  out="${{inp%.inp}}.out"
  if [ -f "$out" ] && grep -q "TERMINATED NORMALLY" "$out"; then continue; fi
  n=$((n+1))
  "$ORCA/orca" "$inp" > "$out" 2>&1 </dev/null
  if grep -q "TERMINATED NORMALLY" "$out"; then
    ok=$((ok+1))
  else
    bad=$((bad+1)); echo "  WARN did not terminate: $inp"
  fi
  if [ $((n % 50)) -eq 0 ]; then echo "  ... $n run, $ok ok, $bad failed  $(date +%H:%M)"; fi
done
echo
echo "run=$n ok=$ok failed=$bad"
echo "end=$(date)"
echo "now run, from ~/system_development:"
echo "  python3 s10_cheap_method_benchmark.py analyse {work.resolve()} \\"
echo "      05_qmmm/19_ensemble_barriers/ensemble_barriers.tsv"
""")
    pbs.chmod(0o755)
    print(f"wrote {len(jobs)} inputs and {pbs}")
    print(f"  {len(frames)} frames x 2 states x {len(METHODS)} methods x 2 field settings")
    print(f"  probe: {PROBE_Q:+.1f} at {BURSCHOWSKY_R} A from O3, outward along centroid->O3")
    print()
    print("  ORCA is not on PATH on this cluster and cannot run on a login node.")
    print(f"  Submit with:  qsub {pbs.resolve()}")


def final_energy(path):
    """Bare single point: the plain FINAL SINGLE POINT ENERGY line, no (MM) or (QM/MM)."""
    t = Path(path).read_text(errors="replace")
    if "ORCA TERMINATED NORMALLY" not in t:
        return None
    hits = re.findall(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", t)
    return float(hits[-1]) if hits else None


def analyse(work, tsv):
    import statistics as st
    work = Path(work)
    dft = {}
    hdr = None
    for l in Path(tsv).read_text().splitlines():
        if l.startswith("#"):
            continue
        f = l.rstrip("\n").split("\t")
        if f[0] == "frame":
            hdr = f
            continue
        r = dict(zip(hdr, f))
        try:
            dft[r["frame"]] = float(r["barrier_vac"])
        except (ValueError, KeyError):
            pass

    E = {}
    for out in work.glob("*.out"):
        m = re.match(r"(\d+)_(R|TS)_(.+?)_(bare|field)\.out$", out.name)
        if not m:
            continue
        E[(m.group(1), m.group(2), m.group(3), m.group(4))] = final_energy(out)

    frames = sorted({k[0] for k in E})
    print(f"analysing {len(frames)} frames\n")

    def barrier(fr, meth, field):
        r = E.get((fr, "R", meth, field))
        ts = E.get((fr, "TS", meth, field))
        return None if r is None or ts is None else (ts - r) * HARTREE

    def pearson(a, b):
        ma, mb = st.mean(a), st.mean(b)
        num = sum((x - ma) * (y - mb) for x, y in zip(a, b))
        da = math.sqrt(sum((x - ma) ** 2 for x in a))
        db = math.sqrt(sum((y - mb) ** 2 for y in b))
        return num / (da * db) if da and db else float("nan")

    print("TEST A - bare barrier against the committed DFT barrier_vac")
    print(f"  {'method':<12}{'n':>4}{'mean':>9}{'DFT mean':>10}{'offset':>9}"
          f"{'r':>8}{'max dev':>9}  verdict")
    refbar = {}
    for meth in METHODS:
        pts = [(fr, barrier(fr, meth, "bare")) for fr in frames]
        pts = [(fr, b) for fr, b in pts if b is not None and fr in dft]
        if len(pts) < 3:
            print(f"  {meth:<12}{len(pts):>4}   insufficient completed jobs")
            continue
        mine = [b for _, b in pts]
        theirs = [dft[fr] for fr, _ in pts]
        refbar[meth] = dict(pts)
        r = pearson(mine, theirs)
        off = st.mean(mine) - st.mean(theirs)
        dev = max(abs(m - t - off) for m, t in zip(mine, theirs))
        ok = "USABLE" if r >= 0.95 else "FAILS A"
        print(f"  {meth:<12}{len(pts):>4}{st.mean(mine):>9.2f}{st.mean(theirs):>10.2f}"
              f"{off:>+9.2f}{r:>8.3f}{dev:>9.2f}  {ok}")

    print("\nTEST B - response to a +1 charge 3.2 A from the ether oxygen")
    print("  (Burschowsky et al. PNAS 2014: Arg90->citrulline, NH2 3.2 A from the ether")
    print("   oxygen, worth 5.9 kcal/mol of TS stabilisation)")
    ref_resp = {}
    for fr in frames:
        b = barrier(fr, "ref", "bare")
        f_ = barrier(fr, "ref", "field")
        if b is not None and f_ is not None:
            ref_resp[fr] = f_ - b
    if ref_resp:
        v = list(ref_resp.values())
        print(f"\n  DFT response: mean {st.mean(v):+.2f} kcal/mol, "
              f"range {min(v):+.2f} to {max(v):+.2f}, n={len(v)}")
        print("  A negative response means the charge LOWERS the barrier, i.e. catalyses.")
    print(f"\n  {'method':<12}{'n':>4}{'mean resp':>11}{'vs DFT':>9}"
          f"{'worst rel err':>15}{'sign ok':>9}  verdict")
    for meth in METHODS:
        if meth == "ref":
            continue
        pts = []
        for fr in frames:
            if fr not in ref_resp:
                continue
            b = barrier(fr, meth, "bare")
            f_ = barrier(fr, meth, "field")
            if b is not None and f_ is not None:
                pts.append((ref_resp[fr], f_ - b))
        if len(pts) < 3:
            print(f"  {meth:<12}{len(pts):>4}   insufficient completed jobs")
            continue
        mine = [m for _, m in pts]
        theirs = [t for t, _ in pts]
        rel = max(abs(m - t) / abs(t) for t, m in pts if abs(t) > 1e-6)
        sign = all((m < 0) == (t < 0) for t, m in pts)
        ok = "USABLE" if (rel <= 0.25 and sign) else "SCREEN ONLY"
        print(f"  {meth:<12}{len(pts):>4}{st.mean(mine):>+11.2f}"
              f"{st.mean(mine)-st.mean(theirs):>+9.2f}{100*rel:>14.0f}%"
              f"{str(sign):>9}  {ok}")

    print("\nRULE, fixed before running:")
    print("  Test A r >= 0.95                      -> usable as a screen")
    print("  Test B within 25% and same sign       -> usable for scoring a design")
    print("  Passing A but failing B               -> pre-filter only, never scoring")
    print("\nNote: this benchmark EVALUATES on fixed DFT geometries. AM1's documented")
    print("failure for this reaction is geometric (substrate compression, Claeyssens 2011")
    print("p.11). Nothing here licenses letting a cheap method OPTIMISE geometries.")


if __name__ == "__main__":
    if len(sys.argv) < 2 or sys.argv[1] not in ("generate", "analyse"):
        sys.exit(__doc__)
    if sys.argv[1] == "generate":
        generate(sys.argv[2], sys.argv[3])
    else:
        analyse(sys.argv[2], sys.argv[3])
