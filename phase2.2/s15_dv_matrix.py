#!/usr/bin/env python3
"""
s15_dv_matrix.py - build the multi-frame difference-potential matrix that s14 optimises.

Reworked from phase2b_charge_design/03_dvpot/s3_dv_on_grid_v2.py, which evaluated the
difference potential for ONE frame. Four changes, each for a stated reason.

  1. MULTI-FRAME. s14 needs a_if for every site i and every frame f, not a single
     vector. The output is an (S, F) matrix in the row order of the grid file, which is
     what s14's solve mode reads.

  2. POST-CUT FRAMES ONLY. The committed dv_grid_ensemble_mean.tsv states in its own
     header that it was built from ten frames, all below the 20000 ps equilibration cut
     and all at the 400 kJ spring. Phase 1's finalised position is that the pilot set
     must not be used for reported results, so the frame list here defaults to the
     post-cut frames that carry an in vacuo reference.

  3. IN VACUO DENSITIES, NOT QM/MM. The Phase 2.2 design environment is bare substrate
     plus charges, with no protein and no continuum ([D2018] never runs charges and
     COSMO together; [B2021] uses explicit molecules and no continuum). So the
     potentials must come from the BARE substrate densities at the QM/MM-optimised
     geometries - exactly the in vacuo single points of step 19 - and not from the
     QM/MM densities, which already contain the protein field the design is meant to
     replace.

  4. UNITS AND SIGN MADE EXPLICIT IN THE OUTPUT. The matrix is written in
     kcal/mol per unit charge so that s14's objective is directly in kcal/mol.

THE SIGN CONVENTION, AND WHY IT IS THE ONE s14 NEEDS
Dv_i = V_TS(r_i) - V_R(r_i). The interaction of a charge q at r_i with a fixed density
is qV(r_i), so the differential stabilisation is

    ddE_i(q) = [E_TS(q) - E_TS(0)] - [E_R(q) - E_R(0)] = q * Dv_i

which is exactly s14's linear-response model with a_i = Dv_i. Dv < 0 therefore means a
POSITIVE charge at that site lowers the barrier. This is the convention the existing
s1_signcheck.py validates with both a positive and a negative control, and it is
preserved unchanged here.

NO NEW SCF IS PERFORMED. orca_vpot only evaluates the potential of an existing density
at new positions, so the map is directly comparable to the barriers those same densities
define.

WHAT THIS DOES NOT PRODUCE
The second-order coefficient b_i. The measured b is the differential polarisation and
needs SCF at two or more charge magnitudes per site, which orca_vpot cannot give. s14
runs on the linear model and says so; see the plan's section on the Tier 2 extension.

USAGE, on the HPC where the in vacuo .gbw files live
    python3 s15_dv_matrix.py prepare <grid.tsv> <workdir>
    qsub <workdir>/s15_dv.pbs
    python3 s15_dv_matrix.py assemble <grid.tsv> <workdir> <out_amatrix.tsv>
"""
import sys
import math
from pathlib import Path

ANG2BOHR = 1.8897259886
HARTREE2KCAL = 627.5094740631

# The frame set is READ FROM THE GRID'S OWN HEADER, not hardcoded. s12_build_grid.py
# writes "# frames enclosed (N): 20000(2) 21634(2) ..." and that is by construction the
# set the surface covers, so the potentials must be evaluated on exactly those frames.
#
# An earlier version hardcoded a 22-frame list. That is the same staleness that left
# s8_invacuo_new.pbs computing in vacuo references for five frames while the ensemble
# had grown to thirty: a list written once and never revisited. Reading the grid header
# means the two can never disagree.
FALLBACK_FRAMES = ["20000", "21634", "23268", "24883", "26495", "33320", "34991",
                   "36665", "38344", "40485", "41786", "42436", "43087", "43738",
                   "45688", "46990", "47641", "48292", "49592", "50892", "52192",
                   "53493"]


def frames_from_grid(path):
    """Read the enclosed-frame list that s12_build_grid.py records in its header."""
    import re
    for l in Path(path).read_text().splitlines():
        if l.startswith("# frames enclosed"):
            return re.findall(r"(\d+)\(\d+\)", l)
        if not l.startswith("#"):
            break
    return None


def read_grid(path):
    """Grid written by s12_build_grid.py: idx, x, y, z, nearest_atom_A."""
    sites = []
    for l in Path(path).read_text().splitlines():
        if l.startswith("#") or l.startswith("idx"):
            continue
        f = l.split("\t")
        if len(f) >= 4:
            sites.append((int(f[0]), float(f[1]), float(f[2]), float(f[3])))
    return sites


def prepare(grid_path, work):
    work = Path(work)
    work.mkdir(parents=True, exist_ok=True)
    sites = read_grid(grid_path)
    if not sites:
        sys.exit(f"no sites parsed from {grid_path}")

    # orca_vpot returns values POSITIONALLY with no site identifiers. If the row order
    # were lost the whole map would be silently scrambled, so it is written separately
    # and checked on assembly.
    with open(work / "points_bohr.xyz", "w") as fh:
        fh.write(f"{len(sites)}\n")
        for _, x, y, z in sites:
            fh.write(f"{x*ANG2BOHR:.10f} {y*ANG2BOHR:.10f} {z*ANG2BOHR:.10f}\n")
    with open(work / "row_order.tsv", "w") as fh:
        fh.write("row\tidx\tx\ty\tz\n")
        for r, (i, x, y, z) in enumerate(sites):
            fh.write(f"{r}\t{i}\t{x:.6f}\t{y:.6f}\t{z:.6f}\n")

    frames = frames_from_grid(grid_path)
    if frames:
        print(f"frame set read from the grid header: {len(frames)} frames")
    else:
        frames = FALLBACK_FRAMES
        print(f"  WARNING: no frame list in the grid header; falling back to the")
        print(f"  hardcoded {len(frames)}-frame list, which may be stale. Rebuild the")
        print(f"  grid with s12_build_grid.py so the two cannot disagree.")
    (work / "frames.txt").write_text("\n".join(frames) + "\n")
    pbs = work / "s15_dv.pbs"
    pbs.write_text(f"""#!/usr/bin/env bash
#PBS -N cm_s15dv
#PBS -l select=1:ncpus=1:mem=8gb
#PBS -l walltime=04:00:00
#PBS -m ae
#PBS -M 18660916@sun.ac.za
#PBS -j oe
#PBS -o {work.resolve()}/s15_dv.pbs.out
#
# orca_vpot on the committed in vacuo densities. No SCF is run.
# One call per frame per state; each is seconds, so this is serial and unparallelised.

set -uo pipefail
ORCA=/home/apps2/ORCA/6.0.1
VAC=/home/18660916/system_development/05_qmmm/20_invacuo
WORK={work.resolve()}

export PATH="/apps/openmpi/4.1.1/bin:$PATH"
export LD_LIBRARY_PATH="$ORCA/lib:/apps/openmpi/4.1.1/lib:/apps/mambaforge/envs/medaka/lib:${{LD_LIBRARY_PATH:-}}"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1

cd "$WORK" || exit 1
if ! mkdir .running 2>/dev/null; then echo "FAIL: another job holds .running"; exit 1; fi
trap 'rmdir "$WORK/.running" 2>/dev/null' EXIT
echo "host=$(hostname) start=$(date)"
echo "frames: {' '.join(frames)}"
echo

miss=0; ok=0
for f in {' '.join(frames)}; do
  for tag in R TS; do
    gbw="$VAC/sp_${{f}}_${{tag}}.gbw"
    out="vpot_${{f}}_${{tag}}.out"
    if [ ! -s "$gbw" ]; then
      echo "  MISSING $gbw"; miss=$((miss+1)); continue
    fi
    [ -s "$out" ] && {{ ok=$((ok+1)); continue; }}
    # ORCA 6 stores the SCF density INSIDE the .densities container rather than as a
    # loose .scfp file, and orca_vpot's second argument is the density's NAME in that
    # container, not a path. Verified 2026-09-30 on sp_24883_R: the call returns in
    # 0.4 s and writes four columns with the potential last. There is therefore no
    # .scfp on disk and none is needed; testing for one rejects every frame.
    # The manual's note that a mismatched basename must be passed as a FIFTH argument
    # is the tell: the second argument is a container entry.
    "$ORCA/orca_vpot" "$gbw" "sp_${{f}}_${{tag}}.scfp" \\
        points_bohr.xyz "$out" "$VAC/sp_${{f}}_${{tag}}" > /dev/null 2>&1
    if [ -s "$out" ]; then ok=$((ok+1)); else echo "  FAILED $f $tag"; miss=$((miss+1)); fi
  done
done
echo
echo "evaluated=$ok missing_or_failed=$miss"
echo "end=$(date)"
echo "now: python3 s15_dv_matrix.py assemble <grid.tsv> {work.resolve()} amatrix.tsv"
""")
    pbs.chmod(0o755)
    print(f"{len(sites)} sites written to {work/'points_bohr.xyz'} (Bohr)")
    print(f"row order recorded in {work/'row_order.tsv'}")
    print(f"{len(frames)} frames, 2 states each = {2*len(frames)} orca_vpot calls")
    print(f"\nsubmit:  qsub {pbs.resolve()}")
    print("\nNOTE: no .scfp file is needed. ORCA 6 keeps the SCF density inside the")
    print(".densities container and orca_vpot takes the density NAME as its second")
    print("argument, with the container basename as an optional fifth. The committed")
    print("in vacuo single points already carry everything required.")


def read_vpot(path, n):
    """orca_vpot output: a count line, then one potential per point in Hartree/e."""
    lines = [l for l in Path(path).read_text().splitlines() if l.strip()]
    vals = []
    for l in lines:
        f = l.split()
        if len(f) == 1 and f[0].isdigit():
            continue
        try:
            vals.append(float(f[-1]))
        except ValueError:
            continue
    if len(vals) != n:
        raise ValueError(f"{path}: {len(vals)} values, expected {n}")
    return vals


def assemble(grid_path, work, out_path):
    work, out_path = Path(work), Path(out_path)
    sites = read_grid(grid_path)
    n = len(sites)

    order = [l.split("\t") for l in
             (work / "row_order.tsv").read_text().splitlines()[1:] if l.strip()]
    if len(order) != n:
        sys.exit(f"row_order has {len(order)} rows but the grid has {n}. The grid has "
                 f"changed since prepare was run; rerun prepare.")
    for r, (i, x, y, z) in enumerate(sites):
        if int(order[r][1]) != i or abs(float(order[r][2]) - x) > 1e-6:
            sys.exit(f"row {r} does not match row_order.tsv. orca_vpot returns values "
                     f"positionally, so the map would be scrambled. Rerun prepare.")

    fl = work / "frames.txt"
    if fl.exists():
        want = [x for x in fl.read_text().split() if x]
    else:
        want = frames_from_grid(grid_path) or FALLBACK_FRAMES
    frames, cols, missing = [], [], []
    for fr in want:
        pr, pts = work / f"vpot_{fr}_R.out", work / f"vpot_{fr}_TS.out"
        if not (pr.exists() and pts.exists()):
            missing.append(fr)
            continue
        vR, vTS = read_vpot(pr, n), read_vpot(pts, n)
        cols.append([(vTS[k] - vR[k]) * HARTREE2KCAL for k in range(n)])
        frames.append(fr)

    if not cols:
        sys.exit("no frame produced both potentials; nothing to assemble")
    F = len(cols)

    hdr = ["# a-matrix for s14_design_milp.py, produced by s15_dv_matrix.py",
           "# a_if = [V_TS(r_i) - V_R(r_i)] for frame f, in kcal/mol per unit charge.",
           "# ddE_f(q) = sum_i a_if q_i, so a < 0 means a POSITIVE charge there lowers",
           "#   the barrier. Convention validated by 03_dvpot/s1_signcheck.py, which",
           "#   uses both a positive and a negative control.",
           "# Densities: IN VACUO single points at the QM/MM-optimised geometries, i.e.",
           "#   bare substrate with no protein and no continuum, matching the Phase 2.2",
           "#   design environment. NOT the QM/MM densities, which already contain the",
           "#   protein field the design is meant to replace.",
           "# No SCF was run; orca_vpot evaluates existing densities at new positions.",
           "#",
           f"# sites:  {n}   (row order matches the grid file exactly)",
           f"# frames: {F}   {' '.join(frames)}",
           f"# grid:   {Path(grid_path).name}"]
    if missing:
        hdr.append(f"# MISSING, excluded: {' '.join(missing)}")
    hdr.append("# ONLY the linear coefficient is here. The measured second-order term")
    hdr.append("#   b (differential polarisation) is NOT included; s14 runs the linear")
    hdr.append("#   model and reports that scope with every result.")

    with open(out_path, "w") as fh:
        fh.write("\n".join(hdr) + "\n")
        for k in range(n):
            fh.write("\t".join(f"{cols[f][k]:.6f}" for f in range(F)) + "\n")

    flat = [v for c in cols for v in c]
    stab = sum(1 for v in flat if v < 0)
    print(f"wrote {out_path}: {n} sites x {F} frames")
    if missing:
        print(f"  EXCLUDED {len(missing)} frames with no potentials: {' '.join(missing)}")
    print(f"  range  {min(flat):+.4f} to {max(flat):+.4f} kcal/mol per unit charge")
    print(f"  stabilising for a POSITIVE charge (a < 0): {stab}/{len(flat)} "
          f"= {100*stab/len(flat):.1f}%")
    print(f"  a value near 50% is expected: the difference potential changes sign across")
    print(f"  space, so roughly half the grid favours a negative charge instead.")

    per = [sum(1 for v in c if v < 0) / n for c in cols]
    print(f"  per-frame stabilising fraction: min {min(per)*100:.0f}%, "
          f"max {max(per)*100:.0f}%")
    # site-level agreement across frames: the question the clustering test asks
    agree = sum(1 for k in range(n)
                if all(cols[f][k] < 0 for f in range(F))
                or all(cols[f][k] > 0 for f in range(F)))
    print(f"  sites where ALL frames agree on sign: {agree}/{n} = {100*agree/n:.1f}%")
    print(f"  Sites where frames disagree on the SIGN cannot be used by any single")
    print(f"  design; the measured sign reversal at one frame beyond 4 A is this effect.")


def main():
    if len(sys.argv) < 4 or sys.argv[1] not in ("prepare", "assemble"):
        sys.exit(__doc__)
    if sys.argv[1] == "prepare":
        prepare(sys.argv[2], sys.argv[3])
    else:
        if len(sys.argv) < 5:
            sys.exit(__doc__)
        assemble(sys.argv[2], sys.argv[3], sys.argv[4])
    return 0


if __name__ == "__main__":
    sys.exit(main())
