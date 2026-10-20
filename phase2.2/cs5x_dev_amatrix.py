#!/usr/bin/env python3
"""
cs5x_dev_amatrix.py - EXPLORATORY (development, not Stage 5 production): the first-order coefficient a of every site of
the development grid (phase2.2/cs_stage5/exploratory/dev_sites.tsv) in each of the 30 Phase 1 band frames, from the
existing in vacuo densities of 05_qmmm/20_invacuo (B3LYP-D3BJ/def2-SVP, the densities A_v2 was built from). No SCF is
run: orca_vpot evaluates each frame's R and TS density at the sites, taken back to that frame's own coordinates by the
inverse of its alignment transform (phase2.2/aligned/frame_<F>/transform.txt), as s15_dv_matrix.py did for A_v2.
a = (V_TS - V_R) x 627.5094740631, kcal/mol per e; a < 0: a positive charge there lowers the barrier.

Purpose: a realistic 30-frame test set (sites, L in exploratory/dev_L.tsv, and this A) for building and testing the
optimiser's code while the SP package runs. No optimiser setting may be chosen by the barriers its designs reach here.

  python3 cs5x_dev_amatrix.py prepare  <results repo> <OUT> [HPC_DIR]
  qsub OUT/dev_A.pbs                                   (hpc1, with OUT at HPC_DIR)
  python3 cs5x_dev_amatrix.py assemble <results repo> <WORK> <OUT_RESULTS>

Rows: the 2769 development sites, then grid_v2's 233 sites as a positive control. Checks in assemble (all must pass,
or nothing is written): the inputs are what prepare writes (byte for byte); each vpot.out holds one value per point at
the points written (1e-6 Bohr); the control rows reproduce A_v2 within 1e-4 kcal/mol per e (A_v2 converted Angstrom to
Bohr with 1.8897259886, this script with ORCA's 1/0.5291772083: s15 estimated the effect at no more than 3.2e-5).
prepare also checks that every frame's 20_invacuo R and TS geometry, moved by its transform, is the aligned geometry
(1e-6 A), so the densities belong to the frames the grid was built on, and writes densities.sha256: the committed sha256
(05_qmmm/20_invacuo/CHECKSUMS_binaries.sha256) of the 120 .gbw and .densities files the job reads. The job checks them on
hpc1 before any orca_vpot call and stops if one differs; assemble requires its record (densities_check.txt) to show all
120 as OK. Python 3.6+, numpy.
"""
import hashlib, os, shutil, sys, tempfile
from pathlib import Path
import numpy as np

H = 627.5094740631
ANG2BOHR = 1.0 / 0.5291772083                    # ORCA 6.0.1's factor (s15_dv_matrix.py's note)
DEV = "phase2.2/cs_stage5/exploratory/dev_sites.tsv"
V2, AV2 = "phase2.2/grid_v2.tsv", "phase2.2/A_v2.tsv"
ALN, VAC = "phase2.2/aligned", "05_qmmm/20_invacuo"
HPC_DEFAULT = "/home/18660916/system_development/phase2.2/cs5x_dev_A"
VAC_HPC = "/home/18660916/system_development/05_qmmm/20_invacuo"

def stop(m): sys.stderr.write("STOP: " + m + "\n"); sys.exit(1)
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def xyz(p):
    L = Path(p).read_text().splitlines(); n = int(L[0].split()[0])
    return np.array([[float(v) for v in l.split()[1:4]] for l in L[2:2 + n]])

def sites(p, idcol):
    rows, h = [], None
    for l in Path(p).read_text().splitlines():
        if not l.strip() or l.startswith("#"): continue
        t = l.split("\t")
        if t[0] == idcol: h = t; continue
        rows.append([float(t[h.index(k)]) for k in ("x", "y", "z")])
    return np.array(rows)

def frames_of(res):
    f = sorted((d.name[6:] for d in (Path(res) / ALN).iterdir() if d.name.startswith("frame_")), key=int)
    A2f = [l for l in (Path(res) / AV2).read_text().splitlines() if l.startswith("# frames")][0].split()[3:]
    if f != A2f: stop("the aligned frames differ from A_v2's")
    return f

def build(res, out, hpc):
    """write every input into out (which must not exist); deterministic"""
    res, out = Path(res), Path(out)
    frames = frames_of(res)
    Q = np.vstack([sites(res / DEV, "site"), sites(res / V2, "idx")])
    nd = len(Q) - 233
    if nd != 2769: stop("dev_sites.tsv does not hold the 2769 development sites")
    out.mkdir(parents=True)
    for f in frames:
        M = np.loadtxt(res / ALN / ("frame_" + f) / "transform.txt"); R, t = M[:3], M[3]
        for s, n in (("R", "reactant_qm"), ("TS", "transition_state_qm")):
            if np.abs(xyz(res / VAC / ("%s_%s.xyz" % (f, s))) @ R + t - xyz(res / ALN / ("frame_" + f) / (n + ".xyz"))).max() > 1e-6:
                stop("frame %s %s: the 20_invacuo geometry is not the aligned geometry moved back" % (f, s))
        D = Q - t                                            # (Q - t) @ R.T, written out element by element: matrix
        P = np.stack([(D[:, 0] * R[k, 0] + D[:, 1] * R[k, 1]) + D[:, 2] * R[k, 2] for k in range(3)], 1) * ANG2BOHR
        #   products go through BLAS, whose rounding differs between machines; these operations do not
        with open(out / ("points_%s.xyz" % f), "w") as fh:
            fh.write("%d\n" % len(P)); fh.write("".join("%.10f %.10f %.10f\n" % tuple(p) for p in P))
    (out / "frames.txt").write_text("\n".join(frames) + "\n")
    ck = {}
    for l in (res / VAC / "CHECKSUMS_binaries.sha256").read_text().splitlines():
        h, n = l.split("  ", 1); ck[n] = h
    want = ["sp_%s_%s.%s" % (f, s, e) for f in frames for s in ("R", "TS") for e in ("gbw", "densities")]
    if any(n not in ck for n in want): stop("CHECKSUMS_binaries.sha256 lacks a density of the development frames")
    (out / "densities.sha256").write_text("".join("%s  %s\n" % (ck[n], n) for n in want))
    (out / "row_order.txt").write_text("# rows 0-%d: %s (site 0-%d); rows %d-%d: %s (idx 0-232), control\n# sha256 %s %s; %s %s\n" % (
        nd - 1, DEV, nd - 1, nd, nd + 232, V2, DEV, sha(res / DEV), V2, sha(res / V2)))
    pbs = """#!/bin/bash
#PBS -N cm_s5xdevA
#PBS -l select=1:ncpus=1:mem=8gb
#PBS -l walltime=04:00:00
#PBS -m ae
#PBS -M 18660916@sun.ac.za
#PBS -j oe
#PBS -o {W}/dev_A.pbs.out
# EXPLORATORY development A matrix (cs5x_dev_amatrix.py): orca_vpot on the committed 20_invacuo densities, no SCF.
# The orca_vpot call is s15_dv_matrix.py's (verified 30 September 2026). Restartable.
set -uo pipefail
ORCA=/home/apps2/ORCA/6.0.1
VAC={VAC}
W={W}
export PATH="/apps/openmpi/4.1.1/bin:$PATH"
export LD_LIBRARY_PATH="$ORCA/lib:/apps/openmpi/4.1.1/lib:/apps/mambaforge/envs/medaka/lib:${{LD_LIBRARY_PATH:-}}"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
cd "$W" || exit 1
if ! mkdir .running 2>/dev/null; then echo "FAIL: another job holds .running"; exit 1; fi
trap 'rmdir "$W/.running" 2>/dev/null' EXIT
echo "host=$(hostname) start=$(date)"
# the densities must be the committed ones (05_qmmm/20_invacuo/CHECKSUMS_binaries.sha256) before anything is evaluated
( cd "$VAC" && sha256sum -c "$W/densities.sha256" ) > densities_check.txt 2>&1 || {{ echo "FAIL: a density on hpc1 is not the committed one (densities_check.txt)"; exit 1; }}
echo "densities: $(grep -c ': OK$' densities_check.txt) of 120 match the committed sha256"
ok=0; miss=0
for f in $(cat frames.txt); do
  for tag in R TS; do
    gbw="$VAC/sp_${{f}}_${{tag}}.gbw"; out="vpot_${{f}}_${{tag}}.out"
    [ -s "$gbw" ] || {{ echo "MISSING $gbw"; miss=$((miss+1)); continue; }}
    [ -s "$out" ] && {{ ok=$((ok+1)); continue; }}
    "$ORCA/orca_vpot" "$gbw" "sp_${{f}}_${{tag}}.scfp" "points_${{f}}.xyz" "$out" "$VAC/sp_${{f}}_${{tag}}" > /dev/null 2>&1
    if [ -s "$out" ]; then ok=$((ok+1)); else echo "FAILED $f $tag"; miss=$((miss+1)); fi
  done
done
echo "evaluated=$ok missing_or_failed=$miss"
echo "end=$(date)"
""".format(W=hpc, VAC=VAC_HPC)
    (out / "dev_A.pbs").write_text(pbs)
    files = sorted(p for p in out.rglob("*") if p.is_file())
    (out / "INPUT_CHECKSUMS.sha256").write_text("".join("%s  %s\n" % (sha(p), p.relative_to(out)) for p in files))
    return frames, nd

def prepare(res, out, hpc=HPC_DEFAULT):
    if os.path.exists(out): stop("%s exists" % out)
    frames, nd = build(res, out, hpc)
    print("prepared %s: %d frames x R, TS; %d development sites + 233 control sites; submit: qsub %s/dev_A.pbs" % (out, len(frames), nd, hpc))

def vpot(p, pts):
    P = np.array([[float(v) for v in l.split()] for l in Path(pts).read_text().splitlines()[1:] if l.strip()])
    L = [l.split() for l in Path(p).read_text().splitlines() if l.strip()]
    if not L or len(L[0]) != 1 or int(L[0][0]) != len(P): stop("%s: count line" % p)
    V = np.array([[float(v) for v in l] for l in L[1:]])
    if V.shape != (len(P), 4) or np.abs(V[:, :3] - P).max() > 1e-6: stop("%s: points differ from %s" % (p, Path(pts).name))
    return V[:, 3]

def assemble(res, work, out):
    if os.path.exists(out): stop("%s exists" % out)
    work = Path(work)
    hpc = [l for l in (work / "dev_A.pbs").read_text().splitlines() if l.startswith("W=")][0][2:]
    tmp = tempfile.mkdtemp()
    try:
        ref = Path(tmp) / "inputs"
        frames, nd = build(res, ref, hpc)                    # the inputs, made again from the committed files
        for l in (ref / "INPUT_CHECKSUMS.sha256").read_text().splitlines():
            h, f = l.split("  ", 1)
            if not (work / f).is_file() or sha(work / f) != h: stop("%s in %s is not what prepare writes" % (f, work))
        if (work / "INPUT_CHECKSUMS.sha256").read_bytes() != (ref / "INPUT_CHECKSUMS.sha256").read_bytes(): stop("INPUT_CHECKSUMS.sha256 differs")
    finally:
        shutil.rmtree(tmp)
    want = [l.split("  ", 1)[1] for l in (work / "densities.sha256").read_text().splitlines()]
    got = (work / "densities_check.txt").read_text().splitlines() if (work / "densities_check.txt").is_file() else []
    if sorted(got) != sorted(n + ": OK" for n in want): stop("densities_check.txt does not show the 120 densities as the committed ones")
    for f in frames:
        for s in ("R", "TS"):
            if not (work / ("vpot_%s_%s.out" % (f, s))).is_file(): stop("vpot_%s_%s.out is missing: see dev_A.pbs.out" % (f, s))
    cols = []
    for f in frames:
        pts = work / ("points_%s.xyz" % f)
        cols.append((vpot(work / ("vpot_%s_TS.out" % f), pts) - vpot(work / ("vpot_%s_R.out" % f), pts)) * H)
    A = np.stack(cols, 1)
    A2 = np.loadtxt(Path(res) / AV2, comments="#", ndmin=2)
    d = np.abs(A[nd:] - A2).max()
    if d > 1e-4: stop("the control rows do not reproduce A_v2: max |difference| %.2e kcal/mol per e" % d)
    Ad = A[:nd]
    rows = [l.split("\t") for l in (Path(res) / DEV).read_text().splitlines() if not l.startswith("#")][1:]
    lay = np.array([r[1] for r in rows])
    os.makedirs(out)
    head = ["# EXPLORATORY development A matrix, cs5x_dev_amatrix.py (sha256 %s): NOT Stage 5 production" % sha(__file__),
            "# a_if = V_TS,f(r_i) - V_R,f(r_i), kcal/mol per e, from the committed 05_qmmm/20_invacuo densities (B3LYP-D3BJ/def2-SVP)",
            "#   of the Phase 1 band frames; a < 0: a positive charge at site i lowers frame f's barrier (first order)",
            "# rows: %s (row i = site i; sha256 %s)" % (DEV, sha(Path(res) / DEV)),
            "# frames (%d): %s" % (len(frames), " ".join(frames)),
            "# control: grid_v2's 233 sites reproduce A_v2 within %.1e kcal/mol per e" % d]
    with open(Path(out) / "A_dev.tsv", "x") as fh:
        fh.write("\n".join(head) + "\n")
        for r in Ad: fh.write("\t".join(("%.6f" % v).replace("-0.000000", "0.000000") for v in r) + "\n")
    rep = head[:5] + ["",
           "control rows against A_v2: max |difference| %.2e kcal/mol per e (limit 1e-4)" % d,
           "range %.3f to %.3f kcal/mol per e; a < 0 in %.1f%% of entries (sites x frames)" % (Ad.min(), Ad.max(), 100 * (Ad < 0).mean())]
    for L_ in ("anion", "cation"):
        m = lay == L_; Am = Ad[m]
        rep.append("%-6s %4d sites: per-site mean a %.3f to %.3f; every frame a < 0 at %d sites, a > 0 at %d, mixed at %d" % (
            L_, m.sum(), Am.mean(1).min(), Am.mean(1).max(), (Am < 0).all(1).sum(), (Am > 0).all(1).sum(),
            m.sum() - (Am < 0).all(1).sum() - (Am > 0).all(1).sum()))
    Path(out, "DEV_A_REPORT.txt").write_text("\n".join(rep) + "\n")
    print("\n".join(rep))

if __name__ == "__main__":
    a = sys.argv[1:]
    if a[:1] == ["prepare"] and len(a) in (3, 4): prepare(*a[1:])
    elif a[:1] == ["assemble"] and len(a) == 4: assemble(*a[1:])
    else: stop(__doc__.split("\n\n")[2])
