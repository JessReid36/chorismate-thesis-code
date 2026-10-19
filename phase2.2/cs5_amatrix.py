#!/usr/bin/env python3
"""
cs5_amatrix.py - Stage 5 (phase2.2, STAGE5_CRITERIA.txt): the first-order coefficient a of every grid site in every
frame, a_if = V_TS,f(r_i) - V_R,f(r_i) in kcal/mol per e, from in vacuo B3LYP-D3BJ/def2-SVP densities of the bare
substrate at the frame's TS and R (R_IRC), as s15_dv_matrix.py built A_v2.tsv (PHASE22_PLAN_UPDATE_20261009.md
section 2; Stage 4e kept this level: SIGN ROBUST). Sokalski (1985), Eq. 19 (p. 404): for a set of point charges
q_C, the barrier change is sum_C q_C (V_TC - V_SC); "the difference of molecular potentials ... represents the
catalytic activity of a single unitary positive point charge" and serves "for any fractional point charge, including
negative charges" (p. 405). So a < 0: a positive charge there lowers the barrier.

  python3 cs5_amatrix.py prepare  <aligned_dir> <sites.tsv> <OUT> [HPC_DIR]   (PC; writes the inputs, to commit)
  qsub OUT/stage5_a.pbs                                                        (hpc1, with OUT copied to HPC_DIR)
  python3 cs5_amatrix.py assemble <aligned_dir> <sites.tsv> <WORK> <OUT_RESULTS>   (PC, on the fetched outputs)

Differences from s15, each for a stated reason:
  - new single points, run on the ALIGNED geometries (cs5_align.py), so the grid needs no transform back (s15 reused
    densities computed in each frame's own coordinates); the in vacuo energy and potential do not depend on where the
    molecule sits, and check (b) below measures that on one frame;
  - points converted to Bohr with ORCA 6.0.1's own factor, 1/0.5291772083 (s15 kept 1.8897259886 for A_v2's
    reproducibility; the difference, 7.7e-8 relative, moves a point at most 6.4e-6 A);
  - V at P is also evaluated (AP = V_P - V_R), for the optimiser stage's open question on path terms; Stage 5 uses
    only A;
  - the sites file is cs5_grid.py's (two layers), read by its header.
Input of every single point (as 05_qmmm/20_invacuo and s17_verify_amatrix.sh):
  ! B3LYP D3BJ def2-SVP def2/J RIJCOSX TightSCF / %maxcore 3000 / %pal nprocs 8 end / %scf MaxIter 300 end /
  * xyzfile -2 1 substrate.xyz
CHECKS (assemble; all must pass, or A.tsv is not written):
  (a) every run: ORCA TERMINATED NORMALLY, SCF CONVERGED, 118 electrons, the echoed coordinates equal substrate.xyz
      (1e-5 A); every vpot.out has one value per site at the points written (1e-6 Bohr), so the row order holds;
  (b) rigid-motion invariance: the alignment frame with the largest rotation is also run in its own (original)
      coordinates, with the sites taken back by its inverse transform; max |V_original - V_aligned| <= 1e-4 Eh/e
      (0.06 kcal/mol per e; the DFT and COSX integration grids are fixed in space, so not exactly invariant);
  (c) finite difference: at the first frame, three sites (those nearest O3, O5 and O1 at R), charges of +-0.05 e as
      ORCA point charges at R and TS: a_fd = [E_TS(+) - E_TS(-) - E_R(+) + E_R(-)] / (2 x 0.05) must equal a within
      0.01 kcal/mol per e + 1% of |a| (Hellmann-Feynman: dE/dq = V for the variational SCF energy; the odd third-order
      term is about 1e-4 of a at 0.05 e, from Stage 4e's 3-4% at 1 e);
  (d) every input is unchanged since prepare (INPUT_CHECKSUMS.sha256); runs.tsv and row_order.tsv are what prepare
      would write for these geometries and sites; each run's substrate.xyz is the aligned geometry (byte for byte; the
      invariance run's within 1e-8 A of the inverse transform), and the sites file was built on these geometries.
A_REPORT.txt also tabulates each frame's in vacuo energies at R, TS and P and the in vacuo barrier E(TS) - E(R).
"""
import hashlib, math, os, re, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cs5_grid as G

H = 627.5094740631                     # kcal/mol per Eh (as s15)
ANG2BOHR = 1.0 / 0.5291772083          # ORCA 6.0.1's factor (s15_dv_matrix.py's note, fitted to 3168 coordinate pairs)
DELTA = 0.05                           # e, finite-difference charge (check c)
HEAD = "! B3LYP D3BJ def2-SVP def2/J RIJCOSX TightSCF\n%maxcore 3000\n%pal nprocs 8 end\n%scf MaxIter 300 end\n"
HPC_DEFAULT = "/home/18660916/system_development/phase2.2/cs_stage5/a_matrix"
STATES = ("R", "TS", "P")
FD_ATOMS = ("O3", "O5", "O1")

def stop(m): sys.stderr.write("STOP: " + m + "\n"); sys.exit(1)
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def frames_and_sites(adir, sites_p):
    frames = G.frame_names(adir)
    rows, com = G.read_sites(sites_p)
    if not rows or list(rows[0]) != G.HEAD: stop("%s is not a cs5_grid.py sites file" % sites_p)
    if G.frames_of(com) != frames: stop("the sites file's frames differ from %s's" % adir)
    G.same_geometries(com, adir, [(f, s) for f in frames for s in STATES])
    return frames, rows, G.sites_xyz(rows)

def transform(adir, f):
    M = np.loadtxt(Path(adir) / ("frame_" + f) / "transform.txt"); return M[:3], M[3]

def rot_angle(R): return math.degrees(math.acos(max(-1.0, min(1.0, (np.trace(R) - 1) / 2))))

def write_points(p, X):
    with open(p, "w") as fh:
        fh.write("%d\n" % len(X))
        for x, y, z in X * ANG2BOHR: fh.write("%.10f %.10f %.10f\n" % (x, y, z))

def xyz_text(X, title):
    return "24\n%s\n" % title + "".join("%-3s%18.10f%18.10f%18.10f\n" % (n[0], x[0], x[1], x[2]) for n, x in zip(G.NAMES, X))

def plan(adir, frames, Q):
    """the runs: (name, kind, frame, state, site, charge, points file)"""
    runs = [("sp_%s_%s" % (f, s), "sp", f, s, "", "", "points.xyz") for f in frames for s in STATES]
    fi = max(frames, key=lambda f: (rot_angle(transform(adir, f)[0]), -int(f)))
    runs.append(("inv_%s_R" % fi, "inv", fi, "R", "", "", "points_inv.xyz"))
    f0 = frames[0]
    _, X0 = G.read_xyz(Path(adir) / ("frame_" + f0) / G.STATES["R"])
    fd_sites = []
    for a in FD_ATOMS:
        d = np.linalg.norm(Q - X0[G.NAMES.index(a)], axis=1); i = int(d.argmin())
        if i not in fd_sites: fd_sites.append(i)
    for i in fd_sites:
        for s in ("R", "TS"):
            for tag, q in (("p", DELTA), ("m", -DELTA)):
                runs.append(("fd_%s_%s_%d_%s" % (f0, s, i, tag), "fd", f0, s, str(i), "%+.2f" % q, "-"))
    return runs, fi, f0, fd_sites

def prepare(adir, sites_p, out, hpc=HPC_DEFAULT):
    if os.path.exists(out): stop("%s exists" % out)
    frames, rows, Q = frames_and_sites(adir, sites_p)
    runs, fi, f0, fd_sites = plan(adir, frames, Q)
    out = Path(out); out.mkdir(parents=True)
    write_points(out / "points.xyz", Q)
    Ri, ti = transform(adir, fi)
    write_points(out / "points_inv.xyz", (Q - ti) @ Ri.T)
    with open(out / "row_order.tsv", "w") as fh:
        fh.write("# sites %s sha256 %s\nrow\tsite\tlayer\tx\ty\tz\n" % (sites_p, sha(sites_p)))
        for k, r in enumerate(rows): fh.write("%d\t%s\t%s\t%s\t%s\t%s\n" % (k, r["site"], r["layer"], r["x"], r["y"], r["z"]))
    for name, kind, f, s, site, q, pts in runs:
        d = out / name; d.mkdir()
        src = Path(adir) / ("frame_" + f) / G.STATES[s]
        if kind == "inv":
            _, X = G.read_xyz(src)
            (d / "substrate.xyz").write_text(xyz_text((X - ti) @ Ri.T, "frame %s %s in its original coordinates (inverse of transform.txt)" % (f, s)))
        else:
            (d / "substrate.xyz").write_text(src.read_text())
        inp = HEAD
        if kind == "fd":
            p = Q[int(site)]
            (d / "pc.pc").write_text("1\n%10.6f %18.10f %18.10f %18.10f\n" % (float(q), p[0], p[1], p[2]))
            inp += '%pointcharges "pc.pc"\n'
        (d / "job.inp").write_text(inp + "* xyzfile -2 1 substrate.xyz\n")
    (out / "runs.tsv").write_text("# run\tkind\tframe\tstate\tsite\tcharge_e\tpoints\n" + "".join("\t".join(r) + "\n" for r in runs))
    pbs = ["#!/bin/bash", "#PBS -N cm_s5a", "#PBS -l select=1:ncpus=8:mem=24gb", "#PBS -l walltime=24:00:00",
           "#PBS -m ae", "#PBS -M 18660916@sun.ac.za", "#PBS -j oe", "#PBS -o %s/stage5_a.pbs.out" % hpc,
           "# Stage 5 A matrix (cs5_amatrix.py, STAGE5_CRITERIA.txt): %d in vacuo single points, orca_vpot on the sp and inv runs." % len(runs),
           "# Restartable: a run whose job.out ends normally is not repeated; its vpot.out is made if missing.",
           'export PATH="/apps/openmpi/4.1.1/bin:$PATH"',
           'export LD_LIBRARY_PATH="/home/apps2/ORCA/6.0.1/lib:/apps/openmpi/4.1.1/lib:/apps/mambaforge/envs/medaka/lib:${LD_LIBRARY_PATH:-}"',
           "export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1", "ORCA=/home/apps2/ORCA/6.0.1", "W=%s" % hpc,
           'cd "$W" || exit 1', 'if ! mkdir .running 2>/dev/null; then echo "FAIL: another job holds .running"; exit 1; fi',
           "trap 'rmdir \"$W/.running\" 2>/dev/null' EXIT", 'echo "host=$(hostname) start=$(date)"',
           "run() {",
           '  d="$W/$1"',
           '  if [ -s "$d/job.out" ] && grep -q "ORCA TERMINATED NORMALLY" "$d/job.out"; then echo "$1 done before"',
           '  else ( cd "$d" && "$ORCA/orca" job.inp > job.out 2>&1 </dev/null ); echo "$1 orca_exit=$? normal=$(grep -c \'ORCA TERMINATED NORMALLY\' "$d/job.out") $(date)"; fi',
           '  if [ "$2" != "-" ] && [ ! -s "$d/vpot.out" ] && grep -q "ORCA TERMINATED NORMALLY" "$d/job.out"; then',
           '    ( cd "$d" && "$ORCA/orca_vpot" job.gbw job.scfp "$W/$2" vpot.out job > vpot.log 2>&1 ); echo "$1 vpot_exit=$? lines=$(wc -l < "$d/vpot.out")"',
           "  fi",
           "}"]
    pbs += ["run %s %s" % (r[0], r[6]) for r in runs] + ['echo "end=$(date)"']
    (out / "stage5_a.pbs").write_text("\n".join(pbs) + "\n")
    files = sorted(p for p in out.rglob("*") if p.is_file())
    (out / "INPUT_CHECKSUMS.sha256").write_text("".join("%s  %s\n" % (sha(p), p.relative_to(out)) for p in files))
    print("prepared %d runs in %s: %d frames x R, TS, P; invariance on frame %s (rotation %.1f deg); finite differences at frame %s,"
          " sites %s (nearest %s at R), +-%.2f e" % (len(runs), out, len(frames), fi, rot_angle(Ri), f0,
                                                    " ".join(map(str, fd_sites)), ", ".join(FD_ATOMS), DELTA))

# ------------------------------------------------------------------ assemble
def orca_run(d):
    """(energy Eh, coordinates echoed) after the checks of (a); a string on failure"""
    p = Path(d) / "job.out"
    if not p.is_file(): return "no job.out"
    t = p.read_text(errors="replace")
    if "****ORCA TERMINATED NORMALLY****" not in t: return "ORCA did not terminate normally"
    if not re.search(r"SCF CONVERGED AFTER\s+\d+\s+CYCLES", t) or "SCF NOT CONVERGED" in t: return "SCF not converged"
    m = re.search(r"Number of Electrons\s+NEL\s+\.\.\.\.\s+(\d+)", t)
    if not m or int(m.group(1)) != 118: return "not 118 electrons"
    e = re.findall(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", t)
    if len(e) != 1: return "%d final energies" % len(e)
    k = t.find("CARTESIAN COORDINATES (ANGSTROEM)")
    X = np.array([[float(v) for v in l.split()[1:4]] for l in t[k:].splitlines()[2:26]])
    _, S = G.read_xyz(Path(d) / "substrate.xyz")
    if X.shape != (24, 3) or np.abs(X - S).max() > 1e-5: return "echoed coordinates differ from substrate.xyz"
    if (Path(d) / "pc.pc").is_file() and "Reading point charge file                   ... ok (1 point charges)" not in t:
        return "the point charge was not read"
    return float(e[0])

def vpot(d, pts):
    """potentials (Eh/e) in point order, after checking the echoed points"""
    P = np.array([[float(v) for v in l.split()] for l in Path(pts).read_text().splitlines()[1:] if l.strip()])
    p = Path(d) / "vpot.out"
    if not p.is_file(): return "no vpot.out"
    L = [l.split() for l in p.read_text().splitlines() if l.strip()]
    if not L or len(L[0]) != 1 or int(L[0][0]) != len(P): return "vpot.out count line"
    V = np.array([[float(v) for v in l] for l in L[1:]])
    if V.shape != (len(P), 4) or np.abs(V[:, :3] - P).max() > 1e-6: return "vpot.out points differ from %s" % Path(pts).name
    return V[:, 3]

def assemble(adir, sites_p, work, out):
    if os.path.exists(out): stop("%s exists" % out)
    work = Path(work)
    frames, rows, Q = frames_and_sites(adir, sites_p)
    bad, rep = [], []
    for l in (work / "INPUT_CHECKSUMS.sha256").read_text().splitlines():          # (d)
        h, f = l.split("  ", 1)
        if not (work / f).is_file() or sha(work / f) != h: bad.append("(d) %s changed or missing since prepare" % f)
    runs = [l.split("\t") for l in (work / "runs.tsv").read_text().splitlines()[1:] if l.strip()]
    want, fi, f0, fd_sites = plan(adir, frames, Q)
    if [tuple(r) for r in runs] != want: stop("runs.tsv differs from the plan for these geometries and sites")
    ro = [l.split("\t") for l in (work / "row_order.tsv").read_text().splitlines()[2:] if l.strip()]
    if [(r[1], r[3], r[4], r[5]) for r in ro] != [(r["site"], r["x"], r["y"], r["z"]) for r in rows]: stop("row_order.tsv differs from the sites file")
    E, V = {}, {}
    for name, kind, f, s, site, q, pts in runs:
        src = Path(adir) / ("frame_" + f) / G.STATES[s]
        if kind != "inv" and (work / name / "substrate.xyz").read_bytes() != src.read_bytes():
            bad.append("(d) %s: substrate.xyz is not %s" % (name, src)); continue
        if kind == "inv":
            Ri, ti = transform(adir, f); _, X = G.read_xyz(src); _, Y = G.read_xyz(work / name / "substrate.xyz")
            if np.abs((X - ti) @ Ri.T - Y).max() > 1e-8: bad.append("(d) %s: substrate.xyz is not frame %s's R in its own coordinates" % (name, f)); continue
        e = orca_run(work / name)
        if isinstance(e, str): bad.append("(a) %s: %s" % (name, e)); continue
        E[name] = e
        if pts != "-":
            v = vpot(work / name, work / pts)
            if isinstance(v, str): bad.append("(a) %s: %s" % (name, v)); continue
            V[name] = v
    rep.append("(a) %d runs, %d energies and %d potential sets read; %d problems" % (len(runs), len(E), len(V), sum(1 for b in bad if b.startswith("(a)"))))
    if not bad:
        A = np.stack([(V["sp_%s_TS" % f] - V["sp_%s_R" % f]) * H for f in frames], 1)
        AP = np.stack([(V["sp_%s_P" % f] - V["sp_%s_R" % f]) * H for f in frames], 1)
        dv = np.abs(V["inv_%s_R" % fi] - V["sp_%s_R" % fi]).max(); de = abs(E["inv_%s_R" % fi] - E["sp_%s_R" % fi])
        ok_b = dv <= 1e-4
        rep.append("(b) frame %s in its own coordinates (rotation %.1f deg): max |dV| %.2e Eh/e (limit 1e-4), |dE| %.2e Eh -> %s"
                   % (fi, rot_angle(transform(adir, fi)[0]), dv, de, "PASS" if ok_b else "FAIL"))
        if not ok_b: bad.append("(b) invariance")
        k0 = frames.index(f0)
        for i in fd_sites:
            g = lambda s, t: E["fd_%s_%s_%d_%s" % (f0, s, i, t)]
            afd = (g("TS", "p") - g("TS", "m") - g("R", "p") + g("R", "m")) / (2 * DELTA) * H
            a = A[i, k0]; tol = 0.01 + 0.01 * abs(a); ok = abs(afd - a) <= tol
            rep.append("(c) frame %s site %d (%s layer): a %+.4f, a_fd %+.4f, diff %+.4f (limit %.4f) kcal/mol per e -> %s"
                       % (f0, i, rows[i]["layer"], a, afd, afd - a, tol, "PASS" if ok else "FAIL"))
            if not ok: bad.append("(c) site %d" % i)
    os.makedirs(out)
    verdict = "PASS" if not bad else "FAIL"
    head = ["# Stage 5 A matrix, cs5_amatrix.py assemble (STAGE5_CRITERIA.txt); cs5_amatrix.py sha256 %s" % sha(__file__),
            "# sites %s sha256 %s (row i = site i); aligned geometries %s; work %s" % (sites_p, sha(sites_p), adir, work),
            "# frames (%d): %s" % (len(frames), " ".join(frames))]
    en = []
    if all("sp_%s_%s" % (f, s) in E for f in frames for s in STATES):
        en = ["", "in vacuo energies (Eh) and barrier E(TS) - E(R) (kcal/mol), B3LYP-D3BJ/def2-SVP at the aligned geometries:",
              "frame\tE_R\tE_TS\tE_P\tbarrier_invacuo"] + [
              "%s\t%.9f\t%.9f\t%.9f\t%.4f" % (f, E["sp_%s_R" % f], E["sp_%s_TS" % f], E["sp_%s_P" % f], (E["sp_%s_TS" % f] - E["sp_%s_R" % f]) * H)
              for f in frames]
    Path(out, "A_REPORT.txt").write_text("\n".join(head + ["# checks: " + verdict] + rep + ["PROBLEM: " + b for b in bad] + en) + "\n")
    print("\n".join(rep + ["PROBLEM: " + b for b in bad] + ["checks " + verdict]))
    if bad: sys.exit(1)
    for name, M, what in (("A.tsv", A, "a_if = V_TS,f(r_i) - V_R,f(r_i)"), ("AP.tsv", AP, "ap_if = V_P,f(r_i) - V_R,f(r_i) (not used by Stage 5)")):
        with open(Path(out) / name, "x") as fh:
            fh.write("\n".join(head + ["# %s, kcal/mol per e; in vacuo B3LYP-D3BJ/def2-SVP densities at the aligned geometries" % what,
                                       "# a < 0: a positive charge at site i lowers frame f's barrier (first order); ddE_f = sum_i a_if q_i",
                                       "# range %.4f to %.4f" % (M.min(), M.max())]) + "\n")
            for r in M: fh.write("\t".join("%.6f" % v for v in r) + "\n")
    print("wrote %s/A.tsv and AP.tsv: %d sites x %d frames" % (out, A.shape[0], A.shape[1]))

if __name__ == "__main__":
    a = sys.argv[1:]
    if a[:1] == ["prepare"] and len(a) in (4, 5): prepare(*a[1:])
    elif a[:1] == ["assemble"] and len(a) == 5: assemble(*a[1:])
    else: stop(__doc__.split("\n\n")[1])
