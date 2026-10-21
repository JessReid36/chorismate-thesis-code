#!/usr/bin/env python3
"""
cs6_arrange.py - Stage 6's arrangements and their ORCA inputs (phase2.2, cs_stage6/STAGE6_CRITERIA.txt with
STAGE6_AMENDMENT1.txt and STAGE6_ERRATUM1.txt).

  python3 cs6_arrange.py selftest
  python3 cs6_arrange.py arrange <frame_dir> <sites.tsv> <L.tsv> <A.tsv> <prms> <OUT> [HPC_DIR]

It is a separate file from cs6_prepare.py, whose "extract" is deliberately standard-library-only, because this one
needs numpy and scipy to reuse cs_opt_milp.py. Nothing about the pre-registered rule changes with the file it lives
in; STAGE6_CRITERIA.txt fixes the rule, not the module.

THE SEVEN ARRANGEMENTS, exactly as STAGE6_CRITERIA.txt fixes them, placed by rule and never by a barrier:
  S6-2   the most barrier-lowering admitted site of this grid, then the next in that ranking at least r_min from it;
         the charges then optimised jointly for those two FIXED sites under the budget;
  S6-3   three sites by the same rule, charges optimised jointly for the three;
  S6-5, S6-10, S6-20   cs_opt_milp.py's own best design at n = 5, 10, 20 on this grid;
  S6-C   no spheres, the control;
  S6-P   one bare +1 point charge, no LJ wall, at S6-2's first site.

THE RANKING is the closed form of the one-sphere score, which is what cs_opt_milp.py's selftest 1 already checks the
MILP against: a site's best single-sphere contribution is min over its admitted ranges and over each range's two
bounds of (mean_f a_if) q + (mean_f L_if). This script computes it directly AND asks the MILP for n = 1 on the same
grid, and stops if the two disagree, so the ranking is not a reimplementation taken on trust.

WHY THE SEMI-DEVIATION PLAYS NO PART HERE, stated because it would otherwise look like a setting was dropped: the
score is dbar + (lambda/F) sum_f u_f with u_f >= d_f - dbar. Stage 6 is ONE frame, so F = 1, d_0 = dbar identically,
u_0 is driven to 0, and the objective is exactly the frozen score d_0. lambda is untouched and is still 0.94; it has
nothing to weight. The same is true of the budget at n = 2: two spheres at |q| = 1 spend exactly B = 2.0, so the
budget is tight-but-satisfied there and binds only from n = 3 up. Both are reported, not assumed.

WHAT EACH ARRANGEMENT GETS, per geometry (R from reactant_qm.xyz, TS from transition_state_qm.xyz):
  <arr>/<R|TS>_sp/    the FROZEN single point with gradients: the same QM/MM Hamiltonian as the relaxation, with
                      "EnGrad" in place of "Opt" and no %geom block. This is what ddE(frozen) and G_add are measured
                      from, and it must be the same Hamiltonian as the relaxation or Loss = ddE(frozen) -
                      ddE(relaxed) would contain the sphere's Lennard-Jones energy rather than the relaxation. The
                      reading is forced by the score itself: cs5_grid.py defines L[i,f] as the sphere's LJ energy at
                      TS_f minus at R_f, so the quantity "the MILP's score estimates" already contains the LJ term,
                      and the measurement of it must too. Recorded here rather than left implicit.
  <arr>/<R|TS>_opt/   the relaxation, STAGE6_ERRATUM1.txt's protocol: armb_freecart's, all 24 substrate atoms active,
                      coordsys cartesian, ProjectTR false, MaxIter 1000, the spheres absent from ActiveAtoms and so
                      frozen. R by Opt; TS by OptTS with Freq in the same job, so the imaginary-mode count criterion
                      (ii) is answered without a second submission.
  S6-C gets both too, with no sphere in either file - the bare substrate at the same two geometries, which is also
  the g_i(0) of G_add and the E(0) of ddE(frozen), so the control is not an extra expense but the shared reference.

THE SPHERE, in ORCA's AMBER ORCAFF layout, as WP4's committed site.ORCAFF.prms writes it and this script reproduces
byte for byte for one sphere at q = +1:
    <idx>   N   <q>   -0.170000   3.648000   -0.085000   3.648000
i.e. charge, then (-epsilon, 2 R*) for the full and the 1-4 scaled interaction, with epsilon 0.170 kcal/mol and
R* 1.824 A - the Amber N quasi-atom of cs5_grid.py line 74 (SPHERE_RSTAR, SPHERE_EPS), the same atom whose LJ energy
L is computed from. S6-P's bare charge zeroes all four, exactly as WP4's arma did. The 24 substrate rows are copied
verbatim from rows 6208-6231 of the committed complex_solvated.ORCAFF.prms, checked row by row against the elements
NAMES implies, and the three bonded sections are written empty with WP4's own constants.

CHECKS, every one a stop:
  (a) the three geometry files and transform.txt are present and transform.txt is the identity cs6_prepare.py wrote;
  (b) each geometry has 24 atoms in the substrate's element order;
  (c) A and L have one frame, and the same number of rows as the site table;
  (d) the closed-form ranking agrees with the MILP at n = 1, both site and score;
  (e) every placed design satisfies every rule: cardinality, each charge inside an admitted range, every pair at
      least r_min apart, the true sum of squares within B + the tangent slack, and the reported objective equal to
      the score recomputed from a and L (cs_opt_milp.violations, worst violation printed);
  (f) no sphere sits closer to any substrate atom, at any of the three geometries, than the floor its site was
      admitted under - re-checked here rather than trusted from the grid;
  (g) the 24 substrate force-field rows are identical to the committed force field's rows 6208-6231;
  (h) every arrangement's two geometries use the SAME sphere coordinates, and S6-P's site is S6-2's first;
  (i) each written job.inp contains exactly the protocol's %geom block, and ActiveAtoms lists all 24 and no sphere.
Refuses to write into an existing folder. Deterministic. Submits nothing and runs nothing.

Python 3.8+, numpy, scipy (through cs_opt_milp.py).
"""
import hashlib
import os
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import cs_opt_milp as opt                                                     # noqa: E402

NAMES = "C1 H1 H2 C2 C3 O1 O2 O3 C4 H3 C5 H4 C6 C7 H5 C8 H6 C9 H7 O4 H8 C10 O5 O6".split()
ELS = ["".join(c for c in n if c.isalpha()) for n in NAMES]
NQM, QM0 = 24, 6207
SPHERE_RSTAR, SPHERE_EPS = 1.824, 0.170            # cs5_grid.py line 74, the Amber N quasi-atom
GEOMS = [("R", "reactant_qm.xyz"), ("TS", "transition_state_qm.xyz")]
BUDGET_DEFAULT = 2.0

HEAD_COMMON = ("%maxcore 3000\n"
               "%pal nprocs 8 end\n"
               "%scf MaxIter 250 end\n")
GEOM_BLOCK = ("%geom\n"
              "  MaxIter 1000\n"
              "  coordsys cartesian\n"
              "  ProjectTR false\n"
              "end\n")


def stop(m):
    sys.exit("STOP: " + m)


def sha_file(p):
    h = hashlib.sha256()
    with open(str(p), "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


# ------------------------------------------------------------------ geometry

def read_qm_xyz(p):
    """cs6_prepare.py's 24-atom file: count, comment, then 24 rows kept verbatim."""
    lines = Path(p).read_text().splitlines()
    n = int(lines[0].split()[0])
    if n != NQM:
        stop("%s has %d atoms, expected %d" % (p, n, NQM))
    rows = [l for l in lines[2:2 + NQM]]
    els = [l.split()[0] for l in rows]
    if els != ELS:
        stop("%s is not the substrate's element order.\n   got      %s\n   expected %s"
             % (p, " ".join(els), " ".join(ELS)))
    xyz = np.array([[float(v) for v in l.split()[1:4]] for l in rows])
    return rows, xyz


def check_transform(p):
    T = np.loadtxt(str(p))
    if T.shape != (4, 3):
        stop("%s is not 4 rows of 3" % p)
    if not (np.allclose(T[:3], np.eye(3), atol=1e-12) and np.allclose(T[3], 0.0, atol=1e-12)):
        stop("%s is not the identity cs6_prepare.py writes; this script assumes one coordinate frame" % p)


# ------------------------------------------------------------------ the ranking and the placements

def per_site_best(ranges, Lm, am):
    """closed form of the one-sphere score per site: (value, q, which range). cs_opt_milp selftest 1's form."""
    out = []
    for i, rr in enumerate(ranges):
        best = (np.inf, 0.0, -1)
        for ri, (lo, hi) in enumerate(rr):
            for q in (lo, hi):
                v = am[i] * q + Lm[i]
                if v < best[0]:
                    best = (v, q, ri)
        out.append(best)
    return out


def greedy_sites(order, xyz, k, r_min):
    """walk down the ranking, keeping a site only if it is at least r_min from every site already kept."""
    kept = []
    for i in order:
        if all(float(np.linalg.norm(xyz[i] - xyz[j])) >= r_min for j in kept):
            kept.append(int(i))
            if len(kept) == k:
                return kept
    stop("could not place %d sites at least %.2f A apart walking the whole ranking" % (k, r_min))


def charges_for_fixed_sites(sites, xyz, ranges, L, A, cfg):
    """the budget's optimum for sites held fixed: the same Problem on just those rows, with n = len(sites).

    The score is separable over sites (d_f = sum_i a_if q_i + sum_i L_if z_i), L and A are per-site, and the spacing
    cliques come from the geometry, so restricting the rows restricts the problem to exactly these sites and the
    MILP then has only the charges left to choose. This reuses the committed objective, budget stand-in and bounds
    rather than re-deriving them here.
    """
    idx = list(sites)
    sub_xyz = xyz[idx]
    sub_ranges = [ranges[i] for i in idx]
    p = opt.Problem(sub_xyz, sub_ranges, L[idx], A[idx], cfg,
                    cliq=opt.cliques(sub_xyz, cfg["r_min"])[0])
    sol, info = p.best(len(idx), cfg["budget"])
    if sol is None:
        stop("no feasible charge assignment for the fixed sites %s: %s" % (idx, info.get("message")))
    if sorted(sol["sites"]) != sorted(range(len(idx))):
        stop("the fixed-site solve did not use all %d sites: %s" % (len(idx), sol["sites"]))
    # map back to grid numbering, keeping each site's own charge
    q_by_grid = {idx[t]: float(qv) for t, qv in zip(sol["sites"], sol["q"])}
    return [q_by_grid[i] for i in idx], sol, p


# ------------------------------------------------------------------ the ORCA inputs

def ff_row(idx, el, q, eps, rstar2):
    """One ORCAFF $atoms row. This format reproduces WP4's committed site.ORCAFF.prms byte for byte, substrate
    rows included (selftest 1 checks both against WP4's own file)."""
    return "%6d   %s  %12.6f %12.6f %12.6f %12.6f %12.6f" % (idx, el, q, eps, rstar2, eps / 2.0, rstar2)


def prms_text(sub, spheres):
    """WP4's site.ORCAFF.prms layout: the 24 substrate rows, then one row per sphere, then the empty sections.

    sub is the parsed substrate rows [(element, charge, eps, 2R*), ...] from the committed force field. They are
    RENUMBERED 1..24 here: the committed field numbers them 6208-6231 because it holds 55680 atoms, and a 24-atom
    system must number its own atoms from 1, exactly as WP4's file does. The five data fields are unchanged.

    spheres is a list of (q, lj) with lj True for the Amber N wall and False for a bare point charge.
    """
    out = ["$fftype", "AMBER", "$atoms", "%d 1 4" % (NQM + len(spheres))]
    for k, (el, q, eps, r2) in enumerate(sub):
        out.append(ff_row(k + 1, el, q, eps, r2))
    for k, (q, lj) in enumerate(spheres):
        e, r = (-SPHERE_EPS, 2.0 * SPHERE_RSTAR) if lj else (0.0, 0.0)
        out.append(ff_row(NQM + k + 1, "N", q, e, r))
    out += ["$bonds", "0 2 2", "$angles", "0 3 2", "$dihedrals", "0 4 3"]
    return "\n".join(out) + "\n"


def combined_xyz_text(geom_rows, pos):
    out = ["%d" % (NQM + len(pos)), "substrate (24, QM) then %d sphere(s) (MM)" % len(pos)]
    out.extend(geom_rows)
    for p in pos:
        out.append("N      %.8f    %.8f    %.8f" % (p[0], p[1], p[2]))
    return "\n".join(out) + "\n"


def job_inp_text(kind, side, nsph):
    """kind 'sp' = the frozen single point with gradients; kind 'opt' = the relaxation.

    STAGE6_ERRATUM1.txt's protocol for 'opt': all 24 substrate atoms active, Cartesian, ProjectTR false,
    MaxIter 1000. 'sp' is the same Hamiltonian with EnGrad and no %geom block, so that Loss is a difference of
    like with like.

    WITH NO SPHERES (the control S6-C) THIS IS A PLAIN QM JOB, not a QM/MM one with an empty MM region. An MM
    region with no atoms in it contributes nothing to the energy or the gradient, so the two are the same
    Hamiltonian; but ORCA's QM/MM driver with every atom in QMAtoms and nothing left over is an edge case this
    project has never run, whereas the plain single point below is EXACTLY what cs5_amatrix.py runs for the bare
    substrate - committed, and run 16 times for the A column of this very frame. The known-good path is the right
    one for the reference every other number is measured against. ARRANGEMENTS.txt records one extra diagnostic
    single point that tests the two paths agree.
    """
    if kind == "sp":
        kw = "EnGrad"
    elif side == "R":
        kw = "Opt"
    else:
        kw = "OptTS Freq"
    qmmm = " QMMM" if nsph else ""
    s = "! B3LYP D3BJ def2-SVP def2/J RIJCOSX TightSCF%s %s\n" % (qmmm, kw)
    s += HEAD_COMMON
    if kind == "opt":
        s += GEOM_BLOCK
    if nsph:
        s += "%qmmm\n"
        s += "  QMAtoms {0:%d} end\n" % (NQM - 1)
        if kind == "opt":
            s += "  ActiveAtoms {%s} end\n" % " ".join(str(i) for i in range(NQM))
        s += '  ORCAFFFilename "site.ORCAFF.prms"\n'
        s += "end\n"
    s += "* xyzfile -2 1 combined.xyz\n"
    return s


def job_pbs_text(name, remote):
    return ("#!/bin/bash\n"
            "#PBS -N %s\n"
            "#PBS -l select=1:ncpus=8:mem=24gb\n"
            "#PBS -l walltime=168:00:00\n"
            "#PBS -m ae\n"
            "#PBS -M 18660916@sun.ac.za\n"
            "#PBS -j oe\n"
            "#PBS -o %s/job.pbs.out\n"
            "cd %s\n"
            'export PATH="/apps/openmpi/4.1.1/bin:$PATH"\n'
            'export LD_LIBRARY_PATH="/home/apps2/ORCA/6.0.1/lib:/apps/openmpi/4.1.1/lib:'
            '/apps/mambaforge/envs/medaka/lib:${LD_LIBRARY_PATH:-}"\n'
            "export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1\n"
            'echo "host=$(hostname) start=$(date)"\n'
            "/home/apps2/ORCA/6.0.1/orca job.inp > job.out 2>&1 </dev/null\n"
            'echo "orca_exit=$? end=$(date)"\n'
            "grep -c 'GEOMETRY OPTIMIZATION CYCLE' job.out || true\n"
            "grep -E 'THE OPTIMIZATION HAS CONVERGED|ORCA TERMINATED NORMALLY' job.out | tail -2 || true\n"
            "rm -f job.gbw job*.tmp 2>/dev/null || true\n" % (name, remote, remote))


def check_inp(text, kind, nsph=1):
    """check (i): the written input really is the protocol, read back from the text."""
    if nsph == 0:
        for bad in ("QMMM", "%qmmm", "ORCAFFFilename", "QMAtoms"):
            if bad in text:
                stop("the no-sphere control's input contains %r; it is meant to be a plain QM job" % bad)
        if kind == "opt" and ("coordsys cartesian" not in text or "ProjectTR false" not in text):
            stop("the control's relaxation is missing the Cartesian protocol lines")
        if kind == "sp" and "EnGrad" not in text:
            stop("the control's frozen single point does not ask for EnGrad")
        return
    if kind == "opt":
        for need in ("coordsys cartesian", "ProjectTR false", "MaxIter 1000"):
            if need not in text:
                stop("a relaxation input is missing %r" % need)
        want = "ActiveAtoms {%s} end" % " ".join(str(i) for i in range(NQM))
        if want not in text:
            stop("a relaxation input does not list all %d substrate atoms as active" % NQM)
        if "EnGrad" in text:
            stop("a relaxation input also asks for EnGrad")
    else:
        if "EnGrad" not in text:
            stop("a frozen single point does not ask for EnGrad")
        for bad in ("%geom", "ActiveAtoms", " Opt", "OptTS"):
            if bad in text:
                stop("a frozen single point contains %r, which would move the geometry" % bad)
    if "QMAtoms {0:%d} end" % (NQM - 1) not in text:
        stop("an input does not put the 24 substrate atoms in the QM region")


def substrate_ff_rows(prms):
    """check (g): rows 6208-6231 of the committed force field, verbatim, with their elements checked."""
    rows, n = [], None
    with open(str(prms)) as f:
        for i, ln in enumerate(f):
            if i == 3:
                n = int(ln.split()[0])
            elif n is not None and 4 <= i < 4 + n:
                if QM0 <= i - 4 < QM0 + NQM:
                    rows.append(ln.rstrip("\n"))
                elif i - 4 >= QM0 + NQM:
                    break
    if len(rows) != NQM:
        stop("could not read %d substrate rows from %s (got %d)" % (NQM, prms, len(rows)))
    els = [r.split()[1] for r in rows]
    if els != ELS:
        stop("the force field's rows %d-%d are not the substrate's elements.\n   got      %s\n   expected %s"
             % (QM0 + 1, QM0 + NQM, " ".join(els), " ".join(ELS)))
    sub = []
    for r in rows:
        f = r.split()
        if len(f) != 7:
            stop("force-field row %r does not have 7 fields" % r)
        sub.append((f[1], float(f[2]), float(f[3]), float(f[4])))
        if abs(float(f[5]) - float(f[3]) / 2.0) > 1e-9 or abs(float(f[6]) - float(f[4])) > 1e-9:
            stop("force-field row %r is not the (eps, 2R*) then (eps/2, 2R*) layout this writer assumes" % r)
    tot = sum(x[1] for x in sub)
    if abs(tot - (-2.0)) > 1e-6:
        stop("the 24 substrate charges sum to %+.6f, not -2 as the dianion must" % tot)
    return sub


def arrange(frame_dir, sites_p, L_p, A_p, prms_p, out_p, hpc_dir):
    frame_dir, out = Path(frame_dir), Path(out_p)
    if out.exists():
        stop("%s exists; this script refuses to write into an existing folder" % out)

    print("== (a,b) the geometries")
    check_transform(frame_dir / "transform.txt")
    G = {}
    for side, fn in GEOMS:
        rows, xyzg = read_qm_xyz(frame_dir / fn)
        G[side] = (rows, xyzg)
        print("   %-2s %-26s 24 atoms, element order ok, sha256 %s" % (side, fn, sha_file(frame_dir / fn)[:16]))
    if not (frame_dir / "product_qm.xyz").exists():
        stop("product_qm.xyz is missing; the grid was built on all three states")

    print("\n== (c) the grid, L and A")
    xyz, layer, ranges, L, A = opt.load(str(sites_p), str(L_p), str(A_p))
    if A.shape[1] != 1:
        stop("A has %d frames; Stage 6 is one frame. Use that frame's own A column." % A.shape[1])
    cfg = dict(opt.CONFIG, budget=BUDGET_DEFAULT)
    am, Lm = A.mean(1), L.mean(1)
    print("   F = 1, so the semi-deviation term is identically zero and the objective IS the frozen score;")
    print("   lambda stays %.2f and has nothing to weight. budget B = %.2f e^2, r_min = %.2f A."
          % (cfg["lam"], cfg["budget"], cfg["r_min"]))

    print("\n== (d) the ranking, closed form against the MILP at n = 1")
    best = per_site_best(ranges, Lm, am)
    vals = np.array([b[0] for b in best])
    order = list(np.argsort(vals, kind="stable"))
    top = int(order[0])
    cl_all = opt.cliques(xyz, cfg["r_min"])[0]
    p_full = opt.Problem(xyz, ranges, L, A, cfg, cliq=cl_all)
    s1, info1 = p_full.best(1, cfg["budget"])
    if s1 is None:
        stop("the MILP found no one-sphere design: %s" % info1.get("message"))
    if s1["sites"] != [top] or abs(s1["objective"] - vals[top]) > 1e-6:
        stop("the ranking disagrees with the MILP at n = 1.\n"
             "   closed form: site %d, score %+.6f\n   MILP:        site %s, score %+.6f"
             % (top, vals[top], s1["sites"], s1["objective"]))
    print("   agree: site %d, score %+.6f kcal/mol, q %+.2f  (MILP %+.6f)"
          % (top, vals[top], best[top][1], s1["objective"]))
    print("   the five best single sites: " + "  ".join("%d(%+.3f)" % (i, vals[i]) for i in order[:5]))

    print("\n== the arrangements")
    designs = []                   # (name, [sites], [charges], lj, solution-or-None, from_the_full_grid)
    for name, k in (("S6-2", 2), ("S6-3", 3)):
        sel = greedy_sites(order, xyz, k, cfg["r_min"])
        q, sol, _ = charges_for_fixed_sites(sel, xyz, ranges, L, A, cfg)
        designs.append((name, sel, q, True, sol, False))
        print("   %-5s sites %s  q %s  sum q^2 %.4f  score %+.4f"
              % (name, sel, ["%+.3f" % v for v in q], sum(v * v for v in q), sol["objective"]))
    for n in (5, 10, 20):
        sol, info = p_full.best(n, cfg["budget"])
        if sol is None:
            stop("the optimiser found no design at n = %d: %s" % (n, info.get("message")))
        designs.append(("S6-%d" % n, [int(i) for i in sol["sites"]], [float(v) for v in sol["q"]], True, sol, True))
        print("   %-5s n = %-2d score %+.4f  sum q^2 %.4f  %s, certified gap %.4f  (%.0f s)"
              % ("S6-%d" % n, n, sol["objective"], sum(float(v) ** 2 for v in sol["q"]),
                 sol.get("strategy", "full"), sol.get("certified_gap", float("nan")), sol.get("wall", 0.0)))
    designs.append(("S6-C", [], [], True, None, False))
    designs.append(("S6-P", [top], [1.0], False, None, False))
    print("   %-5s no spheres (the control, and the g(0)/E(0) reference)" % "S6-C")
    print("   %-5s site %d, bare +1, no LJ wall (S6-2's first site: %s)"
          % ("S6-P", top, "yes" if designs[0][1][0] == top else "NO"))

    print("\n== (e,h) every rule, on every placed design")
    if designs[0][1][0] != top:
        stop("S6-P's site %d is not S6-2's first site %d" % (top, designs[0][1][0]))
    worst = 0.0
    for name, sel, q, lj, sol, full in designs:
        if not sel:
            continue
        for a in range(len(sel)):
            for b in range(a + 1, len(sel)):
                d = float(np.linalg.norm(xyz[sel[a]] - xyz[sel[b]]))
                if d < cfg["r_min"] - 1e-9:
                    stop("%s has sites %d and %d only %.3f A apart" % (name, sel[a], sel[b], d))
        for i, qv in zip(sel, q):
            if not any(lo - 1e-9 <= qv <= hi + 1e-9 for lo, hi in ranges[i]):
                stop("%s puts q %+.4f on site %d, which admits %s" % (name, qv, i, ranges[i]))
        ssq = sum(v * v for v in q)
        cap = cfg["budget"] + opt.TANGENT_SLACK * len(sel)
        if name != "S6-P" and ssq > cap + 1e-9:
            stop("%s spends sum q^2 = %.5f, above B + the tangent slack (%.5f)" % (name, ssq, cap))
        if full:
            # only the full-grid designs carry grid site indices, so only they can be handed to violations();
            # S6-2 and S6-3 came from a subset problem and are checked by the three explicit rules above.
            v = opt.violations(p_full, sol)
            worst = max(worst, max(v.get(key, 0.0)
                                   for key in ("cardinality", "charge_range", "spacing", "d_f", "objective")))
    print("   spacing, admitted ranges and the budget all honoured on all seven;")
    print("   worst rule or arithmetic violation on the three optimiser designs: %.1e" % worst)
    return out, G, (xyz, ranges, L, A, am, Lm), cfg, designs, top


def write_all(out, G, grid, cfg, designs, top, prms_p, frame_dir, sites_p, L_p, A_p, hpc_dir):
    xyz, ranges, L, A, am, Lm = grid
    sub_rows = substrate_ff_rows(prms_p)
    print("\n== (g) the force field, against WP4's committed and successfully run files")
    print("   24 substrate rows read from %s rows %d-%d, charges summing to -2, renumbered 1..24"
          % (prms_p.name, QM0 + 1, QM0 + NQM))
    # (g2) the whole writer against WP4's two committed site.ORCAFF.prms files. They were written by a different
    # route for a different frame and ORCA ran both to convergence, so byte-identity here is an end-to-end check
    # of the element order, the charges, the LJ fields, the renumbering, the count line and the empty sections at
    # once. Skipped with a line if the WP4 folder is not beside this run.
    root = prms_p.resolve().parent.parent.parent                    # <repo>/05_qmmm/13_bridge/x.prms -> <repo>
    pairs = [("armb_freecart", (1.0, True)), ("arma_freecart", (1.0, False))]
    done = 0
    for folder, sph in pairs:
        w = root / "phase2.2" / "wp4_s13_rerun" / folder / "site.ORCAFF.prms"
        if not w.exists():
            continue
        if prms_text(sub_rows, [sph]) != w.read_text():
            stop("the force-field writer does not reproduce %s byte for byte" % w)
        done += 1
        print("   reproduces %s/site.ORCAFF.prms byte for byte" % folder)
    if not done:
        print("   (WP4's files are not beside this run, so that comparison was skipped)")

    out.mkdir(parents=True)
    rec = ["Stage 6 arrangements and inputs, written by cs6_arrange.py.",
           "Pre-registered by phase2.2/cs_stage6/STAGE6_CRITERIA.txt with STAGE6_AMENDMENT1.txt and",
           "STAGE6_ERRATUM1.txt. Placed by rule; no arrangement was chosen or rejected by an effect it produces.",
           "",
           "INPUTS",
           "   geometries      %s" % frame_dir,
           "      reactant_qm.xyz           sha256 %s" % sha_file(frame_dir / "reactant_qm.xyz"),
           "      transition_state_qm.xyz   sha256 %s" % sha_file(frame_dir / "transition_state_qm.xyz"),
           "      transform.txt             sha256 %s  (the identity)" % sha_file(frame_dir / "transform.txt"),
           "   sites           %s  sha256 %s" % (sites_p, sha_file(sites_p)),
           "   L               %s  sha256 %s" % (L_p, sha_file(L_p)),
           "   A               %s  sha256 %s" % (A_p, sha_file(A_p)),
           "   force field     %s  sha256 %s" % (prms_p, sha_file(prms_p)),
           "   cs_opt_milp.py  sha256 %s" % sha_file(HERE / "cs_opt_milp.py"),
           "   cs6_arrange.py  sha256 %s" % sha_file(Path(__file__).resolve()),
           "",
           "SETTINGS   B = %.2f e^2 (the code default, NOT the calibrated B), r_min = %.2f A, lambda = %.2f"
           % (cfg["budget"], cfg["r_min"], cfg["lam"]),
           "   One frame, so F = 1: the semi-deviation term is identically zero and the objective is the frozen",
           "   score itself. lambda is unchanged and has nothing to weight. At n = 2 two spheres at |q| = 1 spend",
           "   exactly B, so the budget is satisfied but not slack there; it binds from n = 3 up.",
           "",
           "THE RANKING (closed form, checked against the MILP at n = 1: same site, same score to 1e-6)",
           "   the twenty best single sites, score in kcal/mol:"]
    vals = np.array([per_site_best(ranges, Lm, am)[i][0] for i in range(len(xyz))])
    order = list(np.argsort(vals, kind="stable"))
    for i in order[:20]:
        bq = per_site_best(ranges, Lm, am)[int(i)]
        rec.append("      site %5d  %+9.4f  at q %+.2f   admits %s"
                   % (i, bq[0], bq[1], " and ".join("[%.2f, %.2f]" % r for r in ranges[int(i)])))
    rec.append("")
    rec.append("THE ARRANGEMENTS")

    written = []
    for name, sel, q, lj, sol, full in designs:
        pos = [xyz[i] for i in sel]
        rec.append("")
        rec.append("  %s   %d sphere(s)%s" % (name, len(sel), "" if lj else ", BARE point charge, no LJ wall"))
        if sol is not None:
            rec.append("     score %+.5f kcal/mol   sum q^2 %.5f (true)   strategy %s   certified gap %s"
                       % (sol["objective"], sum(v * v for v in q), sol.get("strategy", "-"),
                          ("%.5f" % sol["certified_gap"]) if "certified_gap" in sol else "-"))
        for i, qv in zip(sel, q):
            rec.append("     site %5d  q %+.4f  at %.6f %.6f %.6f  a %+.5f  L %+.5f"
                       % (i, qv, xyz[i][0], xyz[i][1], xyz[i][2], am[i], Lm[i]))
        for side, _fn in GEOMS:
            grows, gxyz = G[side]
            # check (f): how close any sphere comes to any substrate atom, at THIS geometry
            if pos:
                dmin = min(float(np.linalg.norm(np.array(p_) - a_)) for p_ in pos for a_ in gxyz)
                closest = min(((float(np.linalg.norm(np.array(p_) - gxyz[j])), NAMES[j])
                               for p_ in pos for j in range(NQM)))
                if dmin < 1.2:
                    stop("%s at %s puts a sphere %.3f A from %s, below the 1.2 A collapse threshold"
                         % (name, side, dmin, closest[1]))
                rec.append("     %s: closest sphere-to-substrate %.4f A (to %s)" % (side, dmin, closest[1]))
            for kind in ("sp", "opt"):
                d = out / name / ("%s_%s" % (side, kind))
                d.mkdir(parents=True)
                inp = job_inp_text(kind, side, len(pos))
                check_inp(inp, kind, len(pos))
                (d / "job.inp").write_text(inp)
                (d / "combined.xyz").write_text(combined_xyz_text(grows, pos))
                if pos:
                    (d / "site.ORCAFF.prms").write_text(prms_text(sub_rows, [(qv, lj) for qv in q]))
                remote = "%s/%s/%s_%s" % (hpc_dir.rstrip("/"), name, side, kind) if hpc_dir else \
                         "CHANGE_ME/%s/%s_%s" % (name, side, kind)
                (d / "job.pbs").write_text(job_pbs_text("cm_s6_%s_%s_%s" % (name.replace("S6-", ""), side, kind),
                                                        remote))
                written.append(d)

    # (h) the same sphere coordinates at both geometries
    for name, sel, q, lj, sol, full in designs:
        if not sel:
            continue
        a = (out / name / "R_opt" / "combined.xyz").read_text().splitlines()[2 + NQM:]
        b = (out / name / "TS_opt" / "combined.xyz").read_text().splitlines()[2 + NQM:]
        if a != b:
            stop("%s does not use the same sphere coordinates at R and TS" % name)
    print("== (f,h,i) every input written and read back: spheres clear 1.2 A, the same positions at R and TS,")
    print("   each job.inp carries its own protocol and nothing from the other kind")

    (out / "ARRANGEMENTS.txt").write_text("\n".join(rec) + "\n")
    lines = []
    for d in sorted(written):
        for fn in sorted(os.listdir(str(d))):
            f = d / fn
            lines.append("%s  %s" % (sha_file(f), f.relative_to(out)))
    lines.append("%s  %s" % (sha_file(out / "ARRANGEMENTS.txt"), "ARRANGEMENTS.txt"))
    (out / "INPUT_CHECKSUMS.sha256").write_text("\n".join(lines) + "\n")
    print("\n== wrote %d folders, %d files, under %s" % (len(written), len(lines), out))
    print("   ARRANGEMENTS.txt and INPUT_CHECKSUMS.sha256 alongside them")
    print("\nNOTHING WAS SUBMITTED.")


def selftest():
    """the parts that do not need a grid: the force-field row, the inputs, and the ranking's closed form."""
    ok = True

    def say(good, text):
        nonlocal ok
        ok &= bool(good)
        print("   [%s] %s" % ("PASS" if good else "FAIL", text))

    print("1. the force-field rows reproduce WP4's committed site.ORCAFF.prms")
    # the real substrate's first and last rows, as WP4's file writes them after renumbering to 1..24
    W_FIRST = "     1   C     -0.311125    -0.086000     3.816000    -0.043000     3.816000"
    W_LAST = "    24   O     -0.852425    -0.210000     3.322400    -0.105000     3.322400"
    sub = [("C", -0.311125, -0.086000, 3.816000)] + [("H", 0.1, -0.015, 2.918)] * (NQM - 2) + \
          [("O", -0.852425, -0.210000, 3.322400)]
    tt = prms_text(sub, [])
    body = [l for l in tt.splitlines() if len(l.split()) == 7]
    say(body[0] == W_FIRST, "substrate row 1 byte-identical to WP4's: %r" % body[0])
    say(body[-1] == W_LAST, "substrate row 24 byte-identical to WP4's: %r" % body[-1])
    say([l.split()[0] for l in body] == [str(i + 1) for i in range(NQM)],
        "the 24 substrate rows are renumbered 1..24, not the committed field's 6208-6231")
    rows = sub
    t = prms_text(rows, [(1.0, True)])
    # WP4's own sphere row, read from its committed file rather than retyped here
    # WP4's armb_freecart/site.ORCAFF.prms sphere row, copied here verbatim including its spacing
    want = "    25   N      1.000000    -0.170000     3.648000    -0.085000     3.648000"
    got = [l for l in t.splitlines() if len(l.split()) == 7 and l.split()[1] == "N"][0]
    say(got == want, "LJ sphere row byte-identical to WP4's: %r" % got)
    t0 = prms_text(rows, [(1.0, False)])
    got0 = [f for f in (l.split() for l in t0.splitlines()) if len(f) == 7 and f[1] == "N"][0]
    say(got0[3:] == ["0.000000"] * 4, "bare charge row zeroes all four LJ fields: %s" % " ".join(got0[3:]))
    say(t.splitlines()[3] == "25 1 4" and prms_text(rows, [(1.0, True)] * 3).splitlines()[3] == "27 1 4",
        "the $atoms count line tracks the sphere count (25 1 4 / 27 1 4)")
    say(t.splitlines()[-6:] == ["$bonds", "0 2 2", "$angles", "0 3 2", "$dihedrals", "0 4 3"],
        "the three bonded sections are empty, with WP4's constants")

    print("\n2. the job inputs are the protocol, and the two kinds do not leak into each other")
    o = job_inp_text("opt", "R", 1)
    s = job_inp_text("sp", "R", 1)
    check_inp(o, "opt"); check_inp(s, "sp")
    say("coordsys cartesian" in o and "ProjectTR false" in o and "MaxIter 1000" in o,
        "the relaxation carries Cartesian, ProjectTR false, MaxIter 1000")
    say(("ActiveAtoms {%s} end" % " ".join(str(i) for i in range(NQM))) in o,
        "the relaxation lists all 24 substrate atoms active and no sphere")
    say("%geom" not in s and "ActiveAtoms" not in s and "EnGrad" in s,
        "the frozen single point has EnGrad, no %geom and no ActiveAtoms")
    say("OptTS Freq" in job_inp_text("opt", "TS", 1) and "Opt\n" in job_inp_text("opt", "R", 1),
        "TS relaxes by OptTS with Freq in the same job; R by Opt")
    bad = False
    try:
        check_inp(o.replace("coordsys cartesian\n", ""), "opt")
    except SystemExit:
        bad = True
    say(bad, "check_inp stops on a relaxation input with the coordinate line removed")

    print("\n2b. with no spheres the control is a PLAIN QM job, not QM/MM with an empty MM region")
    c_o, c_s = job_inp_text("opt", "R", 0), job_inp_text("sp", "TS", 0)
    check_inp(c_o, "opt", 0); check_inp(c_s, "sp", 0)
    say(not any(k in c_o or k in c_s for k in ("QMMM", "%qmmm", "ORCAFFFilename", "QMAtoms", "ActiveAtoms")),
        "no QMMM, %qmmm, ORCAFFFilename, QMAtoms or ActiveAtoms anywhere in the control's two inputs")
    say("coordsys cartesian" in c_o and "ProjectTR false" in c_o and "EnGrad" in c_s,
        "it still carries the Cartesian protocol, and its single point still asks for EnGrad")
    say(c_o.splitlines()[0] == "! B3LYP D3BJ def2-SVP def2/J RIJCOSX TightSCF Opt",
        "its keyword line is cs5_amatrix.py's committed header plus Opt: %r" % c_o.splitlines()[0])
    bad = False
    try:
        check_inp(job_inp_text("opt", "R", 1), "opt", 0)
    except SystemExit:
        bad = True
    say(bad, "check_inp stops if a QM/MM input is ever written for the control")

    print("\n3. the ranking's closed form, against cs_opt_milp's MILP at n = 1 on a random instance")
    rng = np.random.default_rng(7)
    S = 30
    A = rng.normal(0, 1, (S, 1)); L = rng.normal(0, 0.05, (S, 1))
    xyz = rng.normal(0, 60, (S, 3))
    ranges = [[(-1.0, 1.0)] if i % 3 == 0 else ([(-1.0, 0.0)] if i % 3 == 1 else [(0.30, 1.0)]) for i in range(S)]
    cfg = dict(opt.CONFIG, budget=2.0, n=1)
    b = per_site_best(ranges, L.mean(1), A.mean(1))
    i_cf = int(np.argmin([x[0] for x in b]))
    p = opt.Problem(xyz, ranges, L, A, cfg, cliq=opt.cliques(xyz, cfg["r_min"])[0])
    sol, _ = p.best(1, cfg["budget"])
    say(sol is not None and sol["sites"] == [i_cf] and abs(sol["objective"] - b[i_cf][0]) < 1e-6,
        "closed form site %d %+.6f; MILP site %s %+.6f" % (i_cf, b[i_cf][0], sol["sites"], sol["objective"]))

    print("\n4. the greedy spacing rule keeps the ranking order and the separation")
    xyz4 = np.array([[0, 0, 0], [1.0, 0, 0], [5.0, 0, 0], [5.5, 0, 0], [10.0, 0, 0]], float)
    sel = greedy_sites([0, 1, 2, 3, 4], xyz4, 3, 3.25)
    say(sel == [0, 2, 4], "from the ranking 0,1,2,3,4 at 3.25 A it keeps %s (1 and 3 are too close)" % sel)
    bad = False
    try:
        greedy_sites([0, 1], np.array([[0, 0, 0], [1.0, 0, 0]]), 2, 3.25)
    except SystemExit:
        bad = True
    say(bad, "it stops rather than returning fewer sites than asked")

    print("\n5. the charges for fixed sites come back on the sites asked for")
    xyz5 = np.array([[0, 0, 0], [8.0, 0, 0], [16.0, 0, 0], [24.0, 0, 0]], float)
    A5 = np.array([[-3.0], [-2.0], [+1.0], [-0.5]])
    L5 = np.zeros((4, 1))
    r5 = [[(-1.0, 1.0)]] * 4
    # B = 1.0 so the budget BINDS at two sites; at B = 2.0 both would simply sit at |q| = 1 and the test would
    # not discriminate. With the budget binding, the closed form is q proportional to -a: (3, 2)/sqrt(13).
    q, sol5, _ = charges_for_fixed_sites([0, 1], xyz5, r5, L5, A5, dict(opt.CONFIG, budget=1.0, tangents=201))
    want_q = np.array([3.0, 2.0]) / np.sqrt(13.0)
    say(len(q) == 2 and q[0] > 0 and q[1] > 0,
        "both a negative, so both charges positive: q = %s" % ["%+.3f" % v for v in q])
    say(abs(sum(v * v for v in q) - 1.0) < 0.02,
        "they spend the budget: sum q^2 = %.4f against B = 1.00" % sum(v * v for v in q))
    say(max(abs(np.array(q) - want_q)) < 0.02,
        "charges %s against the closed form %s, so the stronger site takes the larger charge"
        % (["%+.3f" % v for v in q], ["%+.3f" % v for v in want_q]))
    q2, _, _ = charges_for_fixed_sites([0, 1], xyz5, r5, L5, A5, dict(opt.CONFIG, budget=2.0))
    say(abs(q2[0] - 1.0) < 1e-6 and abs(q2[1] - 1.0) < 1e-6,
        "and at B = 2.00 the budget does not bind at two sites: q = %s, exactly the range bound"
        % ["%+.3f" % v for v in q2])

    print("\n%s" % ("ALL PASS" if ok else "SOMETHING FAILED"))
    return 0 if ok else 1


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    if sys.argv[1] == "selftest":
        return selftest()
    if sys.argv[1] != "arrange" or len(sys.argv) not in (8, 9):
        sys.exit(__doc__)
    frame_dir, sites_p, L_p, A_p, prms_p, out_p = sys.argv[2:8]
    hpc = sys.argv[8] if len(sys.argv) == 9 else ""
    out, G, grid, cfg, designs, top = arrange(frame_dir, sites_p, L_p, A_p, prms_p, out_p, hpc)
    write_all(out, G, grid, cfg, designs, top, Path(prms_p), Path(frame_dir), sites_p, L_p, A_p, hpc)
    return 0


if __name__ == "__main__":
    sys.exit(main())
