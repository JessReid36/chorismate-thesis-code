#!/usr/bin/env python3
"""Certified sparse charge optimiser -- Option A: NEUTRALISED substrate + path-wide force integrity.

Our optimiser (certified convex LP), not GOCAT's GA. Extends the path-wide force-integrity model
with a fixed NEUTRALISING layer that mimics the enzyme's arginine salt bridges: fixed +NEUT_Q point
charges are placed just beyond each carboxylate's oxygens (on the outside), locally cancelling the
-1 carboxylate "handles" the catalytic field would otherwise grab and tear off. These neutralisers
are NOT optimised -- they are part of every design (they enter the oracle field and the force model);
the certified optimiser then designs the catalytic charges around the neutralised substrate.

NOTE: Option B (explicit guanidinium groups instead of point charges) is the planned physical-
fidelity follow-up. This abstract point-charge version is the fast decisive test and the baseline
that B will be compared against.

Objective:  minimise  sum_i q_i * dV_i  +  lam1 * sum_i |q_i|   (catalytic charges only)
Constraints: box |q_i|<=QMAX ; net-neutral sum q_i = 0 (catalytic layer) ;
  for each geometry g in {R,TS,P} and each guarded bond, the TOTAL differential field-force
  (catalytic charges + fixed neutralisers) stays within +-FMAX:  |g.q + c_fixed| <= FMAX.

Certified LP via native HiGHS. Runs on the PC.
"""
import os
import numpy as np
import highspy

HERE    = os.path.dirname(os.path.abspath(__file__))
INP     = os.path.join(HERE, "..", "inputs")
GRIDS   = ["dv_grid_ext.tsv", "dv_grid_9to15.tsv"]
REACT   = "reactant.xyz"
GEOMS   = [("reactant.xyz", "substrate_charges.txt"),
           ("ts.xyz",       "substrate_charges_ts.txt"),
           ("product.xyz",  "substrate_charges_product.txt")]
CHARGED = True
QMAX    = 1.0
PRUNE   = 1.0e-3
NLAM    = 6
LAM_LO  = 0.05
LAM_HI  = 0.95
SHELL_MIN = 0.0
FMAX    = 0.010
USE_NEUTRALISERS = True
NEUT_Q  = 1.0
NEUT_D  = 2.0
CARBOXYLATES = [(4, 5, 6), (21, 22, 23)]
COV = {"C":0.77,"O":0.66,"N":0.70,"H":0.31}
BOND_TOL = 1.30
BOHR = 0.52917721


def load_grid():
    xyz, shell, dV, seen = [], [], [], set()
    for gf in GRIDS:
        for ln in open(os.path.join(INP, gf)).read().splitlines()[1:]:
            p = ln.split("\t")
            if len(p) < 8:
                continue
            x, y, z = float(p[1]), float(p[2]), float(p[3])
            key = (round(x, 3), round(y, 3), round(z, 3))
            if key in seen:
                continue
            seen.add(key)
            xyz.append((x, y, z)); shell.append(float(p[4])); dV.append(float(p[7]))
    return np.array(xyz), np.array(shell), np.array(dV)


def load_geom(xyzf, qf):
    L = open(os.path.join(INP, xyzf)).read().splitlines()[2:26]
    els = [l.split()[0] for l in L]
    R = np.array([[float(v) for v in l.split()[1:4]] for l in L])
    Q = np.array([float(x) for x in open(os.path.join(INP, qf)).read().replace(",", " ").split()])
    assert len(Q) == 24, "%s must have 24 values, got %d" % (qf, len(Q))
    return els, R, Q


def neutraliser_sites(R):
    sites = []
    for (cC, oA, oB) in CARBOXYLATES:
        b = (R[oA] - R[cC]) + (R[oB] - R[cC]); b = b / np.linalg.norm(b)
        sites.append(0.5 * (R[oA] + R[oB]) + NEUT_D * b)
    return np.array(sites)


def guarded_bonds(els, R):
    b = []
    for a in range(len(els)):
        for c in range(a + 1, len(els)):
            if {els[a], els[c]} == {"C", "O"} and \
               np.linalg.norm(R[a] - R[c]) < BOND_TOL * (COV[els[a]] + COV[els[c]]):
                b.append((a, c))
    return b


def force_rows(grid_xyz, R, Q, bonds, fixed_sites=None, fixed_q=None):
    Rb = R / BOHR; gb = grid_xyz / BOHR
    fs = (fixed_sites / BOHR) if fixed_sites is not None else None
    out = []
    for (a, b) in bonds:
        u = (R[b] - R[a]); u = u / np.linalg.norm(u)
        da = Rb[a] - gb; db = Rb[b] - gb
        fa = Q[a] * da / (np.linalg.norm(da, axis=1) ** 3)[:, None]
        fbb = Q[b] * db / (np.linalg.norm(db, axis=1) ** 3)[:, None]
        g = (fa - fbb) @ u
        c = 0.0
        if fs is not None:
            for s_pos, sq in zip(fs, fixed_q):
                dsa = Rb[a] - s_pos; dsb = Rb[b] - s_pos
                Fa = Q[a] * sq * dsa / (np.linalg.norm(dsa) ** 3)
                Fb = Q[b] * sq * dsb / (np.linalg.norm(dsb) ** 3)
                c += float((Fa - Fb) @ u)
        out.append((g, c))
    return out


def solve_l1(dV, lam1, qmax, charged, frows):
    N = dV.size; M = 2 * N
    cost = np.concatenate([dV + lam1, -dV + lam1]).astype(np.float64)
    h = highspy.Highs(); h.setOptionValue("output_flag", False)
    h.addVars(M, np.zeros(M), np.full(M, qmax))
    h.changeColsCost(M, np.arange(M, dtype=np.int32), cost)
    idx = np.arange(M, dtype=np.int32)
    if charged:
        h.addRow(0.0, 0.0, M, idx, np.concatenate([np.ones(N), -np.ones(N)]).astype(np.float64))
    for (g, c) in frows:
        h.addRow(-FMAX - c, FMAX - c, M, idx, np.concatenate([g, -g]).astype(np.float64))
    h.run()
    st = h.getModelStatus()
    if st != highspy.HighsModelStatus.kOptimal:
        raise RuntimeError("HiGHS: %s" % h.modelStatusToString(st))
    x = np.array(h.getSolution().col_value)
    return x[:N] - x[N:], float(h.getInfo().objective_function_value)


def write_design(tag, sub_lines, neut_sites, xyz, shell, dV, q):
    act = np.where(np.abs(q) >= PRUNE)[0]
    nneu = len(neut_sites) if neut_sites is not None else 0
    n = 24 + nneu + len(act)
    with open(os.path.join(HERE, "explore_%s_coords.xyz" % tag), "w") as f:
        f.write("%d\nexplore %s: 24 substrate + %d neutralisers + %d catalytic\n" % (n, tag, nneu, len(act)))
        for l in sub_lines:
            f.write(l.rstrip() + "\n")
        for s_pos in (neut_sites if nneu else []):
            f.write("C   %14.8f %14.8f %14.8f\n" % (s_pos[0], s_pos[1], s_pos[2]))
        for i in act:
            f.write("C   %14.8f %14.8f %14.8f\n" % (xyz[i, 0], xyz[i, 1], xyz[i, 2]))
    with open(os.path.join(HERE, "explore_%s_charges.txt" % tag), "w") as f:
        vals = ["0.0000"] * 24 + ["%.4f" % NEUT_Q] * nneu + ["%.4f" % q[i] for i in act]
        f.write(",".join(vals))
    with open(os.path.join(HERE, "explore_%s_report.tsv" % tag), "w") as f:
        f.write("grid_idx\tshell\tdV\tq\n")
        for i in act:
            f.write("%d\t%.1f\t%+.6f\t%+.4f\n" % (i, shell[i], dV[i], q[i]))
    return act


def shell_hist(shell, act):
    if len(act) == 0:
        return ""
    v, c = np.unique(shell[act].astype(int), return_counts=True)
    return " ".join("%dA:%d" % (a, b) for a, b in zip(v, c))


def main():
    xyz, shell, dV = load_grid()
    if SHELL_MIN > 0.0:
        keep = shell >= SHELL_MIN
        xyz, shell, dV = xyz[keep], shell[keep], dV[keep]
        print("SHELL_MIN=%.1f A -> %d candidate sites" % (SHELL_MIN, keep.sum()))
    sub_lines = open(os.path.join(INP, REACT)).read().splitlines()[2:26]

    els0, R0, Q0 = load_geom(*GEOMS[0])
    neut_sites = neutraliser_sites(R0) if USE_NEUTRALISERS else None
    fixed_q = np.full(len(neut_sites), NEUT_Q) if USE_NEUTRALISERS else None
    if USE_NEUTRALISERS:
        print("neutralisers: %d x +%.1f at carboxylates %s (fixed, not optimised)"
              % (len(neut_sites), NEUT_Q, [c[0] for c in CARBOXYLATES]))
        for k, (cC, oA, oB) in enumerate(CARBOXYLATES):
            dmin = min(np.linalg.norm(neut_sites[k] - R0[oA]), np.linalg.norm(neut_sites[k] - R0[oB]))
            print("  site %d at %.2f A from nearest carboxylate O" % (k, dmin))

    # show what the neutralisers do to the fragmenting force on the carboxylate bonds (reactant)
    b0 = guarded_bonds(els0, R0)
    fr_off = force_rows(xyz, R0, Q0, b0)
    fr_on  = force_rows(xyz, R0, Q0, b0, fixed_sites=neut_sites, fixed_q=fixed_q) if USE_NEUTRALISERS else fr_off
    print("  carboxylate-bond constant force (reactant), neutralisers off -> on:")
    for k, (a, b) in enumerate(b0):
        if a in [c for tup in CARBOXYLATES for c in tup] or b in [c for tup in CARBOXYLATES for c in tup]:
            print("    bond %s: %+.4e -> %+.4e" % ((a, b), fr_off[k][1], fr_on[k][1]))

    frows = []
    for xyzf, qf in GEOMS:
        els, R, Q = load_geom(xyzf, qf)
        bonds = guarded_bonds(els, R)
        frows += force_rows(xyz, R, Q, bonds, fixed_sites=neut_sites, fixed_q=fixed_q)
        print("  %-13s net %+.2f | %d guarded C-O bonds" % (xyzf, Q.sum(), len(bonds)))
    print("force-integrity at %d geometries, %d bond constraints; neutralisers %s (FMAX=%.4f)"
          % (len(GEOMS), len(frows), "ON" if USE_NEUTRALISERS else "OFF", FMAX))

    dVmax = abs(dV).max()
    lambdas = np.geomspace(LAM_LO, LAM_HI, NLAM) * dVmax
    print("grid: %d sites | shells %.1f-%.1f A | |dV|max %.5f | QMAX=%.2f\n"
          % (len(dV), shell.min(), shell.max(), dVmax, QMAX))
    print("tag  lam1/|dV|max  obj           net_q     #active  max_bond_force  shell histogram")

    summary = open(os.path.join(HERE, "explore_summary.tsv"), "w")
    summary.write("tag\tlam1_frac\tlam1\tobj\tnet_q\tn_active\tmax_bond_force\tshell_hist\n")
    base = ("m%d_ge%d" if SHELL_MIN > 0 else "s%d")
    for k, lam1 in enumerate(lambdas):
        tag = (base % (k, int(SHELL_MIN))) if SHELL_MIN > 0 else (base % k)
        q, obj = solve_l1(dV, lam1, QMAX, CHARGED, frows)
        act = write_design(tag, sub_lines, neut_sites, xyz, shell, dV, q)
        maxf = max((abs(float(g @ q) + c) for (g, c) in frows), default=0.0)
        hist = shell_hist(shell, act)
        print("%-6s %.3f         %+.5e  %+.1e  %6d   %.4e     %s"
              % (tag, lam1 / dVmax, obj, q.sum(), len(act), maxf, hist))
        summary.write("%s\t%.4f\t%.6e\t%.6e\t%.2e\t%d\t%.4e\t%s\n"
                      % (tag, lam1 / dVmax, lam1, obj, q.sum(), len(act), maxf, hist))
    summary.close()
    print("\nOption A: designs include fixed neutralisers + optimised catalytic charges.")
    print("wrote explore_* + explore_summary.tsv")


if __name__ == "__main__":
    main()
