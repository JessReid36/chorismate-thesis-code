#!/usr/bin/env python3
"""
cs5_grid.py - Stage 5 of the charged-sphere development (phase2.2, STAGE5_CRITERIA.txt): the candidate-site grid in
two layers with charge-dependent floors, and the spheres' Lennard-Jones term L for every site and frame.

  python3 cs5_grid.py grid    <aligned_dir> <prms> <floors.tsv> <out_sites.tsv> [key=value ...]
  python3 cs5_grid.py lj      <aligned_dir> <prms> <sites.tsv> <out_L.tsv>
  python3 cs5_grid.py check   <aligned_dir> <prms> <floors.tsv> <sites.tsv>
  python3 cs5_grid.py ljcheck <results repo>

aligned_dir: frame_<F>/{reactant_qm,transition_state_qm,product_qm}.xyz, the 24 substrate atoms of each frame aligned
onto one reference (cs5_align.py; s16_align_frames.py for the development frames). prms: the committed force field,
05_qmmm/13_bridge/complex_solvated.ORCAFF.prms (substrate rows 6208-6231). floors.tsv: cs5_floors.py's floors.

THE SURFACE. Dittner's thesis (2019), Eq. 2.16 (p. 18): each charge lies on the sphere of radius d_J about one atom J
and outside every other, over "a common vdW surface exposed by all atoms of all frames" (p. 10). Dittner's d_J are van
der Waals radii. Here d_J is the FLOOR of atom J for the sphere's charge (cs5_floors.py): how close a sphere of that
charge comes to rest (Stage 2's minima), or its Lennard-Jones contact (the LJ pair minimum, R*_sphere + R*_J) where it
has no resting place. A sized, charged sphere rests at a distance that depends on its charge (Behrens' thesis (2024),
section 5.2, p. 60: quasi-atoms are "locally optimized to their equilibrium distance"), so the grid has two layers, each
built as one Dittner surface over every atom of every geometry (states R, TS and P of every frame):
  ANION/NEUTRAL layer  radii = the q <= 0 floors: the LJ contact (Stage 2: no negative sphere has a resting place);
  CATION layer         radii = the q = +1 floors: no nearer to atom J than any resting +1 sphere came (its closest
                       approach where that is nearer than LJ contact - every atom but H3 and H7, since the spheres
                       rest at the oxygens, close to their neighbours - and 3.68 A from the hydroxyl hydrogen H8).
                       A floor is one atom's closest approach applied around that atom in every direction, so the
                       layer is an envelope, not a set of resting places; CHECK tests it radially (below).
Per layer: candidate points on a Fibonacci lattice (density points per A^2) on every atom's sphere in every geometry,
kept if outside every sphere of every geometry; enclosed pockets removed by a flood fill from the box boundary through
the voxels (spacing vox) that no sphere covers (a point is kept if its voxel, 1.5 voxels out along its sphere's normal,
is reached); then deterministic greedy Poisson-disk thinning at the resolution res, in a fixed order (geometry - frame
by number, then state R, TS, P -, atom, lattice index), so a rerun reproduces the grid exactly.
A cation-layer site that lies ON the anion/neutral surface (its clearance to the q <= 0 floors is zero within 1e-6 A)
is dropped when an anion/neutral site within res admits every charge it admits: there the two surfaces coincide and
the anion/neutral layer holds that place at the same resolution. Cation-layer sites inside that surface (a cation
resting nearer than LJ contact) and outside it (beside H8, which a cation does not approach as near as the q <= 0
floor) are kept, as are on-surface sites whose admitted charges no anion/neutral site within res offers.

CHARGES A SITE ADMITS: every charge whose floors the site clears at every geometry. q in [-1, 0] if it clears the
q <= 0 floors (q_min_neg -1; otherwise q_min_neg 0 and no charge <= 0 at that site); q in [q_min_pos, +1] for the
smallest positive level whose floors it clears (0 = any q > 0; 0.25, 0.50, 0.75, 1.00; "none" = no positive charge).
A charge between two levels takes the floors of the level below it (cs5_floors.py). An empty site (no sphere, no L)
is always allowed. The spacing between spheres (>= 3.25 A) is the optimiser's constraint, over both layers.

THE LJ TERM. For site i and frame f, L[i,f] = E_LJ(sphere at i, substrate at TS_f) - E_LJ(sphere at i, substrate at
R_f), kcal/mol: the Amber N sphere (R* 1.824 A, eps 0.170 kcal/mol) with the substrate's committed LJ parameters,
combined as ORCA does (Stage 1 S1.4), R_min,ij = (R_min,sphere + R_min,J)/2 on the full R_min, eps_ij = sqrt(eps_s eps_J),
E = eps_ij [(R_min,ij/r)^12 - 2 (R_min,ij/r)^6]. It does not depend on the sphere's charge and applies when a sphere
occupies the site. Sphere-sphere terms are the same at R and TS (spheres fixed) and cancel from the barrier.

CHECK re-reads a sites file and verifies, separately from the build, each site's admitted charges against the floors
at every geometry, its layer, on_atom, on_state, clearance and distance columns, the frame set and the sha256 of the
geometries, floors and force field. It also reports a physical test of the positive admitted charges with Stage 2's
rigid model (Coulomb in vacuum with the substrate's charges, plus the LJ above; wp5_rigid_lj_scan.py): at every site
that admits a positive charge, for the lowest such charge (0.05 e where any q > 0 is admitted) and for +1, the force
on the sphere along the site's outward normal (from the atom it lies on) at every geometry. An inward force means a
sphere held on that line would come at least as near as the site; the report counts the sites where it is pushed
outward and how far out its radial minimum then lies. LJCHECK compares the LJ formula with ORCA's own MM energies in
Stage 4c's L1 runs (dE_LJ(L1), STAGE4C_REPORT.txt, 3 decimals).

Python 3.6+, numpy. Parameters, the sha256 of every input and of this file are written into every output header.
"""
import hashlib, math, re, sys
from pathlib import Path
import numpy as np

NAMES = "C1 H1 H2 C2 C3 O1 O2 O3 C4 H3 C5 H4 C6 C7 H5 C8 H6 C9 H7 O4 H8 C10 O5 O6".split()
QM_FIRST, QM_LAST = 6208, 6231                       # substrate rows of the committed force field (1-based)
STATES = {"R": "reactant_qm.xyz", "TS": "transition_state_qm.xyz", "P": "product_qm.xyz"}
SPHERE_RSTAR, SPHERE_EPS = 1.824, 0.170              # Amber N (approach document; Stages 0-3)
DEFAULTS = {"states": "R,TS,P", "density": "48.0", "res": "0.5", "vox": "0.25", "margin": "6.0"}
COLS = ["neg", "pos0", "p025", "p050", "p075", "p100"]
POS = ["pos0", "p025", "p050", "p075", "p100"]
LEVEL_LO = {"pos0": 0.0, "p025": 0.25, "p050": 0.5, "p075": 0.75, "p100": 1.0}   # smallest positive q each level admits
LAYER_COL = {"cation": "p100", "anion": "neg"}
TOL = 1e-6                                           # A: clearance tolerance (floating point only)
HEAD = ["site", "layer", "x", "y", "z", "on_atom", "on_state", "clearance_A", "q_min_neg", "q_min_pos"] + ["d_" + n for n in NAMES]

def die(m): sys.stderr.write("STOP: " + m + "\n"); sys.exit(1)
def fmt(v, nd):
    """v to nd decimals, with a value that rounds to zero written without a sign: floating-point noise of order 1e-16
    gives -0.0000 on one machine and 0.0000 on another, and the file must not depend on the machine"""
    t = "%.*f" % (nd, v)
    return t[1:] if t.startswith("-") and float(t) == 0.0 else t
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def cfg_of(args):
    c = dict(DEFAULTS)
    for a in args:
        k, _, v = a.partition("=")
        if k not in c: die("unknown parameter " + k)
        c[k] = v
    return c

def read_xyz(p):
    L = Path(p).read_text().splitlines()
    n = int(L[0].split()[0])
    if n != 24: die("%s has %d atoms, not 24" % (p, n))
    el = [l.split()[0] for l in L[2:2 + n]]
    X = np.array([[float(v) for v in l.split()[1:4]] for l in L[2:2 + n]])
    return el, X

def frame_names(adir):
    """the frame_<F> folders of adir, in frame-number order"""
    f = [d.name[6:] for d in Path(adir).iterdir() if d.is_dir() and d.name.startswith("frame_")]
    if not f or not all(x.isdigit() for x in f): die("%s: frame_<number> folders expected" % adir)
    return sorted(f, key=int)

def load_frames(adir, states):
    adir = Path(adir)
    frames = frame_names(adir)
    G = {}
    for f in frames:
        for s in states:
            el, X = read_xyz(adir / ("frame_" + f) / STATES[s])
            if [e[0] for e in el] != [n[0] for n in NAMES]: die("frame %s %s: element order differs from the substrate's" % (f, s))
            G[(f, s)] = X
    return frames, G

def geometry(adir, states):
    """every geometry's 24 atoms stacked, in a fixed order (frame, then state as given)"""
    frames, G = load_frames(adir, states)
    keys = [(f, s) for f in frames for s in states]
    C = np.concatenate([G[k] for k in keys])
    return frames, G, keys, C, np.tile(np.arange(24), len(keys)), np.repeat(np.arange(len(keys)), 24)

def load_lj(prms, charges=False):
    L = Path(prms).read_text().splitlines()
    n = int(L[3].split()[0])
    rows = [l.split() for l in L[4:4 + n] if QM_FIRST <= int(l.split()[0]) <= QM_LAST]
    if len(rows) != 24 or [r[1][0] for r in rows] != [x[0] for x in NAMES]: die("the force field file does not hold the 24 substrate rows")
    out = (np.array([abs(float(r[3])) for r in rows]), np.array([float(r[4]) for r in rows]))
    return out + (np.array([float(r[2]) for r in rows]),) if charges else out

def geometry_sha(adir, keys):
    return ["%s:%s=%s" % (f, s, sha(Path(adir) / ("frame_" + f) / STATES[s])[:16]) for f, s in keys]

def same_geometries(com, adir, keys):
    """stop unless the sites file was built on exactly these aligned geometries (its header's sha256 line)"""
    line = [l for l in com if l.startswith("#   geometry sha256: ")]
    if not line: die("the sites file records no geometry sha256")
    want = dict(x.split("=") for x in line[0].split(": ", 1)[1].split())
    got = dict(x.split("=") for x in geometry_sha(adir, keys))
    if any(want.get(k) != v for k, v in got.items()): die("the aligned geometries in %s differ from those the grid was built on" % adir)

def read_floors(p):
    """floors by charge level (cs5_floors.py): {column: array over NAMES}"""
    F, hdr = {}, False
    for l in Path(p).read_text().splitlines():
        if not l.strip() or l.startswith("#"): continue
        t = l.split("\t")
        if t[0] == "atom":
            if t[1:7] != COLS: die("floors.tsv columns are not " + " ".join(COLS))
            hdr = True; continue
        F[t[0]] = [float(v) for v in t[1:7]]
    if not hdr or sorted(F) != sorted(NAMES): die("floors.tsv must give the floors of each of the 24 substrate atoms")
    out = {c: np.array([F[n][k] for n in NAMES]) for k, c in enumerate(COLS)}
    for a, b in zip(POS[:-1], POS[1:]):                    # a site cleared for one level must be cleared for every higher one
        if (out[b] > out[a] + 1e-12).any(): die("floors.tsv: the %s floors are not monotone in q (%s above %s)" % (b, b, a))
    return out

def fib_sphere(n):
    i = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * i / n); th = math.pi * (1 + 5 ** 0.5) * i
    return np.stack([np.cos(th) * np.sin(phi), np.sin(th) * np.sin(phi), np.cos(phi)], 1)

def pair_dist(A, B, chunk=20000):
    """distances between the rows of A and B, in chunks of A"""
    return np.concatenate([np.sqrt(((A[c:c + chunk, None, :] - B[None, :, :]) ** 2).sum(-1)) for c in range(0, len(A), chunk)]) \
        if len(A) else np.zeros((0, len(B)))

def admitted(dist, F, ngeo):
    """per site: clearance to each level's floors (min over atoms and geometries), q_min_neg, q_min_pos (None = none)"""
    clear = {c: (dist - np.tile(F[c], ngeo)[None, :]).min(1) for c in COLS}
    qneg = np.where(clear["neg"] >= -TOL, -1.0, 0.0)
    qpos = [next((LEVEL_LO[c] for c in POS if clear[c][i] >= -TOL), None) for i in range(dist.shape[0])]
    return clear, qneg, qpos

# ------------------------------------------------------------------ grid
def surface(C, R, atom_of, geo_of, ngeo, cfg):
    """sites on the union of spheres (centres C, radii R): outside every sphere, open to the outside, thinned to res.
    Returns the sites, the open candidates (for coverage) and the counts (lattice, outside, open)."""
    dens, res, vox, margin = float(cfg["density"]), float(cfg["res"]), float(cfg["vox"]), float(cfg["margin"])
    pts, src = [], []
    for a in range(len(C)):
        n = max(12, int(math.ceil(dens * 4 * math.pi * R[a] ** 2)))
        pts.append(C[a] + R[a] * fib_sphere(n)); src.append(np.full(n, a))
    P = np.concatenate(pts); S = np.concatenate(src)
    n_lattice = len(P)
    for g in range(ngeo):                                          # outside every sphere, geometry by geometry
        Cg, Rg = C[geo_of == g], R[geo_of == g]
        bad = np.zeros(len(P), bool)
        for c0 in range(0, len(P), 400000):                       # in chunks, to bound memory
            dist = np.sqrt(((P[c0:c0 + 400000, None, :] - Cg[None, :, :]) ** 2).sum(-1))
            bad[c0:c0 + 400000] = (dist < Rg[None, :] - 1e-9).any(1)
        P, S = P[~bad], S[~bad]
    n_cand = len(P)
    lo = C.min(0) - R.max() - margin; hi = C.max(0) + R.max() + margin   # enclosed pockets: flood fill from the box boundary
    shape = np.ceil((hi - lo) / vox).astype(int) + 1
    gx = [lo[k] + vox * np.arange(shape[k]) for k in range(3)]
    covered = np.zeros(shape, bool)
    for a in range(len(C)):
        i0 = np.maximum(((C[a] - R[a] - lo) / vox).astype(int) - 1, 0); i1 = np.minimum(((C[a] + R[a] - lo) / vox).astype(int) + 2, shape)
        sx, sy, sz = (gx[k][i0[k]:i1[k]] for k in range(3))
        dd = (sx[:, None, None] - C[a, 0]) ** 2 + (sy[None, :, None] - C[a, 1]) ** 2 + (sz[None, None, :] - C[a, 2]) ** 2
        covered[i0[0]:i1[0], i0[1]:i1[1], i0[2]:i1[2]] |= dd < (R[a] - 1e-6) ** 2
    out_r = np.zeros(shape, bool)
    out_r[0, :, :] = out_r[-1, :, :] = out_r[:, 0, :] = out_r[:, -1, :] = out_r[:, :, 0] = out_r[:, :, -1] = True
    out_r &= ~covered
    while True:
        g = out_r.copy()
        g[1:] |= out_r[:-1]; g[:-1] |= out_r[1:]; g[:, 1:] |= out_r[:, :-1]; g[:, :-1] |= out_r[:, 1:]
        g[:, :, 1:] |= out_r[:, :, :-1]; g[:, :, :-1] |= out_r[:, :, 1:]
        g &= ~covered
        if (g == out_r).all(): break
        out_r = g
    nrm = (P - C[S]) / R[S][:, None]
    idx = np.clip(np.round((P + 1.5 * vox * nrm - lo) / vox).astype(int), 0, shape - 1)
    outside = out_r[idx[:, 0], idx[:, 1], idx[:, 2]]
    P, S = P[outside], S[outside]
    n_open = len(P)
    order = np.lexsort((np.arange(len(P)), atom_of[S], geo_of[S]))   # deterministic greedy Poisson-disk thinning
    cell = {}; keep = []
    for j in order:
        p = P[j]; c = tuple(np.floor(p / res).astype(int)); bad = False
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    for k in cell.get((c[0] + dx, c[1] + dy, c[2] + dz), ()):
                        if ((P[k] - p) ** 2).sum() < res * res: bad = True; break
                    if bad: break
                if bad: break
            if bad: break
        if not bad: keep.append(j); cell.setdefault(c, []).append(j)
    return P[np.array(keep)], P, (n_lattice, n_cand, n_open)

def held(p, qpos, anion, res):
    """True if an anion/neutral-layer site within res of p admits every positive charge that qpos admits (every
    anion/neutral site admits q <= 0)"""
    Qa, qa = anion
    near = np.where(np.sqrt(((Qa - p) ** 2).sum(1)) < res)[0]
    return any(qa[j] is not None and (qpos is None or qa[j] <= qpos + 1e-9) for j in near)

def build_grid(adir, prms, floors_p, out, args):
    cfg = cfg_of(args)
    if Path(out).exists(): die("%s exists; give a new file name" % out)
    states = cfg["states"].split(",")
    frames, G, keys, C, atom_of, geo_of = geometry(adir, states)
    ngeo = len(keys)
    F = read_floors(floors_p)
    load_lj(prms)                                                 # the force field must hold the substrate rows
    rows, stats, anion = [], {}, None
    for lay in ("anion", "cation"):
        Q, cand, st = surface(C, np.tile(F[LAYER_COL[lay]], ngeo), atom_of, geo_of, ngeo, cfg)
        cov = pair_dist(cand, Q).min(1).max() if len(Q) else float("nan")     # every open candidate lies within res of a site
        dist = pair_dist(Q, C)
        clear, qneg, qpos = admitted(dist, F, ngeo)
        own = clear[LAYER_COL[lay]]
        if (own < -TOL).any(): die("internal: a %s site lies inside its own layer's floors" % lay)
        on = (dist - np.tile(F[LAYER_COL[lay]], ngeo)[None, :]).argmin(1)
        dmin = np.stack([dist[:, atom_of == a].min(1) for a in range(24)], 1)
        kept = dropped = 0
        for i in range(len(Q)):
            if lay == "cation" and abs(clear["neg"][i]) <= TOL and held(Q[i], qpos[i], anion, float(cfg["res"])):
                dropped += 1; continue            # on the anion/neutral surface, and that layer holds this place
            rows.append((lay, Q[i], NAMES[atom_of[on[i]]], "%s:%s" % keys[geo_of[on[i]]], own[i], qneg[i], qpos[i], dmin[i]))
            kept += 1
        if lay == "anion": anion = (Q, qpos)
        stats[lay] = st + (len(Q), dropped, kept, cov)
    allQ = np.array([r[1] for r in rows])
    nn = pair_dist(allQ, allQ); np.fill_diagonal(nn, np.inf); nn = nn.min(1)
    na = sum(1 for r in rows if r[0] == "anion"); nc = len(rows) - na
    adm = {}
    for r in rows:
        k = ("q<=0 " if r[5] < 0 else "") + ("q>0" if r[6] == 0.0 else "q>=%.2f" % r[6] if r[6] is not None else "")
        adm[(r[0], k.strip() or "nothing")] = adm.get((r[0], k.strip() or "nothing"), 0) + 1
    if any(k[1] == "nothing" for k in adm): die("internal: a site admits no charge")
    hdr = ["# candidate-site grid, cs5_grid.py grid (Stage 5, STAGE5_CRITERIA.txt); cs5_grid.py sha256 %s" % sha(__file__),
           "# two layers, each built as Dittner's thesis Eq. 2.16 (union over every atom of every geometry of atom-centred spheres),",
           "#   radii = the floors of floors.tsv: ANION/NEUTRAL layer on the q <= 0 floors, CATION layer on the q = +1 floors;",
           "#   open to the outside (flood fill); deterministic greedy Poisson-disk thinning at res within each layer; cation-layer",
           "#   sites ON the anion/neutral surface (clearance to the q <= 0 floors zero within %g A) dropped where an anion/neutral" % TOL,
           "#   site within res admits every charge they admit: that layer holds the place",
           "# each site admits every charge whose floors it clears at every geometry: q in [q_min_neg, 0] if q_min_neg = -1;",
           "#   q in [q_min_pos, +1] (q_min_pos 0.00: any q > 0; none: no positive charge); an empty site is always allowed",
           "# frames (%d): %s" % (len(frames), " ".join(frames)),
           "# states: %s; aligned geometries %s" % (",".join(states), adir),
           "#   geometry sha256: " + " ".join(geometry_sha(adir, keys)),
           "# force field %s sha256 %s" % (prms, sha(prms)), "# floors %s sha256 %s" % (floors_p, sha(floors_p)),
           "# parameters: " + " ".join("%s=%s" % kv for kv in sorted(cfg.items())),
           "# anion/neutral layer: lattice %d, outside %d, open %d, thinned %d, dropped %d, sites %d; coverage max %.4f A" % stats["anion"],
           "# cation layer:        lattice %d, outside %d, open %d, thinned %d, dropped %d, sites %d; coverage max %.4f A" % stats["cation"],
           "#   (coverage: largest distance from an open candidate point to its layer's nearest thinned site, below res by construction;",
           "#   a dropped cation site has an anion/neutral site within res admitting its charges, so every point is within 2 res)",
           "# sites by admitted charges: " + "; ".join("%s %s %d" % (k[0], k[1], v) for k, v in sorted(adm.items())),
           "# nearest neighbour over both layers: min %.3f median %.3f max %.3f A" % (nn.min(), float(np.median(nn)), nn.max()),
           "\t".join(HEAD)]
    with open(out, "x") as o:
        o.write("\n".join(hdr) + "\n")
        for i, r in enumerate(rows):
            o.write("\t".join([str(i), r[0]] + [fmt(v, 6) for v in r[1]] + [r[2], r[3], fmt(r[4], 4), fmt(r[5], 2),
                              "none" if r[6] is None else fmt(r[6], 2)] + [fmt(v, 4) for v in r[7]]) + "\n")
    print("grid: %d frames x %d states; anion/neutral layer %d sites, cation layer %d sites (%d dropped on the anion/neutral surface);"
          " resolution %s A; nearest neighbour %.3f-%.3f A" % (len(frames), len(states), na, nc, stats["cation"][4], cfg["res"], nn.min(), nn.max()))

# ------------------------------------------------------------------ reading a sites file
def read_sites(p):
    """(rows as dicts, header comment lines); accepts this script's format and grid_v2's (idx x y z ...)"""
    rows, h, com = [], None, []
    for l in Path(p).read_text().splitlines():
        if not l.strip(): continue
        if l.startswith("#"): com.append(l); continue
        t = l.split("\t")
        if t[0] in ("site", "idx"): h = t; continue
        if h is None: die("%s: no column header" % p)
        rows.append(dict(zip(h, t)))
    return rows, com

def sites_xyz(rows): return np.array([[float(r["x"]), float(r["y"]), float(r["z"])] for r in rows])
def frames_of(com):
    m = [l for l in com if l.startswith("# frames (")]
    return m[0].split(":", 1)[1].split() if m else None

# ------------------------------------------------------------------ check
def check(adir, prms, floors_p, sites_p, args):
    rows, com = read_sites(sites_p)
    if not rows or list(rows[0]) != HEAD: die("%s is not a cs5_grid.py sites file" % sites_p)
    states = [l for l in com if l.startswith("# states: ")][0].split()[2].rstrip(";").split(",")
    frames, G, keys, C, atom_of, geo_of = geometry(adir, states)
    ngeo = len(keys)
    if frames_of(com) != frames: die("the sites file's frames differ from %s's" % adir)
    same_geometries(com, adir, keys)
    want = [l for l in com if l.startswith("#   geometry sha256: ")][0].split(": ", 1)[1].split()
    if want != geometry_sha(adir, keys): die("the sites file's geometries are not exactly these, in this order")
    F = read_floors(floors_p)
    if not any(l.startswith("# floors ") and l.endswith(sha(floors_p)) for l in com): die("floors.tsv differs from the one the grid was built on")
    if not any(l.startswith("# force field ") and l.endswith(sha(prms)) for l in com): die("the force field differs from the one the grid was built on")
    Q = sites_xyz(rows)
    dist = pair_dist(Q, C)
    clear, qneg, qpos = admitted(dist, F, ngeo)
    res = float([l for l in com if l.startswith("# parameters:")][0].split("res=")[1].split()[0])
    ma = np.array([r["layer"] == "anion" for r in rows])
    anion = (Q[ma], [qpos[i] for i in np.where(ma)[0]])
    bad = []
    for i, r in enumerate(rows):
        if int(r["site"]) != i: bad.append("row %d: site number %s" % (i, r["site"]))
        lay = r["layer"]
        if lay not in LAYER_COL: bad.append("site %d: layer %s" % (i, lay)); continue
        own = clear[LAYER_COL[lay]][i]
        if own < -TOL: bad.append("site %d: inside its own layer's floors by %.6f A" % (i, -own))
        if abs(own - float(r["clearance_A"])) > 6e-5: bad.append("site %d: clearance %.6f, file %s" % (i, own, r["clearance_A"]))
        if lay == "cation" and abs(clear["neg"][i]) <= TOL and held(Q[i], qpos[i], anion, res):
            bad.append("site %d: a cation site on the anion/neutral surface that an anion/neutral site holds" % i)
        if abs(qneg[i] - float(r["q_min_neg"])) > 1e-9: bad.append("site %d: q_min_neg %s, should be %.2f" % (i, r["q_min_neg"], qneg[i]))
        want = "none" if qpos[i] is None else "%.2f" % qpos[i]
        if r["q_min_pos"] != want: bad.append("site %d: q_min_pos %s, should be %s" % (i, r["q_min_pos"], want))
        if qneg[i] == 0.0 and qpos[i] is None: bad.append("site %d admits no charge" % i)
        Rl = np.tile(F[LAYER_COL[lay]], ngeo)
        j = int((dist[i] - Rl).argmin())
        if r["on_atom"] != NAMES[atom_of[j]] or r["on_state"] != "%s:%s" % keys[geo_of[j]]:
            bad.append("site %d: on %s %s, file %s %s" % (i, NAMES[atom_of[j]], "%s:%s" % keys[geo_of[j]], r["on_atom"], r["on_state"]))
        for a, n in enumerate(NAMES):
            if abs(dist[i, atom_of == a].min() - float(r["d_" + n])) > 6e-5: bad.append("site %d: d_%s" % (i, n))
        # every admitted charge clears its own floors at every atom of every geometry
        if qneg[i] < 0 and (dist[i] < np.tile(F["neg"], ngeo) - TOL).any(): bad.append("site %d: q <= 0 inside a floor" % i)
        if qpos[i] is not None:
            for c in POS:
                if LEVEL_LO[c] >= qpos[i] and (dist[i] < np.tile(F[c], ngeo) - TOL).any(): bad.append("site %d: q >= %.2f inside a floor" % (i, LEVEL_LO[c]))
    nn = pair_dist(Q, Q); np.fill_diagonal(nn, np.inf)
    if nn.min() < 1e-6: bad.append("two sites coincide")
    lay = [r["layer"] for r in rows]
    for L_ in LAYER_COL:                                      # within a layer, the thinning keeps sites at least res apart
        m = np.array([x == L_ for x in lay])
        if m.sum() > 1 and nn[np.ix_(m, m)].min() < res - 1e-9: bad.append("%s layer: two sites closer than res" % L_)
    if bad: die("check FAILED (%d):\n  " % len(bad) + "\n  ".join(bad[:40]))
    for line in radial(rows, Q, G, keys, prms): print(line)
    n_c = sum(1 for x in lay if x == "cation")
    print("check PASS: %d sites (anion/neutral %d, cation %d) over %d frames x %s; every admitted charge clears its floors at"
          " every geometry; columns reproduced" % (len(rows), len(rows) - n_c, n_c, len(frames), ",".join(states)))

def radial(rows, Q, G, keys, prms):
    """Stage 2's rigid model at each site that admits a positive charge: the outward force at every geometry for the
    lowest admitted positive charge (0.05 e for any q > 0) and for +1; informative lines"""
    eps, rmin, qs = load_lj(prms, charges=True)
    K = 332.0637                                                  # kcal A / (mol e^2), as wp5_rigid_lj_scan.py
    Rij = 0.5 * (2 * SPHERE_RSTAR + rmin); eij = np.sqrt(SPHERE_EPS * eps)
    idx = [i for i, r in enumerate(rows) if r["q_min_pos"] != "none"]
    nrm = np.zeros((len(idx), 3))
    for k, i in enumerate(idx):
        f, s = rows[i]["on_state"].split(":"); v = Q[i] - G[(f, s)][NAMES.index(rows[i]["on_atom"])]; nrm[k] = v / np.linalg.norm(v)
    P = Q[idx]
    def force_out(qq, P):                                         # outward force on a sphere of charge qq, per geometry
        F = []
        for key in keys:
            d = P[:, None, :] - G[key][None, :, :]; r = np.sqrt((d * d).sum(-1)); s6 = (Rij[None, :] / r) ** 6
            dEdr = -K * qq[:, None] * qs[None, :] / r ** 2 + eij[None, :] * (-12.0 * s6 * s6 + 12.0 * s6) / r
            g = (dEdr[:, :, None] * d / r[:, :, None]).sum(1)
            F.append(-(g * nrm[:len(P)]).sum(1))
        return np.array(F)                                        # geometries x sites
    out = []
    for label, qq in (("lowest admitted q", np.array([max(float(rows[i]["q_min_pos"]), 0.05) for i in idx])), ("q = +1", np.ones(len(idx)))):
        F = force_out(qq, P); pushed = (F > 0).any(0)
        far = 0.0
        for k in np.where(pushed)[0]:                             # radial minimum along the normal, worst geometry
            t = np.arange(0.0, 1.5001, 0.01); X = P[k] + t[:, None] * nrm[k]
            for key in keys:
                d = X[:, None, :] - G[key][None, :, :]; r = np.sqrt((d * d).sum(-1)); s6 = (Rij[None, :] / r) ** 6
                E = (K * qq[k] * qs[None, :] / r + eij[None, :] * (s6 * s6 - 2 * s6)).sum(1)
                far = max(far, t[E.argmin()])
        out.append("radial (rigid model, informative), %s: pushed outward at some geometry at %d of %d sites admitting q > 0;"
                   " largest outward force %.2f kcal/mol/A; radial minimum at most %.2f A farther out" % (
                       label, int(pushed.sum()), len(idx), F.max(), far))
    return out

# ------------------------------------------------------------------ LJ term
def lj_energy(Q, X, eps, rmin):
    """E_LJ (kcal/mol) of Amber N spheres at the rows of Q with the 24 substrate atoms X; one value per sphere"""
    Rij = 0.5 * (2 * SPHERE_RSTAR + rmin); eij = np.sqrt(SPHERE_EPS * eps)
    r = pair_dist(Q, X)
    s6 = (Rij[None, :] / r) ** 6
    return (eij[None, :] * (s6 * s6 - 2 * s6)).sum(1)

def lj_term(adir, prms, sites_p, out, args):
    if args: die("lj takes no parameters")
    if Path(out).exists(): die("%s exists; give a new file name" % out)
    frames, G = load_frames(adir, ["R", "TS"])
    rows, com = read_sites(sites_p)
    if frames_of(com) != frames: die("the sites file's frames differ from %s's" % adir)
    same_geometries(com, adir, [(f, s) for f in frames for s in ("R", "TS")])
    eps, rmin = load_lj(prms)
    Q = sites_xyz(rows)
    L = np.stack([lj_energy(Q, G[(f, "TS")], eps, rmin) - lj_energy(Q, G[(f, "R")], eps, rmin) for f in frames], 1)
    Rij = 0.5 * (2 * SPHERE_RSTAR + rmin); eij = np.sqrt(SPHERE_EPS * eps); nfar, far_site = 0, np.zeros(len(Q))
    for f in frames:                                   # pairs at 10 A or more, where ORCA's force switch acts
        for s in ("R", "TS"):
            r = pair_dist(Q, G[(f, s)]); s6 = (Rij[None, :] / r) ** 6; m = r >= 10.0
            nfar += int(m.sum()); far_site = np.maximum(far_site, (np.abs(eij[None, :] * (s6 * s6 - 2 * s6)) * m).sum(1))
    hdr = ["# L matrix, cs5_grid.py lj (Stage 5, STAGE5_CRITERIA.txt); cs5_grid.py sha256 %s" % sha(__file__),
           "# L[i,f] = E_LJ(Amber N sphere at site i, substrate at TS_f) - E_LJ(same, at R_f), kcal/mol; independent of q;",
           "#   applies when a sphere occupies site i",
           "# sphere R* %.3f A eps %.3f kcal/mol; substrate LJ from %s (sha256 %s); combining R_min,ij = (R_min,s + R_min,J)/2,"
           % (SPHERE_RSTAR, SPHERE_EPS, prms, sha(prms)),
           "#   eps_ij = sqrt(eps_s eps_J), as ORCA (Stage 1 S1.4)",
           "# ORCA force-switches LJ from 10 A; below 10 A that is one constant per pair (WP1), which cancels in TS - R. Site-atom",
           "#   pairs at 10 A or more: %d of %d; their summed |E_LJ| is at most %.4f kcal/mol at any site and geometry, the order of" % (
               nfar, len(Q) * 24 * 2 * len(frames), far_site.max()),
           "#   what the switch changes in L (ljcheck against ORCA: within 0.0005 kcal/mol on Stage 4c's arrangements)",
           "# sites: %d from %s (sha256 %s); row i = site i" % (len(Q), sites_p, sha(sites_p)),
           "# frames (%d): %s" % (len(frames), " ".join(frames)),
           "# aligned geometries %s" % adir,
           "# range %.3f to %.3f kcal/mol; median %.3f" % (L.min(), L.max(), float(np.median(L)))]
    with open(out, "x") as o:
        o.write("\n".join(hdr) + "\n")
        for row in L: o.write("\t".join(fmt(v, 6) for v in row) + "\n")
    print("L: %d sites x %d frames; range %.3f to %.3f kcal/mol" % (L.shape[0], L.shape[1], L.min(), L.max()))

def ljcheck(res_repo):
    """the LJ formula against ORCA: Stage 4c's L1 runs (bare Amber N spheres, QM/MM, frame 41786's fixed R and TS)"""
    D = Path(res_repo) / "phase2.2/cs_stage4c"
    eps, rmin = load_lj(Path(res_repo) / "05_qmmm/13_bridge/complex_solvated.ORCAFF.prms")
    ref = {}
    for l in (D / "STAGE4C_REPORT.txt").read_text().splitlines():
        m = re.match(r"^\s+(E\d)\s+\d+\s+(\S+\s+){5}([+-]\d+\.\d{3})\s*$", l)
        if m: ref[m.group(1)] = float(m.group(3))
    if sorted(ref) != ["E%d" % k for k in range(1, 9)]: die("STAGE4C_REPORT.txt: the dE_LJ(L1) column was not found for E1-E8")
    worst = 0.0
    for e in sorted(ref):
        E = {}
        for g in ("R", "TS"):
            d = D / ("L1_%s_%s" % (e, g))
            L = (d / "combined.xyz").read_text().splitlines()
            n = int(L[0]); body = [l.split() for l in L[2:2 + n]]
            X = np.array([[float(v) for v in b[1:4]] for b in body[:24]]); Q = np.array([[float(v) for v in b[1:4]] for b in body[24:]])
            P = [l.split() for l in (d / "site.ORCAFF.prms").read_text().splitlines()[4:4 + n]]
            if [b[0][0] for b in body[:24]] != [x[0] for x in NAMES]: die("%s: substrate order" % d)
            if not np.allclose([abs(float(p[3])) for p in P[:24]], eps, atol=1e-9) or not np.allclose([float(p[4]) for p in P[:24]], rmin, atol=1e-9):
                die("%s: substrate LJ differs from the committed force field" % d)
            if any(p[1] != "N" or abs(abs(float(p[3])) - SPHERE_EPS) > 1e-9 or abs(float(p[4]) - 2 * SPHERE_RSTAR) > 1e-9 for p in P[24:]):
                die("%s: sphere LJ is not Amber N" % d)
            E[g] = (lj_energy(Q, X, eps, rmin).sum(), Q)
        if not np.array_equal(E["R"][1], E["TS"][1]): die("%s: spheres differ between R and TS" % e)
        dE = E["TS"][0] - E["R"][0]; worst = max(worst, abs(dE - ref[e]))
        print("  %s  %d sphere(s)  cs5 %+.4f  ORCA %+.3f  diff %+.4f" % (e, len(E["R"][1]), dE, ref[e], dE - ref[e]))
    ok = worst <= 0.0015
    print("ljcheck %s: largest |difference| %.4f kcal/mol (tolerance 0.0015: the report's 3-decimal rounding plus margin)"
          % ("PASS" if ok else "FAIL", worst))
    sys.exit(0 if ok else 1)

if __name__ == "__main__":
    a = sys.argv[1:]
    if a[:1] == ["ljcheck"] and len(a) == 2: ljcheck(a[1])
    elif len(a) >= 5 and a[0] in ("grid", "lj", "check"):
        {"grid": build_grid, "lj": lj_term, "check": check}[a[0]](*a[1:5], a[5:])
    else: die(__doc__.split("\n\n")[1])
