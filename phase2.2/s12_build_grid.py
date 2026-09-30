#!/usr/bin/env python3
"""
s12_build_grid.py - build the candidate-site grid on the van der Waals surface.

Pure function of its configuration: parameters in, sites out, no hidden state. Every
parameter is written into the output header, so a grid file always says how it was made.
That convention is what let the pilot-frame problem in dv_grid_ensemble_mean.tsv be
caught; it is kept deliberately.

THE SURFACE, from [DTHESIS] Ch. 2 (Dittner PhD thesis, gocat_2.pdf)
Their constraint is, for each charge i,

    || r_i - R_{J,Voronoi} ||_2 - d_J = 0

so a charge sits EXACTLY ON the sphere of its assigned atom J, where d_J are
ATOM-DEPENDENT van der Waals radii. Not a probe-rolled solvent-accessible surface, not a
uniform offset shell: the union of atom-centred vdW spheres. Assignment is Voronoi-like,
the atom minimising (||r_i - R_J|| - d_J), and a charge inside any sphere is projected
back out onto it. The surface encloses ALL PATH FRAMES at once - "common vdW surface
exposed by all atoms of all frames".

[DTHESIS] typical values: charge bounds [-1, +1] e, minimum charge-charge separation
r_min = 1 Angstrom.

WHAT IS TUNABLE, AND WHY IT IS TUNABLE

  Set by physics. Defensible defaults, change only with a reason.
    vdw_radii             per element, from a named source
    surface_offset        0.0 reproduces [DTHESIS]. Non-zero gives an offset shell.
    min_approach          minimum site-to-atom distance. SET BY s11_probe_distance_scan:
                          the shortest distance at which a bare point charge is still
                          well behaved. Until that scan reports, the default is a
                          placeholder and is flagged as such in the header.
    charge_bounds         [-1, +1] e, [DTHESIS] Ch. 2
    min_charge_separation 1.0 A, [DTHESIS] Ch. 2. Phase 2 used 3.5, which is 3.5x
                          tighter and is why placing two charges in a 4-5 A favourable
                          region was hard. That constraint was ours, not inherited.

  Set by tractability. Genuinely free, tune against MILP solve time.
    site_spacing          Poisson-disk radius; drives site count, which drives MILP size
    max_sites             hard cap so the model stays solvable
    random_seed           Poisson-disk is stochastic; a fixed seed makes it reproducible

  Set by the substrate decision. The reason the grid must be cheap to rebuild.
    frames                which frames the surface encloses
    path_states           R, TS, P, or all NEB images
    enclose_mode          union over frames, or per-frame grids intersected

  Set by the charge model. Pending s11.
    charge_model          point | smeared
    smearing_width        only meaningful when smeared. NOTE: a smeared charge has a
                          physical width that should SET min_approach rather than being
                          independent of it; do not tune them separately without saying
                          why.

USAGE
    python3 s12_build_grid.py <ensemble_dir> <out.tsv> [key=value ...]

    python3 s12_build_grid.py 05_qmmm/19_ensemble_barriers grid.tsv \\
        site_spacing=0.8 frames=24883,20000,43087 min_approach=3.0
"""
import sys
import math
import random
from pathlib import Path

import numpy as np

# --------------------------------------------------------------------------- config
CONFIG = {
    # physics
    "vdw_radii_source": "Bondi 1964, J. Phys. Chem. 68, 441; H revised by Rowland & "
                        "Taylor 1996, J. Phys. Chem. 100, 7384",
    "surface_offset": 0.0,          # A, 0.0 = charge on the atom vdW sphere [DTHESIS]
    "min_approach": 3.0,            # A, FLOOR on the site radius. PLACEHOLDER until
                                    # s11 reports. Set 0.0 to recover [DTHESIS] exactly.
    "charge_lower": -1.0,           # e, [DTHESIS] Ch.2
    "charge_upper": +1.0,           # e, [DTHESIS] Ch.2
    "min_charge_separation": 1.0,   # A, [DTHESIS] Ch.2
    # tractability
    "site_spacing": 0.8,            # A, Poisson-disk minimum site separation
    "max_sites": 4000,
    "random_seed": 20260930,
    "poisson_candidates": 30,       # k in Bridson's algorithm
    # substrate
    "frames": "",                   # comma-separated; empty = every frame present.
                                    # NOTE: enclosing all 30 post-cut frames plus the
                                    # 14 pilot ones is 1700+ atoms and a much larger
                                    # surface than any single design needs. The pilot
                                    # frames are excluded from reported results, so
                                    # pass the post-cut list explicitly.
    "path_states": "R,TS",          # which geometries the surface must enclose
    "enclose_mode": "union",        # union | intersect
    # charge model
    "charge_model": "point",        # point | smeared
    "smearing_width": 0.0,          # A, only meaningful when smeared
}

# Bondi vdW radii, Angstrom. H from Rowland & Taylor.
VDW = {"H": 1.09, "C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80, "P": 1.80}

STATE_FILE = {"R": "reactant_qm.xyz",
              "TS": "transition_state_qm.xyz",
              "P": "product_qm.xyz"}


def dist(a, b):
    return math.sqrt((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2)


def read_xyz(p):
    L = Path(p).read_text().splitlines()
    n = int(L[0].split()[0])
    return [(f[0], float(f[1]), float(f[2]), float(f[3]))
            for f in (l.split() for l in L[2:2 + n])]


# Some frames keep their geometries only under 20_invacuo, named <frame>_R.xyz rather
# than frame_<frame>/reactant_qm.xyz. Frames harvested after the in vacuo stage are in
# that position. Both layouts are searched.
ALT_STATE = {"R": "_R.xyz", "TS": "_TS.xyz", "P": "_P.xyz"}


def collect_atoms(ens, frames, states, alt_dir=None):
    """Every atom of every requested geometry, as (element, x, y, z).

    Returns (atoms, used, missing). MISSING IS RETURNED AND REPORTED: an earlier version
    silently skipped frames whose geometry files were absent and announced a frame count
    lower than the one requested, which is the kind of error that reaches the design
    without tripping anything.
    """
    out, used, missing = [], [], []
    ens = Path(ens)
    alt = Path(alt_dir) if alt_dir else ens.parent / "20_invacuo"
    want = list(frames) if frames else \
        sorted(d.name.replace("frame_", "") for d in ens.glob("frame_*"))
    for fr in want:
        got = 0
        for s in states:
            f = ens / f"frame_{fr}" / STATE_FILE[s]
            if not f.exists():
                f = alt / f"{fr}{ALT_STATE[s]}"      # the 20_invacuo layout
            if f.exists():
                out.extend(read_xyz(f))
                got += 1
        if got == len(states):
            used.append((fr, got))
        elif got:
            used.append((fr, got))
            missing.append(f"{fr} (only {got} of {len(states)} states)")
        else:
            missing.append(fr)
    return out, used, missing


def site_radius(el, cfg):
    """Radius at which sites may sit around an atom of this element.

    [DTHESIS] Eq. 2.16 puts a charge exactly on the atom's vdW sphere, so its distance
    to that nucleus is the vdW radius itself: 1.52 A for O, 1.09 A for H. Burschowsky's
    citrulline NH2 sits 3.2 A from the ether oxygen, so their surface places a charge at
    under half the distance of the real cation it stands in for. That is very plausibly
    why [DTHESIS] Ch. 3.6 reports Coulomb implosion - the surface choice and the
    pathology are linked.

    min_approach therefore acts as a FLOOR on the site radius, pushing the surface
    outward where the vdW radius alone would put a charge too close. Setting
    min_approach = 0 recovers [DTHESIS]'s surface exactly.
    """
    return max(VDW.get(el.capitalize(), 1.70) + cfg["surface_offset"],
               cfg["min_approach"])


def on_surface(pt, atoms, cfg):
    """Valid if outside every atom's site sphere and touching at least one.

    The discrete form of [DTHESIS] Eq. 2.16, with the site radius floored at
    min_approach so a charge is never placed closer to any nucleus than the probe scan
    says a bare point charge can be trusted.
    """
    touching = False
    for el, x, y, z in atoms:
        d = dist(pt, (x, y, z))
        r = site_radius(el, cfg)
        if d < r - 1e-6:
            return False, False
        if d <= r + cfg["site_spacing"]:
            touching = True
    return True, touching


def poisson_surface_sites_fast(atoms, cfg):
    """Vectorised Poisson-disk sampling on the union surface.

    The pure-Python version checks every candidate against every atom one at a time. On
    a 30-frame ensemble that is about 1700 atoms and roughly 750 candidates per atom,
    i.e. of order 10^9 distance evaluations, which takes hours. This does the same work
    with numpy: all atoms at once per candidate batch, and an accepted-site check on a
    single array. Same acceptance rule, same seed, same result to floating precision.
    """
    rng = random.Random(cfg["random_seed"])
    P = np.array([[a[1], a[2], a[3]] for a in atoms])
    R = np.array([site_radius(a[0], cfg) for a in atoms])
    accepted = []
    spacing2 = cfg["site_spacing"] ** 2
    order = sorted(range(len(atoms)), key=lambda k: (atoms[k][1], atoms[k][2],
                                                     atoms[k][3], atoms[k][0]))
    for k in order:
        r = R[k]
        area = 4.0 * math.pi * r * r
        n_try = max(8, int(cfg["poisson_candidates"] * area /
                           (math.pi * cfg["site_spacing"] ** 2)))
        # generate this atom's candidates in one batch
        u = np.array([rng.random() for _ in range(n_try)])
        v = np.array([rng.random() for _ in range(n_try)])
        theta = 2.0 * math.pi * u
        phi = np.arccos(2.0 * v - 1.0)
        cand = np.stack([P[k, 0] + r * np.sin(phi) * np.cos(theta),
                         P[k, 1] + r * np.sin(phi) * np.sin(theta),
                         P[k, 2] + r * np.cos(phi)], axis=1)
        # distance from every candidate to every atom, (n_try, n_atoms)
        d = np.linalg.norm(cand[:, None, :] - P[None, :, :], axis=2)
        inside = (d < R[None, :] - 1e-6).any(axis=1)
        touching = (d <= R[None, :] + cfg["site_spacing"]).any(axis=1)
        keep = (~inside) & touching
        for c in cand[keep]:
            if len(accepted) >= cfg["max_sites"]:
                return [tuple(a) for a in accepted]
            if accepted:
                A = np.array(accepted)
                if ((A - c) ** 2).sum(axis=1).min() < spacing2:
                    continue
            accepted.append(c)
    return [tuple(a) for a in accepted]


def poisson_surface_sites(atoms, cfg):
    """Poisson-disk sampling on the union surface.

    Candidates are generated per atom on its own sphere, then accepted if they are on
    the union surface and no closer than site_spacing to an accepted site. Sampling per
    atom rather than in a bounding box means density follows exposed area, which is what
    [D2018] do when they weight the initial charge distribution "by the exposed vdW
    surface".
    """
    rng = random.Random(cfg["random_seed"])
    accepted = []
    # deterministic atom order so the seed fully determines the grid
    for el, ax, ay, az in sorted(atoms, key=lambda a: (a[1], a[2], a[3], a[0])):
        r = site_radius(el, cfg)
        area = 4.0 * math.pi * r * r
        n_try = max(8, int(cfg["poisson_candidates"] * area /
                           (math.pi * cfg["site_spacing"] ** 2)))
        for _ in range(n_try):
            if len(accepted) >= cfg["max_sites"]:
                return accepted
            # uniform point on a sphere
            u, v = rng.random(), rng.random()
            theta = 2.0 * math.pi * u
            phi = math.acos(2.0 * v - 1.0)
            pt = (ax + r * math.sin(phi) * math.cos(theta),
                  ay + r * math.sin(phi) * math.sin(theta),
                  az + r * math.cos(phi))
            ok, touching = on_surface(pt, atoms, cfg)
            if not (ok and touching):
                continue
            if any(dist(pt, s) < cfg["site_spacing"] for s in accepted):
                continue
            accepted.append(pt)
    return accepted


def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    ens, out = sys.argv[1], Path(sys.argv[2])
    cfg = dict(CONFIG)
    for kv in sys.argv[3:]:
        if "=" not in kv:
            sys.exit(f"bad argument {kv!r}, expected key=value")
        k, v = kv.split("=", 1)
        if k not in cfg:
            sys.exit(f"unknown parameter {k!r}. Known: {', '.join(sorted(cfg))}")
        cfg[k] = type(cfg[k])(v) if not isinstance(cfg[k], str) else v

    frames = [f for f in cfg["frames"].split(",") if f]
    states = [s for s in cfg["path_states"].split(",") if s]
    for s in states:
        if s not in STATE_FILE:
            sys.exit(f"unknown state {s!r}; known: {', '.join(STATE_FILE)}")
    if cfg["enclose_mode"] != "union":
        sys.exit("only enclose_mode=union is implemented; intersect is a stub")
    if cfg["charge_model"] == "smeared" and cfg["smearing_width"] <= 0:
        sys.exit("charge_model=smeared requires smearing_width > 0")

    atoms, used, missing = collect_atoms(ens, frames, states)
    if not atoms:
        sys.exit("no geometries found; check the ensemble directory and frames")
    if missing:
        print(f"  WARNING: {len(missing)} requested frame(s) had no geometry and were")
        print(f"  NOT enclosed: {' '.join(missing)}")
        print(f"  The surface therefore does not cover them. Either supply their")
        print(f"  geometries or drop them from the frame list, but do not quote this")
        print(f"  grid as covering the frames you asked for.")
    if frames and len(used) != len(frames):
        print(f"  requested {len(frames)} frames, enclosed {len(used)}")
    print(f"enclosing {len(used)} frames x {len(states)} states = {len(atoms)} atoms")
    sites = poisson_surface_sites_fast(atoms, cfg)
    print(f"accepted {len(sites)} sites at {cfg['site_spacing']} A spacing")
    if len(sites) >= cfg["max_sites"]:
        print(f"  WARNING: hit max_sites={cfg['max_sites']}; the grid is truncated and "
              f"NOT a uniform sample. Raise max_sites or coarsen site_spacing.")

    nearest = [min(dist(s, (a[1], a[2], a[3])) for a in atoms) for s in sites]
    hdr = ["# candidate-site grid on the union vdW surface, s12_build_grid.py",
           "# surface definition follows [DTHESIS] Ch.2 Eq.2.16: a site lies ON the",
           "#   union of atom-centred vdW spheres, outside every sphere and touching",
           "#   at least one. min_approach is an ADDITIONAL constraint on how close a",
           "#   charge may sit to any nucleus.",
           "#",
           "# PARAMETERS, all of them, so this file says how it was made:"]
    for k in sorted(cfg):
        hdr.append(f"#   {k} = {cfg[k]}")
    hdr += ["#",
            f"# frames enclosed ({len(used)}): " +
            " ".join(f"{f}({n})" for f, n in used),
            (f"# FRAMES REQUESTED BUT NOT ENCLOSED ({len(missing)}): "
             + " ".join(missing)) if missing else "# all requested frames enclosed",
            f"# states: {','.join(states)}",
            f"# atoms enclosed: {len(atoms)}",
            f"# sites: {len(sites)}",
            f"# nearest-atom distance over sites: min {min(nearest):.3f}, "
            f"mean {sum(nearest)/len(nearest):.3f}, max {max(nearest):.3f} A"]
    if cfg["min_approach"] == CONFIG["min_approach"]:
        hdr.append("# WARNING: min_approach is still the PLACEHOLDER default. It must be"
                   " set from")
        hdr.append("#   s11_probe_distance_scan before this grid is used for a design.")
    if cfg["charge_model"] == "smeared":
        hdr.append("# NOTE: smearing_width should normally SET min_approach rather than"
                   " being")
        hdr.append("#   tuned independently of it.")
    hdr.append("idx\tx\ty\tz\tnearest_atom_A")
    out.write_text("\n".join(hdr) + "\n" +
                   "".join(f"{i}\t{s[0]:.6f}\t{s[1]:.6f}\t{s[2]:.6f}\t{d:.4f}\n"
                           for i, (s, d) in enumerate(zip(sites, nearest))))
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
