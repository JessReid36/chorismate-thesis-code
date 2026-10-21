"""prep_charges.py - the optimal-charges map's data.
Usage: python3 prep_charges.py <results_repo> <out_dir>   (writes <out_dir>/data.json)

Reads only committed Stage 6 files: the frame 42436 substrate geometries cs6_prepare.py extracted, the one-frame
grid and its L, the A column cs5_amatrix.py assembled, and relax/ARRANGEMENTS.txt for the seven placed designs.
Nothing is recomputed that those files already fix; the per-sphere contribution a*q + L is recomputed here and
checked against each design's reported objective."""
import hashlib, json, os, re, sys
import numpy as np

ROOT, OUT = sys.argv[1], sys.argv[2]
T6 = os.path.join(ROOT, "phase2.2", "cs_stage6")
NAMES = "C1 H1 H2 C2 C3 O1 O2 O3 C4 H3 C5 H4 C6 C7 H5 C8 H6 C9 H7 O4 H8 C10 O5 O6".split()
EL = [n.rstrip("0123456789") for n in NAMES]
COV = {"C": 0.76, "H": 0.31, "O": 0.66}


def sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


def xyz(p):
    L = open(p).read().splitlines()
    n = int(L[0].split()[0])
    el = [l.split()[0] for l in L[2:2 + n]]
    X = np.array([[float(v) for v in l.split()[1:4]] for l in L[2:2 + n]])
    return el, X


el, XR = xyz(os.path.join(T6, "frame_42436", "reactant_qm.xyz"))
_, XT = xyz(os.path.join(T6, "frame_42436", "transition_state_qm.xyz"))
_, XP = xyz(os.path.join(T6, "frame_42436", "product_qm.xyz"))
assert el == EL, el


def bonds(X):
    return {(i, j) for i in range(24) for j in range(i + 1, 24)
            if np.linalg.norm(X[i] - X[j]) < COV[el[i]] + COV[el[j]] + 0.40}


bR, bT, bP = bonds(XR), bonds(XT), bonds(XP)
breaking, forming = sorted(bR - bP), sorted(bP - bR)
print("breaking", [(NAMES[i], NAMES[j]) for i, j in breaking])
print("forming ", [(NAMES[i], NAMES[j]) for i, j in forming])

# the grid, as context
rows = [l.rstrip("\n").split("\t") for l in open(os.path.join(T6, "sites_f42436.tsv")) if not l.startswith("#")]
h, rows = rows[0], rows[1:]
ci = {k: h.index(k) for k in ("site", "layer", "x", "y", "z", "q_min_neg", "q_min_pos")}
G = np.array([[float(r[ci["x"]]), float(r[ci["y"]]), float(r[ci["z"]])] for r in rows])
layer = [1 if r[ci["layer"]] == "cation" else 0 for r in rows]
A = np.loadtxt(os.path.join(T6, "A_f42436.tsv"), comments="#").reshape(-1)
L = np.loadtxt(os.path.join(T6, "L_f42436.tsv"), comments="#").reshape(-1)
assert len(A) == len(rows) == len(L), (len(A), len(rows), len(L))
print("grid %d sites, a from %+.3f to %+.3f, L from %+.3f to %+.3f" % (len(A), A.min(), A.max(), L.min(), L.max()))

# the seven designs, from the committed record
rec = open(os.path.join(T6, "relax", "ARRANGEMENTS.txt")).read()
blocks = re.split(r"\n  (S6-\S+)   ", rec)
designs = []
for k in range(1, len(blocks), 2):
    name, body = blocks[k], blocks[k + 1]
    sites = [(int(a), float(b), float(c), float(d)) for a, b, c, d in re.findall(
        r"site\s+(\d+)\s+q\s+([+-][\d.]+)\s+at\s+\S+\s+\S+\s+\S+\s+a\s+([+-][\d.]+)\s+L\s+([+-][\d.]+)", body)]
    m = re.search(r"score ([+-][\d.]+)", body)
    score = float(m.group(1)) if m else None
    bare = "BARE point charge" in body
    close = {g: (float(x), y) for g, x, y in re.findall(r"(R|TS): closest sphere-to-substrate ([\d.]+) A \(to (\w+)\)", body)}
    # every a and L in the record must be the committed matrices' own values
    for s, q, a_, l_ in sites:
        assert abs(a_ - A[s]) < 5e-5, (name, s, a_, A[s])
        assert abs(l_ - L[s]) < 5e-5, (name, s, l_, L[s])
    # A BARE point charge has no Lennard-Jones wall, so the site's L is NOT part of its effect. ARRANGEMENTS.txt
    # prints the site's L in the S6-P block all the same, which is true of the site and not of the arrangement;
    # nothing computed from it anywhere, since no score is reported for S6-P and its frozen effect is measured from
    # a single point whose sphere has all four LJ fields zeroed. Here L is dropped for a bare design so that the
    # per-sphere contribution shown is the one that exists.
    lj = 0.0 if bare else 1.0
    tot = sum(q * a_ + lj * l_ for _, q, a_, l_ in sites)
    if score is not None:
        assert abs(tot - score) < 2e-3, (name, tot, score)
    nearest = {}
    for g, X in (("R", XR), ("TS", XT)):
        nearest[g] = []
        for s, q, a_, l_ in sites:
            d = np.linalg.norm(X - G[s], axis=1)
            nearest[g].append([round(float(d.min()), 3), NAMES[int(d.argmin())]])
    designs.append({"name": name, "bare": bare, "score": score, "lj": lj,
                    "aq": round(sum(q * a_ for _, q, a_, _ in sites), 4),
                    "sl": round(lj * sum(l_ for *_, l_ in sites), 4),
                    "ssq": round(sum(q * q for _, q, _, _ in sites), 5),
                    "close": {g: [v[0], v[1]] for g, v in close.items()},
                    "near": nearest,
                    "sph": [[s, q, round(a_, 4), round(lj * l_, 4), round(q * a_ + lj * l_, 4)]
                            for s, q, a_, l_ in sites]})
    print("%-6s n=%2d  sum a*q %+8.3f  sum L %+7.3f  total %+8.3f  reported %s  OK%s"
          % (name, len(sites), designs[-1]["aq"], designs[-1]["sl"], tot,
             ("%.3f" % score) if score is not None else "-", "   (bare: L dropped)" if bare else ""))

# the six charged residues, from residues42436.py, and how each design's spheres sit against their charge centres.
# This is the one comparison on the page that is not internal to the optimiser, so it is computed here and checked:
# the nearest residue charge centre to every sphere, and the nearest sphere to every residue centre, at BOTH
# geometries, with whether the signs agree.
RJ = json.load(open(os.path.join(OUT, "residues.json")))
assert RJ["frame"] == "42436", RJ["frame"]
cen = {st: {r["lab"]: np.array(r[st][r["ci"]]) for r in RJ["residues"]} for st in ("R", "TS")}
for d in designs:
    d["vsres"] = {}
    for st in ("R", "TS"):
        per_sphere, per_res = [], {}
        for s in d["sph"]:
            g = G[s[0]]
            lab = min(cen[st], key=lambda k: float(np.linalg.norm(g - cen[st][k])))
            per_sphere.append([lab, round(float(np.linalg.norm(g - cen[st][lab])), 2)])
        for r in RJ["residues"]:
            if not d["sph"]:
                per_res[r["lab"]] = None
                continue
            best = min(d["sph"], key=lambda s: float(np.linalg.norm(G[s[0]] - cen[st][r["lab"]])))
            per_res[r["lab"]] = [best[0], best[1],
                                 round(float(np.linalg.norm(G[best[0]] - cen[st][r["lab"]])), 2),
                                 1 if (best[1] > 0) == (r["q"] > 0) else 0]
        d["vsres"][st] = {"sphere": per_sphere, "res": per_res}
print()
print("nearest sphere to each residue charge centre, at R:")
for r in RJ["residues"]:
    for d in designs:
        if d["name"] != "S6-20":
            continue
        v = d["vsres"]["R"]["res"][r["lab"]]
        print("   %-8s rank %d q %+d   site %-5d q %+.2f at %.2f A   %s"
              % (r["lab"], r["rank"], r["q"], v[0], v[1], v[2], "same sign" if v[3] else "opposite"))

data = {
    "res": RJ["residues"],
    "names": NAMES, "el": el,
    "R": np.round(XR, 3).tolist(), "TS": np.round(XT, 3).tolist(),
    "bR": sorted(bR), "bT": sorted(bT & bR & bP), "breaking": breaking, "forming": forming,
    "grid": [[round(float(G[i, 0]), 2), round(float(G[i, 1]), 2), round(float(G[i, 2]), 2),
              layer[i], round(float(A[i]), 2)] for i in range(len(rows))],
    "designs": designs,
    "budget": 2.0, "rmin": 3.25, "frame": "42436",
    "sha": {"A": sha(os.path.join(T6, "A_f42436.tsv"))[:16],
            "sites": sha(os.path.join(T6, "sites_f42436.tsv"))[:16],
            "L": sha(os.path.join(T6, "L_f42436.tsv"))[:16],
            "arr": sha(os.path.join(T6, "relax", "ARRANGEMENTS.txt"))[:16]},
}
os.makedirs(OUT, exist_ok=True)
json.dump(data, open(os.path.join(OUT, "data.json"), "w"), separators=(",", ":"))
print("json bytes", os.path.getsize(os.path.join(OUT, "data.json")))
