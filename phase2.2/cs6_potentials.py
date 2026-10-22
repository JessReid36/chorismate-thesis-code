#!/usr/bin/env python3
"""
cs6_potentials.py - keep the absolute electrostatic potentials the A matrix throws away.

  python3 cs6_potentials.py <results_repo> <OUT>

cs5_amatrix.py runs orca_vpot at every grid site for R, TS and P and then takes differences: A = V_TS - V_R and
AP = V_P - V_R. The absolute potentials are discarded, and with them the only information in the whole pipeline
about ABSOLUTE electrostatic stabilisation rather than the barrier difference. That matters because the score
a_i q_i is built from a difference, and a difference cannot see what holds a dianion together: STAGE6_FINDINGS
_20261022.md Finding 8 shows the bare dianion is not electronically bound in vacuo and that the optimiser's designs
make it less bound, because nothing in the score rewards keeping it bound.

This script reads them back out of the committed a_matrix outputs - no new calculation of any kind - and writes
V_f42436.tsv, one row per site: V_R, V_TS, V_P in kcal/mol per e.

Checks, each a stop:
  (a) every vpot.out read out of the archive matches its sha256 in the committed
      a_matrix_outputs.MANIFEST.sha256, and holds one value per site at the points the grid records, so the row
      order is the site order;
  (b) (V_TS - V_R) * H reproduces the committed A_f42436.tsv to 1e-5 kcal/mol per e, and (V_P - V_R) * H reproduces
      AP_f42436.tsv - which is what proves these are the same numbers the A column was built from;
  (c) V_R is negative at every site, as the potential of a dianion must be.
"""
import hashlib, os, sys, tarfile
from pathlib import Path
import numpy as np

H = 627.5094740631
ROOT, OUT = sys.argv[1], sys.argv[2]
T6 = Path(ROOT) / "phase2.2" / "cs_stage6"
STATES = ("R", "TS", "P")


def stop(m):
    sys.exit("STOP: " + m)


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


arch = T6 / "a_matrix_outputs.tar.xz"
man = T6 / "a_matrix_outputs.MANIFEST.sha256"
if not arch.is_file():
    stop("%s is missing" % arch)
if not man.is_file():
    stop("%s is missing" % man)
MAN = {}
for l in man.read_text().splitlines():
    if l.strip():
        h, n = l.split(None, 1)
        MAN[n.strip()] = h
# (a) the sites the grid records, and the potentials the job wrote
rows = [l.rstrip("\n").split("\t") for l in (T6 / "sites_f42436.tsv").read_text().splitlines() if not l.startswith("#")]
head, rows = rows[0], rows[1:]
ci = {k: head.index(k) for k in ("site", "x", "y", "z")}
P = np.array([[float(r[ci["x"]]), float(r[ci["y"]]), float(r[ci["z"]])] for r in rows])
BOHR = 0.5291772083
V = {}
with tarfile.open(str(arch)) as t:
    names = t.getnames()
    for st in STATES:
        hit = [n for n in names if n.endswith("sp_42436_%s/vpot.out" % st)]
        if len(hit) != 1:
            stop("the archive holds %d vpot.out for state %s" % (len(hit), st))
        raw = t.extractfile(hit[0]).read()
        key = hit[0] if hit[0].startswith("./") else "./" + hit[0].lstrip("./")
        if key not in MAN:
            stop("%s is not listed in %s" % (key, man.name))
        got = hashlib.sha256(raw).hexdigest()
        if got != MAN[key]:
            stop("%s sha256 %s does not match %s's %s" % (key, got[:16], man.name, MAN[key][:16]))
        L = [l.split() for l in raw.decode().splitlines() if l.strip()]
        if int(L[0][0]) != len(P):
            stop("%s declares %s points, the grid has %d sites" % (hit[0], L[0][0], len(P)))
        M = np.array([[float(x) for x in r] for r in L[1:]], dtype=float)
        if M.shape != (len(P), 4):
            stop("%s is %s, expected (%d, 4)" % (hit[0], M.shape, len(P)))
        d = float(np.abs(M[:, :3] * BOHR - P).max())
        if d > 1e-6:
            stop("%s: its points differ from the grid's sites by %.1e A, so the row order is not the site order"
                 % (hit[0], d))
        V[st] = M[:, 3]
        print("   %-2s %s  sha256 matches the manifest, %d sites, points match the grid to %.1e A"
              % (st, os.path.basename(os.path.dirname(hit[0])), len(P), d))

# (b) these are the numbers the A column was built from
for name, got, ref in (("A", (V["TS"] - V["R"]) * H, T6 / "A_f42436.tsv"),
                       ("AP", (V["P"] - V["R"]) * H, T6 / "AP_f42436.tsv")):
    if not ref.is_file():
        stop("%s is missing" % ref)
    want = np.loadtxt(str(ref), comments="#").reshape(-1)
    if want.shape != got.shape:
        stop("%s has %s rows, the potentials give %s" % (ref.name, want.shape, got.shape))
    d = float(np.abs(got - want).max())
    if d > 1e-5:
        stop("%s differs from the committed %s by %.1e kcal/mol per e" % (name, ref.name, d))
    print("   (V_%s - V_R) * H reproduces the committed %s to %.1e kcal/mol per e"
          % ("TS" if name == "A" else "P", ref.name, d))

# (c) a dianion's potential is negative everywhere outside it
if V["R"].max() >= 0:
    stop("V_R is not negative at every site: its maximum is %+.6f Eh/e" % V["R"].max())
print("   V_R from %+.4f to %+.4f Eh/e, negative at every site as a dianion's must be"
      % (V["R"].min(), V["R"].max()))

os.makedirs(OUT, exist_ok=True)
out = Path(OUT) / "V_f42436.tsv"
with open(out, "w") as fh:
    fh.write("# absolute electrostatic potential of the bare substrate at every grid site, kcal/mol per e\n")
    fh.write("# read back from the committed a_matrix outputs by cs6_potentials.py; NO new calculation\n")
    fh.write("# source %s sha256 %s\n" % (arch.name, sha(arch)))
    fh.write("# sites %s sha256 %s (row i = site i)\n" % ("sites_f42436.tsv", sha(T6 / "sites_f42436.tsv")))
    fh.write("# (V_TS - V_R) reproduces A_f42436.tsv and (V_P - V_R) reproduces AP_f42436.tsv, both checked\n")
    fh.write("# the interaction energy of a charge q at site i with state X is q * V_X(i)\n")
    fh.write("site\tV_R\tV_TS\tV_P\n")
    for i in range(len(P)):
        fh.write("%s\t%.6f\t%.6f\t%.6f\n" % (rows[i][ci["site"]], V["R"][i] * H, V["TS"][i] * H, V["P"][i] * H))
print("   wrote %s" % out)
