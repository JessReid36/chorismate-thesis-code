#!/usr/bin/env bash
# wp4a_mgn_am1bcc.sh - WP4, validation group (CHARGE_MODEL_DECISIONS.md D4): AM1-BCC charges for
# methylguanidinium (CH3-NH-C(NH2)2+, the arginine guanidinium group capped with a methyl at CD).
#
# Modelled on Phase 1's step08a_am1bcc_charges.sh, the protocol used for the substrate:
#   antechamber -at gaff -c bcc -nc <net charge>, charges extracted in input order and renormalised
#   to the exact integer. Same reason for running here: AM1-BCC needs sqm, which is broken on hpc1.
#   The hpc1 half (GAFF typing in read-charge mode, parmchk2, the ORCA force-field file) is wp4b.
#
# GEOMETRY SOURCE: the enzyme's own Arg90 (residue 217) in the committed QM/MM reactant,
#   chorismate-thesis-results/05_qmmm/18e_reactant_reduced/reactant_reduced.pdb.
#   Atoms CD..HH22 are taken as they are; CG is replaced by a hydrogen (HD1) 1.09 A from CD along
#   CD->CG, which caps the group as a methyl. antechamber's AM1-BCC calls sqm; its input (sqm.in,
#   recorded below) shows whether sqm optimised the geometry before the charges were taken.
#
# CHECKS (the script stops if any fails): 13 atoms and 12 bonds with expected bond lengths before
#   antechamber; afterwards 13 atoms, no Tripos dot-types, raw charge sum within 1e-3 of +1.
#   Reported, not enforced: equality of charges on topologically equivalent atoms (antechamber
#   equalises them by atomic paths by default for -c bcc) and the sign pattern against Amber's own
#   arginine charges in the committed force field (CZ +0.81, NE -0.53, NH -0.86, H +0.35/+0.45).
#
# RUN ON YOUR PC (needs the micromamba "ambertools" environment set up for step 08a).
#   bash wp4a_mgn_am1bcc.sh
set -euo pipefail

MAMBA="${MAMBA:-$HOME/bin/micromamba}"
export MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$HOME/micromamba}"
ENV="${ENV:-ambertools}"
REPO="${REPO:-$HOME/Desktop/chorismate-thesis-results}"
PDB="$REPO/05_qmmm/18e_reactant_reduced/reactant_reduced.pdb"
WD="$HOME/system_dev_offline/wp4a_mgn"

command -v "$MAMBA" >/dev/null 2>&1 || { echo "STOP: needs the conda AmberTools set up for step 08a ($MAMBA). This runs on your PC, not hpc1."; exit 1; }
[ -s "$PDB" ] || { echo "STOP: missing $PDB (is the results repo at $REPO?)"; exit 1; }
[ -e "$WD" ] && { echo "STOP: $WD exists; remove it first"; exit 1; }
mkdir -p "$WD"; cd "$WD"

echo "=== wp4a: build methylguanidinium from Arg90 (residue 217) of the committed QM/MM reactant ==="
python3 - "$PDB" mgn_input.mol2 <<'PY'
import sys, math
pdb, out = sys.argv[1], sys.argv[2]
keep = ["CD", "HD2", "HD3", "NE", "HE", "CZ", "NH1", "HH11", "HH12", "NH2", "HH21", "HH22"]
xyz = {}
for l in open(pdb):
    if l[:6] in ("ATOM  ", "HETATM") and l[17:20] == "ARG" and int(l[22:26]) == 217:
        xyz[l[12:16].strip()] = [float(l[30:38]), float(l[38:46]), float(l[46:54])]
miss = [a for a in keep + ["CG"] if a not in xyz]
if miss: sys.exit("FAIL atoms missing from residue 217: {}".format(miss))
cd, cg = xyz["CD"], xyz["CG"]
v = [cg[k] - cd[k] for k in range(3)]; n = math.sqrt(sum(x * x for x in v))
xyz["HD1"] = [cd[k] + 1.09 * v[k] / n for k in range(3)]       # methyl cap replaces CG
atoms = ["CD", "HD1", "HD2", "HD3", "NE", "HE", "CZ", "NH1", "HH11", "HH12", "NH2", "HH21", "HH22"]
sybyl = {"CD": "C.3", "NE": "N.pl3", "CZ": "C.cat", "NH1": "N.pl3", "NH2": "N.pl3"}
bonds = [("CD", "HD1", "1"), ("CD", "HD2", "1"), ("CD", "HD3", "1"), ("CD", "NE", "1"), ("NE", "HE", "1"),
         ("NE", "CZ", "1"), ("CZ", "NH1", "2"), ("CZ", "NH2", "1"), ("NH1", "HH11", "1"), ("NH1", "HH12", "1"),
         ("NH2", "HH21", "1"), ("NH2", "HH22", "1")]
# bond-length sanity before anything runs
lim = {"C-H": (1.05, 1.13), "H-N": (0.98, 1.06), "C-N": (1.28, 1.52)}   # keys in sorted element order
for a, b, _ in bonds:
    d = math.sqrt(sum((xyz[a][k] - xyz[b][k]) ** 2 for k in range(3)))
    kind = "-".join(sorted([a[0], b[0]]))
    lo, hi = lim[kind]
    if not lo <= d <= hi: sys.exit("FAIL bond {}-{} = {:.3f} A outside {}-{}".format(a, b, d, lo, hi))
idx = {a: i + 1 for i, a in enumerate(atoms)}
with open(out, "w") as fh:
    fh.write("@<TRIPOS>MOLECULE\nMGN\n{} {} 1 0 0\nSMALL\nUSER_CHARGES\n\n@<TRIPOS>ATOM\n".format(len(atoms), len(bonds)))
    for a in atoms:
        x, y, z = xyz[a]
        fh.write("{:>7d} {:<4} {:>10.4f} {:>10.4f} {:>10.4f} {:<6} 1 MGN  0.0000\n".format(idx[a], a, x, y, z, sybyl.get(a, "H")))
    fh.write("@<TRIPOS>BOND\n")
    for i, (a, b, t) in enumerate(bonds, 1):
        fh.write("{:>6d} {:>4d} {:>4d} {}\n".format(i, idx[a], idx[b], t))
print("PASS mgn_input.mol2: {} atoms, {} bonds, every bond length in range; HD1 placed on CD->CG at 1.09 A".format(len(atoms), len(bonds)))
PY

echo
echo "=== wp4a: AM1-BCC via antechamber (calls sqm) ==="
if ! "$MAMBA" run -n "$ENV" antechamber -i mgn_input.mol2 -fi mol2 -o mgn_am1bcc.mol2 -fo mol2 \
      -at gaff -c bcc -nc 1 -rn MGN > antechamber_am1bcc.log 2>&1; then
  echo "FAIL antechamber -c bcc failed"; tail -20 antechamber_am1bcc.log; exit 1
fi
[ -s mgn_am1bcc.mol2 ] || { echo "FAIL no mgn_am1bcc.mol2"; tail -20 antechamber_am1bcc.log; exit 1; }
echo "PASS AM1-BCC produced mgn_am1bcc.mol2"
[ -s sqm.in ] || { echo "FAIL no sqm.in - cannot record the AM1 settings"; exit 1; }
{ echo "sqm.in as antechamber wrote it (AM1 keywords; maxcyc absent means sqm's default, a full optimisation):"; cat sqm.in; } > sqm_keywords.txt
grep -i -E "qm_theory|maxcyc|grms_tol" sqm.in | head -3

echo
echo "=== wp4a: extract charges, renormalise to exactly +1, report ==="
amber_ver="$("$MAMBA" list -n "$ENV" ambertools 2>/dev/null | awk '$1=="ambertools"{print $2}' | head -1)"
python3 - mgn_input.mol2 mgn_am1bcc.mol2 "${amber_ver:-unknown}" "$PDB" <<'PY'
import sys, hashlib
from pathlib import Path
inp, bcc, ver, pdb = sys.argv[1:5]
def read(p):
    atoms, sec = [], None
    for l in open(p):
        if l.startswith("@<TRIPOS>"): sec = l.strip(); continue
        if sec == "@<TRIPOS>ATOM" and len(l.split()) >= 9:
            f = l.split(); atoms.append((f[1], f[5], float(f[8])))
    return atoms
src, out = read(inp), read(bcc)
if len(out) != 13: sys.exit("FAIL AM1-BCC mol2 has {} atoms, expected 13".format(len(out)))
if [a[0] for a in out] != [a[0] for a in src]: sys.exit("FAIL atom order changed")
dot = sorted({a[1] for a in out if "." in a[1]})
if dot: sys.exit("FAIL Tripos dot-types remain: {}".format(dot))
q = [a[2] for a in out]; raw = sum(q)
if abs(raw - 1.0) > 1e-3: sys.exit("FAIL raw AM1-BCC sum {:+.4f}, expected +1".format(raw))
res = 1.0 - raw; q = [x + res / len(q) for x in q]
Path("charges_mgn_am1bcc.dat").write_text("".join("{: .8f}\n".format(x) for x in q))
name = [a[0] for a in out]; typ = [a[1] for a in out]
groups = {"methyl H": ["HD1", "HD2", "HD3"], "terminal N": ["NH1", "NH2"], "terminal N-H": ["HH11", "HH12", "HH21", "HH22"]}
lines = []
for g, mem in groups.items():
    v = [q[name.index(m)] for m in mem]
    lines.append("  {:<14} spread {:.4f} e  ({})".format(g, max(v) - min(v), ", ".join("{} {:+.4f}".format(m, q[name.index(m)]) for m in mem)))
amber = {"CZ": 0.8076, "NE": -0.5295, "NH1": -0.8627, "NH2": -0.8627, "HE": 0.3456, "HH11": 0.4478, "HH12": 0.4478, "HH21": 0.4478, "HH22": 0.4478}
signs = [n for n in amber if (q[name.index(n)] > 0) != (amber[n] > 0)]
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
tab = "\n".join("  {:<5} {:<4} {:+.6f}{}".format(n, t, x, "   (Amber Arg {:+.4f})".format(amber[n]) if n in amber else "") for n, t, x in zip(name, typ, q))
prov = """wp4a - AM1-BCC charges for methylguanidinium (WP4 validation group, CHARGE_MODEL_DECISIONS.md D4)

Geometry: Arg90 (residue 217) of {pdb}, atoms CD..HH22; CG replaced by HD1 at 1.09 A along CD->CG.
Method (Phase 1 step08a protocol):
  antechamber -i mgn_input.mol2 -fi mol2 -o mgn_am1bcc.mol2 -fo mol2 -at gaff -c bcc -nc 1 -rn MGN
  AmberTools {ver}; charges in input order, renormalised to exactly +1 by spreading the {res:+.6f} e
  rounding residual over 13 atoms.

Charges (GAFF type, AM1-BCC charge; Amber's own arginine charge from the committed force field for comparison):
{tab}
  raw sum {raw:+.6f}, renormalised sum {s:+.6f}

Equivalent atoms (antechamber equalises by atomic paths by default for -c bcc):
{eq}
Sign pattern against Amber's arginine: {sg}

sha256:
  mgn_input.mol2         {h1}
  mgn_am1bcc.mol2        {h2}
  charges_mgn_am1bcc.dat {h3}
Next: wp4b on hpc1 - GAFF typing in read-charge mode, parmchk2, and the ORCA force-field file.
""".format(pdb=pdb, ver=ver, res=res, tab=tab, raw=raw, s=sum(q), eq="\n".join(lines),
           sg="all agree" if not signs else "DIFFERS at " + ", ".join(signs),
           h1=sha(inp), h2=sha(bcc), h3=sha("charges_mgn_am1bcc.dat"))
Path("wp4a_provenance.txt").write_text(prov)
print(prov)
PY

sha256sum mgn_input.mol2 mgn_am1bcc.mol2 charges_mgn_am1bcc.dat sqm.in sqm.out > sha256_wp4a.txt
echo "outputs in $WD: mgn_input.mol2 mgn_am1bcc.mol2 charges_mgn_am1bcc.dat wp4a_provenance.txt sqm_keywords.txt sqm.in sqm.out antechamber_am1bcc.log sha256_wp4a.txt"
echo "WP4A DONE - paste the provenance printed above"
