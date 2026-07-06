#!/usr/bin/env bash
# Step 04 - place chorismate into the repaired A/B/C active sites.
#
# Method: each chorismate ligand is positioned by least-squares (Kabsch)
# superposition onto the crystallographic transition-state-analogue (TSA) pose
# of its target A/B/C active site. The CP2K ligand coordinates are registered to
# the J/K/L copy of 2CHT (see step 01b); since all 12 chains are crystallo-
# graphically equivalent copies, superposing from the J/K/L TSA onto the A/B/C
# TSA introduces no bias and places chorismate in the A/B/C site by the observed
# analogue geometry. Placement is validated by TSA-fit RMSD, ligand-centre offset
# to the target site, a protein-ligand heavy-atom contact audit, and recorded
# closest-approach distances to the defined active-site residues.
#
# NOTE: the combined PDB is written WITHOUT bond/charge records. The placed .mol2
# files remain the authoritative source of ligand connectivity and formal charge
# (chorismate net charge = -2), re-established at parameterisation (step ~12).
#
# Catalytic-contact distances are CLOSEST HEAVY-ATOM approaches on unrelaxed
# placeholder coordinates - they indicate whether each active-site residue is in
# the expected neighbourhood, not precise interaction geometry. Reported, not gated.
#
# Residue set: Arg90, Glu78, Arg7, Arg116, Tyr108. Arg63 was measured but EXCLUDED:
# it lies ~21 A from chorismate in all three sites, i.e. it is not a substrate-
# contact residue in this system (dropped on the basis of measurement, not omitted).
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

root="$HOME/system_development"
raw2cht="$root/01_inputs/structures/2cht_raw.pdb"
protein="$root/02_preparation/protein_only/abc_repaired_clean.pdb"
ligdir="$root/01_inputs/ligands"
outdir="$root/02_preparation/ligand_placement"
admin="$root/00_admin"
mkdir -p "$outdir" "$admin"

python3 - "$raw2cht" "$protein" "$ligdir" "$outdir" "$admin" <<'PY'
import sys, math
from pathlib import Path
from collections import defaultdict
import numpy as np

raw2cht, protein, ligdir, outdir, admin = (Path(sys.argv[1]), Path(sys.argv[2]),
    Path(sys.argv[3]), Path(sys.argv[4]), Path(sys.argv[5]))

transform_rep = admin  / "step04_placement_transform_report.tsv"
contact_rep   = admin  / "step04_placement_contact_report.txt"
audit_rep     = admin  / "step04_placement_audit.txt"
catalytic_rep = admin  / "step04_catalytic_contact_report.tsv"
combined_pdb  = outdir / "abc_with_chorismate_unprotonated.pdb"

# Substrate-contact residues (Arg63 excluded: measured ~21 A, not a contact).
catalytic_residues = [90, 78, 7, 116, 108]  # Arg90, Glu78, Arg7, Arg116, Tyr108

placements = [
    dict(lig="liga.mol2", out="cha_a_placed.mol2", src=("K",210), tgt=("A",203),
         lig_chain="A", lig_resid=201, resname="CHA", label="cha_a"),
    dict(lig="ligb.mol2", out="cha_b_placed.mol2", src=("J",212), tgt=("B",201),
         lig_chain="B", lig_resid=201, resname="CHA", label="cha_b"),
    dict(lig="ligc.mol2", out="cha_c_placed.mol2", src=("L",211), tgt=("C",202),
         lig_chain="C", lig_resid=201, resname="CHA", label="cha_c"),
]

def read_tsa(path):
    g = defaultdict(dict)
    for line in path.read_text().splitlines():
        if line.startswith("HETATM") and line[17:20].strip() == "TSA":
            key = (line[21], int(line[22:26])); atom = line[12:16].strip()
            g[key][atom] = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return g

def read_mol2_atoms(path):
    atoms, sec = [], None
    for line in path.read_text().splitlines():
        if line.startswith("@<TRIPOS>"): sec = line.strip(); continue
        if sec == "@<TRIPOS>ATOM" and len(line.split()) >= 6:
            p = line.split()
            atoms.append(dict(id=int(p[0]), name=p[1], x=float(p[2]), y=float(p[3]),
                              z=float(p[4]), type=p[5],
                              charge=float(p[8]) if len(p) > 8 else 0.0))
    return atoms

def read_mol2_bonds(path):
    bonds, sec = [], None
    for line in path.read_text().splitlines():
        if line.startswith("@<TRIPOS>"): sec = line.strip(); continue
        if sec == "@<TRIPOS>BOND" and line.strip():
            bonds.append(line)
    return bonds

def kabsch(P, Q):
    Pc, Qc = P.mean(0), Q.mean(0)
    C = (P - Pc).T @ (Q - Qc)
    V, S, Wt = np.linalg.svd(C)
    d = np.sign(np.linalg.det(V @ Wt))
    U = V @ np.diag([1, 1, d]) @ Wt
    rmsd = math.sqrt((((P - Pc) @ U + Qc - Q) ** 2).sum() / len(P))
    return Pc, Qc, U, rmsd

def elem_of(a):
    e = a["type"].split(".")[0]
    e = "".join(c for c in e if c.isalpha())
    return (e[:2].capitalize() if len(e) > 1 else e.upper()) if e else "X"

def is_h(name, elem): return elem.upper() == "H" or name.strip().startswith("H")
def dist(a, b): return math.sqrt(float(((a - b) ** 2).sum()))

def read_protein(path):
    out = []
    for line in path.read_text().splitlines():
        if line.startswith("ATOM  "):
            out.append(dict(atom=line[12:16].strip(), resname=line[17:20].strip(),
                chain=line[21], resid=int(line[22:26]),
                xyz=np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])]),
                element=line[76:78].strip() if len(line) >= 78 else ""))
    return out

tsa = read_tsa(raw2cht)
prot = read_protein(protein)
assert prot, "no protein atoms"

placed_all, trows, clines, alines, catrows = [], [], [], [], []

for pl in placements:
    if pl["src"] not in tsa: sys.exit(f"FAIL missing source TSA {pl['src']}")
    if pl["tgt"] not in tsa: sys.exit(f"FAIL missing target TSA {pl['tgt']}")
    src, tgt = tsa[pl["src"]], tsa[pl["tgt"]]
    common = sorted(set(src) & set(tgt))
    if len(common) < 6: sys.exit(f"FAIL too few common TSA atoms {pl['src']}->{pl['tgt']}: {len(common)}")

    P = np.array([src[n] for n in common]); Q = np.array([tgt[n] for n in common])
    Pc, Qc, U, rmsd = kabsch(P, Q)

    atoms = read_mol2_atoms(ligdir / pl["lig"])
    if not atoms: sys.exit(f"FAIL no atoms in {pl['lig']}")
    placed = []
    for a in atoms:
        q = (np.array([a["x"], a["y"], a["z"]]) - Pc) @ U + Qc
        b = dict(a); b["x"], b["y"], b["z"] = map(float, q); placed.append(b)

    bonds = read_mol2_bonds(ligdir / pl["lig"])
    with (outdir / pl["out"]).open("w") as f:
        f.write("@<TRIPOS>MOLECULE\n" + pl["resname"] + "\n")
        f.write(f"{len(placed):5d} {len(bonds):5d}     1     0     0\nSMALL\nUSER_CHARGES\n\n")
        f.write("@<TRIPOS>ATOM\n")
        for a in placed:
            f.write(f"{a['id']:7d} {a['name']:<8s} {a['x']:10.4f} {a['y']:10.4f} {a['z']:10.4f} "
                    f"{a['type']:<8s} {1:4d} {pl['resname']:<8s} {a['charge']:10.6f}\n")
        f.write("@<TRIPOS>BOND\n")
        for bl in bonds: f.write(bl + "\n")

    centre = np.array([[a["x"], a["y"], a["z"]] for a in placed]).mean(0)
    tgt_centre = np.array(list(tgt.values())).mean(0)
    offset = dist(centre, tgt_centre)

    heavy = [a for a in placed if not is_h(a["name"], elem_of(a))]
    heavy_xyz = [np.array([a["x"], a["y"], a["z"]]) for a in heavy]

    close, min_d, min_pair = {}, None, None
    for L in heavy_xyz:
        for pa in prot:
            if is_h(pa["atom"], pa["element"]): continue
            dd = dist(L, pa["xyz"])
            if min_d is None or dd < min_d:
                min_d = dd; min_pair = (pa["chain"], pa["resid"], pa["resname"], pa["atom"], dd)
            if dd <= 4.0:
                k = (pa["chain"], pa["resid"], pa["resname"])
                if k not in close or dd < close[k]: close[k] = dd

    for resid in catalytic_residues:
        res_atoms = [pa for pa in prot
                     if pa["chain"] == pl["lig_chain"] and pa["resid"] == resid
                     and not is_h(pa["atom"], pa["element"])]
        if not res_atoms:
            catrows.append((pl["label"], pl["lig_chain"], resid, "NA", "absent")); continue
        best, best_pair = None, None
        for L in heavy_xyz:
            for pa in res_atoms:
                dd = dist(L, pa["xyz"])
                if best is None or dd < best:
                    best = dd; best_pair = (pa["resname"], pa["atom"])
        catrows.append((pl["label"], pl["lig_chain"], resid, f"{best:.3f}",
                        f"{best_pair[0]}:{best_pair[1]}"))

    trows.append((pl["label"], pl["lig"], f"{pl['src'][0]}{pl['src'][1]}",
                  f"{pl['tgt'][0]}{pl['tgt'][1]}", len(common), rmsd, offset))
    clines.append(f"contacts within 4.0 A of {pl['label']}")
    clines.append("-" * 50)
    for (ch, r, rn), dd in sorted(close.items()):
        clines.append(f"{ch} {r:4d} {rn:3s}  min_dist={dd:6.3f}")
    clines.append("")
    alines += [f"{pl['label']}_atoms\t{len(placed)}",
               f"{pl['label']}_heavy\t{len(heavy)}",
               f"{pl['label']}_tsa_fit_rmsd_A\t{rmsd:.4f}",
               f"{pl['label']}_centre_offset_A\t{offset:.3f}",
               f"{pl['label']}_min_prot_lig_heavy_A\t{min_d:.3f}"]
    if min_pair:
        pc, pr, prn, pan, dd = min_pair
        alines.append(f"{pl['label']}_min_pair\tprot:{pc}:{pr}:{prn}:{pan}:{dd:.3f}")
    if min_d is not None and min_d < 0.80:
        sys.exit(f"FAIL severe overlap for {pl['label']}: {min_d:.3f} A")
    placed_all.append((pl, placed))

with transform_rep.open("w") as f:
    f.write("label\tinput\tsource_tsa\ttarget_tsa\tcommon_atoms\tfit_rmsd_A\tcentre_offset_A\n")
    for lab, lig, s, t, n, r, o in trows:
        f.write(f"{lab}\t{lig}\t{s}\t{t}\t{n}\t{r:.4f}\t{o:.3f}\n")
contact_rep.write_text("\n".join(clines) + "\n")
audit_rep.write_text("step04_placement_audit\n" + "\n".join(alines) + "\n")

with catalytic_rep.open("w") as f:
    f.write("# Arg63 excluded from set: measured ~21 A from chorismate (not a contact residue)\n")
    f.write("ligand\tchain\tresid\tclosest_heavy_A\tligand_partner_atom\n")
    for lab, ch, resid, d, partner in catrows:
        f.write(f"{lab}\t{ch}\t{resid}\t{d}\t{partner}\n")

with combined_pdb.open("w") as f:
    for line in protein.read_text().splitlines():
        if line.startswith("ATOM  ") or line.startswith("TER"):
            f.write(line + "\n")
    serial = len(prot) + 1
    for pl, placed in placed_all:
        for a in placed:
            f.write(f"HETATM{serial:5d} {a['name'][:4]:<4s} {pl['resname']:>3s} "
                    f"{pl['lig_chain']:1s}{pl['lig_resid']:4d}    "
                    f"{a['x']:8.3f}{a['y']:8.3f}{a['z']:8.3f}"
                    f"{1.00:6.2f}{0.00:6.2f}          {elem_of(a):>2s}\n")
            serial += 1
        f.write("TER\n")
    f.write("END\n")

print("placement report:")
for lab, lig, s, t, n, r, o in trows:
    print(f"  {lab}: {lig} {s}->{t} | fit RMSD {r:.4f} A | centre offset {o:.3f} A")

resname_of = {90:"ARG",78:"GLU",7:"ARG",116:"ARG",108:"TYR"}
print("\ncatalytic-residue closest heavy-atom approach (A):")
print(f"  {'site':6s} " + " ".join(f"{resname_of[r]}{r:<3d}" for r in catalytic_residues))
for pl in placements:
    row = {resid: d for lab, ch, resid, d, p in catrows if lab == pl["label"]}
    cells = " ".join(f"{row.get(r,'NA'):>6s}" for r in catalytic_residues)
    print(f"  {pl['lig_chain']:6s} {cells}")
print("\n  (Arg63 excluded: ~21 A in all sites, not a substrate-contact residue)")
print(f"\nWROTE {catalytic_rep}")
print(f"WROTE {combined_pdb}")
PY
echo "STEP 04 DONE"
