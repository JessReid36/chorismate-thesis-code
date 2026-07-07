#!/usr/bin/env bash
# Step 07b - audit and repair the H++ server output, and accept it for tleap.
#
# This combines what attempt_3 did across two steps (step 15 audit + step 15b
# chain-restore/histidine-rename) into one streamlined step, with two upgrades:
#   * expected residue/heavy-atom counts are derived from the step-07a input
#     (the template), not hard-coded, so the check is self-validating;
#   * heavy-atom drift is split into BACKBONE vs SIDE-CHAIN, because H++ performs
#     legitimate His/Asn/Gln orientation flips that move side-chain heavy atoms by
#     ~2 A. Backbone atoms must not move (that is the real integrity check); large
#     side-chain drift localised to a flipped ring/amide is expected and benign.
#
# What it does:
#   1. Restore chain IDs + canonical 1-127 numbering. The H++ output has the same
#      residues in the same order as the template; we walk both in lockstep and
#      copy the template's chain + residue number onto each H++ residue. Robust to
#      H++ stripping chains AND renumbering at once.
#   2. Rename histidines from the placed titratable ring protons: HD1 -> HID,
#      HE2 -> HIE, both -> HIP (neither -> HIE, and flagged).
#   3. Audit: chain ranges, residue/atom counts vs template, backbone/side-chain
#      drift, missing/extra heavy atoms, peptide C-N links (watch the step-02 graft
#      joins), duplicate atom names, CYS disulfide proximity, non-ff14SB resnames.
#   4. Write a clean accepted PDB (chains, canonical numbering, HID/HIE/HIP, TER)
#      for the tleap build in step 09.
#
# Hard failure (blocks tleap) only for structural-identity problems: chain range
# failure, residue-count mismatch, no hydrogens, duplicate atom names, a missing
# backbone C/N, or a missing heavy atom. Geometry issues (drift, odd C-N length,
# disulfide proximity, non-standard names) are REVIEW flags, not blocks - matching
# attempt_3's treatment of the 2.245 A flip it accepted with review.
#
# H++ run settings for the record: pH 7.0, salinity 0.15 M, internal dielectric
# 10, external dielectric 80.
set -euo pipefail

export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
       NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1

root="$HOME/system_development"
hppdir="$root/03_amber/protonation_hpp"
admin="$root/00_admin"
template="$hppdir/hpp_input_protein_only_with_TER.pdb"   # 07a output: chains + 1-127 numbering
hppout="$hppdir/hpp_output_from_server.pdb"              # raw H++ return
accepted="$hppdir/abc_protonated_hpp_accepted.pdb"       # <- written for tleap
mkdir -p "$admin"

python3 - "$template" "$hppout" "$accepted" "$admin" <<'PY'
import sys, math
from collections import defaultdict, Counter
from pathlib import Path

template, hppout, accepted, admin = (Path(a) for a in sys.argv[1:5])
his_tsv   = admin / "step07b_his_assignment.tsv"
drift_tsv = admin / "step07b_heavy_atom_drift.tsv"
summ_tsv  = admin / "step07b_audit_summary.tsv"

BACKBONE = {"N", "CA", "C", "O", "OXT"}
BB_DRIFT_FLAG = 0.02          # backbone should be ~0
SC_REPORT_TOP = 12            # how many top side-chain drifts to list
FF14SB_STD = {
    "ALA","ARG","ASN","ASP","ASH","CYS","CYX","CYM","GLN","GLU","GLH","GLY",
    "HID","HIE","HIP","ILE","LEU","LYS","LYN","MET","PHE","PRO","SER","THR",
    "TRP","TYR","VAL",
}
HIS_RING = {"ND1", "NE2", "CE1", "CD2", "CG"}   # imidazole ring => histidine

def die(msg):
    sys.exit(f"FAIL {msg}")

for p in (template, hppout):
    if not p.is_file():
        die(f"missing input: {p}")

# ---- fixed-column PDB parsing (never whitespace-split) ----------------------
def is_atom(l):   return l.startswith("ATOM  ") or l.startswith("HETATM")
def rec(l):       return l[0:6]
def aname(l):     return l[12:16].strip()
def resname(l):   return l[17:20].strip()
def chain(l):     return l[21]
def resseq(l):    return l[22:26].strip()
def icode(l):     return l[26]
def xyz(l):       return (float(l[30:38]), float(l[38:46]), float(l[46:54]))
def element(l):   return l[76:78].strip()

def is_h(l):
    e = element(l)
    if e:
        return e == "H"
    nm = aname(l).lstrip("0123456789")
    return nm[:1] == "H"

def dist(p, q):
    return math.sqrt((p[0]-q[0])**2 + (p[1]-q[1])**2 + (p[2]-q[2])**2)  # manual, pre-3.8 safe

def read_residues(path, use_chain):
    residues, cur, cur_key, ter = [], None, None, 0
    for l in path.read_text().splitlines():
        if l.startswith("TER"):
            ter += 1; cur = None; cur_key = None; continue
        if not is_atom(l):
            continue
        key = (chain(l), resseq(l), icode(l)) if use_chain else (resseq(l), icode(l), resname(l))
        if key != cur_key:
            cur = {"chain": chain(l), "resseq": resseq(l), "resname": resname(l), "lines": []}
            residues.append(cur); cur_key = key
        cur["lines"].append(l)
    return residues, ter

tmpl_res, _        = read_residues(template, use_chain=True)
hpp_res, hpp_ter   = read_residues(hppout,  use_chain=False)

if not tmpl_res: die("template has no residues")
if not hpp_res:  die("H++ output has no residues")
if len(tmpl_res) != len(hpp_res):
    die(f"residue-count mismatch: template {len(tmpl_res)} vs H++ {len(hpp_res)} "
        f"(H++ must not add/drop residues)")

def is_histidine(reslines):
    return HIS_RING.issubset({aname(l) for l in reslines})

# ---- lockstep: restore chain+number, rename His, collect restored atoms -----
restored = []                 # (chain, resid_int, resname_final, line)
his_rows = []                 # (chain, resid, input_resname, HD1, HE2, suggested)
resname_diffs = []            # (chain, resid, template_name, hpp_name)
serial = 0
for t, h in zip(tmpl_res, hpp_res):
    ch = t["chain"]; rs = int(t["resseq"])
    rn = h["resname"] or t["resname"]          # keep H++ protonation-state name for non-His

    if is_histidine(h["lines"]) or rn in {"HIS","HID","HIE","HIP","HDP"}:
        names = {aname(l) for l in h["lines"]}
        hd1, he2 = "HD1" in names, "HE2" in names
        if rn in {"HID","HIE","HIP"}:      suggested = rn
        elif hd1 and he2:                   suggested = "HIP"
        elif hd1:                           suggested = "HID"
        elif he2:                           suggested = "HIE"
        else:                               suggested = "HIE"   # neither: default, flagged below
        his_rows.append((ch, rs, h["resname"], "yes" if hd1 else "no",
                         "yes" if he2 else "no", suggested))
        rn = suggested
    else:
        if t["resname"] != h["resname"] and h["resname"]:
            resname_diffs.append((ch, rs, t["resname"], h["resname"]))

    for l in h["lines"]:
        serial += 1
        newline = (f"{rec(l)}{serial:5d} {l[12:16]}{l[16]}"
                   f"{rn:>3s} {ch}{rs:4d}{l[26]}{l[27:]}")
        restored.append((ch, rs, rn, newline))

# ---- write accepted PDB (TER between chains) --------------------------------
with accepted.open("w") as f:
    prev = None
    for ch, rs, rn, line in restored:
        if prev is not None and ch != prev:
            f.write("TER\n")
        f.write(line + "\n"); prev = ch
    f.write("TER\nEND\n")

# ---- indexes for audits -----------------------------------------------------
tmpl_heavy = {}     # (chain,resid,atom) -> xyz
for t in tmpl_res:
    for l in t["lines"]:
        if not is_h(l):
            tmpl_heavy[(t["chain"], int(t["resseq"]), aname(l))] = xyz(l)

restored_heavy = {}
restored_atoms_by_res = defaultdict(list)
xyz_by_res = defaultdict(dict)
resname_by_res = {}
restored_res_ids = defaultdict(set)
for ch, rs, rn, line in restored:
    a = aname(line)
    restored_atoms_by_res[(ch, rs)].append(a)
    resname_by_res[(ch, rs)] = rn
    restored_res_ids[ch].add(rs)
    if not is_h(line):
        restored_heavy[(ch, rs, a)] = xyz(line)
        xyz_by_res[(ch, rs)][a] = xyz(line)

# chain ranges (expected span/count derived from template)
tmpl_span = {}
for t in tmpl_res:
    ch = t["chain"]; r = int(t["resseq"])
    lo, hi = tmpl_span.get(ch, (10**9, -1))
    tmpl_span[ch] = (min(lo, r), max(hi, r))
tmpl_chain_res = Counter(t["chain"] for t in tmpl_res)

chain_rows, chain_fail = [], []
for ch in sorted(tmpl_span):
    ids = sorted(restored_res_ids.get(ch, []))
    if not ids:
        chain_rows.append((ch, "MISSING", "MISSING", 0, "FAIL")); chain_fail.append(f"chain {ch} missing"); continue
    lo, hi, n = min(ids), max(ids), len(ids)
    exp_lo, exp_hi = tmpl_span[ch]; exp_n = tmpl_chain_res[ch]
    status = "PASS" if (lo == exp_lo and hi == exp_hi and n == exp_n) else "FAIL"
    chain_rows.append((ch, lo, hi, n, status))
    if status != "PASS":
        chain_fail.append(f"chain {ch} range/count abnormal: min={lo} max={hi} count={n} (expected {exp_lo}-{exp_hi}, {exp_n})")

# drift split backbone vs side-chain, matched by (chain,resid,atom)
common = sorted(set(tmpl_heavy) & set(restored_heavy))
missing_heavy = sorted(set(tmpl_heavy) - set(restored_heavy))
extra_heavy   = sorted(set(restored_heavy) - set(tmpl_heavy))
bb_drifts, sc_drifts, sc_detail = [], [], []
for key in common:
    d = dist(tmpl_heavy[key], restored_heavy[key])
    if key[2] in BACKBONE:
        bb_drifts.append(d)
    else:
        sc_drifts.append(d); sc_detail.append((d, key[0], key[1], resname_by_res.get((key[0],key[1]),"?"), key[2]))
max_bb = max(bb_drifts) if bb_drifts else 0.0
max_sc = max(sc_drifts) if sc_drifts else 0.0
rms_bb = math.sqrt(sum(d*d for d in bb_drifts)/len(bb_drifts)) if bb_drifts else 0.0
rms_sc = math.sqrt(sum(d*d for d in sc_drifts)/len(sc_drifts)) if sc_drifts else 0.0
sc_detail.sort(reverse=True)

# peptide C-N links within each chain
bad_links, checked = [], 0
for ch in sorted(restored_res_ids):
    ids = sorted(restored_res_ids[ch])
    idset = restored_res_ids[ch]
    for r1 in ids:
        r2 = r1 + 1
        if r2 not in idset:
            continue
        C = xyz_by_res[(ch, r1)].get("C"); N = xyz_by_res[(ch, r2)].get("N")
        if C is None or N is None:
            bad_links.append((ch, r1, r2, "missing_C_or_N", None)); continue
        checked += 1
        d = dist(C, N)
        if d < 1.15 or d > 1.70:
            bad_links.append((ch, r1, r2, "abnormal_C_N_distance", d))
missing_CN  = [b for b in bad_links if b[3] == "missing_C_or_N"]
abnormal_CN = [b for b in bad_links if b[3] == "abnormal_C_N_distance"]

# duplicate atom names per residue
dup_rows = []
for key, names in restored_atoms_by_res.items():
    dups = [n for n, c in Counter(names).items() if c > 1]
    if dups:
        dup_rows.append((key, dups))

# CYS SG proximity (possible disulfide)
sg = [(ch, rs, xyz_by_res[(ch,rs)]["SG"]) for (ch,rs) in xyz_by_res
      if "SG" in xyz_by_res[(ch,rs)] and resname_by_res.get((ch,rs)) in {"CYS","CYX","CYM"}]
ss_pairs = []
for i in range(len(sg)):
    for j in range(i+1, len(sg)):
        d = dist(sg[i][2], sg[j][2])
        if d < 2.5:
            ss_pairs.append((sg[i][0], sg[i][1], sg[j][0], sg[j][1], d))

# counts
tmpl_heavy_n = len(tmpl_heavy)
restored_heavy_n = len(restored_heavy)
restored_h_n = sum(1 for ch, rs, rn, line in restored if is_h(line))
resname_counts = Counter(rn for (ch, rs), rn in resname_by_res.items())
non_std = sorted(n for n in resname_counts if n not in FF14SB_STD)
his_default_hie = [r for r in his_rows if r[2] not in {"HID","HIE","HIP"} and r[3]=="no" and r[4]=="no"]

# ---- write reports ----------------------------------------------------------
with his_tsv.open("w") as f:
    f.write("chain\tresid\tinput_resname\tHD1_placed\tHE2_placed\tsuggested_amber_resname\n")
    for r in his_rows:
        f.write("\t".join(map(str, r)) + "\n")

with drift_tsv.open("w") as f:
    f.write("# backbone drift should be ~0; large SIDE-CHAIN drift localised to a residue = H++ flip (benign)\n")
    f.write(f"summary\tmax_backbone_A={max_bb:.3f}\trms_backbone_A={rms_bb:.3f}\t"
            f"max_sidechain_A={max_sc:.3f}\trms_sidechain_A={rms_sc:.3f}\n")
    f.write("rank\tdrift_A\tchain\tresid\tresname\tatom\n")
    for i, (d, ch, rs, rn, at) in enumerate(sc_detail[:SC_REPORT_TOP], 1):
        f.write(f"{i}\t{d:.3f}\t{ch}\t{rs}\t{rn}\t{at}\n")

with summ_tsv.open("w") as f:
    f.write("metric\tvalue\n")
    for k, v in [
        ("template_residues", len(tmpl_res)), ("hpp_residues", len(hpp_res)),
        ("template_heavy_atoms", tmpl_heavy_n), ("restored_heavy_atoms", restored_heavy_n),
        ("hydrogens_added", restored_h_n), ("hpp_TER_count", hpp_ter),
        ("histidines", len(his_rows)),
        ("suggested_his", ",".join(f"{a}:{b}" for a,b in sorted(Counter(r[5] for r in his_rows).items())) or "none"),
        ("resname_counts", ",".join(f"{a}:{b}" for a,b in sorted(resname_counts.items()))),
        ("non_ff14sb_resnames", ",".join(non_std) or "none"),
        ("max_backbone_drift_A", f"{max_bb:.4f}"), ("rms_backbone_drift_A", f"{rms_bb:.4f}"),
        ("max_sidechain_drift_A", f"{max_sc:.4f}"), ("rms_sidechain_drift_A", f"{rms_sc:.4f}"),
        ("missing_heavy_atoms", len(missing_heavy)), ("extra_heavy_atoms", len(extra_heavy)),
        ("peptide_links_checked", checked), ("missing_CN_links", len(missing_CN)),
        ("abnormal_CN_links", len(abnormal_CN)), ("duplicate_atom_residues", len(dup_rows)),
        ("possible_disulfides_under_2.5A", len(ss_pairs)),
    ]:
        f.write(f"{k}\t{v}\n")

# ---- hard failures (structural identity) vs review flags --------------------
hard = []
if chain_fail:        hard.append("chain range/identity failure")
if restored_h_n == 0: hard.append("no hydrogens present after H++")
if dup_rows:          hard.append("duplicate atom names in a residue")
if missing_CN:        hard.append(f"{len(missing_CN)} peptide link(s) missing backbone C or N")
if missing_heavy:     hard.append(f"{len(missing_heavy)} original heavy atom(s) missing from H++ output")

review = []
if extra_heavy:                          review.append(f"{len(extra_heavy)} extra heavy atom(s) vs template")
if max_bb > BB_DRIFT_FLAG:               review.append(f"BACKBONE moved > {BB_DRIFT_FLAG} A (max {max_bb:.3f}) - investigate, backbone should be fixed")
if max_sc > 0.5:                         review.append(f"side-chain drift up to {max_sc:.3f} A - expected if an His/Asn/Gln was flipped (see drift table)")
if abnormal_CN:                          review.append(f"{len(abnormal_CN)} peptide C-N distance(s) out of 1.15-1.70 A (watch step-02 graft joins; minimisation may relax)")
if ss_pairs:                             review.append(f"{len(ss_pairs)} CYS SG pair(s) < 2.5 A - confirm CYS vs CYX before tleap")
if non_std:                              review.append(f"non-ff14SB resname(s): {','.join(non_std)} - map before tleap")
if any(r[2] == "HIS" for r in his_rows): review.append("H++ left histidines as HIS; accepted PDB uses the HID/HIE/HIP names")
if his_default_hie:                      review.append(f"{len(his_default_hie)} histidine(s) had neither HD1 nor HE2 - defaulted to HIE, inspect")

# ---- stdout summary ---------------------------------------------------------
print("STEP 07b - H++ output audit + repair")
print(f"  residues: template {len(tmpl_res)}  H++ {len(hpp_res)}  (matched)")
print(f"  heavy atoms: template {tmpl_heavy_n}  restored {restored_heavy_n}   + {restored_h_n} H   (TER={hpp_ter})")
print("  chain ranges:")
for ch, lo, hi, n, st in chain_rows:
    print(f"    chain {ch}: {lo}-{hi} ({n} res) {st}")
print(f"  histidines ({len(his_rows)}):")
for ch, rs, inn, hd1, he2, sug in his_rows:
    print(f"    {ch}/{rs}: {inn or '-'} -> {sug}  (HD1={hd1}, HE2={he2})")
print(f"  drift: backbone max {max_bb:.3f} A (rms {rms_bb:.3f}) | side-chain max {max_sc:.3f} A (rms {rms_sc:.3f})")
if sc_detail and sc_detail[0][0] > 0.5:
    print("    top side-chain movers (likely a flip):")
    for d, ch, rs, rn, at in sc_detail[:5]:
        print(f"      {d:.3f} A  {ch}/{rs} {rn} {at}")
print(f"  peptide links checked {checked}, missing C/N {len(missing_CN)}, abnormal C-N {len(abnormal_CN)}")
print(f"  cysteines: {resname_counts.get('CYS',0)} CYS / {resname_counts.get('CYX',0)} CYX; SG pairs <2.5A: {len(ss_pairs)}")
print(f"  resnames: {', '.join(f'{a}:{b}' for a,b in sorted(resname_counts.items()))}")

print()
if hard:
    print("  HARD FAILURES (block tleap):")
    for x in hard: print(f"    FAIL {x}")
if review:
    print("  REVIEW (non-blocking):")
    for x in review: print(f"    REVIEW {x}")
if not hard and not review:
    print("  PASS - no issues")

print(f"\n  WROTE {accepted}")
print(f"  WROTE {his_tsv}")
print(f"  WROTE {drift_tsv}")
print(f"  WROTE {summ_tsv}")

if hard:
    sys.exit("\nRESULT: HARD FAIL - do not proceed to tleap until resolved")
print("\nRESULT: PASS" + (" (with review notes above)" if review else "") + " - accepted structure ready for step 08/09")
PY
echo "STEP 07b DONE"
