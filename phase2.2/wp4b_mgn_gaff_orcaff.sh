#!/usr/bin/env bash
# wp4b_mgn_gaff_orcaff.sh - RUN ON hpc1 (login node; antechamber in read-charge mode, parmchk2,
# tleap and orca_mm -convff are seconds each and run no QM).
#
# WP4 validation group (CHARGE_MODEL_DECISIONS.md D4), second half. Follows Phase 1's
# step08b_ligand_gaff.sh (GAFF typing with the off-HPC AM1-BCC charges, -c rc, parmchk2 -s gaff),
# step09b_tleap_build.sh (tleap: source leaprc.gaff, loadmol2, loadamberparams, saveamberparm) and
# step13a_bridge.sh (orca_mm -convff -AMBER <prmtop>). Heredocs inside a script run fine here
# (step08b uses them); only pasting them at the prompt fails.
#
# Inputs, scp'd from the PC (wp4a): ~/mgn_input.mol2 and ~/charges_mgn_am1bcc.dat, checked against
# the checksums wp4a printed. Output folder: ~/system_development/phase2.2/wp4_site_models/.
#
# AUDIT (stops on failure): 13 atoms, 12 bonds, net +1, the AM1-BCC charges preserved to 1e-4, the
# same GAFF types wp4a reported (c3 h1 nh hn cz), parmchk2 needing no guessed parameters (no ATTN),
# and every atom of the ORCA force-field file carrying the charge from wp4a and the LJ of its GAFF
# type in gaff.dat (length column = 2 R*, ORCA's convention measured in WP1).
# REPORTED: the GAFF LJ of each type next to the enzyme's own arginine (Amber N R* 1.824, eps 0.170;
# polar H 0.6, 0.0157), the D4 statement that only the charges differ.
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
ORCA=/home/apps2/ORCA/6.0.1
W=$HOME/system_development/phase2.2/wp4_site_models
SHA_MOL2=c6e426d36193edeb2829cabeb075bd51cb48bdeecc4a431c8d7aee2cf5856299
SHA_CHG=87448b9fd7d7c9b8c26eb45262a460b5e0e52f2b48302ad56530a94d5e09f2ac
stop() { echo "STOP: $*"; exit 1; }

[ -e "$W" ] && stop "$W exists; remove it first"
for f in mgn_input.mol2 charges_mgn_am1bcc.dat; do [ -s "$HOME/$f" ] || stop "~/$f missing - scp it from the PC first"; done
[ "$(sha256sum "$HOME/mgn_input.mol2" | cut -c1-64)" = "$SHA_MOL2" ] || stop "mgn_input.mol2 differs from wp4a's"
[ "$(sha256sum "$HOME/charges_mgn_am1bcc.dat" | cut -c1-64)" = "$SHA_CHG" ] || stop "charges file differs from wp4a's"
mkdir -p "$W"; cp -p "$HOME/mgn_input.mol2" "$HOME/charges_mgn_am1bcc.dat" "$W/"; cd "$W"
echo "inputs match wp4a's checksums"

set +u; export PERL5LIB="${PERL5LIB:-}"; module load app/amber22/22; set -u
for exe in antechamber parmchk2 tleap; do command -v "$exe" >/dev/null || stop "missing $exe"; done
echo "AMBERHOME=$AMBERHOME"

echo "=== GAFF typing with the AM1-BCC charges (read-charge mode, no sqm) ==="
antechamber -i mgn_input.mol2 -fi mol2 -o mgn_gaff.mol2 -fo mol2 -at gaff -c rc -cf charges_mgn_am1bcc.dat \
            -nc 1 -rn MGN > antechamber_mgn_gaff.log 2>&1 || { tail -20 antechamber_mgn_gaff.log; stop "antechamber failed"; }
parmchk2 -i mgn_gaff.mol2 -f mol2 -o mgn.frcmod -s gaff > parmchk2_mgn.log 2>&1 || { tail -20 parmchk2_mgn.log; stop "parmchk2 failed"; }
[ -s mgn.frcmod ] || stop "parmchk2 wrote no frcmod"

echo "=== tleap: standalone group topology ==="
{ echo "source leaprc.gaff"
  echo "loadamberparams mgn.frcmod"
  echo "MGN = loadmol2 mgn_gaff.mol2"
  echo "check MGN"
  echo "charge MGN"
  echo "saveamberparm MGN mgn.prmtop mgn.inpcrd"
  echo "quit"; } > tleap_mgn.in
tleap -f tleap_mgn.in > tleap_mgn.log 2>&1 || { tail -20 tleap_mgn.log; stop "tleap failed"; }
[ -s mgn.prmtop ] || { tail -20 tleap_mgn.log; stop "tleap wrote no prmtop"; }

echo "=== ORCA force-field file ==="
export LD_LIBRARY_PATH="$ORCA/lib:${LD_LIBRARY_PATH:-}"
"$ORCA/orca_mm" -convff -AMBER mgn.prmtop > convff_mgn.log 2>&1 || { tail -20 convff_mgn.log; stop "orca_mm -convff failed"; }
[ -s mgn.ORCAFF.prms ] || { tail -20 convff_mgn.log; stop "no mgn.ORCAFF.prms"; }

echo "=== audit ==="
python3 - mgn_input.mol2 mgn_gaff.mol2 charges_mgn_am1bcc.dat mgn.frcmod mgn.ORCAFF.prms "$AMBERHOME/dat/leap/parm/gaff.dat" tleap_mgn.log <<'PY'
import sys
inp, gaff, chg, frc, prms, gaffdat, tlog = sys.argv[1:8]
def mol2(p):
    atoms, nb, sec = [], 0, None
    for l in open(p):
        if l.startswith("@<TRIPOS>"): sec = l.strip(); continue
        if sec == "@<TRIPOS>ATOM" and len(l.split()) >= 9:
            f = l.split(); atoms.append((f[1], f[5], float(f[8])))
        elif sec == "@<TRIPOS>BOND" and len(l.split()) >= 4: nb += 1
    return atoms, nb
fails = []
src, _ = mol2(inp); out, nb = mol2(gaff)
q = [float(x) for x in open(chg) if x.strip()]
if len(out) != 13 or nb != 12: fails.append("mgn_gaff.mol2 has {} atoms, {} bonds".format(len(out), nb))
if [a[0] for a in out] != [a[0] for a in src]: fails.append("atom order changed")
if abs(sum(a[2] for a in out) - 1.0) > 1e-4: fails.append("net charge {:+.5f}".format(sum(a[2] for a in out)))
dq = max(abs(a[2] - x) for a, x in zip(out, q))
if dq > 1e-4: fails.append("charges not preserved (max diff {:.1e})".format(dq))
types = [a[1] for a in out]
want = {"CD": "c3", "HD1": "h1", "HD2": "h1", "HD3": "h1", "NE": "nh", "HE": "hn", "CZ": "cz",
        "NH1": "nh", "HH11": "hn", "HH12": "hn", "NH2": "nh", "HH21": "hn", "HH22": "hn"}
bad = [a[0] for a in out if want.get(a[0]) != a[1]]
if bad: fails.append("GAFF types differ from wp4a's at {}".format(bad))
ftxt = open(frc, errors="replace").read()
if "ATTN" in ftxt: fails.append("parmchk2 flagged parameters needing attention (ATTN)")
# GAFF LJ from gaff.dat (MOD4 RE block)
lj, inmod = {}, False
for l in open(gaffdat, errors="replace"):
    if l.startswith("MOD4"): inmod = True; continue
    if inmod:
        if not l.strip() or l.startswith("END"): break
        f = l.split()
        if len(f) >= 3:
            try: lj[f[0]] = (float(f[1]), float(f[2]))
            except ValueError: pass
for t in set(types):
    if t not in lj: fails.append("GAFF type {} has no LJ entry in gaff.dat".format(t))
# ORCA force-field file
L = open(prms).read().splitlines()
i = L.index("$atoms"); n = int(L[i + 1].split()[0]); rows = [l.split() for l in L[i + 2:i + 2 + n]]
if n != 13: fails.append("ORCAFF file has {} atoms".format(n))
worst_q = worst_r = worst_e = 0.0
for (name, t, _), x, r in zip(out, q, rows):
    worst_q = max(worst_q, abs(float(r[2]) - x))
    rs, ep = lj[t]
    worst_r = max(worst_r, abs(float(r[4]) - 2 * rs)); worst_e = max(worst_e, abs(abs(float(r[3])) - ep))
if worst_q > 1e-4: fails.append("ORCAFF charges differ from wp4a's by up to {:.1e}".format(worst_q))
if worst_r > 1e-4 or worst_e > 1e-4: fails.append("ORCAFF LJ differs from gaff.dat (R {:.1e}, eps {:.1e})".format(worst_r, worst_e))
tl = open(tlog, errors="replace").read()
print("mgn_gaff.mol2: {} atoms, {} bonds, net {:+.5f}; charges preserved to {:.1e}; types {}".format(
    len(out), nb, sum(a[2] for a in out), dq, " ".join(sorted(set(types)))))
print("parmchk2: {} ATTN; frcmod {} lines".format(ftxt.count("ATTN"), len(ftxt.splitlines())))
print("tleap: {}".format("; ".join(l.strip() for l in tl.splitlines() if "Total unperturbed charge" in l or "Checking" in l or "Unit is OK" in l or "Warning" in l)[:300] or "see tleap_mgn.log"))
print("ORCAFF file: 13 atoms; charges match wp4a to {:.1e}; LJ match gaff.dat to {:.1e} (R) and {:.1e} (eps)".format(worst_q, worst_r, worst_e))
print("\nGAFF Lennard-Jones by type (R* A, eps kcal/mol), against the enzyme's own arginine (Amber):")
amber = {"nh": (1.824, 0.170, "Amber N"), "hn": (0.600, 0.0157, "Amber H on N"), "cz": (1.908, 0.086, "Amber CZ (C)"),
         "c3": (1.908, 0.1094, "Amber CD (CT)"), "h1": (1.387, 0.0157, "Amber HD (H1)")}
for t in ("nh", "hn", "cz", "c3", "h1"):
    a = amber[t]
    print("  {:<3} GAFF {:.4f} {:.4f}   {:<14} {:.4f} {:.4f}   {}".format(t, lj[t][0], lj[t][1], a[2], a[0], a[1],
          "same" if abs(lj[t][0] - a[0]) < 1e-3 and abs(lj[t][1] - a[1]) < 1e-4 else "DIFFERENT"))
print("\nRESULT: " + ("all checks PASS" if not fails else "FAIL - " + "; ".join(fails)))
sys.exit(1 if fails else 0)
PY

sha256sum mgn_input.mol2 charges_mgn_am1bcc.dat mgn_gaff.mol2 mgn.frcmod mgn.prmtop mgn.inpcrd mgn.ORCAFF.prms > sha256_wp4b.txt
echo "outputs in $W:"; ls "$W"
echo "WP4B DONE - paste everything from '=== audit ===' down"
