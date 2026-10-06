#!/usr/bin/env bash
# wp4c_prepare.sh - RUN ON hpc1 (login node: builds inputs only, runs no QM). WP4, the s13 rerun.
#
# WHAT IT ASKS: where does a +1 site stop against the substrate, for three site models, when it starts
# where the enzyme's own Arg90 sits? (charge-model plan WP4; CHARGE_MODEL_DECISIONS.md D2, D4, D-WP4)
#   arm a  bare +1 site (QM/MM site with zero LJ)          - the control; expected to collapse
#   arm b  +1 site with Amber's N Lennard-Jones (R* 1.824, eps 0.170) - a model property
#   arm c  methylguanidinium, all-atom MM (wp4a/wp4b)       - the validation group
# Arms a and b start at Arg90's CZ. CZ is checked, not assumed, to be the group's charge centre: with the
# AM1-BCC charges the charge centre sum(q r)/Q lies 0.11 A from CZ (Amber's own arginine charges: 0.19 A).
#
# FREEDOM - two variants, applied identically to all three arms (the s13 test2 lesson):
#   free  every substrate atom moves.
#   anch  the four ring carbons farthest from Arg90's CZ (C5 C6 C7 C8, 5.4-6.4 A; the next nearest
#         heavy atom is C9 at 4.1 A) are held, standing in for the pocket that holds the substrate in
#         the enzyme. Every atom of the contact region - O3, C2, C4, C9 and both carboxylates - is free.
#         (A first draft held s16's alignment core instead. That core contains O3, the contact atom, so
#         the O3 contact would have been fixed by construction; it was replaced.)
#   Holding substrate atoms is a modelling choice, not a literature protocol: the cluster literature
#   holds the PROTEIN and leaves chorismate free (Agbaglo et al. 2024 freeze Calpha and some Cbeta atoms
#   to "mimic the semi-rigid character of the protein tertiary structure"), because many residues hold
#   the substrate there. With one group and no pocket, the rigid model turns a free substrate away from
#   the Arg90 pose (WP4_CRITERIA.txt), so only the anchored variant can test the Arg90-like contact.
# THE SITE OR GROUP: the point sites of arms a and b are held (a design site is a fixed point). The
#   methylguanidinium of arm c is held only at its attachment - CD and HD1, which lies along Arg90's
#   CD->CG bond - and is otherwise free, as the literature holds a side chain at its backbone atoms and
#   lets it move (Behrens & Hartke 2021 hold alpha carbons; Agbaglo et al. 2024 Calpha/Cbeta).
# EMBEDDING CAVEAT: MM charges close to QM density can overpolarise it (Senn & Thiel 2009), which can
#   shorten contacts between a QM acceptor and an MM donor. Phase 1 treated Arg90 exactly this way (MM,
#   electrostatic embedding, the same LJ), so the comparison with the enzyme's contacts is like for like.
# SOURCE GEOMETRY: 05_qmmm/19_ensemble/frame_41786/reactant.pdb (NEB image 0, the optimised reactant;
#   frame 41786 has no reactant-end dip and is untouched by C4 and C5). Its 24 substrate atoms are
#   checked against the committed 05_qmmm/20_invacuo/41786_R.xyz before anything is written.
# METHOD: as s13 - B3LYP D3BJ def2-SVP def2/J RIJCOSX TightSCF QMMM Opt, QM = substrate (charge -2,
#   singlet), electrostatic embedding, ORCA 6.0.1; force field rows of the substrate copied verbatim from
#   05_qmmm/13_bridge/complex_solvated.ORCAFF.prms (rows 6208-6231), the group from wp4b's
#   mgn.ORCAFF.prms with its bonded terms renumbered.
# Heredocs run fine inside scripts on hpc1 (s13 and step08b use them).
set -euo pipefail
ROOT=${ROOT:-/home/18660916/system_development}
ORCA=/home/apps2/ORCA/6.0.1
W=$ROOT/phase2.2/wp4_s13_rerun
PDB=$ROOT/05_qmmm/19_ensemble/frame_41786/reactant.pdb
RXYZ=$ROOT/05_qmmm/20_invacuo/41786_R.xyz
PRMS=$ROOT/05_qmmm/13_bridge/complex_solvated.ORCAFF.prms
MGN=$ROOT/phase2.2/wp4_site_models/mgn.ORCAFF.prms
MGN_SHA=$(grep ' mgn.ORCAFF.prms$' "$ROOT/phase2.2/wp4_site_models/sha256_wp4b.txt" | cut -c1-64)
stop() { echo "STOP: $*"; exit 1; }
[ -e "$W" ] && stop "$W exists; remove it first"
for f in "$PDB" "$RXYZ" "$PRMS" "$MGN"; do [ -s "$f" ] || stop "missing $f"; done
[ "$(sha256sum "$MGN" | cut -c1-64)" = "$MGN_SHA" ] || stop "mgn.ORCAFF.prms differs from wp4b's checksum"
mkdir -p "$W"

python3 - "$PDB" "$RXYZ" "$PRMS" "$MGN" "$W" <<'PY'
import sys, math, os
pdb, rxyz, prms, mgnprms, W = sys.argv[1:6]
# ---- source geometry
sub, arg = [], {}
for l in open(pdb):
    if l[:6] not in ("ATOM  ", "HETATM"):
        continue
    xyz = [float(l[30:38]), float(l[38:46]), float(l[46:54])]
    if l[17:20] == "CHA" and int(l[22:26]) == 383:
        sub.append((l[76:78].strip() or l[12:16].strip()[0], l[12:16].strip(), xyz))
    if l[17:20] == "ARG" and int(l[22:26]) == 217:
        arg[l[12:16].strip()] = xyz
if len(sub) != 24:
    sys.exit("FAIL found {} substrate atoms in {}".format(len(sub), pdb))
R = [l.split() for l in open(rxyz).read().splitlines()[2:26]]
dmax = max(abs(float(R[i][k + 1]) - sub[i][2][k]) for i in range(24) for k in range(3))
if dmax > 5e-4 or [r[0] for r in R] != [s[0] for s in sub]:
    sys.exit("FAIL reactant.pdb substrate differs from 20_invacuo/41786_R.xyz (max {:.1e} A)".format(dmax))
print("substrate: 24 atoms, identical to the committed 41786_R.xyz (max |dx| {:.1e} A)".format(dmax))
# methylguanidinium on Arg90, atom order of mgn_gaff.mol2 / mgn.ORCAFF.prms
order = ["CD", "HD1", "HD2", "HD3", "NE", "HE", "CZ", "NH1", "HH11", "HH12", "NH2", "HH21", "HH22"]
cd, cg = arg["CD"], arg["CG"]
v = [cg[k] - cd[k] for k in range(3)]; n = math.sqrt(sum(x * x for x in v))
arg["HD1"] = [cd[k] + 1.09 * v[k] / n for k in range(3)]
mgn = [(o[0], o, arg[o]) for o in order]
def dist(a, b): return math.sqrt(sum((a[k] - b[k]) ** 2 for k in range(3)))
close = min(dist(a[2], b[2]) for a in mgn for b in sub)
print("methylguanidinium on Arg90: closest approach to the substrate {:.2f} A".format(close))
if close < 1.5:
    sys.exit("FAIL group overlaps the substrate")
site = arg["CZ"]
# ---- force-field rows
L = open(prms).read().splitlines()
assert L[2].strip() == "$atoms"
srows = [l.split() for l in L[4:4 + int(L[3].split()[0])] if 6208 <= int(l.split()[0]) <= 6231]
assert len(srows) == 24 and [r[1] for r in srows] == [s[0] for s in sub], "prms rows do not match the substrate"
M = open(mgnprms).read().splitlines()
sec, msec = None, {}
for l in M:
    if l.startswith("$"):
        sec = l.strip(); msec[sec] = []; continue
    msec[sec].append(l)
mrows = [l.split() for l in msec["$atoms"][1:1 + int(msec["$atoms"][0].split()[0])]]
assert [r[1] for r in mrows] == [m[0] for m in mgn], "mgn.ORCAFF.prms element order differs from the built group"
def row(i, f):
    return "%6d   %-3s%12.6f%13.6f%13.6f%13.6f%13.6f\n" % (i, f[0], f[1], f[2], f[3], f[4], f[5])
def write_prms(path, extra_atoms, bonded):
    with open(path, "w") as fh:
        fh.write("$fftype\nAMBER\n$atoms\n%d 1 4\n" % (24 + len(extra_atoms)))
        for i, r in enumerate(srows, 1):
            fh.write(row(i, (r[1], float(r[2]), float(r[3]), float(r[4]), float(r[5]), float(r[6]))))
        for j, f in enumerate(extra_atoms, 25):
            fh.write(row(j, f))
        for name, nidx, ncol in (("$bonds", 2, 2), ("$angles", 3, 2), ("$dihedrals", 4, 3)):
            rows = bonded.get(name, [])
            fh.write("%s\n%d %d %d\n" % (name, len(rows), nidx, ncol))
            for r in rows:
                fh.write("".join("%6d   " % (int(x) + 24) for x in r[:nidx]) + "   ".join(r[nidx:]) + "\n")
def mgn_bonded():
    out = {}
    for name in ("$bonds", "$angles", "$dihedrals"):
        blk = msec.get(name, ["0"])
        k = int(blk[0].split()[0])
        out[name] = [l.split() for l in blk[1:1 + k]]
    return out
SITE = {"a": ("N", 1.0, 0.0, 0.0, 0.0, 0.0), "b": ("N", 1.0, -0.170, 3.648, -0.085, 3.648)}
# anchor: the ring carbons more than 5 A from Arg90's CZ (a rule, then checked against the expected set)
names = [s_[1] for s_ in sub]
ring = ["C4", "C5", "C6", "C7", "C8", "C9"]
ANCHOR = [names.index(r) for r in ring if dist(sub[names.index(r)][2], site) > 5.0]
print("ring carbons from CZ (A): " + ", ".join("{} {:.2f}".format(r, dist(sub[names.index(r)][2], site)) for r in ring))
if [names[i] for i in ANCHOR] != ["C5", "C6", "C7", "C8"]:
    sys.exit("FAIL anchor rule gave {} instead of C5 C6 C7 C8".format([names[i] for i in ANCHOR]))
# charge centre of the group with its AM1-BCC charges, read from mgn.ORCAFF.prms
qg = [float(r[2]) for r in mrows]
cc = [sum(q * m[2][k] for q, m in zip(qg, mgn)) / sum(qg) for k in range(3)]
print("group charge centre (AM1-BCC) is {:.2f} A from CZ".format(dist(cc, site)))
sub_free = list(range(24))
sub_anch = [i for i in range(24) if i not in ANCHOR]
grp = list(range(26, 37))                       # methylguanidinium minus CD (24) and HD1 (25)
def braces(ix): return "{" + " ".join(str(i) for i in ix) + "}"
active = {("ab", "free"): braces(sub_free), ("ab", "anch"): braces(sub_anch),
          ("c", "free"): braces(sub_free + grp), ("c", "anch"): braces(sub_anch + grp)}
for arm in ("a", "b", "c"):
    for var in ("free", "anch"):
        d = os.path.join(W, "arm{}_{}".format(arm, var)); os.makedirs(d)
        with open(os.path.join(d, "combined.xyz"), "w") as fh:
            extra = [("N", site)] if arm != "c" else [(m[0], m[2]) for m in mgn]
            fh.write("{}\nframe 41786 substrate + arm {} ({})\n".format(24 + len(extra), arm, var))
            for e, nm, x in sub:
                fh.write("{:<3}{:15.8f}{:15.8f}{:15.8f}\n".format(e, *x))
            for e, x in extra:
                fh.write("{:<3}{:15.8f}{:15.8f}{:15.8f}\n".format(e, *x))
        if arm == "c":
            write_prms(os.path.join(d, "site.ORCAFF.prms"),
                       [(r[1], float(r[2]), float(r[3]), float(r[4]), float(r[5]), float(r[6])) for r in mrows], mgn_bonded())
        else:
            write_prms(os.path.join(d, "site.ORCAFF.prms"), [SITE[arm]], {})
        with open(os.path.join(d, "job.inp"), "w") as fh:
            fh.write("! B3LYP D3BJ def2-SVP def2/J RIJCOSX TightSCF QMMM Opt\n%maxcore 3000\n%pal nprocs 8 end\n"
                     "%scf MaxIter 250 end\n%geom MaxIter 500 end\n%qmmm\n  QMAtoms {0:23} end\n"
                     "  ActiveAtoms " + active[("c" if arm == "c" else "ab", var)] + " end\n  ORCAFFFilename \"site.ORCAFF.prms\"\nend\n"
                     "* xyzfile -2 1 combined.xyz\n")
open(os.path.join(W, "start_site_xyz.txt"), "w").write("Arg90 CZ (arms a, b start): {:.3f} {:.3f} {:.3f}\n".format(*site))
print("written: 6 folders (arms a, b, c x free, anch); anchor = substrate indices {} ({}); arm c group held at CD, HD1".format(
    ANCHOR, " ".join(names[i] for i in ANCHOR)))
PY

for d in "$W"/arm*; do
  n=$(basename "$d")
  { echo "#!/bin/bash"
    echo "#PBS -N cm24_wp4_$n"
    echo "#PBS -l select=1:ncpus=8:mem=24gb"
    echo "#PBS -l walltime=168:00:00"
    echo "#PBS -m ae"
    echo "#PBS -M 18660916@sun.ac.za"
    echo "#PBS -j oe"
    echo "#PBS -o $d/job.pbs.out"
    echo "cd $d"
    echo "export PATH=\"/apps/openmpi/4.1.1/bin:\$PATH\""
    echo "export LD_LIBRARY_PATH=\"$ORCA/lib:/apps/openmpi/4.1.1/lib:/apps/mambaforge/envs/medaka/lib:\${LD_LIBRARY_PATH:-}\""
    echo "export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1"
    echo "echo \"host=\$(hostname) start=\$(date)\""
    echo "$ORCA/orca job.inp > job.out 2>&1 </dev/null"
    echo "echo \"orca_exit=\$? end=\$(date)\""
    echo "grep -c 'GEOMETRY OPTIMIZATION CYCLE' job.out; grep -E 'THE OPTIMIZATION HAS CONVERGED|ORCA TERMINATED NORMALLY' job.out | tail -2"
    echo "rm -f job.gbw job*.tmp 2>/dev/null"
  } > "$d/job.pbs"
done

cat > "$W/WP4_CRITERIA.txt" <<'EOF'
WP4 criteria, fixed before submission.

WHAT WP4 IS FOR (the charge-model plan): (1) whether a +1 site carrying the LJ wall of a real atom stops
where the cheap rigid model says, so the rigid model can set the grid floor (WP9) across the grid;
(2) whether a real group placed at a design site's position reproduces the enzyme's Arg90 contact, i.e.
whether a site can be realised (D2, D4); (3) whether the bare +1 collapses under the same freedom as the
other arms (the s13 test2 lesson).

START, frame 41786, as the enzyme holds Arg90 after QM/MM optimisation: N(NE)..O2 2.71 A,
N(NH2)..O3 2.83 A, CZ..O3 3.55 A - Arg90 bridges the ether O3 and the carboxylate O2. These are the
enzyme's own QM/MM contacts at this frame, the like-for-like reference for arm c.

Every run: converged (THE OPTIMIZATION HAS CONVERGED), or, if stopped, the contact distance moves by at
most 0.05 A over the last third of the cycles (s13's rule).

arm a, bare +1, both variants (control): expected to collapse onto a substrate oxygen (closest O below
  1.2 A of the site), as s13_lj_test2's bare arm did (1.009 A, C-O split 1.345/1.215 A).
arm b, +1 with Amber N LJ:
  free - rigid-model expectation 2.584 A to O2, with O3 at 2.78 A (computed from this start, substrate
         rigid, site free; a freely translating substrate is the same motion). Expected 2.28-2.88 A.
  anch - the substrate cannot translate; only a carboxylate can turn towards the site. Expected between
         the LJ stop (about 2.6 A) and the start (3.55 A to O3, 3.64 A to O2). Reported, with the motion.
arm c, methylguanidinium (held at CD and HD1, otherwise free):
  anch - THE VALIDATION TEST. Pass if the bridge is kept (one N-H to O3, one to a carboxylate oxygen) and
         nearest N..O is within 2.62-3.11 A and CZ..nearest O within 3.18-3.86 A (Arg90 MD mean +- 2 sd,
         CHARGE_MODEL_DECISIONS.md D-WP4). Also reported: each contact's change from the enzyme's QM/MM
         value above; how far the group's charge centre moved; and the closest approach of the charge
         centre to ANY substrate atom, hydrogens included - the quantity the grid floor (WP9) constrains.
  free - informative. With the group FULLY frozen, the rigid model turns the substrate by about 230
         degrees and 5.5 A to a double-carboxylate bridge (N..O1 2.77, N..O6 2.81, CZ..O6 3.49 A;
         37 kcal/mol lower in MM). The group here can turn about its attachment, so the outcome may
         differ. Reported: binding mode and contacts against the same windows; no pass or fail.
Substrate integrity, arms b and c: no bond breaks (every bonded distance within 0.15 A of the start).
Limitation: one frame. The rigid stop varies by 0.002-0.004 A across the 30 frames (WP5), but Arg90's
contacts spread by 0.12-0.17 A (MD sd); a single frame tests the mechanism, not that spread.
EOF

echo "== frame 41786 inputs written to $W"
ls "$W"
echo "== arm c (anch) input:"; cat "$W/armc_anch/job.inp"
echo "== arm c force field, header lines:"; grep -n '^\$' -A1 "$W/armc_anch/site.ORCAFF.prms"
echo "== arm b site row:"; sed -n '29p' "$W/armb_free/site.ORCAFF.prms"
echo "To submit all six after looking these over:"
echo "  for d in $W/arm*; do qsub \$d/job.pbs; done"
