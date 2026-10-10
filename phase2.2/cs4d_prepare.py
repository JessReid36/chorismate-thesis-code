#!/usr/bin/env python3
"""
cs4d_prepare.py - Stage 4d inputs (positive and negative controls; STAGE4D_CRITERIA.txt). RUN ON hpc1 (login node;
builds inputs only, runs no QM). Python 3.6, standard library.

Reads, each checked first against the sha256 listed in STAGING/cs4d_sources.sha256 (built from the results repo):
  05_qmmm/19_ensemble/frame_41786/reactant.pdb               NEB image 0, full system (env_R)
  05_qmmm/19_ensemble/frame_41786/neb_NEB-CI_converged.xyz   NEB climbing image, full system (env_TS)
  05_qmmm/19_ensemble/frame_41786/complex_solvated.ORCAFF.prms   Phase 1 charges and LJ parameters
  STAGING/41786_R.xyz, 41786_TS.xyz   the substrate geometries (copies of 05_qmmm/20_invacuo's, Phase 1's in vacuo inputs)
  phase2.2/cs_stage4c/L1_E*_{R,TS}/{combined.xyz,site.ORCAFF.prms,job.out}   for the LJ self-test
  STAGING/grid_v2.tsv, STAGING/transform_41786.txt           grid site 37, taken back to frame 41786's own coordinates
Stops on any failed check. Writes OUT (default ROOT/phase2.2/cs_stage4d, which must not exist):
  <run>/job.inp, <run>/substrate.xyz (42 runs), pc/<env>.pc, runs.tsv, lj.tsv, residues.tsv, batch_{a,b,c}.pbs,
  GENERATOR_LOG.txt. Output is deterministic (no dates), so a second run into another folder must be identical.
USAGE  python3 cs4d_prepare.py STAGING [OUT]       (add --check to verify the sources and stop before writing)
"""
import hashlib, lzma, math, os, re, sys
ROOT = os.environ.get("ROOT", "/home/18660916/system_development")
ARGS = [a for a in sys.argv[1:] if a != "--check"]
STAGING = ARGS[0]
OUT = ARGS[1] if len(ARGS) > 1 else os.path.join(ROOT, "phase2.2/cs_stage4d")
F41 = "05_qmmm/19_ensemble/frame_41786/"
QM0, NQM, NATOMS = 6207, 24, 55680               # QMAtoms {6207:6230} in every Phase 1 input (0-based)
ORCA_KCAL = 627.508107                           # ORCA's MM energy unit, WP1_REPORT.txt
LOG = []
def say(s): LOG.append(s); print(s)
def stop(s): print("STOP: " + s); sys.exit(1)
def dist(a, b): return math.sqrt(sum((a[k] - b[k]) ** 2 for k in range(3)))
def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""): h.update(b)
    return h.hexdigest()
# ---- 1. sources
if os.path.exists(OUT) and "--check" not in sys.argv: stop(OUT + " exists")
EXP = {}
for l in open(os.path.join(STAGING, "cs4d_sources.sha256")):
    if l.strip(): EXP[l[66:].strip()] = l[:64]
def src(rel):
    p = os.path.join(STAGING, rel[8:]) if rel.startswith("STAGING/") else os.path.join(ROOT, rel)
    if rel not in EXP: stop("no expected checksum for " + rel)
    if not os.path.isfile(p): stop("missing " + p)
    if sha(p) != EXP[rel]: stop("checksum differs from the results repo's record: " + rel)
    return p
for rel in sorted(EXP): src(rel)
say("sources: %d files, each matches its sha256 in the results repo" % len(EXP))
# ---- 2. force field
def prms_atoms(p):
    L = open(p).read().splitlines(); i = [k for k, l in enumerate(L) if l.strip() == "$atoms"][0]
    n = int(L[i + 1].split()[0]); A = []
    for l in L[i + 2:i + 2 + n]:
        f = l.split()
        A.append({"el": f[1], "q": float(f[2]), "eps": abs(float(f[3])), "R": float(f[4])})
    return A
FF = prms_atoms(src(F41 + "complex_solvated.ORCAFF.prms"))
if len(FF) != NATOMS: stop("force field has %d atoms, expected %d" % (len(FF), NATOMS))
# ---- 3. geometries
PDB = []
for l in open(src(F41 + "reactant.pdb")):
    if l[:6] in ("ATOM  ", "HETATM"):
        PDB.append({"name": l[12:16].strip(), "res": l[17:20].strip(), "seq": int(l[22:26]),
                    "x": (float(l[30:38]), float(l[38:46]), float(l[46:54]))})
L = open(src(F41 + "neb_NEB-CI_converged.xyz")).read().splitlines()
if int(L[0]) != NATOMS: stop("CI xyz has %s atoms" % L[0])
CI = [(l.split()[0], tuple(float(v) for v in l.split()[1:4])) for l in L[2:2 + NATOMS]]
if len(PDB) != NATOMS: stop("reactant.pdb has %d atoms" % len(PDB))
bad = [i for i in range(NATOMS) if FF[i]["el"].lower() != CI[i][0].lower()]
if bad: stop("%d elements differ between the force field and the CI xyz (first: atom %d)" % (len(bad), bad[0]))
same = sum(1 for i in range(NATOMS) if max(abs(PDB[i]["x"][k] - CI[i][1][k]) for k in range(3)) <= 6e-4)
if same < 50000: stop("only %d atoms coincide between reactant.pdb and the CI xyz; atom order in doubt" % same)
say("force field, reactant.pdb, CI xyz: %d atoms each; elements agree; %d atoms coincide (frozen region)" % (NATOMS, same))
def xyzfile(rel):
    L = open(src(rel)).read().splitlines()
    return [(l.split()[0], tuple(float(v) for v in l.split()[1:4])) for l in L[2:2 + NQM]]
SUB = {"R": xyzfile("STAGING/41786_R.xyz"), "TS": xyzfile("STAGING/41786_TS.xyz")}
qmi = list(range(QM0, QM0 + NQM))
if any(PDB[i]["res"] != "CHA" or PDB[i]["seq"] != 383 for i in qmi): stop("atoms 6207-6230 are not CHA 383")
if [FF[i]["el"].lower() for i in qmi] != [a[0].lower() for a in SUB["R"]]: stop("substrate elements differ")
dR = max(abs(PDB[i]["x"][k] - SUB["R"][j][1][k]) for j, i in enumerate(qmi) for k in range(3))
dT = max(abs(CI[i][1][k] - SUB["TS"][j][1][k]) for j, i in enumerate(qmi) for k in range(3))
if dR > 6e-4 or dT > 1e-6: stop("substrate differs: R %.1e A, TS %.1e A" % (dR, dT))
say("substrate = atoms 6207-6230 (CHA 383): reactant.pdb matches 41786_R.xyz (%.1e A), CI xyz matches 41786_TS.xyz (%.1e A)" % (dR, dT))
SLJ = [(FF[i]["eps"], FF[i]["R"]) for i in qmi]
ENV = [i for i in range(NATOMS) if not QM0 <= i < QM0 + NQM]
# ---- 4. LJ (ORCA's rule, WP1) and its self-test against Stage 4c's L1 MM energies
def lj(sub, others):
    """sub: [(x, eps, R)] substrate; others: [(x, eps, R)]; kcal/mol, every pair, no cutoff"""
    e = 0.0
    for (xa, ea, ra) in sub:
        for (xb, eb, rb) in others:
            if ea == 0.0 or eb == 0.0: continue
            r2 = (xa[0] - xb[0]) ** 2 + (xa[1] - xb[1]) ** 2 + (xa[2] - xb[2]) ** 2
            s6 = (((ra + rb) / 2.0) ** 2 / r2) ** 3
            e += math.sqrt(ea * eb) * (s6 * s6 - 2.0 * s6)
    return e
def text(p):
    return open(p, errors="replace").read() if os.path.exists(p) else lzma.open(p + ".xz", "rt", errors="replace").read()
worst = 0.0
for k in range(1, 9):
    d = {}
    for g in ("R", "TS"):
        run = "phase2.2/cs_stage4c/L1_E%d_%s/" % (k, g)
        A = prms_atoms(src(run + "site.ORCAFF.prms"))
        X = [tuple(float(v) for v in l.split()[1:4]) for l in open(src(run + "combined.xyz")).read().splitlines()[2:2 + len(A)]]
        mm = float(re.findall(r"FINAL SINGLE POINT ENERGY \(MM\)\s+(-?\d+\.\d+)", text(src(run + "job.out")))[-1])
        d[g] = (lj([(X[i], A[i]["eps"], A[i]["R"]) for i in range(24)],
                   [(X[i], A[i]["eps"], A[i]["R"]) for i in range(24, len(A))]), mm)
        if [(A[i]["eps"], A[i]["R"]) for i in range(24)] != SLJ: stop("Stage 4c's substrate LJ parameters differ from Phase 1's")
    mine = d["TS"][0] - d["R"][0]; orca = (d["TS"][1] - d["R"][1]) * ORCA_KCAL
    worst = max(worst, abs(mine - orca))
if worst > 1e-4: stop("LJ self-test: this script's dE_LJ differs from ORCA's L1 MM energies by %.2e kcal/mol" % worst)
say("LJ self-test: Stage 4c L1 E1-E8, dE_LJ from this script against ORCA's MM energies: max |diff| %.1e kcal/mol" % worst)
# ---- 5. residues (Szefczyk et al. 2004's charged residues at this active site; system numbering, step04 report)
RES = [("R90", 217, "ARG", 1), ("R7", 134, "ARG", 1), ("E78", 205, "GLU", -1), ("R116", 243, "ARG", 1),
       ("R63p", 63, "ARG", 1), ("K60p", 60, "LYS", 1)]
RIDX, CZ = {}, {}
rows = []
for lab, seq, name, z in RES:
    # match number and name: water residue numbers can wrap past 9999 and repeat a protein residue's number
    idx = [i for i in range(NATOMS) if PDB[i]["seq"] == seq and PDB[i]["res"] == name]
    if not idx or max(idx) - min(idx) + 1 != len(idx) or idx[-1] >= QM0:
        stop("residue %s %d not found as one contiguous protein residue" % (name, seq))
    qn = sum(FF[i]["q"] for i in idx)
    if abs(qn - z) > 2e-3: stop("residue %s %d has net charge %.4f, expected %+d" % (name, seq, qn, z))
    dmin = min(dist(PDB[i]["x"], a[1]) for i in idx for a in SUB["R"])
    if dmin > 10.0: stop("%s %d is %.1f A from the substrate; residue numbering in doubt" % (name, seq, dmin))
    RIDX[lab] = idx
    cz = [i for i in idx if PDB[i]["name"] == "CZ"]
    if name == "ARG":
        if len(cz) != 1: stop("no single CZ in %s" % lab)
        CZ[lab] = PDB[cz[0]]["x"]
    rows.append("%s\t%d\t%s\t%d\t%+.4f\t%.3f\t%s" % (lab, seq, name, len(idx), qn, dmin,
                "%.3f %.3f %.3f" % CZ[lab] if lab in CZ else "-"))
if max(abs(CZ["R90"][k] - (40.832, 18.750, 42.551)[k]) for k in range(3)) > 2e-3:
    stop("Arg90 CZ is not where WP4 and Stage 3 placed it (40.832, 18.750, 42.551)")
say("residues: names and whole-residue charges as expected; Arg90 CZ at WP4's position")
# ---- 6. grid site 37 in frame 41786's own coordinates (aligned = original @ R + t)
T = [[float(v) for v in l.split()] for l in open(src("STAGING/transform_41786.txt")) if not l.startswith("#")]
Rm, t = T[:3], T[3]
G = [l.rstrip("\n").split("\t") for l in open(src("STAGING/grid_v2.tsv")) if not l.startswith("#") and not l.startswith("idx")]
if len(G) != 233 or int(G[37][0]) != 37: stop("grid_v2.tsv row 37 not found")
ga = [float(v) - t[k] for k, v in enumerate(G[37][1:4])]
G37 = tuple(sum(ga[j] * Rm[i][j] for j in range(3)) for i in range(3))     # (aligned - t) @ R.T
d37 = dist(G37, CZ["R90"])
if abs(d37 - 0.648) > 0.002: stop("grid site 37 is %.3f A from Arg90 CZ, expected 0.648" % d37)
say("grid site 37 in frame 41786's coordinates: %.6f %.6f %.6f (%.3f A from Arg90 CZ)" % (G37 + (d37,)))
if "--check" in sys.argv:
    say("check only: sources, atom order, substrate, LJ self-test, residues and grid site verified; nothing written"); sys.exit(0)
# ---- 7. environments (point-charge sets) and LJ sets
pos = {"envR": {i: PDB[i]["x"] for i in ENV}, "envTS": {i: CI[i][1] for i in ENV}}
NLJ = (0.170, 3.648)                                    # Amber N, the sphere's LJ (Stages 1-4c)
def charges(geo, scale=1.0, zero=(), extra=()):
    z = set(zero)
    return [(0.0 if i in z else scale * FF[i]["q"], pos[geo][i]) for i in ENV] + list(extra)
PC = {"NC0": [(0.0, pos["envR"][i]) for i in ENV], "ENZ": charges("envR"), "ENZTS": charges("envTS"),
      "REV": charges("envR", -1.0)}
for lab in RIDX:
    PC["KO_" + lab] = charges("envR", zero=RIDX[lab])
    PC["ALONE_" + lab] = [(FF[i]["q"], pos["envR"][i]) for i in RIDX[lab]]
PC["SPH_R90CZ"] = [(1.0, CZ["R90"])]
PC["SPH_R7CZ"] = [(1.0, CZ["R7"])]
PC["SPH_G37"] = [(1.0, G37)]
PC["SUB_R90"] = charges("envR", zero=RIDX["R90"], extra=[(1.0, CZ["R90"])])
LJSET = {"envR_all": [(pos["envR"][i], FF[i]["eps"], FF[i]["R"]) for i in ENV],
         "envTS_all": [(pos["envTS"][i], FF[i]["eps"], FF[i]["R"]) for i in ENV],
         "sph_R90CZ": [(CZ["R90"],) + NLJ], "sph_R7CZ": [(CZ["R7"],) + NLJ], "sph_G37": [(G37,) + NLJ]}
for lab in RIDX: LJSET["res_" + lab] = [(pos["envR"][i], FF[i]["eps"], FF[i]["R"]) for i in RIDX[lab]]
LJOF = {"NC0": "envR_all", "ENZ": "envR_all", "ENZTS": "envTS_all", "REV": "envR_all", "SUB_R90": "envR_all",
        "SPH_R90CZ": "sph_R90CZ", "SPH_R7CZ": "sph_R7CZ", "SPH_G37": "sph_G37"}
for lab in RIDX: LJOF["KO_" + lab] = "envR_all"; LJOF["ALONE_" + lab] = "res_" + lab
# ---- 8. write
os.makedirs(os.path.join(OUT, "pc"))
for e in sorted(PC):
    with open(os.path.join(OUT, "pc", e + ".pc"), "w") as f:
        f.write("%d\n" % len(PC[e]))
        for q, x in PC[e]: f.write("%10.6f %14.8f %14.8f %14.8f\n" % (q, x[0], x[1], x[2]))
ljrows = []
for s in sorted(LJSET):
    v = {}
    for g in ("R", "TS"):
        v[g] = lj([(SUB[g][j][1],) + SLJ[j] for j in range(NQM)], LJSET[s])
    ljrows.append("%s\t%d\t%.6f\t%.6f\t%.6f" % (s, len(LJSET[s]), v["R"], v["TS"], v["TS"] - v["R"]))
open(os.path.join(OUT, "lj.tsv"), "w").write(
    "# substrate-environment LJ, kcal/mol, ORCA's rule (WP1), every pair, no cutoff; dLJ = TS - R\n"
    "# lj_set\tn_atoms\tLJ_R\tLJ_TS\tdLJ\n" + "\n".join(ljrows) + "\n")
open(os.path.join(OUT, "residues.tsv"), "w").write(
    "# label\tresseq\tresname\tn_atoms\tnet_q\tmin_dist_to_substrate_R_A\tCZ (frame 41786 coordinates)\n" + "\n".join(rows) + "\n")
HEAD = """! B3LYP D3BJ def2-SVP def2/J RIJCOSX TightSCF
%maxcore 3000
%pal nprocs 8 end
%scf MaxIter 300 end
%output
  Print[P_Hirshfeld] 1
end
"""
SRCXYZ = {"R": src("STAGING/41786_R.xyz"), "TS": src("STAGING/41786_TS.xyz")}
ORDER = ["bare", "NC0", "ENZ", "ENZTS", "REV", "SUB_R90"] + ["KO_" + r[0] for r in RES] + \
        ["ALONE_" + r[0] for r in RES] + ["SPH_R90CZ", "SPH_R7CZ", "SPH_G37"]
BATCH = {"bare": "a", "NC0": "a", "ENZ": "a", "ENZTS": "a", "REV": "a", "SUB_R90": "a"}
runrows = []
for e in ORDER:
    for g in ("R", "TS"):
        run = "%s_%s" % (e, g); d = os.path.join(OUT, run); os.makedirs(d)
        open(os.path.join(d, "substrate.xyz"), "wb").write(open(SRCXYZ[g], "rb").read())
        open(os.path.join(d, "job.inp"), "w").write(HEAD + ('%pointcharges "env.pc"\n' if e != "bare" else "") +
                                                     "* xyzfile -2 1 substrate.xyz\n")
        runrows.append("%s\t%s\t%s\t%d\t%s\t%s" % (run, e if e != "bare" else "-", g, len(PC[e]) if e != "bare" else 0,
                       LJOF.get(e, "-"), BATCH.get(e, "b" if e.startswith("KO_") else "c")))
open(os.path.join(OUT, "runs.tsv"), "w").write(
    "# run\tpoint-charge file (pc/<name>.pc, copied to env.pc by the batch job)\tgeometry\tn_point_charges\tlj_set\tbatch\n"
    + "\n".join(runrows) + "\n")
PBS = """#!/bin/bash
#PBS -N cm34d_%s
#PBS -l select=1:ncpus=8:mem=24gb
#PBS -l walltime=168:00:00
#PBS -m ae
#PBS -M 18660916@sun.ac.za
#PBS -j oe
#PBS -o %s/batch_%s.pbs.out
export PATH="/apps/openmpi/4.1.1/bin:$PATH"
export LD_LIBRARY_PATH="/home/apps2/ORCA/6.0.1/lib:/apps/openmpi/4.1.1/lib:/apps/mambaforge/envs/medaka/lib:${LD_LIBRARY_PATH:-}"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
echo "host=$(hostname) start=$(date)"
"""
SDH = "/home/18660916/system_development/phase2.2/cs_stage4d"
for b in "abc":
    s = PBS % (b, SDH, b)
    for r in runrows:
        f = r.split("\t")
        if f[5] != b: continue
        s += "cd %s/%s\n" % (SDH, f[0])
        if f[1] != "-": s += "cp ../pc/%s.pc env.pc\n" % f[1]
        s += "/home/apps2/ORCA/6.0.1/orca job.inp > job.out 2>&1 </dev/null\n"
        s += 'echo "%s orca_exit=$? $(grep -c \'ORCA TERMINATED NORMALLY\' job.out) $(date)"\n' % f[0]
        s += "rm -f env.pc job.gbw job*.tmp 2>/dev/null\n"
    s += 'echo "end=$(date)"\n'
    open(os.path.join(OUT, "batch_%s.pbs" % b), "w").write(s)
say("written: %d runs, %d point-charge files, lj.tsv, residues.tsv, runs.tsv, batch_a/b/c.pbs" % (len(runrows), len(PC)))
open(os.path.join(OUT, "GENERATOR_LOG.txt"), "w").write("cs4d_prepare.py\n" + "\n".join(LOG) + "\n")
