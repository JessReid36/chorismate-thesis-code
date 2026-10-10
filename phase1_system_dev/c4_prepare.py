#!/usr/bin/env python3
"""
c4_prepare.py - C4 inputs (Phase 1 audit item C4; C4_CRITERIA.txt). RUN ON hpc1 (login node; builds inputs only, runs
no QM). Python 3.6, standard library.

Reads, each checked first against the sha256 listed in STAGING/c4_sources.sha256 (built from the results repo only):
  hpc1 files, per frame (ROOT-relative)
    05_qmmm/19_ensemble/frame_N/reactant.pdb           image 0 of the band: the Phase 1 reactant (full system)
    05_qmmm/19_ensemble/frame_N/frame_N_CHA2.pdb        the frame's MD structure (the input of the Phase 1 optimisation)
    05_qmmm/<band>/frame_N/neb_MEP.allxyz               the band as written to disk (10 structures, full system)
    05_qmmm/19_ensemble/frame_41786/complex_solvated.ORCAFF.prms   the force field (one file, identical for every frame)
  staging copies of committed files (STAGING/repo/<repo path>)
    05_qmmm/23_c4_reactant/c4_bands.tsv                 frames, k_low, the written energies, the runs
    05_qmmm/19_ensemble/frame_N/reactant_opt.inp        the Phase 1 reactant optimisation input (the template)
    05_qmmm/<band>/frame_N/neb.inp                      the band's input (active-atom list)
    05_qmmm/<band>/frame_N/neb_MEP.QMRegion_trj.xyz     the written band's substrate path
Writes OUT (default ROOT/05_qmmm/23_c4_reactant, which must not exist), deterministic (no dates):
  <run>/start.xyz, <run>/reopt.inp, <run>/reopt.pbs, <run>/complex_solvated.ORCAFF.prms -> ../ (35 runs)
  complex_solvated.ORCAFF.prms, c4_bands.tsv (copy), active_atoms.txt, runs.tsv,
  frames/frame_N/{image0_active.xyz,md_substrate.xyz},
  GENERATOR_LOG.txt
USAGE  python3 c4_prepare.py STAGING [OUT]       (add --check to verify the sources and stop before writing)
"""
import hashlib, math, os, re, shutil, sys
ROOT = os.environ.get("ROOT", "/home/18660916/system_development")
ARGS = [a for a in sys.argv[1:] if a != "--check"]
STAGING = ARGS[0]
OUT = ARGS[1] if len(ARGS) > 1 else os.path.join(ROOT, "05_qmmm/23_c4_reactant")
SD = "/home/18660916/system_development/05_qmmm/23_c4_reactant"   # where the jobs run (fixed in the job scripts)
QM0, NQM, NATOMS = 6207, 24, 55680
H = 627.5094740631
LOG = []
def say(s): LOG.append(s); print(s)
def stop(s): print("STOP: " + s); sys.exit(1)
def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""): h.update(b)
    return h.hexdigest()
def kabsch_rmsd(A, B):
    """RMSD after optimal superposition (Horn's quaternion method; standard library only)"""
    n = len(A)
    ca = [sum(p[k] for p in A) / n for k in range(3)]; cb = [sum(p[k] for p in B) / n for k in range(3)]
    a = [[p[k] - ca[k] for k in range(3)] for p in A]; b = [[p[k] - cb[k] for k in range(3)] for p in B]
    S = [[sum(a[i][r] * b[i][c] for i in range(n)) for c in range(3)] for r in range(3)]
    (xx, xy, xz), (yx, yy, yz), (zx, zy, zz) = S
    N = [[xx + yy + zz, yz - zy, zx - xz, xy - yx], [yz - zy, xx - yy - zz, xy + yx, zx + xz],
         [zx - xz, xy + yx, -xx + yy - zz, yz + zy], [xy - yx, zx + xz, yz + zy, -xx - yy + zz]]
    lam = max(abs(v) for row in N for v in row) * 4 + 1.0
    v = [1.0, 0.0, 0.0, 0.0]
    for _ in range(500):
        w = [sum((N[i][j] + (lam if i == j else 0.0)) * v[j] for j in range(4)) for i in range(4)]
        s = sum(x * x for x in w) ** 0.5; v = [x / s for x in w]
    emax = sum(v[i] * sum(N[i][j] * v[j] for j in range(4)) for i in range(4))
    ga = sum(x * x for p in a for x in p); gb = sum(x * x for p in b for x in p)
    return max(0.0, (ga + gb - 2.0 * emax) / n) ** 0.5
# ---- 1. sources
if os.path.exists(OUT) and "--check" not in sys.argv: stop(OUT + " exists")
EXP = {}
for l in open(os.path.join(STAGING, "c4_sources.sha256")):
    if l.strip(): EXP[l[66:].rstrip("\n")] = l[:64]
def src(rel):
    p = os.path.join(STAGING, "repo", rel[5:]) if rel.startswith("REPO/") else os.path.join(ROOT, rel)
    if rel not in EXP: stop("no expected checksum for " + rel)
    if not os.path.isfile(p): stop("missing " + p)
    return p
for rel in sorted(EXP):
    if sha(src(rel)) != EXP[rel]: stop("checksum differs from the results repo's record: " + rel)
say("sources: %d files, each matches its sha256 in the results repo" % len(EXP))
# ---- 2. band table, force field
BT = [l.rstrip("\n").split("\t") for l in open(src("REPO/05_qmmm/23_c4_reactant/c4_bands.tsv")) if not l.startswith("#")]
HD, BT = BT[0], BT[1:]
BAND = {r[0]: dict(zip(HD, r)) for r in BT}
FRAMES = [r[0] for r in BT]
if len(FRAMES) != 30: stop("c4_bands.tsv has %d frames" % len(FRAMES))
def prms_elements(p):
    L = open(p).read().splitlines(); i = [k for k, l in enumerate(L) if l.strip() == "$atoms"][0]
    n = int(L[i + 1].split()[0])
    return [l.split()[1] for l in L[i + 2:i + 2 + n]]
FFP = src("05_qmmm/19_ensemble/frame_41786/complex_solvated.ORCAFF.prms")
EL = prms_elements(FFP)
if len(EL) != NATOMS: stop("force field has %d atoms" % len(EL))
# ---- 3. inputs of Phase 1: template and active atoms
def braces(t, key):
    m = re.search(key + r"\s*\{([^}]*)\}", t)
    if not m: return None
    out = []
    for tok in m.group(1).split():
        if ":" in tok:
            a, b = tok.split(":"); out.extend(range(int(a), int(b) + 1))
        else: out.append(int(tok))
    return out
ACTIVE = None; TEMPL = {}
for f in FRAMES:
    t = open(src("REPO/05_qmmm/19_ensemble/frame_%s/reactant_opt.inp" % f)).read()
    tn = open(src("REPO/05_qmmm/%s/frame_%s/neb.inp" % (BAND[f]["band"], f))).read()
    L = t.splitlines()
    want = ["! QMMM B3LYP D3BJ def2-SVP def2/J RIJCOSX L-Opt", "%maxcore 3000", "%pal nprocs 8 end", "%scf MaxIter 200 end",
            "%geom MaxIter 2000 end", "%qmmm", "  QMAtoms {6207:6230} end"]
    if L[:7] != want or L[8:] != ['  ORCAFFFilename "complex_solvated.ORCAFF.prms"', "end",
                                   "*pdbfile -2 1 frame_%s_CHA2.pdb" % f] or not L[7].startswith("  ActiveAtoms {"):
        stop("frame %s: reactant_opt.inp is not the expected Phase 1 template" % f)
    act = braces(t, "ActiveAtoms"); actn = braces(tn, "ActiveAtoms")
    if braces(t, "QMAtoms") != list(range(QM0, QM0 + NQM)) or braces(tn, "QMAtoms") != list(range(QM0, QM0 + NQM)):
        stop("frame %s: QM atoms" % f)
    if act != actn: stop("frame %s: active atoms differ between reactant_opt.inp and the band's neb.inp" % f)
    if ACTIVE is None: ACTIVE = act
    elif act != ACTIVE: stop("frame %s: a different active-atom list" % f)
    TEMPL[f] = L
if len(ACTIVE) != 1959 or len(set(ACTIVE)) != 1959 or ACTIVE != sorted(ACTIVE) or not set(range(QM0, QM0 + NQM)) <= set(ACTIVE):
    stop("active-atom list")
AS = set(ACTIVE)
say("Phase 1 inputs: reactant_opt.inp is the same template in every frame; the same 1959 active atoms (including the"
    " 24 QM atoms) in every reactant_opt.inp and band neb.inp")
# ---- 4. per frame: geometries
def read_pdb(p):
    A = []
    for l in open(p):
        if l[:6] in ("ATOM  ", "HETATM"):
            A.append((l[12:16].strip(), l[17:20].strip(), (float(l[30:38]), float(l[38:46]), float(l[46:54]))))
    return A
def read_allxyz(p):
    """an ORCA .allxyz file: (lines, [(comment line, first atom line, natoms)]); rows are split only when asked for"""
    L = open(p).read().splitlines(); S = []; i = 0
    while i < len(L):
        s = L[i].strip()
        if not s or s.startswith(">"): i += 1; continue
        n = int(s.split()[0])
        if i + 2 + n > len(L): stop("truncated structure in " + p)
        S.append((L[i + 1], i + 2, n)); i += n + 2
    return L, S
def rows_of(L, s):
    R = [l.split() for l in L[s[1]:s[1] + s[2]]]
    if any(len(r) < 4 for r in R): stop("a short atom line")
    return R
def read_xyz_frames(p):
    L = open(p).read().splitlines(); S = []; i = 0
    while i < len(L):
        if not L[i].strip(): i += 1; continue
        n = int(L[i].split()[0]); S.append([[float(v) for v in l.split()[1:4]] for l in L[i + 2:i + 2 + n]]); i += n + 2
    return S
def fl(r): return (float(r[1]), float(r[2]), float(r[3]))
def med(x): x = sorted(x); return (x[len(x) // 2] + x[(len(x) - 1) // 2]) / 2
def dmax(a, b): return max(abs(a[k] - b[k]) for k in range(3))
TXT = {}; info = []
for f in FRAMES:
    b = BAND[f]; band = b["band"]; kl = int(b["k_low"])
    R = read_pdb(src("05_qmmm/19_ensemble/frame_%s/reactant.pdb" % f))
    M = read_pdb(src("05_qmmm/19_ensemble/frame_%s/frame_%s_CHA2.pdb" % (f, f)))
    if len(R) != NATOMS or len(M) != NATOMS: stop("frame %s: reactant.pdb %d atoms, MD structure %d" % (f, len(R), len(M)))
    for i in range(QM0, QM0 + NQM):
        if R[i][1] != "CHA" or M[i][1] != "CHA" or R[i][0] != M[i][0]: stop("frame %s: atom %d is not the same substrate atom" % (f, i))
    XL, XS = read_allxyz(src("05_qmmm/%s/frame_%s/neb_MEP.allxyz" % (band, f)))
    if len(XS) != 10 or any(s[2] != NATOMS for s in XS): stop("frame %s: neb_MEP.allxyz is not 10 structures of %d atoms" % (f, NATOMS))
    X0 = rows_of(XL, XS[0]); K = rows_of(XL, XS[kl]); CI = rows_of(XL, XS[int(b["ci"])])
    if any(X0[i][0].lower() != EL[i].lower() or K[i][0].lower() != EL[i].lower() for i in range(NATOMS)):
        stop("frame %s: allxyz elements differ from the force field" % f)
    d0 = max(dmax(fl(X0[i]), R[i][2]) for i in range(NATOMS))
    if d0 > 6e-4: stop("frame %s: allxyz image 0 differs from reactant.pdb by %.1e A" % (f, d0))
    dfix = max(dmax(fl(K[i]), R[i][2]) for i in range(NATOMS) if i not in AS)
    if dfix > 6e-4: stop("frame %s: an atom outside the active region moved in image %d (%.1e A)" % (f, kl, dfix))
    T = read_xyz_frames(src("REPO/05_qmmm/%s/frame_%s/neb_MEP.QMRegion_trj.xyz" % (band, f)))
    if len(T) != 10: stop("frame %s: written QM-region path" % f)
    dq = max(dmax(fl(K[QM0 + j]), T[kl][j]) for j in range(NQM))
    if dq > 2e-6: stop("frame %s: allxyz image %d's substrate differs from the written path (%.1e A)" % (f, kl, dq))
    m = re.search(r"\sE\s+(-?\d+\.\d+)", XS[kl][0])
    ecom = "%.1e" % abs(float(m.group(1)) - float(b["Ew_klow"])) if m else "none"
    if m and abs(float(m.group(1)) - float(b["Ew_klow"])) > 1e-6:
        stop("frame %s: allxyz image %d's energy (%s) is not the written band's (%s)" % (f, kl, m.group(1), b["Ew_klow"]))
    rmd_r = kabsch_rmsd([R[i][2] for i in range(QM0, QM0 + NQM)], [M[i][2] for i in range(QM0, QM0 + NQM)])
    rmd_k = kabsch_rmsd([fl(K[i]) for i in range(QM0, QM0 + NQM)], [M[i][2] for i in range(QM0, QM0 + NQM)])
    sq = [sum((float(K[i][k + 1]) - float(X0[i][k + 1])) ** 2 for k in range(3)) for i in ACTIVE]
    ENVA = [i for i in ACTIVE if not QM0 <= i < QM0 + NQM]
    def env_rms(A, Bs): return (sum(sum((float(A[i][k + 1]) - float(Bs[i][k + 1])) ** 2 for k in range(3)) for i in ENVA) / len(ENVA)) ** 0.5
    info.append((f, kl, d0, dfix, dq, ecom, rmd_r, rmd_k, (sum(sq) / len(sq)) ** 0.5, max(sq) ** 0.5, env_rms(X0, CI), env_rms(K, CI)))
    starts = {0: X0, kl: K}
    imgs = ([0] if "S0" in b["runs"].split(",") else []) + [kl]
    TXT[f] = {
        "image0_active": "%d\nframe %s image 0 (reactant.pdb as written in neb_MEP.allxyz), active atoms; 5th column: 0-based"
                         " atom index\n" % (len(ACTIVE), f)
                         + "".join("%s %s %s %s %d\n" % (X0[i][0], X0[i][1], X0[i][2], X0[i][3], i) for i in ACTIVE),
        "md_substrate": "%d\nframe %s MD structure (frame_%s_CHA2.pdb), atoms %d-%d\n" % (NQM, f, f, QM0, QM0 + NQM - 1)
                        + "".join("%s %.3f %.3f %.3f\n" % (EL[i], M[i][2][0], M[i][2][1], M[i][2][2]) for i in range(QM0, QM0 + NQM)),
        "start": {img: "%d\nframe %s, image %d of the written band (05_qmmm/%s/frame_%s/neb_MEP.allxyz)\n" % (NATOMS, f, img, band, f)
                       + "".join("%s %s %s %s\n" % (r[0], r[1], r[2], r[3]) for r in starts[img]) for img in imgs}}
    del XL, XS, X0, K, CI, R, M
say("geometry checks (per frame): allxyz image 0 = reactant.pdb; atoms outside the active region unmoved in image k_low;"
    " image k_low's substrate = the committed written path; atoms 6207-6230 are the substrate in both PDB files")
say("frame  k_low  img0-pdb  fixed    substrate  E(comment)  RMSD_MD: Phase1 R  k_low  active RMS k_low-0 (max)"
    "  environment RMS to CI: image 0  k_low")
for r in info:
    say("%s  %d      %.1e  %.1e  %.1e    %-10s  %.3f              %.3f  %.3f (%.3f)              %.3f    %.3f%s" % (
        r + ("   <- start above 0.25 A from its MD structure" if r[7] > 0.25 else "",)))
say("substrate RMSD to its MD structure: Phase 1 reactant median %.3f, max %.3f A; written image k_low %.3f-%.3f A (A2"
    " threshold 0.30); written image k_low against image 0, active-region RMS: max %.3f A (A3 threshold 0.20)" % (
        med(r[6] for r in info), max(r[6] for r in info), min(r[7] for r in info), max(r[7] for r in info), max(r[8] for r in info)))
say("environment (active atoms other than the substrate) RMS displacement to the climbing image: from image 0 median %.3f,"
    " max %.3f A; from k_low median %.3f, max %.3f A" % (med(r[10] for r in info), max(r[10] for r in info),
                                                        med(r[11] for r in info), max(r[11] for r in info)))
if "--check" in sys.argv:
    print("check only: nothing written"); sys.exit(0)
# ---- 5. write
os.makedirs(OUT)
shutil.copyfile(FFP, os.path.join(OUT, "complex_solvated.ORCAFF.prms"))
shutil.copyfile(src("REPO/05_qmmm/23_c4_reactant/c4_bands.tsv"), os.path.join(OUT, "c4_bands.tsv"))
open(os.path.join(OUT, "active_atoms.txt"), "w").write("\n".join(str(i) for i in ACTIVE) + "\n")
PBS = """#!/bin/bash
#PBS -N c4_{run}
#PBS -l select=1:ncpus=8:mem=120gb
#PBS -l walltime=176:00:00
#PBS -m ae
#PBS -M 18660916@sun.ac.za
#PBS -j oe
#PBS -o {sd}/{run}/reopt.pbs.out
cd {sd}/{run}
export PATH="/apps/openmpi/4.1.1/bin:$PATH"
export LD_LIBRARY_PATH="/home/apps2/ORCA/6.0.1/lib:/apps/openmpi/4.1.1/lib:/apps/mambaforge/envs/medaka/lib:${{LD_LIBRARY_PATH:-}}"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
for f in reopt.inp start.xyz complex_solvated.ORCAFF.prms; do [[ -s "$f" ]] || {{ echo "FAIL missing $f"; exit 1; }}; done
[[ -e reopt.out ]] && {{ echo "FAIL reopt.out exists: this run was already started"; exit 1; }}
[[ $(df -Pk . | awk 'NR==2{{print $4}}') -ge 3000000 ]] || {{ echo "FAIL less than 3 GB free here (the trajectory needs up to ~1.3 GB)"; exit 1; }}
echo "host=$(hostname) start=$(date)"
/home/apps2/ORCA/6.0.1/orca reopt.inp > reopt.out 2>&1 </dev/null
echo "orca_exit=$? $(date)"
grep -q "ORCA TERMINATED NORMALLY" reopt.out && grep -q "The minimization has converged" reopt.out && echo C4_REOPT_CONVERGED || echo C4_REOPT_INCOMPLETE
env -u LD_LIBRARY_PATH python3 {sd}/c4_extract.py {sd}/{run} && echo "extract=ok" || echo "extract=FAILED"
# a finished run's trajectory (55680 atoms per step) is not needed: its last geometry is reopt.xyz. A run killed by the
# walltime never reaches this line, so its trajectory stays for the continuation (C4_CRITERIA.txt, M1)
rm -f reopt_trj.dcd reopt*.tmp 2>/dev/null
echo "end=$(date)"
"""
rows = []
for f in FRAMES:
    b = BAND[f]; kl = int(b["k_low"])
    fd = os.path.join(OUT, "frames", "frame_" + f); os.makedirs(fd)
    open(os.path.join(fd, "image0_active.xyz"), "w").write(TXT[f]["image0_active"])
    open(os.path.join(fd, "md_substrate.xyz"), "w").write(TXT[f]["md_substrate"])
    runs = ([("S0", 0)] if "S0" in b["runs"].split(",") else []) + [("S2", kl)]
    for step, img in runs:
        run = "%s_%s" % (step, f); d = os.path.join(OUT, run); os.makedirs(d)
        open(os.path.join(d, "start.xyz"), "w").write(TXT[f]["start"][img])
        L = list(TEMPL[f])
        L[0] = L[0] + " TightSCF"
        L[4:5] = ["%geom", "  MaxIter 2000", "  TolMaxG 1.0e-4", "  TolRMSG 3.0e-5", "end"]
        L[-1] = "*xyzfile -2 1 start.xyz"
        open(os.path.join(d, "reopt.inp"), "w").write("\n".join(L) + "\n")
        open(os.path.join(d, "reopt.pbs"), "w").write(PBS.format(run=run, sd=SD))
        os.symlink("../complex_solvated.ORCAFF.prms", os.path.join(d, "complex_solvated.ORCAFF.prms"))
        ew = b["Ew_0"] if img == 0 else b["Ew_klow"]
        rows.append("%s\t%s\t%s\t%d\t%s\t%s" % (run, f, step, img, ew, "yes" if (step == "S0" or img == 0) else "no"))
open(os.path.join(OUT, "runs.tsv"), "w").write(
    "# run\tframe\tstep\tstart_image\tEw_start (Eh, written band)\tstep0_run\n" + "\n".join(rows) + "\n")
say("written: %d runs (%d S0, %d S2), 30 frames' image-0 active regions and MD substrates, runs.tsv, active_atoms.txt,"
    " one force-field copy (each run links to it)" % (len(rows), sum(1 for r in rows if r.split("\t")[2] == "S0"),
                                                     sum(1 for r in rows if r.split("\t")[2] == "S2")))
open(os.path.join(OUT, "GENERATOR_LOG.txt"), "w").write("\n".join(LOG) + "\n")
