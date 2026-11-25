#!/usr/bin/env python3
"""cs6_d59_prepare.py - form the inputs of STAGE 6, DIAGNOSTIC 59 from committed files.

AMENDMENT 11, A11.5 item 5, and AMENDMENT 13, A13.4 item 1 (b): the ten survivors of D56 with the largest
working-level lowering, run with functions at their positive charges and behind the helium wall on their
positive charges. This program writes the inputs; it computes nothing and reads no output.

  THE FUNCTIONS (D53's test, its forms copied):
     F   the design's charges (D56's pc.pc), and the functions of D50 (set X2, treatment A) at each
         positive charge of the design, as ghost atoms
     N   the same charges and no functions (the partner)
     GF  no charge; the same functions (the control with functions)
     GN  no charge, no functions (the plain control; one for each level and state)
     at CAM (CAM-B3LYP-D3BJ/def2-TZVPD), WXV (wB97X-V/def2-TZVPD) and HF (Hartree-Fock/def2-TZVPD, F, GF
     and GN only), R and TS, the dianion. The heads of the inputs are D53's committed heads, the ghost
     block D53's committed block, word for word.
  THE WALL (D57's W and Z forms, the helium set only):
     WHe  each positive charge a coreless centre carrying it, behind the helium set; each negative
          charge a bare point charge; dianion and monoanion, R and TS, def2-SVPD (SV) and def2-TZVPD (TZ)
     ZHe  the same centres with no charge and no point charge, the dianion, R and TS, SV and TZ

Every geometry and every position of a charge is copied from a committed file: the substrate of
d56/SV_C_ctl_<state>/substrate.xyz (the same file as D53's), the charges of D56's committed pc.pc of each
design, the heads and the ghost block of D53's committed inputs, the helium block of
d23/he_d16_S6-10_wall0_R/job.inp; the method line and the SCF block of the wall runs are D56's.

usage:
   python3 cs6_d59_prepare.py --res <results repository> --out <folder to write d59 into>
The reader cs6_d59.py calls build() and stops unless the committed d59/ is what build() gives.
"""
import argparse, csv, hashlib, math, os, sys

S6REL = "phase2.2/cs_stage6"
DIR = "d59"
HPC = "/home/18660916/system_development/phase2.2/cs_stage6/" + DIR
NDES = 10                                   # A11.5 item 5: at most the ten with the largest lowering
STATES = ("R", "TS")
FLEVELS = ("CAM", "WXV", "HF")
FKINDS = {"CAM": ("F", "N", "GF"), "WXV": ("F", "N", "GF"), "HF": ("F", "GF")}
D53_HEAD = {"CAM": "d53/DI6_CAM_F_x6_R/job.inp", "WXV": "d53/DI6_WXV_F_x6_R/job.inp", "HF": "d53/DI6_HF_F_x6_R/job.inp"}
D53_GHOST = "d53/DI6_CAM_F_x6_R/job.inp"
NB_PLAIN = 682                              # def2-TZVPD on the substrate (D53: N and GN)
NB_GHOST = 17                               # set X2 on one ghost hydrogen (D53: (784 - 682) / 6)
WLEVELS = (("SV", "def2-SVPD", 1800), ("TZ", "def2-TZVPD", 5000))
CAP = {"F": 10200, "GF": 10200, "N": 9000, "GN": 9000, "HF": 4500}   # D53's caps, by kind; HF by level
PER_JOB = {"FG": 8, "NG": 8, "HF": 14}      # runs of one job, so that every job ends within 24 h at its caps
LIMIT_S = 24 * 3600


def rd(p):
    with open(p) as f:
        return f.read()


def designs(S6):
    """A11.5 item 5, by program: the survivors of the committed d56_survivors.tsv, the most negative
    working-level ddE first, at most NDES of them"""
    rows = [r for r in csv.DictReader(open(os.path.join(S6, "d56_survivors.tsv")), delimiter="\t") if r["the_design"] == "SURVIVOR"]
    rows.sort(key=lambda r: float(r["ddE_SV_kcal"]))
    return [r["design"] for r in rows[:NDES]]


def charges_of(S6, design, state="R"):
    """[(q, ('x', 'y', 'z') as written in D56's pc.pc)] of the design's committed point charges"""
    L = [l.split() for l in rd(os.path.join(S6, "d56", "SV_%s_di_%s" % (design, state), "pc.pc")).split("\n")[1:] if l.strip()]
    return [(float(c[0]), (c[1], c[2], c[3])) for c in L]


def d53_head(S6, level):
    """the committed head of D53's F input at this level: every line before %pointcharges"""
    L = rd(os.path.join(S6, D53_HEAD[level])).split("\n")
    return L[:L.index('%pointcharges "pc.pc"')]


def d53_ghost_block(S6):
    """the lines of D53's committed ghost block that follow the line 'H: x y z' of its first ghost"""
    L = rd(os.path.join(S6, D53_GHOST)).split("\n")
    i = [k for k, l in enumerate(L) if l.startswith("H: ")][0]
    j = [k for k, l in enumerate(L) if l.startswith("H: ")][1]
    return L[i + 1:j]


def he_block(S6):
    L = rd(os.path.join(S6, "d23/he_d16_S6-10_wall0_R/job.inp")).split("\n")
    i = [k for k, l in enumerate(L) if l.split()[:1] == ["He>"]][0]
    j = i + 1
    while L[j].strip() != "end":
        j += 1
    return L[i + 1:j + 1]


def wall_head(basis):
    return "! CAM-B3LYP D3BJ %s def2/J RIJCOSX TightSCF\n%%maxcore 3000\n%%pal nprocs 8 end\n%%scf\n  MaxIter 300\n  STABPerform true\nend\n" % basis


def chunks(names, size):
    return [names[k:k + size] for k in range(0, len(names), size)]


def runs(S6):
    """[(run, job, cap, kind)] in the order the jobs make them, and the job names"""
    D = designs(S6)
    jobs = []
    # the functions: with them (controls first), without them (control first), Hartree-Fock
    for lv in ("CAM", "WXV"):
        fg = [("%s_%s_%s_%s" % (lv, k, d, st), CAP[k], k) for st in STATES for d in D for k in ("GF", "F")]
        for c in chunks(fg, PER_JOB["FG"]):
            jobs.append(c)
        ng = [x for st in STATES for x in [("%s_GN_%s" % (lv, st), CAP["GN"], "GN")] + [("%s_N_%s_%s" % (lv, d, st), CAP["N"], "N") for d in D]]
        for c in chunks(ng, PER_JOB["NG"]):
            jobs.append(c)
    hf = [x for st in STATES for x in [("HF_GN_%s" % st, CAP["HF"], "GN")] + [("HF_%s_%s_%s" % (k, d, st), CAP["HF"], k) for d in D for k in ("GF", "F")]]
    for c in chunks(hf, PER_JOB["HF"]):
        jobs.append(c)
    nf = len(jobs)
    # the wall: two designs to a job, each with its walled runs and then its wall-alone controls
    cap = dict((lv, c) for lv, _, c in WLEVELS)
    for k in range(0, len(D), 2):
        rs = []
        for d in D[k:k + 2]:
            rs += [("%s_WHe_%s_%s_%s" % (lv, d, kd, st), cap[lv], "W") for lv, _, _ in WLEVELS for kd in ("di", "mo") for st in STATES]
            rs += [("%s_ZHe_%s_di_%s" % (lv, d, st), cap[lv], "Z") for lv, _, _ in WLEVELS for st in STATES]
        jobs.append(rs)
    names = ["J%02d" % (k + 1) for k in range(nf)] + ["W%d" % (k + 1) for k in range(len(jobs) - nf)]
    out = [(r, names[i], c, kind) for i, js in enumerate(jobs) for r, c, kind in js]
    assert len(set(r for r, _, _, _ in out)) == len(out)
    for j in names:
        assert sum(c for _, jj, c, _ in out if jj == j) <= LIMIT_S, j
    return out, names


def nb_of(S6, run):
    """the number of basis functions a run of the functions asks for"""
    p = run.split("_")
    if p[1] in ("N", "GN"):
        return NB_PLAIN
    d = "_".join(p[2:-1])
    return NB_PLAIN + NB_GHOST * sum(1 for q, _ in charges_of(S6, d) if q > 0)


def build(S6):
    """every file of d59/, as {relative path: text}"""
    R, J = runs(S6)
    sub = {st: rd(os.path.join(S6, "d56/SV_C_ctl_%s/substrate.xyz" % st)) for st in STATES}
    atoms = {st: [l for l in sub[st].split("\n")[2:] if l.strip()] for st in STATES}
    assert all(len(atoms[st]) == 24 for st in atoms)
    for st in STATES:
        assert sub[st] == rd(os.path.join(S6, "d53/DI6_CAM_GN_%s/substrate.xyz" % st)), "D53's substrate is not D56's"
    heads = dict((lv, d53_head(S6, lv)) for lv in FLEVELS)
    ghost = d53_ghost_block(S6)
    he = he_block(S6)
    basis = dict((lv, b) for lv, b, _ in WLEVELS)
    F = {}
    for r, j, cap, kind in R:
        p = r.split("_")
        st = p[-1]
        if p[0] in FLEVELS:
            lv, k = p[0], p[1]
            h = "\n".join(heads[lv]) + "\n"
            if k == "GN":
                F[r + "/job.inp"] = h + "* xyzfile -2 1 substrate.xyz\n"
                F[r + "/substrate.xyz"] = sub[st]
                continue
            d = "_".join(p[2:-1])
            Q = charges_of(S6, d, st)
            pos = [x for q, x in Q if q > 0]
            if k in ("F", "N"):
                F[r + "/pc.pc"] = rd(os.path.join(S6, "d56", "SV_%s_di_%s" % (d, st), "pc.pc"))
                h += '%pointcharges "pc.pc"\n'
            if k == "N":
                F[r + "/job.inp"] = h + "* xyzfile -2 1 substrate.xyz\n"
                F[r + "/substrate.xyz"] = sub[st]
            else:
                F[r + "/job.inp"] = h + "* xyz -2 1\n" + "".join("%s\n" % a for a in atoms[st]) + \
                    "".join("H: %s %s %s\n" % x + "".join("%s\n" % b for b in ghost) for x in pos) + "*\n"
        else:
            lv, t = p[0], p[1]
            d = "_".join(p[2:-2])
            kd = p[-2]
            Q = charges_of(S6, d, st)
            if t == "WHe":
                walls, pcs = [(q, x) for q, x in Q if q > 0], [(q, x) for q, x in Q if q < 0]
            else:
                walls, pcs = [(0.0, x) for q, x in Q if q > 0], []
            chg, mult = (-2, 1) if kd == "di" else (-1, 2)
            txt = wall_head(basis[lv])
            if pcs:
                txt += '%pointcharges "pc.pc"\n'
                F[r + "/pc.pc"] = "%d\n" % len(pcs) + "".join(" %9.6f  %s  %s  %s\n" % (q, x[0], x[1], x[2]) for q, x in pcs)
            F[r + "/job.inp"] = txt + "* xyz %d %d\n" % (chg, mult) + "".join("%s\n" % a for a in atoms[st]) + "".join(
                "  He>  %.6f  %s  %s  %s  NewECP\n" % (q, x[0], x[1], x[2]) + "".join("%s\n" % b for b in he) for q, x in walls) + "*\n"
    F["RUNS.tsv"] = "run\tjob\tcap_s\tkind\tbasis_functions\n" + "".join(
        "%s\t%s\t%d\t%s\t%s\n" % (r, j, c, k, nb_of(S6, r) if r.split("_")[0] in FLEVELS else "") for r, j, c, k in R)
    F["DESIGNS.tsv"] = "design\tpositive_charges\tnegative_charges\n" + "".join(
        "%s\t%d\t%d\n" % (d, sum(1 for q, _ in charges_of(S6, d) if q > 0), sum(1 for q, _ in charges_of(S6, d) if q < 0)) for d in designs(S6))
    for j in J:
        F["run_%s.pbs" % j] = pbs(j, [(r, c) for r, jj, c, _ in R if jj == j], len(J))
    F["INPUT_CHECKSUMS.sha256"] = "".join("%s  %s\n" % (hashlib.sha256(F[p].encode()).hexdigest(), p) for p in sorted(F))
    return F


def pbs(j, mine, njobs):
    t = """#!/bin/bash
#PBS -N cm_d59_%s
#PBS -l select=1:ncpus=8:mem=48gb
#PBS -l walltime=24:00:00
#PBS -m ae
#PBS -M 18660916@sun.ac.za
#PBS -j oe
#PBS -o %s/run_%s.pbs.out
# Stage 6 DIAGNOSTIC 59: the ten survivors of D56 with the largest lowering, with functions at their positive
# charges and behind the helium wall (AMENDMENT 11, A11.5 item 5; AMENDMENT 13, A13.4 item 1 (b)).
# Job %s of %d. It makes %d runs in the order below (name, cap in s); a run that fails does not stop the ones after it.
# No -q: the allocation is automatic and walltime 24:00:00 routes to day. The job script is D56's committed run_C0.pbs
# with its names, runs and caps changed. No run is made twice:
#   submitted again, it makes only the runs whose folder holds no output yet (D44 NOTE 3).
export PATH="/apps/openmpi/4.1.1/bin:$PATH"
export LD_LIBRARY_PATH="/home/apps2/ORCA/6.0.1/lib:/apps/openmpi/4.1.1/lib:/apps/mambaforge/envs/medaka/lib:${LD_LIBRARY_PATH:-}"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
ORCA=/home/apps2/ORCA/6.0.1
W=%s
cd "$W" || exit 1
if ! mkdir .running_%s 2>/dev/null; then echo "FAIL: another job holds .running_%s in $W"; exit 1; fi
trap 'rmdir "$W/.running_%s" 2>/dev/null' EXIT
echo "host=$(hostname) start=$(date)"
echo "free memory at start:"; free -g | head -2
run() {
  d="$W/$1"
  if [ -e "$d/job.out" ]; then echo "$1 has an output already: not run again"
  else ( cd "$d" && timeout -k 10 "$2" "$ORCA/orca" job.inp > job.out 2>&1 </dev/null ); echo "$1 orca_exit=$? normal=$(grep -c 'ORCA TERMINATED NORMALLY' "$d/job.out") $(date)"; fi
}
""" % (j.lower(), HPC, j, j, njobs, len(mine), HPC, j, j, j)
    for r, c in mine:
        t += "run %s %d\n" % (r, c)
    t += 'echo "end=$(date)"\n'
    t += "n=0; for r in %s; do grep -q 'ORCA TERMINATED NORMALLY' \"$W/$r/job.out\" 2>/dev/null && n=$((n+1)); done\n" % " ".join(r for r, _ in mine)
    t += 'echo "finished: $n of %d"\n' % len(mine)
    return t


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--res", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    S6 = os.path.join(a.res, S6REL)
    F = build(S6)
    dst = os.path.join(a.out, DIR)
    if os.path.exists(dst):
        sys.exit("REFUSED: %s exists already." % dst)
    for p in sorted(F):
        os.makedirs(os.path.dirname(os.path.join(dst, p)) or dst, exist_ok=True)
        with open(os.path.join(dst, p), "w") as f:
            f.write(F[p])
    R, J = runs(S6)
    print("wrote %d files into %s: %d runs in %d jobs; designs %s" % (len(F), dst, len(R), len(J), ", ".join(designs(S6))))


if __name__ == "__main__":
    main()
