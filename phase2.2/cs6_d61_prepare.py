#!/usr/bin/env python3
"""cs6_d61_prepare.py - form the inputs of STAGE 6, DIAGNOSTIC 61 from committed files.

D59's test of the functions at the positive charges of the ten survivors (AMENDMENT 11, A11.5 item 5;
AMENDMENT 13, A13.4 item 1 (b)), made again for the runs that carry ghost atoms, because in D59 every run
with the charges and the functions together stopped in ORCA's start-up: cs6_d59_prepare.py put each ghost
centre exactly on its point charge ("Zero distance encountered between atoms 24 and charge 2"; D59 NOTE 1).

  F   the design's charges (D56's committed pc.pc) and the functions of D50 (set X2, treatment A) as a
      ghost hydrogen near each positive charge
  GF  no charge; the same functions at the same centres (the control with functions)
  at CAM (CAM-B3LYP-D3BJ/def2-TZVPD), WXV (wB97X-V/def2-TZVPD) and HF (Hartree-Fock/def2-TZVPD), R and TS,
  the dianion, for the ten designs of D59. N and GN are D59's committed runs and are not made here.

The one change from D59's inputs: each ghost centre stands 0.01 A beyond its positive charge, on the line from
the nearest substrate atom, the rule of the committed cs6_d48_prepare.py (GHOST_OFF, gpos) that D50, D53 and
D54 used; build() re-forms D53's committed ghost centres by this rule and stops unless they are D53's to the
six decimals written. Everything else is D59's: the heads and the ghost block (copied from D53's committed
inputs by the committed cs6_d59_prepare.py), the designs, the charges, the substrate, the caps.

usage:
   python3 cs6_d61_prepare.py --res <results repository> --code <code repository> --out <folder to write d61 into>
The reader cs6_d61.py calls build() and stops unless the committed d61/ is what build() gives.
"""
import argparse, glob, hashlib, importlib.util, os, sys

S6REL = "phase2.2/cs_stage6"
DIR = "d61"
HPC = "/home/18660916/system_development/phase2.2/cs_stage6/" + DIR
STATES = ("R", "TS")
GHOST_OFF = 0.01                 # A, cs6_d48_prepare.py's: "the ghost centre sits this far beyond the charge, on the line from the nearest atom"
NEAR_LO, NEAR_HI = 0.005, 0.02   # A, cs6_d53_prepare.py's check that a ghost centre lies 0.01 A from its charge
D53_CHECK = ("d53/DI6_*_F_*",)   # D53's committed inputs with functions and charges
SMOKE = "HF_F_S_n05_k2_R"        # the first run of the first job: one run with charges and functions, through ORCA first
PER_JOB = {"CAM": 8, "WXV": 8, "HF": 14}   # D59's: every job ends within 24 h at its caps
LIMIT_S = 24 * 3600


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def rd(p):
    with open(p) as f:
        return f.read()


def mods(code):
    P2 = os.path.join(code, "phase2.2")
    return load(os.path.join(P2, "cs6_d59_prepare.py"), "p59_for_d61"), load(os.path.join(P2, "cs6_d47_prepare.py"), "p47_for_d61")


def dist(p, q):
    return ((p[0] - q[0]) ** 2 + (p[1] - q[1]) ** 2 + (p[2] - q[2]) ** 2) ** 0.5


def ghost_rule(S6, p47):
    """cs6_d48_prepare.py's gpos: the charge moved GHOST_OFF away from the nearest atom of R or TS, rounded to 6 decimals"""
    G = dict((st, p47.read_geom(os.path.join(S6, p47.GEOM[st]))) for st in STATES)

    def gpos(xyz):
        d = dict((st, min((dist(xyz, g[1:]), k) for k, g in enumerate(G[st]))) for st in STATES)
        m = min((d[st][0], d[st][1], st) for st in STATES)
        X = G[m[2]][m[1]][1:]
        L = dist(X, xyz)
        return tuple(round(xyz[j] + GHOST_OFF * (xyz[j] - X[j]) / L, 6) for j in range(3))
    return gpos


def check_rule_on_d53(S6, gpos):
    """every ghost centre the rule forms from the positive charges of D53's committed inputs with functions is one
    of that input's committed ghost centres, to the six decimals written; returns the number of centres checked"""
    n = 0
    for pat in D53_CHECK:
        for r in sorted(glob.glob(os.path.join(S6, pat))):
            inp, pc = os.path.join(r, "job.inp"), os.path.join(r, "pc.pc")
            if not (os.path.isfile(inp) and os.path.isfile(pc)):
                continue
            gh = set(tuple(float(v) for v in l.split()[1:4]) for l in rd(inp).split("\n") if l.startswith("H: "))
            pos = [tuple(float(v) for v in l.split()[1:4]) for l in rd(pc).split("\n")[1:] if l.strip() and float(l.split()[0]) > 0]
            for x in pos:
                if gpos(x) not in gh:
                    raise SystemExit("STOP: the rule does not give D53's committed ghost centre for the charge at %s in %s." % (x, os.path.relpath(r, S6)))
                n += 1
    if n < 60:
        raise SystemExit("STOP: only %d of D53's ghost centres were checked against the rule." % n)
    return n


def runs(S6, code):
    """[(run, job, cap, kind)] in the order the jobs make them, and the job names. The first job is Hartree-Fock and
    its first run is SMOKE, so that one run with charges and functions goes through ORCA before the others start."""
    p59, _ = mods(code)
    D = p59.designs(S6)
    jobs = []
    hf = [("HF_%s_%s_%s" % (k, d, st), p59.CAP["HF"], k) for st in STATES for d in reversed(D) for k in ("F", "GF")]
    hf = [x for x in hf if x[0] == SMOKE] + [x for x in hf if x[0] != SMOKE]
    for c in p59.chunks(hf, PER_JOB["HF"]):
        jobs.append(c)
    for lv in ("CAM", "WXV"):
        fg = [("%s_%s_%s_%s" % (lv, k, d, st), p59.CAP[k], k) for st in STATES for d in D for k in ("GF", "F")]
        for c in p59.chunks(fg, PER_JOB[lv]):
            jobs.append(c)
    names = ["J%02d" % (k + 1) for k in range(len(jobs))]
    out = [(r, names[i], c, kind) for i, js in enumerate(jobs) for r, c, kind in js]
    assert len(set(r for r, _, _, _ in out)) == len(out) == 120
    assert out[0][0] == SMOKE
    for j in names:
        assert sum(c for _, jj, c, _ in out if jj == j) <= LIMIT_S, j
    return out, names


def ghosts_of(S6, code, design, st):
    p59, p47 = mods(code)
    gpos = ghost_rule(S6, p47)
    return [gpos(tuple(float(v) for v in x)) for q, x in p59.charges_of(S6, design, st) if q > 0]


def build(S6, code):
    """every file of d61/, as {relative path: text}"""
    p59, p47 = mods(code)
    gpos = ghost_rule(S6, p47)
    nchk = check_rule_on_d53(S6, gpos)
    R, J = runs(S6, code)
    sub = {st: rd(os.path.join(S6, "d56/SV_C_ctl_%s/substrate.xyz" % st)) for st in STATES}
    atoms = {st: [l for l in sub[st].split("\n")[2:] if l.strip()] for st in STATES}
    assert all(len(atoms[st]) == 24 for st in atoms)
    heads = dict((lv, p59.d53_head(S6, lv)) for lv in p59.FLEVELS)
    ghost = p59.d53_ghost_block(S6)
    F = {}
    for r, j, cap, kind in R:
        p = r.split("_")
        lv, k, st = p[0], p[1], p[-1]
        d = "_".join(p[2:-1])
        Q = p59.charges_of(S6, d, st)
        pos = [tuple(float(v) for v in x) for q, x in Q if q > 0]
        gh = [gpos(x) for x in pos]
        for x, g in zip(pos, gh):
            dd = dist(x, g)
            if not NEAR_LO < dd < NEAR_HI:
                raise SystemExit("STOP: a ghost centre of %s is %.6f A from its charge." % (r, dd))
        for g in gh:
            for q, x in Q:
                if dist(g, tuple(float(v) for v in x)) < NEAR_LO:
                    raise SystemExit("STOP: a ghost centre of %s lies within %.3f A of a point charge." % (r, NEAR_LO))
        h = "\n".join(heads[lv]) + "\n"
        if k == "F":
            F[r + "/pc.pc"] = rd(os.path.join(S6, "d56", "SV_%s_di_%s" % (d, st), "pc.pc"))
            h += '%pointcharges "pc.pc"\n'
        F[r + "/job.inp"] = h + "* xyz -2 1\n" + "".join("%s\n" % a for a in atoms[st]) + \
            "".join("H: %.6f %.6f %.6f\n" % g + "".join("%s\n" % b for b in ghost) for g in gh) + "*\n"
    F["RUNS.tsv"] = "run\tjob\tcap_s\tkind\tbasis_functions\n" + "".join(
        "%s\t%s\t%d\t%s\t%d\n" % (r, j, c, k, p59.nb_of(S6, r)) for r, j, c, k in R)
    F["GHOSTS.tsv"] = "design\tstate\tcharge_x\tcharge_y\tcharge_z\tghost_x\tghost_y\tghost_z\tdistance_A\n" + "".join(
        "%s\t%s\t%s\t%s\t%s\t%.6f\t%.6f\t%.6f\t%.6f\n" % ((d, st) + tuple(x) + gpos(tuple(float(v) for v in x)) + (dist(tuple(float(v) for v in x), gpos(tuple(float(v) for v in x))),))
        for d in p59.designs(S6) for st in STATES for q, x in p59.charges_of(S6, d, st) if q > 0)
    F["RULE_CHECK.txt"] = "D53's committed ghost centres re-formed by the rule of cs6_d48_prepare.py (GHOST_OFF %.2f A): %d of %d, to six decimals.\n" % (GHOST_OFF, nchk, nchk)
    for j in J:
        F["run_%s.pbs" % j] = pbs(j, [(r, c) for r, jj, c, _ in R if jj == j], len(J))
    F["INPUT_CHECKSUMS.sha256"] = "".join("%s  %s\n" % (hashlib.sha256(F[p].encode()).hexdigest(), p) for p in sorted(F))
    return F


def pbs(j, mine, njobs):
    t = """#!/bin/bash
#PBS -N cm_d61_%s
#PBS -l select=1:ncpus=8:mem=48gb
#PBS -l walltime=24:00:00
#PBS -m ae
#PBS -M 18660916@sun.ac.za
#PBS -j oe
#PBS -o %s/run_%s.pbs.out
# Stage 6 DIAGNOSTIC 61: D59's test of the functions at the positive charges of the ten survivors, made again for the
# runs that carry ghost atoms, each ghost centre 0.01 A beyond its charge (D59 NOTE 1).
# Job %s of %d. It makes %d runs in the order below (name, cap in s); a run that fails does not stop the ones after it.
# No -q: the allocation is automatic and walltime 24:00:00 routes to day. The job script is D59's committed run_J01.pbs
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
    ap.add_argument("--code", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    S6 = os.path.join(a.res, S6REL)
    F = build(S6, a.code)
    dst = os.path.join(a.out, DIR)
    if os.path.exists(dst):
        sys.exit("REFUSED: %s exists already." % dst)
    for p in sorted(F):
        os.makedirs(os.path.dirname(os.path.join(dst, p)) or dst, exist_ok=True)
        with open(os.path.join(dst, p), "w") as f:
            f.write(F[p])
    R, J = runs(S6, a.code)
    print("wrote %d files into %s: %d runs in %d jobs" % (len(F), dst, len(R), len(J)))


if __name__ == "__main__":
    main()
