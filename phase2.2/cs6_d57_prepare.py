#!/usr/bin/env python3
"""cs6_d57_prepare.py - form the inputs of STAGE 6, DIAGNOSTIC 57 from committed files.

AMENDMENT 12, A12.4. This program writes the inputs; it computes nothing and reads no output.

  W runs: at each number of charges of D56, the best-ranked design that D56's committed result calls
          NOT A SURVIVOR (A12.4 item 1, applied here by program), with each positive charge a coreless
          centre carrying the charge behind the helium set (WHe) or the neon set (WNe), each negative
          charge a bare point charge (item 5); dianion and monoanion, R and TS, both levels (item 2).
  A runs: of those designs, the one with the most charges (ranked better on a tie) with EVERY charge
          walled, both sets (item 6).
  Z runs: the controls of the wall alone: the same centres with no charge and no point charge, dianion,
          R and TS, both levels, for each set and for each design of W and A.
  I runs: the measurement of the wall against real ions at one site (item 4): the substrate alone;
          +1 e as a bare point charge; +1 e behind the helium set, the neon set and the published capping
          block (A12.3 item 5); each of the three with no charge; Na+ and K+ as atoms; their functions
          alone (ghosts); R and TS; both levels. Dianion throughout.

Every geometry and every position of a charge is copied from a committed file: the substrate and the
point charges of D56's committed inputs, the site of item 4 from sites_f42436.tsv, the helium block
from d23/he_d16_S6-10_wall0_R/job.inp and the neon block from additivity/def2-SVP/n10_wall0_R/job.inp.
The method line, the SCF block and its stability analysis are D56's.

usage:
   python3 cs6_d57_prepare.py --res <results repository> --out <folder to write d57 into>
The reader cs6_d57.py calls build() and stops unless the committed d57/ is what build() gives.
"""
import argparse, csv, hashlib, os, sys

S6REL = "phase2.2/cs_stage6"
DIR = "d57"
HPC = "/home/18660916/system_development/phase2.2/cs_stage6/" + DIR
SITE_I = 1690                       # the site of item 4: DIAGNOSTIC 3's site
LEVELS = (("SV", "def2-SVPD", 1800), ("TZ", "def2-TZVPD", 5000))
SETS = ("He", "Ne")
NS = (3, 4, 5, 6, 8, 10, 12, 15, 20)
SHI = ["  N_core 0", "  lmax f", "  s 1", "    1  1.732000000  14.676000000  2", "  p 1", "    1  1.115000000  5.175700000  2",
       "  d 1", "    1  1.203000000  -1.816000000  2", "  f 1", "    1  1.000000000  0.000000000  2", "end"]
IKINDS = ("ctl", "bare", "He", "Ne", "Shi", "He0", "Ne0", "Shi0", "Na", "K", "gNa", "gK")


def rd(p):
    with open(p) as f:
        return f.read()


def walled_designs(S6):
    """A12.4 item 1, by program: at each n, the NOT A SURVIVOR with the smallest rank_at_n"""
    rows = list(csv.DictReader(open(os.path.join(S6, "d56_survivors.tsv")), delimiter="\t"))
    out = []
    for n in NS:
        c = sorted([r for r in rows if int(r["n"]) == n and r["the_design"] != "SURVIVOR"], key=lambda r: int(r["rank_at_n"]))
        if c:
            out.append((c[0]["design"], n, int(c[0]["rank_at_n"])))
    return out


def both_ways(S6):
    """A12.4 item 6: of the walled designs, the one with the most charges, ranked better on a tie"""
    w = walled_designs(S6)
    return sorted(w, key=lambda x: (-x[1], x[2]))[0][0] if w else None


def head(basis):
    return "! CAM-B3LYP D3BJ %s def2/J RIJCOSX TightSCF\n%%maxcore 3000\n%%pal nprocs 8 end\n%%scf\n  MaxIter 300\n  STABPerform true\nend\n" % basis


def block(path, tag):
    L = rd(path).split("\n")
    i = [k for k, l in enumerate(L) if l.split()[:1] == [tag]][0]
    j = i + 1
    while L[j].strip() != "end":
        j += 1
    return L[i + 1:j + 1]


def charges_of(S6, design):
    """[(q, 'x y z' as written in D56's pc.pc)] of the design's committed point charges"""
    L = [l.split() for l in rd(os.path.join(S6, "d56", "SV_%s_di_R" % design, "pc.pc")).split("\n")[1:] if l.strip()]
    return [(float(c[0]), (c[1], c[2], c[3])) for c in L]


def runs(S6):
    """[(run, job, cap)] in the order the jobs make them"""
    W = walled_designs(S6)
    A = both_ways(S6)
    jobs = []
    for k, (d, n, rk) in enumerate(W):
        jobs.append(("W%d" % (k + 1), ["%s_W%s_%s_%s_%s" % (lv, s, d, kd, st) for s in SETS for lv, _, _ in LEVELS for kd in ("di", "mo") for st in ("R", "TS")]))
    if A:
        jobs.append(("W%d" % (len(W) + 1), ["%s_A%s_%s_%s_%s" % (lv, s, A, kd, st) for s in SETS for lv, _, _ in LEVELS for kd in ("di", "mo") for st in ("R", "TS")]))
    variants = [("Z", d) for d, _, _ in W] + ([("Y", A)] if A else [])
    for lv, _, _ in LEVELS:
        jobs.append(("W%d" % (len(jobs) + 1), ["%s_%s%s_%s_di_%s" % (lv, t, s, d, st) for t, d in variants for s in SETS for st in ("R", "TS")]))
    ion = lambda lv: ["%s_I_%s_%s" % (lv, k, st) for k in IKINDS for st in ("R", "TS")]
    jobs.append(("I1", ion("SV")))
    jobs.append(("I2", ion("TZ")[:12]))
    jobs.append(("I3", ion("TZ")[12:]))
    cap = dict((lv, c) for lv, _, c in LEVELS)
    out = [(r, j, cap[r[:2]]) for j, rs in jobs for r in rs]
    assert len(set(r for r, _, _ in out)) == len(out)
    return out, [j for j, _ in jobs]


def build(S6):
    """every file of d57/, as {relative path: text}"""
    R, J = runs(S6)
    sub = {st: rd(os.path.join(S6, "d56/SV_C_ctl_%s/substrate.xyz" % st)) for st in ("R", "TS")}
    atoms = {st: [l for l in sub[st].split("\n")[2:] if l.strip()] for st in ("R", "TS")}
    assert all(len(atoms[st]) == 24 for st in atoms)
    blk = {"He": ("He>", block(os.path.join(S6, "d23/he_d16_S6-10_wall0_R/job.inp"), "He>")),
           "Ne": ("Ne>", block(os.path.join(S6, "additivity/def2-SVP/n10_wall0_R/job.inp"), "Ne>")),
           "Shi": ("Na>", SHI)}
    site = None
    for l in rd(os.path.join(S6, "sites_f42436.tsv")).split("\n"):
        c = l.split("\t")
        if c[0] == str(SITE_I):
            site = tuple("%.10f" % float(v) for v in c[2:5])
    assert site is not None
    basis = dict((lv, b) for lv, b, _ in LEVELS)
    F = {}

    def centres(st, items):
        return "".join("%s\n" % a for a in atoms[st]) + "".join(
            "  %s  %.6f  %s  %s  %s  NewECP\n" % (blk[s][0], q, x[0], x[1], x[2]) + "".join("%s\n" % b for b in blk[s][1]) for s, q, x in items)

    for r, j, cap in R:
        lv, t = r[:2], r.split("_")[1]
        st = r.rsplit("_", 1)[1]
        h = head(basis[lv])
        if t[0] in "WAZY":
            s = t[1:]
            d = r[len(lv) + 1 + len(t) + 1:].rsplit("_", 2)[0]
            kd = r.rsplit("_", 2)[1]
            Q = charges_of(S6, d)
            if t[0] == "W":
                walls, pcs = [(s, q, x) for q, x in Q if q > 0], [(q, x) for q, x in Q if q < 0]
            elif t[0] == "A":
                walls, pcs = [(s, q, x) for q, x in Q], []
            elif t[0] == "Z":
                walls, pcs = [(s, 0.0, x) for q, x in Q if q > 0], []
            else:
                walls, pcs = [(s, 0.0, x) for q, x in Q], []
            chg, mult = (-2, 1) if kd == "di" else (-1, 2)
            txt = h
            if pcs:
                txt += '%pointcharges "pc.pc"\n'
                F[r + "/pc.pc"] = "%d\n" % len(pcs) + "".join(" %9.6f  %s  %s  %s\n" % (q, x[0], x[1], x[2]) for q, x in pcs)
            F[r + "/job.inp"] = txt + "* xyz %d %d\n" % (chg, mult) + centres(st, walls) + "*\n"
        else:
            k = r.split("_")[2]
            if k == "ctl":
                F[r + "/job.inp"] = h + "* xyzfile -2 1 substrate.xyz\n"
                F[r + "/substrate.xyz"] = sub[st]
            elif k == "bare":
                F[r + "/job.inp"] = h + '%pointcharges "pc.pc"\n* xyzfile -2 1 substrate.xyz\n'
                F[r + "/pc.pc"] = "1\n %9.6f  %s  %s  %s\n" % (1.0, site[0], site[1], site[2])
                F[r + "/substrate.xyz"] = sub[st]
            elif k in ("He", "Ne", "Shi", "He0", "Ne0", "Shi0"):
                F[r + "/job.inp"] = h + "* xyz -2 1\n" + centres(st, [(k.rstrip("0"), 0.0 if k.endswith("0") else 1.0, site)]) + "*\n"
            else:
                sym, chg = {"Na": ("Na", -1), "K": ("K", -1), "gNa": ("Na:", -2), "gK": ("K:", -2)}[k]
                F[r + "/job.inp"] = h + "* xyzfile %d 1 system.xyz\n" % chg
                F[r + "/system.xyz"] = "25\nsubstrate of frame 42436 at %s (d56/SV_C_ctl_%s), then %s at site %d\n" % (st, st, sym, SITE_I) + \
                    "".join("%s\n" % a for a in atoms[st]) + "  %-3s %s  %s  %s\n" % (sym, site[0], site[1], site[2])
    F["RUNS.tsv"] = "run\tjob\tcap_s\n" + "".join("%s\t%s\t%d\n" % x for x in R)
    for j in J:
        F["run_%s.pbs" % j] = pbs(j, [(r, c) for r, jj, c in R if jj == j], len(J))
    F["INPUT_CHECKSUMS.sha256"] = "".join("%s  %s\n" % (hashlib.sha256(F[p].encode()).hexdigest(), p) for p in sorted(F))
    return F


def pbs(j, mine, njobs):
    t = """#!/bin/bash
#PBS -N cm_d57_%s
#PBS -l select=1:ncpus=8:mem=48gb
#PBS -l walltime=24:00:00
#PBS -m ae
#PBS -M 18660916@sun.ac.za
#PBS -j oe
#PBS -o %s/run_%s.pbs.out
# Stage 6 DIAGNOSTIC 57: D56's best-ranked failed designs behind the wall, and the wall against real ions (AMENDMENT 12).
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
    print("wrote %d files into %s: %d runs in %d jobs; walled designs %s; every charge walled: %s" % (
        len(F), dst, len(R), len(J), ", ".join(d for d, _, _ in walled_designs(S6)), both_ways(S6)))


if __name__ == "__main__":
    main()
