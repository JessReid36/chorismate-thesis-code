#!/usr/bin/env python3
"""DIAGNOSTIC 33 stage A: the A matrix at CAM-B3LYP/def2-SVPD.

Every input is LIFTED from the committed Stage 5 A matrix and has exactly one line changed - the
method line. The geometries, the point-charge files, the evaluation points and the run names are
copied byte for byte, so a_i^new and the committed a_i differ by the level of theory and by
nothing else. This refuses to proceed if the committed method line is not what it expects.
"""
import sys, shutil, hashlib
from pathlib import Path

OLD = "! B3LYP D3BJ def2-SVP def2/J RIJCOSX TightSCF"
NEW = "! CAM-B3LYP D3BJ def2-SVPD def2/J RIJCOSX TightSCF"
COPY = ("points.xyz", "points_inv.xyz", "row_order.tsv", "runs.tsv")

def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""): h.update(b)
    return h.hexdigest()

def main(src, dst, hpc):
    src, dst = Path(src), Path(dst)
    if not (src / "runs.tsv").is_file():
        raise SystemExit("no runs.tsv under %s; that is the committed A matrix tree" % src)
    runs = [l.split("\t")[0] for l in (src / "runs.tsv").read_text().splitlines()
            if l.strip() and not l.startswith("#")]
    if len(runs) != 16:
        raise SystemExit("expected 16 committed A-matrix runs, runs.tsv lists %d" % len(runs))
    dst.mkdir(parents=True, exist_ok=True)
    man = []
    for r in runs:
        s, d = src / r, dst / r
        inp = s / "job.inp"
        if not inp.is_file(): raise SystemExit("%s has no job.inp" % s)
        lines = inp.read_text().splitlines(keepends=True)
        hits = [i for i, l in enumerate(lines) if l.rstrip("\n") == OLD]
        if len(hits) != 1:
            raise SystemExit("%s: expected exactly one committed method line, found %d.\n"
                             "       Refusing to guess which line to change." % (inp, len(hits)))
        d.mkdir(parents=True, exist_ok=True)
        lines[hits[0]] = NEW + "\n"
        (d / "job.inp").write_text("".join(lines))
        for extra in ("substrate.xyz", "pc.pc"):
            if (s / extra).is_file():
                shutil.copyfile(s / extra, d / extra)
                man.append((sha(d / extra), "%s/%s" % (r, extra)))
        man.append((sha(d / "job.inp"), "%s/job.inp" % r))
    for f in COPY:
        if not (src / f).is_file(): raise SystemExit("the committed tree has no %s" % f)
        shutil.copyfile(src / f, dst / f)
        man.append((sha(dst / f), f))

    npts = int((dst / "points.xyz").read_text().split("\n", 1)[0].strip())
    sp  = [r for r in runs if r.startswith(("sp_", "inv_"))]
    fd  = [r for r in runs if r.startswith("fd_")]
    pbs = ["#!/bin/bash",
           "#PBS -N cm_d33a",
           "#PBS -l select=1:ncpus=8:mem=24gb",
           "#PBS -l walltime=06:00:00",
           "#PBS -m ae",
           "#PBS -M 18660916@sun.ac.za",
           "#PBS -j oe",
           "#PBS -o %s/run.pbs.out" % hpc,
           "# Stage 6 DIAGNOSTIC 33 stage A: the A matrix at CAM-B3LYP/def2-SVPD. %d runs, lifted" % len(runs),
           "# from the committed Stage 5 A matrix with only the method line changed. No -q: the",
           "# allocation is automatic and walltime 06:00:00 routes through workq to day.",
           "# The charge-free control at CAM-B3LYP/def2-SVPD takes about 5 minutes, so %d runs fit" % len(runs),
           "# inside this window with room to spare. Restartable: resubmit and it skips what finished.",
           'export PATH="/apps/openmpi/4.1.1/bin:$PATH"',
           'export LD_LIBRARY_PATH="/home/apps2/ORCA/6.0.1/lib:/apps/openmpi/4.1.1/lib:/apps/mambaforge/envs/medaka/lib:${LD_LIBRARY_PATH:-}"',
           "export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1",
           "ORCA=/home/apps2/ORCA/6.0.1",
           'W=%s' % hpc,
           'cd "$W" || exit 1',
           'if ! mkdir .running 2>/dev/null; then echo "FAIL: another job holds .running in $W"; exit 1; fi',
           "trap 'rmdir \"$W/.running\" 2>/dev/null' EXIT",
           'echo "host=$(hostname) start=$(date)"',
           "run() {",
           '  d="$W/$1"',
           '  if [ -s "$d/job.out" ] && grep -q "ORCA TERMINATED NORMALLY" "$d/job.out"; then echo "$1 done before"',
           '  else ( cd "$d" && "$ORCA/orca" job.inp > job.out 2>&1 </dev/null ); echo "$1 orca_exit=$? normal=$(grep -c \'ORCA TERMINATED NORMALLY\' "$d/job.out") $(date)"; fi',
           "}",
           "vpot() {",
           '  d="$W/$1"; pts="$W/$2"',
           '  if [ -s "$d/vpot.out" ]; then echo "$1 vpot done before"; return; fi',
           '  if ! grep -q "ORCA TERMINATED NORMALLY" "$d/job.out" 2>/dev/null; then echo "$1 vpot skipped, run not normal"; return; fi',
           '  ( cd "$d" && "$ORCA/orca_vpot" job.gbw job.scfp "$pts" vpot.out job > vpot.log 2>&1 )',
           '  echo "$1 vpot_exit=$? lines=$(wc -l < "$d/vpot.out" 2>/dev/null || echo 0) expected=%d"' % (npts + 1),
           "}",
           "# the three bare single points and the inverse-transform check first: they carry the potentials",
           ]
    for r in sp: pbs.append("run %s" % r)
    pbs.append("# then the finite-difference pairs that validate a_i against a_fd")
    for r in fd: pbs.append("run %s" % r)
    pbs.append("# orca_vpot at every site, on the runs that carry potentials")
    for l in (src / "runs.tsv").read_text().splitlines():
        if not l.strip() or l.startswith("#"): continue
        c = l.split("\t")
        if len(c) >= 7 and c[6].strip() and c[6].strip() != "-":
            pbs.append("vpot %s %s" % (c[0], c[6].strip()))
    pbs += ['echo "end=$(date)"',
            'echo "finished: $(grep -lc \'ORCA TERMINATED NORMALLY\' */job.out 2>/dev/null | wc -l) of %d"' % len(runs)]
    (dst / "run.pbs").write_text("\n".join(pbs) + "\n")
    man.append((sha(dst / "run.pbs"), "run.pbs"))
    (dst / "INPUT_CHECKSUMS.sha256").write_text(
        "".join("%s  %s\n" % (h, p) for h, p in sorted(man, key=lambda x: x[1])))
    print("%d runs prepared under %s" % (len(runs), dst))
    print("  method line: %s" % NEW)
    print("  %d evaluation points, %d potential-carrying runs, %d finite-difference runs"
          % (npts, len(sp), len(fd)))
    print("  %d files checksummed in INPUT_CHECKSUMS.sha256" % len(man))

if __name__ == "__main__":
    if len(sys.argv) != 4:
        raise SystemExit("usage: cs6_d33a_prepare.py <committed a_matrix dir> <out dir> <hpc path>")
    main(*sys.argv[1:])
