#!/usr/bin/env python3
"""DIAGNOSTIC 33 stage B: the inputs for the designs stage A selected.

Bare arms take the d31 route: %pointcharges with a pc.pc. Walled arms take the d32 route: inline
geometry, the committed DIAGNOSTIC 23 argon block at every centre, and no pc.pc, because ORCA will
not take both. The argon block is lifted verbatim from a committed d23 input by the same check
cs6_d32_prepare.py uses, and this refuses to proceed if the centres there do not share one block.
"""
import re, sys, hashlib
from pathlib import Path

HEAD = ["! CAM-B3LYP D3BJ def2-SVPD def2/J RIJCOSX TightSCF",
        "%maxcore 3000", "%pal nprocs 8 end", "%scf", "  MaxIter 300", "  STABPerform true", "end"]

def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""): h.update(b)
    return h.hexdigest()

def ecp_block(d23_inp):
    """the argon parameter block, verbatim from a committed DIAGNOSTIC 23 input"""
    t = Path(d23_inp).read_text()
    blocks = re.findall(r"^  Ar>\s+-?\d+\.\d+\s+\S+\s+\S+\s+\S+\s+NewECP\n((?:.*\n)*?)end\n", t, re.M)
    if not blocks: raise SystemExit("no Ar> NewECP block found in %s" % d23_inp)
    if len(set(blocks)) != 1:
        raise SystemExit("the %d centres in %s do not share one parameter block" % (len(blocks), d23_inp))
    b = next(iter(set(blocks)))
    if "N_core 18" not in b: raise SystemExit("the block is not a full argon core (N_core 18)")
    return b.rstrip("\n")

def sites_xyz(p, code):
    """site coordinates from the COMMITTED optimiser's own reader, not a parser written here"""
    import importlib.util
    spec = importlib.util.spec_from_file_location("opt", str(Path(code) / "phase2.2" / "cs_opt_milp.py"))
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
    xyz, _layer, _ranges = m.read_sites(str(p))
    return {i: (float(x), float(y), float(z)) for i, (x, y, z) in enumerate(xyz)}

def main(designs, sites_p, geom_dir, d23_inp, code, out, hpc):
    out = Path(out); out.mkdir(parents=True, exist_ok=True)
    S = sites_xyz(sites_p, code)
    ECP = ecp_block(d23_inp)
    rows = [l.split("\t") for l in Path(designs).read_text().splitlines()
            if l.strip() and not l.startswith("#") and not l.startswith("name")]
    if not rows: raise SystemExit("no designs in %s" % designs)
    geom = {st: Path(geom_dir) / ("ctl_%s" % st) / "substrate.xyz" for st in ("R", "TS")}
    for st, g in geom.items():
        if not g.is_file(): raise SystemExit("no committed geometry at %s" % g)
    man, made = [], []

    def write_run(name, st, pairs, walled):
        d = out / ("%s_%s" % (name, st)); d.mkdir(parents=True, exist_ok=True)
        lines = list(HEAD)
        if walled:
            xyz = geom[st].read_text().splitlines()
            n = int(xyz[0].split()[0])
            atoms = [l for l in xyz[2:2 + n] if l.strip()]
            if len(atoms) != n: raise SystemExit("%s: expected %d atoms" % (geom[st], n))
            lines.append("* xyz -2 1")
            lines += ["  " + " ".join(l.split()) for l in atoms]
            for q, (x, y, z) in pairs:
                lines.append("  Ar>  %9.6f      %.10f        %.10f        %.10f  NewECP" % (q, x, y, z))
                lines.append(ECP)
                lines.append("end")
            lines.append("*")
        else:
            lines.append('%pointcharges "pc.pc"')
            lines.append("* xyzfile -2 1 substrate.xyz")
            pc = ["%d" % len(pairs)]
            for q, (x, y, z) in pairs:
                pc.append(" %9.6f      %.10f      %.10f      %.10f" % (q, x, y, z))
            (d / "pc.pc").write_text("\n".join(pc) + "\n")
            man.append((sha(d / "pc.pc"), "%s_%s/pc.pc" % (name, st)))
        (d / "job.inp").write_text("\n".join(lines) + "\n")
        man.append((sha(d / "job.inp"), "%s_%s/job.inp" % (name, st)))
        if not walled:
            (d / "substrate.xyz").write_bytes(geom[st].read_bytes())
            man.append((sha(d / "substrate.xyz"), "%s_%s/substrate.xyz" % (name, st)))
        made.append("%s_%s" % (name, st))

    wall_ns = {}
    for r in rows:
        name, n, q, arm = r[0], int(r[1]), float(r[2]), r[3]
        sites = [int(x) for x in r[7].split(",")]
        chg = [float(x) for x in r[8].split(",")]
        if len(sites) != n or len(chg) != n:
            raise SystemExit("%s: %d sites and %d charges for n = %d" % (name, len(sites), len(chg), n))
        if any(abs(abs(c) - q) > 1e-9 for c in chg):
            raise SystemExit("%s: a charge is not |q| = %.2f exactly" % (name, q))
        pairs = [(c, S[s]) for c, s in zip(chg, sites)]
        for st in ("R", "TS"): write_run(name, st, pairs, arm == "wall")
        if arm == "wall": wall_ns[n] = [S[s] for s in sites]
    for n, pts in sorted(wall_ns.items()):
        for st in ("R", "TS"):
            write_run("wall0_n%02d" % n, st, [(0.0, p) for p in pts], True)

    pbs = ["#!/bin/bash", "#PBS -N cm_d33b", "#PBS -l select=1:ncpus=8:mem=24gb",
           "#PBS -l walltime=12:00:00", "#PBS -m ae", "#PBS -M 18660916@sun.ac.za",
           "#PBS -j oe", "#PBS -o %s/run.pbs.out" % hpc,
           "# Stage 6 DIAGNOSTIC 33 stage B: %d runs, CAM-B3LYP/def2-SVPD, STABPerform on every one." % len(made),
           "# Ordered ASCENDING in Sigma|q|, so the designs inside D30's confirmed window land first",
           "# and an overrun loses only the two that are outside it. No -q: the allocation is",
           "# automatic and walltime 12:00:00 routes through workq to day. Restartable.",
           'export PATH="/apps/openmpi/4.1.1/bin:$PATH"',
           'export LD_LIBRARY_PATH="/home/apps2/ORCA/6.0.1/lib:/apps/openmpi/4.1.1/lib:/apps/mambaforge/envs/medaka/lib:${LD_LIBRARY_PATH:-}"',
           "export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1",
           "ORCA=/home/apps2/ORCA/6.0.1", "W=%s" % hpc, 'cd "$W" || exit 1',
           'if ! mkdir .running 2>/dev/null; then echo "FAIL: another job holds .running in $W"; exit 1; fi',
           "trap 'rmdir \"$W/.running\" 2>/dev/null' EXIT",
           'echo "host=$(hostname) start=$(date)"',
           "run() {", '  d="$W/$1"',
           '  if [ -s "$d/job.out" ] && grep -q "ORCA TERMINATED NORMALLY" "$d/job.out"; then echo "$1 done before"',
           '  else ( cd "$d" && "$ORCA/orca" job.inp > job.out 2>&1 </dev/null ); echo "$1 orca_exit=$? normal=$(grep -c \'ORCA TERMINATED NORMALLY\' "$d/job.out") $(date)"; fi',
           "}"]
    order = {r[0]: float(r[4]) for r in rows}
    seq = sorted(made, key=lambda nm: (order.get(nm.rsplit("_", 1)[0], -1.0), nm))
    for nm in seq: pbs.append("run %s" % nm)
    pbs += ['echo "end=$(date)"',
            'echo "finished: $(grep -l \'ORCA TERMINATED NORMALLY\' */job.out 2>/dev/null | wc -l) of %d"' % len(made)]
    (out / "run.pbs").write_text("\n".join(pbs) + "\n")
    man.append((sha(out / "run.pbs"), "run.pbs"))
    Path(out, "INPUT_CHECKSUMS.sha256").write_text(
        "".join("%s  %s\n" % (h, p) for h, p in sorted(man, key=lambda x: x[1])))
    print("%d runs prepared under %s" % (len(made), out))
    print("  %d designs, %d wall-only controls at n = %s"
          % (len(rows), len(wall_ns), ", ".join(str(n) for n in sorted(wall_ns))))
    print("  %d files checksummed" % len(man))

if __name__ == "__main__":
    if len(sys.argv) != 8:
        raise SystemExit("usage: cs6_d33b_prepare.py <designs.tsv> <sites.tsv> <geom dir> <d23 inp> <code dir> <out> <hpc>")
    main(*sys.argv[1:])
