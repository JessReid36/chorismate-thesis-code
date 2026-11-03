#!/usr/bin/env python3
# cs6_relaxA.py -- the relaxation inputs for one readable arrangement.
#
#   python3 cs6_relaxA.py --results <results-repo> --stem n10_bare_A --tree additivity \
#                         --template S6-10 --out <directory>
#
# Stage 6's committed relaxation inputs (relax/S6-5, S6-10, S6-20) carry the
# optimiser's best designs at full field, which DIAGNOSTIC 20 later placed at or
# above the measured instability onset and DIAGNOSTIC 17 placed outside the
# readable window. STAGE6_AMENDMENT3_20261103.txt adds one readable arrangement.
# This script writes its four inputs by taking the committed template of the same
# sphere count and substituting ONLY the ten sphere charges and the ten sphere
# positions, so that every other byte - the protocol, the substrate geometry, the
# force-field parameters, the sphere's Lennard-Jones terms - is the committed one.

import argparse
import hashlib
import sys
from pathlib import Path

JOBS = ("R_sp", "TS_sp", "R_opt", "TS_opt")
STATE = {"R_sp": "R", "R_opt": "R", "TS_sp": "TS", "TS_opt": "TS"}
NSUB = 24                     # the substrate atoms, QM in every job
SCALE = 0.5                   # the field factor between the template and this arrangement


def die(msg):
    sys.stderr.write("cs6_relaxA: %s\n" % msg)
    sys.exit(2)


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def read_pc(p):
    """a pc.pc as ORCA reads it: a count, then q x y z"""
    lines = [l for l in Path(p).read_text().splitlines() if l.strip()]
    n = int(lines[0].split()[0])
    rows = [l.split() for l in lines[1:]]
    if len(rows) != n or any(len(r) != 4 for r in rows):
        die("%s is not a well-formed pc.pc" % p)
    return rows                       # strings, so nothing is re-formatted


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", required=True)
    ap.add_argument("--stem", required=True)        # n10_bare_A
    ap.add_argument("--tree", required=True)        # additivity
    ap.add_argument("--template", required=True)    # S6-10
    ap.add_argument("--name", required=True)        # n10_A
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    res, out = Path(a.results), Path(a.out)
    T6 = res / "phase2.2/cs_stage6"
    tmpl = T6 / "relax" / a.template

    for j in JOBS:
        for f in ("job.inp", "job.pbs", "combined.xyz", "site.ORCAFF.prms"):
            if not (tmpl / j / f).is_file():
                die("the template %s/%s has no %s" % (a.template, j, f))

    # the design, from its own committed frozen inputs, one pc.pc per state
    pcs = {}
    for st in ("R", "TS"):
        p = T6 / a.tree / "def2-SVPD" / ("%s_%s" % (a.stem, st)) / "pc.pc"
        if not p.is_file():
            die("no committed pc.pc for %s at %s" % (a.stem, st))
        pcs[st] = read_pc(p)
    if [r[0] for r in pcs["R"]] != [r[0] for r in pcs["TS"]]:
        die("%s carries different charges at R and TS" % a.stem)
    nsph = len(pcs["R"])

    written = []
    for j in JOBS:
        st = STATE[j]
        src = tmpl / j
        dst = out / j
        dst.mkdir(parents=True, exist_ok=True)

        # combined.xyz: byte for byte the committed one. The sites are the same sites;
        # only the charges differ, and charges do not live in the xyz.
        raw = (src / "combined.xyz").read_text()
        lines = raw.splitlines()
        if int(lines[0].split()[0]) != NSUB + nsph:
            die("%s/%s has %s atoms; this design needs %d spheres on %d substrate atoms"
                % (a.template, j, lines[0].split()[0], nsph, NSUB))
        for k, r in enumerate(pcs[st]):
            f = lines[2 + NSUB + k].split()
            if any(abs(float(f[1 + c]) - float(r[1 + c])) > 1e-9 for c in range(3)):
                die("sphere %d of %s/%s is not where %s's own pc.pc puts it"
                    % (k + 1, a.template, j, a.stem))
        (dst / "combined.xyz").write_text(raw)

        # the force field: the committed file with only the sphere charges substituted
        prms = (src / "site.ORCAFF.prms").read_text().splitlines()
        i = prms.index("$atoms")
        n = int(prms[i + 1].split()[0])
        rows = prms[i + 2:i + 2 + n]
        if n != NSUB + nsph:
            die("%s/%s's force field has %d atoms, not %d" % (a.template, j, n, NSUB + nsph))
        new = []
        for k, row in enumerate(rows):
            if k < NSUB:
                new.append(row)
                continue
            f = row.split()
            if f[1] != "N":
                die("sphere %d of the template is %r, not an Amber N sphere" % (k - NSUB + 1, f[1]))
            q = float(f[2]) * SCALE
            want = float(pcs[st][k - NSUB][0])
            if abs(q - want) > 1e-12:
                die("sphere %d: the template's %.6f at %g is %.6f, and %s's pc.pc says %.6f"
                    % (k - NSUB + 1, float(f[2]), SCALE, q, a.stem, want))
            new.append("%6s   %-3s %12.6f %12.6f %12.6f %12.6f %12.6f"
                       % (f[0], f[1], q, float(f[3]), float(f[4]), float(f[5]), float(f[6])))
        (dst / "site.ORCAFF.prms").write_text(
            "\n".join(prms[:i + 2] + new + prms[i + 2 + n:]) + "\n")

        # the protocol, byte for byte
        (dst / "job.inp").write_text((src / "job.inp").read_text())

        # the submission script, with only the job name and the two paths changed
        pbs = (src / "job.pbs").read_text()
        pbs = pbs.replace("cm_s6_%s_%s" % (a.template.replace("S6-", ""), j),
                          "cm_s6_%s_%s" % (a.name.replace("-", ""), j))
        pbs = pbs.replace("relax/%s/%s" % (a.template, j), "relax/%s/%s" % (a.name, j))
        (dst / "job.pbs").write_text(pbs)

        for f in ("job.inp", "job.pbs", "combined.xyz", "site.ORCAFF.prms"):
            written.append(("relax/%s/%s/%s" % (a.name, j, f), sha(dst / f)))

    sys.stdout.write("   %d spheres, 4 jobs, %d files\n" % (nsph, len(written)))
    (out / "SHA256").write_text("\n".join("%-44s  %s" % (p, s) for p, s in written) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
