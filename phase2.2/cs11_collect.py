#!/usr/bin/env python3
"""Stage a finished run's job.out files into its committed input directory.

The inputs were committed before the runs; this copies ONLY the outputs back beside them, after
checking byte for byte that what ran is what was committed. Nothing else is copied: ORCA's scratch
(job.cpcm, job.property.txt, job.bibtex, job.densitiesinfo) is hundreds of kilobytes per run and is
not evidence of anything.

usage: cs11_collect.py <run-dir> <committed-dir> <result-file>
"""
import sys, hashlib
from pathlib import Path


def stop(m):
    print("STOP: %s" % m, file=sys.stderr)
    raise SystemExit(1)


def main(run, dest, result):
    run, dest, result = Path(run), Path(dest), Path(result)
    man = dest / "MANIFEST.sha256"
    if not man.is_file():
        stop("no committed MANIFEST.sha256 in %s" % dest)
    bad, n = [], 0
    for l in man.read_text().splitlines():
        if not l.strip():
            continue
        want, rel = l.split(None, 1)
        f = run / rel.strip()
        if not f.is_file():
            bad.append("%s is missing from the run directory" % rel.strip())
        elif hashlib.sha256(f.read_bytes()).hexdigest() != want:
            bad.append("%s does not match the committed input" % rel.strip())
        n += 1
    if bad:
        stop("the runs are not the committed inputs:\n   " + "\n   ".join(bad))

    lines, got = [], 0
    for d in sorted(p for p in dest.iterdir() if p.is_dir()):
        src = run / d.name / "job.out"
        if not src.is_file():
            stop("missing output: %s" % src)
        if "ORCA TERMINATED NORMALLY" not in src.read_text(errors="replace"):
            stop("%s did not terminate normally" % src)
        (d / "job.out").write_bytes(src.read_bytes())
        lines.append("%s  %s/job.out"
                     % (hashlib.sha256((d / "job.out").read_bytes()).hexdigest(), d.name))
        got += 1
    if not result.is_file():
        stop("result file not found: %s" % result)
    (dest / result.name).write_bytes(result.read_bytes())
    lines.append("%s  %s" % (hashlib.sha256(result.read_bytes()).hexdigest(), result.name))
    (dest / "OUTPUTS.sha256").write_text("\n".join(lines) + "\n")
    print("   %-14s %2d inputs verified, %2d job.out staged, %s copied"
          % (dest.name, n, got, result.name))


if __name__ == "__main__":
    if len(sys.argv) != 4:
        stop(__doc__.strip().splitlines()[-1])
    main(*sys.argv[1:])
