#!/usr/bin/env python3
"""check_input_echo.py FOLDER - for every run folder in FOLDER that has job.inp and job.out (or job.out.xz), check that
the input ORCA echoed at the top of job.out is the committed job.inp, line for line. Proves each output came from the
committed input (for the zero-charge runs: every sphere charge 0.000000). Prints the count; exit 1 on any mismatch."""
import glob, lzma, os, re, sys
D = sys.argv[1]
def echo(t):
    i = t.find("INPUT FILE"); j = t.find("****END OF INPUT****", i)
    if i < 0 or j < 0: return None
    e = [m.group(1).rstrip() for m in re.finditer(r"^\|\s*\d+> ?(.*)$", t[i:j], re.M)]
    while e and not e[-1]: e.pop()
    return e
n = bad = 0
for d in sorted(glob.glob(os.path.join(D, "*", "job.inp"))):
    r = os.path.dirname(d); p = os.path.join(r, "job.out")
    if os.path.exists(p): t = open(p, errors="replace").read()
    elif os.path.exists(p + ".xz"): t = lzma.open(p + ".xz", "rt", errors="replace").read()
    else: continue
    inp = [l.rstrip() for l in open(d).read().splitlines()]
    while inp and not inp[-1]: inp.pop()
    n += 1
    if echo(t) != inp:
        bad += 1; print("input echoed in the output differs from the committed job.inp: " + os.path.basename(r))
print("%d outputs checked against their committed inputs, %d mismatches" % (n, bad))
sys.exit(1 if bad else 0)
