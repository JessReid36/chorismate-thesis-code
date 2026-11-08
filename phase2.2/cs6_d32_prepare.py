#!/usr/bin/env python3
"""DIAGNOSTIC 32: generate the walled (n, q) grid. A PAIRED design.

Every condition reuses DIAGNOSTIC 31's committed sites and charges exactly, with one argon ECP
centre placed at each charge. Nothing about the design is re-derived: the charges and coordinates
are read out of d31's committed pc.pc files, the substrate geometry is copied line for line from
the committed substrate.xyz, and the ECP parameter block is lifted verbatim from a committed
DIAGNOSTIC 23 input. So the only difference from d31 is the presence of the wall.

A walled arm carries its charges on the Ar> lines, not in a pc.pc, which is why the geometry has
to be inline: the ECP centres sit in the same coordinate block as the substrate.
"""
import argparse, hashlib, re, sys
from pathlib import Path

NS = (1, 2, 3, 5, 10, 20)
LADDER = (0.25, 0.50, 0.75, 1.00)
WALL0_NS = (1, 5, 20)
KEYWORD = "! CAM-B3LYP D3BJ def2-SVPD def2/J RIJCOSX TightSCF"


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""): h.update(b)
    return h.hexdigest()


def ecp_block(d23_inp):
    """the argon parameter block, verbatim from a committed DIAGNOSTIC 23 input.
    All ten centres in that file carry the identical block; this refuses to proceed otherwise."""
    t = Path(d23_inp).read_text()
    blocks = re.findall(r"^  Ar>\s+-?\d+\.\d+\s+\S+\s+\S+\s+\S+\s+NewECP\n((?:.*\n)*?)end\n", t, re.M)
    if not blocks: raise SystemExit("no Ar> NewECP block found in %s" % d23_inp)
    uniq = set(blocks)
    if len(uniq) != 1:
        raise SystemExit("the %d centres in %s do not share one parameter block" % (len(blocks), d23_inp))
    b = next(iter(uniq))
    if "N_core 18" not in b: raise SystemExit("the block is not a full argon core (N_core 18)")
    return b.rstrip("\n")


def read_pc(p):
    L = Path(p).read_text().splitlines()
    rows = []
    for ln in L[1:]:
        f = ln.split()
        if len(f) >= 4: rows.append((float(f[0]), f[1], f[2], f[3]))
    if int(L[0]) != len(rows): raise SystemExit("%s: header says %s, has %d" % (p, L[0], len(rows)))
    return rows


def job_inp(geom_lines, rows, block, zero=False):
    L = [KEYWORD, "%maxcore 3000", "%pal nprocs 8 end", "%scf", "  MaxIter 300",
         "  STABPerform true", "end", "* xyz -2 1"]
    L += geom_lines
    for q, x, y, z in rows:
        L.append("  Ar>  %f      %s        %s        %s  NewECP" % (0.0 if zero else q, x, y, z))
        L.append(block)
        L.append("end")
    L += ["*", ""]
    return "\n".join(L)


def main():
    a = argparse.ArgumentParser()
    a.add_argument("--d31", required=True, help="the committed d31 tree, for the sites and charges")
    a.add_argument("--geom", required=True, help="the committed substrate.xyz pair's parent directory")
    a.add_argument("--d23", required=True, help="a committed DIAGNOSTIC 23 argon input, for the ECP block")
    a.add_argument("--out", required=True)
    a = a.parse_args()
    D31, GE, OUT = Path(a.d31), Path(a.geom), Path(a.out)
    block = ecp_block(a.d23)
    print("ECP block taken verbatim from %s" % Path(a.d23).name)
    print("   sha256 of that input: %s" % sha(a.d23))
    print("   %d lines, full argon core" % len(block.splitlines()))
    geom = {}
    for st in ("R", "TS"):
        p = GE / ("ctl_%s" % st) / "substrate.xyz"
        lines = p.read_text().splitlines()
        n = int(lines[0])
        geom[st] = lines[2:2 + n]
        print("   geometry %-2s %d atoms, copied line for line from %s (sha256 %s)"
              % (st, n, p.name, sha(p)[:16]))

    plan = ([("A/n%02d" % n, "A/n%02d" % n, False) for n in NS]
            + [("B/n%02d_q%.2f" % (n, q), "B/n%02d_q%.2f" % (n, q), False)
               for n in NS for q in LADDER]
            + [("wall0/n%02d" % n, "A/n%02d" % n, True) for n in WALL0_NS])
    made, probs, table = 0, [], []
    for name, src, zero in plan:
        for st in ("R", "TS"):
            spc = D31 / ("%s_%s" % (src, st)) / "pc.pc"
            if not spc.exists(): probs.append("%s: no committed d31 pc.pc at %s" % (name, spc)); continue
            rows = read_pc(spc)
            d = OUT / ("%s_%s" % (name, st)); d.mkdir(parents=True, exist_ok=True)
            (d / "job.inp").write_text(job_inp(geom[st], rows, block, zero))
            made += 1
        # the charges and coordinates must be d31's, unchanged
        r_, t_ = read_pc(D31 / ("%s_R" % src) / "pc.pc"), read_pc(D31 / ("%s_TS" % src) / "pc.pc")
        if r_ != t_: probs.append("%s: d31's R and TS charges differ" % name)
        table.append((name, src, zero, len(r_),
                      sum(abs(q) for q, *_ in r_), max((abs(q) for q, *_ in r_), default=0.0)))
    if probs:
        for p in probs[:10]: print("   PROBLEM " + p)
        raise SystemExit("refusing to write: %d problem(s)" % len(probs))

    # read back what was written and check it against d31 and against the geometry
    bad = 0
    for name, src, zero, *_ in table:
        for st in ("R", "TS"):
            t = (OUT / ("%s_%s" % (name, st)) / "job.inp").read_text()
            cen = re.findall(r"^  Ar>\s+(-?\d+\.\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+NewECP$", t, re.M)
            rows = read_pc(D31 / ("%s_%s" % (src, st)) / "pc.pc")
            if len(cen) != len(rows): print("   MISMATCH %s_%s centre count" % (name, st)); bad += 1; continue
            for (cq, cx, cy, cz), (q, x, y, z) in zip(cen, rows):
                want = 0.0 if zero else q
                if abs(float(cq) - want) > 1e-9 or (cx, cy, cz) != (x, y, z):
                    print("   MISMATCH %s_%s a centre does not match d31" % (name, st)); bad += 1; break
            gl = [l for l in t.splitlines() if re.fullmatch(r"  [A-Z][a-z]? .*", l)]
            if gl != geom[st]: print("   MISMATCH %s_%s geometry" % (name, st)); bad += 1
            if t.count("N_core 18") != len(rows):
                print("   MISMATCH %s_%s one ECP block per centre" % (name, st)); bad += 1
    print("\n%d job directories written, %d mismatches" % (made, bad))
    print("\n  condition          paired with        n    Sigma|q|   max|q|")
    for name, src, zero, n, s, m in table:
        print("  %-18s %-18s %-4d %8.4f %8.4f%s" % (name, src, n, 0.0 if zero else s,
                                                     0.0 if zero else m, "   (wall only, no charge)" if zero else ""))
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
