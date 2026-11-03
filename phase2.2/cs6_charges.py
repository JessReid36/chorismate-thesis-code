#!/usr/bin/env python3
# cs6_charges.py -- write the Stage 6 per-site charge vectors to a data file.
#
#   python3 cs6_charges.py --results <results-repo> --out <directory>
#
# The nine whole-design conditions of Stage 6 each carry a set of point charges.
# Until now those vectors existed only inside the committed ORCA inputs (each
# bare arm's pc.pc) and, in prose, inside relax/ARRANGEMENTS.txt. A max|q| column
# populated from memory rather than from pc.pc is what
# RESULTS_charge_representation_domain_DRAFT_ERRATUM1_20261102.md corrected. This
# script makes the vectors a file, so no figure derived from them has to be
# remembered again.
#
# The source of truth is each bare arm's committed pc.pc: charge, then x, y, z,
# in the order ORCA read them. Site identity, layer, host atom and clearance come
# from sites_f42436.tsv by exact coordinate match. Nothing is rounded on the way
# in; the TSV carries the bytes the inputs carry.

import argparse
import hashlib
import math
import sys
from pathlib import Path

SITES = "phase2.2/cs_stage6/sites_f42436.tsv"
TREES = ["additivity/def2-SVPD", "negatives/def2-SVPD", "negatives_armA/def2-SVPD",
         "d16/def2-SVPD", "d17/def2-SVPD"]
# label, bare-arm stem, the diagnostic the condition belongs to, readable under D6.3
NINE = [("S6-2 at 1.3000 e",   "S6-2_bare_F1",   "D16", "no"),
        ("S6-2 at 1.0777 e",   "S6-2_bare_M",    "D17", "yes"),
        ("S6-3 at 1.3000 e",   "S6-3_bare_F1",   "D16", "yes"),
        ("S6-5 at 1.3000 e",   "S6-5_bare_F1",   "D16", "yes"),
        ("S6-5 at 1.6166 e",   "S6-5_bare_M",    "D17", "no"),
        ("S6-10 at 1.3000 e",  "S6-10_bare_F1",  "D16", "yes"),
        ("S6-10 at 2.1900 e",  "n10_bare_A",     "D15", "yes"),
        ("S6-20 at 1.3000 e",  "S6-20_bare_F1",  "D16", "yes"),
        ("S6-20 at 2.1900 e",  "S6-20_bare_F2",  "D16", "yes")]
DETAIL = "CHARGE_VECTORS_20261103.tsv"
SUMMARY = "CHARGE_VECTORS_SUMMARY_20261103.tsv"
TOL = 1e-6          # a pc.pc coordinate and its site must agree to this, in angstrom


def die(msg):
    sys.stderr.write("cs6_charges: %s\n" % msg)
    sys.exit(2)


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def read_pc(p):
    """A pc.pc as ORCA reads it: a count, then 'q x y z' per line."""
    lines = [l for l in Path(p).read_text().splitlines() if l.strip()]
    n = int(lines[0].split()[0])
    rows = []
    for l in lines[1:]:
        f = l.split()
        if len(f) != 4:
            die("%s has a point with %d fields, not 4" % (p, len(f)))
        rows.append((l.split()[0], float(f[0]), float(f[1]), float(f[2]), float(f[3])))
    if len(rows) != n:
        die("%s says %d points and carries %d" % (p, n, len(rows)))
    return rows


def read_sites(p):
    out = []
    for l in Path(p).read_text().splitlines():
        if not l.strip() or l.startswith("#"):
            continue
        f = l.split("\t")
        if f[0] == "site":
            hdr = f
            continue
        out.append({hdr[i]: f[i] for i in range(len(f))})
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    res, out = Path(a.results), Path(a.out)
    T6 = res / "phase2.2/cs_stage6"
    sites = read_sites(res / SITES)
    ssha = sha(res / SITES)

    det = ["\t".join(("condition", "diagnostic", "stem", "tree", "point", "site", "layer",
                      "on_atom", "x", "y", "z", "q", "clearance_A", "q_min_neg", "q_min_pos"))]
    summ = ["\t".join(("condition", "diagnostic", "stem", "tree", "n", "sum_abs_q", "max_abs_q",
                       "min_q", "max_q", "n_positive", "n_negative", "readable_D6_3",
                       "source", "source_sha256"))]
    for label, stem, diag, readable in NINE:
        tree = next((t for t in TREES if (T6 / t / (stem + "_R") / "pc.pc").is_file()), None)
        if tree is None:
            die("no committed tree carries %s" % stem)
        pr, pt = T6 / tree / (stem + "_R") / "pc.pc", T6 / tree / (stem + "_TS") / "pc.pc"
        if not pt.is_file():
            die("%s has a reactant arm but no transition-state arm" % stem)
        if pr.read_bytes() != pt.read_bytes():
            die("%s carries different charges at R and TS; one vector per condition would be a\n"
                "            fiction. Write both arms before using this file." % stem)
        rows = read_pc(pr)
        qs = [r[1] for r in rows]
        for i, (qtxt, q, x, y, z) in enumerate(rows, 1):
            near = [s for s in sites
                    if math.dist((float(s["x"]), float(s["y"]), float(s["z"])), (x, y, z)) <= TOL]
            if len(near) != 1:
                die("point %d of %s matches %d sites of %s within %g A; a charge vector is only\n"
                    "            as good as the site it sits on." % (i, stem, len(near), SITES, TOL))
            s = near[0]
            det.append("\t".join((label, diag, stem, tree.split("/")[0], str(i), s["site"],
                                  s["layer"], s["on_atom"],
                                  "%.6f" % x, "%.6f" % y, "%.6f" % z, qtxt,
                                  s["clearance_A"], s["q_min_neg"], s["q_min_pos"])))
        summ.append("\t".join((label, diag, stem, tree.split("/")[0], str(len(qs)),
                               "%.4f" % sum(abs(q) for q in qs),
                               "%.4f" % max(abs(q) for q in qs),
                               "%.6f" % min(qs), "%.6f" % max(qs),
                               str(sum(1 for q in qs if q > 0)),
                               str(sum(1 for q in qs if q < 0)), readable,
                               str(pr.relative_to(res)), sha(pr))))
    head = ("# Stage 6 per-site charge vectors, written by cs6_charges.py from each bare arm's\n"
            "# committed pc.pc. sites_f42436.tsv sha256 %s. Coordinates in angstrom,\n"
            "# charges in e, exactly as the ORCA inputs carry them. Every figure derived from these\n"
            "# vectors is read from this file or from pc.pc, never from a summary table or memory.\n"
            % ssha)
    out.mkdir(parents=True, exist_ok=True)
    (out / DETAIL).write_text(head + "\n".join(det) + "\n")
    (out / SUMMARY).write_text(head + "\n".join(summ) + "\n")
    sys.stdout.write("   %s  %d rows\n   %s  %d rows\n"
                     % (DETAIL, len(det) - 1, SUMMARY, len(summ) - 1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
