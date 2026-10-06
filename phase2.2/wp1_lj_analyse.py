#!/usr/bin/env python3
"""
wp1_lj_analyse.py - WP1: which Lennard-Jones combining rule ORCA 6.0.1 applies to ORCAFF.prms
parameters, and what the length column means. Reads the 648 two-atom MM single points written by
wp1_lj.pbs and tests them against the criteria that job wrote before it ran (WP1_CRITERIA.txt).

Format read (seen in the job's output before this was written):
    FINAL SINGLE POINT ENERGY (MM)       -0.000214922924      (Eh; the last such line is used)

Candidate rules, Rc = length column, e = |epsilon column| (kcal/mol), E(r) = e_ij[(R_ij/r)^12 - 2(R_ij/r)^6]:
    A  R_ij = (Rc_i + Rc_j)/2              (column = full R_min; Amber's R*_i + R*_j)
    B  R_ij = sqrt(Rc_i Rc_j)              (geometric mean of R_min)
    C  R_ij = Rc_i + Rc_j                  (column = R_min/2)
    D  R_ij = 2^(1/6) (Rc_i + Rc_j)/2      (column = sigma)
    all with e_ij = sqrt(e_i e_j)

Criteria (WP1_CRITERIA.txt):
    1. switched off: exactly one rule reproduces every energy of all four pairs within 1e-8 Eh;
       every other rule misses by more than 1e-4 Eh somewhere in 2-6 A
    2. default switching: energies = switched-off energies + one constant per pair, within 1e-8 Eh
    3. all 648 single points terminate normally

AMENDMENT to criterion 1, made before the full analysis ran. The one energy seen first (Na-Oc,
3.00 A, switching off) matches rule A to 4.7e-10 Eh, a relative 2.2e-6. A relative difference of that
size is what a different kcal/mol-to-Eh constant inside ORCA would give, and at 2.0 A, where the
energies reach ~0.2 Eh, it alone would exceed the absolute 1e-8 Eh bound even for the right rule. The
bound mixed a units question into a shape question. Criterion 1 is therefore tested as:
    1a. per rule, fit ONE scale factor s over all switched-off energies of all four pairs;
        the right rule leaves residuals within 1e-6 |E| + 1e-10 Eh at every point, and s itself is
        reported (it measures ORCA's energy-unit constant against 627.5094740631);
    1b. every other rule leaves a residual above 1e-3 |E| at some point with |E| > 1e-5 Eh.
The original absolute test is still computed and printed for the record.

Writes wp1_energies.tsv and WP1_REPORT.txt in the WP1 folder; exits non-zero if a criterion fails.

USAGE   python3 wp1_lj_analyse.py [wp1_folder]
Python 3.6 compatible, standard library only (runs on the hpc1 login node).
"""
import math
import re
import sys
from pathlib import Path

H2KCAL = 627.5094740631
W = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("/home/18660916/system_development/phase2.2/wp1_lj_combining")
RULES = {
    "A": lambda a, b: (a + b) / 2.0,
    "B": lambda a, b: math.sqrt(a * b),
    "C": lambda a, b: a + b,
    "D": lambda a, b: 2.0 ** (1.0 / 6.0) * (a + b) / 2.0,
}
DESC = {"A": "R_ij = (Rc_i + Rc_j)/2 (column = full R_min)", "B": "R_ij = sqrt(Rc_i Rc_j)",
        "C": "R_ij = Rc_i + Rc_j (column = R_min/2)", "D": "column = sigma"}


def energy(p):
    t = p.read_text(errors="replace")
    if "ORCA TERMINATED NORMALLY" not in t:
        return None
    m = re.findall(r"FINAL SINGLE POINT ENERGY \(MM\)\s+(-?\d+\.\d+)", t)
    return float(m[-1]) if m else None


def lj(e, R, r):
    x = (R / r) ** 6
    return e * (x * x - 2.0 * x) / H2KCAL


def main():
    pairs = [l.split() for l in (W / "pairs.txt").read_text().splitlines() if l.strip()]
    dists = [l.strip() for l in (W / "distances.txt").read_text().splitlines() if l.strip()]
    out, fails, rep = [], [], []
    say = rep.append
    E = {}
    missing = 0
    for sw in ("default", "noswitch"):
        for name, e1, p1, r1, e2, p2, r2, mult in pairs:
            for r in dists:
                v = energy(W / sw / name / "sp_{}.out".format(r))
                if v is None:
                    missing += 1
                E[(sw, name, r)] = v
                out.append("{}\t{}\t{}\t{}".format(sw, name, r, "" if v is None else "{:.12f}".format(v)))
    (W / "wp1_energies.tsv").write_text("switching\tpair\tr_A\tE_MM_Eh\n" + "\n".join(out) + "\n")
    say("WP1 report - ORCA 6.0.1 Lennard-Jones combining rule (wp1_lj_analyse.py)")
    say("energies read: {} of {}".format(len(E) - missing, len(E)))
    if missing:
        fails.append("criterion 3: {} single points missing or not terminated normally".format(missing))

    say("\nCriterion 1 as first written, switching off: largest |E_ORCA - E_rule| over 2-6 A, Eh (for the record)")
    say("  pair   " + "".join("{:>14}".format(k) for k in RULES))
    worst = {k: 0.0 for k in RULES}
    pts = {k: [] for k in RULES}
    for name, e1, p1, r1, e2, p2, r2, mult in pairs:
        e = math.sqrt(abs(float(p1)) * abs(float(p2)))
        row = []
        for k, f in RULES.items():
            R = f(float(r1), float(r2))
            dd = []
            for r in dists:
                eo = E[("noswitch", name, r)]
                if eo is None:
                    continue
                er = lj(e, R, float(r))
                pts[k].append((eo, er))
                dd.append(abs(eo - er))
            worst[k] = max(worst[k], max(dd))
            row.append(max(dd))
        say("  {:<6} ".format(name) + "".join("{:>14.2e}".format(x) for x in row))

    say("\nCriterion 1 as amended: one fitted scale factor per rule, all pairs together")
    say("  rule   scale s            max residual (Eh)   max residual/|E| (|E|>1e-5)   within 1e-6|E|+1e-10")
    winners, res_rel = [], {}
    for k in RULES:
        P = pts[k]
        sc = sum(a * b for a, b in P) / sum(b * b for a, b in P)
        rabs = max(abs(a - sc * b) for a, b in P)
        rrel = max(abs(a - sc * b) / abs(a) for a, b in P if abs(a) > 1e-5)
        ok = all(abs(a - sc * b) <= 1e-6 * abs(a) + 1e-10 for a, b in P)
        res_rel[k] = rrel
        if ok:
            winners.append(k)
        say("  {:<6} {:.10f}   {:>12.2e}       {:>12.2e}                  {}".format(k, sc, rabs, rrel, "yes" if ok else "no"))
        if ok:
            say("         -> ORCA's energy constant implied: {:.6f} kcal/mol per Eh (this script uses {})".format(H2KCAL / sc, H2KCAL))
    losers_ok = all(res_rel[k] > 1e-3 for k in RULES if k not in winners)
    if len(winners) == 1 and losers_ok:
        k = winners[0]
        say("  PASS: rule {} - {} - is the only rule that fits; every other rule misses by >1e-3 |E|".format(k, DESC[k]))
    else:
        fails.append("criterion 1 (amended): rules that fit: {}; all others miss by >1e-3: {}".format(winners or "none", losers_ok))
        say("  FAIL: " + fails[-1])

    say("\nCriterion 2, default switching = switched off + one constant per pair")
    for name, e1, p1, r1, e2, p2, r2, mult in pairs:
        d = [E[("default", name, r)] - E[("noswitch", name, r)] for r in dists
             if E[("default", name, r)] is not None and E[("noswitch", name, r)] is not None]
        spread = max(d) - min(d)
        ok = spread <= 1e-8
        say("  {:<6} constant {:+.6e} Eh ({:+.6f} kcal/mol); spread over 2-6 A {:.1e} Eh  {}".format(
            name, sum(d) / len(d), sum(d) / len(d) * H2KCAL, spread, "PASS" if ok else "FAIL"))
        if not ok:
            fails.append("criterion 2: {} switching offset varies by {:.1e} Eh".format(name, spread))

    if winners:
        k = winners[0]
        say("\nPair minima under rule {} (the values every LJ distance in the plan should use):".format(k))
        for name, e1, p1, r1, e2, p2, r2, mult in pairs:
            R = RULES[k](float(r1), float(r2))
            say("  {:<6} {}-{}: R_min {:.4f} A, well depth {:.6f} kcal/mol".format(
                name, e1, e2, R, math.sqrt(abs(float(p1)) * abs(float(p2)))))
    say("\nRESULT: " + ("all criteria PASS" if not fails else "FAIL - " + "; ".join(fails)))
    (W / "WP1_REPORT.txt").write_text("\n".join(rep) + "\n")
    print("\n".join(rep))
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
