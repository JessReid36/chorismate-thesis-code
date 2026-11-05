#!/usr/bin/env python3
"""Stage 6 DIAGNOSTIC 27's reader.

Reads the 24 committed ORCA outputs, computes the three control barriers, the nine ddE values and
the two gap sequences, applies the decision rules D27.4 fixed BEFORE the run, and renders the
result document. It computes nothing that the pre-registration did not name, and it chooses no
threshold: TOL comes from D27.4, which took it from criterion D5.3(d).

The gates are D27.3's, in its order. Gate 7 - that every figure in the rendered document matches
the run - is NOT checked here, on purpose: it is checked by verify.py, which parses the outputs
with its own code so that the generator does not grade itself.
"""
import argparse, hashlib, json, re, sys
from pathlib import Path

H = 627.5094740631              # Eh -> kcal/mol, this project's constant
TOL = 1.0                       # kcal/mol, from D27.4, which took it from D5.3(d)
GATE_KCAL = 1.0e-4              # D27.3 GATE 3/4
GATE_EH = 1.0e-6                # D27.3 GATE 4, the house figure
PROP_EH = 1.0e-9                # job.out vs job.property.txt, two prints of one number

BASES = ["def2-SVP", "def2-SVPD", "def2-TZVPD"]
NBASIS = {"def2-SVP": 264, "def2-SVPD": 402}        # ORCA's own, from committed runs
NBASIS_PREDICTED = {"def2-TZVPD": 682}              # pyscf's count; recorded, not asserted
VECTORS = ["1.3000e_n10", "2.8800e_n10", "5.6433e_n20"]
SUMQ = {"1.3000e_n10": 1.3000, "2.8800e_n10": 2.8800, "5.6433e_n20": 5.6433}
NCH = {"1.3000e_n10": 10, "2.8800e_n10": 10, "5.6433e_n20": 20}
# the committed def2-SVPD anchors this run must reproduce
D16_CTL_BARRIER = 18.332256
D16_DDE_1300    = -6.134121
D16_CTL_R, D16_CTL_TS = -836.249865855127, -836.220651545264


class Gates:
    def __init__(self): self.rows, self.failed = [], []
    def __call__(self, ok, gate, msg):
        self.rows.append("   %s  %-7s %s" % ("pass" if ok else "FAIL", gate, msg))
        if not ok: self.failed.append("%s  %s" % (gate, msg))
        return ok
    def report(self): return "\n".join(self.rows)


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""): h.update(b)
    return h.hexdigest()


def read_out(d):
    """Everything this diagnostic needs from one run, and nothing else."""
    t = (d / "job.out").read_text(errors="replace")
    fsp = re.findall(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", t)
    nb = re.findall(r"Number of basis functions\s*\.*\s*(\d+)", t)
    nel = re.findall(r"Number of Electrons\s+NEL\s*\.*\s*(\d+)", t)
    scf = re.search(r"SCF CONVERGED AFTER\s+(\d+)\s+CYCLES", t)
    prop = None
    pf = d / "job.property.txt"
    if pf.exists():
        m = re.findall(r"&TOTALENERGY\s*\[&Type \"Double\"\]\s*(-?[\d.]+e[+-]\d+)", pf.read_text(errors="replace"), re.I)
        if m: prop = float(m[-1])
    return {"fsp": float(fsp[-1]) if fsp else None,
            "nbasis": sorted(set(int(x) for x in nb)),
            "nel": sorted(set(int(x) for x in nel)),
            "cycles": int(scf.group(1)) if scf else None,
            "normal": "****ORCA TERMINATED NORMALLY****" in t,
            "prop": prop}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tree", required=True, help="the d27 tree carrying the 24 run directories")
    ap.add_argument("--out", required=True, help="where to write D27_REPORT.txt, d27.tsv and the result document")
    ap.add_argument("--result-name", required=True)
    a = ap.parse_args()
    tree, out = Path(a.tree), Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    g = Gates()
    names = ["ctl_R", "ctl_TS"] + ["%s_%s" % (v, s) for v in VECTORS for s in ("R", "TS")]

    # ---- read everything first; read nothing into a conclusion yet --------------
    R = {}
    for b in BASES:
        for n in names:
            d = tree / b / n
            if not (d / "job.out").exists():
                g(False, "GATE 5", "%s/%s has no job.out - the job did not produce it" % (b, n)); continue
            R[(b, n)] = read_out(d)

    # ---- GATE 0  provenance ----------------------------------------------------
    man = tree / "HPC_MANIFEST.sha256"
    inp = tree / "INPUT_CHECKSUMS.sha256"
    if g(man.exists(), "GATE 0", "hpc1's own HPC_MANIFEST.sha256 is present") and \
       g(inp.exists(), "GATE 0", "the committed INPUT_CHECKSUMS.sha256 is present"):
        mh = dict((n, h) for h, n in (l.split(None, 1) for l in man.read_text().split("\n") if l.strip()))
        mh = dict((v.strip(), k) for k, v in mh.items()) if False else \
             dict((l.split(None, 1)[1].strip().lstrip("./"), l.split(None, 1)[0])
                  for l in man.read_text().split("\n") if l.strip())
        bad = []
        for l in inp.read_text().split("\n"):
            if not l.strip(): continue
            h, n = l.split(None, 1); n = n.strip()
            if n == "INPUT_CHECKSUMS.sha256": continue
            if mh.get(n) != h: bad.append(n)
        n_checked = sum(1 for l in inp.read_text().split("\n")
                        if l.strip() and l.split(None, 1)[1].strip() != "INPUT_CHECKSUMS.sha256")
        g(not bad, "GATE 0", "all %d committed inputs are byte-identical to what hpc1 ran%s"
          % (n_checked, "" if not bad else "; differs: " + ", ".join(bad[:3])))

    # ---- GATE 0b  hpc1 held nothing the record does not account for ------------
    if man.exists():
        KNOWN_IN = {"job.inp", "substrate.xyz", "pc.pc"}
        KNOWN_OUT = {"job.out", "job.property.txt", "job.gbw", "job.densities",
                     "job.densitiesinfo", "job.bibtex", "job.cpcm", "job.cpcm_corr"}
        TOP = {"run.pbs", "run.pbs.out", "INPUT_CHECKSUMS.sha256", "HPC_MANIFEST.sha256"}
        stray, ctl_pc = [], []
        for l in man.read_text().split("\n"):
            if not l.strip(): continue
            n = l.split(None, 1)[1].strip().lstrip("./")
            parts = n.split("/")
            if len(parts) == 1:
                if parts[0] not in TOP: stray.append(n)
                continue
            if len(parts) != 3: stray.append(n); continue
            b, run, f = parts
            if b not in BASES or run not in names: stray.append(n); continue
            if not (f in KNOWN_IN or f in KNOWN_OUT or f.endswith(".tmp")): stray.append(n); continue
            if run.startswith("ctl_") and f == "pc.pc": ctl_pc.append(n)
        g(not stray, "GATE 0b", "every file hpc1 held is an expected input, output or scratch file%s"
          % ("" if not stray else "; unexpected: " + ", ".join(stray[:4])))
        g(not ctl_pc, "GATE 0b", "no control directory on hpc1 held a pc.pc%s"
          % ("" if not ctl_pc else "; found: " + ", ".join(ctl_pc)))

    # ---- GATE 1  the groups differ only in the basis ---------------------------
    diffs = []
    for n in names:
        ref = (tree / BASES[0] / n / "job.inp").read_text().replace(BASES[0], "<B>")
        for b in BASES[1:]:
            if (tree / b / n / "job.inp").read_text().replace(b, "<B>") != ref: diffs.append("%s/%s" % (b, n))
    g(not diffs, "GATE 1", "all 24 inputs are identical once the basis token is removed%s"
      % ("" if not diffs else "; differ: " + ", ".join(diffs[:3])))

    # ---- GATE 2  control vs charged --------------------------------------------
    bad = []
    for b in BASES:
        c = (tree / b / "ctl_R/job.inp").read_text()
        for v in VECTORS:
            for s in ("R", "TS"):
                ref = (tree / b / ("ctl_%s" % s) / "job.inp").read_text()
                if (tree / b / ("%s_%s" % (v, s)) / "job.inp").read_text().replace(
                        '%pointcharges "pc.pc"\n', "") != ref: bad.append("%s/%s_%s" % (b, v, s))
        if (tree / b / "ctl_R/pc.pc").exists() or (tree / b / "ctl_TS/pc.pc").exists():
            bad.append("%s control has a pc.pc" % b)
    g(not bad, "GATE 2", "each charged input is its control plus exactly the %%pointcharges line%s"
      % ("" if not bad else "; wrong: " + ", ".join(bad[:3])))

    # ---- the quantities --------------------------------------------------------
    ctlb, barr, ddE = {}, {}, {}
    for b in BASES:
        if (b, "ctl_R") not in R or (b, "ctl_TS") not in R: continue
        ctlb[b] = (R[(b, "ctl_TS")]["fsp"] - R[(b, "ctl_R")]["fsp"]) * H
        for v in VECTORS:
            if (b, v + "_R") not in R or (b, v + "_TS") not in R: continue
            barr[(b, v)] = (R[(b, v + "_TS")]["fsp"] - R[(b, v + "_R")]["fsp"]) * H
            ddE[(b, v)] = barr[(b, v)] - ctlb[b]

    # ---- GATE 3  the committed production number reproduces --------------------
    have = ("def2-SVPD", "1.3000e_n10") in ddE
    if g(have, "GATE 3", "the def2-SVPD 1.3000e arm was read"):
        d = ddE[("def2-SVPD", "1.3000e_n10")]
        g(abs(d - D16_DDE_1300) <= GATE_KCAL, "GATE 3",
          "def2-SVPD 1.3000e ddE %.6f reproduces d16's committed %.6f (|d| %.2e <= %.0e)"
          % (d, D16_DDE_1300, abs(d - D16_DDE_1300), GATE_KCAL))

    # ---- GATE 4  the control barrier reproduces --------------------------------
    if "def2-SVPD" in ctlb:
        g(abs(ctlb["def2-SVPD"] - D16_CTL_BARRIER) <= GATE_KCAL, "GATE 4",
          "def2-SVPD control barrier %.6f reproduces d16's committed %.6f" % (ctlb["def2-SVPD"], D16_CTL_BARRIER))
        for n, ref in (("ctl_R", D16_CTL_R), ("ctl_TS", D16_CTL_TS)):
            v = R[("def2-SVPD", n)]["fsp"]
            g(abs(v - ref) <= GATE_EH, "GATE 4",
              "def2-SVPD %-7s %.12f reproduces d16's %.12f to %.0e Eh" % (n, v, ref, GATE_EH))

    # ---- GATE 5  termination and convergence -----------------------------------
    nn = [k for k, v in R.items() if not v["normal"]]
    nc = [k for k, v in R.items() if v["cycles"] is None]
    g(len(R) == 24, "GATE 5", "all 24 runs produced an output (%d read)" % len(R))
    g(not nn, "GATE 5", "every run reports ORCA TERMINATED NORMALLY%s"
      % ("" if not nn else "; did not: " + ", ".join("/".join(k) for k in nn[:3])))
    g(not nc, "GATE 5", "every run reports SCF CONVERGED AFTER n CYCLES%s"
      % ("" if not nc else "; did not: " + ", ".join("/".join(k) for k in nc[:3])))

    # ---- GATE 6  the basis is the one asked for --------------------------------
    badnel = [k for k, v in R.items() if v["nel"] != [118]]
    g(not badnel, "GATE 6", "every run reports Number of Electrons NEL = 118%s"
      % ("" if not badnel else "; not: " + ", ".join("/".join(k) for k in badnel[:3])))
    tz = set()
    for b in BASES:
        got = set(tuple(v["nbasis"]) for k, v in R.items() if k[0] == b)
        if b in NBASIS:
            g(got == {(NBASIS[b],)}, "GATE 6", "%-11s every run reports %d basis functions" % (b, NBASIS[b]))
        else:
            g(len(got) == 1, "GATE 6", "%-11s all eight runs agree on the basis count" % b)
            if len(got) == 1:
                tz = list(got)[0]
                g(tz[0] > NBASIS["def2-SVPD"], "GATE 6",
                  "%-11s RECORDED at %d basis functions, above def2-SVPD's 402%s" % (
                      b, tz[0], "" if tz[0] == NBASIS_PREDICTED[b]
                      else "; the document PREDICTED %d, so the prediction was wrong - this changes no quantity"
                           % NBASIS_PREDICTED[b]))

    # ---- GATE 8  two independent prints of each energy -------------------------
    pb = [("/".join(k), abs(v["fsp"] - v["prop"])) for k, v in R.items()
          if v["prop"] is not None and abs(v["fsp"] - v["prop"]) > PROP_EH]
    npr = [k for k, v in R.items() if v["prop"] is None]
    g(not npr, "GATE 8", "every run has a job.property.txt to check job.out against%s"
      % ("" if not npr else "; missing: %d" % len(npr)))
    g(not pb, "GATE 8", "job.out and job.property.txt agree to %.0e Eh in all %d runs%s"
      % (PROP_EH, len(R) - len(npr), "" if not pb else "; worst " + str(pb[0])))

    # ---- the decision, by D27.4's rules, which were fixed before the run -------
    rows, verdicts = [], {}
    for v in VECTORS:
        if any((b, v) not in ddE for b in BASES): continue
        a, bb, c = ddE[(BASES[0], v)], ddE[(BASES[1], v)], ddE[(BASES[2], v)]
        gsd, gdt = abs(bb - a), abs(c - bb)
        not_converging = gdt >= gsd
        not_converged = gdt > TOL
        verdict = "NOT CONVERGING" if not_converging else ("NOT CONVERGED" if not_converged else "CONVERGED")
        verdicts[v] = verdict
        rows.append({"vector": v, "n": NCH[v], "sum_abs_q": SUMQ[v],
                     "ddE_def2_SVP": a, "ddE_def2_SVPD": bb, "ddE_def2_TZVPD": c,
                     "gap_SD": gsd, "gap_DT": gdt,
                     "shrinking": not not_converging, "within_TOL": not not_converged,
                     "verdict": verdict})
    # the third vector's arbitration, stated in D27.4 before the run
    arb = None
    v3 = VECTORS[2]
    if all((b, v3) in ddE for b in BASES):
        a, bb, c = ddE[(BASES[0], v3)], ddE[(BASES[1], v3)], ddE[(BASES[2], v3)]
        if abs(c - bb) <= TOL:   arb = "def2-SVP IS THE ARTEFACT"
        elif abs(c - a) <= TOL:  arb = "def2-SVPD IS THE ARTEFACT"
        else:                    arb = "BOTH ARE ARTEFACTS"
    overall = (list(verdicts.values())[0] if verdicts and len(set(verdicts.values())) == 1
               else "PER VECTOR; NO SINGLE VERDICT IS FORMED")

    # ---- write the machine outputs --------------------------------------------
    hdr = ["vector", "n", "sum_abs_q", "ddE_def2_SVP", "ddE_def2_SVPD", "ddE_def2_TZVPD",
           "gap_SD", "gap_DT", "shrinking", "within_TOL", "verdict"]
    with open(out / "d27.tsv", "w") as f:
        f.write("\t".join(hdr) + "\n")
        for r in rows:
            f.write("\t".join(("%.6f" % r[k] if isinstance(r[k], float) else str(r[k])) for k in hdr) + "\n")
    payload = {"H": H, "TOL": TOL, "control_barrier": ctlb, "nbasis_tzvpd": (tz[0] if tz else None),
               "ddE": dict(("%s|%s" % k, v) for k, v in ddE.items()),
               "fsp": dict(("%s|%s" % k, v["fsp"]) for k, v in R.items()),
               "cycles": dict(("%s|%s" % k, v["cycles"]) for k, v in R.items()),
               "rows": rows, "arbitration": arb, "overall": overall,
               "gates_failed": g.failed}
    (out / "d27_values.json").write_text(json.dumps(payload, indent=1, sort_keys=True))
    with open(out / "D27_REPORT.txt", "w") as f:
        f.write("DIAGNOSTIC 27: the basis test, read from the 24 committed outputs\n")
        f.write("TOL = %.1f kcal/mol, from D27.4, which took it from criterion D5.3(d)\n\n" % TOL)
        f.write("GATES\n%s\n\n" % g.report())
        f.write("CONTROL BARRIERS (bare substrate, kcal/mol)\n")
        for b in BASES:
            if b in ctlb: f.write("   %-11s %12.6f\n" % (b, ctlb[b]))
        f.write("\nddE BY VECTOR AND BASIS (kcal/mol)\n")
        f.write("   %-12s %5s %9s %12s %12s %12s %10s %10s  %s\n" % (
            "vector", "n", "sum|q|", "def2-SVP", "def2-SVPD", "def2-TZVPD", "gap_SD", "gap_DT", "verdict"))
        for r in rows:
            f.write("   %-12s %5d %9.4f %12.6f %12.6f %12.6f %10.6f %10.6f  %s\n" % (
                r["vector"], r["n"], r["sum_abs_q"], r["ddE_def2_SVP"], r["ddE_def2_SVPD"],
                r["ddE_def2_TZVPD"], r["gap_SD"], r["gap_DT"], r["verdict"]))
        f.write("\nTHIRD VECTOR ARBITRATION (D27.4): %s\n" % arb)
        f.write("OVERALL: %s\n" % overall)
        if g.failed:
            f.write("\nFAILED GATES, committed as failures:\n")
            for x in g.failed: f.write("   %s\n" % x)
    print(g.report()); print()
    print("control barriers: " + ", ".join("%s %.6f" % (b, ctlb[b]) for b in BASES if b in ctlb))
    for r in rows:
        print("   %-12s SVP %10.6f  SVPD %10.6f  TZVPD %10.6f   gap_SD %8.6f  gap_DT %8.6f  %s"
              % (r["vector"], r["ddE_def2_SVP"], r["ddE_def2_SVPD"], r["ddE_def2_TZVPD"],
                 r["gap_SD"], r["gap_DT"], r["verdict"]))
    print("   third vector: %s" % arb)
    print("   OVERALL: %s" % overall)
    return 2 if g.failed else 0


if __name__ == "__main__":
    sys.exit(main())
