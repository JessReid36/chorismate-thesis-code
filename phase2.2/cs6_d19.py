#!/usr/bin/env python3
"""
cs6_d19.py - STAGE6_DIAGNOSTIC19_20261030.txt.

Is the near-degeneracy at full field multireference, or a small gap B3LYP can still be trusted
with? Two indicators: FOD (TPSS/def2-TZVP, 5000 K) and SCF stability at the production method.

  python3 cs6_d19.py prepare <results_repo> <OUT> <HPC_DIR>
  python3 cs6_d19.py analyse <results_repo> <RUNS> <OUT_DIR>

NOTHING ALREADY COMMITTED IS RETYPED OR RE-RUN.
  cs6_transfer    design, head, read_xyz
  cs6_additivity  members
  cs13_qscale     FRAME_SHA, GEOM, read_xyz_raw
The point-charge files are COPIED from the committed field tree, not regenerated, so this ladder
is the same ladder the gap was measured on.

THE STABILITY PARSER IS WRITTEN AGAINST DOCUMENTATION THAT DOES NOT QUOTE ORCA'S OUTPUT.
If it cannot find the block it STOPS and prints what it did find, rather than guessing a verdict.
Re-reading costs no HPC time.
"""
import hashlib, re, sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import cs13_qscale as Q13
import cs6_transfer as TR
import cs6_additivity as AD

PREREG = "STAGE6_DIAGNOSTIC19_20261030.txt"
NAME, NCH = "S6-20", 20
LAMS = (0.00, 0.25, 0.50, 0.75, 1.00)          # D19.2 A
STAB_LAMS = (0.00, 1.00)                        # D19.2 B
PROD_BASIS = "def2-SVP"
TAG = {0.00: "lam000", 0.25: "lam025", 0.50: "lam050", 0.75: "lam075", 1.00: "lam100"}
FIELD_PC = "field/def2-SVP/%s_%s/pc.pc"
R_RISE, R_FLAT = 2.0, 1.25                      # D19.3 conventions, fixed before any run
NELEC = 118


def stop(m):
    print("STOP: " + m); sys.stdout.flush(); sys.exit(1)


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def design_q(T6, lam):
    d = TR.design(T6, NAME, NCH)
    return AD.members(d, d["site"], lam)


def gate_charges(T6, lam, st):
    """GATE 2: the committed field tree's charges ARE this ladder's charges."""
    p = T6 / (FIELD_PC % (TAG[lam], st))
    if not p.is_file():
        stop("the committed charge file is missing: %s" % p)
    rows = [l.split() for l in p.read_text().strip().splitlines()]
    if len(rows[0]) == 1:
        rows = rows[1:]
    mem = design_q(T6, lam)
    if len(rows) != NCH or len(mem) != NCH:
        stop("%s carries %d charges and the design has %d, expected %d"
             % (p, len(rows), len(mem), NCH))
    # Compared AT THE PRECISION THE FILE STORES (%.6f for charges, %.10f for coordinates),
    # not by float equality. lambda 0.75 is not an exact binary fraction, so 0.15 * 0.75 is
    # 0.11249999999999999 while the file's "0.112500" parses to exactly 0.1125 - a difference of
    # one unit in the last place. Bit equality would reject a file that is correct to every digit
    # it records. Anything genuinely different at 5e-7 or more is still caught.
    for (_, q, X), r in zip(mem, rows):
        if "%.6f" % q != "%.6f" % float(r[0]) or \
           any("%.10f" % a != "%.10f" % float(b) for a, b in zip(X, r[1:4])):
            stop("%s does not match the design at lambda %.2f, at the precision the file stores. "
                 "This ladder would not be the ladder the gap was measured on." % (p, lam))
    return p


def fod_head():
    return "! FOD\n%maxcore 3000\n%pal nprocs 8 end\n"


def stab_head():
    """The committed production head, with the stability request added. Not retyped."""
    h = TR.head(PROD_BASIS)
    old = "%scf MaxIter 300 end\n"
    if old not in h:
        stop("cs6_transfer.head() no longer carries %r; the stability head cannot be derived "
             "from it without retyping, which this does not do." % old)
    return h.replace(old, "%scf\n  MaxIter 300\n  STABPerform true\nend\n")


def write_run(d, head, geo_text, pc_src):
    d.mkdir(parents=True)
    (d / "substrate.xyz").write_text(geo_text)
    body = head
    if pc_src is not None:
        (d / "pc.pc").write_bytes(Path(pc_src).read_bytes())   # copied, never regenerated
        body += '%pointcharges "pc.pc"\n'
    body += "* xyzfile -2 1 substrate.xyz\n"
    (d / "job.inp").write_text(body)


def rn(kind, lam, st):
    return "%s_%s_%s" % (kind, TAG[lam], st)


def prepare(root, out, hpc):
    T6 = Path(root) / "phase2.2" / "cs_stage6"
    out = Path(out)
    if out.exists():
        stop("%s already exists; this never overwrites" % out)
    geo = {}
    for st, fn in (("R", "reactant_qm.xyz"), ("TS", "transition_state_qm.xyz")):
        p = T6 / "frame_42436" / fn
        if sha(p) != Q13.FRAME_SHA[fn]:
            stop("%s is not the committed frame" % p)
        geo[st] = p.read_text()

    runs = []
    for lam in LAMS:
        for st in Q13.GEOM:
            pc = None if lam == 0.00 else gate_charges(T6, lam, st)
            write_run(out / rn("fod", lam, st), fod_head(), geo[st], pc)
            runs.append(rn("fod", lam, st))
    for lam in STAB_LAMS:
        for st in Q13.GEOM:
            pc = None if lam == 0.00 else gate_charges(T6, lam, st)
            write_run(out / rn("stab", lam, st), stab_head(), geo[st], pc)
            runs.append(rn("stab", lam, st))
    if len(runs) != 14:
        stop("%d runs were built, expected the 14 D19.2 states" % len(runs))

    pbs = AD.PBS % dict(tag="19", hpc=hpc, basis=".", n=len(runs),
                        runs="\n".join('run "%s"' % r for r in runs))
    pbs = pbs.replace("cm_d1419", "cm_d19") \
             .replace("DIAGNOSTIC 14 (cs6_additivity.py)", "DIAGNOSTIC 19 (cs6_d19.py)") \
             .replace("%d single points at .." % len(runs),
                      "%d single points: 10 FOD (TPSS/def2-TZVP, 5000 K) and 4 stability "
                      "(%s)." % (len(runs), PROD_BASIS)) \
             .replace("%s/./" % hpc, "%s/" % hpc).replace("W=%s/." % hpc, "W=%s" % hpc)
    if re.search(r"^#PBS([ \t].*)?[ \t]-q([ \t]|$)", pbs, re.M):
        stop("the PBS script requests a queue; the allocation is automatic")
    (out / "run.pbs").write_text(pbs)

    print("DIAGNOSTIC 19 prepared at %s" % out)
    print("   FOD    TPSS/def2-TZVP, 5000 K, lambda %s, R and TS   (10 runs)"
          % ", ".join("%.2f" % x for x in LAMS))
    print("   STAB   %s + STABPerform, lambda %s, R and TS        (4 runs)"
          % (PROD_BASIS, ", ".join("%.2f" % x for x in STAB_LAMS)))
    print("   gate   every pc.pc is COPIED from the committed field tree, byte for byte, and")
    print("          matches the design at the precision the file stores")
    print("   runs   %s" % ", ".join(runs))
    return 0


def read_fod(p):
    t = Path(p).read_text(errors="replace")
    if "ORCA TERMINATED NORMALLY" not in t:
        stop("%s: ORCA did not terminate normally" % p)
    m = re.findall(r"N_FOD\s*[=:]\s*(-?\d+\.\d+)", t)
    if not m:
        stop("%s: no 'N_FOD =' line. The FOD analysis did not run or ORCA's label has changed."
             % p)
    return float(m[-1]), t


def read_stab(p):
    """ORCA's printed wording is not quoted in the documentation this was written against.
    Several spellings are tried; if none matches, this STOPS and shows what is there."""
    t = Path(p).read_text(errors="replace")
    if "ORCA TERMINATED NORMALLY" not in t:
        stop("%s: ORCA did not terminate normally" % p)
    low = t.lower()
    if "stability analysis" not in low:
        stop("%s: no stability-analysis section. Was STABPerform honoured?" % p)
    k = low.find("stability analysis")
    blk = t[k:k + 4000]
    unstable = bool(re.search(r"\bis unstable\b|\bunstable\b|instabilit", blk, re.I))
    stable = bool(re.search(r"\bis stable\b|\bwavefunction is stable\b", blk, re.I))
    eig = [float(x) for x in re.findall(r"(-?\d+\.\d{3,})", blk)]
    if not (unstable or stable or eig):
        stop("%s: the stability section was found but neither a verdict nor an eigenvalue could "
             "be read from it. THIS IS A PARSER PROBLEM, NOT A RUN PROBLEM - the run is kept and "
             "re-reading costs nothing. The section reads:\n%s" % (p, blk[:1500]))
    lowest = min(eig) if eig else None
    if unstable and not stable:
        v = "UNSTABLE"
    elif stable and not unstable:
        v = "STABLE"
    elif lowest is not None:
        v = "UNSTABLE" if lowest < 0 else "STABLE"
    else:
        stop("%s: the stability section is ambiguous and carries no eigenvalue:\n%s"
             % (p, blk[:1500]))
    return v, lowest


def analyse(root, runs, outdir):
    T6 = Path(root) / "phase2.2" / "cs_stage6"
    runs = Path(runs)
    for lam in LAMS:
        for st in Q13.GEOM:
            if lam > 0:
                gate_charges(T6, lam, st)

    F, S = {}, {}
    for lam in LAMS:
        for st in Q13.GEOM:
            d = runs / rn("fod", lam, st)
            n, t = read_fod(d / "job.out")
            nel = re.search(r"Number of Electrons\s+NEL\s+\.+\s+(\d+)", t)
            if not nel or int(nel.group(1)) != NELEC:
                stop("%s: %s electrons, expected %d" % (d, nel.group(1) if nel else "none", NELEC))
            if lam > 0 and not re.search(r"ok \(%d point charges\)" % NCH, t):
                stop("GATE 3: %s does not report ORCA reading %d point charges. A flat ladder "
                     "here would be flat for the wrong reason." % (d, NCH))
            if lam == 0 and re.search(r"Reading point charge file", t):
                stop("%s: the lambda 0 run reports a point charge file and must not have one" % d)
            F[(lam, st)] = n
    for lam in STAB_LAMS:
        for st in Q13.GEOM:
            S[(lam, st)] = read_stab(runs / rn("stab", lam, st) / "job.out")

    L = []; w = L.append
    w("Stage 6, DIAGNOSTIC 19 RESULT.")
    w("Read by cs6_d19.py against STAGE6_DIAGNOSTIC19_20261030.txt. Criteria as pre-registered.")
    w("")
    w("GATE 1, the zero-field stability controls")
    g1 = all(S[(0.00, st)][0] == "STABLE" for st in Q13.GEOM)
    for st in Q13.GEOM:
        v, e = S[(0.00, st)]
        w("   stab_lam000_%-2s  %-9s lowest eigenvalue %s" % (st, v, "%.6f" % e if e is not None else "n/a"))
    if not g1:
        w("   GATE 1 FAILED. A charge-free closed-shell substrate reports an instability, so the")
        w("   indicator is not measuring what it is meant to and nothing below is readable.")
        Path(outdir).mkdir(parents=True, exist_ok=True)
        (Path(outdir) / "DIAGNOSTIC19_RESULT.txt").write_text("\n".join(L) + "\n")
        print("\n".join(L)); return 1
    w("   GATE 1 passed.")
    w("")
    w("READING 1, does FOD rise with the field?   (TPSS/def2-TZVP, 5000 K; all systems 118 e)")
    w("   %-8s %12s %12s" % ("lambda", "N_FOD R", "N_FOD TS"))
    for lam in LAMS:
        w("   %-8.2f %12.6f %12.6f" % (lam, F[(lam, "R")], F[(lam, "TS")]))
    ratios = {}
    for st in Q13.GEOM:
        z = F[(0.00, st)]
        ratios[st] = (F[(1.00, st)] / z) if z else None
    mono = {st: all(F[(LAMS[i + 1], st)] >= F[(LAMS[i], st)] for i in range(len(LAMS) - 1))
            for st in Q13.GEOM}
    for st in Q13.GEOM:
        w("   R_FOD %-3s %s   monotonic in lambda: %s"
          % (st, ("%.3f" % ratios[st]) if ratios[st] is not None else "undefined (N_FOD at zero field is 0)",
             "yes" if mono[st] else "no"))
    if all(r is not None and r >= R_RISE for r in ratios.values()):
        v1 = "RISES WITH FIELD"
    elif all(r is not None and r <= R_FLAT for r in ratios.values()):
        v1 = "FLAT"
    else:
        v1 = "WEAK OR SPLIT"
    w("   thresholds, conventions fixed in D19.3: rises >= %.2f, flat <= %.2f" % (R_RISE, R_FLAT))
    w("   VERDICT: %s" % v1)
    w("")
    w("READING 2, is the production wavefunction stable at full field?")
    for st in Q13.GEOM:
        v, e = S[(1.00, st)]
        w("   stab_lam100_%-2s  %-9s lowest eigenvalue %s" % (st, v, "%.6f" % e if e is not None else "n/a"))
    v2 = "UNSTABLE AT FULL FIELD" if any(S[(1.00, st)][0] == "UNSTABLE" for st in Q13.GEOM) \
         else "STABLE AT FULL FIELD"
    w("   VERDICT: %s" % v2)
    w("")
    w("READING 3, do the two indicators agree?")
    if v1 == "RISES WITH FIELD" and v2 == "UNSTABLE AT FULL FIELD":
        v3 = "AGREE, MULTIREFERENCE"
    elif v1 == "FLAT" and v2 == "STABLE AT FULL FIELD":
        v3 = "AGREE, SINGLE REFERENCE"
    else:
        v3 = "DISAGREE"
    w("   VERDICT: %s" % v3)
    w("")
    w("D19.4 stands whatever was returned. ORCA's own documentation states that FOD is not a")
    w("fully reliable diagnostic for multireference character and should be a first indication")
    w("combined with more sophisticated diagnostics. A positive result here licenses a CASSCF or")
    w("NEVPT2 calculation; it does not substitute for one. No energy in this diagnostic enters")
    w("any Dset, d_shield or barrier, and no D6.3 class is assigned to a def2-TZVP run.")

    outdir = Path(outdir); outdir.mkdir(parents=True, exist_ok=True)
    p = outdir / "DIAGNOSTIC19_RESULT.txt"
    p.write_text("\n".join(L) + "\n")
    print("\n".join(L)); print("\nwritten: %s" % p)
    return 0


if __name__ == "__main__":
    if len(sys.argv) == 5 and sys.argv[1] == "prepare":
        raise SystemExit(prepare(sys.argv[2], sys.argv[3], sys.argv[4]))
    if len(sys.argv) == 5 and sys.argv[1] == "analyse":
        raise SystemExit(analyse(sys.argv[2], sys.argv[3], sys.argv[4]))
    print(__doc__); raise SystemExit(2)
