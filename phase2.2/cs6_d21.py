#!/usr/bin/env python3
"""
cs6_d21.py - STAGE6_DIAGNOSTIC21_20261031.txt.

Are the designs that carry every readable result sitting on stable solutions?

  python3 cs6_d21.py prepare <results_repo> <OUT> <HPC_DIR>
  python3 cs6_d21.py analyse <results_repo> <RUNS> <OUT_DIR>

NOTHING ALREADY COMMITTED IS RETYPED OR RE-RUN.
  cs6_transfer    head, H
  cs6_additivity  PBS
  cs13_qscale     FRAME_SHA, GEOM, SHIELD_REL, SHIELD_ABS
  cs6_d19         read_stab                      (as corrected by D19 NOTE 1)
Every pc.pc is COPIED from the committed run it pairs with, never regenerated, and every input
is gated to reduce to that run's committed input when the stability head is swapped back.
"""
import hashlib, re, sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import cs13_qscale as Q13
import cs6_transfer as TR
import cs6_additivity as AD
import cs6_d19 as D19

PREREG = "STAGE6_DIAGNOSTIC21_20261031.txt"
BASIS = "def2-SVPD"
GATE_EH = 1e-6
BUDGET = 1.067

# D21.2's table: name, the committed bare arm it pairs with, n, sum|q|, and GATE 2's targets.
SET = (
    ("S6-2_D17",    "d17/def2-SVPD/S6-2_bare_M",     2, 1.0777, -836.258538470, -836.241335755),
    ("S6-2_D16",    "d16/def2-SVPD/S6-2_bare_F1",    2, 1.3000, -836.262982927, -836.248538981),
    ("S6-3_D16",    "d16/def2-SVPD/S6-3_bare_F1",    3, 1.3000, -836.215800869, -836.200566235),
    ("S6-5_D16",    "d16/def2-SVPD/S6-5_bare_F1",    5, 1.3000, -836.272747122, -836.256686962),
    ("S6-5_D17",    "d17/def2-SVPD/S6-5_bare_M",     5, 1.6166, -836.281860904, -836.269391654),
    ("S6-10_D16",   "d16/def2-SVPD/S6-10_bare_F1",  10, 1.3000, -836.197018417, -836.177579450),
    ("S6-10_D15",   "transfer/def2-SVPD/lam050",    10, 2.1900, -836.168363751, -836.156155383),
    ("S6-20_D16F1", "d16/def2-SVPD/S6-20_bare_F1",  20, 1.3000, -836.236651180, -836.215399314),
    ("S6-20_D16F2", "d16/def2-SVPD/S6-20_bare_F2",  20, 2.1900, -836.232942143, -836.217532568),
)
CTL = {"R": -836.249865855128, "TS": -836.220651545262}
# D21.2 GATE 1: D20's committed control verdict, reused rather than repeated.
D20_RESULT = "d20/DIAGNOSTIC20_RESULT.txt"
D20_CTL_E0 = {"R": 0.05655489, "TS": 0.01882783}
# D20's measured stable ceiling, for READING 3.
LADDER_STABLE_TO = 2.88


def stop(m):
    print("STOP: " + m); sys.stdout.flush(); sys.exit(1)


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def stab_head():
    h = TR.head(BASIS)
    old = "%scf MaxIter 300 end\n"
    if old not in h:
        stop("cs6_transfer.head() no longer carries %r; the stability head cannot be derived "
             "from it without retyping, which this does not do." % old)
    return h.replace(old, "%scf\n  MaxIter 300\n  STABPerform true\nend\n")


def fsp_last(text):
    e = re.findall(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", text)
    if not e:
        stop("an output carries no FINAL SINGLE POINT ENERGY")
    return float(e[-1])


def gate1(T6):
    """D21.2 GATE 1: D20 measured this basis's control with this head and found it stable."""
    p = T6 / D20_RESULT
    if not p.is_file():
        stop("%s is not there; GATE 1 reuses DIAGNOSTIC 20's control measurement and cannot "
             "verify it without the committed result." % p)
    t = p.read_text(errors="replace")
    for st in Q13.GEOM:
        m = re.search(r"%s ctl_%-2s\s+(\w+)\s+E\(0\.\.2\)\s+([-+]\d+\.\d+)" % (BASIS, st), t)
        if not m:
            m = re.search(r"ctl_%s\s+(\w+)\s+E\(0\.\.2\)\s+([-+]\d+\.\d+)" % st, t)
        if not m:
            stop("could not find DIAGNOSTIC 20's %s control verdict for %s in %s"
                 % (BASIS, st, p))
        if m.group(1).upper() != "STABLE":
            stop("GATE 1: DIAGNOSTIC 20 reports the %s control %s as %s, not STABLE. The "
                 "acceptance test fails on the control, so nothing below is readable."
                 % (BASIS, st, m.group(1)))
        if abs(float(m.group(2)) - D20_CTL_E0[st]) > 5e-8:
            stop("GATE 1: DIAGNOSTIC 20's %s control %s has E(0) = %s, not the %.8f D21.2 states"
                 % (BASIS, st, m.group(2), D20_CTL_E0[st]))
    return True


def gate_targets(T6):
    """Every GATE 2 target is read from the committed file and checked against D21.2."""
    out = {}
    for nm, src, n, f, eR, eTS in SET:
        for st, want in (("R", eR), ("TS", eTS)):
            p = T6 / ("%s_%s" % (src, st)) / "job.out"
            if not p.is_file():
                stop("the committed source arm is missing: %s" % p)
            e = fsp_last(p.read_text(errors="replace"))
            if abs(e - want) > 5e-9:
                stop("the committed energy for %s %s is %.9f, not the %.9f D21.2 states"
                     % (nm, st, e, want))
            out[(nm, st)] = e
    return out


def gate_pc(T6, src, st, n):
    p = T6 / ("%s_%s" % (src, st)) / "pc.pc"
    if not p.is_file():
        stop("the committed charge file is missing: %s" % p)
    rows = [l.split() for l in p.read_text().strip().splitlines()]
    if len(rows[0]) == 1:
        if int(rows[0][0]) != n:
            stop("%s declares %s charges, D21.2 says %d" % (p, rows[0][0], n))
        rows = rows[1:]
    if len(rows) != n:
        stop("%s carries %d charges, D21.2 says %d" % (p, len(rows), n))
    return p


def gate_reduces(mine, committed):
    a = Path(mine).read_text(); b = Path(committed).read_text()
    h = stab_head()
    if not a.startswith(h):
        stop("%s does not begin with the stability head" % mine)
    if TR.head(BASIS) + a[len(h):] != b:
        stop("%s does not reduce to the committed input %s when its head is swapped back.\n"
             "   The only difference between them must be STABPerform." % (mine, committed))


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
    gate1(T6)
    gate_targets(T6)

    runs = []
    for nm, src, n, f, _, _ in SET:
        for st in Q13.GEOM:
            pc = gate_pc(T6, src, st, n)
            d = out / ("%s_%s" % (nm, st))
            d.mkdir(parents=True)
            (d / "substrate.xyz").write_text(geo[st])
            (d / "pc.pc").write_bytes(pc.read_bytes())
            (d / "job.inp").write_text(stab_head() + '%pointcharges "pc.pc"\n'
                                       '* xyzfile -2 1 substrate.xyz\n')
            gate_reduces(d / "job.inp", T6 / ("%s_%s" % (src, st)) / "job.inp")
            runs.append("%s_%s" % (nm, st))
    if len(runs) != 18:
        stop("%d runs were built, expected the 18 D21.2 states" % len(runs))

    pbs = AD.PBS % dict(tag="21", hpc=hpc, basis=".", n=len(runs),
                        runs="\n".join('run "%s"' % r for r in runs))
    pbs = pbs.replace("cm_d1421", "cm_d21") \
             .replace("DIAGNOSTIC 14 (cs6_additivity.py)", "DIAGNOSTIC 21 (cs6_d21.py)") \
             .replace("%d single points at .." % len(runs),
                      "%d single points with STABPerform, at %s." % (len(runs), BASIS)) \
             .replace("%s/./" % hpc, "%s/" % hpc).replace("W=%s/." % hpc, "W=%s" % hpc)
    if re.search(r"^#PBS([ \t].*)?[ \t]-q([ \t]|$)", pbs, re.M):
        stop("the PBS script requests a queue; the allocation is automatic")
    (out / "run.pbs").write_text(pbs)

    print("DIAGNOSTIC 21 prepared at %s" % out)
    print("   %d runs, %s, bare arm, the nine readable whole-design conditions" % (len(runs), BASIS))
    for nm, src, n, f, _, _ in SET:
        print("      %-13s n=%-3d sum|q| %.4f e   from %s" % (nm, n, f, src))
    print("   gate 1  DIAGNOSTIC 20's %s control is committed and STABLE" % BASIS)
    print("   gate 2  every target read from the committed files agrees with D21.2")
    print("   gate 3  every input reduces to its source's committed input: only STABPerform differs")
    return 0


def analyse(root, runs, outdir):
    T6 = Path(root) / "phase2.2" / "cs_stage6"
    runs = Path(runs)
    gate1(T6)
    T = gate_targets(T6)
    H = TR.H

    S = {}
    for nm, src, n, f, _, _ in SET:
        for st in Q13.GEOM:
            S[(nm, st)] = D19.read_stab(runs / ("%s_%s" % (nm, st)) / "job.out")

    L = []; w = L.append
    w("Stage 6, DIAGNOSTIC 21 RESULT.")
    w("Read by cs6_d21.py against STAGE6_DIAGNOSTIC21_20261031.txt. Criteria as pre-registered.")
    w("")
    w("GATE 1, the %s control, reused from DIAGNOSTIC 20 and not repeated" % BASIS)
    w("   ctl_R   STABLE  E(0) %+.8f      ctl_TS  STABLE  E(0) %+.8f"
      % (D20_CTL_E0["R"], D20_CTL_E0["TS"]))
    w("   GATE 1 passed: verified present and STABLE in the committed %s." % D20_RESULT)
    w("")
    w("GATE 2, each run's FIRST converged energy against its source's committed energy")
    w("   %-13s %-3s %19s %19s %11s" % ("design", "st", "first converged", "committed", "diff Eh"))
    g2 = {}
    for nm, src, n, f, _, _ in SET:
        for st in Q13.GEOM:
            an, tot, fsp = S[(nm, st)]
            disp = fsp[-1] - tot[-1]
            first = tot[0] + disp
            ok = abs(first - T[(nm, st)]) <= GATE_EH
            g2[(nm, st)] = (ok, first)
            w("   %-13s %-3s %19.9f %19.9f %11.2e%s"
              % (nm, st, first, T[(nm, st)], first - T[(nm, st)], "" if ok else "   FAILED"))
    nbad = sum(1 for ok, _ in g2.values() if not ok)
    if nbad:
        w("   GATE 2 FAILED on %d run(s): STABPerform changed their SCF, so they are not the" % nbad)
        w("   production solution and READING 2 is withheld for those designs (D21.2).")
    else:
        w("   GATE 2 passed: the first solution IS the production solution, everywhere.")
    w("")

    w("READING 1, are the readable designs stable?")
    w("   %-13s %-3s %8s %-9s %s" % ("design", "st", "sum|q|", "verdict", "E(0..2)"))
    unstable = []
    for nm, src, n, f, _, _ in SET:
        for st in Q13.GEOM:
            an, tot, fsp = S[(nm, st)]
            v, eig = an[0]
            if v == "UNSTABLE":
                unstable.append((nm, st, f))
            w("   %-13s %-3s %8.4f %-9s %s" % (nm, st, f, v,
              "  ".join("%+.8f" % x for x in eig[:3])))
            if len(an) > 1:
                w("      ORCA reconverged; second analysis %s, E(0) %+.8f, total energy "
                  "%+.4f kcal/mol" % (an[-1][0], an[-1][1][0], (tot[-1] - tot[0]) * H))
    v1 = "ALL STABLE" if not unstable else "SOME UNSTABLE"
    w("   VERDICT: %s" % v1)
    if unstable:
        w("   unstable: %s" % ", ".join("%s %s (%.4f e)" % u for u in unstable))
    w("")

    w("READING 2, what does any instability cost?")
    corr = {}
    for nm, src, n, f, _, _ in SET:
        if not any(S[(nm, st)][0][0][0] == "UNSTABLE" for st in Q13.GEOM):
            continue
        if not all(g2[(nm, st)][0] for st in Q13.GEOM):
            w("   %s: withheld, GATE 2 failed" % nm); continue
        anR, totR, fspR = S[(nm, "R")]; anT, totT, fspT = S[(nm, "TS")]
        dR, dT = fspR[-1] - totR[-1], fspT[-1] - totT[-1]
        first = ((totT[0] + dT) - (totR[0] + dR)) - (CTL["TS"] - CTL["R"])
        stab = ((totT[-1] + dT) - (totR[-1] + dR)) - (CTL["TS"] - CTL["R"])
        corr[nm] = (first * H, stab * H, (stab - first) * H, f)
    if not corr:
        w("   No design is unstable, so there is nothing to correct.")
        v2 = "NOTHING TO CORRECT"
    else:
        tol = max(Q13.SHIELD_REL * max(abs(c) for _, _, c, _ in corr.values()), Q13.SHIELD_ABS)
        w("   %-13s %8s %13s %13s %13s %9s" % ("design", "sum|q|", "ddE(first)",
                                               "ddE(stable)", "correction", "% budget"))
        for nm, (a, b, c, f) in sorted(corr.items()):
            w("   %-13s %8.4f %13.4f %13.4f %13.4f %9.1f"
              % (nm, f, a, b, c, 100.0 * abs(c) / BUDGET))
        w("   tol = max(%.2f x largest |correction|, %.2f) = %.4f kcal/mol   (one basis; "
          "nothing is pooled)" % (Q13.SHIELD_REL, Q13.SHIELD_ABS, tol))
        v2 = "CORRECTION MATERIAL" if any(abs(c) > tol for _, _, c, _ in corr.values()) \
             else "CORRECTION NEGLIGIBLE"
    w("   VERDICT: %s" % v2)
    over = [(nm, c, f) for nm, (_, _, c, f) in corr.items() if abs(c) > BUDGET]
    if over:
        w("")
        w("   D21.3's stated consequence applies. The Dset and the barrier these designs carry")
        w("   are PROVISIONAL until recomputed on a stable solution, and the readable window")
        w("   cannot be claimed on their numbers:")
        for nm, c, f in sorted(over):
            w("      %-13s sum|q| %.4f e   correction %+.4f kcal/mol" % (nm, f, c))
    w("")

    w("READING 3, does DIAGNOSTIC 20's onset transfer across arrangements?")
    w("   D20 measured S6-20's ladder stable to sum|q| = %.2f e. Every design here is at or"
      % LADDER_STABLE_TO)
    w("   below %.4f e, so the ladder predicts all eighteen runs stable."
      % max(f for _, _, _, f, _, _ in SET))
    below = [(nm, st, f) for nm, st, f in unstable if f < LADDER_STABLE_TO]
    if not unstable:
        v3 = "ONSET TRANSFERS"
        w("   Nine arrangements at four charge counts (n = 2, 3, 5, 10, 20) and three fields")
        w("   are all stable, as one arrangement's onset predicted.")
    else:
        v3 = "ONSET DOES NOT TRANSFER" if below else "ONSET TRANSFERS"
        if below:
            w("   These are unstable BELOW the ladder's measured stable ceiling:")
            for nm, st, f in below:
                w("      %-13s %-3s sum|q| %.4f e" % (nm, st, f))
            w("   The onset is therefore a property of the arrangement and not of the field")
            w("   alone. D20's boundary cannot be read as a field threshold.")
    w("   VERDICT: %s" % v3)
    w("")
    w("D21.4 stands whatever was returned. Stability is NECESSARY AND NOT SUFFICIENT: a stable")
    w("restricted solution is the lowest single-reference solution ORCA's rotations could find,")
    w("and it does not make the system single-reference. DIAGNOSTIC 19's READING 3 returned")
    w("DISAGREE on multireference character and that remains open at every field here. ALL")
    w("STABLE licenses no claim about it. ORCA's restart is a bounded heuristic, so any")
    w("correction above is a LOWER BOUND. The bare arm only: nothing here is a d_shield, which")
    w("is a difference of two arms.")

    outdir = Path(outdir); outdir.mkdir(parents=True, exist_ok=True)
    p = outdir / "DIAGNOSTIC21_RESULT.txt"
    p.write_text("\n".join(L) + "\n")
    print("\n".join(L)); print("\nwritten: %s" % p)
    return 1 if nbad else 0


if __name__ == "__main__":
    if len(sys.argv) == 5 and sys.argv[1] == "prepare":
        raise SystemExit(prepare(sys.argv[2], sys.argv[3], sys.argv[4]))
    if len(sys.argv) == 5 and sys.argv[1] == "analyse":
        raise SystemExit(analyse(sys.argv[2], sys.argv[3], sys.argv[4]))
    print(__doc__); raise SystemExit(2)
