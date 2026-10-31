#!/usr/bin/env python3
"""
cs6_d20.py - STAGE6_DIAGNOSTIC20_20261031.txt.

Where along the field does the production wavefunction stop being a minimum, and what does the
stable solution do to the barrier?

  python3 cs6_d20.py prepare <results_repo> <OUT> <HPC_DIR>
  python3 cs6_d20.py analyse <results_repo> <RUNS> <OUT_DIR>

NOTHING ALREADY COMMITTED IS RETYPED OR RE-RUN.
  cs6_transfer    design, head, read_xyz, H
  cs6_additivity  members
  cs13_qscale     FRAME_SHA, GEOM, read_xyz_raw, ECP_TAIL, SHIELD_REL, SHIELD_ABS
  cs6_d19         read_stab                              (as corrected by D19 NOTE 1)
Point-charge files are COPIED from the committed field tree, never regenerated.

THE INPUT GATES ARE EXACT. Every input this writes must reduce to the committed production
input, byte for byte, when its head is swapped back for cs6_transfer.head(). That is what makes
GATE 3's energy comparison meaningful: the only difference is STABPerform.
"""
import hashlib, re, sys, tarfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import cs13_qscale as Q13
import cs6_transfer as TR
import cs6_additivity as AD
import cs6_d19 as D19

PREREG = "STAGE6_DIAGNOSTIC20_20261031.txt"
NAME, NCH = "S6-20", 20
SVP, SVPD = "def2-SVP", "def2-SVPD"
TAG = {0.00: "lam000", 0.25: "lam025", 0.50: "lam050", 0.75: "lam075", 1.00: "lam100"}
GATE_EH = 1e-6

# D20.2: the runs, as the pre-registration tables them.
BARE = {SVP: (0.25, 0.50, 0.75), SVPD: (0.00, 0.25, 0.50, 0.75, 1.00)}
WALLED = (SVP, 1.00)

# D20.2 GATE 3's targets, as stated in the pre-registration. They are also READ from the
# committed files; a disagreement between the two is a stop, not an adjustment.
PROD = {
    (SVP, 0.25, "R"): -836.126173661, (SVP, 0.25, "TS"): -836.104876539,
    (SVP, 0.50, "R"): -836.120938885, (SVP, 0.50, "TS"): -836.108502797,
    (SVP, 0.75, "R"): -836.124785698, (SVP, 0.75, "TS"): -836.122550043,
    (SVPD, 0.25, "R"): -836.235775313, (SVPD, 0.25, "TS"): -836.215421621,
    (SVPD, 0.50, "R"): -836.233234395, (SVPD, 0.50, "TS"): -836.222558303,
    (SVPD, 0.75, "R"): -836.260454738, (SVPD, 0.75, "TS"): -836.247550257,
    (SVPD, 1.00, "R"): -836.346305931, (SVPD, 1.00, "TS"): -836.326879768,
}
CTL = {(SVP, "R"): -836.139795518709, (SVP, "TS"): -836.110376019855,
       (SVPD, "R"): -836.249865855128, (SVPD, "TS"): -836.220651545262}
WALL = {("wall", "R"): -835.801435357405, ("wall", "TS"): -835.811665184024,
        ("wall0", "R"): -836.138274953040, ("wall0", "TS"): -836.109033316175}


def stop(m):
    print("STOP: " + m); sys.stdout.flush(); sys.exit(1)


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def stab_head(basis):
    """The committed production head with the stability request added. Not retyped."""
    h = TR.head(basis)
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


def field_out(T6, basis, lam, st):
    a = T6 / "field" / "outputs" / "field_outputs.tar.xz"
    if not a.is_file():
        stop("the committed field archive is missing: %s" % a)
    t = tarfile.open(str(a))
    key = "./%s/%s_%s/job.out" % (basis, TAG[lam], st)
    names = {m.name: m for m in t.getmembers()}
    if key not in names:
        stop("%s does not contain %s" % (a, key))
    return t.extractfile(names[key]).read().decode(errors="replace")


def committed_targets(T6):
    """Read every GATE 3 target from the committed files and check it against the value the
    pre-registration states. Both must agree; neither is trusted alone."""
    got = {}
    for (basis, lam, st), want in sorted(PROD.items()):
        e = fsp_last(field_out(T6, basis, lam, st))
        if abs(e - want) > 5e-9:
            stop("the committed field energy for %s %s %s is %.9f, not the %.9f D20.2 states"
                 % (basis, TAG[lam], st, e, want))
        got[(basis, lam, st)] = e
    for (basis, st), want in sorted(CTL.items()):
        p = T6 / "additivity" / basis / ("ctl_%s" % st) / "job.out"
        if not p.is_file():
            stop("the committed control is missing: %s" % p)
        e = fsp_last(p.read_text(errors="replace"))
        if abs(e - want) > 5e-9:
            stop("the committed %s control %s is %.12f, not the %.12f D20.2 states"
                 % (basis, st, e, want))
        got[(basis, "ctl", st)] = e
    for (arm, st), want in sorted(WALL.items()):
        p = T6 / "d18" / SVP / ("%s_%s" % (NAME + "_" + arm if arm == "wall0" else NAME + "_wall", st))
        p = (T6 / "d18" / SVP / ("%s_%s_%s" % (NAME, arm, st))) if not p.is_dir() else p
        q = p / "job.out"
        if not q.is_file():
            stop("the committed d18 arm is missing: %s" % q)
        e = fsp_last(q.read_text(errors="replace"))
        if abs(e - want) > 5e-9:
            stop("the committed d18 %s %s is %.12f, not the %.12f D20.2 states"
                 % (arm, st, e, want))
        got[("d18", arm, st)] = e
    return got


def charges(T6, lam):
    d = TR.design(T6, NAME, NCH)
    return AD.members(d, d["site"], lam)


def gate_pc(T6, basis, lam, st):
    """GATE 2, compared at the precision the file stores (see cs6_d19 for why not bit equality)."""
    p = T6 / "field" / basis / ("%s_%s" % (TAG[lam], st)) / "pc.pc"
    if not p.is_file():
        stop("the committed charge file is missing: %s" % p)
    rows = [l.split() for l in p.read_text().strip().splitlines()]
    if len(rows[0]) == 1:
        rows = rows[1:]
    mem = charges(T6, lam)
    if len(rows) != NCH or len(mem) != NCH:
        stop("%s carries %d charges, the design has %d" % (p, len(rows), len(mem)))
    for (_, q, X), r in zip(mem, rows):
        if "%.6f" % q != "%.6f" % float(r[0]) or \
           any("%.10f" % a != "%.10f" % float(b) for a, b in zip(X, r[1:4])):
            stop("%s does not match the design at lambda %.2f" % (p, lam))
    return p


def write_bare(d, basis, geo_text, pc_src):
    d.mkdir(parents=True)
    (d / "substrate.xyz").write_text(geo_text)
    body = stab_head(basis)
    if pc_src is not None:
        (d / "pc.pc").write_bytes(Path(pc_src).read_bytes())
        body += '%pointcharges "pc.pc"\n'
    body += "* xyzfile -2 1 substrate.xyz\n"
    (d / "job.inp").write_text(body)


def write_walled(d, basis, geo_text, raw, mem, zero):
    """cs6_additivity.write_run's walled arm, with the stability head. The gate below checks the
    result reduces to d18's committed input byte for byte."""
    d.mkdir(parents=True)
    (d / "substrate.xyz").write_text(geo_text)
    body = ["* xyz -2 1"] + list(raw)
    for _, q, X in mem:
        body.append("  Ne>  %.6f  %18.8f %18.8f %18.8f  NewECP"
                    % (0.0 if zero else q, X[0], X[1], X[2]))
        body.append(Q13.ECP_TAIL)
    body.append("*")
    (d / "job.inp").write_text(stab_head(basis) + "\n".join(body) + "\n")


def gate_reduces(mine, committed, basis):
    """The input differs from the committed production input in the head and nowhere else."""
    a = Path(mine).read_text(); b = Path(committed).read_text()
    if not a.startswith(stab_head(basis)):
        stop("%s does not begin with the stability head" % mine)
    red = TR.head(basis) + a[len(stab_head(basis)):]
    if red != b:
        stop("%s does not reduce to the committed input %s when its head is swapped back.\n"
             "   The only difference between them must be STABPerform." % (mine, committed))


def rn(kind, lam, st):
    if kind == "wall0":
        return "%s_wall0_%s" % (NAME, st)
    if kind == "wall":
        return "%s_wall_%s_%s" % (NAME, TAG[lam], st)
    return "bare_%s_%s" % (TAG[lam], st)


def prepare(root, out, hpc):
    T6 = Path(root) / "phase2.2" / "cs_stage6"
    out = Path(out)
    if out.exists():
        stop("%s already exists; this never overwrites" % out)
    geo, raw = {}, {}
    for st, fn in (("R", "reactant_qm.xyz"), ("TS", "transition_state_qm.xyz")):
        p = T6 / "frame_42436" / fn
        if sha(p) != Q13.FRAME_SHA[fn]:
            stop("%s is not the committed frame" % p)
        geo[st], raw[st] = p.read_text(), Q13.read_xyz_raw(p)
    committed_targets(T6)

    runs = []
    for basis in (SVP, SVPD):
        for lam in BARE[basis]:
            for st in Q13.GEOM:
                pc = None if lam == 0.00 else gate_pc(T6, basis, lam, st)
                d = out / basis / rn("bare", lam, st)
                write_bare(d, basis, geo[st], pc)
                if lam > 0:
                    gate_reduces(d / "job.inp",
                                 T6 / "field" / basis / ("%s_%s" % (TAG[lam], st)) / "job.inp",
                                 basis)
                else:
                    gate_reduces(d / "job.inp",
                                 T6 / "additivity" / basis / ("ctl_%s" % st) / "job.inp", basis)
                runs.append("%s/%s" % (basis, rn("bare", lam, st)))
    wb, wl = WALLED
    mem = charges(T6, wl)
    for kind, zero in (("wall", False), ("wall0", True)):
        for st in Q13.GEOM:
            d = out / wb / rn(kind, wl, st)
            write_walled(d, wb, geo[st], raw[st], mem, zero)
            gate_reduces(d / "job.inp", T6 / "d18" / wb / ("%s_%s_%s" % (NAME, kind, st))
                         / "job.inp", wb)
            runs.append("%s/%s" % (wb, rn(kind, wl, st)))
    if len(runs) != 20:
        stop("%d runs were built, expected the 20 D20.2 states" % len(runs))

    pbs = AD.PBS % dict(tag="20", hpc=hpc, basis=".", n=len(runs),
                        runs="\n".join('run "%s"' % r for r in runs))
    pbs = pbs.replace("cm_d1420", "cm_d20") \
             .replace("DIAGNOSTIC 14 (cs6_additivity.py)", "DIAGNOSTIC 20 (cs6_d20.py)") \
             .replace("%d single points at .." % len(runs),
                      "%d single points with STABPerform, at %s and %s." % (len(runs), SVP, SVPD)) \
             .replace("%s/./" % hpc, "%s/" % hpc).replace("W=%s/." % hpc, "W=%s" % hpc)
    if re.search(r"^#PBS([ \t].*)?[ \t]-q([ \t]|$)", pbs, re.M):
        stop("the PBS script requests a queue; the allocation is automatic")
    (out / "run.pbs").write_text(pbs)

    print("DIAGNOSTIC 20 prepared at %s" % out)
    print("   %-10s bare lambda %s, R and TS" % (SVP, ", ".join("%.2f" % x for x in BARE[SVP])))
    print("   %-10s bare lambda %s, R and TS" % (SVPD, ", ".join("%.2f" % x for x in BARE[SVPD])))
    print("   %-10s walled and wall-at-zero at lambda %.2f, R and TS" % (wb, wl))
    print("   gate   every GATE 3 target read from the committed files agrees with D20.2")
    print("   gate   every input reduces to the committed production input when the head is")
    print("          swapped back: the only difference is STABPerform")
    print("   runs   %d" % len(runs))
    return 0


def analyse(root, runs, outdir):
    T6 = Path(root) / "phase2.2" / "cs_stage6"
    runs = Path(runs)
    T = committed_targets(T6)
    H = TR.H

    S = {}
    for basis in (SVP, SVPD):
        for lam in BARE[basis]:
            for st in Q13.GEOM:
                if lam > 0:
                    gate_pc(T6, basis, lam, st)
                S[(basis, "bare", lam, st)] = D19.read_stab(
                    runs / basis / rn("bare", lam, st) / "job.out")
    wb, wl = WALLED
    for kind in ("wall", "wall0"):
        for st in Q13.GEOM:
            S[(wb, kind, wl, st)] = D19.read_stab(runs / wb / rn(kind, wl, st) / "job.out")

    def first_energy(k):
        an, tot, fsp = S[k]
        return fsp[-1] if len(an) == 1 else None, tot

    L = []; w = L.append
    w("Stage 6, DIAGNOSTIC 20 RESULT.")
    w("Read by cs6_d20.py against STAGE6_DIAGNOSTIC20_20261031.txt. Criteria as pre-registered.")
    w("")

    # ---- GATE 3 -------------------------------------------------------------
    w("GATE 3, each run's FIRST converged energy against the committed production energy")
    w("   %-10s %-8s %-3s %18s %18s %11s" % ("basis", "arm", "st", "first converged",
                                             "committed", "diff Eh"))
    g3 = {}
    for basis in (SVP, SVPD):
        for lam in BARE[basis]:
            for st in Q13.GEOM:
                an, tot, fsp = S[(basis, "bare", lam, st)]
                disp = fsp[-1] - tot[-1]            # constant at fixed geometry
                first = tot[0] + disp
                want = T[(basis, lam, st)] if lam > 0 else T[(basis, "ctl", st)]
                ok = abs(first - want) <= GATE_EH
                g3[(basis, "bare", lam, st)] = ok
                w("   %-10s %-8s %-3s %18.9f %18.9f %11.2e%s"
                  % (basis, TAG[lam], st, first, want, first - want, "" if ok else "   FAILED"))
    for kind in ("wall", "wall0"):
        for st in Q13.GEOM:
            an, tot, fsp = S[(wb, kind, wl, st)]
            disp = fsp[-1] - tot[-1]
            first = tot[0] + disp
            want = T[("d18", kind, st)]
            ok = abs(first - want) <= GATE_EH
            g3[(wb, kind, wl, st)] = ok
            w("   %-10s %-8s %-3s %18.9f %18.9f %11.2e%s"
              % (wb, kind, st, first, want, first - want, "" if ok else "   FAILED"))
    nbad = sum(1 for v in g3.values() if not v)
    if nbad:
        w("   GATE 3 FAILED on %d run(s). STABPerform changed the SCF for those arms, so the" % nbad)
        w("   comparison READING 2 makes is meaningless for them and is withheld (D20.2).")
    else:
        w("   GATE 3 passed: the first solution IS the production solution, everywhere.")
    w("")

    # ---- GATE 1 -------------------------------------------------------------
    w("GATE 1, the def2-SVPD zero-field controls")
    g1 = True
    for st in Q13.GEOM:
        an, tot, fsp = S[(SVPD, "bare", 0.00, st)]
        v, eig = an[0]
        g1 = g1 and v == "STABLE"
        w("   %s ctl_%-2s  %-8s  E(0..2) %s" % (SVPD, st, v,
          "  ".join("%+.8f" % x for x in eig[:3])))
    if not g1:
        w("   GATE 1 FAILED: a charge-free closed-shell substrate reports an instability at")
        w("   %s, so the indicator is not measuring what it is meant to and nothing below" % SVPD)
        w("   is readable.")
        Path(outdir).mkdir(parents=True, exist_ok=True)
        (Path(outdir) / "DIAGNOSTIC20_RESULT.txt").write_text("\n".join(L) + "\n")
        print("\n".join(L)); return 1
    w("   GATE 1 passed.")
    w("")

    # ---- READING 1 ----------------------------------------------------------
    w("READING 1, where does the instability begin?   (the two bases are never pooled)")
    onset = {}
    for basis in (SVP, SVPD):
        w("   %s" % basis)
        seq = []
        for lam in BARE[basis]:
            row = []
            for st in Q13.GEOM:
                an, tot, fsp = S[(basis, "bare", lam, st)]
                v, eig = an[0]
                row.append((st, v, eig))
            unst = any(v == "UNSTABLE" for _, v, _ in row)
            seq.append((lam, unst))
            for st, v, eig in row:
                w("      lambda %.2f %-3s %-8s  E(0..2) %s" % (lam, st, v,
                  "  ".join("%+.8f" % x for x in eig[:3])))
        if basis == SVP:
            w("      lambda 0.00 and 1.00 at this basis were measured by DIAGNOSTIC 19:")
            w("         0.00  STABLE   (E(0) +0.05617271 R, +0.01728066 TS)")
            w("         1.00  UNSTABLE (E(0) -0.01083954 R, -0.03473872 TS)")
            seq = [(0.00, False)] + seq + [(1.00, True)]
        bad = [l for l, u in seq if u]
        if not bad:
            onset[basis] = "STABLE THROUGHOUT"
        else:
            first_bad = min(bad)
            tail = [u for l, u in seq if l >= first_bad]
            onset[basis] = ("ONSET AT lambda %.2f" % first_bad) if all(tail) else "NOT MONOTONIC"
            if all(tail) and first_bad == min(l for l, _ in seq if l > 0):
                onset[basis] = "UNSTABLE THROUGHOUT"
        w("      VERDICT: %s" % onset[basis])
    w("")

    # ---- READING 2 ----------------------------------------------------------
    w("READING 2, what does the stable solution do to the barrier?")
    corr = {}
    for basis in (SVP, SVPD):
        for lam in BARE[basis]:
            if lam == 0.00:
                continue
            kR, kT = (basis, "bare", lam, "R"), (basis, "bare", lam, "TS")
            if not (g3[kR] and g3[kT]):
                w("   %s lambda %.2f: withheld, GATE 3 failed" % (basis, lam)); continue
            anR, totR, fspR = S[kR]; anT, totT, fspT = S[kT]
            if len(anR) == 1 and len(anT) == 1:
                continue
            dR, dT = fspR[-1] - totR[-1], fspT[-1] - totT[-1]
            first = ((totT[0] + dT) - (totR[0] + dR)) - (T[(basis, "ctl", "TS")] - T[(basis, "ctl", "R")])
            stab = ((totT[-1] + dT) - (totR[-1] + dR)) - (T[(basis, "ctl", "TS")] - T[(basis, "ctl", "R")])
            corr[(basis, lam)] = (first * H, stab * H, (stab - first) * H)
    if not corr:
        w("   No run found an instability, so there is no correction to measure.")
        v2 = "CORRECTION NEGLIGIBLE"
    else:
        w("   %-10s %-7s %14s %14s %14s %9s" % ("basis", "lambda", "ddE(first)", "ddE(stable)",
                                                "correction", "% budget"))
        for (basis, lam), (a, b, c) in sorted(corr.items()):
            w("   %-10s %-7.2f %14.4f %14.4f %14.4f %9.1f"
              % (basis, lam, a, b, c, 100.0 * abs(c) / 1.067))
        tol = max(Q13.SHIELD_REL * max(abs(c) for _, _, c in corr.values()), Q13.SHIELD_ABS)
        w("   tol = max(%.2f x largest |correction|, %.2f) = %.4f kcal/mol"
          % (Q13.SHIELD_REL, Q13.SHIELD_ABS, tol))
        v2 = "CORRECTION MATERIAL" if any(abs(c) > tol for _, _, c in corr.values()) \
             else "CORRECTION NEGLIGIBLE"
    w("   VERDICT: %s" % v2)
    over = [(b, l, c) for (b, l), (_, _, c) in corr.items() if abs(c) > 1.067]
    if over:
        w("")
        w("   D20.3's stated consequence applies. These corrections exceed DIAGNOSTIC 11.3's")
        w("   budget of 1.067 kcal/mol, so every committed barrier, Dset and d_shield at or")
        w("   above that lambda at that basis is PROVISIONAL until recomputed on a stable")
        w("   solution, and must be quoted as provisional or not quoted. DIAGNOSTIC 18's")
        w("   ddE(bare) = -24.6898 kcal/mol is among them.")
        for b, l, c in sorted(over):
            w("      %-10s lambda %.2f   correction %+.4f kcal/mol" % (b, l, c))
    w("")

    # ---- READING 3 ----------------------------------------------------------
    w("READING 3, does the Pauli wall change stability?   (%s, lambda %.2f)" % (wb, wl))
    bare_unstable = True      # DIAGNOSTIC 19, both states
    w("   bare   UNSTABLE in both states (DIAGNOSTIC 19, not re-run here)")
    wv = []
    for kind in ("wall", "wall0"):
        for st in Q13.GEOM:
            an, tot, fsp = S[(wb, kind, wl, st)]
            v, eig = an[0]
            if kind == "wall":
                wv.append(v)
            w("   %-6s %-3s %-8s  E(0..2) %s" % (kind, st, v,
              "  ".join("%+.8f" % x for x in eig[:3])))
    wall_unstable = any(v == "UNSTABLE" for v in wv)
    if wall_unstable == bare_unstable:
        v3 = "WALL DOES NOT CHANGE STABILITY"
    elif bare_unstable and not wall_unstable:
        v3 = "WALL STABILISES"
    else:
        v3 = "WALL DESTABILISES"
    w("   VERDICT: %s" % v3)
    w("")
    w("D20.4 stands whatever was returned. The reconverged solution is another RESTRICTED")
    w("solution reached by rotating orbitals: it is not a multireference treatment and it does")
    w("not make the system single-reference. DIAGNOSTIC 19's READING 3 returned DISAGREE and")
    w("that stands. ORCA's restart is a bounded heuristic, so every correction above is a LOWER")
    w("BOUND on the correction, not the correction. The numbers here are for S6-20 at these")
    w("lambdas and do not transfer to DIAGNOSTIC 16's or 17's designs without measuring them.")

    outdir = Path(outdir); outdir.mkdir(parents=True, exist_ok=True)
    p = outdir / "DIAGNOSTIC20_RESULT.txt"
    p.write_text("\n".join(L) + "\n")
    print("\n".join(L)); print("\nwritten: %s" % p)
    return 1 if nbad else 0


if __name__ == "__main__":
    if len(sys.argv) == 5 and sys.argv[1] == "prepare":
        raise SystemExit(prepare(sys.argv[2], sys.argv[3], sys.argv[4]))
    if len(sys.argv) == 5 and sys.argv[1] == "analyse":
        raise SystemExit(analyse(sys.argv[2], sys.argv[3], sys.argv[4]))
    print(__doc__); raise SystemExit(2)
