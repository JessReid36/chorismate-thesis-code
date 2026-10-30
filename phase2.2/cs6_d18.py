#!/usr/bin/env python3
"""
cs6_d18.py - STAGE6_DIAGNOSTIC18_20261030.txt (experiment E3 of ROUTE_ASSESSMENT_20261028).

Does a Pauli wall open the HOMO-LUMO gap of S6-20 at full strength at def2-SVP?

  python3 cs6_d18.py prepare <results_repo> <OUT> <HPC_DIR>
  python3 cs6_d18.py analyse <results_repo> <RUNS> <OUT_DIR>

NOTHING ALREADY COMMITTED IS RETYPED OR RE-RUN.
  cs6_additivity  write_run, members, PBS, one0
  cs6_negatives   zeroed, com, corrected_dip, check_correction        (AMENDMENT 1 A1.2)
  cs6_transfer    design, head, controls, verdict, one, H, read_xyz
  cs13_qscale     FRAME_SHA, GEOM, read_xyz_raw, SHIELD_REL, SHIELD_ABS
  cs6_audit       D_LW, R_DIP, GAP_ABS, R_GAP, R_CYC                  (D6.3's constants)

THE BARE ARM IS REUSED, NOT RE-RUN (D18.2). It is field/def2-SVP/lam100_{R,TS}, which lives
inside field/outputs/field_outputs.tar.xz. analyse() extracts it and puts it through the SAME
committed validator, cs6_transfer.one(), that a fresh run would go through.
"""
import hashlib, re, shutil, sys, tarfile, tempfile
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import cs13_qscale as Q13
import cs6_transfer as TR
import cs6_additivity as AD
import cs6_negatives as NG
import cs6_audit as AUD

PREREG = "STAGE6_DIAGNOSTIC18_20261030.txt"
NAME, NCH = "S6-20", 20
BASIS = "def2-SVP"
LAM = 1.0                     # D18.2: the value at which a committed bare arm exists
CTL_EH = 1e-6                 # D18.2 GATE

# D18.2's gate targets: the committed def2-SVP control energies, same geometry, same method line.
CTL_REF = {"R": -836.139795518709, "TS": -836.110376019855}
CTL_SRC = "additivity/def2-SVP/ctl_%s"

# D18.2: the reused bare arm, and the archive it lives in.
BARE_ARC = "field/outputs/field_outputs.tar.xz"
BARE_MEM = "./def2-SVP/lam100_%s/job.out"
BARE_PC = "field/def2-SVP/lam100_%s/pc.pc"

# D18.1, measured and stated before any run of this diagnostic exists.
BARE_GAP = {"R": 0.6099, "TS": 0.9236}
# D18.2, the committed charge set this must reproduce exactly.
SUMABS, MAXQ, SUMQ2 = 5.7600, 0.5500, 2.0476


def stop(m):
    print("STOP: " + m); sys.stdout.flush(); sys.exit(1)


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def charges(T6):
    """S6-20 at lambda = 1.00, from the committed design. Not retyped."""
    d = TR.design(T6, NAME, NCH)
    mem = AD.members(d, d["site"], LAM)
    if len(mem) != NCH:
        stop("the design has %d charges, expected %d" % (len(mem), NCH))
    sa = sum(abs(q) for _, q, _ in mem)
    mq = max(abs(q) for _, q, _ in mem)
    s2 = sum(q * q for _, q, _ in mem)
    for got, want, nm, tol in ((sa, SUMABS, "sum|q|", 5e-5), (mq, MAXQ, "max|q|", 5e-5),
                               (s2, SUMQ2, "sum q^2", 5e-5)):
        if abs(got - want) > tol:
            stop("%s is %.6f, not the %.4f D18.2 states" % (nm, got, want))
    return mem


def gate_pairs_with_bare(T6, mem):
    """D18.2's load-bearing gate. The bare arm is reused, so the charges these runs carry must
    be the committed bare arm's charges EXACTLY - not to a tolerance, exactly."""
    for st in Q13.GEOM:
        p = T6 / (BARE_PC % st)
        if not p.is_file():
            stop("the reused bare arm's charge file is missing: %s" % p)
        rows = [l.split() for l in p.read_text().strip().splitlines()]
        if len(rows[0]) == 1:
            if int(rows[0][0]) != NCH:
                stop("%s declares %s charges, expected %d" % (p, rows[0][0], NCH))
            rows = rows[1:]
        if len(rows) != NCH:
            stop("%s carries %d charges, expected %d" % (p, len(rows), NCH))
        for (_, q, X), r in zip(mem, rows):
            if float(r[0]) != q:
                stop("%s: charge %r is not the design's %r. The walled arm would not pair "
                     "with the reused bare arm." % (p, float(r[0]), q))
            for a, b in zip(X, r[1:4]):
                if float(b) != float(a):
                    stop("%s: a coordinate differs from the design's. The walled arm would not "
                         "pair with the reused bare arm." % p)


def gate_same_head(T6):
    """D18.2: the generated method line must be the reused bare arm's, byte for byte."""
    mine = TR.head(BASIS)
    for st in Q13.GEOM:
        p = T6 / ("field/def2-SVP/lam100_%s/job.inp" % st)
        theirs = p.read_text().split("%pointcharges")[0]
        if theirs != mine:
            stop("the method line this generator writes is not the reused bare arm's.\n"
                 "  generated: %r\n  committed: %r" % (mine, theirs))


def gate_controls(T6):
    """D18.2 GATE: the committed control energies these runs must reproduce."""
    out = {}
    for st in Q13.GEOM:
        p = T6 / (CTL_SRC % st) / "job.out"
        if not p.is_file():
            stop("the committed def2-SVP control is missing: %s" % p)
        e = re.findall(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", p.read_text(errors="replace"))
        if len(e) != 1:
            stop("%s: %d final energies, expected 1" % (p, len(e)))
        if abs(float(e[0]) - CTL_REF[st]) > 1e-9:
            stop("the committed control %s is %.12f, not the %.12f D18.2 states"
                 % (st, float(e[0]), CTL_REF[st]))
        out[st] = float(e[0])
    return out


def rn(arm, st):
    return "ctl_%s" % st if arm == "none" else "%s_%s_%s" % (NAME, arm, st)


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
    NG.check_correction(T6)

    mem = charges(T6)
    gate_pairs_with_bare(T6, mem)
    gate_same_head(T6)
    ctl = gate_controls(T6)

    runs = []
    for st in Q13.GEOM:
        AD.write_run(out / BASIS / rn("none", st), BASIS, geo[st], raw[st], [], "none")
        runs.append(rn("none", st))
    for arm in ("wall0", "wall"):
        for st in Q13.GEOM:
            AD.write_run(out / BASIS / rn(arm, st), BASIS, geo[st], raw[st], mem, arm)
            runs.append(rn(arm, st))
    if len(runs) != 6:
        stop("%d runs were built, expected the 6 D18.2 states" % len(runs))

    pbs = AD.PBS % dict(tag="18", hpc=hpc, basis=BASIS, n=len(runs),
                        runs="\n".join('run "%s"' % r for r in runs))
    pbs = pbs.replace("cm_d1418", "cm_d18").replace("DIAGNOSTIC 14 (cs6_additivity.py)",
                                                    "DIAGNOSTIC 18 (cs6_d18.py)")
    if re.search(r"^#PBS([ \t].*)?[ \t]-q([ \t]|$)", pbs, re.M):
        stop("the PBS script requests a queue; the allocation is automatic")
    (out / BASIS / "run.pbs").write_text(pbs)

    print("DIAGNOSTIC 18 prepared at %s" % out)
    print("   design      %s, n=%d, lambda=%.6f at %s" % (NAME, NCH, LAM, BASIS))
    print("   charges     sum|q| %.4f e   max|q| %.4f e   sum q^2 %.4f e^2  (2.4%% over the "
          "2.0 budget, declared in D18.2)" % (SUMABS, MAXQ, SUMQ2))
    print("   gate        charges are byte-exact against the reused bare arm, both states")
    print("   gate        method line is byte-exact against the reused bare arm")
    print("   gate        committed controls %.12f / %.12f Eh" % (ctl["R"], ctl["TS"]))
    print("   runs        %s" % ", ".join(runs))
    print("   bare arm    REUSED from %s, not re-run" % (BARE_ARC,))
    return 0


def bare_runs(T6, tmp):
    """Extract the reused bare arm so the committed validator can be applied to it unchanged."""
    arc = T6 / BARE_ARC
    if not arc.is_file():
        stop("the reused bare arm's archive is missing: %s" % arc)
    t = tarfile.open(str(arc))
    names = {m.name: m for m in t.getmembers()}
    out = {}
    for st in Q13.GEOM:
        key = BARE_MEM % st
        if key not in names:
            stop("%s does not contain %s" % (arc, key))
        d = Path(tmp) / ("bare_%s" % st)
        d.mkdir(parents=True)
        (d / "job.out").write_bytes(t.extractfile(names[key]).read())
        out[st] = d
    return out


def analyse(root, runs, outdir):
    T6 = Path(root) / "phase2.2" / "cs_stage6"
    runs = Path(runs)
    C = TR.controls(T6)
    GEO = {st: TR.read_xyz(T6 / "frame_42436" / fn) for st, fn in
           (("R", "reactant_qm.xyz"), ("TS", "transition_state_qm.xyz"))}
    O = {st: NG.com(T6, st) for st in Q13.GEOM}
    mem = charges(T6)
    gate_pairs_with_bare(T6, mem)
    zer = NG.zeroed(mem)

    E, D, CLS, WHY = {}, {}, {}, {}

    # the six new runs
    for st in Q13.GEOM:
        d = runs / BASIS / rn("none", st)
        e, nb, dd = AD.one0(d, GEO[st], 0, False)
        E[("ctl", st)], D[("ctl", st)] = e, dd
        CLS[("ctl", st)], WHY[("ctl", st)] = TR.verdict(dd, C)
        if abs(e - CTL_REF[st]) > CTL_EH:
            stop("GATE: control %s is %.12f Eh, %.2e from the committed %.12f. D18.2 requires "
                 "1e-6." % (st, e, abs(e - CTL_REF[st]), CTL_REF[st]))
    for arm, mm in (("wall0", zer), ("wall", mem)):
        for st in Q13.GEOM:
            d = runs / BASIS / rn(arm, st)
            e, nb, dd = AD.one0(d, GEO[st], len(mm), True)
            dd = dict(dd)
            dd["dip"] = NG.corrected_dip(d, mm, O[st])      # A1.2, as D14-D17 apply it
            E[(arm, st)], D[(arm, st)] = e, dd
            CLS[(arm, st)], WHY[(arm, st)] = TR.verdict(dd, C)

    # the reused bare arm, through the committed validator
    tmp = tempfile.mkdtemp(prefix="d18bare_")
    try:
        bd = bare_runs(T6, tmp)
        for st in Q13.GEOM:
            e, nb, dd = TR.one(bd[st], GEO[st], NCH)
            E[("bare", st)], D[("bare", st)] = e, dd
            CLS[("bare", st)], WHY[("bare", st)] = TR.verdict(dd, C)
            if abs(dd["gap"] - BARE_GAP[st]) > 5e-4:
                stop("the reused bare arm's %s gap reads %.4f eV, not the %.4f D18.1 states"
                     % (st, dd["gap"], BARE_GAP[st]))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    H = TR.H
    ddE_bare = ((E[("bare", "TS")] - E[("bare", "R")])
                - (E[("ctl", "TS")] - E[("ctl", "R")])) * H
    ddE_wall = ((E[("wall", "TS")] - E[("wall", "R")])
                - (E[("wall0", "TS")] - E[("wall0", "R")])) * H
    dset = ddE_bare - ddE_wall

    L = []
    w = L.append
    w("Stage 6, DIAGNOSTIC 18 RESULT (experiment E3).")
    w("Read by cs6_d18.py against STAGE6_DIAGNOSTIC18_20261030.txt. Criteria as pre-registered.")
    w("")
    w("D18.2  THE RUNS")
    w("   %s at %s, lambda = %.6f;  sum|q| %.4f e, max|q| %.4f e, sum q^2 %.4f e^2"
      % (NAME, BASIS, LAM, SUMABS, MAXQ, SUMQ2))
    w("   GATE  controls reproduce the committed def2-SVP energies to 1e-6 Eh:")
    for st in Q13.GEOM:
        w("         ctl_%-2s  %.12f Eh   committed %.12f   diff %.2e"
          % (st, E[("ctl", st)], CTL_REF[st], abs(E[("ctl", st)] - CTL_REF[st])))
    w("   The bare arm is REUSED from %s and was not re-run." % BARE_ARC)
    w("")
    w("   %-16s %-10s %9s %6s %9s %9s" % ("arm", "D6.3", "gap eV", "cyc", "maxlw", "dip"))
    for arm in ("ctl", "bare", "wall0", "wall"):
        for st in Q13.GEOM:
            d = D[(arm, st)]
            w("   %-16s %-10s %9.4f %6d %9.4f %9.3f"
              % ("%s_%s" % (arm, st), CLS[(arm, st)], d["gap"], d["cyc"], d["maxlw"], d["dip"]))
            if WHY[(arm, st)]:
                w("      %s" % "; ".join(WHY[(arm, st)]))
    w("")
    w("READING 1, does the wall open the gap?")
    gaps = {st: D[("wall", st)]["gap"] for st in Q13.GEOM}
    opened = all(gaps[st] >= AUD.GAP_ABS for st in Q13.GEOM)
    v1 = "WALL OPENS THE GAP" if opened else "WALL DOES NOT OPEN THE GAP"
    for st in Q13.GEOM:
        w("   %-3s walled %.4f eV   bare %.4f eV   GAP_ABS %.2f eV"
          % (st, gaps[st], D[("bare", st)]["gap"], AUD.GAP_ABS))
    w("   VERDICT: %s" % v1)
    w("")
    w("READING 2, d_shield for a whole design at full strength at %s." % BASIS)
    tol = max(Q13.SHIELD_REL * abs(dset), Q13.SHIELD_ABS)
    w("   ddE(bare) %+.4f   ddE(wall) %+.4f   Dset %+.4f kcal/mol   tol %.4f"
      % (ddE_bare, ddE_wall, dset, tol))
    v2 = "MEASURABLE" if abs(dset) > tol else "NO ARTEFACT TO MEASURE"
    w("   VERDICT: %s" % v2)
    if abs(dset) > tol:
        w("   %.1f%% of DIAGNOSTIC 11.3's 1.067 kcal/mol budget" % (100.0 * abs(dset) / 1.067))
    unread = [k for k in (("bare", "R"), ("bare", "TS"), ("wall", "R"), ("wall", "TS"))
              if CLS[k] not in ("CLEAN", "POLARISED")]
    if unread:
        w("   NOT READABLE under A1.3: %s. Dset is reported and no claim rests on it (D18.3)."
          % ", ".join("%s_%s (%s)" % (a, s, CLS[(a, s)]) for a, s in unread))
    else:
        w("   Both arms are CLEAN or POLARISED, so Dset is readable under A1.3.")
    w("")
    w("READING 3, does the density move at %s?" % BASIS)
    w("   %-12s %10s %10s %10s %10s" % ("arm", "maxlw", "excess", "dip", "dip/ctl"))
    trip = {}
    for arm in ("bare", "wall"):
        t3 = False
        for st in Q13.GEOM:
            d = D[(arm, st)]; c = C[d["nb"]]
            ex = d["maxlw"] - c["maxlw"]; dr = d["dip"] / c["dip"]
            t3 = t3 or ex >= AUD.D_LW or dr >= AUD.R_DIP
            w("   %-12s %10.4f %10.4f %10.3f %10.3f"
              % ("%s_%s" % (arm, st), d["maxlw"], ex, d["dip"], dr))
        trip[arm] = t3
    w("   D_LW %.2f   R_DIP %.1f" % (AUD.D_LW, AUD.R_DIP))
    if not trip["bare"] and not trip["wall"]:
        v3 = "NEITHER ARM MOVES"
    elif trip["bare"] and not trip["wall"]:
        v3 = "THE WALL REMOVES IT"
    elif trip["wall"] and not trip["bare"]:
        v3 = "THE WALL ADDS IT"
    else:
        v3 = "BOTH ARMS MOVE"
    w("   VERDICT: %s" % v3)
    w("")
    w("D18.4's limitation stands whatever Reading 1 returned: a gap that stays closed with the")
    w("wall shows the near-degeneracy is not produced by the representation. It does not show")
    w("it is physical. At this gap single-reference B3LYP is outside its range and ddE(bare)")
    w("= %+.4f kcal/mol is not a number the method is entitled to report." % ddE_bare)

    outdir = Path(outdir); outdir.mkdir(parents=True, exist_ok=True)
    p = outdir / "DIAGNOSTIC18_RESULT.txt"
    p.write_text("\n".join(L) + "\n")
    print("\n".join(L))
    print("\nwritten: %s" % p)
    return 0


if __name__ == "__main__":
    if len(sys.argv) == 5 and sys.argv[1] == "prepare":
        raise SystemExit(prepare(sys.argv[2], sys.argv[3], sys.argv[4]))
    if len(sys.argv) == 5 and sys.argv[1] == "analyse":
        raise SystemExit(analyse(sys.argv[2], sys.argv[3], sys.argv[4]))
    print(__doc__); raise SystemExit(2)
