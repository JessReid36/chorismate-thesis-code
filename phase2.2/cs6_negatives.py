#!/usr/bin/env python3
"""
cs6_negatives.py - STAGE6_DIAGNOSTIC14_AMENDMENT1_20261028.txt, A1.4 to A1.6.

Why are S6-10's six negative charges inert behind the Pauli wall at def2-SVP and not at
def2-SVPD?  Each of the six alone, at arm B, bare and walled, at both bases.

  python3 cs6_negatives.py prepare <results_repo> <OUT> <HPC_DIR>
  python3 cs6_negatives.py analyse <results_repo> <RUNS> <OUT_DIR>

NOTHING ALREADY COMMITTED IS RETYPED.
  cs6_additivity  write_run(), members(), rname(), one0(), the PBS template, ARM_LAM
  cs6_transfer    design(), lambdas(), head(), controls(), verdict(), one()
  cs13_qscale     ECP_TAIL, FRAME_SHA, GEOM, read_xyz_raw, descriptors, SHIELD_REL, SHIELD_ABS
Only what Amendment 1 adds is here: the dipole correction of A1.2, and A1.5's three readings.
"""
import hashlib, re, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import cs13_qscale as Q13
import cs6_transfer as TR
import cs6_additivity as AD

AMEND = "STAGE6_DIAGNOSTIC14_AMENDMENT1_20261028.txt"
AMEND2 = "STAGE6_DIAGNOSTIC14_AMENDMENT2_20261028.txt"
DESIGN, NSPH = "S6-10", 10
ARM = "B"                                   # A1.4: arm B only, the larger of the two effects
BASES = ("def2-SVPD", "def2-SVP")           # A1.4: both, because the contrast IS the question
DEBYE = 2.541746473                         # a.u. -> Debye, ORCA's own factor
BOHR = Q13.BOHR
# A1.2's origin: the mass-weighted centroid of the 24 substrate atoms, no Ne> contribution.
# These masses reproduce ORCA's printed nuclear dipole exactly on DIAGNOSTIC 14's committed
# outputs; prepare() checks that before it writes anything.
MASS = {"H": 1.008, "C": 12.011, "O": 15.999, "N": 14.007}
CTL_EH = 1e-6                               # A1.4's gate against D14's committed controls
N6_TOL = 1e-4                               # and against its committed n6neg d_shield


def stop(m):
    print("STOP: " + m); sys.stdout.flush(); sys.exit(1)


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def zeroed(mem):
    """The wall-at-zero arm carries NO charge on its centres - cs6_additivity.write_run forces
    0.0 there whatever is passed. A1.2's correction must therefore subtract nothing for it, so
    the member list it is given has to say zero too. Passing the lambda-1 charges here subtracts
    a dipole that is not in the calculation; the rehearsal caught that."""
    return [(s, 0.0, X) for s, _, X in mem]


def negatives(D):
    s = [t for i, t in enumerate(D["site"]) if D["q"][i] < 0]
    if len(s) != 6:
        stop("expected six negative sites in %s, parsed %d" % (DESIGN, len(s)))
    return s


def com(T6, st):
    fn = "reactant_qm.xyz" if st == "R" else "transition_state_qm.xyz"
    L = [l.split() for l in (T6 / "frame_42436" / fn).read_text().splitlines()[2:26]]
    X = np.array([[float(v) for v in f[1:4]] for f in L])
    m = np.array([MASS[f[0]] for f in L])
    return (m[:, None] * X).sum(0) / m.sum()


def qdip(mem, O):
    """A1.2: the classical dipole of the charges about O, in atomic units."""
    return sum(q * (X - O) for _, q, X in mem) / BOHR


def vecs(d):
    t = (Path(d) / "job.out").read_text(errors="replace")
    out = []
    for k in ("Electronic contribution", "Nuclear contribution"):
        m = re.search(k + r"\s*:\s*(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)", t)
        if not m:
            stop("%s: no %s in the dipole block" % (d, k))
        out.append(np.array([float(x) for x in m.groups()]))
    return out


def corrected_dip(d, mem, O):
    """The printed total dipole with the charges' own classical term removed (A1.2)."""
    e, n = vecs(d)
    return float(np.linalg.norm(e + n - qdip(mem, O))) * DEBYE


# ---------------------------------------------------------------- prepare

def check_correction(T6):
    """A1.2 must reproduce DIAGNOSTIC 14's committed numbers before it is used on anything new."""
    A = T6 / "additivity"
    if not (A / "DIAGNOSTIC14_RESULT.txt").is_file():
        stop("DIAGNOSTIC 14's result is not committed; the correction cannot be checked")
    D = TR.design(T6, DESIGN, NSPH)
    REF = TR.design(T6, TR.REF_DESIGN, TR.REF_N)
    lam = TR.lambdas(D, REF)[AD.ARM_LAM["A"]]
    worst = 0.0
    for st in Q13.GEOM:
        mem = AD.members(D, negatives(D), lam)
        O = com(T6, st)
        obs = vecs(A / "def2-SVPD" / ("n6neg_wall_A_%s" % st))[1] \
            - vecs(A / "def2-SVPD" / ("ctl_%s" % st))[1]
        worst = max(worst, float(np.abs(obs - qdip(mem, O)).max()))
    if worst > 1e-6:
        stop("A1.2's correction does not reproduce DIAGNOSTIC 14's printed nuclear dipole "
             "(worst component off by %.2e a.u.). The origin or the masses are wrong; nothing "
             "is written." % worst)
    print("   A1.2's correction reproduces D14's nuclear dipole to %.1e a.u. (both states)" % worst)


def prepare(root, out, hpc):
    T6 = Path(root) / "phase2.2" / "cs_stage6"
    out = Path(out)
    if out.exists():
        stop("%s already exists - this never overwrites" % out)
    for f in (AMEND, AMEND2):
        if not (T6 / f).is_file():
            stop("%s is not committed. The amendments go in before any input exists." % f)
    print("   %s is committed (sha256 %s)" % (AMEND, sha(T6 / AMEND)[:16]))
    geo, raw = {}, {}
    for st, fn in (("R", "reactant_qm.xyz"), ("TS", "transition_state_qm.xyz")):
        p = T6 / "frame_42436" / fn
        if sha(p) != Q13.FRAME_SHA[fn]:
            stop("%s does not match the sha256 cs13_qscale records" % fn)
        geo[st], raw[st] = p.read_text(), Q13.read_xyz_raw(p)
    print("   both geometries match the sha256 the committed scripts record")
    check_correction(T6)

    D = TR.design(T6, DESIGN, NSPH)
    REF = TR.design(T6, TR.REF_DESIGN, TR.REF_N)
    lam = TR.lambdas(D, REF)[AD.ARM_LAM[ARM]]
    neg = negatives(D)
    i = {s: k for k, s in enumerate(D["site"])}
    print("   arm %s, lambda %.6f; the six negatives of %s:" % (ARM, lam, DESIGN))
    for s in neg:
        dn = Q13.descriptors(D["X"][i[s]], Q13.read_xyz(T6 / "frame_42436" / "reactant_qm.xyz"))["dn"]
        print("      site %-5d q %+.4f -> %+.4f e   nearest atom %.3f A" % (s, D["q"][i[s]],
                                                                            lam * D["q"][i[s]], dn))
    C = TR.controls(T6)
    for nb in (264, 402):
        if nb not in C:
            stop("D6's control table has no nb %d row" % nb)

    out.mkdir(parents=True)
    plan = []
    for basis in BASES:
        rs = []
        for s in neg:
            nm = "s%d" % s
            mem = AD.members(D, [s], lam)
            mem0 = AD.members(D, [s], 1.0)
            for st in Q13.GEOM:
                for arm, m in (("wall0", mem0), ("bare", mem), ("wall", mem)):
                    rn = "%s_%s" % (AD.rname(nm, arm, ARM), st)
                    AD.write_run(out / basis / rn, basis, geo[st], raw[st], m, arm)
                    rs.append(rn)
                    plan.append((basis, nm, arm, st, abs(m[0][1]) if arm != "wall0" else 0.0))
        for st in Q13.GEOM:
            rn = "ctl_%s" % st
            AD.write_run(out / basis / rn, basis, geo[st], raw[st], [], "none")
            rs.append(rn); plan.append((basis, "ctl", "none", st, 0.0))
        (out / basis / "run.pbs").write_text(
            AD.PBS % dict(tag=("n" if basis.endswith("D") else "s"), hpc=hpc, basis=basis,
                          n=len(rs), runs="\n".join("run %s" % r for r in rs)))
    with open(out / "PLAN.tsv", "w") as fh:
        fh.write("# %s A1.4; cs6_negatives.py prepare\n" % AMEND)
        fh.write("# %s negatives at arm %s, lambda %.6f\n" % (DESIGN, ARM, lam))
        fh.write("basis\tsite\tarm\tstate\tabs_q\n")
        for r in plan:
            fh.write("%s\t%s\t%s\t%s\t%.4f\n" % r)
    files = sorted(p for p in out.rglob("*") if p.is_file())
    (out / "INPUT_CHECKSUMS.sha256").write_text(
        "".join("%s  %s\n" % (sha(p), p.relative_to(out)) for p in files))
    print("   %d run directories, %d files including the manifest" % (len(plan), len(files) + 1))
    for basis in BASES:
        print("      %-10s %3d runs" % (basis, sum(1 for r in plan if r[0] == basis)))
    print("   NOTHING WAS SUBMITTED BY THIS BUNDLE.")


# ---------------------------------------------------------------- analyse

def d14_n6neg(T6):
    """The n6neg d_shield DIAGNOSTIC 14 printed, read out of its own committed result."""
    out = {}
    for l in (T6 / "additivity" / "DIAGNOSTIC14_RESULT.txt").read_text().splitlines():
        m = re.match(r"\s*(def2-\S+)\s+arm\s+(\S)\s+d_shield\s+([+-]\d+\.\d+)", l)
        if m:
            out[(m.group(1), m.group(2))] = float(m.group(3))
    if not out:
        stop("could not read DIAGNOSTIC 14's n6neg d_shield out of its committed result")
    return out


def analyse(root, runs, outdir):
    T6, R, outdir = Path(root) / "phase2.2" / "cs_stage6", Path(runs), Path(outdir)
    for f in (AMEND, AMEND2):
        if not (T6 / f).is_file():
            stop("%s is not committed" % f)
    A = T6 / "additivity"
    D = TR.design(T6, DESIGN, NSPH); REF = TR.design(T6, TR.REF_DESIGN, TR.REF_N)
    lam = TR.lambdas(D, REF)[AD.ARM_LAM[ARM]]
    neg = negatives(D); C = TR.controls(T6)
    GEO = {st: TR.read_xyz(T6 / "frame_42436" / fn) for st, fn in
           (("R", "reactant_qm.xyz"), ("TS", "transition_state_qm.xyz"))}
    # descriptors() wants cs13_qscale's labelled form, not cs6_transfer's array
    GLBL = Q13.read_xyz(T6 / "frame_42436" / "reactant_qm.xyz")
    O = {st: com(T6, st) for st in Q13.GEOM}
    idx = {s: k for k, s in enumerate(D["site"])}

    cand = [m for m in R.rglob("INPUT_CHECKSUMS.sha256") if ".git" not in m.parts]
    if len(cand) != 1:
        stop("expected exactly one INPUT_CHECKSUMS.sha256 under %s, found %d" % (R, len(cand)))
    base, bad = cand[0].parent, []
    for line in cand[0].read_text().splitlines():
        want, rel = line.split("  ", 1)
        f = base / rel
        if not f.is_file():
            bad.append("missing %s" % rel)
        elif sha(f) != want:
            bad.append("changed %s" % rel)
    if bad:
        stop("GATE 0 failed, %d input files do not match the manifest:\n   %s"
             % (len(bad), "\n   ".join(bad[:10])))

    E, CL, WHY = {}, {}, {}
    def take(tree, basis, nm, arm, st, mem, walled, key):
        d = tree / basis / (("ctl_%s" % st) if nm == "ctl"
                            else ("%s_%s" % (AD.rname(nm, arm, ARM if arm != "wall0" else None), st)))
        if arm == "bare":
            e, nb, dd = TR.one(d, GEO[st], len(mem))
        else:
            e, nb, dd = AD.one0(d, GEO[st], len(mem), walled)
        if walled:                       # A1.2, and only here
            dd = dict(dd); dd["dip"] = corrected_dip(d, mem, O[st])
        E[key], (CL[key], WHY[key]) = e, TR.verdict(dd, C)

    for basis in BASES:
        for st in Q13.GEOM:
            take(R, basis, "ctl", "none", st, [], False, (basis, "ctl", "none", st))
        for s in neg:
            nm = "s%d" % s
            mem, mem0 = AD.members(D, [s], lam), zeroed(AD.members(D, [s], 1.0))
            for st in Q13.GEOM:
                take(R, basis, nm, "wall0", st, mem0, True, (basis, nm, "wall0", st))
                take(R, basis, nm, "bare", st, mem, False, (basis, nm, "bare", st))
                take(R, basis, nm, "wall", st, mem, True, (basis, nm, "wall", st))
        # DIAGNOSTIC 14's own n6neg runs, committed, re-read here with A1.2 applied
        mem, mem0 = AD.members(D, neg, lam), zeroed(AD.members(D, neg, 1.0))
        for st in Q13.GEOM:
            take(A, basis, "ctl", "none", st, [], False, (basis, "D14ctl", "none", st))
            take(A, basis, "n6neg", "wall0", st, mem0, True, (basis, "n6neg", "wall0", st))
            take(A, basis, "n6neg", "bare", st, mem, False, (basis, "n6neg", "bare", st))
            take(A, basis, "n6neg", "wall", st, mem, True, (basis, "n6neg", "wall", st))

    H = TR.H
    def ddE(basis, nm, arm, ctl):
        c = {g: E[(basis, ctl, "none", g)] for g in Q13.GEOM}
        z = {g: (E[(basis, nm, "wall0", g)] - c[g]) if arm == "wall" else 0.0 for g in Q13.GEOM}
        d = {g: E[(basis, nm, arm, g)] - c[g] - z[g] for g in Q13.GEOM}
        return (d["TS"] - d["R"]) * H
    def dsh(basis, nm, ctl="ctl"):
        return ddE(basis, nm, "bare", ctl) - ddE(basis, nm, "wall", ctl)
    def well(basis, nm):
        return all(CL[(basis, nm, a, g)] in ("CLEAN", "POLARISED")
                   for a in ("wall0", "bare", "wall") for g in Q13.GEOM)

    O_ = ["STAGE 6 - DIAGNOSTIC 14, AMENDMENT 1 (A1.4): where does the negatives' non-inertness",
          "come from, and why only at def2-SVPD?", "",
          "Measured by cs6_negatives.py against %s, committed before" % AMEND,
          "any input of A1.4 existed, and applied unchanged. Arm %s only, lambda %.6f. No geometry" % (ARM, lam),
          "moved. In vacuo. D6.3's classification carries A1.2's correction: on a run whose charges sit",
          "on inline ECP centres the printed dipole has sum q_i (r_i - R_com) removed before the rule is",
          "applied. No tolerance changed.", ""]

    # GATE 1
    O_ += ["GATE 1: THESE CONTROLS MUST REPRODUCE DIAGNOSTIC 14'S, TO %g Eh" % CTL_EH, ""]
    g1 = True
    for basis in BASES:
        for st in Q13.GEOM:
            got, want = E[(basis, "ctl", "none", st)], E[(basis, "D14ctl", "none", st)]
            ok = abs(got - want) <= CTL_EH; g1 &= ok
            O_.append("   %-10s %-2s %18.10f   D14 %18.10f   %s"
                      % (basis, st, got, want, "pass" if ok else "FAIL"))
    O_.append("")
    if not g1:
        return finish(O_, outdir, "GATE 1 failed: the harness has drifted. Nothing below is read.")

    # GATE 2
    W = d14_n6neg(T6)
    O_ += ["GATE 2: THE n6neg d_shield MUST REPRODUCE DIAGNOSTIC 14'S PRINTED VALUE", "",
           "   A1.2 changes the classification, not the energies, so these must be unchanged.", ""]
    g2 = True
    Dset = {}
    for basis in BASES:
        got = dsh(basis, "n6neg", "D14ctl"); want = W[(basis, ARM)]
        Dset[basis] = abs(got)
        ok = abs(got - want) <= N6_TOL; g2 &= ok
        O_.append("   %-10s arm %s   %+8.4f   D14 %+8.4f   %s"
                  % (basis, ARM, got, want, "pass" if ok else "FAIL"))
    O_.append("")
    if not g2:
        return finish(O_, outdir, "GATE 2 failed: this harness disagrees with DIAGNOSTIC 14's own "
                                  "number. Nothing below is read.")

    # the classification, before any reading
    O_ += ["EVERY RUN CLASSIFIED BY D6.3 WITH A1.2'S CORRECTION, READ BEFORE ANY NUMBER", "",
           "   run                              class       note"]
    for basis in BASES:
        for k in sorted([k for k in CL if k[0] == basis], key=lambda k: (k[1], k[2], k[3])):
            O_.append("   %-32s %-10s  %s" % ("%s/%s_%s_%s" % (basis, k[1], k[2], k[3]),
                                              CL[k], "; ".join(WHY[k])))
    O_.append("")

    V1, V2, SH = {}, {}, {}
    for basis in BASES:
        dn = {}
        O_ += ["%s" % basis, "",
               "   site   nearest atom   d_shield   class"]
        rows = []
        for s in neg:
            nm = "s%d" % s
            dn[s] = Q13.descriptors(D["X"][idx[s]], GLBL)["dn"]
            if not well(basis, nm):
                O_.append("   %-5d  %8.3f A    NOT READABLE" % (s, dn[s]))
                rows.append((s, None)); continue
            v = dsh(basis, nm)
            rows.append((s, v))
            O_.append("   %-5d  %8.3f A   %+8.4f   %s" % (s, dn[s], v, CL[(basis, nm, "wall", "R")]))
        SH[basis] = dict(rows)
        unread = [s for s, v in rows if v is None]
        O_.append("")
        if unread or not well(basis, "n6neg"):
            V1[basis] = V2[basis] = "NOT READABLE"
            O_ += ["   NOT READABLE: %s" % (("sites " + ", ".join(str(s) for s in unread))
                                            if unread else "the committed n6neg runs"), ""]
            continue
        Dsum = sum(abs(v) for _, v in rows)
        tol = max(Q13.SHIELD_REL * Dsum, Q13.SHIELD_ABS)
        big = max(rows, key=lambda r: abs(r[1]))
        O_ += ["   Dsum %.4f   Dset %.4f   tol %.4f   Dset-Dsum %+.4f"
               % (Dsum, Dset[basis], tol, Dset[basis] - Dsum), ""]
        if Dset[basis] <= tol:                      # A2.2: nothing to locate, no site is named
            V1[basis] = "NO EFFECT TO LOCATE (Dset %.4f <= tol %.4f)" % (Dset[basis], tol)
        elif abs(abs(big[1]) - Dset[basis]) <= tol:
            V1[basis] = "ONE SITE (%d)" % big[0]
        elif Dsum + tol >= Dset[basis]:
            V1[basis] = "SPREAD"
        else:
            V1[basis] = "NEITHER"
        V2[basis] = "BOUNDED" if Dset[basis] <= Dsum + tol else "NOT BOUNDED"
        triv = "  (trivially satisfied: nothing here exceeds the tolerance)" \
               if Dset[basis] <= tol else ""        # A2.3 disclosure, not a criterion
        O_ += ["   READING 1: %s" % V1[basis], "   READING 2: %s%s" % (V2[basis], triv), ""]

    # READING 3
    O_ += ["READING 3: IS THE BASIS DIFFERENCE PER-SITE OR COLLECTIVE?", ""]
    a, b = BASES
    if "NOT READABLE" in (V1.get(a), V1.get(b)):
        v3 = "NOT READABLE"
    else:
        Dsum = sum(abs(SH[a][s]) for s in neg)
        tol = max(Q13.SHIELD_REL * Dsum, Q13.SHIELD_ABS)
        per = [s for s in neg if abs(SH[a][s] - SH[b][s]) > tol]
        setd = abs(Dset[a] - Dset[b])
        for s in neg:
            O_.append("   site %-5d  %-10s %+8.4f   %-10s %+8.4f   difference %+8.4f"
                      % (s, a, SH[a][s], b, SH[b][s], SH[a][s] - SH[b][s]))
        O_ += ["", "   tol %.4f;  |Dset(%s) - Dset(%s)| = %.4f" % (tol, a, b, setd), ""]
        v3 = ("PER-SITE (%s)" % ", ".join(str(s) for s in per)) if per else \
             ("COLLECTIVE" if setd > tol else "NEITHER")
    O_ += ["   READING 3: %s" % v3, "",
           "d_near is reported and nothing is concluded from it: ERRATUM3 E3.1 records that distance",
           "and magnitude are not yet separated, and this amendment does not separate them.", ""]
    return finish(O_, outdir, None)


def finish(O, outdir, tail):
    if tail:
        O += [tail, ""]
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "DIAGNOSTIC14_AMENDMENT1_RESULT.txt").write_text("\n".join(O) + "\n")
    print("\n".join(O))
    print("\n   wrote %s" % (outdir / "DIAGNOSTIC14_AMENDMENT1_RESULT.txt"))


if __name__ == "__main__":
    A_ = sys.argv[1:]
    if len(A_) == 4 and A_[0] == "prepare":
        prepare(A_[1], A_[2], A_[3])
    elif len(A_) == 4 and A_[0] == "analyse":
        analyse(A_[1], A_[2], A_[3])
    else:
        stop("usage: cs6_negatives.py prepare <results_repo> <OUT> <HPC_DIR>\n"
             "       cs6_negatives.py analyse <results_repo> <RUNS> <OUT_DIR>")
