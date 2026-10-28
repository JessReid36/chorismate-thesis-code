#!/usr/bin/env python3
"""
cs6_additivity.py - STAGE6_DIAGNOSTIC14_20261028.txt.

Does the per-site artefact bound a many-charge design?  S6-10's own sites, nested, at two
scalings, bare point charges against the same charges behind a Pauli wall.

  python3 cs6_additivity.py prepare <results_repo> <OUT> <HPC_DIR>
  python3 cs6_additivity.py analyse <results_repo> <RUNS> <OUT_DIR>

NOTHING ALREADY COMMITTED IS RETYPED HERE.
  cs6_transfer.py  design()   the design, its sites, charges, a and L, parsed out of the committed
                              relax/ARRANGEMENTS.txt and checked against that file's own recorded
                              score and sum q^2
                   lambdas()  arm B is S6-20's matched lambda, computed, not typed
                   head()     the two headers, already checked byte for byte against the committed
                              cs5_amatrix HEAD and the committed qscan C5 header
                   controls() D6's published control medians per basis-set size
                   verdict()  D6.3's rule, thresholds from cs6_audit
                   one()      the output reader, for the runs that carry point charges
  cs13_qscale.py   ECP_TAIL, FRAME_SHA, GEOM, read_xyz_raw, descriptors, SHIELD_REL, SHIELD_ABS
  cs6_audit.py     read()     reached through cs6_transfer

Only what DIAGNOSTIC 14 adds is written here: the nine site sets, and D14.3's reading.
"""
import hashlib, re, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import cs13_qscale as Q13
import cs6_transfer as TR

PREREG = "STAGE6_DIAGNOSTIC14_20261028.txt"
DESIGN, NSPH = "S6-10", 10
ARM_LAM = {"A": "lam050", "B": "lamM2"}     # D14.2's two arms, named by D7's own lambda tags
POS = (1695, 2152, 1089, 1234)              # the positives, largest |a| first; the nesting order
ARM_BASES = {"def2-SVPD": ("A", "B"), "def2-SVP": ("B",)}   # D14.2: def2-SVP is a control, arm B only
RATIO = ("n2", "n3", "n4")
REPRO_TOL = 1e-4                            # kcal/mol, D14.3's three reproduction checks
CTL_EH = 1e-6                               # the def2-SVP control against A_REPORT.txt
N6_ABS = 0.05                               # D14.3: n6neg must be inert to this


def stop(m):
    print("STOP: " + m); sys.stdout.flush(); sys.exit(1)


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def sets_of(D):
    """D14.2's nine site sets, built from the committed design. Only POS is typed."""
    all_sites = list(D["site"])
    for s in POS:
        if s not in all_sites:
            stop("site %d is not in the committed %s" % (s, DESIGN))
        if D["q"][all_sites.index(s)] <= 0:
            stop("site %d is not positive in the committed %s" % (s, DESIGN))
    neg = [s for s in all_sites if D["q"][all_sites.index(s)] < 0]
    if len(neg) != 6:
        stop("expected six negative sites in %s, parsed %d" % (DESIGN, len(neg)))
    S = {"s%d" % s: [s] for s in POS}
    S["n2"], S["n3"], S["n4"] = list(POS[:2]), list(POS[:3]), list(POS[:4])
    S["n6neg"], S["n10"] = neg, all_sites
    return S


def members(D, sites, lam):
    i = {s: k for k, s in enumerate(D["site"])}
    return [(s, lam * float(D["q"][i[s]]), D["X"][i[s]]) for s in sites]


def rname(name, arm, lam):
    return "%s_%s" % (name, arm) if arm == "wall0" else "%s_%s_%s" % (name, arm, lam)


# ---------------------------------------------------------------- prepare

def write_run(d, basis, geo_text, raw, mem, arm):
    """One run directory. The bare arm's three files are byte-identical to D7's, which the
    n10 gate below checks against the committed transfer/ tree."""
    d.mkdir(parents=True)
    (d / "substrate.xyz").write_text(geo_text)     # the committed geometry file, copied verbatim
    if arm == "none":
        (d / "job.inp").write_text(TR.head(basis) + "* xyzfile -2 1 substrate.xyz\n")
        return
    if arm == "bare":
        lines = ["%d" % len(mem)]
        for _, q, X in mem:
            lines.append("%10.6f %18.10f %18.10f %18.10f" % (q, X[0], X[1], X[2]))
        (d / "pc.pc").write_text("\n".join(lines) + "\n")
        (d / "job.inp").write_text(TR.head(basis) + '%pointcharges "pc.pc"\n'
                                                    '* xyzfile -2 1 substrate.xyz\n')
        return
    body = ["* xyz -2 1"] + list(raw)
    for _, q, X in mem:
        body.append("  Ne>  %.6f  %18.8f %18.8f %18.8f  NewECP" % (0.0 if arm == "wall0" else q,
                                                                   X[0], X[1], X[2]))
        body.append(Q13.ECP_TAIL)
    body.append("*")
    (d / "job.inp").write_text(TR.head(basis) + "\n".join(body) + "\n")


PBS = """#!/bin/bash
#PBS -N cm_d14%(tag)s
#PBS -l select=1:ncpus=8:mem=24gb
#PBS -l walltime=24:00:00
#PBS -m ae
#PBS -M 18660916@sun.ac.za
#PBS -j oe
#PBS -o %(hpc)s/%(basis)s/run.pbs.out
# Stage 6 DIAGNOSTIC 14 (cs6_additivity.py): %(n)d single points at %(basis)s. Restartable.
# No -q is given: the allocation is automatic, and walltime 24:00:00 routes through workq to day.
export PATH="/apps/openmpi/4.1.1/bin:$PATH"
export LD_LIBRARY_PATH="/home/apps2/ORCA/6.0.1/lib:/apps/openmpi/4.1.1/lib:/apps/mambaforge/envs/medaka/lib:${LD_LIBRARY_PATH:-}"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
ORCA=/home/apps2/ORCA/6.0.1
W=%(hpc)s/%(basis)s
cd "$W" || exit 1
if ! mkdir .running 2>/dev/null; then echo "FAIL: another job holds .running in $W"; exit 1; fi
trap 'rmdir "$W/.running" 2>/dev/null' EXIT
echo "host=$(hostname) start=$(date)"
run() {
  d="$W/$1"
  if [ -s "$d/job.out" ] && grep -q "ORCA TERMINATED NORMALLY" "$d/job.out"; then echo "$1 done before"
  else ( cd "$d" && "$ORCA/orca" job.inp > job.out 2>&1 </dev/null ); echo "$1 orca_exit=$? normal=$(grep -c 'ORCA TERMINATED NORMALLY' "$d/job.out") $(date)"; fi
}
%(runs)s
echo "end=$(date)"
"""


def prepare(root, out, hpc):
    T6 = Path(root) / "phase2.2" / "cs_stage6"
    out = Path(out)
    if out.exists():
        stop("%s already exists - this never overwrites" % out)
    if not (T6 / PREREG).is_file():
        stop("%s is not committed. The pre-registration goes in before any input exists." % PREREG)
    print("   %s is committed (sha256 %s)" % (PREREG, sha(T6 / PREREG)[:16]))

    geo, raw = {}, {}
    for st, fn in (("R", "reactant_qm.xyz"), ("TS", "transition_state_qm.xyz")):
        p = T6 / "frame_42436" / fn
        if sha(p) != Q13.FRAME_SHA[fn]:
            stop("%s does not match the sha256 cs13_qscale records" % fn)
        geo[st], raw[st] = p.read_text(), Q13.read_xyz_raw(p)
    print("   both geometries match the sha256 the committed scripts record")

    D, REF = TR.design(T6, DESIGN, NSPH), TR.design(T6, TR.REF_DESIGN, TR.REF_N)
    LAM = TR.lambdas(D, REF)
    for a, t in ARM_LAM.items():
        if t not in LAM:
            stop("cs6_transfer.lambdas() does not produce %s, which arm %s is defined as" % (t, a))
    print("   %s: sum a q %+.5f  sum q^2 %.5f  sum |q| %.4f e  score %+.5f"
          % (DESIGN, D["saq"], D["sq2"], D["sabs"], D["score"]))
    for a in ("A", "B"):
        lam = LAM[ARM_LAM[a]]
        print("   ARM %s  lambda %.6f  sum |q| %.4f e  largest |q| %.4f e"
              % (a, lam, lam * D["sabs"], lam * D["maxq"]))
    # the two headers, checked against their committed counterparts exactly as cs6_transfer does
    if TR.head("def2-SVP") != ("! B3LYP D3BJ def2-SVP def2/J RIJCOSX TightSCF\n%maxcore 3000\n"
                               "%pal nprocs 8 end\n%scf MaxIter 300 end\n"):
        stop("the def2-SVP header is not cs5_amatrix.py's HEAD byte for byte")
    ref = T6 / "qscan" / "C5" / "qp1000_R" / "job.inp"
    if TR.head("def2-SVPD") != ref.read_text().split("%pointcharges")[0]:
        stop("the def2-SVPD header is not the committed qscan C5 header byte for byte")
    print("   both headers match their committed counterparts byte for byte")
    S = sets_of(D)
    C = TR.controls(T6)
    for nb in (264, 402):
        if nb not in C:
            stop("D6's control table has no nb %d row; the classification could not be anchored" % nb)
    print("   D6's control medians are present for nb 264 and nb 402")

    out.mkdir(parents=True)
    plan, bybasis = [], {}
    for basis, arms in ARM_BASES.items():
        rs = []
        for name in sorted(S):
            sites = S[name]
            for st in Q13.GEOM:
                mem0 = members(D, sites, 1.0)
                rn = "%s_%s" % (rname(name, "wall0", None), st)
                write_run(out / basis / rn, basis, geo[st], raw[st], mem0, "wall0")
                rs.append(rn); plan.append((basis, name, "wall0", "-", st, len(sites), 0.0)); 
                for a in arms:
                    lam = LAM[ARM_LAM[a]]
                    mem = members(D, sites, lam)
                    for arm in ("bare", "wall"):
                        rn = "%s_%s" % (rname(name, arm, a), st)
                        write_run(out / basis / rn, basis, geo[st], raw[st], mem, arm)
                        rs.append(rn)
                        plan.append((basis, name, arm, a, st, len(sites),
                                     sum(abs(q) for _, q, _ in mem)))
        for st in Q13.GEOM:
            rn = "ctl_%s" % st
            write_run(out / basis / rn, basis, geo[st], raw[st], [], "none")
            rs.append(rn); plan.append((basis, "ctl", "none", "-", st, 0, 0.0))
        bybasis[basis] = rs
        (out / basis / "run.pbs").write_text(
            PBS % dict(tag=("d" if basis.endswith("D") else ""), hpc=hpc, basis=basis,
                       n=len(rs), runs="\n".join("run %s" % r for r in rs)))

    # THE GATE THAT MAKES THIS WORTH DOING: n10's bare arm must be D7's runs, byte for byte.
    print("   checking n10's bare runs against the committed DIAGNOSTIC 7 inputs, byte for byte")
    nchk = 0
    for basis, arms in ARM_BASES.items():
        for a in arms:
            for st in Q13.GEOM:
                mine = out / basis / ("n10_bare_%s_%s" % (a, st))
                d7 = T6 / "transfer" / basis / ("%s_%s" % (ARM_LAM[a], st))
                if not d7.is_dir():
                    stop("the committed D7 run %s is missing; byte-identity cannot be checked" % d7)
                for f in ("job.inp", "pc.pc", "substrate.xyz"):
                    if sha(mine / f) != sha(d7 / f):
                        stop("n10 arm %s %s %s/%s is NOT byte-identical to the committed D7 input "
                             "at %s. D14.2 claims it is; either this generator or that claim is "
                             "wrong, and nothing is submitted until that is resolved."
                             % (a, basis, mine.name, f, d7))
                    nchk += 1
    print("   PASS: %d files across %d n10 bare runs are byte-identical to DIAGNOSTIC 7's"
          % (nchk, nchk // 3))

    with open(out / "PLAN.tsv", "w") as fh:
        fh.write("# %s; cs6_additivity.py prepare\n" % PREREG)
        fh.write("# %s from relax/ARRANGEMENTS.txt: sum a q %+.5f  sum q^2 %.5f  sum |q| %.4f e\n"
                 % (DESIGN, D["saq"], D["sq2"], D["sabs"]))
        fh.write("# arm A lambda %.6f (%s); arm B lambda %.6f (%s), computed by cs6_transfer.lambdas()\n"
                 % (LAM[ARM_LAM["A"]], ARM_LAM["A"], LAM[ARM_LAM["B"]], ARM_LAM["B"]))
        fh.write("basis\tset\tarm\tlambda_arm\tstate\tn\tsum_abs_q\n")
        for r in plan:
            fh.write("%s\t%s\t%s\t%s\t%s\t%d\t%.4f\n" % r)
    files = sorted(p for p in out.rglob("*") if p.is_file())
    (out / "INPUT_CHECKSUMS.sha256").write_text(
        "".join("%s  %s\n" % (sha(p), p.relative_to(out)) for p in files))
    print("   %d run directories, %d files including the manifest" % (len(plan), len(files) + 1))
    for b, rs in bybasis.items():
        print("      %-10s %3d runs, one restartable PBS job" % (b, len(rs)))
    print("   NOTHING WAS SUBMITTED BY THIS BUNDLE.")


# ---------------------------------------------------------------- analyse

def one0(d, ref_xyz, nsph, walled):
    """cs6_transfer.one() for a run that carries NO point charge file.

    A structural transcription of that committed function, which cannot be reused here because its
    point-charge assertion is unconditional. Every check below is its check, in its order, with two
    changes, both of them additions: the output must NOT report a point charge file, and a walled
    run must echo exactly one Ne> centre per sphere. Do not re-tune anything here.
    """
    p = Path(d) / "job.out"
    if not p.is_file():
        stop("%s: no job.out - has the batch run?" % d)
    t = p.read_text(errors="replace")
    if "****ORCA TERMINATED NORMALLY****" not in t:
        stop("%s: ORCA did not terminate normally" % d)
    if not re.search(r"SCF CONVERGED AFTER\s+\d+\s+CYCLES", t) or "SCF NOT CONVERGED" in t:
        stop("%s: the SCF did not converge" % d)
    m = re.search(r"Number of Electrons\s+NEL\s+\.+\s+(\d+)", t)
    if not m or int(m.group(1)) != TR.NEL:
        stop("%s: %s electrons, expected %d" % (d, m.group(1) if m else "no count", TR.NEL))
    e = re.findall(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", t)
    if len(e) != 1:
        stop("%s: %d final energies, expected 1" % (d, len(e)))
    if re.search(r"Reading point charge file", t):
        stop("%s: this run reports a point charge file and must not have one" % d)
    k = t.find("CARTESIAN COORDINATES (ANGSTROEM)")
    if k < 0:
        stop("%s: no echoed coordinates" % d)
    block = []
    for l in t[k:].splitlines()[2:]:
        if not l.strip():
            break
        block.append(l)
    X = np.array([[float(v) for v in l.split()[1:4]] for l in block[:len(ref_xyz)]])
    if X.shape != ref_xyz.shape or np.abs(X - ref_xyz).max() > 1e-5:
        stop("%s: the echoed coordinates differ from the committed geometry" % d)
    ne = sum(1 for l in block if l.split() and l.split()[0] == "Ne>")
    if ne != (nsph if walled else 0):
        stop("%s: %d Ne> centres echoed, expected %d" % (d, ne, nsph if walled else 0))
    if len(block) != len(ref_xyz) + ne:
        stop("%s: %d centres echoed, expected %d" % (d, len(block), len(ref_xyz) + ne))
    nb = re.search(r"Number of basis functions\s+\.+\s+(\d+)", t)
    if not nb:
        stop("%s: no basis-function count" % d)
    dd = TR.AUD.read(t)
    if dd is None:
        stop("%s: cs6_audit.read() does not accept this output" % d)
    return float(e[0]), int(nb.group(1)), dd


def d7_table(T6):
    """D7's committed ddE per basis at each lambda tag, read out of its own result table."""
    p = T6 / "transfer" / "DIAGNOSTIC7_RESULT.txt"
    if not p.is_file():
        stop("%s is missing: the three reproduction checks have nothing to check against" % p)
    out = {}
    for l in p.read_text().splitlines():
        m = re.match(r"\s*(lam\S+)\s+(\d\.\d{6})\s+\d+\.\d+\s+\d+\.\d+\s+[+-]?\d+\.\d+"
                     r"\s+([+-]?\d+\.\d+)\s+\d+\.\d+%\s+([+-]?\d+\.\d+)\s+\d+\.\d+%"
                     r"\s+[+-]?\d+\.\d+\s*$", l)
        if m:
            out[("def2-SVP", m.group(1))] = (float(m.group(2)), float(m.group(3)))
            out[("def2-SVPD", m.group(1))] = (float(m.group(2)), float(m.group(4)))
    if not out:
        stop("could not parse DIAGNOSTIC7_RESULT.txt's lambda table")
    return out


def a_report(T6):
    for l in (T6 / "A_REPORT.txt").read_text().splitlines():
        f = l.split("\t")
        if f and f[0].strip() == "42436" and len(f) >= 3:
            return {"R": float(f[1]), "TS": float(f[2])}
    stop("A_REPORT.txt has no frame 42436 row to check the control against")


def analyse(root, runs, outdir):
    T6, R, outdir = Path(root) / "phase2.2" / "cs_stage6", Path(runs), Path(outdir)
    if not (T6 / PREREG).is_file():
        stop("%s is not committed" % PREREG)
    D, REF = TR.design(T6, DESIGN, NSPH), TR.design(T6, TR.REF_DESIGN, TR.REF_N)
    LAM, S, C = TR.lambdas(D, REF), sets_of(D), TR.controls(T6)
    GEO = {st: TR.read_xyz(T6 / "frame_42436" / fn) for st, fn in
           (("R", "reactant_qm.xyz"), ("TS", "transition_state_qm.xyz"))}

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

    E, NB, CL, WHY = {}, {}, {}, {}
    def take(basis, name, arm, a, st, sites, walled):
        d = base / basis / (("ctl_%s" % st) if name == "ctl"
                            else ("%s_%s" % (rname(name, arm, a), st)))
        if arm == "bare":
            e, nb, dd = TR.one(d, GEO[st], len(sites))
        else:
            e, nb, dd = one0(d, GEO[st], len(sites), walled)
        k = (basis, name, arm, a, st)
        E[k], NB[k] = e, nb
        CL[k], WHY[k] = TR.verdict(dd, C)
    for basis, arms in ARM_BASES.items():
        for st in Q13.GEOM:
            take(basis, "ctl", "none", "-", st, [], False)
        for name in sorted(S):
            for st in Q13.GEOM:
                take(basis, name, "wall0", "-", st, S[name], True)
                for a in arms:
                    take(basis, name, "bare", a, st, S[name], False)
                    take(basis, name, "wall", a, st, S[name], True)

    H = TR.H
    def ddE(basis, name, arm, a):
        c = {g: E[(basis, "ctl", "none", "-", g)] for g in Q13.GEOM}
        z = {g: (E[(basis, name, "wall0", "-", g)] - c[g]) if arm == "wall" else 0.0
             for g in Q13.GEOM}     # Stage 4c Amendment 1: ddE_corr = ddE - ddE(zero charge)
        d = {g: E[(basis, name, arm, a, g)] - c[g] - z[g] for g in Q13.GEOM}
        return (d["TS"] - d["R"]) * H
    def dshield(basis, name, a):
        return ddE(basis, name, "bare", a) - ddE(basis, name, "wall", a)
    def well(basis, name, a):
        ks = [(basis, name, "wall0", "-", g) for g in Q13.GEOM]
        for arm in ("bare", "wall"):
            ks += [(basis, name, arm, a, g) for g in Q13.GEOM]
        return all(CL[k] in ("CLEAN", "POLARISED") for k in ks)

    O = ["STAGE 6 - DIAGNOSTIC 14 RESULT: does the per-site artefact bound a design?", "",
         "Measured by cs6_additivity.py against %s, committed before any input" % PREREG,
         "existed, and applied unchanged. No geometry moved: single points on Diagnostic 1's",
         "committed R and TS, frame 42436. In vacuo, no continuum (STAGE6_AMENDMENT2_20261027).",
         "",
         "The design, its sites and charges come from the committed relax/ARRANGEMENTS.txt through",
         "cs6_transfer.design(); arm B's lambda from cs6_transfer.lambdas(); D6.3's classification",
         "from cs6_transfer.verdict() on cs6_audit.read(), against the committed",
         "DIAGNOSTIC6_RESULT.txt control medians. ddE is corrected by the zero-charge wall term as",
         "STAGE4C_AMENDMENT1 requires (ddE_corr = ddE - ddE(zero charge)); the bare arm needs no",
         "such correction and gets none.",
         "",
         "   ARM A  lambda %.6f   sum |q| %.4f e   largest |q| %.4f e"
         % (LAM[ARM_LAM["A"]], LAM[ARM_LAM["A"]] * D["sabs"], LAM[ARM_LAM["A"]] * D["maxq"]),
         "   ARM B  lambda %.6f   sum |q| %.4f e   largest |q| %.4f e"
         % (LAM[ARM_LAM["B"]], LAM[ARM_LAM["B"]] * D["sabs"], LAM[ARM_LAM["B"]] * D["maxq"]), ""]

    # GATE 1 - the control
    AR = a_report(T6)
    O += ["GATE 1: THE def2-SVP CONTROL MUST REPRODUCE A_REPORT.TXT'S FRAME 42436 ENERGIES", ""]
    g1 = True
    for g in Q13.GEOM:
        got, want = E[("def2-SVP", "ctl", "none", "-", g)], AR[g]
        ok = abs(got - want) <= CTL_EH
        g1 &= ok
        O.append("   %-3s %18.10f   committed %18.10f   %s" % (g, got, want, "pass" if ok else "FAIL"))
    O.append("")
    if not g1:
        return finish(O, outdir, "GATE 1 failed: the control is not the committed one. Nothing below is read.")

    # GATE 2 - the three reproduction checks
    D7 = d7_table(T6)
    O += ["GATE 2: THE THREE n10 REPRODUCTION CHECKS (D14.3)", "",
          "   n10's bare arm is byte-identical to DIAGNOSTIC 7's runs, which cs6_additivity.py",
          "   prepare checked file by file. These check the energies agree to %g kcal/mol." % REPRO_TOL, ""]
    g2 = True
    for basis, arms in ARM_BASES.items():
        for a in arms:
            lam7, want = D7[(basis, ARM_LAM[a])]
            if abs(lam7 - LAM[ARM_LAM[a]]) > 1e-6:
                stop("D7's %s lambda is %.6f, cs6_transfer.lambdas() gives %.6f"
                     % (ARM_LAM[a], lam7, LAM[ARM_LAM[a]]))
            got = ddE(basis, "n10", "bare", a)
            ok = abs(got - want) <= REPRO_TOL
            g2 &= ok
            O.append("   %-10s arm %s (%-6s)  ddE %+9.4f   D7 %+9.4f   %s"
                     % (basis, a, ARM_LAM[a], got, want, "pass" if ok else "FAIL"))
    O.append("")
    if not g2:
        return finish(O, outdir, "GATE 2 failed: n10 does not reproduce DIAGNOSTIC 7. Nothing below is read.")

    # the classification, read before any energy
    O += ["EVERY RUN CLASSIFIED BY D6.3, READ BEFORE ANY RATIO", "",
          "   run                              nb   class       note"]
    for basis in ARM_BASES:
        for k in sorted([k for k in CL if k[0] == basis], key=lambda k: (k[1], k[3], k[2], k[4])):
            O.append("   %-32s %3d   %-10s  %s"
                     % ("%s/%s" % (basis, ("ctl_%s" % k[4]) if k[1] == "ctl"
                                   else "%s_%s" % (rname(k[1], k[2], k[3]), k[4])),
                        NB[k], CL[k], "; ".join(WHY[k])))
    O.append("")

    # CONTROL - n6neg
    O += ["CONTROL: n6neg MUST BE INERT (D14.3)", "",
          "   the six negatives together, |d_shield| <= %.2f kcal/mol, CLEAN at both bases" % N6_ABS, ""]
    g3 = True
    for basis, arms in ARM_BASES.items():
        for a in arms:
            ks = [(basis, "n6neg", "wall0", "-", g) for g in Q13.GEOM]
            for arm in ("bare", "wall"):
                ks += [(basis, "n6neg", arm, a, g) for g in Q13.GEOM]
            cls = sorted({CL[k] for k in ks})
            ds = dshield(basis, "n6neg", a)
            ok = (cls == ["CLEAN"]) and abs(ds) <= N6_ABS
            g3 &= ok
            O.append("   %-10s arm %s   d_shield %+8.4f   %-28s %s"
                     % (basis, a, ds, "/".join(cls), "pass" if ok else "FAIL"))
    O.append("")
    if not g3:
        return finish(O, outdir,
                      "THE n6neg CONTROL FAILED. D14.3: if the negative set is not inert, the sign "
                      "asymmetry is wrong and this diagnostic stops without a verdict.")

    # READING 1
    O += ["READING 1: IS THE SET BOUNDED BY THE SUM OF ITS MEMBERS?", "",
          "   Dsum = sum over members of |d_shield| measured alone; Dset = |d_shield| of the set",
          "   together; tol = max(0.10 x Dsum, %.2f), DIAGNOSTIC 11.3's tolerance, unchanged." % Q13.SHIELD_ABS,
          "",
          "   A ratio point is READABLE only when the set's six runs AND every member single's six",
          "   runs are CLEAN or POLARISED. Requiring the singles too is this harness's reading of",
          "   D14.3's 'any run in a subset', and it is the conservative one: Dsum is built from",
          "   those singles and an unwell single would carry into it. That choice is disclosed here.",
          ""]
    V1, TAB = {}, {}
    for basis, arms in ARM_BASES.items():
        for a in arms:
            O.append("   %s, arm %s (lambda %.6f)" % (basis, a, LAM[ARM_LAM[a]]))
            O.append("      set  Dsum      Dset      tol      Dset-Dsum   ratio    verdict")
            res = []
            for nm in RATIO:
                mem = S[nm]
                bad = [x for x in ["s%d" % s for s in mem] + [nm] if not well(basis, x, a)]
                if bad:
                    O.append("      %-4s  %s" % (nm, "NOT READABLE: %s unwell" % ", ".join(bad)))
                    continue
                dsum = sum(abs(dshield(basis, "s%d" % s, a)) for s in mem)
                dset = abs(dshield(basis, nm, a))
                tol = max(Q13.SHIELD_REL * dsum, Q13.SHIELD_ABS)
                ok = dset <= dsum + tol
                res.append(ok)
                TAB[(basis, a, nm)] = (dsum, dset, tol, ok)
                O.append("      %-4s  %8.4f  %8.4f  %7.4f  %+9.4f   %6.3f   %s"
                         % (nm, dsum, dset, tol, dset - dsum,
                            (dset / dsum if dsum else float("nan")),
                            "bounded" if ok else "NOT BOUNDED"))
            if len(res) < 2:
                V1[(basis, a)] = "NO VERDICT"
                O.append("      fewer than two readable ratio points -> NO VERDICT")
            else:
                V1[(basis, a)] = "BOUNDED" if all(res) else "NOT BOUNDED"
                O.append("      verdict: %s" % V1[(basis, a)])
            O.append("")

    # READING 2
    O += ["READING 2: DOES THE ANSWER CHANGE BETWEEN THE TWO TOTAL FIELDS?", ""]
    b = "def2-SVPD"
    both = [nm for nm in RATIO if (b, "A", nm) in TAB and (b, "B", nm) in TAB]
    va, vb = V1.get((b, "A")), V1.get((b, "B"))
    if len(both) < 2:
        v2 = "NOT COMPARABLE"
    elif va == vb:
        v2 = "CONSTANT"
    elif va == "BOUNDED" and vb == "NOT BOUNDED":
        v2 = "TURNS ON"
    elif va == "NOT BOUNDED" and vb == "BOUNDED":
        v2 = "TURNS OFF"
    else:
        v2 = "NOT COMPARABLE"
    O += ["   at %s: arm A %s, arm B %s, over %d ratio point(s) readable in both" % (b, va, vb, len(both)),
          "", "   READING 2: %s" % v2, ""]
    if "def2-SVP" in ARM_BASES:
        O += ["   def2-SVP is a control and carries arm B only: %s" % V1.get(("def2-SVP", "B")), ""]
    return finish(O, outdir, None)


def finish(O, outdir, tail):
    if tail:
        O += [tail, ""]
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "DIAGNOSTIC14_RESULT.txt").write_text("\n".join(O) + "\n")
    print("\n".join(O))
    print("\n   wrote %s" % (outdir / "DIAGNOSTIC14_RESULT.txt"))


if __name__ == "__main__":
    A = sys.argv[1:]
    if len(A) == 4 and A[0] == "prepare":
        prepare(A[1], A[2], A[3])
    elif len(A) == 4 and A[0] == "analyse":
        analyse(A[1], A[2], A[3])
    else:
        stop("usage: cs6_additivity.py prepare <results_repo> <OUT> <HPC_DIR>\n"
             "       cs6_additivity.py analyse <results_repo> <RUNS> <OUT_DIR>")
