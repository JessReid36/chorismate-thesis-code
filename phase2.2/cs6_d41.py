#!/usr/bin/env python3
"""DIAGNOSTIC 41's reader: the DLPNO-CCSD(T) barrier at def2-SVPD on frame 42436.

One new run, T, the transition state. The reactant leg is d38/N/job.out as committed and is
required by its sha256 (d41/REACTANT_LEG.sha256). Both legs, and D36's two legs at
def2-TZVPD (d36a2/B_R and B_TS), are read by D38's committed reader, cs6_d38.py read_out,
imported by path, so the reactant is read here exactly as D38 read it.

Nothing is typed in except committed reference values, each named with its source, and
each is checked against the committed file it comes from before it is used.

Any gate failure stops the reader and nothing is written.

usage: cs6_d41.py --res <results-repo> --code <code-repo> --prereg <file> --result <file> --tsv <file>
"""
import argparse, hashlib, importlib.util, os, re, subprocess, sys, textwrap

H2K = 627.5094740631
S = "phase2.2/cs_stage6"
NEL = 118
N_CCSDT = -834.533005232                    # G41.3, d38/N, STAGE6_DIAGNOSTIC41_20261115.txt
TOL_ID = 2e-9                               # G41.5
D36_RESULT = "STAGE6_DIAGNOSTIC36_RESULT_20261114.txt"
D36_BARRIER = "23.8603"                     # D36's result, R36.1
D36_T1 = (("R", "0.015649748"), ("TS", "0.018809319"))      # D36's result, R36.4
D38N_T1, D38N_HOMO = "0.016726026", "-2.0444"               # d38/N, D38's result
D36_PHRASE = "a basis uncertainty of roughly one kcal/mol"
SVPD = [("B3LYP-D3BJ", "18.3323"), ("CAM-B3LYP-D3BJ", "26.8248"), ("wB97X-D3", "29.0441"),
        ("wB97M-V", "28.6953"), ("wB97X-V", "30.7479"), ("PW6B95-D3BJ", "22.1571"),
        ("M06-2X-D3(0)", "29.4427"), ("PWPB95-D3BJ", "22.2031")]   # R41.3, as the pre-registration lists them
RANGE = ("-0.8624", "-1.3028")              # the eight functionals' basis shifts, pre-registration section 1
THREE = ("B3LYP-D3BJ", "CAM-B3LYP-D3BJ", "wB97X-D3")           # D36's pre-registered three
GUARDED = ("d33a", "d36", "d36a2", "d38")
COMP = (("E(0), Hartree-Fock", "e0"), ("E(CORR)(strong-pairs)", "strong"), ("E(CORR)(weak-pairs)", "weak"),
        ("E(CORR)(corrected)", "corr"), ("E(CCSD)", "ccsd"), ("Triples Correction (T)", "tcorr"),
        ("E(CCSD(T))", "ccsdt"))


def stop(m):
    sys.stderr.write("STOP  %s\n" % m)
    raise SystemExit(1)


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 16), b""):
            h.update(b)
    return h.hexdigest()


def orca_words(t):
    L = t.splitlines()
    err = [l.strip() for l in L if re.search(r"error|aborting|not converged|No space left|cannot", l, re.I)]
    if err:
        return err[-1][:200]
    return "no error line found; the output ends: " + (L[-1].strip()[:120] if L else "(empty)")


def job_facts(p):
    """R41.6, from the job script's own scheduler output."""
    if not os.path.isfile(p):
        return {}
    t = open(p, errors="replace").read()
    m = re.search(r"R41\.6 scratch=(\S+) avail=(\d+)G margin=(\S+)", t)
    k = re.search(r"R41\.6 sampled maximum while running: (\d+) MB", t)
    z = re.search(r"R41\.6 largest residue after the run:\s+(\d+) MB", t)
    h = re.search(r"host=(\S+)", t)
    return dict(scratch=m.group(1) if m else None, avail=int(m.group(2)) if m else None,
                margin=m.group(3) if m else None, peak=int(k.group(1)) if k else None,
                residue=int(z.group(1)) if z else None, host=h.group(1) if h else None,
                refused=re.findall(r"^   candidate (\S+)\s+(.*)$", t, re.M))


def main():
    ap = argparse.ArgumentParser()
    for a in ("res", "code", "prereg", "result", "tsv"):
        ap.add_argument("--" + a, required=True)
    a = ap.parse_args()
    S6 = os.path.join(a.res, S)
    D = os.path.join(S6, "d41")
    p38 = os.path.join(a.code, "phase2.2", "cs6_d38.py")
    sp = importlib.util.spec_from_file_location("d41_cs6_d38", p38)
    d38 = importlib.util.module_from_spec(sp)
    sp.loader.exec_module(d38)

    rep, gates = [], []
    say = rep.append
    say("STAGE 6, DIAGNOSTIC 41 - RESULT")
    say("The correlated barrier at the production basis: DLPNO-CCSD(T)/def2-SVPD, frame 42436.")
    say("")
    say("Read by phase2.2/cs6_d41.py               sha256 %s" % sha(__file__))
    say("   every output by cs6_d38.py read_out     sha256 %s" % sha(p38))
    say("   pre-registration %-22s sha256 %s" % (os.path.basename(a.prereg), sha(a.prereg)))
    say("")
    say("One new single point, the transition state (d41/T). The reactant leg is d38/N as")
    say("committed and was not run again. No charges, no MM region, no wall.")
    say("")

    # ------------------------------------------------------------------ gates
    geo = {leg: os.path.join(S6, "d33a", "sp_42436_%s" % leg, "substrate.xyz") for leg in ("R", "TS")}
    for lab, p, leg in (("T/substrate.xyz", os.path.join(D, "T", "substrate.xyz"), "TS"),
                        ("d38/N/substrate.xyz", os.path.join(S6, "d38", "N", "substrate.xyz"), "R")):
        if not (os.path.isfile(p) and os.path.isfile(geo[leg])) or sha(p) != sha(geo[leg]):
            gates.append("G41.1 %s is not d33a/sp_42436_%s/substrate.xyz" % (lab, leg))
    x, y = os.path.join(D, "T", "job.inp"), os.path.join(S6, "d38", "N", "job.inp")
    if not (os.path.isfile(x) and os.path.isfile(y)) or sha(x) != sha(y):
        gates.append("G41.2 T/job.inp is not byte-identical to d38/N/job.inp")

    N = d38.read_out(os.path.join(S6, "d38", "N"))
    pr = os.path.join(D, "REACTANT_LEG.sha256")
    want = open(pr).read().split() if os.path.isfile(pr) else []
    if len(want) != 2 or want[1] != "d38/N/job.out":
        gates.append("G41.3 d41/REACTANT_LEG.sha256 is missing or does not name d38/N/job.out")
    elif N is None or N["sha"] != want[0]:
        gates.append("G41.3 d38/N/job.out is not the file whose sha256 d41/REACTANT_LEG.sha256 records")
    if N is None or N.get("ccsdt") is None or abs(N["ccsdt"] - N_CCSDT) > 5e-10:
        gates.append("G41.3 d38/N/job.out does not print E(CCSD(T)) = %.9f" % N_CCSDT)

    T = d38.read_out(os.path.join(D, "T"))
    if T is None:
        gates.append("G41.4 T produced no output")
    else:
        words = orca_words(open(T["path"], errors="replace").read())
        if not T["norm"]:
            gates.append("G41.4 T did not terminate normally. ORCA's words: %s" % words)
        if T["notconv"]:
            gates.append("G41.4 T did not converge its SCF")
        if T["nel"] != NEL:
            gates.append("G41.4 T reports %s electrons, not %d" % (T["nel"], NEL))
        if T["ccsdt"] is None:
            gates.append("G41.4 T prints no CCSD(T) energy")
        elif T["E"] is None or abs(T["E"] - T["ccsdt"]) > 1e-9:
            gates.append("G41.4 T's last FINAL SINGLE POINT ENERGY, %s, is not its E(CCSD(T)), %.9f"
                         % ("none" if T["E"] is None else "%.9f" % T["E"], T["ccsdt"]))
        for b in T["inconsistent"]:
            gates.append("G41.4 T: %s" % b)
        if T["atoms"] and T["atoms"] != d38.xyz(geo["TS"]):
            gates.append("G41.1 T RAN on coordinates that are not frame 42436's transition state")
    for lab, r in (("d38/N", N), ("T", T)):
        if r is None:
            continue
        for b in (r["inconsistent"] if lab != "T" else []):
            gates.append("G41.4 %s: %s" % (lab, b))
        for what, p, q, s in (("E(0) + E(CORR)(corrected) = E(CCSD)", r.get("e0"), r.get("corr"), r.get("ccsd")),
                              ("E(CCSD) + (T) = E(CCSD(T))", r.get("ccsd"), r.get("tcorr"), r.get("ccsdt"))):
            if None in (p, q, s):
                gates.append("G41.5 %s cannot check %s: a term is not printed" % (lab, what))
            elif abs(p + q - s) > TOL_ID:
                gates.append("G41.5 %s fails %s by %.2e Eh" % (lab, what, abs(p + q - s)))
    try:
        q = subprocess.run(["git", "-C", a.res, "diff", "--quiet", "HEAD", "--"] + ["%s/%s" % (S, g) for g in GUARDED])
        if q.returncode != 0:
            gates.append("G41.6 a committed file under %s differs from HEAD" % ", ".join(GUARDED))
    except OSError:
        gates.append("G41.6 could not be checked: git did not run")

    if gates:
        say("GATES")
        for g in gates:
            say("   FAILED  %s" % g)
        say("")
        say("ANY FAILURE STOPS THE DIAGNOSTIC. No reading is made.")
        print("\n".join(rep))
        stop("%d gate(s) failed. Nothing is written." % len(gates))

    # ---------------------------------------------- the committed references, each checked
    Z = {leg: d38.read_out(os.path.join(S6, "d36a2", "B_%s" % leg)) for leg in ("R", "TS")}
    pd = os.path.join(S6, D36_RESULT)
    if not os.path.isfile(pd):
        stop("D36's result document %s is missing; R41.2 cannot be read." % D36_RESULT)
    t36 = open(pd, errors="replace").read()
    if any(Z[leg] is None or Z[leg].get(k) is None for leg in ("R", "TS") for _, k in COMP):
        stop("d36a2/B_R or B_TS does not print every component; R41.2 cannot be read.")
    b36 = (Z["TS"]["ccsdt"] - Z["R"]["ccsdt"]) * H2K
    if "%.4f" % b36 != D36_BARRIER or ("THE BARRIER IS %s kcal/mol." % D36_BARRIER) not in t36:
        stop("D36's barrier recomputed from d36a2 is %.4f; it and %s must both be %s."
             % (b36, D36_RESULT, D36_BARRIER))
    p37 = os.path.join(S6, "d37.tsv")
    if not os.path.isfile(p37):
        stop("d37.tsv is missing; R41.2 and R41.3 cannot be read.")
    F = {}
    for l in open(p37).read().splitlines()[1:]:
        f = l.split("\t")
        if len(f) >= 9 and f[5]:
            F[(f[0], f[2])] = dict(er=float(f[3]), ets=float(f[4]), bar=f[5])
    for name, b in SVPD:
        if (name, "def2-SVPD") not in F or (name, "def2-TZVPD") not in F:
            stop("d37.tsv does not carry %s at both basis sets." % name)
        if F[(name, "def2-SVPD")]["bar"] != b:
            stop("the pre-registration lists %s at %s; d37.tsv says %s." % (name, b, F[(name, "def2-SVPD")]["bar"]))
    rows36 = re.findall(r"^   (\S+)\s+(\d+\.\d{4})\s+([+-]\d+\.\d{4})\s*$", t36, re.M)
    T36 = {}
    for name, bar, dist in rows36:
        T36[name] = (bar, dist)         # the eight-row table comes last and so wins
    order36 = []
    for name, bar, dist in rows36[-8:]:
        order36.append(name)
    if sorted(order36) != sorted(n for n, _ in SVPD):
        stop("the last eight rows of D36's table of distances are not the eight functionals.")

    # ------------------------------------------------------------------ R41.1
    say("=" * 78)
    say("GATES G41.1 TO G41.6: ALL PASS")
    say("=" * 78)
    say("")
    say("   G41.1  both geometries are frame 42436's (d33a), and T ran on the transition state's.")
    say("   G41.2  T/job.inp is d38/N/job.inp byte for byte.")
    say("   G41.3  d38/N/job.out is the file recorded in d41/REACTANT_LEG.sha256 and prints")
    say("          E(CCSD(T)) = %.9f." % N_CCSDT)
    say("   G41.4  T terminated normally, converged its SCF, reports %d electrons, and its last" % NEL)
    say("          FINAL SINGLE POINT ENERGY is its E(CCSD(T)).")
    say("   G41.5  in both legs E(0) + E(CORR)(corrected) = E(CCSD) and E(CCSD) + (T) = E(CCSD(T)),")
    say("          each to within %.0e Eh." % TOL_ID)
    say("   G41.6  no committed file under %s differs from HEAD." % ", ".join(GUARDED))
    say("")
    say("=" * 78)
    say("R41.1  THE DLPNO-CCSD(T) BARRIER AT FRAME 42436, def2-SVPD")
    say("=" * 78)
    say("")
    say("   component                 R (d38/N)          TS (d41/T)     barrier kcal/mol")
    part = {}
    for lab, key in COMP:
        x, y = N.get(key), T.get(key)
        if None in (x, y):
            say("   %-24s %16s %17s" % (lab, "?", "?"))
            continue
        part[key] = (y - x) * H2K
        say("   %-24s %16.9f %17.9f %18.4f" % (lab, x, y, part[key]))
    bb = part["ccsdt"]
    say("")
    say("   THE BARRIER AT def2-SVPD IS %.4f kcal/mol." % bb)
    say("")
    say("   Hartree-Fock alone gives %+.4f. Correlation at the CCSD level takes %.4f off it," % (
        part["e0"], part["e0"] - part["ccsd"]))
    say("   leaving %+.4f, and the triples take a further %.4f off, leaving %+.4f." % (
        part["ccsd"], part["ccsd"] - bb, bb))
    say("   D41.1: this is not a converged barrier and must not be quoted as the in vacuo")
    say("   barrier. D41.2: it is not the experimental number.")
    say("")

    # ------------------------------------------------------------------ R41.2
    say("=" * 78)
    say("R41.2  THE METHOD'S OWN BASIS SHIFT, def2-TZVPD minus def2-SVPD")
    say("=" * 78)
    say("")
    say("   component                 def2-TZVPD (D36)   def2-SVPD (here)        shift")
    shift = {}
    for lab, key in COMP:
        if key in ("strong", "weak"):
            continue
        z = (Z["TS"][key] - Z["R"][key]) * H2K
        shift[key] = z - part[key]
        say("   %-24s %16.4f %18.4f %+12.4f" % (lab, z, part[key], shift[key]))
    sh = shift["ccsdt"]
    say("")
    say("   THE CORRELATED BARRIER SHIFTS BY %+.4f kcal/mol from def2-SVPD to def2-TZVPD." % sh)
    say("   That is D36's %s minus R41.1. D36's number is recomputed here from" % D36_BARRIER)
    say("   d36a2/B_{R,TS}/job.out and is the one %s prints." % D36_RESULT)
    say("")
    say("   The eight functionals' shifts, from the energies of d37.tsv:")
    fs = []
    for name, _ in SVPD:
        lo, hi = F[(name, "def2-SVPD")], F[(name, "def2-TZVPD")]
        v = ((hi["ets"] - hi["er"]) - (lo["ets"] - lo["er"])) * H2K
        fs.append(v)
        say("      %-16s %+.4f" % (name, v))
    ext = ("%.4f" % max(fs), "%.4f" % min(fs))
    say("   They run from %s to %s. The pre-registration states %s to %s:" % (ext[0], ext[1], RANGE[0], RANGE[1]))
    say("   %s." % ("reproduced" if ext == RANGE else "NOT REPRODUCED"))
    say("")
    quoted = D36_PHRASE in " ".join(t36.split())
    say("   D36 treated a def2-SVPD barrier as carrying %s%s%s" % ('"' if quoted else "", D36_PHRASE, '"' if quoted else ""))
    say("   on the functionals' evidence. The rule fixed in advance: if the correlated shift")
    say("   differs from -1.0 by more than 1.0 kcal/mol either way, that statement is not")
    say("   supported for the benchmark itself.")
    if abs(sh + 1.0) > 1.0:
        say("   The correlated shift, %+.4f, differs from -1.0 by %.4f. D36'S STATEMENT IS NOT" % (sh, abs(sh + 1.0)))
        say("   SUPPORTED FOR THE BENCHMARK ITSELF.")
    else:
        say("   The correlated shift, %+.4f, differs from -1.0 by %.4f, which is not more than" % (sh, abs(sh + 1.0)))
        say("   1.0. The sentence reserved for a larger difference is therefore not written.")
    say("   D41.5: two basis sets do not make an extrapolation, and none is attempted.")
    say("")

    # ------------------------------------------------------------------ R41.3
    say("=" * 78)
    say("R41.3  THE EIGHT FUNCTIONALS' COMMITTED def2-SVPD BARRIERS MINUS R41.1")
    say("=" * 78)
    say("")
    say("   THE EIGHT BARRIERS WERE KNOWN AND THE BENCHMARK WAS NOT when this comparison was")
    say("   fixed. It is still not the comparison D36 pre-registered, which named three")
    say("   functionals at def2-TZVPD.")
    say("")
    say("   functional         def2-SVPD   minus R41.1     at def2-TZVPD, D36's result")
    dist = {}
    for name, b in SVPD:
        dist[name] = float(b) - bb
        z = T36.get(name)
        chk = ""
        if z is not None:
            mine = "%+.4f" % (float(F[(name, "def2-TZVPD")]["bar"]) - b36)
            if mine != z[1]:
                chk = "   (recomputed here: %s)" % mine
        say("   %-16s %10s %+13.4f     %s%s" % (name, b, dist[name], "%s" % z[1] if z else "not printed", chk))
    order = sorted((n for n, _ in SVPD), key=lambda n: abs(dist[n]))
    say("")
    say("   By distance from the benchmark, nearest first:")
    for lab, o in (("def2-SVPD, here", order), ("def2-TZVPD, the order of D36's table", order36)):
        say("      %s:" % lab)
        for l in textwrap.wrap(", ".join(o), 66):
            say("         %s" % l)
    say("   %s" % ("THE ORDER IS THE SAME at the two basis sets." if order == order36
                   else "THE ORDER DIFFERS between the two basis sets."))
    near = min(THREE, key=lambda n: abs(dist[n]))
    say("   Of the three functionals D36 pre-registered, the nearest at def2-SVPD is %s," % near)
    say("   at %+.4f kcal/mol." % dist[near])
    say("   D41.7: this adopts no production level of theory.")
    say("")

    # ------------------------------------------------------------------ R41.4, R41.5
    say("=" * 78)
    say("R41.4  THE T1 DIAGNOSTIC")
    say("=" * 78)
    say("")
    say("   %-38s %s" % ("d41/T, transition state, def2-SVPD", "%.9f" % T["t1"] if T["t1"] is not None else "not printed"))
    n1 = "%.9f" % N["t1"] if N["t1"] is not None else "not printed"
    say("   %-38s %s%s" % ("d38/N, reactant, def2-SVPD", n1, "" if n1 == D38N_T1 else "   (the pre-registration states %s)" % D38N_T1))
    for leg, ref in D36_T1:
        z1 = "%.9f" % Z[leg]["t1"] if Z[leg]["t1"] is not None else "not printed"
        say("   %-38s %s%s" % ("D36, %s, def2-TZVPD" % ("reactant" if leg == "R" else "transition state"), z1,
                                 "" if z1 == ref else "   (the pre-registration states %s)" % ref))
    say("")
    say("=" * 78)
    say("R41.5  THE HIGHEST OCCUPIED ORBITAL OF THE HARTREE-FOCK REFERENCE")
    say("=" * 78)
    say("")
    say("   By the committed cs10_cpcm.py convention, as cs6_d38.py reads it.")
    say("   d41/T, transition state   %s eV" % ("%+.4f" % T["homo"] if T["homo"] is not None else "not read"))
    nh = "%+.4f" % N["homo"] if N["homo"] is not None else "not read"
    say("   d38/N, reactant           %s eV%s" % (nh, "" if nh == D38N_HOMO else "   (the pre-registration states %s)" % D38N_HOMO))
    say("")

    # ------------------------------------------------------------------ R41.6
    jf = job_facts(os.path.join(D, "run.pbs.out"))
    say("=" * 78)
    say("R41.6  WHERE IT RAN, AND WHAT IT USED")
    say("=" * 78)
    say("")
    say("   node                               %s" % (jf.get("host") or "not reported"))
    say("   scratch location                   %s" % (jf.get("scratch") or "not reported"))
    say("   free when chosen                   %s" % ("%d GB" % jf["avail"] if jf.get("avail") is not None else "not reported"))
    say("   margin                             %s" % (jf.get("margin") or "not reported"))
    say("   maximum sampled while ORCA ran     %s" % ("%d MB, every 60 s; a lower bound on the true peak" % jf["peak"]
                                                     if jf.get("peak") is not None else "not reported"))
    say("   residue after the run              %s" % ("%d MB; not the peak" % jf["residue"]
                                                     if jf.get("residue") is not None else "not reported"))
    if jf.get("refused"):
        say("   candidates refused before it settled:")
        for c, why in jf["refused"]:
            say("      %-40s %s" % (c, why))
    say("")

    say("=" * 78)
    say("WHAT THIS CANNOT SAY (section 5 of the pre-registration, unchanged)")
    say("=" * 78)
    say("")
    for l in ("D41.1  It is not a converged barrier and must not be quoted as the in vacuo barrier.",
              "D41.2  It is not the experimental number; E7.11 applies.",
              "D41.3  One frame, frozen geometries, in vacuo.",
              "D41.4  DLPNO is an approximation; TightPNO's residual error is not measured. The two",
              "       legs were run in different jobs.",
              "D41.5  Two basis sets do not make an extrapolation, and none is attempted.",
              "D41.6  No charges, no design, no reliability claim, nothing about r_safe.",
              "D41.7  It adopts no production level of theory.",
              "D41.8  Nothing here is shaped toward an outcome.",
              "There is no independent control for the transition state at this basis."):
        say("   %s" % l)

    L = ["component\tR_d38N_Eh\tTS_d41T_Eh\tbarrier_SVPD_kcal\tbarrier_TZVPD_kcal\tshift_kcal"]
    for lab, key in COMP:
        if key not in part:
            continue
        z = (Z["TS"][key] - Z["R"][key]) * H2K
        L.append("%s\t%.9f\t%.9f\t%.4f\t%.4f\t%.4f" % (lab, N[key], T[key], part[key], z, z - part[key]))
    open(a.tsv, "w").write("\n".join(L) + "\n")
    open(a.result, "w").write("\n".join(rep) + "\n")
    print("\n".join(rep))
    return 0


if __name__ == "__main__":
    sys.exit(main())
