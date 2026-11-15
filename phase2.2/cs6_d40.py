#!/usr/bin/env python3
"""DIAGNOSTIC 40's reader: vertical detachment energies at the density-functional level.

VDE = E(monoanion) - E(dianion), same geometry, same external charges, positive means bound.
That is Probe 2's convention and the conversion is the probes' own (cs11_design_vde.py).
The dianion half of every difference is the committed output of the run the monoanion input
was copied from. It is not run again and it is read here by the same rule as the monoanion.

Energies. A plain or point-charge run: the last FINAL SINGLE POINT ENERGY. A QM/MM run: the
rule of the committed cs6_analyse.py read_run, which is exactly one bare, one (MM) and one
(QM/MM) final energy, bare plus (MM) equal to (QM/MM), and the (QM/MM) total taken.
read_run itself cannot be called on a monoanion because it requires 118 electrons.

Stability. The committed parser cs6_d33b.py stability() is imported by path and used first.
No open-shell output with a stability analysis exists in the record, so its eigenvalue
pattern is untested on one. If it does not read an output, ORCA's own verdict sentences
are read instead, and if those cannot be matched to the analyses the cell says NOT READ
and ORCA's lines are quoted. Nothing is guessed.

TRAH and level shift are detected by the markers of cs6_d28.py and cs6_d33b.py. The
committed dianion's highest occupied orbital is read by cs6_d28.py read_out under D28's
rule, imported by path.

Nothing is typed in except committed reference values, each named with its source.
R40.3 is NOT EVALUATED here: the pre-registration leaves it to a dated addendum.

usage: cs6_d40.py --res <results-repo> --code <code-repo> --prereg <file> --result <file> --tsv <file>
"""
import argparse, hashlib, importlib.util, os, re, subprocess, sys, textwrap

EV = 27.211386245988            # eV per Eh, vde_designs/cs11_design_vde.py
H2K = 627.5094740631
S = "phase2.2/cs_stage6"
NEL_MONO, NEL_DI = 117, 118
CTL_E, TOL_EH = -836.31486034, 1e-6     # G40.4: D28, CAM-B3LYP-D3BJ/def2-SVPD, reactant
P1_VDE, TOL_EV = 0.5013, 0.001          # G40.5: Probe 1's C5, STAGE6_PROBE1_ERRATUM1_20261027.txt
P1_FILE = "STAGE6_PROBE1_ERRATUM1_20261027.txt"
P1_LINE = r"C5\s+def2-SVPD\s+VDE\s+\+0\.5013 eV"
P2_FILE = "vde_probe/PROBE2_RESULT.txt"
P2_SENT = "Any ddE measured at such a condition is not a barrier effect and must not be reported as one."
PREREG = "STAGE6_DIAGNOSTIC40_20261115.txt"
GUARDED = ("d27", "d28", "d31", "d33a", "d33b", "d34", "d37", "qscan")

CONTROL = ("C_CAM_R", "d28/CAM-B3LYP/def2-SVPD/ctl_R")
PART_A = [  # run, committed dianion source, the functional's name in d37.tsv
    ("A_B3LYP", "d27/def2-SVPD/ctl_R", "B3LYP-D3BJ"),
    ("A_CAM", "d28/CAM-B3LYP/def2-SVPD/ctl_R", "CAM-B3LYP-D3BJ"),
    ("A_wB97XD3", "d28/wB97X-D3/def2-SVPD/ctl_R", "wB97X-D3"),
    ("A_wB97MV", "d37/wB97MV_SV_R", "wB97M-V"),
    ("A_wB97XV", "d37/wB97XV_SV_R", "wB97X-V"),
    ("A_PW6B95", "d37/PW6B95_SV_R", "PW6B95-D3BJ"),
    ("A_M062X", "d37/M062X_SV_R", "M06-2X-D3(0)"),
    ("A_PWPB95", "d37/PWPB95_SV_R", "PWPB95-D3BJ"),
]
PART_B = [  # arrangement, source pattern, kind, charge-free pattern, quoted ddE, where it is quoted
    ("bare", "d31/ctl_%s", "plain", None, None, None),
    ("ADD-s540_qm06500", "d34/a2_ADD-s540_qm06500_%s", "qmmm", "d34/a2_S6-C_%s", -3.3985, "D34, R34.9"),
    ("ADD-s678_qm07625", "d34/a2_ADD-s678_qm07625_%s", "qmmm", "d34/a2_S6-C_%s", -4.1567, "D34, R34.9"),
    ("ADD-s678_qm10000", "d34/a2_ADD-s678_qm10000_%s", "qmmm", "d34/a2_S6-C_%s", -5.3701, "D34, R34.9"),
    ("bare_n05_q0.50", "d33b/bare_n05_q0.50_%s", "pc", "d31/ctl_%s", -14.3612, "D33 stage B, with its Erratum 2"),
    ("bare_n11_q0.25", "d33b/bare_n11_q0.25_%s", "pc", "d31/ctl_%s", -13.0007, "D33 stage B, with its Erratum 2"),
]
STATES = ("R", "TS")
FILES = {"plain": ("substrate.xyz",), "pc": ("substrate.xyz", "pc.pc"), "qmmm": ("combined.xyz", "site.ORCAFF.prms")}
GEO = {"plain": "substrate.xyz", "pc": "substrate.xyz", "qmmm": "combined.xyz"}
NQM = 24


def stop(m):
    sys.stderr.write("STOP  %s\n" % m)
    raise SystemExit(1)


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 16), b""):
            h.update(b)
    return h.hexdigest()


def load(path, name):
    sp = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(sp)
    sp.loader.exec_module(m)
    return m


def git(res, *a):
    try:
        q = subprocess.run(["git", "-C", res] + list(a), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except OSError:
        return None, ""
    return q.returncode, q.stdout.decode(errors="replace").strip()


def orca_words(t):
    """ORCA's own words for why a run stopped."""
    L = t.splitlines()
    for i, l in enumerate(L):
        if "INPUT ERROR" in l:
            box = [x.strip() for x in L[i:i + 4] if x.strip() and not set(x.strip()) <= set("!")]
            return " / ".join(box)[:200]
    err = [l.strip() for l in L if re.search(r"error|aborting|not converged|unrecognized", l, re.I)]
    if err:
        return err[-1][:200]
    return "no error line found; the output ends: " + (L[-1].strip()[:120] if L else "(empty)")


def stability(t, d33b):
    """first and last verdict of the stability analyses in one output; never a guess."""
    n = t.count("WAVEFUNCTION STABILITY ANALYSIS")
    sent = [x.upper() for x in re.findall(r"stability analysis[^\n]*?wave ?function is (stable|unstable)", t, re.I)]
    if n == 0 and not sent:
        q = [l.strip()[:150] for l in t.splitlines() if re.search(r"stabil", l, re.I) and not re.search(r"\bsec\b|STABPerform", l)][:4]
        return dict(ran=False, n=0, first=None, last=None, how="no header found", nneg=None, lo=None, quote=q)
    if n == 0:
        return dict(ran=True, n=len(sent), first=sent[0], last=sent[-1], nneg=None, lo=None, quote=[],
                    how="ORCA's verdict sentences; no WAVEFUNCTION STABILITY ANALYSIS header was found")
    c = d33b.stability(t)
    quote = [l.strip()[:150] for l in t.splitlines()
             if re.search(r"stability analysis shows|wavefunction is (un)?stable|UNSTABLE", l) and "WAVEFUNCTION STABILITY ANALYSIS" not in l][:4]
    if c["n"] == n and c["first"] is not None and c["last"] is not None:
        if len(sent) == n and (sent[0] != c["first"] or sent[-1] != c["last"]):
            return dict(ran=True, n=n, first=None, last=None, nneg=None, lo=None, quote=quote,
                        how="NOT READ: cs6_d33b.stability and ORCA's verdict sentences disagree")
        return dict(ran=True, n=n, first=c["first"], last=c["last"], nneg=c["nneg"], lo=c["lo"], quote=quote,
                    how="cs6_d33b.stability")
    if len(sent) == n:
        return dict(ran=True, n=n, first=sent[0], last=sent[-1], nneg=None, lo=None, quote=quote,
                    how="ORCA's verdict sentences; cs6_d33b.stability did not read this output")
    return dict(ran=True, n=n, first=None, last=None, nneg=None, lo=None, quote=quote,
                how="NOT READ: %d analyses and %d verdict sentences" % (n, len(sent)))


def read(d, kind, nel, d33b):
    """One output, by the rule of its kind. 'why' lists every reason it is not readable."""
    p = os.path.join(d, "job.out")
    if not os.path.isfile(p):
        return dict(exists=False, why=["produced no output"], E=None)
    t = open(p, errors="replace").read()
    why = []
    r = dict(exists=True, path=p, text=t, E=None, e_mm=None)
    r["normal"] = "****ORCA TERMINATED NORMALLY****" in t
    if not r["normal"]:
        why.append("did not terminate normally: %s" % orca_words(t))
    r["conv"] = bool(re.search(r"SCF CONVERGED AFTER\s+\d+\s+CYCLES", t)) and "SCF NOT CONVERGED" not in t
    if r["normal"] and not r["conv"]:
        why.append("did not converge its SCF")
    ne = sorted(set(int(x) for x in re.findall(r"Number of Electrons\s+NEL\s*\.*\s*(\d+)", t)))
    r["nel"] = ne
    if ne != [nel]:
        why.append("reports %s electrons, not %d" % (", ".join(map(str, ne)) if ne else "no count of", nel))
    qm = re.findall(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", t)
    mm = re.findall(r"FINAL SINGLE POINT ENERGY \(MM\)\s+(-?\d+\.\d+)", t)
    tot = re.findall(r"FINAL SINGLE POINT ENERGY \(QM/MM\)\s+(-?\d+\.\d+)", t)
    if kind == "qmmm":
        if not (len(qm) == len(mm) == len(tot) == 1):
            why.append("%d bare, %d (MM) and %d (QM/MM) final energies, not one of each" % (len(qm), len(mm), len(tot)))
        elif abs(float(qm[0]) + float(mm[0]) - float(tot[0])) > 1e-9:
            why.append("the QM energy plus the MM energy is not the QM/MM total")
        else:
            r["E"], r["e_mm"] = float(tot[0]), float(mm[0])
        m = re.search(r"Point charges in QM calc\. from MM atoms\.*\s*(\d+)", t)
        nat = None
        try:
            nat = int(open(os.path.join(d, "combined.xyz")).readline().split()[0])
        except (OSError, ValueError, IndexError):
            why.append("combined.xyz could not be read")
        if nat is not None and (not m or int(m.group(1)) != nat - NQM):
            why.append("ORCA read %s MM point charges, the arrangement has %d sphere(s)" % (m.group(1) if m else "no", nat - NQM))
    else:
        if mm or tot:
            why.append("reports MM or QM/MM energies; it is meant to be a plain QM job")
        elif not qm:
            why.append("has no final energy")
        else:
            r["E"] = float(qm[-1])
        m = re.search(r"Reading point charge file\s*\.*\s*ok \((\d+) point charges?\)", t)
        if kind == "pc":
            npc = None
            try:
                npc = int(open(os.path.join(d, "pc.pc")).readline().split()[0])
            except (OSError, ValueError, IndexError):
                why.append("pc.pc could not be read")
            if npc is not None and (not m or int(m.group(1)) != npc):
                why.append("ORCA read %s point charges, pc.pc holds %d" % (m.group(1) if m else "no", npc))
        elif m:
            why.append("read a point-charge file; it is meant to have none")
    s2 = re.findall(r"Expectation value of <S\*\*2>\s*:\s*(-?\d+\.\d+)", t)
    r["s2"] = float(s2[-1]) if s2 else None
    r["trah"] = ("Leaving SCF to start the TRAH-SCF procedure" in t) or ("(TRAH MAcro)" in t)
    r["lshift"] = bool(re.search(r"[Ll]evel ?[Ss]hift set to", t))
    r["stab"] = stability(t, d33b)
    r["why"] = why
    if why:
        r["E"] = None
    return r


def main():
    ap = argparse.ArgumentParser()
    for a in ("res", "code", "prereg", "result", "tsv"):
        ap.add_argument("--" + a, required=True)
    a = ap.parse_args()
    a.res, a.code, a.prereg = (os.path.abspath(x) for x in (a.res, a.code, a.prereg))
    S6 = os.path.join(a.res, S)
    D = os.path.join(S6, "d40")
    p28 = os.path.join(a.code, "phase2.2", "cs6_d28.py")
    p33 = os.path.join(a.code, "phase2.2", "cs6_d33b.py")
    d28 = load(p28, "d40_cs6_d28")
    d33b = load(p33, "d40_cs6_d33b")

    rep, gates = [], []
    say = rep.append
    say("STAGE 6, DIAGNOSTIC 40 - RESULT")
    say("Detachment energies at the density-functional level: eight functionals on the bare")
    say("substrate, and the five arrangements whose lowering is quotable.")
    say("")
    say("Read by phase2.2/cs6_d40.py                sha256 %s" % sha(__file__))
    say("   orbital energies by cs6_d28.py read_out  sha256 %s" % sha(p28))
    say("   stability by cs6_d33b.py stability       sha256 %s" % sha(p33))
    say("   pre-registration %-23s sha256 %s" % (os.path.basename(a.prereg), sha(a.prereg)))
    say("")
    say("VDE = E(monoanion) - E(dianion) at the same geometry with the same external charges,")
    say("in eV at %.12f eV per Eh. Positive means bound. Frame 42436, def2-SVPD." % EV)
    say("The dianion energies are the committed outputs of the runs the monoanion inputs were")
    say("copied from. They are not run again.")
    say("")

    # ------------------------------------------------------------ every run and its source
    runs = [(CONTROL[0], CONTROL[1], "plain", "R")]
    runs += [(r, s, "plain", "R") for r, s, _ in PART_A]
    runs += [("B_%s_%s" % (n, st), pat % st, kind, st) for n, pat, kind, _, _, _ in PART_B for st in STATES]
    if len(runs) != 21:
        stop("the reader's own list holds %d runs, not 21" % len(runs))

    # ------------------------------------------------------------ G40.1 to G40.3
    nfile = 0
    for run, src, kind, st in runs:
        for f in FILES[kind]:
            x, y = os.path.join(D, run, f), os.path.join(S6, src, f)
            if not (os.path.isfile(x) and os.path.isfile(y)) or sha(x) != sha(y):
                gates.append("G40.1 %s/%s is not byte-identical to %s/%s" % (run, f, src, f))
            nfile += 1
        if kind in ("plain", "pc"):
            x, y = os.path.join(D, run, "substrate.xyz"), os.path.join(S6, "d33a", "sp_42436_%s" % st, "substrate.xyz")
            if not (os.path.isfile(x) and os.path.isfile(y)) or sha(x) != sha(y):
                gates.append("G40.1 %s/substrate.xyz is not d33a/sp_42436_%s/substrate.xyz" % (run, st))
        pi, ps = os.path.join(D, run, "job.inp"), os.path.join(S6, src, "job.inp")
        if not (os.path.isfile(pi) and os.path.isfile(ps)):
            gates.append("G40.2 %s or its source has no input" % run)
            continue
        if run == CONTROL[0]:
            if sha(pi) != sha(ps):
                gates.append("G40.2 the control input is not byte-identical to %s/job.inp" % src)
            continue
        li, ls = open(pi).read().split("\n"), open(ps).read().split("\n")
        diff = [i for i in range(min(len(li), len(ls))) if li[i] != ls[i]]
        if len(li) != len(ls) or len(diff) != 1 \
                or ls[diff[0]] != "* xyzfile -2 1 %s" % GEO[kind] or li[diff[0]] != "* xyzfile -1 2 %s" % GEO[kind]:
            gates.append("G40.2 %s does not differ from %s/job.inp in the charge-and-multiplicity line alone" % (run, src))
    x, y = os.path.join(D, "A_B3LYP", "job.inp"), os.path.join(S6, "qscan", "VDE", "C5_mono_R", "job.inp")
    if not (os.path.isfile(x) and os.path.isfile(y)) or sha(x) != sha(y):
        gates.append("G40.3 A_B3LYP/job.inp is not byte-identical to qscan/VDE/C5_mono_R/job.inp")

    # ------------------------------------------------------------ G40.6
    rc, _ = git(a.res, "diff", "--quiet", "HEAD", "--", *["%s/%s" % (S, x) for x in GUARDED])
    if rc is None:
        gates.append("G40.6 could not be checked: git did not run")
    elif rc != 0:
        gates.append("G40.6 a committed file under %s differs from HEAD" % ", ".join(GUARDED))
    rc, pc = git(a.res, "log", "--diff-filter=A", "--format=%H", "-1", "--", "%s/%s" % (S, PREREG))
    dian = sorted(set([s for _, s, _, _ in runs] + [cpat % st for _, _, _, cpat, _, _ in PART_B if cpat for st in STATES]))
    if rc != 0 or not pc:
        gates.append("G40.6 the commit that added %s could not be found" % PREREG)
    else:
        rc1, b1 = git(a.res, "rev-parse", "%s:%s/%s" % (pc, S, PREREG))
        rc2, b2 = git(a.res, "hash-object", a.prereg)
        if rc1 != 0 or rc2 != 0 or b1 != b2:
            gates.append("G40.6 the pre-registration given is not the committed %s" % PREREG)
        for s in dian:
            rc1, b1 = git(a.res, "rev-parse", "%s:%s/%s/job.out" % (pc, S, s))
            rc2, b2 = git(a.res, "hash-object", os.path.join(S6, s, "job.out"))
            if rc1 != 0 or rc2 != 0 or b1 != b2:
                gates.append("G40.6 %s/job.out is not the file at the commit the pre-registration is committed on" % s)

    # ------------------------------------------------------------ the committed references
    pp = os.path.join(S6, P1_FILE)
    if not os.path.isfile(pp) or not re.search(P1_LINE, open(pp).read()):
        gates.append("G40.5 the committed source of +0.5013 eV, %s, does not print it" % P1_FILE)
    pp = os.path.join(S6, P2_FILE)
    if not os.path.isfile(pp) or P2_SENT not in " ".join(open(pp).read().split()):
        gates.append("R40.7 Probe 2's sentence is not in %s" % P2_FILE)

    # ------------------------------------------------------------ read everything
    DI = {}
    for s in dian:
        kind = "qmmm" if s.startswith("d34/") and "S6-C" not in s else ("pc" if s.startswith("d33b/") else "plain")
        r = read(os.path.join(S6, s), kind, NEL_DI, d33b)
        if r["E"] is None:
            gates.append("G40.6 the committed dianion output %s is not readable: %s" % (s, "; ".join(r["why"])))
        DI[s] = r
    MO = {}
    for run, src, kind, st in runs:
        MO[run] = read(os.path.join(D, run), kind, NEL_DI if run == CONTROL[0] else NEL_MONO, d33b)

    def vde(run, src):
        m, d = MO[run], DI[src]
        if m["E"] is None or d["E"] is None:
            return None
        return (m["E"] - d["E"]) * EV

    # ------------------------------------------------------------ G40.4 and G40.5
    c = MO[CONTROL[0]]
    if c["E"] is None:
        gates.append("G40.4 the control %s is not readable: %s" % (CONTROL[0], "; ".join(c["why"])))
    elif abs(c["E"] - CTL_E) > TOL_EH:
        gates.append("G40.4 the control gives %.8f Eh, which differs from D28's %.8f by %.2e, over the %.0e limit"
                     % (c["E"], CTL_E, abs(c["E"] - CTL_E), TOL_EH))
    v0 = vde("A_B3LYP", PART_A[0][1])
    if v0 is None:
        if MO["A_B3LYP"]["E"] is None:
            gates.append("G40.5 A_B3LYP is not readable: %s" % "; ".join(MO["A_B3LYP"]["why"]))
    elif abs(v0 - P1_VDE) > TOL_EV:
        gates.append("G40.5 A_B3LYP's detachment energy is %+.4f eV, not %+.4f within %.3f" % (v0, P1_VDE, TOL_EV))

    if gates:
        say("GATES")
        for g in gates:
            say("   FAILED  %s" % g)
        say("")
        say("A FAILURE OF G40.1 TO G40.6 STOPS THE DIAGNOSTIC. Nothing else is readable.")
        print("\n".join(rep))
        stop("%d gate(s) failed. Nothing is written." % len(gates))

    say("=" * 78)
    say("G40.1 TO G40.6  THE STOPPING GATES")
    say("=" * 78)
    say("")
    say("   G40.1  %d geometry, point-charge and force-field files are byte-identical to the" % nfile)
    say("          committed runs' own, and every bare geometry to d33a/sp_42436_{R,TS}.")
    say("   G40.2  each of the 20 monoanion inputs differs from its committed source in the")
    say("          charge-and-multiplicity line alone; the control input in nothing.")
    say("   G40.3  A_B3LYP/job.inp is byte-identical to qscan/VDE/C5_mono_R/job.inp.")
    say("   G40.6  no committed file under %s" % ", ".join(GUARDED[:-1]))
    say("          or %s differs from HEAD, and the %d committed dianion outputs read are the" % (GUARDED[-1], len(dian)))
    say("          files at the commit the pre-registration was committed on, %s." % pc[:8])
    say("   G40.4 and G40.5 are the two controls, under R40.1.")
    say("")

    # ------------------------------------------------------------ R40.1
    say("=" * 78)
    say("R40.1  THE TWO CONTROLS")
    say("=" * 78)
    say("")
    say("   G40.4  C_CAM_R, D28's own CAM-B3LYP-D3BJ/def2-SVPD dianion input, run again:")
    say("          %.8f Eh   D28 %.8f   differ %.2e   limit %.0e   PASS" % (c["E"], CTL_E, abs(c["E"] - CTL_E), TOL_EH))
    say("   G40.5  A_B3LYP, Probe 1's committed monoanion input, against the committed dianion")
    say("          %s:" % PART_A[0][1])
    say("          VDE %+.4f eV   Probe 1 %+.4f   differ %.4f   limit %.3f   PASS" % (v0, P1_VDE, abs(v0 - P1_VDE), TOL_EV))
    say("          Probe 1's value is the C5 line of %s." % P1_FILE)
    say("")

    # ------------------------------------------------------------ R40.2
    def fs(v, fmt="%.6f"):
        return "-" if v is None else fmt % v

    say("=" * 78)
    say("R40.2  PART A. THE BARE SUBSTRATE AT THE REACTANT, ONE RUN PER FUNCTIONAL")
    say("=" * 78)
    say("")
    t37 = {}
    p37 = os.path.join(S6, "d37.tsv")
    if os.path.isfile(p37):
        for l in open(p37).read().splitlines()[1:]:
            f = l.split("\t")
            if len(f) > 6 and f[2] == "def2-SVPD" and f[6]:
                t37[f[0]] = f[6]
    say("   functional      E(monoanion) Eh  E(dianion) Eh     VDE eV   <S**2>   HOMO eV   HOMO+VDE eV")
    A, notes = {}, []
    for run, src, name in PART_A:
        m, d = MO[run], DI[src]
        v = vde(run, src)
        o = d28.read_out(os.path.join(S6, src, "job.out"))
        h = o.get("homo")
        A[run] = dict(v=v, h=h, m=m, d=d, name=name, src=src)
        if h is None:
            notes.append("%s: the committed dianion output yields no quotable orbital energy under D28's rule." % name)
        elif name not in t37:
            notes.append("%s: d37.tsv has no def2-SVPD orbital energy to compare with." % name)
        elif t37[name] != "%+.4f" % h:
            notes.append("%s: the orbital energy read here, %+.4f, is not d37.tsv's %s." % (name, h, t37[name]))
        if v is None:
            say("   %-15s NOT MEASURED (G40.7)   E(dianion) %.8f" % (name, d["E"]))
            for w in m["why"]:
                say("        %s" % w)
            continue
        say("   %-15s %15.8f %15.8f   %+8.4f   %s   %s   %s"
            % (name, m["E"], d["E"], v, fs(m["s2"]), "%+.4f" % h if h is not None else "      -",
               "%+.4f" % (h + v) if h is not None else "-"))
    say("")
    say("   HOMO is the committed dianion's highest occupied orbital at def2-SVPD and the")
    say("   reactant, read under D28's rule. HOMO+VDE is zero for a functional whose orbital")
    say("   energy equals minus its own detachment energy.")
    if notes:
        for n in notes:
            say("   NOTE  %s" % n)
    else:
        say("   Each orbital energy read here is the one d37.tsv carries.")
    say("   D40.3: these monoanions are not analysed for stability, so each value is an upper")
    say("   bound in that sense. <S**2> is reported and not corrected; the ideal value is 0.75.")
    say("")

    # ------------------------------------------------------------ R40.3
    say("=" * 78)
    say("R40.3  PART A AGAINST THE BENCHMARK")
    say("=" * 78)
    say("")
    say("   NOT EVALUATED HERE. The pre-registration fixes that the comparison with D39's")
    say("   coupled-cluster value is made by a dated addendum to D40 once D39's result is")
    say("   committed, from the number that result prints.")
    say("")

    # ------------------------------------------------------------ R40.4
    say("=" * 78)
    say("R40.4  PART B. CAM-B3LYP-D3BJ/def2-SVPD WITH THE STABILITY ANALYSIS, BOTH STATES")
    say("=" * 78)
    say("")
    say("   BOUND at this condition: VDE above zero. METASTABLE at this condition: VDE at or")
    say("   below zero. Stability is the monoanion's, first analysis then last.")
    say("")
    say("   arrangement       st    VDE eV  word        <S**2>    stability           TRAH  level shift")
    B, tail = {}, []
    for n, pat, kind, _, _, _ in PART_B:
        for st in STATES:
            run, src = "B_%s_%s" % (n, st), pat % st
            m, d = MO[run], DI[src]
            v = vde(run, src)
            why = list(m["why"])
            if v is not None and kind == "qmmm" and abs(m["e_mm"] - d["e_mm"]) > 1e-9:
                why.append("its Lennard-Jones term, %.12f Eh, is not the committed dianion's %.12f" % (m["e_mm"], d["e_mm"]))
                v = None
            word = None if v is None else ("BOUND" if v > 0 else "METASTABLE")
            B[(n, st)] = dict(v=v, word=word, m=m, d=d, why=why, src=src, run=run)
            if v is None:
                say("   %-17s %-2s  NOT MEASURED (G40.7)   E(dianion) %.8f" % (n, st, d["E"]))
                for w in why:
                    say("        %s" % w)
                continue
            sb = m["stab"]
            cell = "NOT FOUND" if not sb["ran"] else ("NOT READ" if sb["first"] is None else "%s/%s" % (sb["first"], sb["last"]))
            say("   %-17s %-2s  %+8.4f  %-10s  %s  %-18s  %-4s  %s"
                % (n, st, v, word, fs(m["s2"]), cell, "yes" if m["trah"] else "no", "yes" if m["lshift"] else "no"))
            if not sb["ran"]:
                tail.append("%s %s: no stability analysis was found in the output, though its input asks for one. Lines of the output that mention stability: %s" % (
                    n, st, " | ".join(sb["quote"]) or "none"))
            elif sb["first"] is None:
                tail.append("%s %s: %s. ORCA's own lines: %s" % (n, st, sb["how"], " | ".join(sb["quote"]) or "(none found)"))
            else:
                if sb["how"] != "cs6_d33b.stability":
                    tail.append("%s %s: read from %s." % (n, st, sb["how"]))
                if sb["first"] == "UNSTABLE":
                    tail.append("%s %s: the first solution was unstable%s; %d analyses were made." % (
                        n, st, "" if sb["nneg"] is None else " (%d negative, lowest %+.6f Eh)" % (sb["nneg"], sb["lo"]), sb["n"]))
    say("")
    say("   Both energies of every row are in d40.tsv.")
    for l in tail:
        for i, w in enumerate(textwrap.wrap(l, 80)):
            say("   %s%s" % ("NOTE  " if i == 0 else "      ", w))
    say("")
    say("   D40.4: where a detachment energy is negative, its sign is the reading and its size")
    say("   is not. D40.7: the transition-state values are the first at that geometry in the")
    say("   record at any level, and have nothing to reproduce.")
    say("")

    # ------------------------------------------------------------ R40.5
    say("=" * 78)
    say("R40.5  A_CAM AGAINST B_bare_R: one species, one level, without and with the analysis")
    say("=" * 78)
    say("")
    x, y = MO["A_CAM"], MO["B_bare_R"]
    if x["E"] is None or y["E"] is None:
        say("   NOT EVALUATED: both monoanion energies are needed, and %s NOT MEASURED."
            % " and ".join(k for k, r in (("A_CAM", x), ("B_bare_R", y)) if r["E"] is None))
    else:
        say("   A_CAM     %.8f Eh   (no stability analysis)" % x["E"])
        say("   B_bare_R  %.8f Eh   (stability analysis: %s)" % (
            y["E"], "NOT FOUND" if not y["stab"]["ran"] else ("NOT READ" if y["stab"]["first"] is None else "%s/%s" % (y["stab"]["first"], y["stab"]["last"]))))
        say("   B_bare_R minus A_CAM  %+.8f Eh, which is %+.4f eV" % (y["E"] - x["E"], (y["E"] - x["E"]) * EV))
    say("")

    # ------------------------------------------------------------ R40.6
    say("=" * 78)
    say("R40.6  CHANGE IN DETACHMENT ENERGY FROM THE BARE SUBSTRATE AT THE SAME STATE")
    say("=" * 78)
    say("")
    say("   arrangement              R eV        TS eV      Reported, not interpreted.")
    for n, _, _, _, _, _ in PART_B[1:]:
        cells = []
        for st in STATES:
            v, b = B[(n, st)]["v"], B[("bare", st)]["v"]
            cells.append("%+10.4f" % (v - b) if v is not None and b is not None else "%10s" % "-")
        say("   %-17s %s  %s" % (n, cells[0], cells[1]))
    if any(B[("bare", st)]["v"] is None for st in STATES):
        say("   A dash means that state, or the bare substrate at that state, is NOT MEASURED.")
    say("")

    # ------------------------------------------------------------ R40.7
    say("=" * 78)
    say("R40.7  WHAT THE RULE OF THE RECORD THEN REACHES")
    say("=" * 78)
    say("")
    say("   Probe 2, %s:" % P2_FILE)
    for l in textwrap.wrap('"%s"' % P2_SENT, 72):
        say("      %s" % l)
    say("   A lowering is a difference between a transition state and a reactant, so the")
    say("   sentence applies to an arrangement's ddE if either of its states is METASTABLE.")
    say("")
    for n, pat, kind, cpat, quoted, where in PART_B[1:]:
        e = {st: DI[pat % st]["E"] for st in STATES}
        c0 = {st: DI[cpat % st]["E"] for st in STATES}
        dde = ((e["TS"] - c0["TS"]) - (e["R"] - c0["R"])) * H2K
        agree = "recomputed here from the committed outputs: %+.4f, %s" % (
            dde, "which agrees" if "%.4f" % dde == "%.4f" % quoted else "WHICH DOES NOT AGREE")
        w = {st: B[(n, st)]["word"] for st in STATES}
        say("   %-17s quoted lowering %+.4f kcal/mol (%s)" % (n, quoted, where))
        say("   %-17s %s" % ("", agree))
        if None in w.values():
            miss = [st for st in STATES if w[st] is None]
            have = ["%s %s, %+.4f eV" % (st, w[st], B[(n, st)]["v"]) for st in STATES if w[st] is not None]
            say("   %-17s NOT CLASSIFIED (G40.8): %s NOT MEASURED." % ("", " and ".join(miss)))
            if have:
                say("   %-17s Measured: %s." % ("", "; ".join(have)))
        elif "METASTABLE" in w.values():
            say("   %-17s R %s, TS %s." % ("", w["R"], w["TS"]))
            say("   %-17s Probe 2's sentence applies to this arrangement's ddE." % "")
            say("   %-17s THE NUMBER IT REACHES IS %+.4f kcal/mol." % ("", quoted))
        else:
            say("   %-17s R BOUND, TS BOUND." % "")
            say("   %-17s This objection does not reach %+.4f kcal/mol at this level and basis." % ("", quoted))
        say("")
    say("   D40 does no more than that. It withdraws nothing and changes no verdict (D40.1):")
    say("   what becomes of a number that R40.7 reaches is a dated erratum and a decision for")
    say("   Jessica.")
    say("")

    # ------------------------------------------------------------ limits
    say("=" * 78)
    say("WHAT THIS CANNOT SAY (section 6 of the pre-registration, unchanged)")
    say("=" * 78)
    say("")
    for l in ("D40.1  It withdraws nothing and changes no verdict.",
              "D40.2  It does not make binding a leg of the reliability criterion.",
              "D40.3  Unrestricted Kohn-Sham differences; part A is not analysed for stability.",
              "D40.4  A METASTABLE species has no bound state to converge to: sign, not size.",
              "D40.5  One basis, def2-SVPD.",
              "D40.6  Vertical, not adiabatic.",
              "D40.7  One frame.",
              "D40.8  All five arrangements sit inside r_safe. D40 says nothing about that.",
              "D40.9  It does not test the functionals against the benchmark: R40.3 waits for D39.",
              "D40.10 It adopts no level of theory.",
              "D40.11 Nothing here is shaped: METASTABLE is reported as plainly as BOUND, in the",
              "       same table, by a rule fixed before any output existed."):
        say("   %s" % l)
    say("")
    nm = sum(1 for x in A.values() if x["v"] is not None) + sum(1 for b in B.values() if b["v"] is not None)
    say("%d of the 20 monoanion runs are measured.%s" % (
        nm, "" if nm == 20 else " %d NOT MEASURED, each listed above with ORCA's own words." % (20 - nm)))

    # ------------------------------------------------------------ the table
    L = ["part\trun\tstate\tsource\tstatus\tE_mono_Eh\tE_di_Eh\tVDE_eV\tS2\tword\tstab_first\tstab_last\tTRAH\tlevel_shift\thomo_di_eV"]

    def row(part, run, st, src, m, d, v, word, h):
        sb = m.get("stab") or {}
        if v is None:
            return "%s\t%s\t%s\t%s\tNOT MEASURED\t\t%.8f\t\t\t\t\t\t\t\t" % (part, run, st, src, d["E"])
        if part == "A":         # R40.2 fixes no stability, TRAH or level-shift reading for part A
            return "%s\t%s\t%s\t%s\tmeasured\t%.8f\t%.8f\t%+.4f\t%s\t\t\t\t\t\t%s" % (
                part, run, st, src, m["E"], d["E"], v, "" if m["s2"] is None else "%.6f" % m["s2"], "" if h is None else "%+.4f" % h)
        return "%s\t%s\t%s\t%s\tmeasured\t%.8f\t%.8f\t%+.4f\t%s\t%s\t%s\t%s\t%s\t%s\t" % (
            part, run, st, src, m["E"], d["E"], v, "" if m["s2"] is None else "%.6f" % m["s2"], word or "",
            sb.get("first") or ("NOT FOUND" if not sb.get("ran") else "NOT READ"), sb.get("last") or ("NOT FOUND" if not sb.get("ran") else "NOT READ"),
            "yes" if m["trah"] else "no", "yes" if m["lshift"] else "no")

    for run, src, name in PART_A:
        x = A[run]
        L.append(row("A", run, "R", src, x["m"], x["d"], x["v"], None, x["h"]))
    for n, pat, kind, _, _, _ in PART_B:
        for st in STATES:
            b = B[(n, st)]
            L.append(row("B", b["run"], st, b["src"], b["m"], b["d"], b["v"], b["word"], None))
    open(a.tsv, "w").write("\n".join(L) + "\n")
    open(a.result, "w").write("\n".join(rep) + "\n")
    print("\n".join(rep))
    return 0


if __name__ == "__main__":
    sys.exit(main())
