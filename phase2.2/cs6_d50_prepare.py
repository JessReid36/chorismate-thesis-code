#!/usr/bin/env python3
"""cs6_d50_prepare.py - Stage 6 DIAGNOSTIC 50: functions at a charge that can hold an electron at 0.25 e.

D49 measured that the def2-TZVPD functions of hydrogen, which D47 and D48 placed at each positive charge, do not hold an
electron at a charge of 0.25 e and fall short of a complete basis by 0.37 eV at 0.50 e. D50 adds diffuse functions to
them and measures what the enlarged sets hold, and what the substrate does when they are there.

  PART C  one electron, a bare charge q and the functions: D49's variant N (D48's input with the dipole moment switched
          off), with the added functions written under the ghost atom. Two sets (X2, X3), ten sizes, Hartree-Fock (HF)
          and the functional of the runs (CAM); CAM with three treatments of the auxiliary basis of the ghost atom
          (A: ORCA's AutoAux for that atom; E: explicit even-tempered functions added; J: def2/J left as it is).
  PART G  the dianion with set X2 at the charge, auxiliary treatments A, E and J, at R and TS: +0.25 e at site 1124, and
          the design B_n11_q0.25 with the functions at its six positive charges; each with a control that holds the
          same functions and no charge, and each job with the plain committed control first.

Every condition, ghost centre, input line and the estimate S come from the committed cs6_d48_prepare.py, which this
script loads by path and whose sha256 it checks. It reads committed files only and edits none.

usage: python3 cs6_d50_prepare.py <results repository> <code repository> <output directory> [<directory of d50 on hpc1>]
"""
import hashlib
import importlib.util
import os
import sys

S6 = "phase2.2/cs_stage6"
P48 = "cs6_d48_prepare.py"
P48_SHA = "c8f88c9367e9b45d60e696aaf40f9c239b29ad373b12dbdf4b8acbe950a1e88b"
HPC = "/home/18660916/system_development/phase2.2/cs_stage6/d50"
ELPROP = ["%elprop", "  Dipole false", "end"]             # D49's variant N: the three lines that let ORCA end normally
GHOST = "H: "                                              # how the committed inputs begin a ghost atom's line
ADD = {"X2": (0.035, 0.012), "X3": (0.035, 0.012, 0.004)}  # exponents added to the ghost atom, an s and a p shell at each
NBF = {"X2": 17, "X3": 21}                                 # 9 of def2-TZVPD hydrogen + 4 for each added exponent
AUX = ("A", "E", "J")                                      # the order in which the rule of R50.3 takes them
AUX_TOP, AUX_RATIO = 0.2, 2.0                              # treatment E: exponents from twice the smallest added one, doubling, up to AUX_TOP
NAUX_J = 11                                                # what def2/J gives a hydrogen atom (D49's outputs)
SIZES = (10, 15, 20, 25, 30, 40, 50, 60, 75, 100)          # hundredths of e
LEVELS = ("HF", "CAM")
TOL_ADEQ = 0.05                                            # eV: a set is ADEQUATE at a size if it holds within this of 13.6057 q^2
G_SET, G_AUX = "X2", ("A", "E", "J")                       # part G: the set and the auxiliary treatments run on the dianion
SITE, Q = 1124, 0.25
A_NAME, B_NAME = "A_s1124_p025", "B_n11_q0.25"
JOBS = [("K1", "cm_d50_k1"), ("G1", "cm_d50_g1"), ("G2", "cm_d50_g2")]
NB_PLAIN = 682                                             # basis functions of the plain committed control (d44/CAMTZ_C_ctl)
# The energy of one electron with the charge AT the centre of the functions, hartree, computed for the pre-registration
# from the s functions alone (exact for one electron; the same calculation gives D49's ten variant-Z energies to
# 1e-10 hartree). What Claude expects of the HF runs; the reader sets the measured energies beside these.
EXPECTED = {
    "X2": {10: -0.0001427555, 15: -0.0090328610, 20: -0.0191908593, 25: -0.0310372767, 30: -0.0449661758, 40: -0.0799939706,
           50: -0.1249880738, 60: -0.1799750261, 75: -0.2811960033, 100: -0.4998217239},
    "X3": {10: -0.0045217594, 15: -0.0111999407, 20: -0.0199985525, 25: -0.0312481639, 30: -0.0449980159, 40: -0.0799944974,
           50: -0.1249882080, 60: -0.1799754827, 75: -0.2811963446, 100: -0.4998227575},
}
EXTRA = ["STAGE6_DIAGNOSTIC49_20261121.txt", "STAGE6_DIAGNOSTIC49_RESULT_20261122.txt", "d49_calibration.tsv", "d49_sprime.tsv", "d47.tsv", "d44.tsv",
         "STAGE6_DIAGNOSTIC48_RESULT_20261121.txt", "d48/CONDITIONS.tsv", "d48/run_E1.pbs", "d48/run_A1.pbs",
         "d49/ONE_HF_p025_N/job.inp", "d49/ONE_CAM_p025_N/job.inp", "d49/ONE_CAM_p025_N/pc.pc"]
COMPARE = ["d48/CAMTZG_G_s1124_di_%s", "d48/CAMTZG_G_n11_q0.25_di_%s", "d48/CAMTZG_B_n11_q0.25_di_%s", "d47/CAMTZ_A_s1124_p025_di_%s",
           "d44/CAMTZ_B_n11_q0.25_di_%s"]


class Refused(Exception):
    pass


def stop(m):
    raise Refused(m)


def sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


def sha_text(t):
    return hashlib.sha256(t.encode()).hexdigest()


def load48(code):
    p = os.path.join(code, "phase2.2", P48)
    if not os.path.isfile(p):
        stop("%s is not in the code repository given." % P48)
    if sha(p) != P48_SHA:
        stop("%s is not the committed script this one was written for." % P48)
    spec = importlib.util.spec_from_file_location("d50_cs6_d48_prepare", p)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def aux_exponents(name):
    e, out = 2.0 * min(ADD[name]), []
    while e <= AUX_TOP + 1e-12:
        out.append(e)
        e *= AUX_RATIO
    return out


def naux(name, aux):
    """the auxiliary functions a ghost hydrogen atom holds: known for J and E, not known in advance for A"""
    if aux == "J":
        return NAUX_J
    if aux == "E":
        return NAUX_J + 9 * len(aux_exponents(name))        # an s, a p and a d shell at each exponent
    return None


def ghost_block(name, aux):
    """the lines written under a ghost atom's line: the added functions, then the treatment of its auxiliary basis.
    aux None: a Hartree-Fock run, which uses no auxiliary basis."""
    L = ["  AddGTO"]
    for sh in ("S", "P"):
        for e in ADD[name]:
            L += ["    %s 1" % sh, "    1 %.6f 1.0" % e]
    L.append("  end")
    if aux == "A":
        L += ["  NewAuxJGTO", '    "AutoAux"', "  end"]
    elif aux == "E":
        L.append("  AddAuxJGTO")
        for sh in ("S", "P", "D"):
            for e in aux_exponents(name):
                L += ["    %s 1" % sh, "    1 %.6f 1.0" % e]
        L.append("  end")
    elif aux not in ("J", None):
        stop("no auxiliary treatment %r." % aux)
    return L


def with_block(text, name, aux, nghost):
    """a committed-form input with the block written under each of its ghost atoms"""
    out, n = [], 0
    for l in text.split("\n"):
        out.append(l)
        if l.startswith(GHOST):
            out += ghost_block(name, aux)
            n += 1
    if n != nghost:
        stop("an input holds %d ghost atoms where %d were expected." % (n, nghost))
    return "\n".join(out)


def one_plain(p48, level):
    """D49's variant N: D48's calibration input with the dipole moment switched off"""
    ls = [dict(p48.ONE)[level]] + list(p48.ONE_BODY)
    k = ls.index(p48.PC_LINE)
    return "\n".join(ls[:k] + ELPROP + ls[k:])


def one_input(p48, level, name, aux):
    return with_block(one_plain(p48, level), name, None if level == "HF" else aux, 1)


def select(s6, p48):
    """the two conditions of part G and their controls, from the committed files through cs6_d48_prepare.select"""
    sel = p48.select(s6)
    for lv in LEVELS:
        if one_plain(p48, lv) != open(os.path.join(s6, "d49", "ONE_%s_p025_N" % lv, "job.inp")).read():
            stop("the calibration input formed here without added functions is not d49's committed variant N.")
    cm = dict((x["name"], x) for x in sel["conds"])
    gm = dict((g["name"], g) for g in sel["gctl"])
    site = sel["sites"][SITE]
    old = cm.get("A_s%04d_p050" % SITE)
    b = cm.get(B_NAME)
    if old is None or b is None or old["gctl"] not in gm or b["gctl"] not in gm:
        stop("cs6_d48_prepare.py does not form the conditions this script builds on.")
    if old["pairs"] != [(0.5, SITE)] or old["pc"] != p48.pc_text([(0.5, site["xyz"])]) or gm[old["gctl"]]["ghosts"] != old["ghosts"]:
        stop("D48's condition at site %d is not formed as this script expects." % SITE)
    ch = [(Q, site["xyz"], SITE)]
    a = dict(name=A_NAME, part="A", pairs=[(Q, SITE)], charges=[(Q, site["xyz"])], rows=[SITE], ghosts=list(old["ghosts"]), gctl=old["gctl"],
             pc=p48.pc_text([(Q, site["xyz"])]), where="site %d, the site every design of D43 shares" % SITE, S={}, names={}, vde={}, pred={})
    for st in p48.STATES:
        s, k, vde = p48.s_min(sel["M"], st, ch)
        a["S"][st], a["names"][st], a["vde"][st], a["pred"][st] = s, SITE, vde, p48.predict(s)
    conds = [a, dict(b)]
    for x in conds:
        if any(q > 0 and abs(q - Q) > 1e-12 for q, i in x["pairs"]):
            stop("%s holds a positive charge that is not %.2f e." % (x["name"], Q))
        if len(x["ghosts"]) != sum(1 for q, i in x["pairs"] if q > 0):
            stop("%s does not hold one ghost centre for each positive charge." % x["name"])
    # the estimate formed here is the committed one
    n47 = 0
    for r in p48.table(os.path.join(s6, "d47.tsv")):
        if r.get("condition") == A_NAME and r.get("kind") == "cond" and r.get("S_CAMTZ_eV"):
            if abs(a["S"][r["state"]] - float(r["S_CAMTZ_eV"])) > 6e-5:
                stop("S of %s at %s formed here is not the one d47.tsv holds." % (A_NAME, r["state"]))
            n47 += 1
    n48 = 0
    for r in p48.table(os.path.join(s6, "d48.tsv")):
        if r.get("condition") == B_NAME and r.get("S_eV"):
            if abs(b["S"][r["state"]] - float(r["S_eV"])) > 6e-5:
                stop("S of %s at %s formed here is not the one d48.tsv holds." % (B_NAME, r["state"]))
            n48 += 1
    if n47 < 2 or n48 != 2:
        stop("d47.tsv and d48.tsv do not hold the estimates this script checks itself against.")
    gctl = [dict(gm[x["gctl"]]) for x in conds]
    return dict(conds=conds, gctl=gctl, sel48=sel, G=sel["G"], M=sel["M"], sites=sel["sites"])


def xname(kind, cond, aux, st):
    """CAMTZX: the committed def2-TZVPD level with the enlarged functions at each positive charge"""
    return "CAMTZX_%s_%s%s_di_%s" % (cond, G_SET, aux, st)


def one_name(level, name, aux, k):
    return "ONE_%s_%s%s_p%03d" % (level, name, "" if level == "HF" else aux, k)


def runs_of(sel, p48):
    out = []

    def add(run, job, level, kind, cond="", state="", ref="", gname="", name="", aux="", size=None):
        out.append(dict(run=run, job=job, level=level, kind=kind, cond=cond, state=state, ref=ref, gname=gname, set=name, aux=aux, size=size))

    for name in sorted(ADD):
        for k in SIZES:
            add(one_name("HF", name, None, k), "K1", "ONE_HF", "one", name=name, size=k)
        for aux in AUX:
            for k in SIZES:
                add(one_name("CAM", name, aux, k), "K1", "ONE_CAM", "one", name=name, aux=aux, size=k)
    for job, x in zip(("G1", "G2"), sel["conds"]):
        add("%s_C_ctl_R" % job, job, "CAMTZ", "ctl", state="R", ref="%s_R (committed)" % p48.CTL)
        for aux in G_AUX:
            for st in p48.STATES:
                add(xname("gctl", x["gctl"], aux, st), job, "CAMTZX", "gctl", cond=x["gctl"], state=st, ref="%s_%s (committed)" % (p48.CTL, st),
                    gname=x["gctl"], name=G_SET, aux=aux)
            for st in p48.STATES:
                add(xname("cond", x["name"], aux, st), job, "CAMTZX", "cond", cond=x["name"], state=st, ref=xname("gctl", x["gctl"], aux, st),
                    gname=x["gctl"], name=G_SET, aux=aux)
    if len(set(r["run"] for r in out)) != len(out):
        stop("two runs carry one name.")
    return out


def input_text(s6, sel, p48, r):
    """the text of one job.inp"""
    if r["kind"] == "one":
        return one_input(p48, r["level"].replace("ONE_", ""), r["set"], r["aux"])
    if r["kind"] == "ctl":
        return p48.input_text(s6, sel["sel48"], "ctl")
    cm = dict((x["name"], x) for x in sel["conds"])
    gm = dict((g["name"], g) for g in sel["gctl"])
    if r["kind"] == "gctl":
        t = p48.input_text(s6, sel["sel48"], "gctl", None, r["state"], r["gname"])
        n = len(gm[r["gname"]]["ghosts"])
    else:
        t = p48.input_text(s6, sel["sel48"], "cond", cm[r["cond"]], r["state"])
        n = len(cm[r["cond"]]["ghosts"])
    return with_block(t, r["set"], r["aux"], n)


def pc_of(sel, p48, r):
    if r["kind"] == "one":
        return p48.one_pc(r["size"])
    if r["kind"] == "cond":
        return dict((x["name"], x) for x in sel["conds"])[r["cond"]]["pc"]
    return None


def nb_of(sel, r):
    """the basis functions a measured run must report"""
    if r["kind"] == "one":
        return NBF[r["set"]]
    if r["kind"] == "ctl":
        return NB_PLAIN
    g = [g for g in sel["gctl"] if g["name"] == r["gname"]][0]
    return NB_PLAIN + NBF[r["set"]] * len(g["ghosts"])


def runs_text(runs):
    L = ["run\tjob\tlevel\tkind\tcondition\tstate\treference_control\tset\taux\tsize_e"]
    for r in runs:
        L.append("\t".join([r["run"], r["job"], r["level"], r["kind"], r["cond"], r["state"], r["ref"], r["set"], r["aux"],
                            "" if r["size"] is None else "%.2f" % (r["size"] / 100.0)]))
    return "\n".join(L) + "\n"


def conditions_text(sel):
    L = ["condition\tpart\twhere\tsites\tcharges_e\tghost_control\tghost_centres\tS_R_eV\tS_TS_eV\tS_names_site_R\tS_names_site_TS\tprediction_R\tprediction_TS"]
    for x in sel["conds"]:
        L.append("\t".join([x["name"], x["part"], x["where"], ",".join(str(i) for q, i in x["pairs"]), ",".join("%+.2f" % q for q, i in x["pairs"]), x["gctl"],
                            ";".join("%.6f %.6f %.6f" % g for g in x["ghosts"]), "%+.4f" % x["S"]["R"], "%+.4f" % x["S"]["TS"],
                            str(x["names"]["R"]), str(x["names"]["TS"]), x["pred"]["R"], x["pred"]["TS"]]))
    return "\n".join(L) + "\n"


def functions_text():
    L = ["set\tbasis\ttreatment\tshell\texponent"]
    for name in sorted(ADD):
        for sh in ("S", "P"):
            for e in ADD[name]:
                L.append("%s\torbital\tall\t%s\t%.6f" % (name, sh, e))
        for sh in ("S", "P", "D"):
            for e in aux_exponents(name):
                L.append("%s\tauxiliary\tE\t%s\t%.6f" % (name, sh, e))
    return "\n".join(L) + "\n"


def expected_text(p48):
    L = ["set\tsize_e\tE_Eh_computed_charge_at_the_centre\tB_eV\texact_eV\tshortfall_eV"]
    for name in sorted(ADD):
        for k in SIZES:
            e, q = EXPECTED[name][k], k / 100.0
            L.append("%s\t%.2f\t%.10f\t%.4f\t%.4f\t%.4f" % (name, q, e, -e * p48.EV, p48.RYD * q * q, p48.RYD * q * q + e * p48.EV))
    return "\n".join(L) + "\n"


def sources(sel, p48):
    out = list(p48.sources(sel["sel48"])) + list(EXTRA)
    for st in p48.STATES:
        out += ["%s/job.out" % (c % st) for c in COMPARE]
    seen, uniq = set(), []
    for f in out:
        if f not in seen:
            seen.add(f)
            uniq.append(f)
    return uniq


def pbs_text(s6, job, runs, hpc):
    """the job script: D48's committed run_E1.pbs (one electron, each run stopped after 900 s) or run_A1.pbs (the
    dianion), with its names and its runs changed and nothing else"""
    jk, pbs = job
    old = "E1" if jk == "K1" else "A1"
    mine = [r["run"] for r in runs if r["job"] == jk]
    ls = open(os.path.join(s6, "d48", "run_%s.pbs" % old)).read().split("\n")
    i0 = next(i for i, l in enumerate(ls) if l.startswith("run "))
    i1 = next(i for i, l in enumerate(ls) if l.startswith('echo "end='))
    head, mid, tail = ls[:i0], ls[i0:i1], ls[i1:]
    names = [l[4:] for l in mid]
    loop = "n=0; for r in %s; do grep -q 'ORCA TERMINATED NORMALLY' \"$W/$r/job.out\" 2>/dev/null && n=$((n+1)); done"
    if any(not l.startswith("run ") for l in mid) or len(tail) != 4 or tail[3] != "" or tail[1] != loop % " ".join(names) \
            or tail[2] != 'echo "finished: $n of %d"' % len(names):
        stop("d48/run_%s.pbs is not the job script this script was written for." % old)
    d48 = "/home/18660916/system_development/phase2.2/cs_stage6/d48"
    rep = [("#PBS -N cm_d48_%s" % old.lower(), "#PBS -N %s" % pbs),
           ("#PBS -o %s/run_%s.pbs.out" % (d48, old), "#PBS -o %s/run_%s.pbs.out" % (hpc, jk)),
           ("# Stage 6 DIAGNOSTIC 48: is the estimate S of DIAGNOSTIC 47 a rule?", "# Stage 6 DIAGNOSTIC 50: functions at a charge that can hold an electron at 0.25 e."),
           ("# This job is one of 7. It makes %d runs in the order below; a run that fails does not stop the ones after it." % len(names),
            "# This job is one of %d. It makes %d runs in the order below; a run that fails does not stop the ones after it." % (len(JOBS), len(mine))),
           ("W=%s" % d48, "W=%s" % hpc),
           ('if ! mkdir .running_%s 2>/dev/null; then echo "FAIL: another job holds .running_%s in $W"; exit 1; fi' % (old, old),
            'if ! mkdir .running_%s 2>/dev/null; then echo "FAIL: another job holds .running_%s in $W"; exit 1; fi' % (jk, jk)),
           ("trap 'rmdir \"$W/.running_%s\" 2>/dev/null' EXIT" % old, "trap 'rmdir \"$W/.running_%s\" 2>/dev/null' EXIT" % jk)]
    for a, b in rep:
        if head.count(a) != 1:
            stop("d48/run_%s.pbs does not hold the line: %s" % (old, a))
        head[head.index(a)] = b
    if any(("d48" in l or old in l or "of 7" in l) for l in head):
        stop("a line of d48/run_%s.pbs still names D48 after the names were changed." % old)
    return "\n".join(head + ["run %s" % n for n in mine] + [tail[0], loop % " ".join(mine), 'echo "finished: $n of %d"' % len(mine), ""])


def build(s6, code, hpc=HPC):
    """every file of d50, as {path under d50: text}; the checksum list last"""
    p48 = load48(code)
    try:
        sel = select(s6, p48)
    except p48.Refused as e:
        stop("cs6_d48_prepare.py: %s" % e)
    runs = runs_of(sel, p48)
    geo = dict((st, open(os.path.join(s6, p48.GEOM[st])).read()) for st in p48.STATES)
    F = {}
    for r in runs:
        F["%s/job.inp" % r["run"]] = input_text(s6, sel, p48, r)
        pc = pc_of(sel, p48, r)
        if pc is not None:
            F["%s/pc.pc" % r["run"]] = pc
        if r["kind"] == "ctl":
            F["%s/substrate.xyz" % r["run"]] = geo[r["state"]]
    F["RUNS.tsv"] = runs_text(runs)
    F["CONDITIONS.tsv"] = conditions_text(sel)
    F["FUNCTIONS.tsv"] = functions_text()
    F["EXPECTED_HF.tsv"] = expected_text(p48)
    F["SOURCES.sha256"] = "".join("%s  %s\n" % (sha(os.path.join(s6, f)), f) for f in sources(sel, p48))
    for job in JOBS:
        F["run_%s.pbs" % job[0]] = pbs_text(s6, job, runs, hpc)
    F["INPUT_CHECKSUMS.sha256"] = "".join("%s  %s\n" % (sha_text(F[k]), k) for k in sorted(F))
    return F, runs, sel, p48


def main():
    if len(sys.argv) not in (4, 5):
        sys.stderr.write("usage: python3 cs6_d50_prepare.py <results repository> <code repository> <output directory> [<directory of d50 on hpc1>]\n")
        return 1
    s6 = os.path.join(sys.argv[1], S6)
    code, out = sys.argv[2], sys.argv[3]
    hpc = sys.argv[4] if len(sys.argv) == 5 else HPC
    if os.path.exists(out):
        sys.stderr.write("STOP  %s exists already; this script overwrites nothing.\n" % out)
        return 1
    try:
        F, runs, sel, p48 = build(s6, code, hpc)
    except Refused as e:
        sys.stderr.write("STOP  %s\n" % e)
        return 1
    for k in sorted(F):
        p = os.path.join(out, k)
        os.makedirs(os.path.dirname(p), exist_ok=True)
        with open(p, "w") as f:
            f.write(F[k])
        if k.endswith(".pbs"):
            os.chmod(p, 0o755)
    print("%d runs (%d of one electron, %d of the dianion), %d files checksummed, %d jobs" % (
        len(runs), sum(1 for r in runs if r["kind"] == "one"), sum(1 for r in runs if r["kind"] != "one"), len(F) - 1, len(JOBS)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
