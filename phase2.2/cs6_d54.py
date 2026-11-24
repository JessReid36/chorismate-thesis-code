#!/usr/bin/env python3
"""Stage 6, DIAGNOSTIC 54: the design at the second level, from two starts. THE READER.

Reads the five runs of phase2.2/cs_stage6/d54 under the pre-registration STAGE6_DIAGNOSTIC54_20261124.txt and writes
the result and three tables. It was committed, with the pre-registration and the inputs, before any job was submitted.

Each run is read by the functions of D53's committed reader (cs6_d53.py), which are D50's: the committed reader of
outputs, the committed reading of stability, the committed classifier and the four legs against a control. The
controls, the partners and the runs of D53 that the new runs are held against are D53's committed outputs, read here
again by the same functions.

Nothing is written unless every gate holds. D54 adopts nothing: it names no design and forms no ddE.

usage: python3 cs6_d54.py --res <results repository> --code <code repository> --prereg <the pre-registration>
                          --result <the result to write> --tsv <the table of runs to write>
"""
import argparse
import importlib.util
import os
import re
import subprocess
import sys
import textwrap

S = "phase2.2/cs_stage6"
PREREG = "STAGE6_DIAGNOSTIC54_20261124.txt"
D53PRE = "STAGE6_DIAGNOSTIC53_20261123.txt"
D53RES = "STAGE6_DIAGNOSTIC53_RESULT_20261124.txt"
D52RES = "STAGE6_DIAGNOSTIC52_RESULT_20261123.txt"
D50PRE = "STAGE6_DIAGNOSTIC50_20261122.txt"
START_FRAC = 0.5    # a run STARTED FROM ITS FILE if its first SCF energy lies less than this fraction as far above its
                    # lowest as in D53's run of the same input, which started from ORCA's own guess
E_TOL = 0.0010      # eV: two runs within this of each other in energy ...
LEFT_TOL = 0.02     # e:  ... and within this in the charge that left the 24 atoms are THE SAME SOLUTION
L_PRINTED, L2_PRINTED = "0.3127", "0.1151"
Q_B = "TOOK ELECTRONS the SCF converged and either a leg fires or the run LEFT THE ATOMS"                                              # D50
Q_C = "RULE C GIVES: L = 0.3127 eV."                                                                                                   # D52's result
Q_D = "If the design is among them at either state, the question of this diagnostic is open, and a dated note says what is done"        # D53
Q_E = "THE SECOND LEVEL IS READ."                                                                                                      # D53's result
GUARDED = ["d53"]
DOCS = ["d52_rules.tsv", "d53.tsv", "d53_words.tsv", "d53_controls.tsv", D53PRE, D53RES, D52RES, D50PRE]
# what D53's committed table of words must hold for the design, for this reading to mean what the pre-registration says
D53_WORDS = {("d11", "R"): dict(CAM_did="TOOK ELECTRONS", CAM_bar="WITHIN", WXV_did="not decided", the_two_levels_together="NOT DECIDED"),
             ("d11", "TS"): dict(CAM_did="DID NOT TAKE", CAM_bar="", WXV_did="not decided", the_two_levels_together="NOT DECIDED")}


def stop(m):
    sys.stderr.write("STOP  %s\n" % m)
    raise SystemExit(1)


def load(path, name):
    sp = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(sp)
    sp.loader.exec_module(m)
    return m


def scf_energies(text):
    """the energies of the SCF iteration tables, in the order printed (the lines: iteration, energy, change, ...)"""
    return [float(m) for m in re.findall(r"^ {1,6}\d+ +(-\d+\.\d{6,}) +[-+]?\d\.\d+e[-+]\d+ ", text, re.M)]


def started(x, guess=None):
    """(STARTED FROM ITS FILE or not or None, first SCF energy, lowest SCF energy); guess: how far above its lowest D53's run began"""
    e = scf_energies(x.get("text") or "")
    if not e:
        return None, None, None
    return (None if guess is None else (e[0] - min(e)) < START_FRAC * guess), e[0], min(e)


def same(x, y, ev):
    """THE SAME SOLUTION: both measured, energies within E_TOL, charge that left the atoms within LEFT_TOL"""
    if x["status"] != "measured" or y["status"] != "measured" or x["left"] is None or y["left"] is None:
        return None
    return abs((x["E"] - y["E"]) * ev) <= E_TOL and abs(x["left"] - y["left"]) <= LEFT_TOL


def seeded_word(s, ref, reproduced, d_lw, ev):
    """R54.4: (the word, E(seeded) - E(reference) in eV or None, why it is not read)"""
    if not reproduced:
        return "NOT READ", None, "the way of starting is not REPRODUCED (R54.2)"
    if s["status"] != "measured":
        return "NOT READ", None, "the run %s" % ("DID NOT CONVERGE" if s["status"] == "DID NOT CONVERGE" else "is NOT MEASURED")
    if s.get("noread"):
        return "NOT READ", None, s["noread"]
    if ref["status"] != "measured" or ref["left"] is None or s["left"] is None:
        return "NOT READ", None, "the run it is held against, %s, is not measured or not read" % ref["run"]
    de = (s["E"] - ref["E"]) * ev
    if same(s, ref, ev):
        return "RETURNED", de, ""
    where = "BELOW" if de < -E_TOL else ("ABOVE" if de > E_TOL else "LEVEL")
    return ("A CAPTURE SOLUTION, %s" if s["left"] >= d_lw else "ANOTHER SOLUTION WITHOUT CAPTURE, %s") % where, de, ""


def second_level_R(o, s, sw, counts_ok):
    """R54.5 at R: (the word, the run in which the electrons were taken or None, why not decided).
    o: the run from its own start; s: the run from the capture solution; sw: its word of R54.4"""
    if not counts_ok:
        return "NOT DECIDED", None, "a count does not hold"
    if o["did"] == "TOOK ELECTRONS":
        return "TAKES ELECTRONS", "own", ""
    if sw[0] in ("A CAPTURE SOLUTION, BELOW", "A CAPTURE SOLUTION, LEVEL"):
        return "TAKES ELECTRONS", "seeded", ""
    other_low = sw[0] in ("ANOTHER SOLUTION WITHOUT CAPTURE, BELOW", "ANOTHER SOLUTION WITHOUT CAPTURE, LEVEL")
    if other_low and s["did"] == "TOOK ELECTRONS":
        return "TAKES ELECTRONS", "seeded", ""
    if sw[0] == "NOT READ":
        return "NOT DECIDED", None, "the run from the capture solution is NOT READ: %s" % sw[2]
    if o["did"] != "DID NOT TAKE":
        return "NOT DECIDED", None, "the run from its own start is not decided: %s" % (o["why"] or "not decided")
    if other_low and s["did"] != "DID NOT TAKE":
        return "NOT DECIDED", None, "the run from the capture solution is another solution, not above, and is not decided by the four legs: %s" % (s["why"] or "not decided")
    return "DOES NOT TAKE", None, ""


def main():
    ap = argparse.ArgumentParser()
    for a in ("res", "code", "prereg", "result", "tsv"):
        ap.add_argument("--" + a, required=True)
    a = ap.parse_args()
    a.res, a.code, a.prereg = (os.path.abspath(x) for x in (a.res, a.code, a.prereg))
    S6 = os.path.join(a.res, S)
    D, D53 = os.path.join(S6, "d54"), os.path.join(S6, "d53")
    paths = dict((k, os.path.join(a.code, "phase2.2", fn)) for k, fn in (
        ("d40", "cs6_d40.py"), ("d33b", "cs6_d33b.py"), ("audit", "cs6_audit.py"), ("r53", "cs6_d53.py"), ("p53", "cs6_d53_prepare.py"), ("p54", "cs6_d54_prepare.py")))
    for k, p in paths.items():
        if not os.path.isfile(p):
            stop("the committed script %s is missing." % p)
    for fn in [a.prereg, os.path.join(D, "RUNS.tsv")] + [os.path.join(S6, x) for x in DOCS]:
        if not os.path.isfile(fn):
            stop("%s is missing." % fn)
    tsv2, tsv3 = os.path.join(os.path.dirname(a.tsv), "d54_words.tsv"), os.path.join(os.path.dirname(a.tsv), "d54_starts.tsv")
    for fn in (a.result, a.tsv, tsv2, tsv3):
        if os.path.exists(fn):
            stop("%s exists already; this reader overwrites nothing." % fn)
    d40, d33b, audit = load(paths["d40"], "d54_cs6_d40"), load(paths["d33b"], "d54_cs6_d33b"), load(paths["audit"], "d54_cs6_audit")
    R, p53, p54 = load(paths["r53"], "d54_cs6_d53"), load(paths["p53"], "d54_cs6_d53_prepare"), load(paths["p54"], "d54_cs6_d54_prepare")
    T = (audit.D_LW, audit.R_DIP, audit.GAP_ABS, audit.R_GAP)
    D_LW, EV, NEL = audit.D_LW, R.EV, R.NEL
    gates = []

    def git(*args):
        try:
            q = subprocess.run(["git", "-C", a.res] + list(args), stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
            return q.returncode, q.stdout.strip()
        except OSError:
            return 127, ""

    # ------------------------------------------------------------ G54.3: the committed files this leans on
    guard = ["%s/%s" % (S, g) for g in GUARDED] + ["%s/%s" % (S, x) for x in DOCS]
    rc, st = git("status", "--porcelain", "--untracked-files=all", "--", *guard)
    if rc != 0:
        gates.append("G54.3 could not be checked: git did not run")
    elif st:
        gates.append("G54.3 a file under d53, or a committed table or document this reader names, is not the committed one")
    rc, pcm = git("log", "--diff-filter=A", "--format=%H", "-1", "--", "%s/%s" % (S, PREREG))
    if rc != 0 or not pcm:
        gates.append("G54.3 the commit that added %s could not be found" % PREREG)
    else:
        rc1, b1 = git("rev-parse", "%s:%s/%s" % (pcm, S, PREREG))
        rc2, b2 = git("hash-object", a.prereg)
        if rc1 != 0 or rc2 != 0 or b1 != b2:
            gates.append("G54.3 the pre-registration given is not the committed %s" % PREREG)

    # ------------------------------------------------------------ G54.1 and G54.2: formed again from the committed files
    try:
        F, runs = p54.build(S6, a.code)
    except Exception as e:
        sys.stdout.write("\n".join("FAILED  " + g for g in gates + ["G54.1 the files of d54 could not be formed again from the committed files: %s: %s" % (type(e).__name__, e)]) + "\n")
        stop("%d gate(s) failed. Nothing was written." % (len(gates) + 1))
    psrc = os.path.join(D, "SOURCES.sha256")
    if not os.path.isfile(psrc) or open(psrc).read() != F["SOURCES.sha256"]:
        gates.append("G54.2 d54/SOURCES.sha256 is not the list formed again from the committed files")
    nsrc = 0
    for l in F["SOURCES.sha256"].splitlines():
        h, fn = l.split(None, 1)
        nsrc += 1
        p = os.path.join(S6, fn)
        if not os.path.isfile(p) or R.sha(p) != h:
            gates.append("G54.2 %s is not the file d54/SOURCES.sha256 lists" % fn)
    nsum = 0
    pic = os.path.join(D, "INPUT_CHECKSUMS.sha256")
    if not os.path.isfile(pic):
        gates.append("G54.1 d54/INPUT_CHECKSUMS.sha256 is missing")
    else:
        if open(pic).read() != F["INPUT_CHECKSUMS.sha256"]:
            gates.append("G54.1 d54/INPUT_CHECKSUMS.sha256 is not the list formed again here")
        for l in open(pic).read().splitlines():
            if l.strip():
                h, fn = l.split(None, 1)
                nsum += 1
                p = os.path.join(D, fn)
                if not os.path.isfile(p) or R.sha(p) != h:
                    gates.append("G54.1 d54/%s is not the file the checksum list holds" % fn)
    for k in sorted(F):
        p = os.path.join(D, k)
        if not os.path.isfile(p) or open(p).read() != F[k]:
            gates.append("G54.1 d54/%s is not the file formed again from the committed files" % k)

    # ------------------------------------------------------------ G54.4: the quoted sentences, the bars, D53's words
    def has(path, sent):
        return " ".join(sent.split()) in " ".join(open(path, errors="replace").read().split())
    for doc, q, what, quoted in ((D50PRE, Q_B, "the reading TOOK ELECTRONS", True), (D53PRE, Q_D, "outcome (f) of D53", True),
                                 (D52RES, Q_C, "RULE C of D52", False), (D53RES, Q_E, "that D53 read the second level", False)):
        if not has(os.path.join(S6, doc), q):
            gates.append("G54.4 %s is not in %s" % (what, doc))
        if quoted and not has(a.prereg, q):
            gates.append("G54.4 the words \"%s\" are not in the pre-registration given" % (q[:60] + " ..."))
    rules = dict(l.split("\t")[:2] for l in open(os.path.join(S6, "d52_rules.tsv")).read().split("\n")[1:] if "\t" in l)
    try:
        L = float(rules["rule_C_L_eV"])
        L2 = 3.0 * float(rules["largest_distance_from_the_limit_WXV_eV"])
    except (KeyError, ValueError):
        L = L2 = None
    if L is None or "%.4f" % L != L_PRINTED or "%.4f" % L2 != L2_PRINTED:
        gates.append("G54.4 d52_rules.tsv does not hold an L and an L2 that print as %s and %s" % (L_PRINTED, L2_PRINTED))
    ls = open(os.path.join(S6, "d53_words.tsv")).read().split("\n")
    w53 = dict(((x["arrangement"], x["state"]), x) for x in (dict(zip(ls[0].split("\t"), l.split("\t"))) for l in ls[1:] if l))
    for k, want in D53_WORDS.items():
        if k not in w53 or any(w53[k].get(c) != v for c, v in want.items()):
            gates.append("G54.4 d53_words.tsv does not hold, for the design at %s, the words the pre-registration states" % k[1])
    ls = open(os.path.join(S6, "d53_controls.tsv")).read().split("\n")
    c53 = dict(((x["level"], x["state"]), x) for x in (dict(zip(ls[0].split("\t"), l.split("\t"))) for l in ls[1:] if l))
    for k in (("WXV", "R"), ("WXV", "TS"), ("CAM", "R"), ("HF", "R")):
        if k not in c53 or c53[k].get("with_functions") != "HOLDS" or c53[k].get("plain_control") != "STABLE":
            gates.append("G54.4 d53_controls.tsv does not say that the control of %s at %s HOLDS" % k)

    # ------------------------------------------------------------ reading every run, and G54.5, G54.6, G54.7
    geom = dict((stt, p53.geometry(S6, stt)) for stt in ("R", "TS"))
    RUN = {}
    for r in runs:
        d = os.path.join(D, r["run"])
        x = R.read_di(d40, d33b, audit, d, "pc")
        x.update(r)
        RUN[r["run"]] = x
        if not x.get("exists"):
            continue
        e = R.echoed(x["text"])
        want = open(os.path.join(d, "job.inp")).read().split("\n") if os.path.isfile(os.path.join(d, "job.inp")) else None
        while want and not want[-1].strip():
            want.pop()
        if e is None:
            if x.get("normal"):
                gates.append("G54.5 %s/job.out ends normally and does not print the input it ran" % r["run"])
        elif want is not None and e != want:
            gates.append("G54.5 the input printed at the top of %s/job.out is not the committed %s/job.inp" % (r["run"], r["run"]))
        if x["status"] == "measured" and x["nbf"] != r["nb"]:
            gates.append("G54.6 %s reports %s basis functions; its input asks for %d" % (r["run"], x["nbf"], r["nb"]))
    pman = os.path.join(D, "HPC_MANIFEST.sha256")
    if not os.path.isfile(pman):
        gates.append("G54.7 d54/HPC_MANIFEST.sha256, which hpc1 writes over its own files at the harvest, is missing: the start orbitals cannot be checked")
    else:
        man = dict((l.split(None, 1)[1].strip(), l.split(None, 1)[0]) for l in open(pman).read().split("\n") if l.strip())
        for l in F["START_ORBITALS.sha256"].splitlines():
            h, fn = l.split(None, 1)
            if man.get(fn.strip()) != h:
                gates.append("G54.7 %s on hpc1 is not the file D53's committed manifest lists for it" % fn.strip())
    # the runs of D53 that the new runs are held against
    REF = {}

    def ref(name, kind):
        if name not in REF:
            x = R.read_di(d40, d33b, audit, os.path.join(D53, name), kind)
            x["run"] = name
            REF[name] = x
        return REF[name]

    if gates:
        sys.stdout.write("\n".join("FAILED  " + g for g in gates) + "\n")
        stop("%d gate(s) failed. Nothing was written." % len(gates))

    # ------------------------------------------------------------ each run against its control of D53
    def count_of(x, level):
        if level == "HF":
            return True, None
        if x["status"] != "measured" or x["count"] is None:
            return False, None
        return abs(x["count"] - NEL) <= R.COUNT_TOL, x["count"] - NEL

    def read_against(x, level, state, guess=None):
        gf_n, gn_n, n_n = p54.REFS[level][state]
        gf, gn = ref(gf_n, "plain"), ref(gn_n, "plain")
        ch = p53.charges_of(S6, p53.PC_SITES, state)
        x.setdefault("why_d40", x.get("why"))
        x["started"], x["e_first"], x["e_low"] = started(x, guess)
        x["noread"] = None if (R.stable(gf) and R.stable(gn)) else "a control of D53 is not measured or not STABLE"
        if x.get("kind") in ("O", "S") and x["status"] == "measured" and x["started"] is not True:
            x["noread"] = "it did not start from its start file" if x["started"] is False else "its first SCF energy could not be read"
        x["z"] = R.legs(x, gf, geom[state][0], ch, T) if not x["noread"] else dict(R.NOZ)
        x["word"] = R.word1(x)
        did, left = R.took(x, D_LW)
        x["did"], x["left"] = did, left
        x["how"] = ("ONTO THE FUNCTIONS" if left >= D_LW else "A LEG ONLY") if did == "TOOK ELECTRONS" else ""
        x["why"] = x["noread"] or (None if x["status"] == "measured" else "the run is %s" % x["status"]) or ("; ".join(x["z"]["held"]) if did == "not decided" else None)
        x["count_ok"], x["count_d"] = count_of(x, level)
        x["gain"] = None
        if n_n is not None:
            n = ref(n_n, "pc")
            four = (x, gf, n, gn)
            if all(v["status"] == "measured" for v in four) and x["count_ok"] and all(count_of(v, level)[0] for v in (gf, n, gn)):
                x["gain"] = ((x["E"] - gf["E"]) - (n["E"] - gn["E"])) * EV
        x["ghost"] = (sum(x["lw"][24:30]) - sum(gf["lw"][24:30])) if (x["status"] == "measured" and x["lw"] and gf.get("lw") and len(x["lw"]) >= 30 and not x["noread"]) else None
        return x

    # D53's own runs of the design, read the same way: what the new runs are held against
    OLD = {}
    for lv, stt in (("WXV", "R"), ("WXV", "TS"), ("CAM", "R"), ("HF", "R")):
        x = ref("DI6_%s_F_d11_%s" % (lv, stt), "pc")
        x.update(level=lv, state=stt)
        OLD[(lv, stt)] = read_against(x, lv, stt)
    for r in runs:
        o = OLD[(r["level"], r["state"])]
        if o["run"] != r["source"]:
            stop("the run of D53 with the input of %s is not the one this reader holds it against." % r["run"])
        read_against(RUN[r["run"]], r["level"], r["state"], (o["e_first"] - o["e_low"]) if o.get("e_first") is not None else None)

    o_r, o_ts = RUN["DI7_WXV_O_d11_R"], RUN["DI7_WXV_O_d11_TS"]
    s_w, s_c, s_h = RUN["DI7_WXV_S_d11_R"], RUN["DI7_CAM_S_d11_R"], RUN["DI7_HF_S_d11_R"]
    # R54.2 the way of starting
    rep = same(s_c, OLD[("CAM", "R")], EV)
    reproduced = s_c["status"] == "measured" and s_c["started"] is True and s_c["count_ok"] and s_c["did"] == "TOOK ELECTRONS" and s_c["left"] is not None and s_c["left"] >= D_LW
    # R54.3 the runs from their own start
    own = dict((stt, dict(x=x, same=same(x, OLD[("WXV", stt)], EV))) for stt, x in (("R", o_r), ("TS", o_ts)))
    # R54.4 the run from the capture solution at the second level, and R54.6 at Hartree-Fock
    ref_w = o_r if o_r["status"] == "measured" and o_r["left"] is not None else OLD[("WXV", "R")]
    sw = seeded_word(s_w, ref_w, reproduced, D_LW, EV)
    sh = seeded_word(s_h, OLD[("HF", "R")], reproduced, D_LW, EV)
    # R54.5 the design at the second level
    wR = second_level_R(o_r, s_w, sw, o_r["count_ok"] and (s_w["count_ok"] or sw[0] == "NOT READ"))
    if not o_ts["count_ok"]:
        wTS = ("NOT DECIDED", None, "a count does not hold")
    elif o_ts["did"] == "TOOK ELECTRONS":
        wTS = ("TAKES ELECTRONS", "own", "")
    elif o_ts["did"] == "DID NOT TAKE":
        wTS = ("DOES NOT TAKE", None, "")
    else:
        wTS = ("NOT DECIDED", None, "the run from its own start is not decided: %s" % (o_ts["why"] or "not decided"))

    def both(state, w):
        cam = w53[("d11", state)]["CAM_did"]
        carrier = {"own": o_r if state == "R" else o_ts, "seeded": s_w}.get(w[1])
        if w[0] == "NOT DECIDED":
            return "NOT DECIDED", "", w[2]
        if w[0] == "TAKES ELECTRONS":
            if cam == "TOOK ELECTRONS":
                if carrier["gain"] is None:
                    return "NOT DECIDED", "", "the gain at the second level is not formed"
                return "TAKES ELECTRONS AT BOTH LEVELS", "%s L2 at WXV" % ("ABOVE" if -carrier["gain"] > L2 else "WITHIN"), ""
            return "TAKEN AT THE SECOND LEVEL ONLY", "", ""
        if cam == "TOOK ELECTRONS":
            return "TAKEN AT THE WORKING FUNCTIONAL ONLY", "%s L" % w53[("d11", state)]["CAM_bar"], ""
        return "HOLDS AT BOTH LEVELS", "", ""
    tR, tTS = both("R", wR), both("TS", wTS)

    # ------------------------------------------------------------ the result
    out = []
    say = out.append

    def f(v, fmt="%+.4f", none="-"):
        return none if v is None else fmt % v

    def para(text, lead="", width=96):
        for l in textwrap.wrap(" ".join(text.split()), width=width - len(lead), break_long_words=False, break_on_hyphens=False):
            say(lead + l)

    def item(label, text, ind=10):
        for i, l in enumerate(textwrap.wrap(" ".join(text.split()), width=96 - ind, break_long_words=False, break_on_hyphens=False)):
            say((("   " + label).ljust(ind) if i == 0 else " " * ind) + l)

    def title(t):
        say("")
        say("=" * 96)
        say(t)
        say("=" * 96)
        say("")
    rc, head = git("rev-parse", "HEAD")
    say("STAGE 6 - DIAGNOSTIC 54 - RESULT")
    say("The design at the second level, from two starts")
    say("=" * 96)
    say("")
    para("Read by cs6_d54.py under %s, which was committed before any job was submitted." % PREREG)
    say("   results repository at %s; the pre-registration was added by %s" % (head[:12] if rc == 0 else "?", pcm[:12]))
    for lab, k in (("reader", None), ("inputs formed again", "p54"), ("reading functions", "r53"), ("outputs read by", "d40"), ("stability by", "d33b"), ("the four legs by", "audit")):
        p = os.path.abspath(__file__) if k is None else paths[k]
        say("   %-20s %s" % (lab, os.path.basename(p)))
        say("   %-20s sha256 %s" % ("", R.sha(p)))
    say("")
    para("""The dianion of frame 42436 with the eleven charges of B_n11_q0.25 and the functions of D50 at the six positive sites,
def2-TZVPD, in vacuo. O: the run of D53 at wB97X-V (WXV) again, started from its own converged orbitals. S: the run at R started
from the solution D53 found at CAM-B3LYP-D3BJ (CAM), in which %s e had left the atoms. Every run is read against D53's committed
control with the same functions and no charge. 'left': the charge that left the 24 atoms, e. THE SAME SOLUTION: energies within
%.4f eV and 'left' within %.2f e. The gain: (run minus control) minus (D53's partner without functions minus the plain control),
eV. L = %.4f eV, L2 = %.4f eV. D54 adopts nothing: it names no design and forms no ddE.""" % (f(OLD[("CAM", "R")]["left"]), E_TOL, LEFT_TOL, L, L2))
    say("")
    say("THE GATES, as they held:")
    item("G54.1", "all %d files of d54/INPUT_CHECKSUMS.sha256 are the ones listed, and every file of d54 is the one formed again here from the committed files." % nsum)
    item("G54.2", "the %d files of d54/SOURCES.sha256 are the ones listed." % nsrc)
    item("G54.3", "no file under d53, and no committed table or document named here, differs from the committed one; the pre-registration given is the committed one.")
    item("G54.4", "the sentences quoted here are in the documents they are quoted from; D52's result holds RULE C and D53's result reads the second level; d52_rules.tsv holds L and L2; D53's committed tables hold the words and the controls the pre-registration states.")
    item("G54.5", "every output that exists prints, as the input it ran, the committed input beside it.")
    item("G54.6", "every measured run reports the number of basis functions its input asks for.")
    item("G54.7", "the manifest hpc1 wrote over d54 gives each start.gbw the sha256 that D53's committed manifest gives the file it is a copy of.")

    title("R54.1  THE RUNS")
    para("""'from file': STARTED FROM ITS FILE, as defined under the table. 'shift set': ORCA printed that it set a level shift, which
withholds the gap leg (the committed rule of D50). 'first / last': the stability analysis, first and last.""")
    say("")
    say("   run                 status           count    from file shift set TRAH first / last")
    for r in runs:
        x = RUN[r["run"]]
        say(("   %-19s %-16s %-8s %-9s %-9s %-4s %s" % (x["run"], x["status"], "-" if (x["level"] == "HF" or x["status"] != "measured") else (f(x["count_d"]) if x["count_ok"] else "NO"),
             {True: "yes", False: "NO", None: "-"}[x["started"]], "yes" if x["lshift"] else "no", "yes" if x["trah"] else "no",
             "%s / %s" % (x["stab"].get("first") or "-", x["stab"].get("last") or "-"))).rstrip())
    bad = [RUN[r["run"]] for r in runs if RUN[r["run"]]["status"] != "measured"]
    say("")
    for r in runs:
        x = RUN[r["run"]]
        if x["level"] != "HF" and x["status"] == "measured" and not x["count_ok"]:
            para("%s: the count of electrons does not hold: it is %s from %d." % (x["run"], "not printed" if x["count_d"] is None else ("%+.4f" % x["count_d"]), NEL), "   ")
    if bad:
        say("   Measured: %d of %d. The others, with ORCA's words:" % (len(runs) - len(bad), len(runs)))
        for x in bad:
            para("%s: %s: %s" % (x["run"], x["status"], "; ".join(x.get("why_d40") or ["no output"])), "   ")
    else:
        say("   Every one of the %d runs is measured." % len(runs))
    say("")
    para("""STARTED FROM ITS FILE: the energy of the first SCF step lies less than half as far above the lowest of the run as in
D53's run of the same input, which started from ORCA's own guess (WXV at R %s, WXV at TS %s, CAM %s, HF %s Eh above). A run that
did not start from its file is NOT READ.""" % tuple("%.3f" % (OLD[k]["e_first"] - OLD[k]["e_low"]) if OLD[k].get("e_first") is not None else "-"
                                                    for k in (("WXV", "R"), ("WXV", "TS"), ("CAM", "R"), ("HF", "R"))), "   ")
    for r in runs:
        x = RUN[r["run"]]
        if x["e_first"] is not None:
            src = OLD.get((r["level"], r["state"])) if r["start"] == r["source"] else None
            extra = ""
            if src is not None and src.get("e_low") is not None:
                extra = "; against the lowest of %s, which it starts from: %+.6f Eh" % (src["run"], x["e_first"] - src["e_low"])
            para("%s: first SCF energy %+.4f Eh above its lowest%s." % (x["run"], x["e_first"] - x["e_low"], extra), "   ")

    title("R54.2  THE WAY OF STARTING: THE CAPTURE SOLUTION, STARTED AGAIN AT ITS OWN LEVEL")
    para("""DI7_CAM_S_d11_R starts CAM from the orbitals of D53's DI6_CAM_F_d11_R. REPRODUCED: it is measured, its count holds, it
STARTED FROM ITS FILE, and it TOOK ELECTRONS with %.2f e or more leaving the atoms. Beside it: whether it is THE SAME SOLUTION as
D53's run.""" % D_LW)
    say("")
    say("   run                 what it did        left      energy against D53's run, eV   <S**2>")
    say("   %-19s %-18s %-9s %-30s %s" % (OLD[("CAM", "R")]["run"], OLD[("CAM", "R")]["did"], f(OLD[("CAM", "R")]["left"]), "-", f(OLD[("CAM", "R")].get("s2"), "%.4f")))
    say("   %-19s %-18s %-9s %-30s %s" % (s_c["run"], s_c["did"], f(s_c["left"]),
        f((s_c["E"] - OLD[("CAM", "R")]["E"]) * EV if s_c["status"] == "measured" and OLD[("CAM", "R")]["status"] == "measured" else None, "%+.6f"), f(s_c.get("s2"), "%.4f")))
    say("")
    say("   THE WAY OF STARTING IS REPRODUCED." if reproduced else "   THE WAY OF STARTING IS NOT REPRODUCED. The runs from the capture solution are NOT READ.")
    say("   Against D53's run it is %s." % {True: "THE SAME SOLUTION", False: "ANOTHER SOLUTION", None: "not compared: a run is not measured or not read"}[rep])

    title("R54.3  THE SECOND LEVEL FROM ITS OWN START")
    para("""Each run against D53's control by the four legs, as D53 read its own: "%s" Beside it: whether it is THE SAME SOLUTION as
D53's run of the same input.""" % Q_B)
    say("")
    say("   state  run                 the run          left      gap, eV   gain      against D53's run")
    for stt in ("R", "TS"):
        x = own[stt]["x"]
        say(("   %-6s %-19s %-16s %-9s %-9s %-9s %s" % (stt, x["run"], x["did"] if not x["noread"] else "NOT READ", f(x["left"]), f(x["z"]["gap"], "%.4f"), f(x["gain"]),
             {True: "THE SAME SOLUTION", False: "ANOTHER SOLUTION", None: "-"}[own[stt]["same"]])).rstrip())
    say("")
    for stt in ("R", "TS"):
        x = own[stt]["x"]
        if x["noread"]:
            para("%s is NOT READ: %s." % (x["run"], x["noread"]), "   ")
        elif x["status"] != "measured":
            para("%s: %s." % (x["run"], x["status"]), "   ")
        else:
            if x["z"]["flags"]:
                para("%s, the legs that fired: %s." % (x["run"], "; ".join(x["z"]["flags"])), "   ")
            if x["z"]["held"]:
                para("%s, the legs withheld: %s." % (x["run"], "; ".join(x["z"]["held"])), "   ")
    for stt in ("R", "TS"):
        x, o = own[stt]["x"], OLD[("WXV", stt)]
        if x["status"] == "measured" and o["status"] == "measured":
            para("%s against D53's %s: energy %+.6f eV, left %s against %s e." % (x["run"], o["run"], (x["E"] - o["E"]) * EV, f(x["left"]), f(o["left"])), "   ")

    title("R54.4  THE SECOND LEVEL FROM THE CAPTURE SOLUTION")
    para("""DI7_WXV_S_d11_R starts WXV from the orbitals of D53's CAM run in which the electrons were taken. RETURNED: it ends as
THE SAME SOLUTION as the run it is held against. Otherwise it is A CAPTURE SOLUTION if %.2f e or more left the atoms and ANOTHER
SOLUTION WITHOUT CAPTURE if not; BELOW, ABOVE or LEVEL by its energy against that run (LEVEL: within %.4f eV).""" % (D_LW, E_TOL))
    say("")
    say("   held against: %s" % ref_w["run"])
    say("   'functions': the charge on the six sets of functions, against the control. 'against it': its")
    say("   energy minus the energy of the run it is held against, eV.")
    say("")
    say("   run                 left      functions  <S**2>   last analysis   against it   gain")
    say(("   %-19s %-9s %-10s %-8s %-15s %-12s %s" % (s_w["run"], f(s_w["left"]), f(s_w["ghost"]), f(s_w.get("s2"), "%.4f"), s_w["stab"].get("last") or "-", f(sw[1], "%+.6f"), f(s_w["gain"]))).rstrip())
    say("")
    say("   %s%s" % (sw[0], (": " + sw[2]) if sw[2] else "."))
    if s_w["status"] == "measured" and not s_w["noread"]:
        para("By the four legs against the control this run reads %s%s." % (
            s_w["did"].upper() if s_w["did"] == "not decided" else s_w["did"],
            (": " + "; ".join(s_w["z"]["flags"] + s_w["z"]["held"])) if (s_w["z"]["flags"] or s_w["z"]["held"]) else ""), "   ")

    title("R54.5  THE DESIGN AT THE SECOND LEVEL, AND THE TWO LEVELS TOGETHER")
    para("""At R the word follows the lowest solution found. TAKES ELECTRONS if the run from its own start TOOK ELECTRONS; or the run
from the capture solution is A CAPTURE SOLUTION that is BELOW or LEVEL; or it is ANOTHER SOLUTION WITHOUT CAPTURE that is BELOW or
LEVEL and itself reads TOOK ELECTRONS by the four legs. DOES NOT TAKE if the run from its own start DID NOT TAKE and the run from
the capture solution RETURNED, or is a solution ABOVE, or is ANOTHER SOLUTION WITHOUT CAPTURE that is BELOW or LEVEL and itself
reads DID NOT TAKE by the four legs. NOT DECIDED otherwise, and if a count does not hold. At TS no capture solution is known at
any level, and the word is that of the run from its own start. Together with the working functional's cell of D53's committed
table, by the words of R53.6.""")
    say("")
    say("   state  CAM (D53)        WXV (D54)         the two levels together")
    for stt, w, t in (("R", wR, tR), ("TS", wTS, tTS)):
        say("   %-6s %-16s %-17s %s%s" % (stt, w53[("d11", stt)]["CAM_did"], w[0], t[0], (" (%s)" % t[1]) if t[1] else ""))
    for stt, w, t in (("R", wR, tR), ("TS", wTS, tTS)):
        if t[0] == "NOT DECIDED":
            para("NOT DECIDED at %s: %s." % (stt, t[2]), "   ")
    say("")
    para("D53's result stands as committed: there the design is NOT DECIDED at both states. This is D54's reading of the same two cells.", "   ")

    title("R54.6  HARTREE-FOCK FROM THE CAPTURE SOLUTION, THE BRACKET")
    para("DI7_HF_S_d11_R against D53's DI6_HF_F_d11_R, by the words of R54.4. Reported; it changes no word above.")
    say("")
    say("   run                 left      <S**2>   last analysis   energy against D53's run, eV")
    say(("   %-19s %-9s %-8s %-15s %s" % (s_h["run"], f(s_h["left"]), f(s_h.get("s2"), "%.4f"), s_h["stab"].get("last") or "-", f(sh[1], "%+.6f"))).rstrip())
    say("")
    say("   %s%s" % (sh[0], (": " + sh[2]) if sh[2] else "."))

    title("WHAT THIS CANNOT SAY (section 6 of the pre-registration, in short; its own words govern)")
    for t in ("One frame, frozen geometries, in vacuo, one old design, charges of one size.",
              "A run started from the capture solution finds the solution nearest to that start. If it RETURNED, no capture solution was found from this start; that is not a proof that none exists.",
              "The stability analysis missed lower solutions in D53 (its Hartree-Fock controls with functions). STABLE is reported and is not a proof.",
              "L and L2 are bars for one whole electron (D52).",
              "D53's result stands as committed. D54 adopts nothing."):
        para(t, "   ")
    say("")
    para("Runs measured: %d of %d. %s The second level from the capture solution: %s. The design, the two levels together: at R, %s%s; at TS, %s%s." % (
        len(runs) - len(bad), len(runs), "THE WAY OF STARTING IS REPRODUCED." if reproduced else "THE WAY OF STARTING IS NOT REPRODUCED.", sw[0],
        tR[0], (" (%s)" % tR[1]) if tR[1] else "", tTS[0], (" (%s)" % tTS[1]) if tTS[1] else ""))

    # ------------------------------------------------------------ the tables
    def tv(v, fmt="%.6f"):
        return "" if v is None else fmt % v
    cols = ["run", "job", "level", "kind", "state", "start_orbitals_of", "status", "E_Eh", "count_minus_118", "count_holds", "basis_functions", "started_from_its_file", "first_SCF_energy_minus_lowest_Eh", "level_shift_set", "TRAH",
            "stability_first", "stability_last", "S2", "lowdin_rise_e", "dipole_ratio", "gap_eV", "flags", "held", "sum_24_minus_control_e", "ghosts_minus_control_e", "word",
            "what_the_run_did", "how_it_took", "gain_eV", "not_read_because"]
    t1 = ["\t".join(cols)]
    for x in [RUN[r["run"]] for r in runs] + [OLD[k] for k in (("WXV", "R"), ("WXV", "TS"), ("CAM", "R"), ("HF", "R"))]:
        z = x["z"]
        t1.append("\t".join([x["run"], x.get("job", "D53"), x["level"], x.get("kind", "F"), x["state"], x.get("start", ""), x["status"], tv(x["E"], "%.10f"), tv(x["count_d"], "%+.6f"),
                             "" if x["level"] == "HF" else ("yes" if x["count_ok"] else "no"), "" if x["nbf"] is None else str(x["nbf"]), {True: "yes", False: "no", None: ""}[x.get("started")] if x.get("kind") else "",
                             tv(x["e_first"] - x["e_low"] if x.get("e_first") is not None else None, "%.6f"), "yes" if x["lshift"] else "no", "yes" if x["trah"] else "no",
                             x["stab"].get("first") or "", x["stab"].get("last") or "", tv(x.get("s2"), "%.6f"), tv(z["rise"], "%+.6f"), tv(z["ratio"], "%.4f"), tv(z["gap"], "%.4f"),
                             "; ".join(z["flags"]), "; ".join(z["held"]), tv(x["left"], "%+.6f"), tv(x["ghost"], "%+.6f"), x["word"], x["did"], x["how"], tv(x["gain"], "%+.6f"), x["noread"] or ""]))
    t2 = ["\t".join(["state", "CAM_did_D53", "CAM_bar_D53", "WXV_own_start_did", "WXV_from_the_capture_solution", "the_design_at_the_second_level", "the_two_levels_together", "mark",
                     "not_decided_because", "L_eV", "L2_eV", "the_way_of_starting"])]
    for stt, o, s, w, t in (("R", o_r, sw[0], wR, tR), ("TS", o_ts, "", wTS, tTS)):
        t2.append("\t".join([stt, w53[("d11", stt)]["CAM_did"], w53[("d11", stt)]["CAM_bar"], o["did"], s, w[0], t[0], t[1], t[2], "%.6f" % L, "%.6f" % L2,
                             "REPRODUCED" if reproduced else "NOT REPRODUCED"]))
    t3 = ["\t".join(["run", "held_against", "word", "E_minus_reference_eV", "left_e", "left_of_reference_e", "S2", "stability_last", "gain_eV", "not_read_because"])]
    for x, rf, w in ((s_c, OLD[("CAM", "R")], ("REPRODUCED" if reproduced else "NOT REPRODUCED", (s_c["E"] - OLD[("CAM", "R")]["E"]) * EV if s_c["status"] == "measured" and OLD[("CAM", "R")]["status"] == "measured" else None, "")),
                     (o_r, OLD[("WXV", "R")], ({True: "THE SAME SOLUTION", False: "ANOTHER SOLUTION", None: "NOT READ"}[own["R"]["same"]], (o_r["E"] - OLD[("WXV", "R")]["E"]) * EV if own["R"]["same"] is not None else None, "")),
                     (o_ts, OLD[("WXV", "TS")], ({True: "THE SAME SOLUTION", False: "ANOTHER SOLUTION", None: "NOT READ"}[own["TS"]["same"]], (o_ts["E"] - OLD[("WXV", "TS")]["E"]) * EV if own["TS"]["same"] is not None else None, "")),
                     (s_w, ref_w, sw), (s_h, OLD[("HF", "R")], sh)):
        t3.append("\t".join([x["run"], rf["run"], w[0], tv(w[1], "%+.6f"), tv(x["left"], "%+.6f"), tv(rf["left"], "%+.6f"), tv(x.get("s2"), "%.6f"), x["stab"].get("last") or "", tv(x["gain"], "%+.6f"), w[2]]))
    text = "\n".join(l.rstrip() for l in out) + "\n"
    for p, t in ((a.result, text), (a.tsv, "\n".join(t1) + "\n"), (tsv2, "\n".join(t2) + "\n"), (tsv3, "\n".join(t3) + "\n")):
        with open(p, "w") as fh:
            fh.write(t)
    sys.stdout.write(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
