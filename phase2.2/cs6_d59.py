#!/usr/bin/env python3
"""cs6_d59.py - the reader of STAGE 6, DIAGNOSTIC 59, committed before any of its outputs exists.

AMENDMENT 11, A11.5 item 5; AMENDMENT 13, A13.4 item 1 (b). It reads, for the ten survivors of D56 with the
largest working-level lowering:
  R59.1 the functions at their positive charges, by D53's committed functions exactly as D53 read its
        arrangements, each design with its own control with functions; the two levels together by D53's
        R53.6; Hartree-Fock beside; and one word for each design;
  R59.2 the helium wall on their positive charges, each run read by D56's committed functions against
        D56's committed controls with the dipole leg withheld as D32.3 withheld it, the word of A9.4 (i)
        formed behind the wall beside D56's;
  R59.3 the wall alone.
It names no design, keeps or drops none, and changes no word of D56.

usage (as cs6_d59_harvest.sh calls it):
   python3 cs6_d59.py --res <results repo> --code <code repo> --prereg <pre-registration>
                      --result <result document to write> --tsv <table of runs to write>
It also writes d59_functions.tsv and d59_wall.tsv beside --tsv. If a gate fails it writes nothing.
"""
import argparse, csv, hashlib, importlib.util, os, re, subprocess, sys, textwrap

PIN_PREREG = "96b9d5340eda65993fcadb9ce941d816f561b51b3ce9a9219ba7af20b2716301"
S6REL = "phase2.2/cs_stage6"
DIR = "d59"
H2K = 627.5094740631
EV = 27.211386245988          # D53's
NEL = 118
COUNT_TOL = 0.02              # e, D53's
L_PRINTED, L2_PRINTED = "0.3127", "0.1151"
TEN = ["S_n10_k1", "S_n10_k3", "S_n10_k2", "S_n08_k1", "S_n08_k3", "S_n08_k2", "S_n06_k1", "S_n06_k2", "S_n06_k3", "S_n05_k2"]  # A13.4 item 1 (b)
W = 96


def sha(p):
    with open(p, "rb") as f:
        return hashlib.sha256(f.read()).hexdigest()


def rd(p):
    with open(p, errors="replace") as f:
        return f.read()


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def stop(why):
    print("THE READER STOPPED AT A GATE: " + why)
    print("Nothing was written.")
    sys.exit(1)


def echoed_input(text):
    out = []
    for l in text.split("\n"):
        if "****END OF INPUT****" in l:
            break
        m = re.match(r"^\|\s*\d+>( (.*))?$", l)
        if m:
            out.append((m.group(2) or "").rstrip())
    return out


def strip(L):
    L = [x.rstrip() for x in L]
    return L[:max([i + 1 for i, x in enumerate(L) if x.strip()] or [0])]


def main():
    ap = argparse.ArgumentParser()
    for a in ("res", "code", "prereg", "result", "tsv"):
        ap.add_argument("--" + a, required=True)
    a = ap.parse_args()
    S6 = os.path.join(a.res, S6REL)
    D = os.path.join(S6, DIR)
    P2 = os.path.join(a.code, "phase2.2")
    NEED = ("cs6_d59_prepare.py", "cs6_d59.py", "cs6_d59_harvest.sh", "cs6_d53.py", "cs6_d56.py", "cs6_d40.py", "cs6_d33b.py", "cs6_audit.py", "cs6_d47_prepare.py")
    # GH.1, GH.2
    if not os.path.isfile(a.prereg) or sha(a.prereg) != PIN_PREREG:
        stop("the pre-registration given is not the committed STAGE6_DIAGNOSTIC59_20261125.txt.")
    for f in NEED:
        q = subprocess.run(["git", "-C", a.code, "diff", "--quiet", "HEAD", "--", "phase2.2/" + f])
        r = subprocess.run(["git", "-C", a.code, "cat-file", "-e", "HEAD:phase2.2/" + f], stderr=subprocess.DEVNULL)
        if q.returncode or r.returncode:
            stop("%s is not committed, or differs from its committed version." % f)
    prep = load(os.path.join(P2, "cs6_d59_prepare.py"), "d59prep")
    R53 = load(os.path.join(P2, "cs6_d53.py"), "r53_59")
    R56 = load(os.path.join(P2, "cs6_d56.py"), "r56_59")
    M = dict((k, load(os.path.join(P2, fn), "m59_" + k)) for k, fn in (
        ("d40", "cs6_d40.py"), ("d33b", "cs6_d33b.py"), ("audit", "cs6_audit.py"), ("r53", "cs6_d53.py"), ("p47", "cs6_d47_prepare.py")))
    geom = dict((st, M["p47"].read_geom(os.path.join(S6, M["p47"].GEOM[st]))) for st in ("R", "TS"))
    T = (M["audit"].D_LW, M["audit"].R_DIP, M["audit"].GAP_ABS, M["audit"].R_GAP)
    D_LW = M["audit"].D_LW
    # GH.3 the inputs, and the designs of A11.5 item 5
    DES = prep.designs(S6)
    if DES != TEN:
        stop("the designs A11.5 item 5 gives from the committed d56_survivors.tsv are not the ten A13.4 item 1 (b) names.")
    F = prep.build(S6)
    tracked = subprocess.run(["git", "-C", a.res, "ls-files", S6REL + "/" + DIR], stdout=subprocess.PIPE, check=True).stdout.decode().split("\n")
    tracked = set(os.path.relpath(t, S6REL + "/" + DIR) for t in tracked if t)
    if not set(F) <= tracked:
        stop("files of d59/ that the program forms are not committed.")
    for p in sorted(F):
        if rd(os.path.join(D, p)) != F[p]:
            stop("d59/%s is not the file the committed program forms from the committed files." % p)
    RUNS, JOBS = prep.runs(S6)
    # GH.4 every output echoes its committed input
    for r, _, _, _ in RUNS:
        o = os.path.join(D, r, "job.out")
        if os.path.isfile(o):
            got = echoed_input(rd(o))
            if got and strip(got) != strip(rd(os.path.join(D, r, "job.inp")).split("\n")):
                stop("%s/job.out does not echo the committed job.inp beside it." % r)
    # GH.6 the rules of D52, as D53 checked them
    rules = dict(l.split("\t")[:2] for l in rd(os.path.join(S6, "d52_rules.tsv")).split("\n")[1:] if "\t" in l)
    try:
        L = float(rules["rule_C_L_eV"])
        L2 = 3.0 * float(rules["largest_distance_from_the_limit_WXV_eV"])
    except (KeyError, ValueError):
        L = L2 = None
    if L is None or "%.4f" % L != L_PRINTED or "%.4f" % L2 != L2_PRINTED or rules.get("rule_B_second_level") != "wB97X-V" \
            or rules.get("rule_A_set") != "X2" or rules.get("rule_A_IntAcc") != "12.0":
        stop("d52_rules.tsv does not hold the rules the pre-registration names (X2, IntAcc 12, wB97X-V, L %s eV, L2 %s eV)." % (L_PRINTED, L2_PRINTED))

    def charges_of(design):
        return [(q, tuple(float(v) for v in x)) for q, x in prep.charges_of(S6, design)]
    npos = dict((d, sum(1 for q, _ in charges_of(d) if q > 0)) for d in DES)

    # D56's controls, and GH.7
    C = {}
    for lv, _, _ in prep.WLEVELS:
        for k in ("ctl", "bare"):
            for st in ("R", "TS"):
                C[(lv, k, st)] = R56.read_run(M, os.path.join(S6, "d56", "%s_C_%s_%s" % (lv, k, st)), k)

    def word_of(get, charges, walled):
        """D56's formation of the word of A9.4 (i), from runs get(lv, kind, state); walled withholds the dipole leg (D57's)"""
        o = {}
        for lv, _, _ in prep.WLEVELS:
            di = dict((st, get(lv, "di", st)) for st in ("R", "TS"))
            mo = dict((st, get(lv, "mo", st)) for st in ("R", "TS"))
            ctl = dict((st, C[(lv, "ctl", st)]) for st in ("R", "TS"))
            if walled:
                z = {}
                for st in ("R", "TS"):
                    ok = ctl[st]["status"] == "measured" and di[st]["status"] in ("measured", "DID NOT CONVERGE")
                    zz = M["r53"].legs(di[st], ctl[st], geom[st], charges, T) if ok else dict(M["r53"].NOZ, flags=[], held=[])
                    zz["flags"] = [f for f in zz["flags"] if not f.startswith("dipole")]
                    zz["held"] = [h for h in zz["held"] if "dipole" not in h]
                    z[st] = zz
                if any(z[st]["flags"] for st in z):
                    four = "UNRELIABLE"
                elif any(di[st]["status"] != "measured" or ctl[st]["status"] != "measured" for st in z):
                    four = "NOT CLASSIFIED"
                elif any(z[st]["held"] for st in z):
                    four = "INCONCLUSIVE"
                else:
                    four = "RELIABLE"
            else:
                four, z = R56.four_legs(M, di, ctl, geom, charges)
            fifth, fw = R56.fifth_leg(di, mo)
            o[lv] = dict(four=four, z=z, fifth=fifth, fw=fw, dde=R56.dde_of(di, ctl), di=di)
        o["agree"] = R56.agreement(o["TZ"]["dde"], o["SV"]["dde"])
        o["at_tz"] = R56.at_level(o["TZ"]["four"], o["TZ"]["fifth"], o["TZ"]["fw"], o["agree"])
        o["work"] = "PASSES THE FIVE LEGS" if o["SV"]["four"] == "RELIABLE" and o["SV"]["fifth"] == "PASSED" else (
            "FAILS" if o["SV"]["four"] == "UNRELIABLE" or any(o["SV"]["fw"][st][1] == "METASTABLE" for st in ("R", "TS")) else "OPEN")
        o["word"] = "SURVIVOR" if o["work"] == "PASSES THE FIVE LEGS" and o["at_tz"] == "PASSES" else (
            "NOT A SURVIVOR: FAILS" if "FAILS" in (o["work"], o["at_tz"]) else "NOT A SURVIVOR: OPEN")
        return o

    surv = dict((r["design"], r) for r in csv.DictReader(open(os.path.join(S6, "d56_survivors.tsv")), delimiter="\t"))
    for d in DES:
        o = word_of(lambda lv, k, st, d=d: R56.read_run(M, os.path.join(S6, "d56", "%s_%s_%s_%s" % (lv, d, k, st)), k), charges_of(d), False)
        if o["word"] != surv[d]["the_design"] or o["word"] != "SURVIVOR":
            stop("D56's functions give %s for %s here, and D56's committed word is %s." % (o["word"], d, surv[d]["the_design"]))

    status = {}
    # ------------------------------------------------------------ R59.1 the functions
    RUN = {}
    nb_want = dict((r, prep.nb_of(S6, r)) for r, _, _, _ in RUNS if r.split("_")[0] in prep.FLEVELS)
    for r, j, cap, kind in RUNS:
        p = r.split("_")
        if p[0] not in prep.FLEVELS:
            continue
        lv, k, st = p[0], p[1], p[-1]
        d = "_".join(p[2:-1]) if k != "GN" else None
        x = R53.read_di(M["d40"], M["d33b"], M["audit"], os.path.join(D, r), "pc" if k in ("F", "N") else "plain")
        x.update(run=r, level=lv, kind=k, design=d, state=st)
        x["control"] = {"F": "%s_GF_%s_%s" % (lv, d, st), "GF": "%s_GN_%s" % (lv, st), "N": "%s_GN_%s" % (lv, st)}.get(k)
        x["partner"] = "%s_N_%s_%s" % (lv, d, st) if k == "F" and lv != "HF" else "-"
        RUN[r] = x
        status[r] = x["status"] if x.get("exists") else "NOT RUN"
        # GH.5 the basis functions the input asks for
        if x["status"] == "measured" and x["nbf"] != nb_want[r]:
            stop("%s reports %s basis functions; its input asks for %d." % (r, x["nbf"], nb_want[r]))
    for x in RUN.values():
        if x["level"] == "HF":
            x["count_ok"], x["count_d"] = True, None
        elif x["status"] != "measured" or x["count"] is None:
            x["count_ok"], x["count_d"] = False, None
        else:
            x["count_d"] = x["count"] - NEL
            x["count_ok"] = abs(x["count_d"]) <= COUNT_TOL
    for x in RUN.values():
        x["noread"], x["z"], x["word"], x["took"], x["how"], x["gain"], x["bar"], x["ghost"], x["holds"], x["nogain"] = None, dict(R53.NOZ), "", ("not decided", None), "", None, "", None, "", ""
    for x in RUN.values():
        if x["kind"] == "GN":
            x["word"] = "STABLE" if R53.stable(x) else ("NOT MEASURED" if x["status"] != "measured" else "NOT STABLE")
    for x in RUN.values():
        if x["kind"] not in ("GF", "N"):
            continue
        c0 = RUN[x["control"]]
        if not R53.stable(c0):
            x["noread"] = "its plain control %s is %s" % (x["control"], c0["word"])
        x["z"] = R53.legs(x, c0, geom[x["state"]], charges_of(x["design"]) if x["kind"] == "N" else [], T) if not x["noread"] else dict(R53.NOZ)
        x["word"] = R53.word1(x)
        if x["kind"] == "GF":
            n = npos[x["design"]]
            x["holds"] = R53.control_holds(x)
            mv = x["z"]["mv"]
            x["left"] = mv[3] if mv else None
            x["ghost"] = sum(x["lw"][24:24 + n]) if x["status"] == "measured" and x["lw"] and len(x["lw"]) >= 24 + n else None
    for x in RUN.values():
        if x["kind"] != "F":
            continue
        n = npos[x["design"]]
        c0 = RUN[x["control"]]
        why = None
        if c0["status"] != "measured":
            why = "its control %s is %s" % (x["control"], c0["status"])
        elif c0["noread"]:
            why = "its control %s is not read: %s" % (x["control"], c0["noread"])
        elif not c0["stab"].get("ran") or c0["stab"].get("last") is None:
            why = "the stability analysis of its control %s was not read" % x["control"]
        elif c0["stab"]["last"] != "STABLE":
            why = "its control %s is UNSTABLE at the last analysis" % x["control"]
        x["noread"] = why
        x["z"] = R53.legs(x, c0, geom[x["state"]], charges_of(x["design"]), T) if not x["noread"] else dict(R53.NOZ)
        x["word"] = R53.word1(x)
        x["took"] = R53.took(x, D_LW)
        x["how"] = ("ONTO THE FUNCTIONS" if x["took"][1] >= D_LW else "A LEG ONLY") if x["took"][0] == "TOOK ELECTRONS" else ""
        if x["status"] == "measured" and x["lw"] and c0.get("lw") and len(x["lw"]) >= 24 + n and len(c0["lw"]) >= 24 + n and not x["noread"]:
            x["ghost"] = sum(x["lw"][24:24 + n]) - sum(c0["lw"][24:24 + n])
        if x["partner"] != "-":
            nn = RUN[x["partner"]]
            gn = RUN[nn["control"]]
            four = (x, c0, nn, gn)
            if all(v["status"] == "measured" for v in four) and all(v["count_ok"] for v in four):
                x["gain"] = ((x["E"] - c0["E"]) - (nn["E"] - gn["E"])) * EV
            else:
                x["nogain"] = "; ".join("%s %s" % (v["run"], v["status"] if v["status"] != "measured" else "has a count that does not hold") for v in four
                                        if v["status"] != "measured" or not v["count_ok"])
        if x["took"][0] == "TOOK ELECTRONS" and x["level"] in ("CAM", "WXV"):
            bar = L if x["level"] == "CAM" else L2
            x["bar"] = "not formed" if x["gain"] is None else ("ABOVE" if -x["gain"] > bar else "WITHIN")

    def cell(level, d, st):
        x = RUN["%s_F_%s_%s" % (level, d, st)]
        c0 = RUN[x["control"]]
        why = x["noread"] or (None if x["status"] == "measured" else ("the run DID NOT CONVERGE" if x["status"] == "DID NOT CONVERGE" else "the run is NOT MEASURED")) or \
            (", ".join(x["z"]["held"]) if x["took"][0] == "not decided" and x["z"]["held"] else None)
        return dict(did=x["took"][0], left=x["took"][1], gain=x["gain"], bar=x["bar"], count_ok=x["count_ok"] and c0["count_ok"], why=why, how=x["how"], x=x,
                    holds=c0["holds"], control=c0["run"], nogain=x["nogain"])

    FW, DW = [], {}
    for d in DES:
        wctl = [RUN["WXV_GN_%s" % st] for st in ("R", "TS")] + [RUN["WXV_GF_%s_%s" % (d, st)] for st in ("R", "TS")]
        second = all(R53.stable(x) for x in wctl)
        words = {}
        for st in ("R", "TS"):
            c, w, h = cell("CAM", d, st), cell("WXV", d, st), cell("HF", d, st)
            word, mark, reason = R53.together(c, w, second)
            FW.append(dict(design=d, state=st, c=c, w=w, h=h, word=word, mark=mark, reason=reason))
            words[st] = word
        if "TAKES ELECTRONS AT BOTH LEVELS" in words.values():
            DW[d] = "TAKES ELECTRONS WITH FUNCTIONS AT ITS POSITIVE CHARGES"
        elif all(v == "HOLDS AT BOTH LEVELS" for v in words.values()):
            DW[d] = "HOLDS WITH FUNCTIONS AT ITS POSITIVE CHARGES"
        else:
            DW[d] = "NOT SETTLED WITH FUNCTIONS"

    # ------------------------------------------------------------ R59.2 the helium wall, R59.3 the wall alone
    def run_of(r, kind):
        d = os.path.join(D, r)
        k = kind if os.path.isfile(os.path.join(d, "pc.pc")) else {"di": "ctl", "mo": "bare"}[kind]
        x = R56.read_run(M, d, k)
        status[r] = x["status"] if os.path.isfile(os.path.join(d, "job.out")) else "NOT RUN"
        return x

    WALL, Z = {}, {}
    for d in DES:
        WALL[d] = word_of(lambda lv, k, st, d=d: run_of("%s_WHe_%s_%s_%s" % (lv, d, k, st), k), charges_of(d), True)
        for lv, _, _ in prep.WLEVELS:
            zz = dict((st, run_of("%s_ZHe_%s_di_%s" % (lv, d, st), "di")) for st in ("R", "TS"))
            ctl = dict((st, C[(lv, "ctl", st)]) for st in ("R", "TS"))
            rise = dict((st, (zz[st]["c"]["maxlw"] - ctl[st]["c"]["maxlw"]) if zz[st]["status"] == "measured" and zz[st].get("c") and ctl[st].get("c") else None) for st in ("R", "TS"))
            Z[(d, lv)] = dict(runs=zz, rise=rise, dde=R56.dde_of(zz, ctl))
            WALL[d][lv]["dde_wall"] = R56.dde_of(WALL[d][lv]["di"], zz) if Z[(d, lv)]["dde"] is not None else None

    # ------------------------------------------------------------ the tables
    base = os.path.dirname(a.tsv)
    f4 = lambda v: "" if v is None else "%.4f" % v
    with open(a.tsv, "w") as f:
        f.write("run\tjob\tstatus\n" + "".join("%s\t%s\t%s\n" % (r, j, status.get(r, "NOT READ")) for r, j, _, _ in RUNS))
    with open(os.path.join(base, "d59_functions.tsv"), "w") as f:
        f.write("design\tstate\tCAM\tCAM_left_e\tCAM_gain_eV\tCAM_bar\tWXV\tWXV_left_e\tWXV_gain_eV\tWXV_bar\tHF\tword\tmark\twhy_not_decided\tdesign_word\n")
        for x in FW:
            f.write("\t".join([x["design"], x["state"], x["c"]["did"], f4(x["c"]["left"]), f4(x["c"]["gain"]), x["c"]["bar"], x["w"]["did"], f4(x["w"]["left"]),
                               f4(x["w"]["gain"]), x["w"]["bar"], x["h"]["did"], x["word"], x["mark"], x["reason"], DW[x["design"]]]) + "\n")
    with open(os.path.join(base, "d59_wall.tsv"), "w") as f:
        f.write("design\tlevel\tfour_legs\tfifth_leg\tddE_kcal\tddE_with_walls_in_place_kcal\twall_alone_ddE_kcal\twall_alone_rise_R_e\twall_alone_rise_TS_e\tagreement\tat_TZ\tword_behind_the_helium_wall\tD56_word\n")
        for d in DES:
            for lv, _, _ in prep.WLEVELS:
                o, z = WALL[d], Z[(d, lv)]
                f.write("\t".join([d, lv, o[lv]["four"], o[lv]["fifth"], f4(o[lv]["dde"]), f4(o[lv]["dde_wall"]), f4(z["dde"]), f4(z["rise"]["R"]), f4(z["rise"]["TS"]),
                                   o["agree"], o["at_tz"], o["word"], surv[d]["the_design"]]) + "\n")

    # ------------------------------------------------------------ the document
    P = lambda s, ind=0: "\n".join(textwrap.wrap(" ".join(s.split()), W, initial_indent=" " * ind, subsequent_indent=" " * ind, break_long_words=False))
    head = subprocess.run(["git", "-C", a.res, "rev-parse", "HEAD"], stdout=subprocess.PIPE).stdout.decode().strip()
    Ls = ["STAGE 6 - DIAGNOSTIC 59 - RESULT", "The ten survivors of D56 with the largest lowering, with functions at their positive charges and behind the helium wall", "=" * W, ""]
    Ls.append(P("Read by cs6_d59.py under STAGE6_DIAGNOSTIC59_20261125.txt, which was committed with the inputs before any job was submitted."))
    Ls.append("   results repository at %s" % head[:12])
    for f in NEED[:2] + NEED[3:]:
        Ls.append("   %-22s sha256 %s" % (f, sha(os.path.join(P2, f))))
    Ls += ["", P("THE GATES, as they held: GH.1 the pre-registration; GH.2 the programs; GH.3 the inputs, and the ten designs of A11.5 item 5 (%s); GH.4 every output echoes its committed input; GH.5 every measured run of the functions reports the basis functions its input asks for; GH.6 d52_rules.tsv holds the rules named (L %.4f eV, L2 %.4f eV); GH.7 D56's functions give D56's word SURVIVOR for each of the ten." % (", ".join(DES), L, L2)), ""]
    cnt = {}
    for r, _, _, _ in RUNS:
        cnt[status.get(r, "NOT READ")] = cnt.get(status.get(r, "NOT READ"), 0) + 1
    Ls += ["   runs: " + "; ".join("%s %d" % kv for kv in sorted(cnt.items())) + " (of %d)" % len(RUNS), ""]
    Ls += ["=" * W, "R59.1  THE FUNCTIONS AT THE POSITIVE CHARGES (D53's reading; LEFT in e, gain in eV)", "=" * W, ""]
    Ls.append("   design     st   CAM                         WXV                         HF")
    for x in FW:
        cc = lambda k: "%-14s %-7s %-6s" % (k["did"][:14], f4(k["gain"]) or "-", (k["bar"] or "")[:6])
        Ls.append("   %-10s %-3s  %s  %s  %s" % (x["design"], x["state"], cc(x["c"]), cc(x["w"]), x["h"]["did"][:14]))
        Ls.append("        -> %s%s%s" % (x["word"], ("  (" + x["mark"] + ")") if x["mark"] else "", ("  [" + x["reason"] + "]") if x["reason"] else ""))
    Ls.append("")
    for d in DES:
        Ls.append("   %-10s %s" % (d, DW[d]))
    Ls += ["", P("Reported beside D56's word and nothing more (A11.5 item 5). A design that takes electrons is reported as such; whether it may then be quoted is for the amendment that follows this result."), ""]
    Ls += ["=" * W, "R59.2  THE HELIUM WALL ON THE POSITIVE CHARGES (the dipole leg withheld, as D32.3; kcal/mol)", "=" * W, ""]
    for d in DES:
        o = WALL[d]
        Ls.append("   %s, behind the helium wall. D56's word: %s" % (d, surv[d]["the_design"]))
        for lv, _, _ in prep.WLEVELS:
            x = o[lv]
            Ls.append("      %s  four legs %-15s fifth %-28s ddE %-9s with walls in place %s" % (lv, x["four"], x["fifth"][:28], f4(x["dde"]) or "-", f4(x["dde_wall"]) or "-"))
        Ls.append("      %s: agreement %s; at def2-TZVPD %s; BEHIND THE HELIUM WALL: %s" % (d, o["agree"], o["at_tz"], o["word"]))
        Ls.append("")
    Ls += [P("Reported beside D56's bare reading and nothing more. No design is named, kept or dropped on it."), ""]
    Ls += ["=" * W, "R59.3  THE WALL ALONE (centres with no charge; kcal/mol, e)", "=" * W, ""]
    for d in DES:
        for lv, _, _ in prep.WLEVELS:
            z = Z[(d, lv)]
            Ls.append("   %-9s %s  barrier change %-9s Loewdin rise R %s TS %s%s" % (d, lv, f4(z["dde"]) or "-", f4(z["rise"]["R"]) or "-", f4(z["rise"]["TS"]) or "-",
                      "   PERTURBS THE DENSITY" if any(v is not None and v >= D_LW for v in z["rise"].values()) else ""))
    Ls += ["", "=" * W, "WHAT THIS CANNOT SAY (section 6 of the pre-registration, in short; its own words govern)", "=" * W, ""]
    Ls.append(P("One frame, frozen geometries, in vacuo, charges of 0.25 e that stand for no atom. Functions at the positive charges only cannot show an electron going to a negative charge. The helium wall on 0.25 e is a model beside the bare one."))
    ntake = sum(1 for d in DES if DW[d].startswith("TAKES"))
    nhold = sum(1 for d in DES if DW[d].startswith("HOLDS"))
    nsurv = sum(1 for d in DES if WALL[d]["word"] == "SURVIVOR")
    Ls += ["", "Runs measured: %d of %d. WITH FUNCTIONS: %d of 10 take electrons, %d hold, %d not settled. BEHIND THE HELIUM WALL: SURVIVOR %d of 10." % (
        cnt.get("measured", 0), len(RUNS), ntake, nhold, 10 - ntake - nhold, nsurv)]
    txt = "\n".join(Ls) + "\n"
    with open(a.result, "w") as f:
        f.write(txt)
    print(txt, end="")


if __name__ == "__main__":
    main()
