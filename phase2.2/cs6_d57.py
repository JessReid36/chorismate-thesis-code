#!/usr/bin/env python3
"""cs6_d57.py - the reader of STAGE 6, DIAGNOSTIC 57, committed before any of its outputs exists.

AMENDMENT 12, A12.4. It reads:
  R57.1 the designs that A12.4 item 1 walls (and item 6's design with every charge walled), each run
        read by D56's committed functions against D56's committed controls, with the dipole leg
        withheld as D32.3 withheld it; the word of A9.4 (i) formed behind each wall, beside D56's;
  R57.2 the wall alone (the Z and Y controls);
  R57.3 the wall against real ions at one site: ion(X) = ddE(X+) - ddE(ghost X), and which of the
        bare and walled +1 e lies nearer the ions, by the rule of the pre-registration.
It names no design, keeps or drops none, and changes no word of D56.

usage (as cs6_d57_harvest.sh calls it):
   python3 cs6_d57.py --res <results repo> --code <code repo> --prereg <pre-registration>
                      --result <result document to write> --tsv <table of runs to write>
It also writes d57_designs.tsv and d57_ions.tsv beside --tsv. If a gate fails it writes nothing.
"""
import argparse, csv, hashlib, importlib.util, os, re, subprocess, sys, textwrap

PIN_PREREG = "84aba779c8324bce2a064b64b19fca6c367724e5a6789927a236130910e1c2e0"
S6REL = "phase2.2/cs_stage6"
DIR = "d57"
H2K = 627.5094740631
TOL_I = 0.30                 # kcal/mol: DIAGNOSTIC 3's limit of an effect too small to count
NEL_I = {"Na": 128, "K": 136}
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
    NEED = ("cs6_d57_prepare.py", "cs6_d57.py", "cs6_d57_harvest.sh", "cs6_d56.py", "cs6_d53.py", "cs6_d40.py", "cs6_d33b.py", "cs6_audit.py", "cs6_d47_prepare.py")
    # GE.1, GE.2
    if not os.path.isfile(a.prereg) or sha(a.prereg) != PIN_PREREG:
        stop("the pre-registration given is not the committed STAGE6_DIAGNOSTIC57_20261125.txt.")
    for f in NEED:
        q = subprocess.run(["git", "-C", a.code, "diff", "--quiet", "HEAD", "--", "phase2.2/" + f])
        r = subprocess.run(["git", "-C", a.code, "cat-file", "-e", "HEAD:phase2.2/" + f], stderr=subprocess.DEVNULL)
        if q.returncode or r.returncode:
            stop("%s is not committed, or differs from its committed version." % f)
    prep = load(os.path.join(P2, "cs6_d57_prepare.py"), "d57prep")
    R56 = load(os.path.join(P2, "cs6_d56.py"), "r56")
    M = dict((k, load(os.path.join(P2, fn), "m57_" + k)) for k, fn in (
        ("d40", "cs6_d40.py"), ("d33b", "cs6_d33b.py"), ("audit", "cs6_audit.py"), ("r53", "cs6_d53.py"), ("p47", "cs6_d47_prepare.py")))
    geom = dict((st, M["p47"].read_geom(os.path.join(S6, M["p47"].GEOM[st]))) for st in ("R", "TS"))
    T = (M["audit"].D_LW, M["audit"].R_DIP, M["audit"].GAP_ABS, M["audit"].R_GAP)
    # GE.3 inputs and the rule of item 1
    F = prep.build(S6)
    tracked = subprocess.run(["git", "-C", a.res, "ls-files", S6REL + "/" + DIR], stdout=subprocess.PIPE, check=True).stdout.decode().split("\n")
    tracked = set(os.path.relpath(t, S6REL + "/" + DIR) for t in tracked if t)
    if not set(F) <= tracked:
        stop("files of d57/ that the program forms are not committed.")
    for p in sorted(F):
        if rd(os.path.join(D, p)) != F[p]:
            stop("d57/%s is not the file the committed program forms from the committed files." % p)
    WD = prep.walled_designs(S6)
    BW = prep.both_ways(S6)
    RUNS, JOBS = prep.runs(S6)
    # GE.4 echo
    for r, _, _ in RUNS:
        o = os.path.join(D, r, "job.out")
        if os.path.isfile(o):
            got = echoed_input(rd(o))
            if got and strip(got) != strip(rd(os.path.join(D, r, "job.inp")).split("\n")):
                stop("%s/job.out does not echo the committed job.inp beside it." % r)

    def charges_of(design):
        L = [l.split() for l in rd(os.path.join(S6, "d56", "SV_%s_di_R" % design, "pc.pc")).split("\n")[1:] if l.strip()]
        return [(float(c[0]), tuple(float(v) for v in c[1:4])) for c in L]

    C = {}
    for lv, _, _ in prep.LEVELS:
        for k in ("ctl", "bare"):
            for st in ("R", "TS"):
                C[(lv, k, st)] = R56.read_run(M, os.path.join(S6, "d56", "%s_C_%s_%s" % (lv, k, st)), k)

    def word_of(get, charges, walled):
        """D56's formation of the word of A9.4 (i), from runs get(lv, kind, state); walled withholds the dipole leg"""
        o = {}
        for lv, _, _ in prep.LEVELS:
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

    # GE.5 D56's functions, on D56's committed outputs of these designs, give D56's committed words
    surv = dict((r["design"], r) for r in csv.DictReader(open(os.path.join(S6, "d56_survivors.tsv")), delimiter="\t"))
    for d, n, rk in WD:
        o = word_of(lambda lv, k, st: R56.read_run(M, os.path.join(S6, "d56", "%s_%s_%s_%s" % (lv, d, k, st)), k), charges_of(d), False)
        if o["word"] != surv[d]["the_design"]:
            stop("D56's functions give %s for %s here, and D56's committed word is %s." % (o["word"], d, surv[d]["the_design"]))

    # the runs of d57
    status = {}

    def run_of(r, kind):
        d = os.path.join(D, r)
        k = kind if os.path.isfile(os.path.join(d, "pc.pc")) else {"di": "ctl", "mo": "bare"}[kind]
        x = R56.read_run(M, d, k)
        status[r] = x["status"] if os.path.isfile(os.path.join(d, "job.out")) else "NOT RUN"
        return x

    DES = []
    for d, n, rk in WD:
        for s in prep.SETS:
            o = word_of(lambda lv, k, st, s=s, d=d: run_of("%s_W%s_%s_%s_%s" % (lv, s, d, k, st), k), charges_of(d), True)
            DES.append(("W", d, s, o))
    if BW:
        for s in prep.SETS:
            o = word_of(lambda lv, k, st, s=s: run_of("%s_A%s_%s_%s_%s" % (lv, s, BW, k, st), k), charges_of(BW), True)
            DES.append(("A", BW, s, o))
    # R57.2 the wall alone, and ddE of the charges with the walls in place
    Z = {}
    for t, d in [("Z", d) for d, _, _ in WD] + ([("Y", BW)] if BW else []):
        for s in prep.SETS:
            for lv, _, _ in prep.LEVELS:
                zz = dict((st, run_of("%s_%s%s_%s_di_%s" % (lv, t, s, d, st), "di")) for st in ("R", "TS"))
                ctl = dict((st, C[(lv, "ctl", st)]) for st in ("R", "TS"))
                rise = dict((st, (zz[st]["c"]["maxlw"] - ctl[st]["c"]["maxlw"]) if zz[st]["status"] == "measured" and zz[st].get("c") and ctl[st].get("c") else None) for st in ("R", "TS"))
                Z[(t, d, s, lv)] = dict(runs=zz, rise=rise, dde=R56.dde_of(zz, ctl))
    for kind, d, s, o in DES:
        for lv, _, _ in prep.LEVELS:
            z = Z[("Z" if kind == "W" else "Y", d, s, lv)]
            o[lv]["dde_wall"] = R56.dde_of(o[lv]["di"], z["runs"]) if z["dde"] is not None else None

    # R57.3 the real ions
    def ion_run(r):
        o = os.path.join(D, r, "job.out")
        t = rd(o) if os.path.isfile(o) else ""
        E = re.findall(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", t)
        nel = re.search(r"Number of Electrons\s+NEL\s+\.\.\.\.\s+(\d+)", t)
        ok = "ORCA TERMINATED NORMALLY" in t and bool(E)
        status[r] = "measured" if ok else ("NOT RUN" if not t else "NOT MEASURED")
        return dict(E=float(E[-1]) if ok else None, nel=int(nel.group(1)) if nel else None)
    IO = {}
    for lv, _, _ in prep.LEVELS:
        X = dict(((k, st), ion_run("%s_I_%s_%s" % (lv, k, st))) for k in prep.IKINDS for st in ("R", "TS"))
        def dd(k):
            if any(X[(kk, st)]["E"] is None for kk in (k, "ctl") for st in ("R", "TS")):
                return None
            return ((X[(k, "TS")]["E"] - X[("ctl", "TS")]["E"]) - (X[(k, "R")]["E"] - X[("ctl", "R")]["E"])) * H2K
        v = dict((k, dd(k)) for k in prep.IKINDS if k != "ctl")
        nelok = all(X[(k, st)]["nel"] in (None, NEL_I.get(k, 118)) for k in prep.IKINDS for st in ("R", "TS"))
        ion = dict((x, None if v[x] is None or v["g" + x] is None else v[x] - v["g" + x]) for x in ("Na", "K"))
        res = dict(dde=v, ion=ion, nelok=nelok, words={})
        if None not in ion.values() and v["bare"] is not None:
            Tm, s = (ion["Na"] + ion["K"]) / 2, abs(ion["Na"] - ion["K"]) / 2
            tol = max(s, TOL_I)
            dist = dict((m, None if v[m] is None else abs(v[m] - Tm)) for m in ("bare", "He", "Ne", "Shi"))
            for m in ("He", "Ne", "Shi"):
                if dist[m] is None:
                    res["words"][m] = "NOT FORMED"
                elif dist[m] + tol < dist["bare"]:
                    res["words"][m] = "NEARER THE IONS THAN THE BARE CHARGE"
                elif dist[m] > dist["bare"] + tol:
                    res["words"][m] = "NOT NEARER"
                else:
                    res["words"][m] = "NOT DISTINGUISHED"
            res.update(T=Tm, s=s, tol=tol, dist=dist)
        IO[lv] = res
    support = "SUPPORTED" if all(IO[lv]["words"].get("He") == "NEARER THE IONS THAN THE BARE CHARGE" for lv in IO) else "NOT SUPPORTED"
    unformed = [lv for lv in IO if "He" not in IO[lv]["words"] or IO[lv]["words"]["He"] == "NOT FORMED"]
    support_line = support + (" (the reading was not formed at %s)" % " and ".join(unformed) if unformed else "")

    # the tables
    base = os.path.dirname(a.tsv)
    f2 = lambda v: "" if v is None else "%.4f" % v
    with open(a.tsv, "w") as f:
        f.write("run\tjob\tstatus\n" + "".join("%s\t%s\t%s\n" % (r, j, status.get(r, "NOT READ")) for r, j, _ in RUNS))
    with open(os.path.join(base, "d57_designs.tsv"), "w") as f:
        f.write("design\twalled\tset\tlevel\tfour_legs\tfifth_leg\tddE_kcal\tddE_with_walls_in_place_kcal\tagreement\tat_TZ\tword_behind_the_wall\tD56_word\n")
        for kind, d, s, o in DES:
            for lv, _, _ in prep.LEVELS:
                f.write("\t".join([d, "every charge" if kind == "A" else "positive charges", s, lv, o[lv]["four"], o[lv]["fifth"], f2(o[lv]["dde"]), f2(o[lv]["dde_wall"]),
                                   o["agree"], o["at_tz"], o["word"], surv[d]["the_design"]]) + "\n")
    with open(os.path.join(base, "d57_ions.tsv"), "w") as f:
        f.write("level\tmodel\tddE_kcal\tdistance_from_ions_kcal\tword\n")
        for lv in IO:
            for m in ("bare", "He", "Ne", "Shi", "He0", "Ne0", "Shi0", "Na", "K", "gNa", "gK"):
                f.write("%s\t%s\t%s\t%s\t%s\n" % (lv, m, f2(IO[lv]["dde"].get(m)), f2(IO[lv].get("dist", {}).get(m)), IO[lv]["words"].get(m, "")))
            for x in ("Na", "K"):
                f.write("%s\tion(%s)\t%s\t\t\n" % (lv, x, f2(IO[lv]["ion"][x])))

    # the document
    P = lambda s, ind=0: "\n".join(textwrap.wrap(" ".join(s.split()), W, initial_indent=" " * ind, subsequent_indent=" " * ind, break_long_words=False))
    head = subprocess.run(["git", "-C", a.res, "rev-parse", "HEAD"], stdout=subprocess.PIPE).stdout.decode().strip()
    L = ["STAGE 6 - DIAGNOSTIC 57 - RESULT", "D56's best-ranked failed designs behind the wall, and the wall measured against real ions", "=" * W, ""]
    L.append(P("Read by cs6_d57.py under STAGE6_DIAGNOSTIC57_20261125.txt, which was committed with the inputs before any job was submitted."))
    L.append("   results repository at %s" % head[:12])
    for f in NEED[:2] + NEED[3:]:
        L.append("   %-22s sha256 %s" % (f, sha(os.path.join(P2, f))))
    L += ["", P("THE GATES, as they held: GE.1 the pre-registration; GE.2 the programs; GE.3 the inputs, and the designs of A12.4 item 1 (%s) and item 6 (%s); GE.4 every output echoes its committed input; GE.5 D56's functions give D56's committed words for these designs." % (", ".join(d for d, _, _ in WD), BW)), ""]
    cnt = {}
    for r, _, _ in RUNS:
        cnt[status.get(r, "NOT READ")] = cnt.get(status.get(r, "NOT READ"), 0) + 1
    L += ["=" * W, "R57.1  THE WALLED DESIGNS (the dipole leg withheld, as D32.3; kcal/mol)", "=" * W, "", "   runs: " + "; ".join("%s %d" % kv for kv in sorted(cnt.items())) + " (of %d)" % len(RUNS), ""]
    for kind, d, s, o in DES:
        L.append("   %s, %s walled, behind the %s wall. D56's word: %s" % (d, "every charge" if kind == "A" else "positive charges", {"He": "helium", "Ne": "neon"}[s], surv[d]["the_design"]))
        for lv, _, _ in prep.LEVELS:
            x = o[lv]
            L.append("      %s  four legs %-15s fifth %-28s ddE %-9s with walls in place %s" % (lv, x["four"], x["fifth"][:28], f2(x["dde"]) or "-", f2(x["dde_wall"]) or "-"))
        L.append("      %s%s: agreement %s; at def2-TZVPD %s; BEHIND THE %s WALL: %s" % (d, " (every charge walled)" if kind == "A" else "", o["agree"], o["at_tz"], {"He": "HELIUM", "Ne": "NEON"}[s], o["word"]))
        L.append("")
    L += [P("Reported beside D56's bare reading and nothing more (A12.4 item 3). No design is named, kept or dropped on it."), ""]
    L += ["=" * W, "R57.2  THE WALL ALONE (centres with no charge; kcal/mol, e)", "=" * W, ""]
    for (t, d, s, lv), z in sorted(Z.items()):
        L.append("   %s %-9s %s %s  barrier change %-9s Loewdin rise R %s TS %s%s" % (t, d, s, lv, f2(z["dde"]) or "-", f2(z["rise"]["R"]) or "-", f2(z["rise"]["TS"]) or "-",
                 "   PERTURBS THE DENSITY" if any(v is not None and v >= M["audit"].D_LW for v in z["rise"].values()) else ""))
    L += ["", "=" * W, "R57.3  THE WALL AGAINST REAL IONS, site %d (kcal/mol)" % prep.SITE_I, "=" * W, ""]
    for lv in IO:
        x = IO[lv]
        L.append("   %s: ddE  bare %s  He %s  Ne %s  Shi %s | walls alone He %s  Ne %s  Shi %s" % (lv, *[f2(x["dde"].get(m)) or "-" for m in ("bare", "He", "Ne", "Shi", "He0", "Ne0", "Shi0")]))
        L.append("       Na+ %s  K+ %s  ghost Na %s  ghost K %s  -> ion(Na) %s  ion(K) %s%s" % (*[f2(x["dde"].get(m)) or "-" for m in ("Na", "K", "gNa", "gK")], f2(x["ion"]["Na"]) or "-", f2(x["ion"]["K"]) or "-",
                 "" if x["nelok"] else "   (an electron count is not the one expected)"))
        if "T" in x:
            L.append("       T %.4f, s %.4f, tolerance %.4f; distances: bare %s, He %s, Ne %s, Shi %s" % (x["T"], x["s"], x["tol"], *[f2(x["dist"][m]) or "-" for m in ("bare", "He", "Ne", "Shi")]))
            for m in ("He", "Ne", "Shi"):
                L.append("       %-3s %s" % (m, x["words"][m]))
        else:
            L.append("       NOT FORMED: a run is missing.")
        L.append("")
    L.append("   THE HELIUM WALL, by section 3 of the pre-registration: %s" % support_line)
    L += ["", "=" * W, "WHAT THIS CANNOT SAY (section 6 of the pre-registration, in short; its own words govern)", "=" * W, ""]
    L.append(P("R57.3 is one site, one frame, a whole charge, frozen geometries, in vacuo, and says nothing of a charge of 0.25 e. The walled designs are a model beside the bare one; no design is named a result."))
    L += ["", "Runs measured: %d of %d. THE HELIUM WALL: %s." % (cnt.get("measured", 0), len(RUNS), support_line)]
    txt = "\n".join(L) + "\n"
    with open(a.result, "w") as f:
        f.write(txt)
    print(txt, end="")


if __name__ == "__main__":
    main()
