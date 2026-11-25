#!/usr/bin/env python3
"""cs6_d60.py - the reader of STAGE 6, DIAGNOSTIC 60, committed before any of its outputs exists.

AMENDMENT 13, A13.4 items 7 and 8, as AMENDMENT 14 has them. It reads:
  R60.1 the control: D33 stage A's runs made again must give D33's committed energies and potentials, and the
        runs of wB97X-V must give D37's committed energies;
  R60.2 the map at wB97X-V/def2-SVPD, V_TS - V_R at the 2626 sites, written as A_f42436_d60_wB97XV.tsv, and how
        it and the committed B3LYP map compare with the working map over the sites that admit a charge;
  R60.3 D56's search made again by its committed program with the working map (which must give D56's 27
        designs), with the map at wB97X-V and with the B3LYP map, the designs side by side, and whether they are
        the same sites and charges;
  R60.4 the random designs of decision 10 under D56's rules, kinds (i) and (ii), 10 000 of each at each number of
        charges, and the place of each of D56's designs among them, at first order on the working map.
It names no design, keeps or drops none, and changes no word of D56.

usage (as cs6_d60_harvest.sh calls it):
   python3 cs6_d60.py --res <results repo> --code <code repo> --prereg <pre-registration>
                      --result <result document to write> --tsv <table of runs to write>
It also writes d60_designs.tsv, d60_null.tsv and A_f42436_d60_wB97XV.tsv beside --tsv. If a gate fails it
writes nothing.
"""
import argparse, contextlib, csv, hashlib, importlib.util, io, os, re, subprocess, sys, textwrap

PIN_PREREG = "847ced6dab913c3973643ab4c37a5961180fe4a188efcd5fee7ef840ed9682aa"
PIN_A14 = "7cedf94b4b74fc34febd900fbedb96d2c2f3ede5ae3215b32bcb73188f69d0ac"
S6REL = "phase2.2/cs_stage6"
DIR = "d60"
TOL_E = 1e-6                  # Eh: an energy made again is the committed one (AMENDMENT 14, item 3)
TOL_V = 1e-4                  # Eh/e: D33's committed limit on a potential (AMENDMENT 14, item 3)
NDRAW = 10000                 # random designs of each kind at each number of charges (decision 10)
SEED = 20261125               # the seed of the random designs, fixed here
NCHECK = 100                  # designs of kind (ii) at each n whose best charges are checked against HiGHS
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


def energy(d):
    """the last final energy of a run that ended normally, else None"""
    p = os.path.join(d, "job.out")
    if not os.path.isfile(p):
        return None
    t = rd(p)
    e = re.findall(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", t)
    return float(e[-1]) if e and "ORCA TERMINATED NORMALLY" in t else None


def main():
    ap = argparse.ArgumentParser()
    for a in ("res", "code", "prereg", "result", "tsv"):
        ap.add_argument("--" + a, required=True)
    a = ap.parse_args()
    import numpy as np
    from scipy.optimize import linprog
    from scipy.stats import spearmanr
    S6 = os.path.join(a.res, S6REL)
    D = os.path.join(S6, DIR)
    P2 = os.path.join(a.code, "phase2.2")
    NEED = ("cs6_d60_prepare.py", "cs6_d60.py", "cs6_d60_harvest.sh", "cs5_amatrix.py", "cs6_d56_prepare.py", "cs_opt_milp.py", "cs6_d47_prepare.py")
    # GV.1 the pre-registration and AMENDMENT 14; GV.2 the programs
    if not os.path.isfile(a.prereg) or sha(a.prereg) != PIN_PREREG:
        stop("the pre-registration given is not the committed STAGE6_DIAGNOSTIC60_20261125.txt.")
    a14 = os.path.join(S6, "STAGE6_AMENDMENT14_20261125.txt")
    if not os.path.isfile(a14) or sha(a14) != PIN_A14:
        stop("STAGE6_AMENDMENT14_20261125.txt is not the committed one.")
    for f in NEED:
        q = subprocess.run(["git", "-C", a.code, "diff", "--quiet", "HEAD", "--", "phase2.2/" + f])
        r = subprocess.run(["git", "-C", a.code, "cat-file", "-e", "HEAD:phase2.2/" + f], stderr=subprocess.DEVNULL)
        if q.returncode or r.returncode:
            stop("%s is not committed, or differs from its committed version." % f)
    prep = load(os.path.join(P2, "cs6_d60_prepare.py"), "d60prep")
    cs5 = load(os.path.join(P2, "cs5_amatrix.py"), "cs5_60")
    d56p = load(os.path.join(P2, "cs6_d56_prepare.py"), "d56p_60")
    opt = load(os.path.join(P2, "cs_opt_milp.py"), "opt_60")
    p47 = load(os.path.join(P2, "cs6_d47_prepare.py"), "p47_60")
    # GV.3 the inputs
    F = prep.build(S6)
    tracked = subprocess.run(["git", "-C", a.res, "ls-files", S6REL + "/" + DIR], stdout=subprocess.PIPE, check=True).stdout.decode().split("\n")
    tracked = set(os.path.relpath(t, S6REL + "/" + DIR) for t in tracked if t)
    if not set(F) <= tracked:
        stop("files of d60/ that the program forms are not committed.")
    for p in sorted(F):
        if rd(os.path.join(D, p)) != F[p]:
            stop("d60/%s is not the file the committed program forms from the committed files." % p)
    RUNS = prep.runs()
    # GV.4 echo
    for r, _ in RUNS:
        o = os.path.join(D, r, "job.out")
        if os.path.isfile(o):
            got = echoed_input(rd(o))
            if got and strip(got) != strip(rd(os.path.join(D, r, "job.inp")).split("\n")):
                stop("%s/job.out does not echo the committed job.inp beside it." % r)
    H = cs5.H
    pts = os.path.join(D, "points.xyz")

    # ------------------------------------------------------------ R60.1 the runs and the control
    status, E, V = {}, {}, {}
    for r, _ in RUNS:
        d = os.path.join(D, r)
        E[r] = energy(d)
        v = cs5.vpot(d, pts) if E[r] is not None else "the run did not end normally"
        V[r] = None if isinstance(v, str) else v
        status[r] = "NOT RUN" if not os.path.isfile(os.path.join(d, "job.out")) else ("measured" if E[r] is not None and V[r] is not None else
                    ("NOT MEASURED: no final energy" if E[r] is None else "NOT MEASURED: %s" % v))
    ref = {}
    for r, _ in RUNS:
        s = os.path.join(S6, prep.source(r))
        ref[r] = dict(E=energy(s), V=(cs5.vpot(s, os.path.join(S6, prep.POINTS)) if r.startswith("CAM") else None))
        if ref[r]["E"] is None or (r.startswith("CAM") and isinstance(ref[r]["V"], str)):
            stop("the committed run %s cannot be read." % prep.source(r))
    ctl = {}
    for r, _ in RUNS:
        x = dict(dE=None, dV=None)
        if E[r] is not None:
            x["dE"] = E[r] - ref[r]["E"]
        if r.startswith("CAM") and V[r] is not None:
            x["dV"] = float(np.abs(V[r] - ref[r]["V"]).max())
        ctl[r] = x
    cam_ok = all(ctl[r]["dE"] is not None and abs(ctl[r]["dE"]) <= TOL_E and ctl[r]["dV"] is not None and ctl[r]["dV"] <= TOL_V for r in ("CAM_R", "CAM_TS"))
    cam_formed = all(status[r] == "measured" for r in ("CAM_R", "CAM_TS"))
    wxv_ok = all(ctl[r]["dE"] is not None and abs(ctl[r]["dE"]) <= TOL_E for r in ("WXV_R", "WXV_TS"))
    wxv_pot = all(V[r] is not None for r in ("WXV_R", "WXV_TS"))
    wxv_formed = all(E[r] is not None for r in ("WXV_R", "WXV_TS"))
    CONTROL = ("THE CONTROL HOLDS" if cam_ok else "THE CONTROL DOES NOT HOLD") if cam_formed else "THE CONTROL IS NOT FORMED"
    WREP = ("REPRODUCES D37" if wxv_ok else "DOES NOT REPRODUCE D37") if wxv_formed else "NOT FORMED"
    USED = cam_ok and wxv_ok and wxv_pot

    # ------------------------------------------------------------ R60.2 the maps
    xyz, layer, ranges = opt.read_sites(os.path.join(S6, p47.SITES))
    nsite = len(xyz)
    a_cam = np.array(p47.read_col(os.path.join(S6, p47.AMAT)), dtype=float)
    a_b3 = np.asarray(opt.read_matrix(os.path.join(S6, "A_f42436.tsv"), nsite, "A")).ravel()
    a_wxv = (V["WXV_TS"] - V["WXV_R"]) * H if USED else None
    cut = d56p.cut_ranges(ranges)
    keep = [i for i, rr in enumerate(cut) if rr]
    if any(len(cut[i]) != 1 for i in keep):
        stop("a site admits two separate ranges of charge after the cut; the random designs assume one.")
    corr = {}
    for name, m in (("wB97X-V", a_wxv), ("B3LYP", a_b3)):
        if m is None:
            continue
        x, y = a_cam[keep], m[keep]
        corr[name] = dict(pearson=float(np.corrcoef(x, y)[0, 1]), spearman=float(spearmanr(x, y)[0]), mad=float(np.abs(y - x).mean()), mx=float(np.abs(y - x).max()))

    # ------------------------------------------------------------ R60.3 the searches
    D56 = list(csv.DictReader(open(os.path.join(S6, "d56", "DESIGNS.tsv")), delimiter="\t"))
    d56 = [dict(name=d["design"], n=int(d["n"]), k=int(d["rank_at_n"]), sites=[int(s) for s in d["sites"].split(",")],
                q=[float(v) for v in d["charges_e"].split(",")], score=float(d["first_order_score_kcal"])) for d in D56]
    for d in d56:
        s = sum(a_cam[i] * v for i, v in zip(d["sites"], d["q"]))
        if abs(s - d["score"]) > 5e-6:
            stop("the first-order score of %s on the working map is %.6f here and %.6f in DESIGNS.tsv." % (d["name"], s, d["score"]))
    key = lambda d: tuple(sorted(zip(d["sites"], [round(v, 4) for v in d["q"]])))
    print("the search with the working map (the check of the search; it may take 10 to 30 minutes) ...", flush=True)
    err = io.StringIO()
    try:
        with contextlib.redirect_stderr(err):
            s_cam, _ = d56p.search(S6, a.code, a_col=list(a_cam))
    except SystemExit:
        stop("D56's committed search program, with the working map, stopped: %s" % (" ".join(err.getvalue().split()) or "no message"))
    if [key(x) for x in s_cam] != [key(d) for d in d56] or [x["name"] for x in s_cam] != [d["name"] for d in d56]:
        stop("D56's committed search program, with the working map, does not give D56's 27 committed designs.")
    SEARCH, FAILED = {}, {}
    for name, m in (("wB97X-V", a_wxv), ("B3LYP", a_b3)):
        if m is None:
            continue
        print("the search with the map at %s (10 to 30 minutes) ..." % name, flush=True)
        err = io.StringIO()
        try:
            with contextlib.redirect_stderr(err):
                SEARCH[name], _ = d56p.search(S6, a.code, a_col=list(m))
        except SystemExit:
            FAILED[name] = " ".join(err.getvalue().split()) or "the search stopped"
    score = lambda d, m: float(sum(m[i] * v for i, v in zip(d["sites"], d["q"])))
    NS = sorted(set(d["n"] for d in d56))
    CMP = {}
    for name, ds in SEARCH.items():
        m = a_wxv if name == "wB97X-V" else a_b3
        rows = {}
        for n in NS:
            mine = [d for d in d56 if d["n"] == n]
            theirs = [d for d in ds if d["n"] == n]
            tk = set(key(d) for d in theirs)
            tsites = set(tuple(sorted(d["sites"])) for d in theirs)
            same = sum(1 for d in mine if key(d) in tk)
            same_sites = sum(1 for d in mine if tuple(sorted(d["sites"])) in tsites)
            share = len(set(key(mine[0])) & set(key(theirs[0]))) / float(n)
            word = "THE SAME DESIGNS" if same == len(mine) else ("PARTLY THE SAME" if same else "DIFFERENT DESIGNS")
            rows[n] = dict(same=same, same_sites=same_sites, share=share, word=word, of=len(mine),
                           best_theirs_own=score(theirs[0], m), best_ours_theirs=score(mine[0], m), best_theirs_cam=score(theirs[0], a_cam))
        CMP[name] = rows

    # ------------------------------------------------------------ R60.4 the random designs
    rng = np.random.Generator(np.random.PCG64(SEED))
    K = np.array(keep)
    P = xyz[K]
    close = np.zeros((len(K), len(K)), dtype=bool)
    for i, j in opt.close_pairs(P, d56p.R_MIN):        # the pairs D56's search forbids, by the committed function
        close[i, j] = close[j, i] = True
    lo = np.array([min(r[0] for r in cut[i]) for i in keep])
    hi = np.array([max(r[1] for r in cut[i]) for i in keep])
    cap = d56p.CAP_Q
    plus_ok = hi >= cap - 1e-9
    minus_ok = lo <= -cap + 1e-9
    ak = a_cam[K]

    def sites_draw(n):
        """random sequential addition over the sites that admit a charge, respecting the spacing"""
        while True:
            chosen, blocked = [], np.zeros(len(K), dtype=bool)
            for j in rng.permutation(len(K)):
                if not blocked[j]:
                    chosen.append(j)
                    blocked |= close[j]
                    if len(chosen) == n:
                        return np.array(chosen)

    def best(ch):
        """the best charges D56's rules allow at these sites: min sum a q, lo <= q <= hi, sum q >= 0 (exact, by order)"""
        a_ = ak[ch]
        q = np.where(a_ > 0, lo[ch], hi[ch]).astype(float)
        deficit = -q.sum()
        if deficit > 0:
            for j in np.argsort(a_):
                room = hi[ch][j] - q[j]
                if room <= 0:
                    continue
                step = min(room, deficit)
                q[j] += step
                deficit -= step
                if deficit <= 1e-15:
                    break
        return q

    NULL, CHECK = {}, []
    for n in NS:
        s1, s2, tries1, odd = [], [], 0, 0
        while len(s1) < NDRAW:
            ch = sites_draw(n)
            tries1 += 1
            sg = np.array([rng.choice([v for v, ok in ((cap, plus_ok[j]), (-cap, minus_ok[j])) if ok]) for j in ch])
            if sg.sum() < -1e-9:
                continue
            s1.append(float((ak[ch] * sg).sum()))
        for t in range(NDRAW):
            ch = sites_draw(n)
            q = best(ch)
            odd += int(np.any(np.abs(np.abs(q) - cap) > 1e-9))
            s2.append(float((ak[ch] * q).sum()))
            if t < NCHECK:
                res = linprog(ak[ch], A_ub=-np.ones((1, n)), b_ub=[0.0], bounds=list(zip(lo[ch], hi[ch])), method="highs")
                CHECK.append(abs(res.fun - s2[-1]) if res.status == 0 else float("inf"))
        NULL[n] = dict(i=np.array(s1), ii=np.array(s2), acc=NDRAW / float(tries1), odd=odd)
    if max(CHECK) > 1e-9:
        stop("the best charges of kind (ii) differ from HiGHS's by %.3g kcal/mol on a checked design." % max(CHECK))
    PCT = {}
    for d in d56:
        z = NULL[d["n"]]
        PCT[d["name"]] = dict(i=100.0 * float((z["i"] > d["score"]).mean()), ii=100.0 * float((z["ii"] > d["score"]).mean()))

    # ------------------------------------------------------------ the tables
    base = os.path.dirname(a.tsv)
    f6 = lambda v: "" if v is None else "%.6f" % v
    with open(a.tsv, "w") as f:
        f.write("run\tjob\tstatus\tenergy_Eh\tcommitted_source\tminus_committed_Eh\tlargest_potential_difference_Eh_per_e\n" + "".join(
            "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (r, prep.JOB, status[r], "" if E[r] is None else "%.9f" % E[r], prep.source(r),
                                              "" if ctl[r]["dE"] is None else "%.3e" % ctl[r]["dE"], "" if ctl[r]["dV"] is None else "%.3e" % ctl[r]["dV"]) for r, _ in RUNS))
    with open(os.path.join(base, "d60_designs.tsv"), "w") as f:
        f.write("map\tdesign\tn\trank_at_n\tsites\tcharges_e\tscore_working_map_kcal\tscore_wB97X-V_map_kcal\tscore_B3LYP_map_kcal\tsame_as_a_D56_design\n")
        allsets = [("working (D56)", s_cam)] + [(k, v) for k, v in SEARCH.items()]
        d56k = set(key(d) for d in d56)
        for name, ds in allsets:
            for d in ds:
                f.write("\t".join([name, d["name"], str(d["n"]), str(d["k"]), ",".join(str(s) for s in d["sites"]), ",".join("%+.4f" % v for v in d["q"]),
                                   f6(score(d, a_cam)), f6(score(d, a_wxv)) if a_wxv is not None else "", f6(score(d, a_b3)), "yes" if key(d) in d56k else "no"]) + "\n")
    with open(os.path.join(base, "d60_null.tsv"), "w") as f:
        f.write("design\tn\tscore_working_map_kcal\tbeats_kind_i_percent\tbeats_kind_ii_percent\tkind_i_median\tkind_i_lowest\tkind_ii_median\tkind_ii_lowest\tkind_i_acceptance\n")
        for d in d56:
            z = NULL[d["n"]]
            f.write("\t".join([d["name"], str(d["n"]), f6(d["score"]), "%.2f" % PCT[d["name"]]["i"], "%.2f" % PCT[d["name"]]["ii"],
                               f6(float(np.median(z["i"]))), f6(float(z["i"].min())), f6(float(np.median(z["ii"]))), f6(float(z["ii"].min())), "%.4f" % z["acc"]]) + "\n")
    if a_wxv is not None:
        with open(os.path.join(base, "A_f42436_d60_wB97XV.tsv"), "w") as f:
            f.write("# a_if = V_TS(r_i) - V_R(r_i), kcal/mol per e, in vacuo wB97X-V/def2-SVPD at the\n# aligned geometries of frame 42436. DIAGNOSTIC 60, cs6_d60.py (AMENDMENT 14).\n")
            f.write("# row i = site i of sites_f42436.tsv (sha256 %s)\n" % sha(os.path.join(S6, p47.SITES)))
            f.write("".join("%.6f\n" % v for v in a_wxv))

    # ------------------------------------------------------------ the document
    Pw = lambda s, ind=0: "\n".join(textwrap.wrap(" ".join(s.split()), W, initial_indent=" " * ind, subsequent_indent=" " * ind, break_long_words=False))
    head = subprocess.run(["git", "-C", a.res, "rev-parse", "HEAD"], stdout=subprocess.PIPE).stdout.decode().strip()
    L = ["STAGE 6 - DIAGNOSTIC 60 - RESULT", "D56's search at a second functional, and against random designs", "=" * W, ""]
    L.append(Pw("Read by cs6_d60.py under STAGE6_DIAGNOSTIC60_20261125.txt and AMENDMENT 14, committed with the inputs before any job was submitted."))
    L.append("   results repository at %s" % head[:12])
    for f in NEED[:2] + NEED[3:]:
        L.append("   %-22s sha256 %s" % (f, sha(os.path.join(P2, f))))
    L += ["", Pw("THE GATES, as they held: GV.1 the pre-registration and AMENDMENT 14 are the committed ones; GV.2 the programs are the committed ones; GV.3 every input of d60/ is the file the committed program forms from the committed files; GV.4 every output echoes its committed input; GV.5 every first-order score of D56's DESIGNS.tsv is the working map's; GV.6 D56's committed search program, with the working map, gives D56's 27 committed designs; GV.7 the best charges of kind (ii) agree with HiGHS on %d designs at each number of charges (largest difference %.1e kcal/mol)." % (NCHECK, max(CHECK))), ""]
    L += ["=" * W, "R60.1  THE RUNS AND THE CONTROL", "=" * W, ""]
    for r, _ in RUNS:
        L.append("   %-7s %-28s E %-18s minus %s: %s%s" % (r, status[r][:28], "-" if E[r] is None else "%.9f" % E[r], prep.source(r),
                 "-" if ctl[r]["dE"] is None else "%+.2e Eh" % ctl[r]["dE"], "" if ctl[r]["dV"] is None else "; largest potential difference %.2e Eh/e" % ctl[r]["dV"]))
    L += ["", "   %s (energies within 1e-6 Eh, potentials within 1e-4 Eh/e of D33's committed ones)." % CONTROL,
          "   THE RUNS OF wB97X-V: %s (energies within 1e-6 Eh of D37's committed ones)." % WREP,
          "   THE MAP AT wB97X-V IS %s." % ("USED" if USED else "NOT USED: " + "; ".join(w for w, ok in (("the control does not hold", cam_ok), ("the runs of wB97X-V do not reproduce D37", wxv_ok), ("a run of wB97X-V has no potentials", wxv_pot)) if not ok)), ""]
    L += ["=" * W, "R60.2  THE MAPS, over the %d sites that admit a charge (kcal/mol per e)" % len(keep), "=" * W, ""]
    for name in ("wB97X-V", "B3LYP"):
        if name in corr:
            c = corr[name]
            L.append("   %-8s against the working map: Pearson %.4f, Spearman %.4f, mean |difference| %.4f, largest %.4f%s" % (
                name, c["pearson"], c["spearman"], c["mad"], c["mx"], "   (def2-SVP, not def2-SVPD)" if name == "B3LYP" else ""))
        else:
            L.append("   %-8s not formed" % name)
    L += ["", "=" * W, "R60.3  THE SEARCH AT EACH MAP, BESIDE D56's (same: the same sites and the same charges)", "=" * W, ""]
    for name in ("wB97X-V", "B3LYP"):
        if name in FAILED:
            L += ["   %s: THE SEARCH STOPPED: %s" % (name, FAILED[name]), ""]
            continue
        if name not in CMP:
            L += ["   %s: not run, as its map is not used." % name, ""]
            continue
        L.append("   %s%s" % (name, " (its basis is def2-SVP, the working map's def2-SVPD; the two maps differ in the basis set too)" if name == "B3LYP" else ""))
        L.append("      n   same   same sites   charges the two best share   its best on its map   D56's best on its map   its best on the working map")
        tot = 0
        for n in NS:
            x = CMP[name][n]
            tot += x["same"]
            L.append("     %2d   %d of %d   %d of %d         %5.2f            %10.4f            %10.4f              %10.4f   %s" % (
                n, x["same"], x["of"], x["same_sites"], x["of"], x["share"], x["best_theirs_own"], x["best_ours_theirs"], x["best_theirs_cam"], x["word"]))
        words = [CMP[name][n]["word"] for n in NS]
        L.append("      %s: THE SAME DESIGNS at %d of %d numbers of charges, PARTLY THE SAME at %d, DIFFERENT DESIGNS at %d; %d of D56's 27 designs found again." % (
            name, words.count("THE SAME DESIGNS"), len(NS), words.count("PARTLY THE SAME"), words.count("DIFFERENT DESIGNS"), tot))
        L.append("")
    L += ["=" * W, "R60.4  D56's DESIGNS AGAINST RANDOM DESIGNS (decision 10; first order, the working map; kcal/mol)", "=" * W, ""]
    L.append(Pw("Kind (i): random sites with the spacing, each charge 0.25 e with a random sign its site admits, drawn again if the charges sum below zero. Kind (ii): random sites with the best charges D56's rules allow there. %s of each at each number of charges; seed %d. 'beats' is the share of random designs whose first-order score is above the design's (a smaller lowering)." % ("{:,}".format(NDRAW).replace(",", " "), SEED)))
    L.append("")
    L.append("   design      score     beats (i)  beats (ii)   (i) median  (i) lowest  (ii) median (ii) lowest")
    for d in d56:
        z = NULL[d["n"]]
        L.append("   %-10s %9.4f   %7.2f%%   %7.2f%%   %9.4f  %9.4f   %9.4f  %9.4f" % (d["name"], d["score"], PCT[d["name"]]["i"], PCT[d["name"]]["ii"],
                 float(np.median(z["i"])), float(z["i"].min()), float(np.median(z["ii"])), float(z["ii"].min())))
    L.append("")
    L.append(Pw("Acceptance of kind (i), the share of drawn designs whose charges sum to zero or more: " + ", ".join("n %d %.3f" % (n, NULL[n]["acc"]) for n in NS) + "."))
    L.append(Pw("Designs of kind (ii) whose best charges include one smaller than 0.25 e in size (zero included), as none of D56's designs holds: " + ", ".join("n %d %d" % (n, NULL[n]["odd"]) for n in NS) + ", of 10 000 at each."))
    L += ["", "=" * W, "WHAT THIS CANNOT SAY (section 5 of the pre-registration, in short; its own words govern)", "=" * W, ""]
    L.append(Pw("First order on one frame, in vacuo. A design found again at another functional is placed by that map, not measured there; a random design is scored, not measured. The 50 full single points of decision 10 are not run. Nothing here names, keeps or drops a design."))
    nb1 = sum(1 for d in d56 if PCT[d["name"]]["i"] >= 99.0)
    nb2 = sum(1 for d in d56 if PCT[d["name"]]["ii"] >= 99.0)
    L += ["", Pw("%s. THE MAP AT wB97X-V IS %s. %s D56's designs beating 99%% or more of the random designs: %d of 27 of kind (i), %d of 27 of kind (ii)." % (
        CONTROL, "USED" if USED else "NOT USED", " ".join("At %s, %d of D56's 27 designs found again." % (k, sum(CMP[k][n]["same"] for n in NS)) for k in ("wB97X-V", "B3LYP") if k in CMP)
        + "".join(" At %s, THE SEARCH STOPPED." % k for k in FAILED), nb1, nb2))]
    txt = "\n".join(L) + "\n"
    with open(a.result, "w") as f:
        f.write(txt)
    print(txt, end="")


if __name__ == "__main__":
    main()
