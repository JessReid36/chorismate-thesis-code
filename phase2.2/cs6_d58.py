#!/usr/bin/env python3
"""cs6_d58.py - the reader of STAGE 6, DIAGNOSTIC 58, committed before any of its outputs exists.

AMENDMENT 13, A13.4 items 2 to 4. It reads:
  1. the gradient each of D56's 18 survivors adds to the substrate at R and at TS (G runs), as
     STAGE6_CRITERIA.txt's G_add (cs6_analyse.g_add), with the largest component, and beside it
     the part the sphere's Lennard-Jones shell would add (cs5_grid's committed function);
  2. the first-order map F of the added gradient, from 144 displaced runs (M runs) and the
     potentials at the 2626 sites, and its check against item 1 at R;
  3. the first-order map E of the change in the energy to remove the last electron, from the
     four controls (E runs) and their potentials, and its check against D56's 54 measured values.
It names no design a result and changes no word of D56. Each map is written to a committed table
with its word; the second search may use a map only if its word is PASSES.

usage (as cs6_d58_harvest.sh calls it):
   python3 cs6_d58.py --res <results repo> --code <code repo> --prereg <pre-registration>
                      --result <result document to write> --tsv <table of runs to write>
It also writes, beside --tsv: d58_forces.tsv, and F_f42436_d58.tsv / E_f42436_d58.tsv when formed.
If a gate fails it prints why and writes nothing.
"""
import argparse, csv, hashlib, importlib.util, os, re, subprocess, sys, textwrap
import numpy as np

PIN_PREREG = "9023bdcaa547af06a130f7beeded1e13bbfa7c2bbd448443cb8f038eebbc17b4"
S6REL = "phase2.2/cs_stage6"
DIR = "d58"
EV = 27.211386245988            # eV per Eh, as cs6_d56.py
T_F = 0.25                      # map F: relative tolerance (pre-registration, section 3)
T_E = 0.15                      # map E: eV (pre-registration, section 3)
TOL_REPRO = 1.0e-6              # Eh: a control reproduces D56's
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
    sys.modules[name] = m
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


def final_energy(t):
    E = re.findall(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", t)
    return float(E[-1]) if E else None


def lj_grad(Q, X, eps, rmin, G):
    """d E_LJ / d X (kcal/mol/A), (24, 3), of spheres at the rows of Q on the substrate atoms X;
    the analytic derivative of cs5_grid.lj_energy"""
    Rij = 0.5 * (2 * G.SPHERE_RSTAR + rmin)
    eij = np.sqrt(G.SPHERE_EPS * eps)
    g = np.zeros_like(X)
    for q in Q:
        d = X - q[None, :]
        r = np.sqrt((d ** 2).sum(1))
        s6 = (Rij / r) ** 6
        dEdr = eij * (-12.0 * s6 * s6 + 12.0 * s6) / r
        g += (dEdr / r)[:, None] * d
    return g


def main():
    ap = argparse.ArgumentParser()
    for a in ("res", "code", "prereg", "result", "tsv"):
        ap.add_argument("--" + a, required=True)
    a = ap.parse_args()
    S6 = os.path.join(a.res, S6REL)
    D = os.path.join(S6, DIR)
    P2 = os.path.join(a.code, "phase2.2")
    NEED = ("cs6_d58_prepare.py", "cs6_d58.py", "cs6_d58_harvest.sh", "cs6_analyse.py", "cs5_grid.py", "cs5_amatrix.py")

    # GD.1 the pre-registration; GD.2 the programs
    if not os.path.isfile(a.prereg) or sha(a.prereg) != PIN_PREREG:
        stop("the pre-registration given is not the committed STAGE6_DIAGNOSTIC58_20261125.txt.")
    for f in NEED:
        q = subprocess.run(["git", "-C", a.code, "diff", "--quiet", "HEAD", "--", "phase2.2/" + f])
        r = subprocess.run(["git", "-C", a.code, "cat-file", "-e", "HEAD:phase2.2/" + f], stderr=subprocess.DEVNULL)
        if q.returncode or r.returncode:
            stop("%s is not committed, or differs from its committed version." % f)
    sys.path.insert(0, P2)
    prep = load(os.path.join(P2, "cs6_d58_prepare.py"), "d58prep")
    G = load(os.path.join(P2, "cs5_grid.py"), "cs5_grid")
    AM = load(os.path.join(P2, "cs5_amatrix.py"), "cs5_amatrix")
    AN = load(os.path.join(P2, "cs6_analyse.py"), "cs6_analyse")
    # GD.3 the committed inputs are what the committed program forms
    F = prep.build(S6)
    tracked = subprocess.run(["git", "-C", a.res, "ls-files", S6REL + "/" + DIR], stdout=subprocess.PIPE, check=True).stdout.decode().split("\n")
    tracked = set(os.path.relpath(t, S6REL + "/" + DIR) for t in tracked if t)
    if not set(F) <= tracked:
        stop("files of d58/ that the program forms are not committed.")
    for p in sorted(F):
        if rd(os.path.join(D, p)) != F[p]:
            stop("d58/%s is not the file the committed program forms from the committed files." % p)
    RUNS, JOBS = prep.runs(S6)
    # GD.4 every output echoes the committed input beside it
    for r, _, _ in RUNS:
        o = os.path.join(D, r, "job.out")
        if os.path.isfile(o):
            got = echoed_input(rd(o))
            if got and strip(got) != strip(rd(os.path.join(D, r, "job.inp")).split("\n")):
                stop("%s/job.out does not echo the committed job.inp beside it." % r)
    # GD.5 the committed Lennard-Jones function, with D56's geometries, gives the committed L_f42436.tsv
    eps, rmin = G.load_lj(os.path.join(a.res, "05_qmmm/13_bridge/complex_solvated.ORCAFF.prms"))
    xyz = lambda p: np.array([[float(v) for v in l.split()[1:4]] for l in rd(p).split("\n")[2:26]])
    X = {"R": xyz(os.path.join(D, "G_C_R/substrate.xyz")), "TS": xyz(os.path.join(D, "G_C_TS/substrate.xyz"))}
    sites = []
    for l in rd(os.path.join(S6, "sites_f42436.tsv")).split("\n"):
        c = l.split("\t")
        if l.startswith("#") or len(c) < 5:
            continue
        try:
            sites.append([float(v) for v in c[2:5]])
        except ValueError:
            pass
    Qs = np.array(sites)
    Lc = np.array([float(l) for l in rd(os.path.join(S6, "L_f42436.tsv")).split("\n") if l.strip() and not l.startswith("#")])
    if len(Qs) != 2626 or np.abs(G.lj_energy(Qs, X["TS"], eps, rmin) - G.lj_energy(Qs, X["R"], eps, rmin) - Lc).max() > 1e-5:
        stop("the committed Lennard-Jones function does not give the committed L_f42436.tsv.")
    # GD.6 the analytic shell gradient is the derivative of that function (central differences, one design)
    des = list(csv.DictReader(open(os.path.join(S6, "d56/DESIGNS.tsv")), delimiter="\t"))
    site_of = {r["design"]: [int(v) for v in r["sites"].split(",")] for r in des}
    q_of = {r["design"]: [float(v) for v in r["charges_e"].split(",")] for r in des}
    Qt = Qs[site_of["S_n10_k1"]]
    ga = lj_grad(Qt, X["R"], eps, rmin, G)
    gn = np.zeros_like(ga)
    for k in range(24):
        for j in range(3):
            Xp, Xm = X["R"].copy(), X["R"].copy()
            Xp[k, j] += 1e-5
            Xm[k, j] -= 1e-5
            gn[k, j] = (G.lj_energy(Qt, Xp, eps, rmin).sum() - G.lj_energy(Qt, Xm, eps, rmin).sum()) / 2e-5
    if np.abs(ga - gn).max() > 1e-4:
        stop("the analytic gradient of the shell is not the derivative of the committed function.")

    # the job logs
    exitc, vex = {}, {}
    for j in JOBS:
        lp = os.path.join(D, "run_%s.pbs.out" % j)
        if os.path.isfile(lp):
            for l in rd(lp).split("\n"):
                m = re.match(r"^(\S+) orca_exit=(\d+) normal=(\d+)", l)
                if m:
                    exitc[m.group(1)] = int(m.group(2))
                m = re.match(r"^(\S+) vpot_exit=(\d+)", l)
                if m:
                    vex[m.group(1)] = int(m.group(2))
    # D56's committed numbers
    d56 = {row["run"]: row for row in csv.DictReader(open(os.path.join(S6, "d56.tsv")), delimiter="\t")}
    dd = {(row["design"], row["level"]): row for row in csv.DictReader(open(os.path.join(S6, "d56_designs.tsv")), delimiter="\t")}

    # R58.1 the runs
    rows, text = [], {}
    for r, j, cap in RUNS:
        o = os.path.join(D, r, "job.out")
        t = rd(o) if os.path.isfile(o) else ""
        text[r] = t
        normal = "ORCA TERMINATED NORMALLY" in t
        st = "NOT RUN" if not t else ("ENDED NORMALLY" if normal else ("STOPPED AT ITS CAP" if exitc.get(r) == 124 else "DID NOT END NORMALLY"))
        E = final_energy(t) if normal else None
        v = ""
        if r[0] in "ME":
            vv = AM.vpot(os.path.join(D, r), os.path.join(D, "points.xyz")) if normal else "not run"
            v = "2626 POTENTIALS" if not isinstance(vv, str) else "NO POTENTIALS (%s)" % vv
        rows.append(dict(run=r, job=j, status=st, orca_exit=exitc.get(r, ""), E_Eh="" if E is None else "%.10f" % E, potentials=v))
    by = {x["run"]: x for x in rows}
    repro = {}
    for r, src in (("G_C_R", "SV_C_ctl_R"), ("G_C_TS", "SV_C_ctl_TS"), ("E_ctl_R", "SV_C_ctl_R"), ("E_ctl_TS", "SV_C_ctl_TS"),
                   ("E_bare_R", "SV_C_bare_R"), ("E_bare_TS", "SV_C_bare_TS")):
        e = by[r]["E_Eh"]
        repro[r] = None if e == "" else float(e) - float(d56[src]["E_Eh"])

    # R58.2 the added gradients
    def eng(r):
        p = os.path.join(D, r, "job.engrad")
        if by[r]["status"] != "ENDED NORMALLY" or not os.path.isfile(p):
            return None
        n, e, g, Z, Xb = AN.read_engrad(p)
        if n != 24 or np.abs(Xb * AN.BOHR - X["R" if r.endswith("_R") else "TS"]).max() > 1e-4:
            return None
        return g
    surv = prep.survivors(S6)
    g0 = {"R": eng("G_C_R"), "TS": eng("G_C_TS")}
    frows, direct = [], {}
    for d in surv:
        for stt in ("R", "TS"):
            g = eng("G_%s_%s" % (d, stt))
            row = dict(design=d, state=stt, G_add_kcal_per_A="", largest_component="", G_add_shell_only="", G_add_with_shell="")
            if g is not None and g0[stt] is not None:
                dg = (g - g0[stt]) * AN.GCONV
                direct[(d, stt)] = g - g0[stt]
                sh = lj_grad(Qs[site_of[d]], X[stt], eps, rmin, G)
                row.update(G_add_kcal_per_A="%.4f" % AN.g_add(g, g0[stt]), largest_component="%.4f" % np.abs(dg).max(),
                           G_add_shell_only="%.4f" % float(np.sqrt((sh ** 2).mean())), G_add_with_shell="%.4f" % float(np.sqrt(((dg + sh) ** 2).mean())))
            frows.append(row)

    # R58.3 the map F
    Fmap, fword, fchk = None, "NOT FORMED", []
    Mruns = [r for r, _, _ in RUNS if r.startswith("M_")]
    if all(by[r]["potentials"] == "2626 POTENTIALS" for r in Mruns):
        h = prep.STEP_BOHR
        Fmap = np.zeros((72, 2626))
        for k in range(24):
            for jx, ax in enumerate("xyz"):
                vp = AM.vpot(os.path.join(D, "M_a%02d_%s_p" % (k, ax)), os.path.join(D, "points.xyz"))
                vm = AM.vpot(os.path.join(D, "M_a%02d_%s_m" % (k, ax)), os.path.join(D, "points.xyz"))
                Fmap[3 * k + jx] = (vp - vm) / (2 * h)              # Eh per e per bohr
        ok = True
        for d in surv:
            if (d, "R") not in direct:
                ok = None
                fchk.append((d, None, None))
                continue
            p = Fmap[:, site_of[d]] @ np.array(q_of[d])
            dg = direct[(d, "R")].reshape(-1)
            rel = float(np.sqrt(((p - dg) ** 2).mean()) / np.sqrt((dg ** 2).mean()))
            mx = float(abs(np.abs(p).max() - np.abs(dg).max()) / np.abs(dg).max())
            fchk.append((d, rel, mx))
            if rel > T_F or mx > T_F:
                ok = False if ok is not None else None
        fword = "NOT FORMED" if ok is None else ("PASSES" if ok else "FAILS")

    # R58.4 the map E
    Emap, eword, echk, vde0 = None, "NOT FORMED", [], {}
    Eruns = ("E_ctl_R", "E_ctl_TS", "E_bare_R", "E_bare_TS")
    if all(by[r]["potentials"] == "2626 POTENTIALS" for r in Eruns):
        if all(repro[r] is not None and abs(repro[r]) <= TOL_REPRO for r in Eruns):
            V = {r: AM.vpot(os.path.join(D, r), os.path.join(D, "points.xyz")) for r in Eruns}
            Emap = {st: (V["E_bare_%s" % st] - V["E_ctl_%s" % st]) * EV for st in ("R", "TS")}      # eV per e
            for st in ("R", "TS"):
                vde0[st] = (float(by["E_bare_%s" % st]["E_Eh"]) - float(by["E_ctl_%s" % st]["E_Eh"])) * EV
            errs = []
            for r in des:
                d = r["design"]
                for st in ("R", "TS"):
                    meas = dd[(d, "SV")]["VDE_%s_eV" % st].strip()
                    if meas in ("", "-"):
                        continue
                    pred = vde0[st] + float(Emap[st][site_of[d]] @ np.array(q_of[d]))
                    echk.append((d, st, pred, float(meas)))
                    errs.append(pred - float(meas))
            eword = "PASSES" if len(errs) == 54 and max(abs(e) for e in errs) <= T_E else ("FAILS" if len(errs) == 54 else "NOT FORMED")
        else:
            eword = "NOT FORMED: A CONTROL DOES NOT REPRODUCE D56"
    d33 = AM.vpot(os.path.join(S6, "d33a/sp_42436_R"), os.path.join(S6, "d33a/points.xyz"))
    vE = AM.vpot(os.path.join(D, "E_ctl_R"), os.path.join(D, "points.xyz")) if by["E_ctl_R"]["potentials"] == "2626 POTENTIALS" else "no"
    d33diff = None if isinstance(vE, str) or isinstance(d33, str) else float(np.abs(vE - d33).max())

    # the tables
    base = os.path.dirname(a.tsv)
    cols = ["run", "job", "status", "orca_exit", "E_Eh", "potentials"]
    with open(a.tsv, "w") as f:
        f.write("\t".join(cols) + "\n" + "".join("\t".join(str(x[c]) for c in cols) + "\n" for x in rows))
    fc = ["design", "state", "G_add_kcal_per_A", "largest_component", "G_add_shell_only", "G_add_with_shell"]
    with open(os.path.join(base, "d58_forces.tsv"), "w") as f:
        f.write("\t".join(fc) + "\n" + "".join("\t".join(x[c] for c in fc) + "\n" for x in frows))
    if Fmap is not None:
        with open(os.path.join(base, "F_f42436_d58.tsv"), "w") as f:
            f.write("# map F, cs6_d58.py (DIAGNOSTIC 58). Row i = site i of sites_f42436.tsv. Column 3k+j = atom k (D56's substrate order), axis j.\n")
            f.write("# F[i,c] = d V(r_i) / d X_c of the substrate without charges at R, frame 42436, CAM-B3LYP-D3BJ/def2-SVPD, by central\n")
            f.write("#   differences of +-%.2f bohr; kcal/mol per A per e. The first-order added GRADIENT of a design is sum_i q_i F[i, :].\n" % prep.STEP_BOHR)
            f.write("# word of the map: %s (tolerance %.2f, pre-registration section 3)\n" % (fword, T_F))
            for i in range(2626):
                f.write("\t".join("%.6f" % v for v in Fmap[:, i] * AN.GCONV) + "\n")
    if Emap is not None:
        with open(os.path.join(base, "E_f42436_d58.tsv"), "w") as f:
            f.write("# map E, cs6_d58.py (DIAGNOSTIC 58). Row i = site i of sites_f42436.tsv.\n")
            f.write("# E[i] = V_one electron removed(r_i) - V_dianion(r_i), frame 42436, CAM-B3LYP-D3BJ/def2-SVPD; eV per e. To first order the energy\n")
            f.write("#   to remove the last electron with a design is VDE0 + sum_i q_i E[i]. VDE0 R %.4f eV, TS %.4f eV (this diagnostic's controls).\n" % (vde0["R"], vde0["TS"]))
            f.write("# word of the map: %s (tolerance %.2f eV, pre-registration section 3)\nE_R_eV_per_e\tE_TS_eV_per_e\n" % (eword, T_E))
            for i in range(2626):
                f.write("%.6f\t%.6f\n" % (Emap["R"][i], Emap["TS"][i]))

    # the document
    P = lambda s, ind=0: "\n".join(textwrap.wrap(" ".join(s.split()), W, initial_indent=" " * ind, subsequent_indent=" " * ind, break_long_words=False))
    head = subprocess.run(["git", "-C", a.res, "rev-parse", "HEAD"], stdout=subprocess.PIPE).stdout.decode().strip()
    L = ["STAGE 6 - DIAGNOSTIC 58 - RESULT", "The force each survivor of D56 adds, and the first-order maps of the force and of the hold on the", "last electron", "=" * W, ""]
    L.append(P("Read by cs6_d58.py under STAGE6_DIAGNOSTIC58_20261125.txt, which was committed with the inputs before any job was submitted."))
    L.append("   results repository at %s" % head[:12])
    for f in NEED[:2] + NEED[3:]:
        L.append("   %-22s sha256 %s" % (f, sha(os.path.join(P2, f))))
    L += ["", P("THE GATES, as they held: GD.1 the pre-registration is the committed one; GD.2 the programs are the committed ones; GD.3 every input of d58/ is the file the committed program forms from the committed files; GD.4 every output that exists echoes the committed input beside it; GD.5 the committed Lennard-Jones function, with D56's geometries, gives the committed L_f42436.tsv; GD.6 the analytic gradient of the shell is the derivative of that function."), ""]
    cnt = {}
    for x in rows:
        cnt[x["status"]] = cnt.get(x["status"], 0) + 1
    L += ["=" * W, "R58.1  THE RUNS", "=" * W, "", "   " + "; ".join("%s %d" % (k, v) for k, v in sorted(cnt.items())) + " (of 186)"]
    bad = [x["run"] for x in rows if x["status"] != "ENDED NORMALLY"]
    if bad:
        L.append(P("not ended normally: " + ", ".join(bad), 3))
    nopot = [x["run"] for x in rows if x["potentials"].startswith("NO")]
    if nopot:
        L.append(P("no potentials: " + ", ".join(nopot), 3))
    L.append("")
    for r in ("G_C_R", "G_C_TS", "E_ctl_R", "E_ctl_TS", "E_bare_R", "E_bare_TS"):
        L.append("   %-10s against D56's committed control: %s" % (r, "not read" if repro[r] is None else "%s (%+.1e Eh)" % ("REPRODUCES" if abs(repro[r]) <= TOL_REPRO else "DIFFERS", repro[r])))
    L.append("   the potentials of E_ctl_R against D33's committed sp_42436_R: %s" % ("not read" if d33diff is None else "largest difference %.1e Eh/e" % d33diff))
    L += ["", "=" * W, "R58.2  THE GRADIENT EACH SURVIVOR ADDS (kcal/mol/A; G_add as STAGE6_CRITERIA.txt defines it)", "=" * W, ""]
    L.append("   design     state  G_add    largest   shell alone  G_add with shell")
    for x in frows:
        L.append("   %-10s %-5s  %-8s %-9s %-12s %s" % (x["design"], x["state"], x["G_add_kcal_per_A"] or "-", x["largest_component"] or "-", x["G_add_shell_only"] or "-", x["G_add_with_shell"] or "-"))
    L += ["", P("Reported; decides nothing (AMENDMENT 13, A13.4 item 2). The shell's part is the gradient the sphere's Lennard-Jones term would add, from the committed parameters; these runs carry bare charges."), ""]
    L += ["=" * W, "R58.3  THE MAP F AND ITS CHECK AT R", "=" * W, ""]
    L.append(P("relative RMS = RMS(map - direct) / RMS(direct) over the 72 components; largest = |largest |map| - largest |direct|| / largest |direct|. Tolerance %.2f for both, for every survivor." % T_F, 3))
    for d, rel, mx in fchk:
        L.append("   %-10s %s" % (d, "not read" if rel is None else "relative RMS %.3f   largest %.3f%s" % (rel, mx, "" if rel <= T_F and mx <= T_F else "   OUTSIDE")))
    L += ["", "   THE MAP F: %s" % fword, ""]
    L += ["=" * W, "R58.4  THE MAP E AND ITS CHECK AGAINST D56", "=" * W, ""]
    if echk:
        e = [p - m for _, _, p, m in echk]
        L.append(P("VDE0 from this diagnostic's controls: R %+.4f eV, TS %+.4f eV. Over D56's %d measured values at the working level: largest |predicted - measured| %.3f eV, mean predicted - measured %+.3f eV. Tolerance %.2f eV." % (vde0["R"], vde0["TS"], len(e), max(abs(x) for x in e), sum(e) / len(e), T_E), 3))
        worst = sorted(echk, key=lambda x: -abs(x[2] - x[3]))[:5]
        for d, st, p, m in worst:
            L.append("   %-10s %-3s predicted %+.3f   measured %+.3f" % (d, st, p, m))
    L += ["", "   THE MAP E: %s" % eword, ""]
    L += ["=" * W, "WHAT THIS CANNOT SAY (section 6 of the pre-registration, in short; its own words govern)", "=" * W, ""]
    L.append(P("One frame, frozen, in vacuo, at the working level. The gradients are of bare charges at fixed geometries; no design has been relaxed here. A map that PASSES agrees with direct runs on D56's designs; it is first order, and the second search's own runs measure every design it places in full."))
    L += ["", "Runs: %d of 186 ended normally. THE MAP F: %s. THE MAP E: %s." % (cnt.get("ENDED NORMALLY", 0), fword, eword)]
    txt = "\n".join(L) + "\n"
    with open(a.result, "w") as f:
        f.write(txt)
    print(txt, end="")


if __name__ == "__main__":
    main()
