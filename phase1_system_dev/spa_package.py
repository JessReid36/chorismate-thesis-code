#!/usr/bin/env python3
"""
spa_package.py - the report of the 30-frame stationary-point package under SP Amendment 2 (results
05_qmmm/24_sp_barriers/SP_AMENDMENT2.txt, with SP_PLAN.txt, SP_AMENDMENT1.txt and SP_ERRATUM1.txt): every frame's state
and values, the attrition rule and the selectivity check (A2.2), the ensemble statistics and the Phase 1 barrier (A2.3),
the surroundings check with its replacement rule (A2.4), the cycle counts (A2.6) and frame 42436's smoke-test-1 reactant
(A2.1). Committed with the amendment, before any package input exists (SP_PLAN.txt, FILES AND CODE). Runs on hpc1 or
on a fetched copy. Python 3.6, standard library; deterministic; reads only, except the per-frame table it writes when
given a TSV path.

  spa_package.py SD BANDS [TSV]
      SD     the folder holding the run folders: f<FRAME>_r8_a1 (S = md, the 30 frames) and f<FRAME>_r8_ci_a1 (S = ci:
             frames 20000, 44388, 59999 and any replacement under A2.4)
      BANDS  05_qmmm/23_c4_reactant/c4_bands.tsv (the 30 frames; barrier_image0 and barrier_klow per frame)
      TSV    optional: write the per-frame table to this new file (an existing file is not overwritten)
  spa_package.py selftest
      the statistics against reference values computed independently (spa_package_stats.py)

Each run folder's own report (spa_report.py) is run first, silently; every record that differs from its ORCA output is
printed as PROBLEM at the end, and the package's numbers are then not to be used until it is resolved. The output is
PROVISIONAL while anything is undecided (a job not finished, a machinery stop not yet re-run, a replacement not yet
run); it is run again after every re-run, and the result is the first output with nothing undecided (A2.9).
Energies in kcal/mol (Eh x 627.5094740631) unless stated.
"""
import contextlib, io, math, os, re, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))   # spa_report.py and spa_package_stats.py beside this file
import spa_report as SR
import spa_package_stats as ST

H = SR.H
RADIUS = 8                          # smoke test 2 met the region rule at 8 A (SP_AMENDMENT1.txt A1.5; f42436_r8_a1)
N_FRAMES = 30
N_MIN = 20                          # SP_PLAN.txt PACKAGE; SP_AMENDMENT2.txt A2.2
SURR_FRAMES = ("20000", "44388", "59999")   # SP_PLAN.txt PACKAGE: the 1st, 15th and 30th by frame number
SURR_DIR = {"20000": 1, "44388": 1, "59999": -1}   # A2.4: a replacement is the nearest frame later (earlier for 59999)
SURR_TOL = 2.0                      # kcal/mol: |B_sp(ci) - B_sp(md)| > 2 exceeds (SP_PLAN.txt PACKAGE; A2.4)
SURR_NFAIL = 2                      # "in at least two of the three"
CHAIN = ("ts", "ts_freq", "irc", "r_irc", "r_irc_freq")    # SP_PLAN.txt QUANTITIES: B_sp is complete when all passed
PAIRS = (("p_irc", "p_irc_freq"), ("r0", "r0_freq"))     # outside the chain: re-run at the job level (A2.2)
ALT = SR.ALT
DH_EXP = (12.7, 0.4)                # BsCM activation enthalpy (Kast, Asif-Ullah & Hilvert 1996, as Ranaghan et al. 2004 and
DG_EXP = 15.4                       #   Claeyssens et al. 2011 quote it), kcal/mol; activation free energy at 25 C
CLAEYSSENS = (11.3, 1.8)            # Claeyssens et al. 2011, Table 2, BsCM, B3LYP/6-31G(d)/CHARMM27: the barrier of the
                                    #   average of 16 potential-energy profiles, with the sd of the 16 pathways
T_RYDE = 300.0                      # K: Ryde 2017's Table 2 is computed at 300 K (A2.3)
STOP1 = "_stop1"

def quiet_report(rd):
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        R = SR.report(rd)
    return R
def job_of(step): return [j for j, steps in SR.JOBS if step in steps][0]
def items(r): return ", ".join("%s=%s" % (k, v) for k, v in r.items() if k not in ("status", "reason", "history")
                               and not k.endswith("sha256"))
def left_queue(rd):
    """jobs of this run folder that qstat no longer lists although their output has no last line (written by the bundle
    that collects the outputs; A2.2)"""
    p = os.path.join(rd, "LEFT_QUEUE.txt")
    return set(l.split()[0] for l in open(p) if l.strip() and not l.startswith("#")) if os.path.isfile(p) else set()
def step_state(rd, R, s):
    """("PASS", ""); ("FAIL", why) - a check failed, or PBS killed the job (walltime); ("INPUT", why) - ORCA reported an
    input error (SP_AMENDMENT1.txt A1.5: corrected in a disclosed amendment); ("STOP", why) - another machinery stop
    (A1.5); ("WAIT", why) - not finished"""
    r = R[s]
    if r.get("status") == "PASS": return "PASS", ""
    job = job_of(s)
    st, use, kb, fails = SR.job_state(rd, job)
    if r.get("status") == "FAIL":
        if SR.disk_failure(rd, s, job): return "STOP", "%s failed and its output reports a full disk or quota" % s
        why = "%s: %s" % (s, r.get("reason", "no reason recorded"))
        if items(r): why += " [values reached: %s]" % items(r)
        if r.get("history"): why += " [cycles: %s]" % r["history"]
        return "FAIL", why
    if st.startswith("killed"): return "FAIL", "%s: %s %s" % (s, job, st)
    if st == "ended":
        w = "%s: %s ended without a record for this step (%s)" % (s, job, "; ".join(fails) or "no FAIL line")
        return ("INPUT" if any("input error" in l.lower() for l in fails) else "STOP"), w
    if st == "running, or ended without its last line" and job in left_queue(rd):
        return "STOP", "%s: %s left the queue without its last line and without a record (LEFT_QUEUE.txt; A2.2)" % (s, job)
    return "WAIT", "%s: %s %s (if qstat no longer lists it, the collection bundle records it in LEFT_QUEUE.txt and it is a"\
                   " machinery stop, A2.2)" % (s, job, st)
def classify(rd, R, s, level):
    """one step's state for the package: pass, failed, stopped (a machinery stop, not decided), not finished.
    level "folder": a chain step, re-run with its whole run folder (A1.5); level "step": a step outside the chain,
    re-run with the job that depends on it (A2.2)"""
    st, w = step_state(rd, R, s)
    if st == "PASS": return "pass", ""
    if st == "FAIL": return "failed", w
    if st == "INPUT":
        return "stopped", ("input error (SP_AMENDMENT1.txt A1.5: the input is corrected in a disclosed amendment, then run"
                           " again): " + w)
    if st == "STOP":
        if level == "folder":
            if os.path.isdir(rd + STOP1): return "failed", "a second machinery stop, after %s%s: %s" % (os.path.basename(rd), STOP1, w)
            return "stopped", ("machinery stop, not decided (A1.5: rename the run folder to %s%s, generate it again, resubmit"
                               " once): %s" % (os.path.basename(rd), STOP1, w))
        if os.path.isdir(os.path.join(rd, s + STOP1)):
            return "failed", "a second machinery stop, after %s%s: %s" % (s, STOP1, w)
        return "stopped", ("machinery stop, not decided (A2.2: rename %s and its job's output with %s, re-run that job and"
                           " the job depending on it once): %s" % (s, STOP1, w))
    return "not finished", w
def cycles_run(rd, s):
    d = os.path.join(rd, s)
    return len([x for x in os.listdir(d) if re.match(r"^c\d\d_env$", x)]) if os.path.isdir(d) else 0
def fnum(x, f="%.2f"): return "-" if x is None else f % x

def frame(rd):
    """everything the package reports for one run folder"""
    out = {"folder": os.path.basename(rd), "state": "absent", "why": "run folder not present", "failed_step": None,
           "B_sp": None, "B_eff": None, "r0_used": None, "r0_note": "", "dR0": None, "ts_minus_r0": None,
           "dE_rxn": None, "imag": None, "lower_bound": None, "pairs": {}, "notes": []}
    for s in ALT: out["cycles_" + s] = None; out["pass_" + s] = False; out["lastdE_" + s] = None
    if not os.path.isdir(rd): return out
    R = quiet_report(rd)
    for s in ALT:
        out["cycles_" + s] = int(R[s]["cycles"]) if R[s].get("cycles") else (cycles_run(rd, s) or None)
        out["pass_" + s] = R[s].get("status") == "PASS"
        if R[s].get("last_core_change_kcal"): out["lastdE_" + s] = float(R[s]["last_core_change_kcal"])
    E = dict((s, float(R[s]["energy_Eh"])) for s in ALT if R[s].get("status") == "PASS" and R[s].get("energy_Eh"))
    state, why, at = "complete", "", None
    for s in CHAIN:
        st, w = classify(rd, R, s, "folder")
        if st == "pass": continue
        state, why, at = st, w, s
        break
    out.update(state=state, why=why, failed_step=at)
    # outside the chain: P_IRC (dE_rxn) and R0 (dR0, B_eff), each with its frequencies
    for a, b in PAIRS:
        sa, wa = classify(rd, R, a, "step")
        if sa == "pass": sb, wb = classify(rd, R, b, "step")
        out["pairs"][a] = (sa, wa) if sa != "pass" else (sb, wb)
        if state == "failed" and out["pairs"][a][0] == "not finished" and "not started" in out["pairs"][a][1]:
            out["pairs"][a] = ("not run", out["pairs"][a][1])
        if out["pairs"][a][0] != "pass" and state == "complete": out["notes"].append(out["pairs"][a][1])
    if R["ts_freq"].get("status") == "PASS": out["imag"] = R["ts_freq"].get("imaginary_counted")
    if "p_irc" in E and "r_irc" in E: out["dE_rxn"] = (E["p_irc"] - E["r_irc"]) * H
    if "r0" in E and "r_irc" in E: out["dR0"] = (E["r0"] - E["r_irc"]) * H
    if "r0" in E and "ts" in E: out["ts_minus_r0"] = (E["ts"] - E["r0"]) * H
    if state == "complete":
        out["B_sp"] = (E["ts"] - E["r_irc"]) * H
        s0, w0 = out["pairs"]["r0"]
        if s0 == "pass":
            out["r0_used"] = True
            out["B_eff"] = (E["ts"] - min(E["r_irc"], E["r0"])) * H
            out["lower_bound"] = E["r0"] < E["r_irc"]
        elif s0 == "failed":
            out["r0_used"] = False; out["B_eff"] = out["B_sp"]; out["lower_bound"] = False
            out["r0_note"] = "R0 failed (no valid minimum): B_eff from R_IRC alone (" + w0 + ")"
        else:
            out["r0_note"] = "R0 not decided yet (" + w0 + ")"
    return out

def read_bands(path):
    rows = [l.rstrip("\n").split("\t") for l in open(path) if l.strip() and not l.startswith("#")]
    head, body = rows[0], rows[1:]
    for c in ("frame", "barrier_image0", "barrier_klow"):
        if c not in head: sys.exit("STOP: %s has no column %s" % (path, c))
    B = dict((r[head.index("frame")], dict(zip(head, r))) for r in body)
    if len(B) != N_FRAMES or len(body) != N_FRAMES: sys.exit("STOP: %s holds %d frames, not %d" % (path, len(B), N_FRAMES))
    for f in SURR_FRAMES + ("42436",):
        if f not in B: sys.exit("STOP: frame %s is not in %s" % (f, path))
    return B

def normal_plot(vals):
    """Ryde 2017 (SI Figure S1) shows histograms and normal plots: here as text"""
    x = sorted(vals); n = len(x)
    if x[-1] - x[0] < 1e-12:
        print("    normal plot and histogram: all values equal"); return
    z = [ST.norm_ppf((i - 0.375) / (n + 0.25)) for i in range(1, n + 1)]
    mx, mz = sum(x) / n, sum(z) / n
    r = sum((a - mx) * (b - mz) for a, b in zip(x, z)) / math.sqrt(sum((a - mx) ** 2 for a in x) * sum((b - mz) ** 2 for b in z))
    print("    normal plot (value against its normal score, Blom positions); correlation r = %.4f" % r)
    cells = ["%7.2f %6.2f" % (a, b) for a, b in zip(x, z)]
    for i in range(0, n, 5): print("      " + "   ".join(cells[i:i + 5]))
    lo = int(math.floor(x[0])); hi = int(math.floor(x[-1]))
    print("    histogram (1 kcal/mol bins):")
    for b in range(lo, hi + 1):
        k = sum(1 for a in x if b <= a < b + 1)
        print("      %5d-%-5d %3d %s" % (b, b + 1, k, "#" * k))
def ryde_rows(sig_kj):
    """Ryde 2017 Table 2 at its own 300 K: sigma scaled by 300/298.15 (the table depends on sigma/RT); the row at or
    above (used: the tabulated numbers rise with sigma) and the row below (shown)"""
    s300 = sig_kj * T_RYDE / ST.T_REF
    above = [r for r in ST.RYDE_T2 if r[0] >= s300]; below = [r for r in ST.RYDE_T2 if r[0] < s300]
    return s300, (above[0] if above else None), (below[-1] if below else None)
def tab(v): return "more than tabulated" if v is None else str(v)

def ensemble(name, vals, provisional):
    print("  %s%s" % (name, "  (PROVISIONAL)" if provisional else ""))
    d = ST.describe(vals)
    if d is None or d["n"] < 2:
        print("    n = %d: too few complete frames for statistics" % (0 if d is None else d["n"])); return d
    n = d["n"]
    print("    n %d   mean %.2f   sd %.2f   sem %.2f   95%% CI of the mean %.2f to %.2f (t = %.3f, %d df)" % (
        n, d["mean"], d["sd"], d["sem"], d["ci95"][0], d["ci95"][1], d["t95"], n - 1))
    print("    median %.2f   range %.2f to %.2f" % (d["median"], d["min"], d["max"]))
    print("    against experiment and the comparator (A2.3): activation enthalpy %.1f +- %.1f; Claeyssens et al. 2011"
          " %.1f +- %.1f (averaged profile; sd of 16 pathways); this mean minus the activation enthalpy %+.2f" % (
              DH_EXP[0], DH_EXP[1], CLAEYSSENS[0], CLAEYSSENS[1], d["mean"] - DH_EXP[0]))
    sig_kj = d["sd"] * ST.KJ
    ea, ca = ST.exp_average(vals), ST.cumulant(vals)
    print("    secondary, estimates of the free-energy barrier at %.2f K, against the activation free energy %.1f:" % (
        ST.T_REF, DG_EXP))
    print("      exponential average %.2f (%+.2f)   cumulant mean - sd^2/2RT %.2f (%+.2f)" % (ea, ea - DG_EXP, ca, ca - DG_EXP))
    s300, up, lo = ryde_rows(sig_kj)
    if up is None:
        print("      sd = %.2f kJ/mol (%.2f at Ryde's 300 K): beyond Ryde 2017's Table 2 (last row sigma = 40); this ensemble"
              " has %d" % (sig_kj, s300, n))
    else:
        print("      sd = %.2f kJ/mol (%.2f at Ryde's 300 K); Ryde 2017 Table 2 (95%% confidence), row sigma = %d: exponential"
              " average %s energies within 10 kJ/mol, %s within 20; cumulant %s within 10, %s within 20; this ensemble has"
              " %d%s" % (sig_kj, s300, up[0], tab(up[2]), tab(up[3]), tab(up[5]), tab(up[6]), n,
                         "" if lo is None else " (row below, sigma = %d: %s, %s; %s, %s)" % (
                             lo[0], tab(lo[2]), tab(lo[3]), tab(lo[5]), tab(lo[6]))))
    if n >= 3:
        W, p = ST.shapiro_wilk(vals); g1, g2 = ST.skew_kurt(vals)
        print("      normality: Shapiro-Wilk W %s, p %s%s; skewness %s, excess kurtosis %s" % (
            fnum(W, "%.4f"), fnum(p, "%.4f"),
            "" if p is None else (" (below 0.05: evidence against a normal distribution, so against the cumulant"
                                  " estimate)" if p < 0.05 else " (no evidence against a normal distribution)"),
            fnum(g1), fnum(g2)))
        normal_plot(vals)
    return d

def surroundings(frames, md, ci):
    """A2.4: per check frame its difference, or its replacement's, or why there is none yet; the verdict"""
    used = set(); rows = []; pend = []; nex = nin = npend = ndead = 0
    for F in SURR_FRAMES:
        a, b = md[F], ci[F]
        res = None
        if a["state"] == "complete" and b["state"] == "complete":
            res = ("result", F, a, b)
        elif a["state"] == "failed" or b["state"] == "failed":
            lost = "; ".join("%s failed (%s)" % (x["folder"], x["why"]) for x in (a, b) if x["state"] == "failed")
            i = frames.index(F); G = None; wait = None
            seq = frames[i + 1:] if SURR_DIR[F] > 0 else frames[:i][::-1]
            for g in seq:
                if g in SURR_FRAMES or g in used: continue
                if md[g]["state"] == "failed": continue
                if md[g]["state"] == "complete": G = g
                else: wait = g
                break
            if G is None and wait is None: res = ("dead", F, "%s; no replacement frame is left" % lost)
            elif G is None: res = ("pending", F, "%s; replacement: waiting for frame %s's md run (%s)" % (lost, wait, md[wait]["state"]))
            else:
                used.add(G); bg = ci[G]
                if bg["state"] == "complete": res = ("result", G, md[G], bg, lost)
                elif bg["state"] == "failed": res = ("dead", F, "%s; its replacement %s failed too (%s)" % (lost, bg["folder"], bg["why"]))
                elif bg["state"] == "absent":
                    res = ("pending", F, "%s; replacement frame %s: generate and run f%s_r%d_ci_a1 (A2.4)" % (lost, G, G, RADIUS))
                else: res = ("pending", F, "%s; replacement %s %s (%s)" % (lost, bg["folder"], bg["state"], bg["why"]))
        else:
            res = ("pending", F, "; ".join("%s %s (%s)" % (x["folder"], x["state"], x["why"]) for x in (a, b) if x["state"] != "complete"))
        if res[0] == "result":
            G, x, y = res[1], res[2], res[3]
            d = y["B_sp"] - x["B_sp"]; big = abs(d) > SURR_TOL
            nex += big; nin += not big
            de = "" if x["B_eff"] is None or y["B_eff"] is None else "; B_eff(ci) - B_eff(md) %+.3f (reported only)" % (y["B_eff"] - x["B_eff"])
            rep = "" if G == F else " (replacing %s: %s)" % (F, res[4])
            rows.append("  %s%s  md %.3f  ci %.3f  |difference| %.3f (%+.3f)  %s%s" % (
                G, rep, x["B_sp"], y["B_sp"], abs(d), d, "exceeds" if big else "within", de))
        elif res[0] == "dead":
            ndead += 1; rows.append("  %s  no difference, finally: %s" % (F, res[2]))
        else:
            npend += 1; pend.append(F); rows.append("  %s  no difference yet: %s" % (F, res[2]))
    if nex >= SURR_NFAIL:
        sv = ("FAILS: the barrier depends on which state's surroundings are frozen (Ryde 2016, item 7). The ensemble is"
              " reported but not used as the Phase 1 barrier; A2.4's remedies follow in a further amendment")
    elif nin >= len(SURR_FRAMES) - SURR_NFAIL + 1:
        sv = "passes"
    elif nex + npend < SURR_NFAIL and nin + npend < len(SURR_FRAMES) - SURR_NFAIL + 1:
        sv = ("CANNOT BE DECIDED: neither outcome can be reached (A2.4: the ensemble is reported but not used as the Phase 1"
              " barrier until a further amendment, written before any further result, decides the check)")
    else:
        sv = "not decided yet"
    return rows, sv, pend

def main(sd, bands_path, tsv=None):
    if tsv and os.path.exists(tsv): sys.exit("STOP: %s exists; give a new file name" % tsv)
    B = read_bands(bands_path)
    frames = sorted(B, key=int)
    md = dict((f, frame(os.path.join(sd, "f%s_r%d_a1" % (f, RADIUS)))) for f in frames)
    ci = dict((f, frame(os.path.join(sd, "f%s_r%d_ci_a1" % (f, RADIUS)))) for f in frames)
    for f in frames: md[f]["band_image0"], md[f]["band_klow"] = float(B[f]["barrier_image0"]), float(B[f]["barrier_klow"])
    print("=" * 120)
    print("SP PACKAGE (SP_AMENDMENT2.txt): %d frames at %d A, S = md; surroundings check frames %s with S = ci" % (
        len(frames), RADIUS, ", ".join(SURR_FRAMES)))
    print("per frame: state; B_sp = E(TS) - E(R_IRC); B_eff = E(TS) - min(E(R_IRC), E(R0)) (R0 only when it passed and its"
          " frequencies show a minimum); dR0 = E(R0) - E(R_IRC); dE_rxn; TS imaginary mode (cm-1); alternation cycles of"
          " TS / R_IRC / P_IRC / R0; the band's barrier_image0 and barrier_klow (c4_bands.tsv, descriptive)")
    print("  %-6s %-12s %7s %7s %7s %7s %9s %11s %7s %7s" % ("frame", "state", "B_sp", "B_eff", "dR0", "dE_rxn", "imag",
                                                              "cycles", "band_i0", "band_kl"))
    for f in frames:
        m = md[f]
        cyc = "/".join("-" if m["cycles_" + s] is None else str(m["cycles_" + s]) for s in ALT)
        print("  %-6s %-12s %7s %7s %7s %7s %9s %11s %7.2f %7.2f%s" % (
            f, m["state"], fnum(m["B_sp"]), fnum(m["B_eff"]), fnum(m["dR0"]), fnum(m["dE_rxn"]), m["imag"] or "-", cyc,
            m["band_image0"], m["band_klow"], "  R0 below R_IRC: B_sp is a lower bound" if m["lower_bound"] else ""))
    for f in frames:
        m = md[f]
        if m["state"] != "complete": print("  %s %s: %s" % (f, m["state"], m["why"]))
        for n_ in m["notes"]: print("  %s (outside the B_sp chain): %s" % (f, n_))
        if m["r0_note"]: print("  %s: %s" % (f, m["r0_note"]))
    comp = [f for f in frames if md[f]["state"] == "complete"]
    failed = [f for f in frames if md[f]["state"] == "failed"]
    open_ = [f for f in frames if md[f]["state"] not in ("complete", "failed")]
    side = ["%s %s" % (f, a) for f in comp for a, b in PAIRS if md[f]["pairs"][a][0] in ("stopped", "not finished")]
    print("=" * 120)
    print("ATTRITION (A2.2): %d of %d complete, %d failed, %d not decided yet (%s)" % (
        len(comp), len(frames), len(failed), len(open_), ", ".join("%s %s" % (f, md[f]["state"]) for f in open_) or "none"))
    if open_ and len(comp) + len(open_) < N_MIN:
        attr = "fewer than %d frames can now be complete" % N_MIN
    elif open_:
        attr = "not decided yet"
    else:
        attr = "at least %d complete" % N_MIN if len(comp) >= N_MIN else "fewer than %d complete" % N_MIN
    print("  rule: at least %d complete frames; else the ensemble is reported but not used as the Phase 1 barrier -> %s" % (N_MIN, attr))
    pv = "  (PROVISIONAL)" if open_ else ""
    if len(failed) >= 2 and len(comp) >= 2:
        a = [md[f]["band_klow"] for f in failed]; b = [md[f]["band_klow"] for f in comp]
        da, db = ST.describe(a), ST.describe(b)
        print("  selectivity check (A2.2): band barrier_klow of the %d failed against the %d complete frames%s" % (len(a), len(b), pv))
        print("    failed: mean %.2f, median %.2f; complete: mean %.2f, median %.2f" % (da["mean"], da["median"], db["mean"], db["median"]))
        if da["sd"] == 0 and db["sd"] == 0:
            print("    -> no spread in either group: no interval")
        else:
            w = ST.welch(a, b)
            print("    difference of means (failed - complete) %.2f, 95%% CI %.2f to %.2f (Welch, %.1f df)" % (
                w["diff"], w["ci95"][0], w["ci95"][1], w["df"]))
            if w["ci95"][0] > 0:
                print("    -> the failed frames' band barriers are higher: the complete frames' mean may be biased low (stated with the results)")
            elif w["ci95"][1] < 0:
                print("    -> the failed frames' band barriers are lower: the complete frames' mean may be biased high (stated with the results)")
            else:
                print("    -> no evidence that failure depends on the band barrier (a test of low power with few failed frames)")
    elif len(failed) == 1 and comp:
        f = failed[0]; v = md[f]["band_klow"]
        print("  selectivity check (A2.2): one failed frame, %s, band barrier_klow %.2f: above %d of the %d complete frames'%s" % (
            f, v, sum(1 for g in comp if md[g]["band_klow"] < v), len(comp), pv))
    else:
        print("  selectivity check: not applicable (%d failed, %d complete)" % (len(failed), len(comp)))
    print("=" * 120)
    print("ENSEMBLE (A2.3; complete frames, S = md)")
    prov = bool(open_)
    dsp = ensemble("B_sp = E(TS) - E(R_IRC), primary", [md[f]["B_sp"] for f in comp], prov)
    pend_eff = [f for f in comp if md[f]["B_eff"] is None]
    ensemble("B_eff = E(TS) - min(E(R_IRC), E(R0))%s" % (
        "" if not pend_eff else ", without frames whose R0 is not decided yet: " + ", ".join(pend_eff)),
        [md[f]["B_eff"] for f in comp if md[f]["B_eff"] is not None], prov or bool(pend_eff))
    lb = [f for f in comp if md[f]["lower_bound"]]
    print("  R0 below R_IRC (B_sp a lower bound) in %d of %d complete frames%s" % (
        len(lb), len(comp), (": " + ", ".join("%s (%.2f)" % (f, -md[f]["dR0"]) for f in lb)) if lb else ""))
    nr0 = [f for f in comp if md[f]["r0_used"] is False]
    if nr0: print("  R0 failed, no valid minimum (B_eff = B_sp) in %d frame(s): %s" % (len(nr0), ", ".join(nr0)))
    print("=" * 120)
    print("CYCLES (A2.6): alternation cycles per stationary point (limit 20; Amendment 1's limit was 8)")
    for s in ALT:
        c = [md[f]["cycles_" + s] for f in comp if md[f]["pass_" + s] and md[f]["cycles_" + s] is not None]
        if c: print("  %-6s passed, in complete frames: %d to %d cycles, median %s" % (s, min(c), max(c), fnum(ST.describe(c)["median"], "%g")))
    over = ["%s %s %d (%s)" % (f, s, md[f]["cycles_" + s], "passed" if md[f]["pass_" + s] else "not passed")
            for f in frames for s in ALT if (md[f]["cycles_" + s] or 0) > 8]
    print("  more than 8 cycles: %s" % (", ".join(over) or "none"))
    for s in ("ts", "r_irc"):
        c = [md[f]["lastdE_" + s] for f in comp if md[f]["lastdE_" + s] is not None]
        if c: print("  %-6s core energy change in the last cycle: %.3f to %.3f (A1.3 b: at most 0.10)" % (s, min(c), max(c)))
    print("=" * 120)
    print("SURROUNDINGS CHECK (SP_PLAN.txt PACKAGE; A2.4): |B_sp(ci) - B_sp(md)| > %.1f exceeds; fails if at least %d of %d"
          " exceed; a check frame whose md or ci run failed is replaced once (A2.4)" % (SURR_TOL, SURR_NFAIL, len(SURR_FRAMES)))
    rows, sv, pend = surroundings(frames, md, ci)
    for l in rows: print(l)
    print("  -> %s" % sv)
    print("=" * 120)
    undecided = [("frame " + f + " " + md[f]["state"]) for f in open_] + side + \
                ["surroundings check place %s" % F for F in pend] + (["records differ from outputs (PROBLEM below):"
                                                                    " resolve and run again"] if SR.PROBLEMS else [])
    stop = []
    if attr != "at least %d complete" % N_MIN and attr != "not decided yet": stop.append("attrition: " + attr)
    if sv.startswith("FAILS") or sv.startswith("CANNOT"): stop.append("surroundings check: " + sv.split(":")[0].lower())
    if stop:
        print("THE PHASE 1 BARRIER (A2.3): not this ensemble (%s)%s" % ("; ".join(stop), "; PROVISIONAL" if undecided else ""))
    elif undecided:
        print("THE PHASE 1 BARRIER (A2.3): not decided yet (PROVISIONAL; undecided: %s)" % ", ".join(undecided))
    else:
        print("THE PHASE 1 BARRIER (A2.3; audit item C3): the conditions are met (at least %d complete frames, the"
              " surroundings check passed, records agree with outputs, nothing undecided): B_sp ensemble mean %.2f kcal/mol,"
              " 95%% CI %.2f to %.2f (n = %d), with B_eff beside it" % (N_MIN, dsp["mean"], dsp["ci95"][0], dsp["ci95"][1], dsp["n"]))
    print("RESULT: %s" % ("PROVISIONAL (run again once nothing is undecided)" if undecided else "FINAL (nothing undecided)"))
    p1 = os.path.join(sd, "f42436_r8", "r0", "STAGE_RESULT.txt")
    m = md["42436"]
    r = SR.result(os.path.join(sd, "f42436_r8"), "r0")
    if not os.path.isfile(p1):
        print("FRAME 42436, smoke test 1's R0 (A2.1): not found (%s)" % p1)
    elif r.get("status") != "PASS" or m["state"] != "complete":
        print("FRAME 42436, smoke test 1's R0 (A2.1): no value (its record: %s; this frame's package run: %s)" % (
            r.get("status", "-"), m["state"]))
    else:
        rf = SR.result(os.path.join(sd, "f42436_r8"), "r0_freq")
        Rm = SR.result(os.path.join(sd, "f42436_r8_a1"), "r_irc"); Tm = SR.result(os.path.join(sd, "f42436_r8_a1"), "ts")
        e1 = float(r["energy_Eh"])
        print("FRAME 42436, smoke test 1's R0 (f42436_r8/r0; SP_PLAN's whole-region optimisation, exploratory; A2.1:"
              " reported as found, not part of the ensemble): E %.6f Eh, %.2f below this frame's R_IRC; E(TS) - E(R0,"
              " smoke test 1) %.2f against B_sp %.2f; its frequencies: %s, counted imaginary modes %s" % (
                  e1, (float(Rm["energy_Eh"]) - e1) * H, (float(Tm["energy_Eh"]) - e1) * H, m["B_sp"],
                  rf.get("status", "not run"), rf.get("imaginary_counted", "-")))
    print("=" * 120)
    print("PROBLEM: " + "\nPROBLEM: ".join(SR.PROBLEMS) if SR.PROBLEMS else "records against outputs (every run folder): no problem")
    if tsv:
        cols = ["frame", "surroundings", "folder", "state", "failed_step", "why", "B_sp", "B_eff", "r0_used", "lower_bound",
                "dR0", "ts_minus_r0", "dE_rxn", "ts_imag_cm1"] + ["cycles_" + s for s in ALT] + \
               ["last_core_change_" + s for s in ALT] + ["p_irc_state", "r0_state", "band_barrier_image0", "band_barrier_klow"]
        with open(tsv, "x") as fh:
            fh.write("# spa_package.py %s %s (SP_AMENDMENT2.txt); kcal/mol; %s\n" % (
                sd, bands_path, "PROVISIONAL" if undecided else "FINAL"))
            fh.write("\t".join(cols) + "\n")
            rowsets = [(f, "md", md[f]) for f in frames] + [(f, "ci", ci[f]) for f in frames
                                                            if f in SURR_FRAMES or ci[f]["state"] != "absent"]
            for f, S, x in rowsets:
                pa = x["pairs"].get("p_irc", ("", "")); pr = x["pairs"].get("r0", ("", ""))
                v = [f, S, x["folder"], x["state"], x.get("failed_step") or "", x["why"].replace("\t", " "),
                     fnum(x["B_sp"], "%.4f"), fnum(x["B_eff"], "%.4f"), "" if x["r0_used"] is None else str(x["r0_used"]),
                     "" if x["lower_bound"] is None else str(x["lower_bound"]), fnum(x["dR0"], "%.4f"),
                     fnum(x["ts_minus_r0"], "%.4f"), fnum(x["dE_rxn"], "%.4f"), x["imag"] or ""] + \
                    ["" if x["cycles_" + s] is None else str(x["cycles_" + s]) for s in ALT] + \
                    [fnum(x["lastdE_" + s], "%.4f") for s in ALT] + \
                    [(pa[0] + (": " + pa[1] if pa[1] else "")).replace("\t", " "),
                     (pr[0] + (": " + pr[1] if pr[1] else "")).replace("\t", " "),
                     B[f]["barrier_image0"], B[f]["barrier_klow"]]
                fh.write("\t".join(v) + "\n")
        print("per-frame table written: %s" % tsv)

if __name__ == "__main__":
    a = sys.argv[1:]
    if a == ["selftest"]: sys.exit(1 if ST.selftest() else 0)
    if len(a) in (2, 3): main(*a)
    else: sys.exit(__doc__)
