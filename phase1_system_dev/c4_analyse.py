#!/usr/bin/env python3
"""
c4_analyse.py - C4 analysis (Phase 1 audit item C4; C4_CRITERIA.txt). Python 3.6, standard library; runs on hpc1 and,
from the committed files, on the PC (byte-identical report). Reads, in the C4 folder:
  c4_bands.tsv, runs.tsv, active_atoms.txt, frames/frame_N/{image0_active.xyz[.xz], md_substrate.xyz}
  per run: reopt.inp, reopt.out[.xz], reopt-minimize-ener.csv, final_active.xyz[.xz]
Writes C4_REPORT.txt, c4_results.tsv (step 2, per frame) and c4_step0.tsv in the same folder.
USAGE  python3 c4_analyse.py [C4 folder]
"""
import lzma, os, re, sys
D = sys.argv[1] if len(sys.argv) > 1 else "/home/18660916/system_development/05_qmmm/23_c4_reactant"
H = 627.5094740631
QM0, NQM = 6207, 24
EH_BOHR = 2625.49964 / 0.529177210903          # kJ/mol/A per Eh/bohr
RT300 = 8.314462618e-3 * 300.0 / 4.184         # kcal/mol, Ryde 2017's T = 300 K
# Ryde 2017, Table 2, cumulant approximation: energies needed for theta = 4, 10, 20 kJ/mol (95% confidence), by sigma
RYDE_CA = {1: (2, 2, 2), 2: (2, 2, 2), 3: (5, 2, 2), 4: (10, 3, 2), 5: (18, 4, 2), 6: (35, 6, 3), 7: (55, 11, 4),
           8: (101, 15, 5), 9: (147, 24, 7), 10: (228, 35, 10), 11: (318, 50, 13), 12: (448, 68, 18), 13: (593, 102, 24),
           14: (809, 127, 32), 15: (1040, 167, 42), 16: (1330, 214, 53), 17: (1680, 266, 66), 18: (2290, 340, 83),
           19: (2650, 407, 106), 20: (3170, 522, 127), 25: (7420, 1320, 292), 30: (15760, 2590, 627),
           35: (30100, 4580, 1120), 40: (51000, 7670, 2070)}
def text(p):
    if os.path.exists(p): return open(p, errors="replace").read()
    if os.path.exists(p + ".xz"): return lzma.open(p + ".xz", "rt", errors="replace").read()
    return None
def tsv(p):
    L = [l.rstrip("\n").split("\t") for l in open(p) if l.strip() and not l.startswith("#")]
    return L
def xyz_idx(p):
    t = text(p)
    if t is None: return None
    L = t.splitlines(); n = int(L[0].split()[0])
    return [(int(l.split()[4]), tuple(float(v) for v in l.split()[1:4])) for l in L[2:2 + n]]
def md_sub(p):
    L = open(p).read().splitlines()
    return [tuple(float(v) for v in l.split()[1:4]) for l in L[2:2 + NQM]]
def kabsch_rmsd(A, B):
    n = len(A)
    ca = [sum(p[k] for p in A) / n for k in range(3)]; cb = [sum(p[k] for p in B) / n for k in range(3)]
    a = [[p[k] - ca[k] for k in range(3)] for p in A]; b = [[p[k] - cb[k] for k in range(3)] for p in B]
    S = [[sum(a[i][r] * b[i][c] for i in range(n)) for c in range(3)] for r in range(3)]
    (xx, xy, xz), (yx, yy, yz), (zx, zy, zz) = S
    N = [[xx + yy + zz, yz - zy, zx - xz, xy - yx], [yz - zy, xx - yy - zz, xy + yx, zx + xz],
         [zx - xz, xy + yx, -xx + yy - zz, yz + zy], [xy - yx, zx + xz, yz + zy, -xx - yy + zz]]
    lam = max(abs(v) for row in N for v in row) * 4 + 1.0
    v = [1.0, 0.0, 0.0, 0.0]
    for _ in range(500):
        w = [sum((N[i][j] + (lam if i == j else 0.0)) * v[j] for j in range(4)) for i in range(4)]
        s = sum(x * x for x in w) ** 0.5; v = [x / s for x in w]
    emax = sum(v[i] * sum(N[i][j] * v[j] for j in range(4)) for i in range(4))
    ga = sum(x * x for p in a for x in p); gb = sum(x * x for p in b for x in p)
    return max(0.0, (ga + gb - 2.0 * emax) / n) ** 0.5
def disp(A, B):
    """active-region displacement, no superposition: (RMS, max, atom index of max)"""
    sq = [(sum((A[i][1][k] - B[i][1][k]) ** 2 for k in range(3)), A[i][0]) for i in range(len(A))]
    m = max(sq)
    return (sum(x[0] for x in sq) / len(sq)) ** 0.5, m[0] ** 0.5, m[1]
def echo(t):
    i = t.find("INPUT FILE"); j = t.find("****END OF INPUT****", i)
    if i < 0 or j < 0: return None
    e = [m.group(1).rstrip() for m in re.finditer(r"^\|\s*\d+> ?(.*)$", t[i:j], re.M)]
    while e and not e[-1]: e.pop()
    return e
ACT = [int(x) for x in open(os.path.join(D, "active_atoms.txt")).read().split()]
BL = tsv(os.path.join(D, "c4_bands.tsv")); BH = BL[0]
BAND = {r[0]: dict(zip(BH, r)) for r in BL[1:]}
FRAMES = [r[0] for r in BL[1:]]
RUNS = tsv(os.path.join(D, "runs.tsv"))
rep = []; say = rep.append
PLATEAU_EH, PLATEAU_N = 2.5e-4, 50      # ORCA's Normal energy criterion, 5e-6 Eh per step, over the final 50 steps
def parse_part(d):
    """one optimisation run folder. m0: problems that make it not assessable (input, settings, records)"""
    t = text(os.path.join(d, "reopt.out"))
    r = {"m0": [], "fin": False, "grad": False, "plateau": False}
    if t is None: r["m0"].append("no output"); return r
    raw = open(os.path.join(d, "reopt.inp")).read()
    inp = [l.rstrip() for l in raw.splitlines()]
    while inp and not inp[-1]: inp.pop()
    if echo(t) != inp: r["m0"].append("echoed input differs from reopt.inp")
    want = {k: float(re.search(k + r"\s+(\S+)", raw).group(1)) * EH_BOHR for k in ("TolMaxG", "TolRMSG")}
    tr = re.findall(r"RMS\(Grad\) Convergence:\s+(\d+\.\d+)", t); tm = re.findall(r"Max\(Grad\) Convergence:\s+(\d+\.\d+)", t)
    if (not tr or not tm or abs(float(tr[0]) - want["TolRMSG"]) > max(6e-4, 2e-3 * want["TolRMSG"])
            or abs(float(tm[0]) - want["TolMaxG"]) > max(6e-4, 2e-3 * want["TolMaxG"])):
        r["m0"].append("thresholds %s/%s" % (tm[0] if tm else "?", tr[0] if tr else "?"))
    cp = os.path.join(d, "reopt-minimize-ener.csv")
    rows = [l.split(";") for l in open(cp) if l.strip() and not l.lstrip().startswith("#")] if os.path.exists(cp) else []
    if not rows: r["m0"].append("no optimiser energies"); return r
    if int(rows[0][0]) != 1: r["m0"].append("energy file does not start at step 1")
    E = [float(x[5]) for x in rows]
    r["steps"] = len(E); r["E_start"] = E[0]; r["E_min"] = E[-1]
    r["slide"] = (E[-1 - PLATEAU_N] - E[-1]) if len(E) > PLATEAU_N else None
    r["grad"] = "ORCA TERMINATED NORMALLY" in t and "The minimization has converged." in t
    r["plateau"] = r["slide"] is None or r["slide"] <= PLATEAU_EH
    r["fin"] = r["grad"] and r["plateau"]
    if not r["grad"]: return r
    nt = re.findall(r"N\(Total\)\s*:\s*(\d+\.\d+)\s+electrons", t)
    if not nt or round(float(nt[-1])) != 118: r["m0"].append("electrons %s" % (nt[-1] if nt else "not found"))
    if not re.search(r"Total Charge:\s+-2\s", t) or not re.search(r"Total Multiplicity:\s+1\s", t): r["m0"].append("charge/multiplicity")
    og = re.findall(r"Max\. Grad\.\s+\|\s+(\d+\.\d+)\s+\|", t)
    r["gmax"] = float(og[-1]) / EH_BOHR if og else None
    sp = re.findall(r"FINAL SINGLE POINT ENERGY \(QM/MM\)\s+(-?\d+\.\d+)", t)
    r["Esp"] = float(sp[-1]) if sp else None
    r["geo"] = xyz_idx(os.path.join(d, "final_active.xyz"))
    if r["geo"] is None or [g[0] for g in r["geo"]] != ACT: r["m0"].append("final_active.xyz missing or not the active atoms (re-run c4_extract.py)")
    return r
def parse_run(run):
    """a run and, if it was continued (criteria M1), its continuation <run>_c: one optimisation, E_start from the first part"""
    r = parse_part(os.path.join(D, run)); r["run"] = run; r["cont"] = False
    dc = os.path.join(D, run + "_c")
    if os.path.isfile(os.path.join(dc, "reopt.inp")) and not r["m0"]:
        c = parse_part(dc)
        r["m0"] = ["continuation: " + x for x in c["m0"]]
        if not c["m0"]:
            for k in ("fin", "grad", "plateau", "slide", "E_min", "gmax", "Esp", "geo"):
                if k in c: r[k] = c[k]
            r["steps"] += c["steps"]
        r["cont"] = True
    return r
P = {}; I0 = {}; MD = {}
for f in FRAMES:
    I0[f] = xyz_idx(os.path.join(D, "frames", "frame_" + f, "image0_active.xyz"))
    MD[f] = md_sub(os.path.join(D, "frames", "frame_" + f, "md_substrate.xyz"))
    if [g[0] for g in I0[f]] != ACT: sys.exit("STOP: frame %s image0_active.xyz is not the active atoms" % f)
for r in RUNS:
    p = parse_run(r[0]); p.update(frame=r[1], step=r[2], img=int(r[3]), Ew=float(r[4])); P[r[0]] = p
say("C4 report (c4_analyse.py) - criteria from C4_CRITERIA.txt; %d runs (runs.tsv), %d frames (c4_bands.tsv)" % (len(RUNS), len(FRAMES)))
say("energies kcal/mol; E_start, E_min: the optimiser's first- and last-step energies; Ew: written band; E: final band")
# ---- machinery
bad = {k: v["m0"] for k, v in P.items() if v["m0"]}
say("\nmachinery")
say("  M0 (echoed input = reopt.inp, thresholds as in reopt.inp, optimiser energies from step 1; once the gradients are met"
    " also 118 electrons, charge -2 singlet and the final active region): %s" % ("PASS, all %d runs" % len(P) if not bad else
    "FAIL (not assessable) for " + "; ".join("%s: %s" % (k, ", ".join(v)) for k, v in sorted(bad.items()))))
ok = {k: v for k, v in P.items() if not v["m0"]}
for v in ok.values(): v["m2"] = (v["E_start"] - v["Ew"]) * H
conv = sorted(k for k, v in ok.items() if v["fin"]); nconv = sorted(k for k, v in ok.items() if not v["fin"])
contd = sorted(k for k, v in ok.items() if v["cont"])
say("  M1 finished (terminated normally, gradients converged, energy over the final %d steps falling by at most %.1e Eh):"
    " %d of %d assessable runs%s%s" % (PLATEAU_N, PLATEAU_EH, len(conv), len(ok),
    ("; continued once: " + " ".join(contd)) if contd else "",
    ("; not finished: " + "; ".join("%s (%s)" % (k, "gradients met, energy still falling %.2f kcal/mol over the final %d steps"
                                                % (ok[k]["slide"] * H, PLATEAU_N) if ok[k]["grad"] else "stopped before the gradients were met")
                                    for k in nconv)) if nconv else ""))
m2 = sorted(v["m2"] for v in ok.values())
m2med = (m2[len(m2) // 2] + m2[(len(m2) - 1) // 2]) / 2 if m2 else 0.0
m2bad = sorted(k for k, v in ok.items() if abs(v["m2"] - m2med) > 0.10)
if m2:
    say("  M2 E_start - Ew(start image): median %+.3f kcal/mol (the offset between the optimiser's energy and the band's),"
        " range %+.3f to %+.3f; every run within 0.10 of the median -> %s" % (m2med, m2[0], m2[-1],
        "PASS" if not m2bad else "FAIL for " + " ".join(m2bad)))
sps = [(v["Esp"] - v["E_min"]) * H for v in ok.values() if v.get("Esp") is not None]
if sps:
    po = [float(BAND[f]["opt_offset"]) for f in FRAMES]
    say("  informative: final single point minus the optimiser's last energy %+.3f to %+.3f kcal/mol (Phase 1's reactant"
        " optimisations: %+.3f to %+.3f)" % (min(sps), max(sps), min(po), max(po)))
good = {k for k, v in ok.items() if v["fin"] and abs(v["m2"] - m2med) <= 0.10}
def B(f, key): return float(BAND[f][key])
# ---- step 0
say("\nstep 0 - which kind of dip (F0 = E_start - E_min of the image-0 run; d = E_0 - E_klow; incomplete if F0 >= d - 1)")
say("  frame  run       d       F0      verdict                    S0 end - S2 end: dE_min  active RMS")
s0rows = []
for f in FRAMES:
    run = "S0_%s" % f if "S0_%s" % f in P else ("S2_%s" % f if int(BAND[f]["k_low"]) == 0 else None)
    if run is None: continue
    if run not in good:
        say("  %s  %-8s  run not assessable or not converged" % (f, run)); continue
    v = P[run]; F0 = (v["E_start"] - v["E_min"]) * H; d = B(f, "depth")
    verdict = "no dip" if int(BAND[f]["k_low"]) == 0 else ("incomplete optimisation" if F0 >= d - 1.0 else "separate local minimum")
    extra = ""
    if run.startswith("S0") and "S2_%s" % f in good:
        w = P["S2_%s" % f]; rm = disp(v["geo"], w["geo"])
        extra = "%+7.2f   %.3f A" % ((v["E_min"] - w["E_min"]) * H, rm[0])
    say("  %s  %-8s  %6.2f  %6.2f   %-25s  %s" % (f, run, d, F0, verdict, extra))
    s0rows.append((f, run, d, F0, verdict))
# ---- step 2
say("\nstep 2 - is k_low a reactant? A1 |dE_reopt| <= 1.0, dE_reopt = (E_min - E_start) + (Ew_klow - E_klow); A2 substrate"
    " RMSD to the MD structure <= 0.30 A; A3 active-region RMS displacement from image 0 <= 0.20 A")
say("  frame  k_low  B_img0  B_klow  dE_reopt A1  RMSD_MD A2  RMS_act A3  max_disp(atom)  steps  verdict")
res = []
for f in FRAMES:
    run = "S2_%s" % f; v = P[run]; b = BAND[f]
    B0 = B(f, "barrier_image0"); Bk = (B(f, "E_ci") - B(f, "E_klow")) * H
    if run not in good:
        if run in ok and not v["fin"] and v["cont"]:
            verdict = "fails step 2 (no minimum after the continuation) -> step 3"
        elif run in ok and not v["fin"]:
            verdict = "waiting: continuation (M1)"
        else:
            verdict = "waiting: " + ("not assessable (M0)" if run not in ok else "M2 failed")
        say("  %s  %-5s  %6.2f  %6.2f  %s" % (f, b["k_low"], B0, Bk, verdict))
        res.append(dict(f=f, k=b["k_low"], B0=B0, Bk=Bk, verdict=verdict))
        continue
    dE = (v["E_min"] - v["E_start"]) * H + (float(b["Ew_klow"]) - B(f, "E_klow")) * H
    sub = [g[1] for g in v["geo"] if QM0 <= g[0] < QM0 + NQM]
    rmd = kabsch_rmsd(sub, MD[f]); act = disp(v["geo"], I0[f])
    a1, a2, a3 = abs(dE) <= 1.0, rmd <= 0.30, act[0] <= 0.20
    verdict = "accepted" if (a1 and a2 and a3) else "fails step 2 (%s) -> step 3" % ",".join(n for n, x in (("A1", a1), ("A2", a2), ("A3", a3)) if not x)
    say("  %s  %-5s  %6.2f  %6.2f  %+6.2f   %-3s %.3f    %-3s %.3f    %-3s %.3f (%d)  %5d  %s" % (
        f, b["k_low"], B0, Bk, dE, "ok" if a1 else "NO", rmd, "ok" if a2 else "NO", act[0], "ok" if a3 else "NO",
        act[1], act[2], v["steps"], verdict))
    res.append(dict(f=f, k=b["k_low"], B0=B0, Bk=Bk, dE=dE, Br=Bk - dE, rmd=rmd, act=act[0], mx=act[1], mxa=act[2],
                    steps=v["steps"], gmax=v["gmax"], verdict=verdict))
acc = [r for r in res if r["verdict"] == "accepted"]
fail = [r["f"] for r in res if r["verdict"].startswith("fails")]
wait = [r["f"] for r in res if r["verdict"].startswith("waiting")]
say("  accepted %d of %d; to step 3: %s; waiting: %s" % (len(acc), len(res), " ".join(fail) or "none",
    "; ".join("%s (%s)" % (r["f"], r["verdict"][9:]) for r in res if r["verdict"].startswith("waiting")) or "none"))
# ---- step 4
def stats(x):
    n = len(x); m = sum(x) / n; s = (sum((v - m) ** 2 for v in x) / (n - 1)) ** 0.5 if n > 1 else float("nan")
    return n, m, s, s / n ** 0.5
final = not fail and not wait
say("\nstep 4 - result%s" % ("" if final else " (PROVISIONAL: over the accepted frames only; final after step 3 and any"
                                            " corrected runs)"))
if len(acc) >= 2:
    n, m, s, e = stats([r["Bk"] for r in acc]); _, m0, s0, e0 = stats([r["B0"] for r in acc]); _, mr, sr, er = stats([r["Br"] for r in acc])
    say("  barrier from k_low (the result): %.2f +- %.2f (sem) kcal/mol, sd %.2f kcal/mol = %.1f kJ/mol, n = %d" % (m, e, s, s * 4.184, n))
    say("  barrier from image 0, same frames: %.2f +- %.2f, sd %.2f (shown once, for comparison)" % (m0, e0, s0))
    say("  barrier from the re-optimised reactants (informative): %.2f +- %.2f, sd %.2f" % (mr, er, sr))
    _, md_, sd_, ed_ = stats([r["dE"] for r in acc])
    say("  dE_reopt over these frames: mean %+.2f +- %.2f (sem): the barrier from k_low is %s than the barrier from the"
        " re-optimised reactant by %.2f kcal/mol on average" % (md_, ed_, "lower" if md_ < 0 else "higher", abs(md_)))
    for lab, mm, ss in (("k_low", m, s), ("image 0", m0, s0)):
        sk = ss * 4.184; keys = sorted(RYDE_CA)
        lo = max(k for k in keys if k <= max(sk, 1)); hi = min((k for k in keys if k >= sk), default=keys[-1])
        n20 = (RYDE_CA[lo][2], RYDE_CA[hi][2]); n10 = (RYDE_CA[lo][1], RYDE_CA[hi][1])
        say("  cumulant (Ryde 2017, T = 300 K), from %s: mean - sd^2/(2RT) = %.2f kcal/mol; sd %.1f kJ/mol needs %d-%d energies for"
            " +-20 kJ/mol and %d-%d for +-10 (Table 2, cumulant); n = %d %s +-20, %s +-10" % (
                lab, mm - ss ** 2 / (2 * RT300), sk, n20[0], n20[1], n10[0], n10[1], n,
                "meets" if n >= n20[1] else ("may meet" if n >= n20[0] else "does not meet"),
                "meets" if n >= n10[1] else ("may meet" if n >= n10[0] else "does not meet")))
else:
    say("  fewer than two accepted frames: no statistics")
cols = ["frame", "k_low", "barrier_image0", "barrier_klow", "dE_reopt", "barrier_reopt", "rmsd_md", "act_rms", "act_max",
        "act_max_atom", "steps", "final_maxgrad_Eh_bohr", "verdict"]
def fmt(r, c):
    x = r.get({"frame": "f", "k_low": "k", "barrier_image0": "B0", "barrier_klow": "Bk", "dE_reopt": "dE", "barrier_reopt": "Br",
               "rmsd_md": "rmd", "act_rms": "act", "act_max": "mx", "act_max_atom": "mxa", "steps": "steps",
               "final_maxgrad_Eh_bohr": "gmax", "verdict": "verdict"}[c])
    if x is None: return ""
    if isinstance(x, float): return "%.2e" % x if c == "final_maxgrad_Eh_bohr" else "%.3f" % x
    return str(x)
open(os.path.join(D, "c4_results.tsv"), "w").write("\t".join(cols) + "\n" + "".join("\t".join(fmt(r, c) for c in cols) + "\n" for r in res))
open(os.path.join(D, "c4_step0.tsv"), "w").write("frame\trun\tdepth\tF0\tverdict\n" + "".join("%s\t%s\t%.3f\t%.3f\t%s\n" % r for r in s0rows))
open(os.path.join(D, "C4_REPORT.txt"), "w").write("\n".join(rep) + "\n")
print("\n".join(rep))
