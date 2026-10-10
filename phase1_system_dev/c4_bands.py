#!/usr/bin/env python3
"""
c4_bands.py - C4 band table (PHASE1_AUDIT_CHECKLIST.md item C4; C4_CRITERIA.txt). Runs on the PC from committed files
only; Python 3.6, standard library. For each of the 30 reported frames it reads the final band:
  05_qmmm/22_c5_nebci/frame_N/   for the eight C5 frames (NEB-CI restarts)
  05_qmmm/19_ensemble/frame_N/   for the other 22 (NEB-CI)
Two sets of image energies exist for every band, and both are used:
  E_k   the final band, the one every tabulated barrier comes from: the last PATH SUMMARY of neb.out (Eh, 5 decimals)
  Ew_k  the band as last written to disk (neb_MEP.QMRegion_trj.xyz and, on hpc1, neb_MEP.allxyz; Eh, 12 decimals):
        ORCA writes these files at the last iteration in neb.NEB.log, which is one before the converged band of the
        PATH SUMMARY, so a re-optimisation started from a written image starts at Ew_k, not E_k; the end points are
        fixed and agree
k_low is the lowest image of the final band at or before the climbing image (image 0 included). Checks (STOP):
  - the frame list is A_v2's 30 frames, which are the full_NAC frames numbered 20000 and above in the harvest
  - neb.out.xz decompresses to the sha256 recorded for the original on hpc1; the band converged; one climbing image
  - Ew_k equal the last iteration of neb.NEB.log (6e-9 Eh), and neb.interp (whose last block is the converged band)
    holds one iteration more than neb.NEB.log; the end points agree between the two bands (6e-6 Eh)
  - image 0's written energy equals the Phase 1 reactant optimisation's final single point (1e-6 Eh)
  - the harvested barrier (ensemble_barriers.tsv; C5 frames: c5_results.tsv) is the PATH SUMMARY's dE of the climbing
    image, and the barrier from image 0 computed here agrees with it within 0.011 kcal/mol (5- and 2-decimal rounding)
Informative column opt_offset: Phase 1's reactant optimisation, final single point minus the optimiser's last-step
energy (kcal/mol).
USAGE  python3 c4_bands.py RESULTS_REPO OUT.tsv
"""
import hashlib, lzma, os, re, sys
RES, OUT = sys.argv[1], sys.argv[2]
Q = os.path.join(RES, "05_qmmm")
H = 627.5094740631
def stop(s): sys.exit("STOP: " + s)
def sums(p, prefix=""):
    d = {}
    for l in open(p):
        if l.strip(): d[prefix + l[66:].split("  (stored")[0].strip()] = l[:64]
    return d
FR = [l for l in open(os.path.join(RES, "phase2.2/A_v2.tsv")) if l.startswith("# frames")][0].split()[3:]
if len(FR) != 30: stop("A_v2 lists %d frames" % len(FR))
HV = {}
for l in open(os.path.join(Q, "19_ensemble_barriers/ensemble_barriers.tsv")):
    f = l.rstrip("\n").split("\t")
    if l.startswith("#") or f[0] == "frame": continue
    HV[f[0]] = f
if sorted(k for k, v in HV.items() if v[1] == "full_NAC" and int(k) >= 20000) != sorted(FR):
    stop("A_v2's frames are not the full_NAC frames numbered 20000 and above")
C5 = sorted(d[6:] for d in os.listdir(os.path.join(Q, "22_c5_nebci")) if d.startswith("frame_"))
if len(C5) != 8: stop("expected 8 C5 frames")
C5NEW = {}
for l in open(os.path.join(Q, "22_c5_nebci/c5_results.tsv")):
    f = l.rstrip("\n").split("\t")
    if f[0] != "frame": C5NEW[f[0]] = float(f[4])
X19 = sums(os.path.join(Q, "19_ensemble/CHECKSUMS_xz_originals.sha256"))
X22 = sums(os.path.join(Q, "22_c5_nebci/OUTPUT_CHECKSUMS.sha256"), "22_c5_nebci/")
def path_summary(t):
    m = list(re.finditer(r"PATH SUMMARY.*?\n(Image.*?\n)(.*?)\n\s*\n", t, re.S))
    if not m: return None
    rows = []
    for l in m[-1].group(2).splitlines():
        g = re.match(r"^\s*(\d+)\s+([-\d.]+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)", l)
        if g: rows.append((int(g.group(1)), float(g.group(2)), float(g.group(3)), "<= CI" in l))
    return rows
def log_last(p):
    """image energies (Eh, 8 decimals) of the last iteration in an ORCA .NEB.log"""
    blocks = [b for b in open(p).read().split(">") if re.search(r"iteration\s*:", b)]
    e = re.search(r"energy\s*:(.*)", blocks[-1])
    return [float(v) for v in e.group(1).split()]
def trj_energies(p):
    L = open(p).read().splitlines(); E = []; i = 0
    while i < len(L):
        if not L[i].strip(): i += 1; continue
        n = int(L[i].split()[0])
        g = re.search(r"\sE\s+(-?\d+\.\d+)\s*$", L[i + 1])
        if not g: stop("no energy in a comment line of " + p)
        E.append(g.group(1)); i += n + 2
    return E
def last_csv_energy(p):
    rows = [l.split(";") for l in open(p) if l.strip() and not l.lstrip().startswith("#")]
    return float(rows[-1][5])
rows = []
for f in FR:
    band = "22_c5_nebci" if f in C5 else "19_ensemble"
    D = os.path.join(Q, band, "frame_" + f)
    raw = lzma.open(os.path.join(D, "neb.out.xz")).read()
    rec = (X22 if f in C5 else X19).get("%s/frame_%s/neb.out" % (band, f))
    if rec is None or hashlib.sha256(raw).hexdigest() != rec: stop("frame %s: neb.out differs from its recorded sha256" % f)
    t = raw.decode("utf-8", "replace")
    if "THE NEB OPTIMIZATION HAS CONVERGED" not in t: stop("frame %s: band not converged" % f)
    ps = path_summary(t)
    WS = trj_energies(os.path.join(D, "neb_MEP.QMRegion_trj.xyz")); W = [float(x) for x in WS]
    LG = log_last(os.path.join(D, "neb.NEB.log"))
    if not ps or len(ps) != 10 or len(W) != 10 or len(LG) != 10: stop("frame %s: expected 10 images" % f)
    E = [r[2] for r in ps]
    if max(abs(W[k] - LG[k]) for k in range(10)) > 6e-9: stop("frame %s: the written band is not neb.NEB.log's last iteration" % f)
    nlog = len([x for x in open(os.path.join(D, "neb.NEB.log")).read().split(">") if re.search(r"iteration\s*:", x)])
    nint = len(re.findall(r"(?m)^Iteration:", open(os.path.join(D, "neb.interp")).read()))
    if nint != nlog + 1: stop("frame %s: neb.interp has %d iterations, neb.NEB.log %d (expected one more)" % (f, nint, nlog))
    if abs(W[0] - E[0]) > 6e-6 or abs(W[9] - E[9]) > 6e-6: stop("frame %s: an end point moved" % f)
    ci = [r[0] for r in ps if r[3]]
    if len(ci) != 1 or not 0 < ci[0] < 9: stop("frame %s: climbing image %s" % (f, ci))
    ci = ci[0]
    klow = min(range(ci + 1), key=lambda k: (E[k], k))
    kwl = min(range(ci + 1), key=lambda k: (W[k], k))
    rt = lzma.open(os.path.join(Q, "19_ensemble/frame_%s/reactant_opt.out.xz" % f), "rt", errors="replace").read()
    sp = float(re.findall(r"FINAL SINGLE POINT ENERGY \(QM/MM\)\s+(-?\d+\.\d+)", rt)[-1])
    if abs(sp - W[0]) > 1e-6: stop("frame %s: image 0 is not the Phase 1 reactant (%.2e Eh)" % (f, sp - W[0]))
    off = (sp - last_csv_energy(os.path.join(Q, "19_ensemble/frame_%s/reactant_opt-minimize-ener.csv" % f))) * H
    b0 = (E[ci] - E[0]) * H; bl = (E[ci] - E[klow]) * H; dep = (E[0] - E[klow]) * H
    ref = C5NEW[f] if f in C5 else float(HV[f][5])
    dci = [l for l in t[t.rfind("PATH SUMMARY"):].splitlines() if "<= CI" in l][0].split()[3]
    if abs(float(dci) - ref) > 1e-9 or abs(b0 - ref) > 0.011:
        stop("frame %s: barrier from image 0 %.3f, PATH SUMMARY %s, harvested %.2f" % (f, b0, dci, ref))
    rows.append(dict(f=f, band=band, ci=ci, klow=klow, kwl=kwl, dist=ps[klow][1], E0=E[0], Ek=E[klow], Eci=E[ci],
                     Ew0=WS[0], Ewk=WS[klow], lag=(W[klow] - E[klow]) * H, dep=dep, b0=b0, bl=bl, off=off))
deep = sorted(rows, key=lambda r: (-r["dep"], r["f"]))[:5]
nodip = [r for r in rows if r["klow"] == 0]
if len(nodip) != 2: stop("expected exactly 2 frames whose lowest image is image 0, found %d" % len(nodip))
S0 = {r["f"] for r in deep}; ND = {r["f"] for r in nodip}
def mean(x): return sum(x) / len(x)
def sd(x):
    m = mean(x); return (sum((v - m) ** 2 for v in x) / (len(x) - 1)) ** 0.5
b0s = [r["b0"] for r in rows]; bls = [r["bl"] for r in rows]; offs = [r["off"] for r in rows]; lags = [r["lag"] for r in rows]
out = ["# C4 band table (c4_bands.py) - the final band of each of the 30 reported frames",
       "# ci: climbing image; k_low: lowest image of the final band at or before it; kw_low: the same on the written band;",
       "#   dist_klow: band distance of k_low (A, all active atoms)",
       "# E_0, E_klow, E_ci: final-band energies (Eh, PATH SUMMARY); Ew_0, Ew_klow: energies of images 0 and k_low as written",
       "#   to disk (Eh, neb_MEP.QMRegion_trj.xyz), the starts of their re-optimisations;",
       "#   lag_klow = (Ew_klow - E_klow) x 627.5094740631 kcal/mol",
       "# depth = E_0 - E_klow; barrier_image0 = E_ci - E_0; barrier_klow = E_ci - E_klow (kcal/mol)",
       "# opt_offset: Phase 1 reactant optimisation, final single point minus the optimiser's last-step energy (kcal/mol)",
       "# runs: S0 = image 0 re-optimised (step 0, the five deepest dips); S2 = k_low re-optimised (step 2, every frame;",
       "#   for the two frames with k_low = 0 it is also their step-0 run)",
       "# frames with depth > 0.5 kcal/mol: %d of %d; barrier from image 0 %.2f +- %.2f, from k_low %.2f +- %.2f (mean +- sd)"
       % (sum(1 for r in rows if r["dep"] > 0.5), len(rows), mean(b0s), sd(b0s), mean(bls), sd(bls)),
       "# opt_offset %.3f to %.3f kcal/mol; lag_klow %.3f to %.3f kcal/mol; k_low differs between the written and the final"
       " band in %d frame(s)" % (min(offs), max(offs), min(lags), max(lags), sum(1 for r in rows if r["klow"] != r["kwl"])),
       "frame\tband\tci\tk_low\tkw_low\tdist_klow\tE_0\tE_klow\tE_ci\tEw_0\tEw_klow\tlag_klow\tdepth\tbarrier_image0\tbarrier_klow"
       "\topt_offset\truns"]
for r in rows:
    runs = "S0,S2" if r["f"] in S0 else ("S2(=S0, no dip)" if r["f"] in ND else "S2")
    out.append("%s\t%s\t%d\t%d\t%d\t%.3f\t%.5f\t%.5f\t%.5f\t%s\t%s\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%s" % (
        r["f"], r["band"], r["ci"], r["klow"], r["kwl"], r["dist"], r["E0"], r["Ek"], r["Eci"], r["Ew0"], r["Ewk"], r["lag"], r["dep"],
        r["b0"], r["bl"], r["off"], runs))
open(OUT, "w").write("\n".join(out) + "\n")
print("c4_bands: 30 frames; step 0: %s (deepest dips) and %s (no dip)" % (" ".join(r["f"] for r in deep), " ".join(sorted(ND))))
