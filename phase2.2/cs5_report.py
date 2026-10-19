#!/usr/bin/env python3
"""
cs5_report.py - Stage 5's descriptive report (STAGE5_CRITERIA.txt section 8). No verdict depends on it; the checks
that must pass are in cs5_align.py, cs5_grid.py check / ljcheck and cs5_amatrix.py assemble.

  python3 cs5_report.py <stage5_dir> <results repo> <out STAGE5_REPORT.txt>

stage5_dir holds extracted/EXTRACT.txt, aligned/ (ALIGNMENT.txt, frame_<F>/...), sites.tsv, L.tsv and, once the A
matrix exists, A.tsv and A_REPORT.txt. From the results repo it reads phase2.2/grid_v2.tsv, A_v2.tsv and the
development alignment phase2.2/aligned/ (for the informative comparison with A_v2: grid_v2's sites are carried into
Stage 5's frame of reference by a core fit of the two references' reactants, and each is paired with the nearest
Stage 5 site within 0.5 A; a is compared over the frames the two share).
"""
import os, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cs5_grid as G
import cs5_align as AL

O_CARB = ("O1", "O2", "O5", "O6")

def matrix(p):
    com = [l for l in Path(p).read_text().splitlines() if l.startswith("#")]
    return np.loadtxt(p, comments="#", ndmin=2), com

def main(sd, res, out):
    sd, res = Path(sd), Path(res)
    if Path(out).exists(): sys.exit("STOP: %s exists" % out)
    W = []
    say = W.append
    say("Stage 5 report (cs5_report.py; descriptive, STAGE5_CRITERIA.txt section 8); cs5_report.py sha256 %s" % G.sha(__file__))
    ex = [l.split("\t") for l in (sd / "extracted/EXTRACT.txt").read_text().splitlines() if not l.startswith("#")][1:]
    frames = G.frame_names(sd / "aligned")
    say("")
    say("FRAMES: %d (%s)" % (len(frames), " ".join(frames)))
    pirc = [r[0] for r in ex if r[2] == "P" and r[3].startswith("irc/")]
    say("  P from the IRC's product end (P_IRC did not pass): %s" % (" ".join(pirc) if pirc else "none"))
    al = (sd / "aligned/ALIGNMENT.txt").read_text().splitlines()
    say("  " + [l for l in al if l.startswith("reference:")][0])
    say("  " + [l for l in al if l.startswith("RMSD on the core after alignment")][0])
    k = al.index("frame\tRMSD_before_A\tRMSD_after_A")
    say("  core RMSD after alignment per frame (A): " + " ".join("%s %s" % tuple(l.split("\t")[0::2]) for l in al[k + 1:] if l.strip()))
    rows, com = G.read_sites(sd / "sites.tsv")
    say("")
    say("GRID (sites.tsv)")
    for l in com:
        if l.startswith(("# anion/neutral layer", "# cation layer", "# sites by admitted", "# nearest neighbour", "# parameters")):
            say("  " + l[2:])
    lay = np.array([r["layer"] for r in rows])
    d = {n: np.array([float(r["d_" + n]) for r in rows]) for n in G.NAMES}
    for L_ in ("anion", "cation"):
        m = lay == L_
        say("  %-6s closest site to O3 %.3f A (%d within 3.6 A); to the carboxylate oxygens %s" % (
            L_, d["O3"][m].min(), int((d["O3"][m] < 3.6).sum()), ", ".join("%s %.3f" % (o, d[o][m].min()) for o in O_CARB)))
    Lm, Lc = matrix(sd / "L.tsv")
    say("")
    say("L (L.tsv), kcal/mol: range %.3f to %.3f" % (Lm.min(), Lm.max()))
    for l in Lc:
        if "summed |E_LJ|" in l: say("  " + l[2:].strip())
    for L_ in ("anion", "cation"):
        v = Lm[lay == L_].mean(1)
        say("  %-6s per-site mean over frames: %.3f to %.3f" % (L_, v.min(), v.max()))
    say("")
    if not (sd / "A.tsv").is_file():
        say("A: not yet assembled")
    else:
        A, Ac = matrix(sd / "A.tsv")
        if A.shape != Lm.shape: sys.exit("STOP: A.tsv and L.tsv differ in shape")
        say("A (A.tsv), kcal/mol per e: range %.3f to %.3f; a < 0 (a cation lowers the barrier) in %.1f%% of entries (sites x frames)" % (
            A.min(), A.max(), 100.0 * (A < 0).mean()))
        say("  " + [l for l in (sd / "A_REPORT.txt").read_text().splitlines() if l.startswith("# checks")][0][2:])
        for L_ in ("anion", "cation"):
            m = lay == L_; Am = A[m]
            agree_neg = int((Am < 0).all(1).sum()); agree_pos = int((Am > 0).all(1).sum())
            say("  %-6s per-site mean a %.3f to %.3f; every frame a < 0 at %d sites, a > 0 at %d, mixed at %d of %d" % (
                L_, Am.mean(1).min(), Am.mean(1).max(), agree_neg, agree_pos, m.sum() - agree_neg - agree_pos, m.sum()))
        # informative: against A_v2 (development grid and alignment)
        _, X5 = AL.read24(str(sd / "aligned" / ("frame_" + frames[0]) / G.STATES["R"]))
        dref = [l for l in (res / "phase2.2/aligned/ALIGNMENT.txt").read_text().splitlines() if l.startswith("reference:")][0].split()[2].rstrip(",")
        _, Xd = AL.read24(str(res / "phase2.2/aligned" / ("frame_" + dref) / G.STATES["R"]))
        Rm, t = AL.kabsch(Xd, X5, AL.CORE)                       # development frame -> Stage 5 frame
        gv, _ = G.read_sites(res / "phase2.2/grid_v2.tsv")
        Q2 = G.sites_xyz(gv) @ Rm + t
        A2 = np.loadtxt(res / "phase2.2/A_v2.tsv", comments="#", ndmin=2)
        f2 = [l for l in (res / "phase2.2/A_v2.tsv").read_text().splitlines() if l.startswith("# frames")][0].split()[3:]
        common = [f for f in frames if f in f2]
        Q5 = G.sites_xyz(rows)
        D = G.pair_dist(Q2, Q5); j = D.argmin(1); ok = D[np.arange(len(Q2)), j] <= 0.5
        if common and ok.any():
            a5 = A[j[ok]][:, [frames.index(f) for f in common]].ravel(); a2 = A2[ok][:, [f2.index(f) for f in common]].ravel()
            say("  informative, against A_v2: %d of %d grid_v2 sites have a Stage 5 site within 0.5 A (development reference frame %s"
                % (ok.sum(), len(Q2), dref))
            say("    carried over by a core fit); over %d frames in common: Pearson r %.3f, same sign %.1f%%, mean |difference| %.3f"
                % (len(common), np.corrcoef(a5, a2)[0, 1], 100.0 * (np.sign(a5) == np.sign(a2)).mean(), np.abs(a5 - a2).mean()))
        else:
            say("  informative, against A_v2: no pairs (frames in common %d, sites paired %d)" % (len(common), int(ok.sum())))
    Path(out).write_text("\n".join(W) + "\n")
    print("\n".join(W))

if __name__ == "__main__":
    if len(sys.argv) != 4: sys.exit(__doc__)
    main(*sys.argv[1:])
