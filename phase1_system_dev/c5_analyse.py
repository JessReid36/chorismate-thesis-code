#!/usr/bin/env python3
"""
c5_analyse.py - PHASE1_AUDIT_CHECKLIST.md item C5: the eight NEB-CI restarts of the post-cut frames that
had run NEB-TS. Applies the criteria in 05_qmmm/22_c5_nebci/C5_CRITERIA.txt (fixed before submission,
committed in 86e8f43) to every frame whose restart has finished; frames still running are reported as
pending, so the script can be rerun as they complete.

Reads, per frame:
  22_c5_nebci/frame_N/   neb.out (last CI-NEB convergence table, last PATH SUMMARY), neb.NEB.log,
                         neb_NEB-CI_converged.QMRegion.xyz (the new climbing image)
  19_ensemble/frame_N/   neb.out (the original run's last PATH SUMMARY, for the end-point energies) and
                         neb_NEB-CI_converged.QMRegion.xyz (the loosely converged climbing image)
  19_ensemble_barriers/ensemble_barriers.tsv  (the tabulated barrier)
The barrier is the '<= CI' row's dE of the LAST path summary in neb.out, as for every tabulated barrier.

Criteria (C5_CRITERIA.txt):
  1. THE NEB OPTIMIZATION HAS CONVERGED, and the final table meets the NEB-CI defaults
     MAX(|FCI|) <= 5.0e-4, RMS(FCI) <= 2.5e-4, MAX(|Fp|) <= 5.0e-3, RMS(Fp) <= 2.5e-3 Eh/bohr
  2. the barrier varies by at most 0.1 kcal/mol over the last 5 iterations of neb.NEB.log
  3. expected: a fall of 0.2-0.8 kcal/mol from the tabulated value; a fall above 2 kcal/mol, or a
     rise, is inspected for environment rearrangement (C4) before it is accepted
  4. the new climbing image within about 0.05 A RMSD (QM region) of the old one
  5. image 0 and the last image keep their energies: within 1e-5 Eh
Writes C5_REPORT.txt and c5_results.tsv into 22_c5_nebci/.

USAGE   python3 c5_analyse.py [system_development folder]
Python 3.6, standard library only (runs on the hpc1 login node).
"""
import math
import re
import sys
from pathlib import Path

SD = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("/home/18660916/system_development")
C = SD / "05_qmmm" / "22_c5_nebci"
E = SD / "05_qmmm" / "19_ensemble"
TSV = SD / "05_qmmm" / "19_ensemble_barriers" / "ensemble_barriers.tsv"
FRAMES = ["20000", "21634", "23268", "24883", "26495", "33320", "34991", "43738"]
H2KCAL = 627.5094740631
TOL = {"MAX(|FCI|)": 5.0e-4, "RMS(FCI)": 2.5e-4, "MAX(|Fp|)": 5.0e-3, "RMS(Fp)": 2.5e-3}


def last_path_summary(t):
    m = list(re.finditer(r"PATH SUMMARY.*?\n(Image.*?\n)(.*?)\n\s*\n", t, re.S))
    if not m:
        return None
    rows = []
    for l in m[-1].group(2).splitlines():
        f = re.match(r"^\s*(\d+)\s+([-\d.]+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)", l)
        if f:
            rows.append((int(f.group(1)), float(f.group(3)), float(f.group(4)), "<= CI" in l))
    return rows


def last_conv_table(t):
    i = t.rfind("CI-NEB convergence")
    if i < 0:
        return {}
    out = {}
    for k, v, tol, yes in re.findall(r"(RMS\(Fp\)|MAX\(\|Fp\|\)|RMS\(FCI\)|MAX\(\|FCI\|\))\s+([\d.Ee+-]+)\s+([\d.Ee+-]+)\s+(\w+)",
                                     t[i:i + 3000])[:4]:
        out[k] = (float(v), float(tol), yes)
    return out


def neb_log_barriers(p):
    out = []
    for b in p.read_text().split("\n>"):
        dd = {}
        for l in b.splitlines():
            if ":" in l:
                k, v = l.split(":", 1)
                dd[k.strip()] = v.strip()
        if "iteration" in dd and "barrier" in dd:
            out.append(float(dd["barrier"].split()[0]) * H2KCAL)
    return out


def read_xyz(p):
    L = p.read_text().splitlines()
    n = int(L[0].split()[0])
    return [[float(v) for v in l.split()[1:4]] for l in L[2:2 + n]]


def rmsd(A, B):
    """RMSD after optimal superposition (Horn's quaternion method)"""
    n = len(A)
    ca = [sum(p[k] for p in A) / n for k in range(3)]
    cb = [sum(p[k] for p in B) / n for k in range(3)]
    a = [[p[k] - ca[k] for k in range(3)] for p in A]
    b = [[p[k] - cb[k] for k in range(3)] for p in B]
    S = [[sum(a[i][r] * b[i][c] for i in range(n)) for c in range(3)] for r in range(3)]
    (xx, xy, xz), (yx, yy, yz), (zx, zy, zz) = S
    N = [[xx + yy + zz, yz - zy, zx - xz, xy - yx], [yz - zy, xx - yy - zz, xy + yx, zx + xz],
         [zx - xz, xy + yx, -xx + yy - zz, yz + zy], [xy - yx, zx + xz, yz + zy, -xx - yy + zz]]
    lam = max(abs(v) for row in N for v in row) * 4 + 1.0
    v = [1.0, 0.0, 0.0, 0.0]
    for _ in range(2000):
        w = [sum((N[i][j] + (lam if i == j else 0.0)) * v[j] for j in range(4)) for i in range(4)]
        s = math.sqrt(sum(x * x for x in w))
        v = [x / s for x in w]
    emax = sum(v[i] * sum(N[i][j] * v[j] for j in range(4)) for i in range(4))
    ga, gb = sum(x * x for p in a for x in p), sum(x * x for p in b for x in p)
    return math.sqrt(max(0.0, (ga + gb - 2.0 * emax) / n))


def tabulated():
    hdr, out = None, {}
    for l in TSV.read_text().splitlines():
        if not l.strip() or l.startswith("#"):
            continue
        f = l.split("\t")
        if hdr is None:
            hdr = f
            continue
        r = dict(zip(hdr, f))
        out[r["frame"]] = float(r["barrier"]) if r.get("barrier") else None
    return out


def main():
    tab = tabulated()
    rep, rows = [], []
    say = rep.append
    say("C5 report (c5_analyse.py) - criteria from 22_c5_nebci/C5_CRITERIA.txt (commit 86e8f43)\n")
    for f in FRAMES:
        F, O = C / "frame_{}".format(f), E / "frame_{}".format(f)
        t = (F / "neb.out").read_text(errors="replace") if (F / "neb.out").exists() else ""
        done = "THE NEB OPTIMIZATION HAS CONVERGED" in t and "ORCA TERMINATED NORMALLY" in t
        if not done:
            ended = (F / "neb.pbs.out").exists() and "orca_exit=" in (F / "neb.pbs.out").read_text(errors="replace")
            say("== frame {}: {}".format(f, "ENDED WITHOUT CONVERGING - inspect neb.out" if ended else "pending (not finished)"))
            rows.append({"frame": f, "status": "failed" if ended else "pending"})
            continue
        ps, ps0 = last_path_summary(t), last_path_summary((O / "neb.out").read_text(errors="replace"))
        ci = [r for r in ps if r[3]]
        ci0 = [r for r in ps0 if r[3]]
        tb = tab.get(f)
        new = ci[0][2]
        conv = last_conv_table(t)
        c1 = len(conv) == 4 and all(conv[k][0] <= TOL[k] and abs(conv[k][1] - TOL[k]) < 1e-12 for k in TOL)
        bl = neb_log_barriers(F / "neb.NEB.log")
        spread = max(bl[-5:]) - min(bl[-5:]) if len(bl) >= 5 else float("nan")
        dB = new - tb
        if dB > 0:
            c3 = "INSPECT (rise)"
        elif dB < -2.0:
            c3 = "INSPECT (fall above 2 kcal/mol; C4)"
        elif -0.8 <= dB <= -0.2:
            c3 = "as expected"
        else:
            c3 = "outside the expected 0.2-0.8 fall (reported)"
        r = rmsd(read_xyz(O / "neb_NEB-CI_converged.QMRegion.xyz"), read_xyz(F / "neb_NEB-CI_converged.QMRegion.xyz"))
        e0 = abs(ps[0][1] - ps0[0][1])
        eN = abs(ps[-1][1] - ps0[-1][1])
        say("== frame {}: converged in {} logged iterations".format(f, len(bl)))
        say("   1. final table: " + ", ".join("{} {:.2e} (tol {:.1e})".format(k, conv[k][0], conv[k][1]) for k in TOL if k in conv)
            + " -> {}".format("PASS" if c1 else "FAIL"))
        say("   2. barrier over the last 5 logged iterations varies by {:.3f} kcal/mol -> {}".format(
            spread, "PASS" if spread <= 0.1 else "FAIL"))
        say("   3. barrier {:.2f} -> {:.2f} kcal/mol ({:+.2f}); climbing image {} -> {} -> {}".format(
            tb, new, dB, ci0[0][0], ci[0][0], c3))
        say("   4. climbing image moved {:.3f} A RMSD (QM region) -> {}".format(r, "PASS" if r <= 0.05 else "beyond about 0.05 A (reported)"))
        say("   5. end points: image 0 {:.1e} Eh, last image {:.1e} Eh -> {}".format(e0, eN, "PASS" if max(e0, eN) <= 1e-5 else "FAIL"))
        rows.append({"frame": f, "status": "converged", "iterations": len(bl), "barrier_tab": tb, "barrier_new": new,
                     "delta": dB, "ci_old": ci0[0][0], "ci_new": ci[0][0], "MAX_FCI": conv.get("MAX(|FCI|)", (float("nan"),))[0],
                     "spread_last5": spread, "ci_rmsd": r, "e0_diff": e0, "eN_diff": eN, "crit1": c1,
                     "crit3": c3})
    done = [x for x in rows if x["status"] == "converged"]
    if done:
        ds = sorted(x["delta"] for x in done)
        say("\nconverged {} of 8; barrier change median {:+.2f}, range {:+.2f} to {:+.2f} kcal/mol "
            "(expected from the 22 NEB-CI frames: median -0.38, range -0.19 to -4.95)".format(
                len(done), ds[len(ds) // 2] if len(ds) % 2 else (ds[len(ds) // 2 - 1] + ds[len(ds) // 2]) / 2, ds[0], ds[-1]))
    keys = []
    for x in rows:
        for k in x:
            if k not in keys:
                keys.append(k)
    lines = ["\t".join(keys)] + ["\t".join(("{:.6g}".format(x[k]) if isinstance(x.get(k), float) else str(x.get(k, "")))
                                           for k in keys) for x in rows]
    (C / "c5_results.tsv").write_text("\n".join(lines) + "\n")
    (C / "C5_REPORT.txt").write_text("\n".join(rep) + "\n")
    print("\n".join(rep))
    return 0


if __name__ == "__main__":
    sys.exit(main())
