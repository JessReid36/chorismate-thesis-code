#!/usr/bin/env python3
"""
step11e_autocorr_nac.py - two read-only analyses of the production trajectory,
both needed before the Phase-1 barrier ensemble can be quoted.

(A) AUTOCORRELATION - are the selected frames independent?

    The barrier ensemble reports a mean and a standard error. A standard error
    computed from correlated samples is too small: it claims precision the data
    does not have. The frames were chosen at ~1670 ps spacing, which is very
    likely enough, but nothing has demonstrated it.

    This computes the integrated autocorrelation time (tau_int) of the two
    quantities that define catalytic competence - the forming C1-C6 distance and
    the Arg90-O13 contact - and compares it with the frame spacing. It also
    reports the effective sample size, which is the honest denominator for the
    standard error.

(B) NAC ORBITAL ALIGNMENT - is the distance criterion the whole criterion?

    Freindorf & Kraka (2019) define chorismate near-attack conformers as
    structures with the forming bond within the van der Waals contact distance
    of 3.7 A AND with the pi-orbitals of the reacting carbons pointing toward
    one another, allowing sufficient overlap for bond formation.

    The step-12a filter tests only the distance. This computes the alignment for
    every selected frame, so the write-up can either state that the frames
    satisfy the full definition or acknowledge that only the distance clause was
    enforced.

    Alignment is measured as the angle between the C1->C6 vector and the local
    pi-axis at each carbon, taken as the normal to the plane of that carbon and
    its two heavy neighbours. Perfect head-on overlap gives 0 deg at both ends.

Usage
-----
  python3 step11e_autocorr_nac.py
  python3 step11e_autocorr_nac.py --stride 1 --max-lag 4000

Read-only. Uses the stdlib NetCDF3 reader pattern from step11c/step12a because
cpptraj is broken cluster-wide. Single-threaded to avoid the login-node
OpenBLAS segfault - set OPENBLAS_NUM_THREADS=1 before running.
"""

import argparse
import os
import re
import struct
import sys
import numpy as np

ARG_N = {"NE", "NH1", "NH2"}


def parse_prmtop(path):
    flags, fmts, cur = {}, {}, None
    for ln in open(path, errors="replace"):
        ln = ln.rstrip("\n")
        if ln.startswith("%FLAG"):
            cur = ln.split()[1]; flags[cur] = []
        elif ln.startswith("%FORMAT"):
            fmts[cur] = ln[ln.find("(") + 1:ln.rfind(")")]
        elif ln.startswith("%"):
            continue
        elif cur is not None:
            flags[cur].append(ln)

    def width(f):
        m = re.match(r"\s*\d*[aAiIeEfFgG](\d+)", f)
        return int(m.group(1)) if m else None

    def toks(name, cast):
        w = width(fmts.get(name, "")); out = []
        for ln in flags.get(name, []):
            if w:
                out += [ln[i:i + w].strip() for i in range(0, len(ln.rstrip()), w)
                        if ln[i:i + w].strip()]
            else:
                out += ln.split()
        return [cast(x) for x in out]

    return (toks("POINTERS", int)[0], toks("ATOM_NAME", str),
            toks("RESIDUE_LABEL", str), toks("RESIDUE_POINTER", int))


def nc_open(path):
    data = open(path, "rb").read(1 << 20)
    fsz = os.path.getsize(path)
    off64 = (data[3] == 2)
    pos = [4]

    def u32():
        v = struct.unpack_from(">I", data, pos[0])[0]; pos[0] += 4; return v

    def off():
        v = struct.unpack_from(">Q" if off64 else ">I", data, pos[0])[0]
        pos[0] += 8 if off64 else 4; return v

    def nm():
        n = u32(); s = data[pos[0]:pos[0] + n].decode("ascii", "replace")
        pos[0] += n + ((4 - (n % 4)) % 4); return s

    numrecs = u32(); dims = []; tag = u32(); ne = u32()
    if tag == 0x0A:
        for _ in range(ne):
            dn = nm(); dl = u32(); dims.append((dn, dl))

    def skip_att():
        t = u32(); n = u32()
        if t == 0x0C:
            for _ in range(n):
                nm(); tp = u32(); k = u32()
                nb = k * {1: 1, 2: 1, 3: 2, 4: 4, 5: 4, 6: 8}.get(tp, 1)
                pos[0] += nb + ((4 - (nb % 4)) % 4)

    skip_att(); rec = []; tag = u32(); nv = u32()
    if tag == 0x0B:
        for _ in range(nv):
            vn = nm(); nd = u32(); dimids = [u32() for _ in range(nd)]
            skip_att(); tp = u32(); vs = u32(); bg = off()
            if nd > 0 and dims[dimids[0]][1] == 0:
                rec.append({"name": vn, "vsize": vs, "begin": bg})
    dimd = {n: l for n, l in dims}
    recsize = sum(v["vsize"] for v in rec)
    frames = numrecs if numrecs != 0xFFFFFFFF else \
        (fsz - min(v["begin"] for v in rec)) // recsize
    return {"frames": frames, "recsize": recsize, "rec": rec, "atom": dimd.get("atom")}


def integrated_tau(x, max_lag):
    """
    Integrated autocorrelation time by the automatic-windowing rule of Sokal:
    sum the normalised autocorrelation until the lag exceeds c*tau, c = 5.
    tau_int = 1 + 2*sum_{k=1..W} rho(k); effective sample size = N / (2*tau_int).
    """
    x = np.asarray(x, float)
    x = x - x.mean()
    n = len(x)
    var = np.dot(x, x) / n
    if var <= 0:
        return float("nan"), float("nan")
    # FFT autocorrelation
    size = 1
    while size < 2 * n:
        size *= 2
    f = np.fft.rfft(x, size)
    acf = np.fft.irfft(f * np.conjugate(f), size)[:n].real
    acf /= acf[0]

    tau = 1.0
    window = min(max_lag, n - 1)
    for k in range(1, window + 1):
        tau += 2.0 * acf[k]
        if k >= 5 * tau:
            window = k
            break
    ess = n / (2.0 * tau) if tau > 0 else float("nan")
    return tau, ess


def pi_axis(p_c, p_a, p_b):
    """Local pi direction at a carbon: normal to the plane it makes with its two
    heavy neighbours. Sign is arbitrary; the caller folds the angle to <= 90."""
    v1 = p_a - p_c
    v2 = p_b - p_c
    n = np.cross(v1, v2)
    nn = np.linalg.norm(n)
    return n / nn if nn > 0 else None


def main():
    ap = argparse.ArgumentParser()
    home = os.path.expanduser("~")
    ap.add_argument("--prmtop",
                    default=f"{home}/system_development/03_amber/tleap_build/complex_solvated.prmtop")
    ap.add_argument("--traj",
                    default=f"{home}/system_development/04_amber_md/10c_production/prod.nc")
    ap.add_argument("--manifest",
                    default=f"{home}/system_development/05_qmmm/12_frame_selection/selection_manifest.tsv")
    ap.add_argument("--stride", type=int, default=1,
                    help="trajectory stride for the autocorrelation (1 = every frame)")
    ap.add_argument("--max-lag", type=int, default=4000)
    ap.add_argument("--site", default="CHA#2")
    args = ap.parse_args()

    natom, names, labels, resptr = parse_prmtop(args.prmtop)
    nres = len(labels)
    starts = [p - 1 for p in resptr] + [natom]
    res_of = np.empty(natom, int)
    for ri in range(nres):
        res_of[starts[ri]:starts[ri + 1]] = ri

    argN = np.array([a for a in range(natom)
                     if labels[res_of[a]] == "ARG" and names[a] in ARG_N])
    cha = []
    for ri in range(nres):
        if labels[ri] == "CHA":
            cha.append((ri, {names[a]: a for a in range(starts[ri], starts[ri + 1])}))
    if not cha:
        sys.exit("FAIL no CHA residues")

    site_idx = int(args.site.split("#")[1]) - 1
    ri, amap = cha[site_idx]
    print(f"site {args.site} = residue {ri+1}\n")

    need = ["C1", "C6", "O3", "C4"]
    for k in need:
        if k not in amap:
            sys.exit(f"FAIL atom {k} not found in {args.site}")

    maxidx = int(max(argN.max(), max(amap.values()))) + 1
    H = nc_open(args.traj)
    if H["atom"] != natom:
        sys.exit(f"FAIL atom count mismatch: traj {H['atom']} vs prmtop {natom}")
    cvar = next(v for v in H["rec"] if v["name"] == "coordinates")
    lvar = next((v for v in H["rec"] if v["name"] == "cell_lengths"), None)

    def frame(i):
        o = cvar["begin"] + i * H["recsize"]
        xyz = np.array(np.memmap(args.traj, dtype=">f4", mode="r", offset=o,
                                 shape=(natom, 3))[:maxidx], float)
        box = np.array([1e9, 1e9, 1e9])
        if lvar is not None:
            bo = lvar["begin"] + i * H["recsize"]
            box = np.array(np.memmap(args.traj, dtype=">f8", mode="r", offset=bo,
                                     shape=(3,)), float)
        return xyz, box

    def mind(a, pts, box):
        d = pts - a
        d -= np.round(d / box) * box
        return float(np.sqrt((d * d).sum(1)).min())

    idxs = list(range(0, H["frames"], args.stride))
    print(f"reading {len(idxs)} of {H['frames']} frames (stride {args.stride}) ...")
    form, contact = [], []
    for i in idxs:
        x, box = frame(i)
        form.append(float(np.linalg.norm(x[amap["C6"]] - x[amap["C1"]])))
        contact.append(mind(x[amap["O3"]], x[argN], box))

    # ---------------------------------------------------------- (A) autocorr
    print("\n(A) AUTOCORRELATION")
    dt_ps = args.stride * 1.0     # 1 ps per written frame (ntwx 500 at dt 2 fs)
    print(f"    frame interval {dt_ps:.0f} ps\n")
    print(f"    {'quantity':<22}{'tau_int (frames)':>18}{'tau_int (ps)':>15}"
          f"{'eff. N':>10}")
    print("    " + "-" * 65)
    taus = {}
    for lbl, series in (("forming C1-C6", form), ("Arg90-O13 contact", contact)):
        tau, ess = integrated_tau(series, args.max_lag)
        taus[lbl] = tau * dt_ps
        print(f"    {lbl:<22}{tau:>18.1f}{tau*dt_ps:>15.1f}{ess:>10.0f}")

    # compare with the selection spacing
    try:
        rows = [l.split("\t") for l in open(args.manifest).read().splitlines()
                if not l.startswith("#") and l.strip()]
        frames_sel = sorted(int(r[1]) for r in rows)
        gaps = np.diff(frames_sel)
        worst = max(taus.values())
        print(f"\n    selected frames: {len(frames_sel)}, "
              f"spacing {gaps.min()}-{gaps.max()} ps (mean {gaps.mean():.0f})")
        print(f"    slowest tau_int: {worst:.1f} ps")
        ratio = gaps.min() / worst if worst > 0 else float("inf")
        print(f"    minimum spacing is {ratio:.1f} x the slowest correlation time")
        if ratio >= 5:
            print("    -> frames are effectively independent; the sem may be "
                  "quoted as sd/sqrt(n)")
        elif ratio >= 2:
            print("    -> marginally independent; state the ratio explicitly and "
                  "treat the sem as a lower bound")
        else:
            print("    -> NOT independent; the sem from sd/sqrt(n) would be "
                  "optimistic. Widen the spacing or inflate the error bar by "
                  "sqrt(2*tau/spacing).")
    except Exception as e:
        print(f"    (manifest comparison skipped: {e})")

    # ------------------------------------------------------------- (B) NAC
    print("\n(B) NAC ORBITAL ALIGNMENT at the selected frames")
    print("    Freindorf & Kraka (2019): a NAC requires the forming bond within")
    print("    3.7 A AND the pi-orbitals of the reacting carbons pointing toward")
    print("    one another. step12a tests the distance only.\n")

    # heavy neighbours of C1 and C6, by connectivity from the first frame
    x0, _ = frame(0)
    heavy = [a for a in amap.values() if not names[a].startswith("H")]

    def neighbours(c, k=2):
        d = [(np.linalg.norm(x0[a] - x0[c]), a) for a in heavy if a != c]
        d.sort()
        return [a for _, a in d[:k]]

    n1 = neighbours(amap["C1"])
    n6 = neighbours(amap["C6"])
    print(f"    pi-axis at C1 from neighbours {[names[a] for a in n1]}")
    print(f"    pi-axis at C6 from neighbours {[names[a] for a in n6]}\n")

    print(f"    {'frame':>7}{'form (A)':>11}{'angle C1':>11}{'angle C6':>11}"
          f"{'verdict':>12}")
    print("    " + "-" * 54)
    try:
        for r in rows:
            fr = int(r[1])
            x, _ = frame(fr)
            v = x[amap["C6"]] - x[amap["C1"]]
            dist = np.linalg.norm(v)
            v = v / dist
            a1 = pi_axis(x[amap["C1"]], x[n1[0]], x[n1[1]])
            a6 = pi_axis(x[amap["C6"]], x[n6[0]], x[n6[1]])
            ang1 = np.degrees(np.arccos(min(1.0, abs(float(np.dot(a1, v))))))
            ang6 = np.degrees(np.arccos(min(1.0, abs(float(np.dot(a6, v))))))
            ok = dist < 3.7 and ang1 < 45 and ang6 < 45
            print(f"    {fr:>7}{dist:>11.3f}{ang1:>11.1f}{ang6:>11.1f}"
                  f"{'NAC' if ok else 'distance only':>12}")
        print("\n    Angles are between the forming-bond vector and the local pi")
        print("    axis at each carbon; 0 deg is perfect head-on overlap. A cutoff")
        print("    of 45 deg is used here as a working threshold - if the frames")
        print("    cluster well below it the distinction is moot, and if they")
        print("    straddle it the threshold must be justified or the orbital")
        print("    clause dropped from the stated criterion with an explanation.")
    except Exception as e:
        print(f"    (skipped: {e})")


if __name__ == "__main__":
    main()
