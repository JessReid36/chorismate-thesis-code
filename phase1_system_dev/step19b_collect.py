#!/usr/bin/env python3
"""
step19b_collect.py - collect the ensemble QM/MM barriers and report the
distribution, so the Phase-1 barrier can be quoted as a mean with a spread
rather than as a single path.

Reads the converged NEB climbing-image barrier from each frame directory under
05_qmmm/19_ensemble/, joins it to the step-12a selection manifest, and reports:

  * mean, standard deviation, and range of the barrier across frames
  * the standard error on the mean, which is what a comparison with a literature
    average or with experiment should be quoted against
  * where the fully characterised frame (820) sits within the distribution
  * correlation of the barrier with the Arg90-O13 contact distance and with the
    near-attack C6-C1 distance recorded at selection time

The last of these is a real test rather than a decoration. The literature
attributes transition-state stabilisation in BsCM largely to Arg90, so if the
barrier varies across frames mainly through electrostatic contact, a negative
correlation between barrier and Arg90 proximity is the expected signature.
Finding it would explain the frame-to-frame spread mechanistically; not finding
it would say the spread comes from somewhere else and is worth understanding
before the mean is quoted.

BARRIER SOURCE
--------------
The converged NEB climbing image, read from the ORCA output. This is the
validated proxy: on frame 820 it gave 15.94 kcal/mol against the fully
characterised reduced-region value of 16.00, agreeing to 0.06. Only converged
bands are counted; an unconverged NEB is reported and excluded rather than
silently averaged in.

Usage
-----
  python3 step19b_collect.py
  python3 step19b_collect.py --ensdir /path/to/19_ensemble --manifest /path/to/selection_manifest.tsv
"""

import argparse
import os
import re
import sys
from pathlib import Path

HARTREE2KCAL = 627.5094740631
REFERENCE_FRAME = 820
REFERENCE_FULL = 16.00     # fully characterised reduced-region barrier, frame 820
REFERENCE_NEB = 15.94      # its NEB-CI proxy, for the validation statement
EXPERIMENT_DH = 12.7       # experimental activation enthalpy, kcal/mol


def read_manifest(path):
    rows = {}
    for line in Path(path).read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        f = line.split("\t")
        rows[int(f[1])] = {
            "ps": float(f[2]),
            "form": float(f[5]),
            "r": float(f[6]),
            "arg90": float(f[7]),
        }
    return rows


def neb_barrier(out_path):
    """
    Converged climbing-image barrier, in kcal/mol.

    Returns (value, 'ok') or (None, reason). The energies are taken from the
    frozen/active energy table that ORCA prints for the converged band; the
    trajectory comment lines are NOT used, because Knarr writes a
    spring-augmented internal energy there that is not the physical barrier.
    """
    if not out_path.exists():
        return None, "no output"
    txt = out_path.read_text(errors="replace")

    if "ORCA TERMINATED NORMALLY" not in txt:
        return None, "did not terminate normally"
    if not re.search(r"THE NEB OPTIMIZATION HAS CONVERGED", txt, re.I):
        return None, "band not converged"

    # last block of image energies: lines ending in '@' carry the HEI column
    hits = re.findall(r"^\s*\d+\s+\S+\s+([-\d.]+)\s+.*@\s*$", txt, re.M)
    if hits:
        try:
            return float(hits[-1]), "ok"
        except ValueError:
            pass

    # fall back to the explicit climbing-image statement
    m = re.findall(r"climbing image.*?([\d.]+)\s*kcal", txt, re.I | re.S)
    if m:
        return float(m[-1]), "ok"

    return None, "converged but no barrier found in output"


def stats(vals):
    n = len(vals)
    mean = sum(vals) / n
    if n > 1:
        sd = (sum((v - mean) ** 2 for v in vals) / (n - 1)) ** 0.5
        sem = sd / n ** 0.5
    else:
        sd = sem = float("nan")
    return mean, sd, sem


def pearson(x, y):
    n = len(x)
    if n < 3:
        return float("nan")
    mx, my = sum(x) / n, sum(y) / n
    sxy = sum((a - mx) * (b - my) for a, b in zip(x, y))
    sxx = sum((a - mx) ** 2 for a in x)
    syy = sum((b - my) ** 2 for b in y)
    if sxx <= 0 or syy <= 0:
        return float("nan")
    return sxy / (sxx * syy) ** 0.5


def main():
    ap = argparse.ArgumentParser()
    home = os.path.expanduser("~")
    ap.add_argument("--ensdir", default=f"{home}/system_development/05_qmmm/19_ensemble")
    ap.add_argument("--manifest",
                    default=f"{home}/system_development/05_qmmm/12_frame_selection/"
                            f"selection_manifest.tsv")
    ap.add_argument("--out", default="ensemble_barriers.tsv")
    args = ap.parse_args()

    man = read_manifest(args.manifest)
    ensdir = Path(args.ensdir)

    results, skipped = [], []

    # frame 820 is carried in from its own directory, by its NEB value, so that
    # every entry in the distribution is the same kind of measurement
    results.append((REFERENCE_FRAME, REFERENCE_NEB, man.get(REFERENCE_FRAME, {})))

    for d in sorted(ensdir.glob("frame_*")):
        try:
            fr = int(d.name.split("_")[1])
        except (IndexError, ValueError):
            continue
        if fr == REFERENCE_FRAME:
            continue
        val, why = neb_barrier(d / "neb.out")
        if val is None:
            skipped.append((fr, why))
        else:
            results.append((fr, val, man.get(fr, {})))

    results.sort()

    print(f"{'frame':>7}{'ps':>9}{'barrier':>10}{'Arg90-O13':>12}"
          f"{'form C6-C1':>12}{'r':>9}")
    print("-" * 59)
    for fr, val, m in results:
        tag = "  <- fully characterised" if fr == REFERENCE_FRAME else ""
        print(f"{fr:>7}{m.get('ps', float('nan')):>9.0f}{val:>10.2f}"
              f"{m.get('arg90', float('nan')):>12.3f}"
              f"{m.get('form', float('nan')):>12.3f}"
              f"{m.get('r', float('nan')):>9.3f}{tag}")

    if skipped:
        print("\nnot included:")
        for fr, why in skipped:
            print(f"  frame {fr}: {why}")

    vals = [v for _, v, _ in results]
    if len(vals) < 2:
        print("\nonly one barrier available - no distribution to report yet")
        sys.exit(0)

    mean, sd, sem = stats(vals)
    print(f"\nENSEMBLE  (n = {len(vals)})")
    print(f"  mean          {mean:.2f} kcal/mol")
    print(f"  sd            {sd:.2f}")
    print(f"  sem           {sem:.2f}   <- quote comparisons against this")
    print(f"  range         {min(vals):.2f} to {max(vals):.2f} "
          f"(spread {max(vals)-min(vals):.2f})")

    ref = [v for f, v, _ in results if f == REFERENCE_FRAME][0]
    below = sum(1 for v in vals if v < ref)
    print(f"\n  frame {REFERENCE_FRAME}: {ref:.2f} by NEB-CI, "
          f"{REFERENCE_FULL:.2f} fully characterised")
    print(f"    sits above {below} of {len(vals)-1} other frames; "
          f"{(ref-mean)/sd:+.2f} sd from the ensemble mean" if sd == sd else "")
    print(f"    proxy validation on this frame: "
          f"{abs(REFERENCE_FULL-REFERENCE_NEB):.2f} kcal/mol")

    print(f"\n  vs experimental dH‡ {EXPERIMENT_DH}: "
          f"ensemble mean is {mean-EXPERIMENT_DH:+.2f} kcal/mol")
    print(f"    (a potential-energy barrier is not dH‡; zero-point and thermal")
    print(f"     corrections are needed before this is a like-for-like comparison)")

    have = [(v, m) for _, v, m in results if m]
    if len(have) >= 3:
        b = [v for v, _ in have]
        print("\n  CORRELATIONS across the ensemble")
        for key, label, expect in [
            ("arg90", "Arg90-O13 contact", "positive: looser contact, higher barrier"),
            ("form", "near-attack C6-C1", "positive: further from attack, higher barrier"),
        ]:
            x = [m[key] for _, m in have]
            r = pearson(x, b)
            print(f"    barrier vs {label:<20} r = {r:+.3f}   (expected {expect})")
        print("\n    A clear positive correlation with the Arg90 contact would say the")
        print("    frame-to-frame spread is electrostatic in origin, consistent with the")
        print("    literature attribution of TS stabilisation to that residue. Absence of")
        print("    one says the spread comes from elsewhere and is worth understanding")
        print("    before the mean is quoted.")

    with open(args.out, "w") as fh:
        fh.write("frame\tprod_ps\tbarrier_kcal\tArg90_O13\tform_C6_C1\tr\n")
        for fr, val, m in results:
            fh.write(f"{fr}\t{m.get('ps','')}\t{val:.3f}\t{m.get('arg90','')}"
                     f"\t{m.get('form','')}\t{m.get('r','')}\n")
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
