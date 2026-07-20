#!/usr/bin/env python3
# Monitor a surrogate relaxation and report the diagnostic that decides the method: how close a
# substrate atom gets to a surrogate atom, tracked across optimisation cycles.
#
# The point-charge run put a substrate oxygen 0.879 A from the designed +1, which is 0.64 A inside
# its own van der Waals radius. Real groups carry core electrons, so the approach should arrest at
# a salt-bridge distance instead: roughly 1.8-2.2 A for H...O, or 2.6-3.0 A nitrogen to oxygen.
# A contact still falling past 1.5 A is still collapse.
#
# Usage: python3 surrmon.py [stem]      default stem surr_ropt

import math
import os
import sys

N_SUB = 24
STEM = sys.argv[1] if len(sys.argv) > 1 else "surr_ropt"
BREAK, FORM = (7, 8), (0, 12)          # ether O3-C4 breaking, C1-C6 forming
REF = {"break": 1.478, "form": 3.181}


def read_frames(path):
    frames, cur = [], []
    for line in open(path, errors="ignore"):
        p = line.split()
        if len(p) == 1 and p[0].isdigit():
            if cur:
                frames.append(cur)
            cur = []
            continue
        if len(p) == 4:
            try:
                cur.append((p[0], float(p[1]), float(p[2]), float(p[3])))
            except ValueError:
                pass
    if cur:
        frames.append(cur)
    return [f for f in frames if len(f) >= N_SUB]


def dist(a, b):
    return math.sqrt((a[1]-b[1])**2 + (a[2]-b[2])**2 + (a[3]-b[3])**2)


def closest(frame):
    sub, surr = frame[:N_SUB], frame[N_SUB:]
    if not surr:
        return None
    best = None
    for s in sub:
        for t in surr:
            d = dist(s, t)
            if best is None or d < best[0]:
                best = (d, s[0], t[0])
    return best


out = STEM + ".out"
trj = STEM + "_trj.xyz"

if not os.path.isfile(out):
    print("%s not started" % STEM)
    sys.exit()

cycles = 0
energy = 0.0
converged = terminated = False
for line in open(out, errors="ignore"):
    if "GEOMETRY OPTIMIZATION CYCLE" in line:
        cycles += 1
    elif "FINAL SINGLE POINT ENERGY" in line:
        try:
            energy = float(line.split()[-1])
        except ValueError:
            pass
    elif "HURRAY" in line or "THE OPTIMIZATION HAS CONVERGED" in line:
        converged = True
    elif "ORCA TERMINATED NORMALLY" in line:
        terminated = True

state = "done" if (terminated and converged) else ("ended unconverged" if terminated else "running")
print("%s: %s, %d cycles, last energy %.6f" % (STEM, state, cycles, energy))

if not os.path.isfile(trj):
    print("no trajectory yet")
    sys.exit()

frames = read_frames(trj)
if not frames:
    print("trajectory empty")
    sys.exit()

print("\ncycle   closest substrate-surrogate   break    form")
step = max(1, len(frames)//12)
for i, f in enumerate(frames):
    if not (i == 0 or i == len(frames)-1 or i % step == 0):
        continue
    c = closest(f)
    if c is None:
        continue
    b = dist(f[BREAK[0]], f[BREAK[1]])
    fm = dist(f[FORM[0]], f[FORM[1]])
    mark = "   <- inside 1.5 A" if c[0] < 1.5 else ""
    print("%5d   %s...%s %6.3f              %6.3f  %6.3f%s" % (i, c[1], c[2], c[0], b, fm, mark))

c = closest(frames[-1])
print("\nlatest: closest %s...%s at %.3f A" % (c[1], c[2], c[0]))
print("  point-charge run collapsed to 0.879 A (O onto +1)")
print("  salt-bridge range is 1.8-2.2 A for H...O, 2.6-3.0 A for N...O")
print("  reference reactant: break %.3f  form %.3f" % (REF["break"], REF["form"]))
