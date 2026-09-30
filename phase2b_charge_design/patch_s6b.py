#!/usr/bin/env python3
# Build s6b_inner_shell.pbs from s6_charge_tests.pbs, changing ONLY:
#   - the site list (constructed inner points instead of grid points)
#   - the tags used in the charge-magnitude scaling loop
#   - WORK, PBS -N, PBS -o
# run_sp, the level of theory, the control, and the entire analysis block are
# copied unchanged, so any difference in the result is due to position alone.
import sys, os
src = "s6_charge_tests.pbs"
dst = "s6b_inner_shell.pbs"
s = open(src).read()

# --- 1. PBS identity and work directory -------------------------------------
a = '#PBS -N '
i = s.index(a); j = s.index('\n', i)
s = s[:i] + '#PBS -N cm_s6b' + s[j:]
assert s.count('06_charge_tests') >= 1
s = s.replace('06_charge_tests', '06b_inner_shell')

# --- 2. the site list --------------------------------------------------------
old_start = "python3 - \"$DV\" > sites.tsv <<'PY'"
old_end   = "echo \"sites selected:\""
i = s.index(old_start); j = s.index(old_end)
new = '''python3 - "$DV" reactant.xyz > sites.tsv <<'PY'
import sys, math
H = 627.5094740631
# Parent: the most stabilising point of the production grid, on the 3.0 shell.
L = [l.split("\\t") for l in open(sys.argv[1]).read().splitlines() if l.strip()]
h = [x.strip() for x in L[0]]
ix, ic = h.index("idx"), h.index("dv_Eh")
xs, ys, zs, sh = h.index("x"), h.index("y"), h.index("z"), h.index("shell")
rows = [(int(r[ix]), float(r[xs]), float(r[ys]), float(r[zs]),
         float(r[sh]), float(r[ic])) for r in L[1:]]
rows.sort(key=lambda r: r[5])
parent = rows[0]

# Substrate geometry, to measure standoff from the vdW surface as step_2_1 does.
RV = {"H":1.20,"C":1.70,"N":1.55,"O":1.52,"S":1.80,"P":1.80}
G = open(sys.argv[2]).read().splitlines()
n = int(G[0].split()[0])
at = [(p[0], float(p[1]), float(p[2]), float(p[3]))
      for p in (l.split() for l in G[2:2+n])]
dist = lambda p,a: math.sqrt(sum((p[i]-a[i+1])**2 for i in range(3)))
standoff = lambda p: min(dist(p,a)-RV[a[0]] for a in at)
nearest  = lambda p: min(at, key=lambda a: dist(p,a)-RV[a[0]])

P = (parent[1], parent[2], parent[3])
a = nearest(P)
s0 = standoff(P)
sys.stderr.write("parent site idx %d: computed standoff %.3f A (production grid says %.1f)\\n"
                 % (parent[0], s0, parent[4]))
sys.stderr.write("  if these disagree, the vdW radii here differ from step_2_1 - STOP\\n")
v = [P[i]-a[i+1] for i in range(3)]
Lv = math.sqrt(sum(c*c for c in v)); u = [c/Lv for c in v]

print("tag\\tidx\\tx\\ty\\tz\\tshell\\tdv_Eh\\tdv_kcal")
print("%s\\t%d\\t%.6f\\t%.6f\\t%.6f\\t%.1f\\t%.9f\\t%.4f"
      % ("parent_3.0", parent[0], P[0], P[1], P[2], parent[4], parent[5], parent[5]*H))
# Inner points: same direction from the same nearest atom, smaller standoff.
# dv_Eh is a 1/r estimate from the parent. The analysis block's "error" column
# then measures how well that extrapolation holds - which is the point of this run.
for target in (2.5, 2.0, 1.5):
    d = s0 - target
    Q = [P[i]-u[i]*d for i in range(3)]
    est = parent[5] * dist(P,a) / dist(Q,a)
    print("%s\\t%d\\t%.6f\\t%.6f\\t%.6f\\t%.1f\\t%.9f\\t%.4f"
          % ("inner_%.1f" % target, -int(target*10), Q[0], Q[1], Q[2],
             target, est, est*H))
PY
'''
s = s[:i] + new + s[j:]

# --- 3. the scaling loop tags ------------------------------------------------
old = 'for tag in most_stabilising most_destabilising shell3_q1; do'
new = 'for tag in parent_3.0 inner_2.0 inner_1.5; do'
assert s.count(old) == 1, s.count(old)
s = s.replace(old, new)

open(dst, "w").write(s)
os.chmod(dst, 0o755)
print("wrote", dst)
