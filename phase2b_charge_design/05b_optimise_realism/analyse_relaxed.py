#!/usr/bin/env python3
# Step 3 analysis (PLAN.md): read a relaxed 24-atom-substrate xyz (or a >24 file with QM substrate first)
# and report the diagnostics that decide whether a preorganised design stayed reactive.
# Reusable across K2/K4/any Step-3 design so the verdict is identical and one command.
#
# Reacting atoms (0-based): C1=0, ether O3=7, C4=8, C6=12.
# Reactive window (near-attack): C1-C6 in [2.6, 3.5] A (TS ~2.5). Bare substrate reactant C1-C6=3.124.
# Reference points (committed):
#   bare-substrate reactant : O3-C4 1.472  C1-C6 3.124
#   K2 mu=0 field-relaxed    : O3-C4 1.434  C1-C6 5.004  (DISTORTED open - the control)
#   TS                       : O3-C4 2.111  C1-C6 2.526

import sys, math

def load(path):
    L=open(path).read().splitlines(); n=int(L[0].split()[0])
    els=[]; xyz=[]
    for i in range(n):
        t=L[2+i].split(); els.append(t[0]); xyz.append([float(v) for v in t[1:4]])
    return els,xyz

def d(a,b): return math.sqrt(sum((a[i]-b[i])**2 for i in range(3)))

# covalent radii (A) for a simple single-molecule connectivity check on the 24-atom substrate
COV={"C":0.76,"N":0.71,"O":0.66,"H":0.31}
def connected(els,xyz,scale=1.3):
    # build adjacency on the FIRST 24 atoms (the QM substrate); BFS from atom 0; all reached => 1 molecule
    N=min(24,len(xyz)); adj={i:[] for i in range(N)}
    for i in range(N):
        for j in range(i+1,N):
            rij=d(xyz[i],xyz[j]); cut=scale*(COV.get(els[i],0.77)+COV.get(els[j],0.77))
            if rij<=cut: adj[i].append(j); adj[j].append(i)
    seen={0}; stack=[0]
    while stack:
        u=stack.pop()
        for v in adj[u]:
            if v not in seen: seen.add(v); stack.append(v)
    return len(seen)==N, len(seen), N

def main(path):
    els,xyz=load(path)
    C1,O3,C4,C6 = xyz[0],xyz[7],xyz[8],xyz[12]
    dC1C6=d(C1,C6); dO3C4=d(O3,C4); dC4C6=d(C4,C6)
    single,reached,N=connected(els,xyz)
    inwin = 2.6<=dC1C6<=3.5
    print("file: %s"%path)
    print("  C1-C6 (forming bond) : %.3f A   %s"%(dC1C6,
          "IN reactive window [2.6-3.5]  <-- preorganised" if inwin else
          ("OPEN (>3.5) - distorted like mu=0 control" if dC1C6>3.5 else "SHORT (<2.6) - past near-attack")))
    print("  O3-C4 (breaking bond): %.3f A   %s"%(dO3C4, "intact ether (reactant-like)" if dO3C4<1.7 else "elongated/broken"))
    print("  C4-C6                : %.3f A"%dC4C6)
    print("  substrate connectivity: %s (%d/%d atoms in one molecule)"%(
          "SINGLE molecule" if single else "FRAGMENTED", reached,N))
    print()
    print("  vs references:  bare reactant C1-C6 3.124 | mu=0 control 5.004 (distorted) | TS 2.526")
    # verdict
    if single and inwin:
        print("  VERDICT: PREORGANISED + intact -> preorg penalty WORKS on this design (a_pot validated).")
    elif single and dC1C6>3.5:
        print("  VERDICT: still OPEN -> penalty did NOT recover near-attack; revisit mode (a_force). [falsifiable arbiter]")
    elif not single:
        print("  VERDICT: FRAGMENTED -> invalid endpoint; do not read barrier. Investigate (over-polarisation?).")
    else:
        print("  VERDICT: C1-C6 short of window - inspect (possible over-compression / product-like).")

if __name__=="__main__":
    if len(sys.argv)<2:
        print("usage: analyse_relaxed.py <relaxed.xyz>"); sys.exit(1)
    main(sys.argv[1])
