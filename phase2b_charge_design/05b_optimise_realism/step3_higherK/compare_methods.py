#!/usr/bin/env python3
# Calibration: compare DFT (B3LYP-D3BJ/def2-SVP) vs r2SCAN-3c relaxed geometries for the SAME designs.
# The two runs differ ONLY in QM method, so any difference in relaxed C1-C6 / connectivity is a pure
# method effect. This is the DESIGN_DECISIONS Tier-2a calibration gate on THIS system.
#   keep r2SCAN-3c as cheap geometry driver if geometries agree (C1-C6 within ~0.2 A, same connectivity)
#   demote to geometry-only / abandon if they disagree (spurious collapse or different near-attack verdict).
#
# Usage: compare_methods.py  [run_dir_root]
#   expects <root>/run_<name>/<name>_relaxed.xyz         (DFT)
#       and <root>/run_<name>_r2scan/<name>_r2scan_relaxed.xyz  (r2SCAN)  [or flat, auto-detected]

import sys, os, math

C1,O3,C4,C6 = 0,7,8,12
def load(path):
    L=open(path).read().splitlines(); n=int(L[0].split()[0])
    els=[]; xyz=[]
    for i in range(n):
        t=L[2+i].split(); els.append(t[0]); xyz.append([float(v) for v in t[1:4]])
    return els,xyz
def d(a,b): return math.sqrt(sum((a[i]-b[i])**2 for i in range(3)))
COV={"C":0.76,"N":0.71,"O":0.66,"H":0.31}
def connected(els,xyz,scale=1.3):
    N=min(24,len(xyz)); adj={i:[] for i in range(N)}
    for i in range(N):
        for j in range(i+1,N):
            if d(xyz[i],xyz[j])<=scale*(COV.get(els[i],.77)+COV.get(els[j],.77)):
                adj[i].append(j); adj[j].append(i)
    seen={0}; st=[0]
    while st:
        u=st.pop()
        for v in adj[u]:
            if v not in seen: seen.add(v); st.append(v)
    return len(seen)==N

def find(root,name,r2):
    # try per-job subdir then flat
    suff="_r2scan_relaxed.xyz" if r2 else "_relaxed.xyz"
    sub =("run_%s_r2scan"%name) if r2 else ("run_%s"%name)
    for p in [os.path.join(root,sub,name+suff), os.path.join(root,name+suff)]:
        if os.path.exists(p): return p
    return None

def geom(path):
    els,xyz=load(path)
    return d(xyz[C1],xyz[C6]), d(xyz[O3],xyz[C4]), connected(els,xyz)

def main(root="."):
    names=["K6_apot","K6_aforce","K10_apot","K10_aforce"]
    print("Calibration: DFT (B3LYP-D3BJ/def2-SVP) vs r2SCAN-3c, same designs, method-only difference")
    print("window = near-attack C1-C6 [2.6-3.5] A ; K2 refs: a_pot 5.28 / a_force 4.06 (DFT)\n")
    print("%-11s | %-22s | %-22s | %s"%("design","DFT  C1-C6 O3-C4 conn","r2SCAN C1-C6 O3-C4 conn","dC1C6  verdict"))
    print("-"*92)
    rows=[]
    for nm in names:
        pd=find(root,nm,False); pr=find(root,nm,True)
        if not pd or not pr:
            print("%-11s | %s"%(nm, "MISSING: DFT=%s r2SCAN=%s"%(bool(pd),bool(pr)))); continue
        c1d,o3d,cd=geom(pd); c1r,o3r,cr=geom(pr)
        dd=abs(c1d-c1r)
        # verdicts: do both agree on in/out of window?
        ind=2.6<=c1d<=3.5; inr=2.6<=c1r<=3.5
        agree_win = "AGREE" if ind==inr else "DISAGREE(window!)"
        agree_geo = "geom OK" if dd<=0.2 else ("geom~%.2f"%dd if dd<=0.5 else "geom DIVERGE %.2f"%dd)
        print("%-11s | %5.3f  %5.3f  %-4s | %5.3f  %5.3f  %-4s | %.3f  %s / %s"%(
            nm,c1d,o3d,"1mol" if cd else "FRAG",c1r,o3r,"1mol" if cr else "FRAG",dd,agree_geo,agree_win))
        rows.append((nm,c1d,c1r,dd,ind,inr,cd,cr))
    if rows:
        mx=max(r[3] for r in rows); mean=sum(r[3] for r in rows)/len(rows)
        winagree=sum(1 for r in rows if r[4]==r[5]); conagree=sum(1 for r in rows if r[6]==r[7])
        print("\nCALIBRATION SUMMARY (%d designs compared):"%len(rows))
        print("  C1-C6 |DFT - r2SCAN|: mean %.3f A, max %.3f A"%(mean,mx))
        print("  near-attack window verdict agreement: %d/%d"%(winagree,len(rows)))
        print("  connectivity agreement:               %d/%d"%(conagree,len(rows)))
        gate = mx<=0.2 and winagree==len(rows) and conagree==len(rows)
        print("  GATE (geometry): %s"%("PASS - r2SCAN-3c usable as cheap geometry driver on this system"
              if gate else "FAIL/PARTIAL - r2SCAN-3c geometry diverges; keep DFT for geometry (see per-design)"))

if __name__=="__main__":
    main(sys.argv[1] if len(sys.argv)>1 else ".")
