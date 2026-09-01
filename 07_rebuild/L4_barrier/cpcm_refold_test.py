#!/usr/bin/env python3
"""Does CPCM eps=4 DRIVE toward the folded reactive geometry, or merely PRESERVE it from a
folded start? Take the vacuum-OPENED structure (C1-C6 ~4.39) and optimize it under CPCM.
- Re-folds to ~3.5  -> CPCM genuinely stabilizes the reactive fold (deep minimum)
- Stays open ~4.4   -> the fold was a shallow basin CPCM only preserves from a folded start
"""
import os, subprocess
import numpy as np
ORCADIR="/home/apps2/ORCA/6.0.1"
BASE=os.path.expanduser("~/07_rebuild")

def read_opt_xyz(path):
    L=open(path).read().splitlines(); n=int(L[0].split()[0]); return [l for l in L[2:2+n]]
def dist(lines,i,j):
    a=np.array([float(x) for x in lines[i].split()[1:4]]); b=np.array([float(x) for x in lines[j].split()[1:4]])
    return np.linalg.norm(a-b)

# the vacuum-opened structure from the control test = ctrl_2_bare_vac.xyz (C1-C6 ~4.39)
src="ctrl_2_bare_vac.xyz"
if not os.path.exists(src):
    print("!! %s not found - need the vacuum-opened geometry from the control test" % src); raise SystemExit
geom=read_opt_xyz(src)
print("starting from vacuum-opened structure: C1-C6=%.3f, O3-C4=%.3f" % (dist(geom,0,12), dist(geom,7,8)))

inp="refold.inp"
with open(inp,"w") as f:
    f.write("! B3LYP D3BJ def2-SVP def2/J RIJCOSX TightSCF Opt\n")
    f.write("%scf MaxIter 300 end\n%cpcm epsilon 4.0 end\n")
    f.write("* xyz -2 1\n")
    for l in geom: f.write(l+"\n")
    f.write("*\n")
subprocess.run(["%s/orca"%ORCADIR, inp], stdout=open("refold.out","w"),
               stderr=subprocess.STDOUT, env=dict(os.environ, PATH=ORCADIR+":"+os.environ["PATH"]))

if os.path.exists("refold.xyz"):
    fin=read_opt_xyz("refold.xyz")
    c=dist(fin,0,12)
    print("FINAL under CPCM: C1-C6=%.3f, O3-C4=%.3f" % (c, dist(fin,7,8)))
    print("="*60)
    if c < 3.7:
        print(">>> RE-FOLDED to %.3f -> CPCM DRIVES toward the reactive fold (deep minimum)." % c)
        print(">>> Your reactive geometry is robustly CPCM-stabilized, not just an inheritance artifact.")
    else:
        print(">>> STAYED OPEN at %.3f -> the fold is a shallow basin CPCM preserves from a folded" % c)
        print(">>> start but does NOT create from open. Your reactive geometry is partly inherited-start-dependent.")
else:
    out=open("refold.out").read()
    print("no refold.xyz - converged:", ("HURRAY" in out or "OPTIMIZATION HAS CONVERGED" in out))
    print("(check refold.out)")
