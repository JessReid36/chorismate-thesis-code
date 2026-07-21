#!/usr/bin/env python3
# Extract the substrate geometry (last frame of a reduced-region QM-region trajectory) into a
# clean .xyz for phase2b/01_geometry, with self-checks. Replaces the old manual 01_substrate swap.
# Self-checks: 24 atoms; elements at the reacting indices (C1=C, O3=O, C4=C, C6=C); and the
# reacting distances match the verified Phase-1 values within tolerance -> guards atom order+frame.
# Pure-Python, pre-3.8. No numpy.
import sys, math
IDX = {"C1":0, "O3":7, "C4":8, "C6":12}
EXPECT_ELEM = {"C1":"C", "O3":"O", "C4":"C", "C6":"C"}
def last_frame(path):
    L=open(path).read().splitlines(); i=0; frames=[]
    while i < len(L):
        try: n=int(L[i].split()[0])
        except (ValueError, IndexError): break
        atoms=[]
        for ln in L[i+2:i+2+n]:
            p=ln.split(); atoms.append((p[0], float(p[1]), float(p[2]), float(p[3])))
        frames.append(atoms); i+=2+n
    if not frames: raise SystemExit("FAIL no frames parsed in %s"%path)
    return frames[-1], len(frames)
def dist(a,b): return math.sqrt(sum((a[k+1]-b[k+1])**2 for k in range(3)))
def main():
    trj, label, out, exp_break, exp_form, tol = sys.argv[1:7]
    exp_break=float(exp_break); exp_form=float(exp_form); tol=float(tol)
    atoms, nframes = last_frame(trj)
    if len(atoms)!=24: raise SystemExit("FAIL %s: %d atoms (expected 24)"%(label,len(atoms)))
    for name,idx in IDX.items():
        if atoms[idx][0]!=EXPECT_ELEM[name]:
            raise SystemExit("FAIL %s: atom %d is %s, expected %s (%s) - order wrong"%(label,idx,atoms[idx][0],EXPECT_ELEM[name],name))
    brk=dist(atoms[IDX["O3"]],atoms[IDX["C4"]]); frm=dist(atoms[IDX["C1"]],atoms[IDX["C6"]])
    db=abs(brk-exp_break); df=abs(frm-exp_form)
    status="OK" if (db<=tol and df<=tol) else "MISMATCH"
    print("%-9s frames=%d  break(O3-C4)=%.3f (exp %.3f, d%.3f)  form(C1-C6)=%.3f (exp %.3f, d%.3f)  [%s]"%(label,nframes,brk,exp_break,db,frm,exp_form,df,status))
    if status!="OK": raise SystemExit("FAIL %s: reacting distances off by > tol (%.3f)"%(label,tol))
    with open(out,"w") as fh:
        fh.write("24\n%s substrate (extracted from %s, last of %d frames)\n"%(label,trj,nframes))
        for el,x,y,z in atoms: fh.write("%-2s %14.8f %14.8f %14.8f\n"%(el,x,y,z))
    print("  wrote %s"%out)
if __name__=="__main__": main()
