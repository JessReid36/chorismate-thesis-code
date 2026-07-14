#!/usr/bin/env python3
# Phase 2 step 6.5 - unconstrained (grid) vs enzyme-constrained (residue positions) design.
#   (1) UNCONSTRAINED : full grid, sign-free -> REUSED from sol_discrete.npz (Stage 3.4), not re-solved.
#   (2) ENZYME sign-LOCKED : candidates = real residue positions, each locked to its real sign.
#   (3) ENZYME sign-FREE   : residue positions, any sign -> isolates positions vs signs as the limit.
# Only (2),(3) solved here (12 candidates) -> instant. unconstrained<<sign-locked => method beats the
# enzyme's arrangement; sign-free<<sign-locked => the enzyme's committed signs (binding cations) cost it.
import sys, numpy as np
from mip import Model, xsum, minimize, BINARY, OptimizationStatus
HART2KCAL = 627.509

def load_residue_dv(npz):
    z=np.load(npz,allow_pickle=True)
    return ([str(n) for n in z["names"]], z["signs"].astype(int), z["centers"], z["dv"])

def unconstrained_from_solution(npz):
    z=np.load(npz,allow_pickle=True); out={}
    for key in z.files:
        if key.endswith("_ddE"):
            tag=key[:-4]; K=int(tag.split("_")[0][1:]); val=float(z[key])
            out.setdefault(K, []).append((tag,val))
    return {K:min(v,key=lambda t:t[1]) for K,v in out.items()}

def miqp_residues(dv, K, sign_lock):
    n=len(dv); md=Model(sense=minimize); md.verbose=0
    if sign_lock is not None:
        a=[md.add_var(var_type=BINARY) for _ in range(n)]
        md += xsum(a)<=K
        md.objective = xsum(dv[i]*int(sign_lock[i])*a[i] for i in range(n))
        st=md.optimize()
        q=np.array([int(sign_lock[i]) if a[i].x and a[i].x>0.5 else 0 for i in range(n)])
    else:
        p=[md.add_var(var_type=BINARY) for _ in range(n)]
        m=[md.add_var(var_type=BINARY) for _ in range(n)]
        for i in range(n): md += p[i]+m[i]<=1
        md += xsum(p[i]+m[i] for i in range(n))<=K
        md.objective = xsum(dv[i]*(p[i]-m[i]) for i in range(n))
        st=md.optimize()
        q=np.array([(1 if p[i].x and p[i].x>0.5 else 0)-(1 if m[i].x and m[i].x>0.5 else 0) for i in range(n)])
    return q, (md.objective_value or 0.0), st

if __name__=="__main__":
    sol_npz=sys.argv[1] if len(sys.argv)>1 else "sol_discrete.npz"
    res_npz=sys.argv[2] if len(sys.argv)>2 else "residue_dv.npz"
    names,signs,centers,dv=load_residue_dv(res_npz)
    uncon=unconstrained_from_solution(sol_npz)
    print(f"enzyme residue candidates: {len(dv)}  "
          f"(cations {int((signs>0).sum())}, anions {int((signs<0).sum())}, neutral {int((signs==0).sum())})\n")
    print(f"{'K':>3} {'unconstrained grid':>26} {'enzyme sign-locked':>19} {'enzyme sign-free':>17}")
    for K in sorted(uncon):
        tag,ug=uncon[K]; Kr=min(K,len(dv))
        ql,ol,stl=miqp_residues(dv,Kr,signs)
        qf,of,stf=miqp_residues(dv,Kr,None)
        s_ok="OPT" if stl==OptimizationStatus.OPTIMAL and stf==OptimizationStatus.OPTIMAL else "chk"
        print(f"{K:3d} {ug*HART2KCAL:14.3f} [{tag:11}] {ol*HART2KCAL:16.3f} {of*HART2KCAL:14.3f}   {s_ok}")
    print("\n(kcal/mol; more negative = better barrier-lowering)")
    Kmax=max(uncon); ql,_,_=miqp_residues(dv,min(Kmax,len(dv)),signs); qf,_,_=miqp_residues(dv,min(Kmax,len(dv)),None)
    print(f"at K={min(Kmax,len(dv))}:")
    print("  sign-locked activates:", [f"{names[i]}({ql[i]:+d})" for i in range(len(dv)) if ql[i]!=0])
    print("  sign-free  activates:", [f"{names[i]}({qf[i]:+d})" for i in range(len(dv)) if qf[i]!=0])
