#!/usr/bin/env python3
# Phase 2 step 6.4 - enzyme spatial overlay. For each charge in a design, find the nearest enzyme
# catalytic residue of the SAME sign (rediscovery) and of ANY sign. Reports per design how much is
# rediscovered enzyme physics vs new space, quantifying whether the optimizer recreates or departs
# from the enzyme's actual electrostatics.
import sys, numpy as np

def load_design(npz, tag):
    z = np.load(npz, allow_pickle=True)
    q = z[f"{tag}_q"]; xyz = z["xyz"]; idx = z["idx"]
    sel = np.where(q!=0)[0]
    return [(int(idx[s]), int(q[s]), xyz[s]) for s in sel]

def load_residues(npz):
    z = np.load(npz, allow_pickle=True)
    return list(zip(z["names"], z["signs"].astype(int), z["centers"]))

def overlay(design, residues, match_radius=4.0):
    rows=[]
    for gidx,q,pos in design:
        same=[(np.linalg.norm(pos-c), nm, s) for nm,s,c in residues if s==q]
        anyr=[(np.linalg.norm(pos-c), nm, s) for nm,s,c in residues]
        same.sort(); anyr.sort()
        nearest_same = same[0] if same else (np.inf,"-",0)
        nearest_any  = anyr[0]
        rows.append(dict(gidx=gidx, q=q, pos=pos,
                         d_same=nearest_same[0], res_same=str(nearest_same[1]),
                         d_any=nearest_any[0], res_any=str(nearest_any[1]), sign_any=nearest_any[2],
                         rediscovered=nearest_same[0] <= match_radius))
    return rows

if __name__ == "__main__":
    npz = sys.argv[1] if len(sys.argv)>1 else "sol_discrete.npz"
    resnpz = sys.argv[2] if len(sys.argv)>2 else "enzyme_residues.npz"
    for match_radius in (4.0, 5.0):
        residues = load_residues(resnpz)
        print(f"\n########## match radius {match_radius} A ##########")
        print(f"enzyme active-site residues (within 12 A): {len(residues)}  "
              f"net sign balance {sum(s for _,s,_ in residues):+d}\n")
        for tag in ("K2_neutral","K2_free","K4_neutral","K6_neutral","K6_free"):
            try:
                design = load_design(npz, tag)
            except KeyError:
                continue
            rows = overlay(design, residues, match_radius)
            nred = sum(r["rediscovered"] for r in rows)
            print(f"=== {tag}: {len(design)} charges ===")
            for r in rows:
                t = "REDISCOVERED" if r["rediscovered"] else "new-space"
                print(f"  q={r['q']:+d} grid {r['gidx']:3d}: nearest same-sign {r['res_same']:>7} "
                      f"@{r['d_same']:5.2f}A | nearest any {r['res_any']:>7}({r['sign_any']:+d})"
                      f"@{r['d_any']:5.2f}A -> {t}")
            print(f"  -> {nred}/{len(rows)} coincide with same-sign residue; {len(rows)-nred} new space\n")
