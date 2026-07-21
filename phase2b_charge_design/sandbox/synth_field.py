#!/usr/bin/env python3
# Synthetic Dv fields for FORMULATION development only (no real data needed).
# Two generators: a single-image grid (endpoint Dv = V_TS - V_R) and a multi-image
# reaction path (per-image Dv_m relative to reactant). Replace with real-data loaders
# (dv_grid_dist.tsv; per-image orca_vpot output) when the consistent-region grid exists.
import math, random

def single_image_grid(n=252, seed=7):
    # points on 3 measured shells; Dv negative (cation-favorable) near an ether-O valley,
    # positive near a peak. Returns list of (id, x, y, z, dv_Ha).
    random.seed(seed); valley=(2.6,0.0,0.0); peak=(-2.2,1.5,0.0); pts=[]
    for i in range(n):
        shell=random.choice([3.2,4.2,5.4]); th,ph=random.uniform(0,math.pi),random.uniform(0,2*math.pi)
        x=shell*math.sin(th)*math.cos(ph); y=shell*math.sin(th)*math.sin(ph); z=shell*math.cos(th)
        dv=-0.022*math.exp(-math.dist((x,y,z),valley)/2.0)+0.018*math.exp(-math.dist((x,y,z),peak)/2.0)
        dv+=random.gauss(0,0.0004); pts.append((i,x,y,z,dv))
    return pts

def path_field(n=300, M=9, seed=5, sharp=True):
    # multi-image path: broad barrier plateau + a migrating (-)-charge valley (sharp field),
    # the regime where whole-path diverges from single-point. Returns (pts, b, dv) where
    # pts=(id,x,y,z), b[m]=bare image barrier (kcal/mol), dv[m][i]=Dv_m at point i (Ha).
    random.seed(seed); pts=[]
    for i in range(n):
        shell=random.choice([3.2,4.2,5.4]); th,ph=random.uniform(0,math.pi),random.uniform(0,2*math.pi)
        pts.append((i, shell*math.sin(th)*math.cos(ph), shell*math.sin(th)*math.sin(ph), shell*math.cos(th)))
    b=[2,8,16,18.0,18.5,18.0,16,8,2][:M]
    lam=0.8 if sharp else 1.8
    dv=[]
    for m in range(M):
        ang=math.pi*(m/(M-1)); valley=(3.0*math.cos(ang), 3.0*math.sin(ang), 0.0)
        dv.append([-0.016*math.exp(-math.dist((x,y,z),valley)/lam) for (_,x,y,z) in pts])
    return pts, b, dv
