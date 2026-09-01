#!/usr/bin/env python3
"""PHASE 1: is a linear-response electrostatic surrogate (DTSS) able to RANK charge designs?
For each design, Delta = sum_a q_a * [V_TS(r_a) - V_S(r_a)], where V is the EXACT DFT
electrostatic potential of the BARE substrate (common reaction path, per DTSS), evaluated at
that design's external-charge positions r_a. Predicted barrier = 17.47 + Delta.
Success = the RANKING/correlation vs the committed barriers holds (NOT absolute match).
This tests the make-or-break premise of the whole deterministic-search plan, on data we ALREADY have.
"""
import os, sys, subprocess, re
import numpy as np
HK = 627.5094740631  # Eh -> kcal/mol
ANG2BOHR = 1.8897259886

ORCADIR = "/home/apps2/ORCA/6.0.1"
SIMPLE = "! B3LYP D3BJ def2-SVP def2/J RIJCOSX TightSCF KeepDens"
BLOCKS = "%cpcm epsilon 4.0 end\n%scf MaxIter 300 end"
NCORES = 8
BARE_BARRIER = 17.47

# --- the labeled designs: (name, coords_file, charges_file, committed_barrier) ---
# coords_file has 24 substrate atoms + N external-charge positions (as XYZ, element is a dummy for the charges)
# charges_file is a comma list: 24 zeros (substrate, not point charges here) + N external charge magnitudes
BASE = os.path.expanduser("~/07_rebuild")
DESIGNS = [
    ("Arg90",     BASE+"/L4_barrier/design_cradleArg90_coords.xyz",
                  BASE+"/L4_barrier/design_cradleArg90_charges.txt", 6.40),
    ("OEEF",      BASE+"/L1_oracle/design_oeefs_q0p50_fb0p013_m0p010_charges.txt", None, 46.37),  # coords resolved below
    ("capacitor", None, None, 29.77),  # resolved below
]

def read_xyz_atoms(path):
    L = open(path).read().splitlines()
    n = int(L[0].split()[0])
    atoms = []
    for ln in L[2:2+n]:
        p = ln.split()
        atoms.append((p[0], float(p[1]), float(p[2]), float(p[3])))
    return atoms

def read_charges(path):
    return [float(x) for x in open(path).read().strip().split(",")]

def bare_geom(which):
    """Return the 24-atom bare substrate geometry lines (element x y z) for reactant or TS."""
    if which == "reactant":
        L = open(BASE+"/L0_groundtruth/L0_react_opt.xyz").read().splitlines()
        n = int(L[0].split()[0]); return [l for l in L[2:2+n]]
    else:  # TS: inline in L0_ts_sp.inp between '* xyz -2 1' and '*'
        L = open(BASE+"/L0_groundtruth/L0_ts_sp.inp").read().splitlines()
        out=[]; grab=False
        for l in L:
            s=l.strip()
            if s.startswith("* xyz"): grab=True; continue
            if grab:
                if s=="*" or s=="": break
                p=s.split()
                if len(p)>=4 and p[0][0].isalpha(): out.append("  ".join(p[:4]))
        return out

def run_scf_and_vpot(tag, geom_lines, points_xyz):
    """SCF on bare substrate (charge -2), then orca_vpot to get V at the given points (Bohr).
    Returns array of V (atomic units) at each point."""
    inp = "scf_%s.inp" % tag
    with open(inp, "w") as f:
        f.write(SIMPLE+"\n"+BLOCKS+"\n%pal nprocs "+str(NCORES)+" end\n")
        f.write("* xyz -2 1\n")
        for l in geom_lines:
            f.write(l+"\n")
        f.write("*\n")
    # run ORCA
    subprocess.run(["%s/orca"%ORCADIR, inp], stdout=open("scf_%s.out"%tag,"w"),
                   stderr=subprocess.STDOUT, env=dict(os.environ, PATH=ORCADIR+":"+os.environ["PATH"]))
    gbw = "scf_%s.gbw"%tag; dens = "scf_%s.densities"%tag
    # points file for orca_vpot: first line = number of points, then x y z in BOHR
    pf = "vpot_%s_points.xyz"%tag
    with open(pf,"w") as f:
        f.write("%d\n"%len(points_xyz))
        for (x,y,z) in points_xyz:
            f.write("%.10f %.10f %.10f\n"%(x*ANG2BOHR, y*ANG2BOHR, z*ANG2BOHR))
    outf = "vpot_%s_out.dat"%tag
    scfp = "scf_%s.scfp" % tag  # member name inside the .densities container (need not exist standalone)
    base = "scf_%s" % tag       # DensContainerBase = 5th arg, lets orca_vpot find the .scfp member
    r = subprocess.run(["%s/orca_vpot"%ORCADIR, gbw, scfp, pf, outf, base],
                   stdout=open("vpot_%s.log"%tag,"w"), stderr=subprocess.STDOUT,
                   env=dict(os.environ, PATH=ORCADIR+":"+os.environ["PATH"]))
    # parse orca_vpot output: it writes V per point (atomic units)
    V=[]
    for line in open(outf):
        s=line.split()
        if len(s) < 4:   # skip the header count line (single number) and blanks
            continue
        try:
            V.append(float(s[-1]))   # V is the last column of "x y z V"
        except ValueError:
            pass
    return np.array(V)

def resolve_design_files():
    """Explicit clean design files (24+N atoms), NOT run_* trajectory dirs.
    Verified atom counts: Arg90 27 (24+3), OEEF 30 (24+6), capacitor 76 (24+52)."""
    return [
        ("Arg90",     BASE+"/L4_barrier/design_cradleArg90_coords.xyz",
                      BASE+"/L4_barrier/design_cradleArg90_charges.txt", 6.40),
        ("OEEF",      BASE+"/L1_oracle/design_oeefs_q0p50_fb0p013_m0p010_coords.xyz",
                      BASE+"/L1_oracle/design_oeefs_q0p50_fb0p013_m0p010_charges.txt", 46.37),
        ("capacitor", BASE+"/L1_oracle/design_cap_f0p020_n5_coords.xyz",
                      BASE+"/L1_oracle/design_cap_f0p020_n5_charges.txt", 29.77),
    ]

def main():
    designs = resolve_design_files()
    print("=== resolved designs ===")
    for name,cx,cq,b in designs:
        print("  %-10s coords=%s" % (name, os.path.basename(cx) if cx else "MISSING"))
    # gather each design's external-charge positions (atoms 24..end) + magnitudes (24..end)
    dz = {}
    for name,cx,cq,b in designs:
        atoms = read_xyz_atoms(cx)
        q = read_charges(cq)
        ext_pos = [(x,y,z) for (_,x,y,z) in atoms[24:]]
        ext_q   = q[24:]
        assert len(ext_pos)==len(ext_q), "%s: %d pos vs %d q"%(name,len(ext_pos),len(ext_q))
        dz[name] = (ext_pos, ext_q, b)
        print("  %-10s : %d external charges, sum q = %+.3f, barrier %.2f" % (name, len(ext_q), sum(ext_q), b))
    # union of ALL charge positions -> one vpot evaluation per state
    all_points=[]; index_map={}
    for name,(pos,q,b) in dz.items():
        index_map[name]=[]
        for p in pos:
            index_map[name].append(len(all_points)); all_points.append(p)
    print("\n=== running bare-substrate SCF + vpot at %d charge points (reactant, then TS) ===" % len(all_points))
    V_S  = run_scf_and_vpot("react", bare_geom("reactant"), all_points)
    V_TS = run_scf_and_vpot("ts",    bare_geom("ts"),       all_points)
    if len(V_S)!=len(all_points) or len(V_TS)!=len(all_points):
        print("!! vpot returned %d/%d values (expected %d) -- parsing needs a look"%(len(V_S),len(V_TS),len(all_points)))
        print("   (inspect vpot_react_out.dat / vpot_ts_out.dat)")
        return
    dV = V_TS - V_S  # exact DFT potential difference at each point (atomic units)
    print("\n=== DTSS surrogate: Delta = sum_a q_a*(V_TS - V_S), predicted barrier = 17.47 + Delta ===")
    print("%-10s %12s %12s %12s   %s" % ("design","Delta(kcal)","pred_barr","true_barr","dir"))
    rows=[]
    for name,(pos,q,b) in dz.items():
        idx = index_map[name]
        delta_au = sum(q[i]*dV[idx[i]] for i in range(len(q)))
        delta_kcal = delta_au*HK
        pred = BARE_BARRIER + delta_kcal
        direction = "LOWER" if delta_kcal<0 else "raise"
        rows.append((name, delta_kcal, pred, b))
        print("%-10s %12.2f %12.2f %12.2f   %s" % (name, delta_kcal, pred, b, direction))
    # add the bare reference (Delta=0 by construction)
    rows.append(("bare", 0.0, BARE_BARRIER, BARE_BARRIER))
    print("%-10s %12.2f %12.2f %12.2f   %s" % ("bare",0.0,BARE_BARRIER,BARE_BARRIER,"ref"))
    # ranking check
    print("\n=== RANKING CHECK ===")
    true_order = sorted(rows, key=lambda r: r[3])
    pred_order = sorted(rows, key=lambda r: r[2])
    print("true ranking (by committed barrier):  ", " < ".join(r[0] for r in true_order))
    print("surrogate ranking (by predicted):     ", " < ".join(r[0] for r in pred_order))
    # Spearman-ish: do the orders match?
    tnames=[r[0] for r in true_order]; pnames=[r[0] for r in pred_order]
    match = tnames==pnames
    # correlation
    tb=np.array([r[3] for r in rows]); pb=np.array([r[2] for r in rows])
    if len(rows)>2 and tb.std()>0 and pb.std()>0:
        corr=np.corrcoef(tb,pb)[0,1]
        print("Pearson correlation (pred vs true barrier): %.3f" % corr)
    print("\n>>> RANK ORDER MATCHES: %s" % match)
    if match:
        print(">>> PASS: linear-response electrostatic surrogate reproduces the ranking.")
        print(">>> Phase 2 (MILP on the bare path) has a sound scoring function.")
    else:
        print(">>> Ranking does NOT fully match. Check WHICH designs misorder:")
        print(">>> If only the violent field designs (capacitor/OEEF) misorder, that is the")
        print(">>> geometric-ceiling effect (fields violate the fixed-path premise) -- still informative.")

if __name__ == "__main__":
    main()
