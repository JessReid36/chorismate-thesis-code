#!/usr/bin/env python3
# Step 0 (PLAN.md, report section B(a)): linear per-site coefficients for the preorganisation term.
#
# PLAN definition (primary): a_s = V_s(C1) - V_s(C6)
#   = the electrostatic POTENTIAL a unit +1 charge at candidate site s induces across the C1<->C6 pair.
#   Same differential-potential construction as Dv = V_TS - V_R, but ACROSS the two reacting carbons
#   instead of across two path states. Linear in the source charge (Coulomb superposition), so the total
#   distorting drive of a design = sum_s q_s * a_s is LINEAR in the placement variables -> drops into the
#   MILP as a hard cap or penalty WITHOUT breaking the certificate.
#
#   Physical reading: V_s(C1) - V_s(C6) is the potential DROP from C1 to C6 due to site s. A large drop
#   means the field pushes (developing) negative charge from one carbon toward the other / separates the
#   charges along the bond - the electrostatic signature of driving the reacting carbons apart rather than
#   holding them in the near-attack geometry.
#
# Cross-check column (secondary): a_force = [F(C6)-F(C1)] . u_hat, the mechanical force-imbalance along the
#   C1->C6 axis (option in the plan's discussion). Sign agreement between a_pot and a_force validates the mode.
#
# k=1, geometry in Angstrom. Only relative magnitudes/signs matter; we normalise max|a_pot|=1.

import sys, math

def load_xyz(path):
    L=open(path).read().splitlines(); n=int(L[0].split()[0])
    return [[float(t) for t in L[2+i].split()[1:4]] for i in range(n)]

def load_grid(path):
    rows=[]
    with open(path) as fh:
        next(fh)
        for ln in fh:
            p=ln.rstrip("\n").split("\t")
            rows.append((int(p[0]), float(p[1]), float(p[2]), float(p[3]), p[4]))  # idx,x,y,z,shell
    return rows

def dist(a,b): return math.sqrt(sum((a[i]-b[i])**2 for i in range(3)))

def main(reactant_xyz, dv_grid, out):
    R=load_xyz(reactant_xyz); C1=R[0]; C6=R[12]
    u=[C6[i]-C1[i] for i in range(3)]; Lax=math.sqrt(sum(c*c for c in u)); u=[c/Lax for c in u]
    grid=load_grid(dv_grid)
    rows=[]
    for idx,x,y,z,shell in grid:
        s=[x,y,z]
        r1=dist(C1,s); r6=dist(C6,s)
        # PRIMARY: potential from a unit +1 at s, evaluated at C1 and C6 (k=1)
        a_pot = (1.0/r1) - (1.0/r6)              # V_s(C1) - V_s(C6)
        # SECONDARY cross-check: force imbalance along bond axis
        def F(at):
            d=[at[i]-s[i] for i in range(3)]; r=math.sqrt(sum(c*c for c in d))
            inv=1.0/(r*r*r) if r>1e-6 else 0.0
            return [d[i]*inv for i in range(3)]
        F1=F(C1); F6=F(C6); a_force=sum((F6[i]-F1[i])*u[i] for i in range(3))
        rows.append((idx,shell,a_pot,a_force))
    amax=max(abs(r[2]) for r in rows) or 1.0
    with open(out,"w") as fh:
        fh.write("idx\tshell\ta_pot\ta_pot_norm\ta_force_xcheck\n")
        for idx,shell,ap,af in rows:
            fh.write("%d\t%s\t%.6e\t%+.4f\t%.6e\n"%(idx,shell,ap,ap/amax,af))
    # summary + sign-agreement cross-check
    agree=sum(1 for _,_,ap,af in rows if (ap>0)==(af>0))
    pos=[ap for _,_,ap,_ in rows if ap>0]; neg=[ap for _,_,ap,_ in rows if ap<0]
    print("sites=%d  C1-C6 axis=%.3f A  (C1 closer=> a_pot>0)"%(len(rows),Lax))
    print("a_pot>0 (drop C1->C6, distorting-signature): %d sites"%len(pos))
    print("a_pot<0 (drop C6->C1):                        %d sites"%len(neg))
    print("sign agreement a_pot vs a_force cross-check: %d/%d (%.0f%%)"%(agree,len(rows),100*agree/len(rows)))
    print("wrote %s"%out)

if __name__=="__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])
